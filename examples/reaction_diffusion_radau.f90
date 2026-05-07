! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot: GPT-5
!
! reaction_diffusion_radau.f90 - Fully implicit Runge-Kutta for a 1-D
! reaction-diffusion travelling wave
!
! Solves
!
!   u_t = u_xx + (1 - u) u^2    on  x in (0, 10),  t in (0, T]
!
! with Dirichlet boundary conditions and initial data imposed from the
! analytical travelling-wave solution
!
!   u(x, t) = 1 / (1 + exp(v (x - v t))),   v = sqrt(1/2).
!
! Time integration uses the 2-stage Radau IIA method:
!
!     [ 5/12  -1/12 ]
! A = [             ],   b = [ 3/4, 1/4 ],   c = [ 1/3, 1 ].
!     [ 3/4    1/4  ]
!
! The stage equations are fully coupled, so every Newton step solves a
! 2-by-2 block sparse linear system with Y12M.  Since the Jacobian pattern is
! constant across all Newton steps and time steps, this example demonstrates
! Case (iii): same sparsity, changing values.
!
! Usage
! -----
!   reaction_diffusion_radau [--help] [N] [nsteps] [T] [output_file]
!
!     N           total grid size (>= 3, default 41)
!     nsteps      number of time steps (default 20)
!     T           final time (default 1.0)
!     output_file output file name (default: reaction_diffusion_radau.dat)
!
module reaction_diffusion_radau_solver
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, real64
   implicit none
   private
   public :: run, dp

   integer, parameter :: dp = real64
   integer, parameter :: nstage = 2
   real(dp), parameter :: x_left = 0.0_dp
   real(dp), parameter :: x_right = 10.0_dp
   real(dp), parameter :: wave_speed = sqrt(0.5_dp)
   real(dp), parameter :: rk_a(nstage, nstage) = reshape([ &
      5.0_dp / 12.0_dp, 3.0_dp / 4.0_dp, &
      -1.0_dp / 12.0_dp, 1.0_dp / 4.0_dp ], [nstage, nstage])
   real(dp), parameter :: rk_b(nstage) = [3.0_dp / 4.0_dp, 1.0_dp / 4.0_dp]
   real(dp), parameter :: rk_c(nstage) = [1.0_dp / 3.0_dp, 1.0_dp]
   real(dp), parameter :: newton_tol = 1.0e-10_dp
   real(dp), parameter :: y12m_growth_limit = 16.0_dp
   integer, parameter :: y12m_workspace_factor = 20
   integer, parameter :: max_newton = 12

   type :: reaction_diffusion_problem
      integer :: ndof = 0
      real(dp) :: h2inv = 0.0_dp
   end type reaction_diffusion_problem

   abstract interface
      subroutine irk_rhs_callback(ctx, t, y, f)
         import :: dp
         class(*), intent(in) :: ctx
         real(dp), intent(in) :: t
         real(dp), intent(in) :: y(:)
         real(dp), intent(out) :: f(:)
      end subroutine irk_rhs_callback

      subroutine irk_jacobian_callback(ctx, t, y, a, row, col, nz)
         import :: dp
         class(*), intent(in) :: ctx
         real(dp), intent(in) :: t
         real(dp), intent(in) :: y(:)
         real(dp), intent(out) :: a(:)
         integer, intent(out) :: row(:), col(:)
         integer, intent(out) :: nz
      end subroutine irk_jacobian_callback
   end interface

contains

   elemental function u_exact(x, t) result(u)
      real(dp), intent(in) :: x, t
      real(dp) :: u
      real(dp) :: z

      z = x - wave_speed * t
      u = 1.0_dp / (1.0_dp + exp(wave_speed * z))
   end function u_exact

   elemental function reaction(u) result(f)
      real(dp), intent(in) :: u
      real(dp) :: f

      f = (1.0_dp - u) * u**2
   end function reaction

   elemental function reaction_prime(u) result(df)
      real(dp), intent(in) :: u
      real(dp) :: df

      df = u * (2.0_dp - 3.0_dp * u)
   end function reaction_prime

   subroutine compute_stage_residual(ndof, dt, t_stage, u_old, y_stage, rhs_stage, &
      residual, fnorm, ctx, rhs_cb)
      integer, intent(in) :: ndof
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: t_stage(nstage)
      real(dp), intent(in) :: u_old(ndof)
      real(dp), intent(in) :: y_stage(ndof, nstage)
      real(dp), intent(out) :: rhs_stage(ndof, nstage)
      real(dp), intent(out) :: residual(ndof, nstage)
      real(dp), intent(out) :: fnorm
      class(*), intent(in) :: ctx
      procedure(irk_rhs_callback) :: rhs_cb

      integer :: i, m, j

      do j = 1, nstage
         call rhs_cb(ctx, t_stage(j), y_stage(:, j), rhs_stage(:, j))
      end do

      fnorm = 0.0_dp
      do m = 1, nstage
         do i = 1, ndof
            residual(i, m) = y_stage(i, m) - u_old(i)
            do j = 1, nstage
               residual(i, m) = residual(i, m) - dt * rk_a(m, j) * rhs_stage(i, j)
            end do
            fnorm = max(fnorm, abs(residual(i, m)))
         end do
      end do
   end subroutine compute_stage_residual

   subroutine assemble_coupled_jacobian(ndof, dt, t_stage, y_stage, jac_val, jac_row, &
      jac_col, jac_nz, a, rnr, snr, nz, ctx, jac_cb)
      integer, intent(in) :: ndof
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: t_stage(nstage)
      real(dp), intent(in) :: y_stage(ndof, nstage)
      real(dp), intent(out) :: jac_val(:, :)
      integer, intent(out) :: jac_row(:, :), jac_col(:, :)
      integer, intent(out) :: jac_nz(nstage)
      real(dp), intent(out) :: a(:)
      integer, intent(out) :: rnr(:), snr(:)
      integer, intent(out) :: nz
      class(*), intent(in) :: ctx
      procedure(irk_jacobian_callback) :: jac_cb

      integer :: j, k, m, row_offset, col_offset
      real(dp) :: coeff, entry

      nz = 0
      do j = 1, nstage
         call jac_cb(ctx, t_stage(j), y_stage(:, j), jac_val(:, j), jac_row(:, j), jac_col(:, j), &
            jac_nz(j))
      end do

      do m = 1, nstage
         row_offset = (m - 1) * ndof
         do j = 1, nstage
            col_offset = (j - 1) * ndof
            coeff = -dt * rk_a(m, j)

            do k = 1, jac_nz(j)
               nz = nz + 1
               rnr(nz) = row_offset + jac_row(k, j)
               snr(nz) = col_offset + jac_col(k, j)
               entry = coeff * jac_val(k, j)
               if (m == j .and. jac_row(k, j) == jac_col(k, j)) entry = entry + 1.0_dp
               a(nz) = entry
            end do
         end do
      end do
   end subroutine assemble_coupled_jacobian

   subroutine solve_newton_system_y12m(nunknown, nz, a, snr, nn, rnr, nn1, pivot, rhs, &
      ha, iha, aflag, iflag, ifail, t_y12mb, t_y12mc, t_y12md)
      use y12m, only: y12mb, y12mc, y12md
      integer, intent(in) :: nunknown, nz, nn, nn1, iha
      real(dp), intent(inout) :: a(nn), pivot(nunknown), rhs(nunknown)
      integer, intent(inout) :: snr(nn), rnr(nn1), ha(iha, 11), iflag(10)
      real(dp), intent(inout) :: aflag(8)
      integer, intent(out) :: ifail
      integer, intent(inout) :: t_y12mb, t_y12mc, t_y12md

      integer :: t0, t1

      call system_clock(t0)
      call y12mb(nunknown, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
      call system_clock(t1)
      t_y12mb = t_y12mb + (t1 - t0)
      if (ifail /= 0) return

      call system_clock(t0)
      call y12mc(nunknown, nz, a, snr, nn, rnr, nn1, pivot, rhs, ha, iha, aflag, iflag, ifail)
      call system_clock(t1)
      t_y12mc = t_y12mc + (t1 - t0)
      if (ifail /= 0) return
      iflag(4) = 2

      call system_clock(t0)
      call y12md(nunknown, a, nn, rhs, pivot, snr, ha, iha, iflag, ifail)
      call system_clock(t1)
      t_y12md = t_y12md + (t1 - t0)
   end subroutine solve_newton_system_y12m

   subroutine integrate_radau_iia_fixed(u, t0, dt, nsteps, jac_nz_max, ctx, rhs_cb, jac_cb, &
      newton_iters, final_residuals, t_y12mb, t_y12mc, t_y12md)
      integer, intent(in) :: nsteps, jac_nz_max
      real(dp), intent(in) :: t0, dt
      real(dp), intent(inout) :: u(:)
      class(*), intent(in) :: ctx
      integer, intent(out) :: newton_iters(nsteps)
      real(dp), intent(out) :: final_residuals(nsteps)
      integer, intent(out) :: t_y12mb, t_y12mc, t_y12md
      procedure(irk_rhs_callback) :: rhs_cb
      procedure(irk_jacobian_callback) :: jac_cb

      integer :: ndof, nunknown, nz_max, nz, nn, nn1, iha
      integer :: iflag(10), ifail
      integer :: step, iter, iter_used
      real(dp) :: aflag(8)
      real(dp) :: t_stage(nstage), fnorm
      real(dp), allocatable :: y_stage(:, :), rhs_stage(:, :), residual(:, :)
      real(dp), allocatable :: jac_val(:, :), a(:), pivot(:), delta(:)
      integer, allocatable :: jac_row(:, :), jac_col(:, :), jac_nz(:)
      integer, allocatable :: snr(:), rnr(:), ha(:, :)
      logical :: converged

      ndof = size(u)
      nunknown = nstage * ndof
      nz_max = nstage * nstage * jac_nz_max
      nn = max(y12m_workspace_factor * nunknown, y12m_workspace_factor * nz_max)
      nn1 = nn
      iha = nunknown

      allocate(y_stage(ndof, nstage), rhs_stage(ndof, nstage), residual(ndof, nstage))
      allocate(jac_val(jac_nz_max, nstage), jac_row(jac_nz_max, nstage), jac_col(jac_nz_max, nstage))
      allocate(jac_nz(nstage), a(nn), pivot(nunknown), delta(nunknown))
      allocate(snr(nn), rnr(nn1), ha(nunknown, 11))

      iflag = 0
      iflag(2) = 3        ! Markowitz search width for pivot selection
      iflag(3) = 1        ! enable column interchanges
      iflag(4) = 1        ! compute ordering on the first Jacobian
      iflag(5) = 1        ! each factorization solves one Newton RHS; no saved L solve phase
      aflag(1) = y12m_growth_limit
      aflag(2) = 1.0e-12_dp
      aflag(3) = 1.0e+16_dp
      aflag(4) = 1.0e-12_dp
      aflag(5:8) = 0.0_dp

      t_y12mb = 0
      t_y12mc = 0
      t_y12md = 0

      do step = 1, nsteps
         t_stage = t0 + (real(step - 1, kind=dp) + rk_c) * dt
         y_stage(:, 1) = u
         y_stage(:, 2) = u
         converged = .false.
         iter_used = 0

         do iter = 1, max_newton
            call compute_stage_residual(ndof, dt, t_stage, u, y_stage, rhs_stage, residual, fnorm, &
               ctx, rhs_cb)
            if (fnorm < newton_tol) then
               converged = .true.
               iter_used = iter - 1
               exit
            end if

            delta = -reshape(residual, [nunknown])
            call assemble_coupled_jacobian(ndof, dt, t_stage, y_stage, jac_val, jac_row, jac_col, &
               jac_nz, a, rnr, snr, nz, ctx, jac_cb)
            call solve_newton_system_y12m(nunknown, nz, a, snr, nn, rnr, nn1, pivot, delta, ha, iha, &
               aflag, iflag, ifail, t_y12mb, t_y12mc, t_y12md)

            if (ifail /= 0) then
               write(error_unit, '(a,i0,a,i0)') &
                  'ERROR: Y12M solve at time step ', step, ' returned IFAIL = ', ifail
               stop 1
            end if

            y_stage = y_stage + reshape(delta, [ndof, nstage])
         end do

         if (.not. converged) then
            call compute_stage_residual(ndof, dt, t_stage, u, y_stage, rhs_stage, residual, fnorm, &
               ctx, rhs_cb)
            if (fnorm < newton_tol) then
               converged = .true.
               iter_used = max_newton
            end if
         end if

         if (.not. converged) then
            write(error_unit, '(a,i0,a,es12.4)') &
               'ERROR: Newton failed at time step ', step, ', residual = ', fnorm
            stop 1
         end if

         newton_iters(step) = iter_used
         final_residuals(step) = fnorm
         u = u + dt * matmul(rhs_stage, rk_b)
      end do
   end subroutine integrate_radau_iia_fixed

   subroutine reaction_diffusion_rhs(ctx, t, y, f)
      class(*), intent(in) :: ctx
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: f(:)

      integer :: k
      real(dp) :: left_bc, right_bc

      select type (problem => ctx)
      type is (reaction_diffusion_problem)
         left_bc = u_exact(x_left, t)
         right_bc = u_exact(x_right, t)

         do k = 1, problem%ndof
            f(k) = (-2.0_dp * y(k)) * problem%h2inv + reaction(y(k))
            if (k > 1) f(k) = f(k) + y(k - 1) * problem%h2inv
            if (k < problem%ndof) f(k) = f(k) + y(k + 1) * problem%h2inv
         end do

         f(1) = f(1) + left_bc * problem%h2inv
         f(problem%ndof) = f(problem%ndof) + right_bc * problem%h2inv
      class default
         write(error_unit, '(a)') 'ERROR: unsupported RHS callback context'
         stop 1
      end select
   end subroutine reaction_diffusion_rhs

   subroutine reaction_diffusion_jacobian(ctx, t, y, a, row, col, nz)
      class(*), intent(in) :: ctx
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: a(:)
      integer, intent(out) :: row(:), col(:)
      integer, intent(out) :: nz

      integer :: k

      select type (problem => ctx)
      type is (reaction_diffusion_problem)
         ! This example has a time-independent semidiscrete Jacobian, but keep a
         ! reference to `t` here so the callback still matches the general
         ! J(t,y) API used by the fixed-step Radau driver.
         associate (unused_t => t)
         end associate
         nz = 0
         do k = 1, problem%ndof
            if (k > 1) then
               nz = nz + 1
               row(nz) = k
               col(nz) = k - 1
               a(nz) = problem%h2inv
            end if

            nz = nz + 1
            row(nz) = k
            col(nz) = k
            a(nz) = -2.0_dp * problem%h2inv + reaction_prime(y(k))

            if (k < problem%ndof) then
               nz = nz + 1
               row(nz) = k
               col(nz) = k + 1
               a(nz) = problem%h2inv
            end if
         end do
      class default
         write(error_unit, '(a)') 'ERROR: unsupported Jacobian callback context'
         stop 1
      end select
   end subroutine reaction_diffusion_jacobian

   subroutine write_solution(N, h, t, u, outfile)
      integer, intent(in) :: N
      real(dp), intent(in) :: h, t
      real(dp), intent(in) :: u(:)
      character(len=*), intent(in) :: outfile

      integer :: i, funit
      real(dp) :: x

      open(newunit=funit, file=outfile, status='unknown', action='write')
      write(funit, '(a)') '# 1D reaction-diffusion equation -- 2-stage Radau IIA'
      write(funit, '(a,es10.3,a,i0)') '# T=', t, '  N=', N
      write(funit, '(a)') '# Columns: x  u_numerical  u_exact'

      x = x_left
      write(funit, '(3(1x,es14.6))') x, u_exact(x, t), u_exact(x, t)

      do i = 1, size(u)
         x = x_left + real(i, kind=dp) * h
         write(funit, '(3(1x,es14.6))') x, u(i), u_exact(x, t)
      end do

      x = x_right
      write(funit, '(3(1x,es14.6))') x, u_exact(x, t), u_exact(x, t)
      close(funit)
   end subroutine write_solution

   subroutine run(N, nsteps, T, outfile)
      integer, intent(in) :: N, nsteps
      real(dp), intent(in) :: T
      character(len=*), intent(in) :: outfile

      integer :: ndof, jac_nz_max
      real(dp) :: h, dt, h2inv
      real(dp), allocatable :: u(:), final_residuals(:)
      integer, allocatable :: newton_iters(:)
      integer :: i, step, clock_rate
      integer :: t_y12mb, t_y12mc, t_y12md
      real(dp) :: s_y12mb, s_y12mc, s_y12md
      real(dp) :: err_max, err_rms, diff
      type(reaction_diffusion_problem) :: problem

      ndof = N - 2
      h = (x_right - x_left) / real(N - 1, kind=dp)
      dt = T / real(nsteps, kind=dp)
      h2inv = 1.0_dp / h**2

      write(output_unit, '(a,i0,a,i0,a)') 'Grid points   : ', N, '  (', ndof, ' interior DOFs)'
      write(output_unit, '(a,es10.3)') 'Domain length : ', x_right - x_left
      write(output_unit, '(a,i0)') 'Time steps    : ', nsteps
      write(output_unit, '(a,es10.3)') 'Final time T  : ', T
      write(output_unit, '(a,es10.3)') 'dt            : ', dt
      write(output_unit, '(a,es10.3)') 'Newton tol    : ', newton_tol
      write(output_unit, '(a,i0)') 'Max Newton it : ', max_newton

      allocate(u(ndof), newton_iters(nsteps), final_residuals(nsteps))

      do i = 1, ndof
         u(i) = u_exact(x_left + real(i, kind=dp) * h, 0.0_dp)
      end do

      problem%ndof = ndof
      problem%h2inv = h2inv
      jac_nz_max = 3 * ndof - 2   ! tridiagonal semidiscrete Jacobian

      call system_clock(count_rate=clock_rate)
      call integrate_radau_iia_fixed(u, 0.0_dp, dt, nsteps, jac_nz_max, problem, &
         reaction_diffusion_rhs, reaction_diffusion_jacobian, newton_iters, final_residuals, &
         t_y12mb, t_y12mc, t_y12md)

      write(output_unit, '(/,a)') 'Per-step Newton iteration counts:'
      write(output_unit, '(a)') '  step  iterations  final_stage_residual'
      do step = 1, nsteps
         write(output_unit, '(2x,i4,2x,i4,4x,es12.4)') step, newton_iters(step), final_residuals(step)
      end do

      err_max = 0.0_dp
      err_rms = 0.0_dp
      do i = 1, ndof
         diff = u(i) - u_exact(x_left + real(i, kind=dp) * h, T)
         err_max = max(err_max, abs(diff))
         err_rms = err_rms + diff**2
      end do
      err_rms = sqrt(err_rms / real(ndof, kind=dp))

      write(output_unit, '(/,a,es12.4)') 'Final max error : ', err_max
      write(output_unit, '(a,es12.4)') 'Final RMS error : ', err_rms
      write(output_unit, '(a,f8.3)') 'Average Newton iterations per step : ', &
         real(sum(newton_iters), kind=dp) / real(nsteps, kind=dp)
      write(output_unit, '(a,i0)') 'Maximum Newton iterations          : ', maxval(newton_iters)

      call write_solution(N, h, T, u, trim(outfile))
      write(output_unit, '(a)') 'Solution written to ' // trim(outfile)

      s_y12mb = real(t_y12mb, kind=dp) / real(clock_rate, kind=dp)
      s_y12mc = real(t_y12mc, kind=dp) / real(clock_rate, kind=dp)
      s_y12md = real(t_y12md, kind=dp) / real(clock_rate, kind=dp)

      write(output_unit, '(/,a)') 'Timing summary:'
      write(output_unit, '(a,g14.6,a)') '  y12mb total           : ', s_y12mb, ' s'
      write(output_unit, '(a,g14.6,a)') '  y12mc total           : ', s_y12mc, ' s'
      write(output_unit, '(a,g14.6,a)') '  y12md total           : ', s_y12md, ' s'
   end subroutine run
end module reaction_diffusion_radau_solver

program reaction_diffusion_radau
   use reaction_diffusion_radau_solver, only: run, dp
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none

   integer :: N, nsteps, ios, iarg, nargs, pos_count
   real(dp) :: T
   character(len=256) :: arg, outfile

   N = 41
   nsteps = 20
   T = 1.0_dp
   outfile = 'reaction_diffusion_radau.dat'

   nargs = command_argument_count()
   iarg = 1
   pos_count = 0

   do while (iarg <= nargs)
      call get_command_argument(iarg, arg)

      if (trim(arg) == '--help' .or. trim(arg) == '-h') then
         write(*, '(a)') 'Usage: reaction_diffusion_radau [--help] [N] [nsteps] [T] [output_file]'
         write(*, '(a)') '  N           total grid size (>= 3, default 41)'
         write(*, '(a)') '  nsteps      number of time steps (default 20)'
         write(*, '(a)') '  T           final time (default 1.0)'
         write(*, '(a)') '  output_file output file (default: reaction_diffusion_radau.dat)'
         stop 0
      else
         pos_count = pos_count + 1
         select case (pos_count)
         case (1)
            read(arg, *, iostat=ios) N
            if (ios /= 0 .or. N < 3) then
               write(error_unit, '(a)') 'ERROR: N must be an integer >= 3'
               stop 1
            end if
         case (2)
            read(arg, *, iostat=ios) nsteps
            if (ios /= 0 .or. nsteps < 1) then
               write(error_unit, '(a)') 'ERROR: nsteps must be a positive integer'
               stop 1
            end if
         case (3)
            read(arg, *, iostat=ios) T
            if (ios /= 0 .or. T <= 0.0_dp) then
               write(error_unit, '(a)') 'ERROR: T must be a positive real'
               stop 1
            end if
         case (4)
            outfile = trim(arg)
         end select
      end if

      iarg = iarg + 1
   end do

   call run(N, nsteps, T, trim(outfile))
end program reaction_diffusion_radau
