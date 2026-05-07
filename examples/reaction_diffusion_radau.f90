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
   use y12m_example_helpers, only: write_three_column_output
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

   type, abstract :: linear_solver
   contains
      procedure(linear_solve_if), deferred :: solve
   end type linear_solver

   abstract interface
      subroutine linear_solve_if(self, n, nz, a, rnr, snr, rhs, istat)
         import :: linear_solver, dp
         class(linear_solver), intent(inout) :: self
         integer, intent(in) :: n, nz
         real(dp), intent(inout) :: a(:), rhs(:)
         integer, intent(inout) :: rnr(:), snr(:)
         integer, intent(inout) :: istat
      end subroutine linear_solve_if

      subroutine irk_rhs_callback(n, t, y, f, ipar, rpar, istat)
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: t
         real(dp), intent(in) :: y(n)
         real(dp), intent(out) :: f(n)
         integer, intent(in) :: ipar(*)
         real(dp), intent(in) :: rpar(*)
         integer, intent(inout) :: istat
      end subroutine irk_rhs_callback

      subroutine irk_jacobian_callback(n, nrow, ncol, t, y, a, lda, row, col, ldc, nz, ipar, rpar, istat)
         import :: dp
         integer, intent(in) :: n, nrow, ncol, lda, ldc
         real(dp), intent(in) :: t
         real(dp), intent(in) :: y(n)
         real(dp), intent(out) :: a(lda)
         integer, intent(out) :: row(ldc), col(ldc)
         integer, intent(out) :: nz
         integer, intent(in) :: ipar(*)
         real(dp), intent(in) :: rpar(*)
         integer, intent(inout) :: istat
      end subroutine irk_jacobian_callback
   end interface

   type, extends(linear_solver) :: y12m_linear_solver
      integer :: n = 0
      integer :: nn = 0
      integer :: nn1 = 0
      integer :: iha = 0
      integer :: iflag(10) = 0
      real(dp) :: aflag(8) = 0.0_dp
      real(dp), allocatable :: pivot(:)
      real(dp), allocatable :: work_a(:)
      integer, allocatable :: work_snr(:), work_rnr(:)
      integer, allocatable :: ha(:, :)
   contains
      procedure :: initialize => y12m_initialize
      procedure :: solve => y12m_solve
   end type y12m_linear_solver

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

   subroutine y12m_initialize(self, n, nz_max)
      class(y12m_linear_solver), intent(inout) :: self
      integer, intent(in) :: n, nz_max

      self%n = n
      self%nn = max(y12m_workspace_factor * n, y12m_workspace_factor * nz_max)
      self%nn1 = self%nn
      self%iha = n

      allocate(self%pivot(n), self%work_a(self%nn), self%work_snr(self%nn), self%work_rnr(self%nn1), &
         self%ha(n, 11))

      self%iflag = 0
      self%iflag(2) = 3
      self%iflag(3) = 1
      self%iflag(4) = 1
      self%iflag(5) = 1
      self%aflag(1) = y12m_growth_limit
      self%aflag(2) = 1.0e-12_dp
      self%aflag(3) = 1.0e+16_dp
      self%aflag(4) = 1.0e-12_dp
      self%aflag(5:8) = 0.0_dp
   end subroutine y12m_initialize

   subroutine y12m_solve(self, n, nz, a, rnr, snr, rhs, istat)
      use y12m, only: y12mb, y12mc, y12md
      class(y12m_linear_solver), intent(inout) :: self
      integer, intent(in) :: n, nz
      real(dp), intent(inout) :: a(:), rhs(:)
      integer, intent(inout) :: rnr(:), snr(:)
      integer, intent(inout) :: istat

      integer :: ifail

      self%work_a(1:nz) = a(1:nz)
      self%work_snr(1:nz) = snr(1:nz)
      self%work_rnr(1:nz) = rnr(1:nz)

      ifail = 0
      call y12mb(n, nz, self%work_a, self%work_snr, self%nn, self%work_rnr, self%nn1, self%ha, self%iha, &
         self%aflag, self%iflag, ifail)
      if (ifail /= 0) then
         istat = ifail
         return
      end if

      call y12mc(n, nz, self%work_a, self%work_snr, self%nn, self%work_rnr, self%nn1, self%pivot, rhs, self%ha, &
         self%iha, &
         self%aflag, self%iflag, ifail)
      if (ifail /= 0) then
         istat = ifail
         return
      end if

      self%iflag(4) = 2
      call y12md(n, self%work_a, self%nn, rhs, self%pivot, self%work_snr, self%ha, self%iha, self%iflag, ifail)
      if (ifail /= 0) then
         istat = ifail
      end if
   end subroutine y12m_solve

   ! Assemble stage residuals for the coupled 2-stage Radau system:
   !
   !   G1 = Y1 - u_n - dt * (a11*F1 + a12*F2)
   !   G2 = Y2 - u_n - dt * (a21*F1 + a22*F2)
   !
   ! where Fj = F(t_n + c_j*dt, Yj).
   subroutine compute_stage_residual(ndof, dt, t_stage, u_old, y_stage, rhs_stage, residual, fnorm, &
      ipar, rpar, rhs_cb, istat)
      integer, intent(in) :: ndof
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: t_stage(nstage)
      real(dp), intent(in) :: u_old(ndof)
      real(dp), intent(in) :: y_stage(ndof, nstage)
      real(dp), intent(out) :: rhs_stage(ndof, nstage)
      real(dp), intent(out) :: residual(ndof, nstage)
      real(dp), intent(out) :: fnorm
      integer, intent(in) :: ipar(*)
      real(dp), intent(in) :: rpar(*)
      procedure(irk_rhs_callback) :: rhs_cb
      integer, intent(inout) :: istat

      call rhs_cb(ndof, t_stage(1), y_stage(:, 1), rhs_stage(:, 1), ipar, rpar, istat)
      if (istat /= 0) return
      call rhs_cb(ndof, t_stage(2), y_stage(:, 2), rhs_stage(:, 2), ipar, rpar, istat)
      if (istat /= 0) return

      residual(:, 1) = y_stage(:, 1) - u_old - dt * (rk_a(1, 1) * rhs_stage(:, 1) + rk_a(1, 2) * rhs_stage(:, 2))
      residual(:, 2) = y_stage(:, 2) - u_old - dt * (rk_a(2, 1) * rhs_stage(:, 1) + rk_a(2, 2) * rhs_stage(:, 2))

      fnorm = max(maxval(abs(residual(:, 1))), maxval(abs(residual(:, 2))))
   end subroutine compute_stage_residual

   ! Assemble the Newton matrix for the coupled stage residuals:
   !
   !   dG_m/dY_j = I*delta_mj - dt * a_mj * J_j,
   !
   ! where J_j = dF/dy(t_n + c_j*dt, Y_j).
   !
   ! This creates a 2x2 block sparse COO matrix with ndof x ndof blocks.
   subroutine assemble_coupled_jacobian(ndof, dt, t_stage, y_stage, jac_nz_max, jac_val, jac_row, &
      jac_col, jac_nz, a, rnr, snr, nz, ipar, rpar, jac_cb, istat)
      integer, intent(in) :: ndof, jac_nz_max
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: t_stage(nstage)
      real(dp), intent(in) :: y_stage(ndof, nstage)
      real(dp), intent(out) :: jac_val(jac_nz_max, nstage)
      integer, intent(out) :: jac_row(jac_nz_max, nstage), jac_col(jac_nz_max, nstage)
      integer, intent(out) :: jac_nz(nstage)
      real(dp), intent(out) :: a(:)
      integer, intent(out) :: rnr(:), snr(:)
      integer, intent(out) :: nz
      integer, intent(in) :: ipar(*)
      real(dp), intent(in) :: rpar(*)
      procedure(irk_jacobian_callback) :: jac_cb
      integer, intent(inout) :: istat

      integer :: j, k, m, row_offset, col_offset
      real(dp) :: coeff, entry

      nz = 0

      call jac_cb(ndof, ndof, ndof, t_stage(1), y_stage(:, 1), jac_val(:, 1), jac_nz_max, jac_row(:, 1), &
         jac_col(:, 1), jac_nz_max, jac_nz(1), ipar, rpar, istat)
      if (istat /= 0) return

      call jac_cb(ndof, ndof, ndof, t_stage(2), y_stage(:, 2), jac_val(:, 2), jac_nz_max, jac_row(:, 2), &
         jac_col(:, 2), jac_nz_max, jac_nz(2), ipar, rpar, istat)
      if (istat /= 0) return

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

   ! Perform one fully implicit 2-stage Radau IIA step.
   !
   ! Newton is full Newton (not modified Newton): the Jacobian is rebuilt and
   ! refactorized every Newton iteration.
   subroutine integrate_radau_iia_step(u, t, dt, jac_nz_max, ipar, rpar, rhs_cb, jac_cb, linsolve, &
      newton_iter, final_residual)
      real(dp), intent(inout) :: u(:)
      real(dp), intent(in) :: t, dt
      integer, intent(in) :: jac_nz_max
      integer, intent(in) :: ipar(*)
      real(dp), intent(in) :: rpar(*)
      procedure(irk_rhs_callback) :: rhs_cb
      procedure(irk_jacobian_callback) :: jac_cb
      class(linear_solver), intent(inout) :: linsolve
      integer, intent(out) :: newton_iter
      real(dp), intent(out) :: final_residual

      integer :: ndof, nunknown, nz_max, nz
      integer :: iter, istat
      real(dp) :: t_stage(nstage), fnorm
      real(dp), allocatable :: y_stage(:, :), rhs_stage(:, :), residual(:, :)
      real(dp), allocatable :: jac_val(:, :), a(:), delta(:)
      integer, allocatable :: jac_row(:, :), jac_col(:, :), jac_nz(:)
      integer, allocatable :: snr(:), rnr(:)
      logical :: converged

      ndof = size(u)
      nunknown = nstage * ndof
      nz_max = nstage * nstage * jac_nz_max

      allocate(y_stage(ndof, nstage), rhs_stage(ndof, nstage), residual(ndof, nstage))
      allocate(jac_val(jac_nz_max, nstage), jac_row(jac_nz_max, nstage), jac_col(jac_nz_max, nstage))
      allocate(jac_nz(nstage), a(nz_max), delta(nunknown))
      allocate(snr(nz_max), rnr(nz_max))

      t_stage = t + rk_c * dt
      y_stage(:, 1) = u
      y_stage(:, 2) = u
      converged = .false.
      newton_iter = 0
      final_residual = huge(1.0_dp)

      do iter = 1, max_newton
         istat = 0
         call compute_stage_residual(ndof, dt, t_stage, u, y_stage, rhs_stage, residual, fnorm, &
            ipar, rpar, rhs_cb, istat)
         if (istat /= 0) then
            write(error_unit, '(a,i0)') 'ERROR: RHS callback returned ISTAT = ', istat
            stop 1
         end if

         if (fnorm < newton_tol) then
            converged = .true.
            newton_iter = iter - 1
            final_residual = fnorm
            exit
         end if

         delta = -reshape(residual, [nunknown])

         istat = 0
         call assemble_coupled_jacobian(ndof, dt, t_stage, y_stage, jac_nz_max, jac_val, jac_row, jac_col, &
            jac_nz, a, rnr, snr, nz, ipar, rpar, jac_cb, istat)
         if (istat /= 0) then
            write(error_unit, '(a,i0)') 'ERROR: Jacobian callback returned ISTAT = ', istat
            stop 1
         end if

         istat = 0
         call linsolve%solve(nunknown, nz, a, rnr, snr, delta, istat)
         if (istat /= 0) then
            write(error_unit, '(a,i0)') 'ERROR: linear solve returned ISTAT = ', istat
            stop 1
         end if

         y_stage = y_stage + reshape(delta, [ndof, nstage])
      end do

      if (.not. converged) then
         istat = 0
         call compute_stage_residual(ndof, dt, t_stage, u, y_stage, rhs_stage, residual, fnorm, &
            ipar, rpar, rhs_cb, istat)
         if (istat /= 0) then
            write(error_unit, '(a,i0)') 'ERROR: RHS callback returned ISTAT = ', istat
            stop 1
         end if
         if (fnorm < newton_tol) then
            converged = .true.
            newton_iter = max_newton
            final_residual = fnorm
         end if
      end if

      if (.not. converged) then
         write(error_unit, '(a,es12.4)') 'ERROR: Newton failed, residual = ', fnorm
         stop 1
      end if

      u = u + dt * (rk_b(1) * rhs_stage(:, 1) + rk_b(2) * rhs_stage(:, 2))
   end subroutine integrate_radau_iia_step

   subroutine reaction_diffusion_rhs(n, t, y, f, ipar, rpar, istat)
      integer, intent(in) :: n
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(n)
      real(dp), intent(out) :: f(n)
      integer, intent(in) :: ipar(*)
      real(dp), intent(in) :: rpar(*)
      integer, intent(inout) :: istat

      integer :: i
      real(dp) :: left_bc, right_bc, h2inv

      if (ipar(1) /= n) then
         istat = 1
         return
      end if

      h2inv = rpar(1)
      left_bc = u_exact(x_left, t)
      right_bc = u_exact(x_right, t)

      f = (-2.0_dp * y) * h2inv + reaction(y)
      if (n >= 2) then
         f(1) = f(1) + y(2) * h2inv
         if (n > 2) then
            do concurrent (i = 2:n - 1)
               f(i) = f(i) + (y(i - 1) + y(i + 1)) * h2inv
            end do
         end if
         f(n) = f(n) + y(n - 1) * h2inv
      end if
      f(1) = f(1) + left_bc * h2inv
      f(n) = f(n) + right_bc * h2inv
   end subroutine reaction_diffusion_rhs

   subroutine reaction_diffusion_jacobian(n, nrow, ncol, t, y, a, lda, row, col, ldc, nz, ipar, rpar, istat)
      integer, intent(in) :: n, nrow, ncol, lda, ldc
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(n)
      real(dp), intent(out) :: a(lda)
      integer, intent(out) :: row(ldc), col(ldc)
      integer, intent(out) :: nz
      integer, intent(in) :: ipar(*)
      real(dp), intent(in) :: rpar(*)
      integer, intent(inout) :: istat

      integer :: k
      real(dp) :: h2inv

      if (ipar(1) /= n) then
         istat = 1
         return
      end if

      if (nrow /= n .or. ncol /= n) then
         istat = 2
         return
      end if

      ! The semidiscrete Jacobian is time-independent for this travelling wave
      ! setup. Keep t in the signature to preserve the generic J(t,y) callback API.
      h2inv = rpar(1)
      nz = 0
      do k = 1, n
         if (k > 1) then
            nz = nz + 1
            if (nz > min(lda, ldc)) then
               istat = 3
               return
            end if
            row(nz) = k
            col(nz) = k - 1
            a(nz) = h2inv
         end if

         nz = nz + 1
         if (nz > min(lda, ldc)) then
            istat = 3
            return
         end if
         row(nz) = k
         col(nz) = k
         a(nz) = -2.0_dp * h2inv + reaction_prime(y(k))

         if (k < n) then
            nz = nz + 1
            if (nz > min(lda, ldc)) then
               istat = 3
               return
            end if
            row(nz) = k
            col(nz) = k + 1
            a(nz) = h2inv
         end if
      end do
   end subroutine reaction_diffusion_jacobian

   subroutine write_solution(N, h, t, u, outfile)
      integer, intent(in) :: N
      real(dp), intent(in) :: h, t
      real(dp), intent(in) :: u(:)
      character(len=*), intent(in) :: outfile

      integer :: i
      real(dp), allocatable :: x(:), u_num(:), u_ref(:)
      character(len=80), dimension(2) :: header

      allocate(x(N), u_num(N), u_ref(N))

      do i = 1, N
         x(i) = x_left + real(i - 1, dp) * h
         u_ref(i) = u_exact(x(i), t)
      end do

      u_num(1) = u_ref(1)
      u_num(N) = u_ref(N)
      if (N > 2) u_num(2:N - 1) = u

      write(header(1), '(a)') '1D reaction-diffusion equation -- 2-stage Radau IIA'
      write(header(2), '(a,es10.3,a,i0)') 'T=', t, '  N=', N

      call write_three_column_output(x, u_num, u_ref, trim(outfile), header, size(header), &
         'Columns: x  u_numerical  u_exact')
   end subroutine write_solution

   subroutine run(N, nsteps, T, outfile)
      integer, intent(in) :: N, nsteps
      real(dp), intent(in) :: T
      character(len=*), intent(in) :: outfile

      integer :: ndof, jac_nz_max
      real(dp) :: h, dt, t_now, h2inv
      real(dp), allocatable :: u(:), final_residuals(:)
      integer, allocatable :: newton_iters(:)
      integer :: i, step
      real(dp) :: err_max, err_rms, diff
      integer :: ipar(1)
      real(dp) :: rpar(1)
      type(y12m_linear_solver) :: y12m_solver

      ndof = N - 2
      h = (x_right - x_left) / real(N - 1, dp)
      dt = T / real(nsteps, dp)
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
         u(i) = u_exact(x_left + real(i, dp) * h, 0.0_dp)
      end do

      ipar(1) = ndof
      rpar(1) = h2inv
      jac_nz_max = 3 * ndof - 2
      call y12m_solver%initialize(2 * ndof, 4 * jac_nz_max)

      t_now = 0.0_dp
      do step = 1, nsteps
         call integrate_radau_iia_step(u, t_now, dt, jac_nz_max, ipar, rpar, reaction_diffusion_rhs, &
            reaction_diffusion_jacobian, y12m_solver, newton_iters(step), final_residuals(step))
         t_now = t_now + dt
      end do

      write(output_unit, '(/,a)') 'Per-step Newton iteration counts:'
      write(output_unit, '(a)') '  step  iterations  final_stage_residual'
      do step = 1, nsteps
         write(output_unit, '(2x,i4,2x,i4,4x,es12.4)') step, newton_iters(step), final_residuals(step)
      end do

      err_max = 0.0_dp
      err_rms = 0.0_dp
      do i = 1, ndof
         diff = u(i) - u_exact(x_left + real(i, dp) * h, T)
         err_max = max(err_max, abs(diff))
         err_rms = err_rms + diff**2
      end do
      err_rms = sqrt(err_rms / real(ndof, dp))

      write(output_unit, '(/,a,es12.4)') 'Final max error : ', err_max
      write(output_unit, '(a,es12.4)') 'Final RMS error : ', err_rms
      write(output_unit, '(a,f8.3)') 'Average Newton iterations per step : ', &
         real(sum(newton_iters), dp) / real(nsteps, dp)
      write(output_unit, '(a,i0)') 'Maximum Newton iterations          : ', maxval(newton_iters)

      call write_solution(N, h, T, u, trim(outfile))
      write(output_unit, '(a)') 'Solution written to ' // trim(outfile)
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
