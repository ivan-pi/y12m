! SPDX-License-Identifier: GPL-2.0-only
!
! heat_implicit.f90 - Implicit time-stepping for the 2-D heat equation
!
! Solves the 2-D heat equation
!
!   u_t = kappa * (u_xx + u_yy)    on  (0,1) x (0,1),  t in (0, T]
!
! with homogeneous Dirichlet boundary conditions
!
!   u = 0   on the boundary
!
! and smooth initial condition
!
!   u(x, y, 0) = sin(pi x) sin(pi y)
!
! Exact solution (separation of variables):
!
!   u(x, y, t) = exp(-2 pi^2 kappa t) sin(pi x) sin(pi y)
!
! Discretisation
! --------------
! Total grid N x N (including boundary nodes), spacing h = 1/(N-1).
! Interior nodes: ndof = (N-2)^2.
! DOF index: k(i,j) = (j-2)*(N-2) + (i-2) + 1.
!
! Time integration: backward (implicit) Euler with step dt = T / nsteps.
!
!   A u^{n+1} = u^n,    A = I + r K,    r = dt kappa / h^2
!
! K is the ndof x ndof positive-definite 5-point stiffness matrix
! (diagonal = 4, off-diagonal = -1 for each interior neighbour).
!
! Solver strategy -- Case (ii): same matrix, many right-hand sides
! ----------------------------------------------------------------
! The coefficient matrix A is CONSTANT across all nsteps time steps.
!
!   * y12mb + y12mc are called ONCE before the time loop to factorize A.
!   * y12md (IFLAG(5) = 2 for step 1, then IFLAG(5) = 3 for all later
!     steps) is called at every time step to perform forward and backward
!     substitution for the new right-hand side.
!
! With an iterative solver the same problem would require re-converging
! the linear solver at every time step; a direct solver amortises the
! factorization cost over all nsteps solves.
!
! Usage
! -----
!   heat_implicit [--help] [N] [nsteps] [T] [kappa] [output_file]
!
!     N           total grid size (>= 3, default 42)
!     nsteps      number of time steps (default 200)
!     T           final time (default 0.1)
!     kappa       thermal diffusivity > 0 (default 1.0)
!     output_file output data file name (default: heat_implicit.dat)
!

! ==============================================================
! Module heat_implicit_solver
! ==============================================================
module heat_implicit_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  !> Exact solution at time t.
  elemental function u_exact(x, y, t, kappa) result(u)
    real(dp), intent(in) :: x, y, t, kappa
    real(dp) :: u
    u = exp(-2.0_dp * pi**2 * kappa * t) * sin(pi * x) * sin(pi * y)
  end function u_exact

  !> Assemble the matrix A = I + r*K in COO format.
  !>
  !> K is the 5-point positive-definite Laplacian stiffness
  !> (diagonal = 4, off-diag = -1 for each interior neighbour).
  !> A has diagonal = 1 + 4*r and off-diagonal = -r.
  !> Entries are emitted in a fixed row-major order so that the
  !> same (rnr, snr) pattern can be reused across time steps.
  subroutine assemble_matrix(N, r, a, rnr, snr, nz)
    integer, intent(in) :: N
    real(dp), intent(in) :: r
    real(dp), intent(out) :: a(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    integer :: ndof_1d, i, j, k_row

    ndof_1d = N - 2
    nz = 0

    do j = 2, N - 1
      do i = 2, N - 1
        k_row = (j - 2) * ndof_1d + (i - 2) + 1

        ! Diagonal: 1 + 4r
        nz = nz + 1
        rnr(nz) = k_row
        snr(nz) = k_row
        a(nz) = 1.0_dp + 4.0_dp * r

        ! West neighbour (i-1, j)
        if (i > 2) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 2) * ndof_1d + (i - 3) + 1
          a(nz) = -r
        end if

        ! East neighbour (i+1, j)
        if (i < N - 1) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 2) * ndof_1d + (i - 1) + 1
          a(nz) = -r
        end if

        ! South neighbour (i, j-1)
        if (j > 2) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 3) * ndof_1d + (i - 2) + 1
          a(nz) = -r
        end if

        ! North neighbour (i, j+1)
        if (j < N - 1) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 1) * ndof_1d + (i - 2) + 1
          a(nz) = -r
        end if

      end do
    end do
  end subroutine assemble_matrix

  ! ---------------------------------------------------------------
  ! Main driver: assemble, factorize once, time-step with LU reuse
  ! ---------------------------------------------------------------
  subroutine run(N, nsteps, T, kappa, outfile)
    use y12m, only: y12mb, y12mc, y12md
    integer, intent(in) :: N, nsteps
    real(dp), intent(in) :: T, kappa
    character(len=*), intent(in) :: outfile

    integer :: ndof_1d, ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h, dt, r
    real(dp), allocatable :: a(:), pivot(:), b(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail
    integer :: i, j, k, step, funit
    integer :: t0, t1, clock_rate
    integer :: t_assemble, t_factor, t_steps
    real(dp) :: s_assemble, s_factor, s_steps
    real(dp) :: err_max

    ndof_1d = N - 2
    ndof = ndof_1d * ndof_1d
    h = 1.0_dp / real(N - 1, dp)
    dt = T / real(nsteps, dp)
    r = dt * kappa / h**2

    write(output_unit, '(a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
    write(output_unit, '(a,i0)')     'Time steps    : ', nsteps
    write(output_unit, '(a,es10.3)') 'Final time T  : ', T
    write(output_unit, '(a,es10.3)') 'dt            : ', dt
    write(output_unit, '(a,es10.3)') 'kappa         : ', kappa
    write(output_unit, '(a,es10.3)') 'r = dt*k/h^2  : ', r

    ! Workspace: 5-point stencil -> at most 5*ndof non-zeros
    nz_max = 5 * ndof
    nn = max(3 * nz_max, 2 * (ndof_1d + 2) * ndof)
    nn1 = nn
    iha = ndof
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))

    call system_clock(count_rate=clock_rate)

    ! =====================================================================
    ! 1. Assemble A and set initial condition b = u^0
    ! =====================================================================
    call system_clock(t0)
    call assemble_matrix(N, r, a, rnr, snr, nz)

    do j = 2, N - 1
      do i = 2, N - 1
        k = (j - 2) * ndof_1d + (i - 2) + 1
        b(k) = u_exact(real(i - 1, dp) * h, real(j - 1, dp) * h, 0.0_dp, kappa)
      end do
    end do
    call system_clock(t1)
    t_assemble = t1 - t0

    ! =====================================================================
    ! 2. Factorize A ONCE with IFLAG(5) = 2 (keep L for subsequent solves)
    ! =====================================================================
    call system_clock(t0)

    iflag = 0
    iflag(2) = 3        ! Markowitz search width
    iflag(3) = 1        ! use column interchanges
    iflag(4) = 0        ! no structural reuse
    iflag(5) = 2        ! keep L after factorization (Case ii)
    aflag(1) = 16.0_dp
    aflag(2) = 1.0e-12_dp
    aflag(3) = 1.0e+16_dp
    aflag(4) = 1.0e-12_dp
    aflag(5:8) = 0.0_dp

    call y12mb(ndof, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12mb returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if

    ! y12mc factorizes A and applies the row permutation + forward
    ! substitution L^{-1} to b.  On exit b holds L^{-1} P b_0.
    call y12mc(ndof, nz, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
        aflag, iflag, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12mc returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if

    call system_clock(t1)
    t_factor = t1 - t0

    ! =====================================================================
    ! 3. Time loop: one y12md call per step
    !    Step 1: IFLAG(5) = 2 -- y12mc already applied L, so only U^{-1}
    !    Steps 2..nsteps: IFLAG(5) = 3 -- apply full L^{-1} U^{-1}
    ! =====================================================================
    call system_clock(t0)

    ! Step 1 (iflag(5) still = 2 from y12mc): backward substitution only
    call y12md(ndof, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12md step 1 returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if

    ! Steps 2..nsteps: b = u^{step-1} is the new RHS -> y12md gives u^{step}
    iflag(5) = 3
    do step = 2, nsteps
      call y12md(ndof, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
      if (ifail /= 0) then
        write(error_unit, '(a,i0,a,i0)') &
            'ERROR: y12md step ', step, ' returned IFAIL = ', ifail
        stop 1, quiet = .true.
      end if
    end do

    call system_clock(t1)
    t_steps = t1 - t0

    ! =====================================================================
    ! 4. Error at t = T
    ! =====================================================================
    err_max = 0.0_dp
    do j = 2, N - 1
      do i = 2, N - 1
        k = (j - 2) * ndof_1d + (i - 2) + 1
        err_max = max(err_max, &
            abs(b(k) - u_exact(real(i - 1, dp) * h, real(j - 1, dp) * h, T, kappa)))
      end do
    end do
    write(output_unit, '(a,es12.4)') 'Max error at T: ', err_max

    ! =====================================================================
    ! 5. Write solution
    ! =====================================================================
    open(newunit=funit, file=outfile, status='unknown', action='write')
    write(funit, '(a)') '# 2D heat equation -- backward Euler implicit time-stepping'
    write(funit, '(a,i0,a,es10.3,a,es10.3)') &
        '# N=', N, '  T=', T, '  kappa=', kappa
    write(funit, '(a)') '# Columns: x  y  u_numerical  u_exact'
    do j = 1, N
      do i = 1, N
        if (i == 1 .or. i == N .or. j == 1 .or. j == N) then
          write(funit, '(4(1x,es14.6))') &
              real(i - 1, dp) * h, real(j - 1, dp) * h, 0.0_dp, &
              u_exact(real(i - 1, dp) * h, real(j - 1, dp) * h, T, kappa)
        else
          k = (j - 2) * ndof_1d + (i - 2) + 1
          write(funit, '(4(1x,es14.6))') &
              real(i - 1, dp) * h, real(j - 1, dp) * h, b(k), &
              u_exact(real(i - 1, dp) * h, real(j - 1, dp) * h, T, kappa)
        end if
      end do
      write(funit, *)
    end do
    close(funit)
    write(output_unit, '(a)') 'Solution written to ' // trim(outfile)

    ! =====================================================================
    ! 6. Timing summary
    ! =====================================================================
    s_assemble = real(t_assemble, dp) / real(clock_rate, dp)
    s_factor   = real(t_factor,   dp) / real(clock_rate, dp)
    s_steps    = real(t_steps,    dp) / real(clock_rate, dp)

    write(output_unit, '(/,a)') 'Timing summary (LU reuse advantage):'
    write(output_unit, '(a,g14.6,a)') '  Assembly + IC         : ', s_assemble, ' s'
    write(output_unit, '(a,g14.6,a)') '  LU factorization      : ', s_factor,   ' s  (done ONCE)'
    write(output_unit, '(a,i0,a,g14.6,a)') &
        '  Time loop (', nsteps, ' solves) : ', s_steps, ' s'
    if (s_steps > 0.0_dp) then
      write(output_unit, '(a,g14.6,a)') &
          '  Per-step solve        : ', s_steps / real(nsteps, dp), ' s/step'
    end if
    write(output_unit, '(a,g14.6,a)') &
        '  Naive estimate        : ', s_factor * real(nsteps, dp), &
        ' s  (if factorized each step)'

  end subroutine run
end module heat_implicit_solver

! ==============================================================
! Main program: parse CLI, call run.
! ==============================================================
program heat_implicit
  use heat_implicit_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N, nsteps, ios, iarg, nargs, pos_count
  real(kind(1.0d0)) :: T, kappa
  character(len=256) :: arg, outfile

  N       = 42
  nsteps  = 200
  T       = 0.1d0
  kappa   = 1.0d0
  outfile = 'heat_implicit.dat'

  nargs     = command_argument_count()
  iarg      = 1
  pos_count = 0

  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: heat_implicit [--help] [N] [nsteps] [T] [kappa] [output_file]'
      write(*, '(a)') '  N           total grid size (>= 3, default 42)'
      write(*, '(a)') '  nsteps      number of time steps (default 200)'
      write(*, '(a)') '  T           final time (default 0.1)'
      write(*, '(a)') '  kappa       thermal diffusivity (default 1.0)'
      write(*, '(a)') '  output_file output file (default: heat_implicit.dat)'
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
        if (ios /= 0 .or. T <= 0.0d0) then
          write(error_unit, '(a)') 'ERROR: T must be a positive real'
          stop 1
        end if
      case (4)
        read(arg, *, iostat=ios) kappa
        if (ios /= 0 .or. kappa <= 0.0d0) then
          write(error_unit, '(a)') 'ERROR: kappa must be a positive real'
          stop 1
        end if
      case (5)
        outfile = trim(arg)
      end select
    end if

    iarg = iarg + 1
  end do

  call run(N, nsteps, T, kappa, trim(outfile))

end program heat_implicit
