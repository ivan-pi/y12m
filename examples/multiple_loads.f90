! SPDX-License-Identifier: GPL-2.0-only
!
! multiple_loads.f90 - Poisson equation with multiple right-hand sides
!
! Solves a family of 2-D Poisson problems on the unit square:
!
!   -nabla^2 u_k = f_k(x, y)    on (0,1) x (0,1),  k = 1, ..., nrhs
!
! with homogeneous Dirichlet boundary conditions  u_k = 0  on the boundary.
!
! Modal load cases (f_k):
!
!   f_k(x, y) = sin(k pi x) sin(pi y),     k = 1, 2, ..., nrhs
!
! Exact solutions (closed form):
!
!   u_k(x, y) = sin(k pi x) sin(pi y) / ((k^2 + 1) pi^2)
!
! Derivation: substituting into the Poisson equation gives
!   -nabla^2 u_k = (k^2 + 1) pi^2 u_k = f_k,  hence u_k = f_k / ((k^2+1) pi^2).
!
! Discretisation
! --------------
! Total grid N x N (including boundary nodes), spacing h = 1/(N-1).
! Interior nodes: ndof = (N-2)^2.
! DOF index: k(i,j) = (j-2)*(N-2) + (i-2) + 1.
!
! The discrete system  K u_k = b_k  uses the standard 5-point stiffness:
!   K diagonal = 4,  off-diagonal = -1  (neighbours).
! The right-hand side is  b_k(i,j) = h^2 f_k(x_i, y_j).
!
! Solver strategy -- Case (ii): same matrix, many right-hand sides
! ----------------------------------------------------------------
! The stiffness matrix K is IDENTICAL for all nrhs load cases; only the
! right-hand side changes.
!
!   * y12mb + y12mc are called ONCE to factorize K.
!   * y12md is called once per load case (IFLAG(5) = 3) to solve the
!     triangular systems with each b_k.
!
! A naive approach would factorize K nrhs times (cost nrhs * O(nnz_LU)).
! Reusing the LU factorization reduces that to 1 * O(nnz_LU) + nrhs *
! O(nnz_LU) for the substitutions -- a substantial saving when nrhs is large.
!
! Usage
! -----
!   multiple_loads [--help] [N] [nrhs] [output_file]
!
!     N           total grid size (>= 3, default 42)
!     nrhs        number of load cases (>= 1, default 8)
!     output_file output data file name (default: multiple_loads.dat)
!

! ==============================================================
! Module multiple_loads_solver
! ==============================================================
module multiple_loads_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  !> Exact solution for load case k.
  elemental function u_exact_k(x, y, k) result(u)
    real(dp), intent(in) :: x, y
    integer, intent(in) :: k
    real(dp) :: u
    u = sin(real(k, dp) * pi * x) * sin(pi * y) &
        / (real(k * k + 1, dp) * pi**2)
  end function u_exact_k

  !> Right-hand side f_k(x,y) = sin(k pi x) sin(pi y).
  elemental function f_load_k(x, y, k) result(f)
    real(dp), intent(in) :: x, y
    integer, intent(in) :: k
    real(dp) :: f
    f = sin(real(k, dp) * pi * x) * sin(pi * y)
  end function f_load_k

  !> Assemble the 5-point Poisson stiffness matrix K in COO format.
  !> K has diagonal = 4 and off-diagonal = -1 for each interior neighbour.
  !> Entries are emitted in a fixed row-major order.
  subroutine assemble_stiffness(N, a, rnr, snr, nz)
    integer, intent(in) :: N
    real(dp), intent(out) :: a(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    integer :: ndof_1d, i, j, k_row

    ndof_1d = N - 2
    nz = 0

    do j = 2, N - 1
      do i = 2, N - 1
        k_row = (j - 2) * ndof_1d + (i - 2) + 1

        ! Diagonal
        nz = nz + 1
        rnr(nz) = k_row
        snr(nz) = k_row
        a(nz) = 4.0_dp

        ! West neighbour (i-1, j)
        if (i > 2) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 2) * ndof_1d + (i - 3) + 1
          a(nz) = -1.0_dp
        end if

        ! East neighbour (i+1, j)
        if (i < N - 1) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 2) * ndof_1d + (i - 1) + 1
          a(nz) = -1.0_dp
        end if

        ! South neighbour (i, j-1)
        if (j > 2) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 3) * ndof_1d + (i - 2) + 1
          a(nz) = -1.0_dp
        end if

        ! North neighbour (i, j+1)
        if (j < N - 1) then
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = (j - 1) * ndof_1d + (i - 2) + 1
          a(nz) = -1.0_dp
        end if

      end do
    end do
  end subroutine assemble_stiffness

  ! ---------------------------------------------------------------
  ! Main driver: factorize K once, solve nrhs times
  ! ---------------------------------------------------------------
  subroutine run(N, nrhs, outfile)
    use y12m, only: y12mb, y12mc, y12md
    integer, intent(in) :: N, nrhs
    character(len=*), intent(in) :: outfile

    integer :: ndof_1d, ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h, h2
    real(dp), allocatable :: a(:), pivot(:), b(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp), allocatable :: b_all(:,:), u_all(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail
    integer :: i, j, k, kload, funit
    integer :: t0, t1, clock_rate
    integer :: t_assemble, t_factor, t_solve
    real(dp) :: s_assemble, s_factor, s_solve
    real(dp) :: err_max, err_k

    ndof_1d = N - 2
    ndof = ndof_1d * ndof_1d
    h = 1.0_dp / real(N - 1, dp)
    h2 = h * h

    write(output_unit, '(a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
    write(output_unit, '(a,i0)') 'Load cases    : ', nrhs

    ! Workspace
    nz_max = 5 * ndof
    nn = max(3 * nz_max, 2 * (ndof_1d + 2) * ndof)
    nn1 = nn
    iha = ndof
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
    allocate(b_all(ndof, nrhs), u_all(ndof, nrhs))

    call system_clock(count_rate=clock_rate)

    ! =====================================================================
    ! 1. Assemble stiffness K and all nrhs right-hand sides
    ! =====================================================================
    call system_clock(t0)
    call assemble_stiffness(N, a, rnr, snr, nz)

    do kload = 1, nrhs
      do j = 2, N - 1
        do i = 2, N - 1
          k = (j - 2) * ndof_1d + (i - 2) + 1
          b_all(k, kload) = h2 * f_load_k(real(i - 1, dp) * h, &
              real(j - 1, dp) * h, kload)
        end do
      end do
    end do
    call system_clock(t1)
    t_assemble = t1 - t0

    ! =====================================================================
    ! 2. Factorize K ONCE with IFLAG(5) = 2 (keep L for reuse)
    !    Load the first RHS into b so y12mc can forward-substitute it.
    ! =====================================================================
    call system_clock(t0)

    b = b_all(:, 1)

    iflag = 0
    iflag(2) = 3        ! Markowitz search width
    iflag(3) = 1        ! use column interchanges
    iflag(4) = 0        ! no structural reuse
    iflag(5) = 2        ! keep L (Case ii)
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

    ! y12mc factorizes K and applies L^{-1} P to the first RHS b.
    call y12mc(ndof, nz, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
        aflag, iflag, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12mc returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if

    call system_clock(t1)
    t_factor = t1 - t0

    ! =====================================================================
    ! 3. Solve for each load case
    !    Load case 1: IFLAG(5) = 2 (y12mc already applied L^{-1})
    !    Load cases 2..nrhs: IFLAG(5) = 3 (fresh RHS each time)
    ! =====================================================================
    call system_clock(t0)

    ! Solve for load case 1
    call y12md(ndof, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12md load 1 returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if
    u_all(:, 1) = b

    ! Solve for load cases 2..nrhs using the stored LU factorization
    iflag(5) = 3
    do kload = 2, nrhs
      b = b_all(:, kload)
      call y12md(ndof, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
      if (ifail /= 0) then
        write(error_unit, '(a,i0,a,i0)') &
            'ERROR: y12md load ', kload, ' returned IFAIL = ', ifail
        stop 1, quiet = .true.
      end if
      u_all(:, kload) = b
    end do

    call system_clock(t1)
    t_solve = t1 - t0

    ! =====================================================================
    ! 4. Error check
    ! =====================================================================
    write(output_unit, '(/,a)') 'Error table:'
    write(output_unit, '(a)') '  k   max|u_num - u_exact|'
    err_max = 0.0_dp
    do kload = 1, nrhs
      err_k = 0.0_dp
      do j = 2, N - 1
        do i = 2, N - 1
          k = (j - 2) * ndof_1d + (i - 2) + 1
          err_k = max(err_k, abs(u_all(k, kload) - &
              u_exact_k(real(i - 1, dp) * h, real(j - 1, dp) * h, kload)))
        end do
      end do
      write(output_unit, '(2x,i3,3x,es12.4)') kload, err_k
      err_max = max(err_max, err_k)
    end do
    write(output_unit, '(a,es12.4)') 'Max over all load cases: ', err_max

    ! =====================================================================
    ! 5. Write last solution to output file
    ! =====================================================================
    open(newunit=funit, file=outfile, status='unknown', action='write')
    write(funit, '(a)') '# Poisson equation -- multiple load cases'
    write(funit, '(a,i0,a,i0)') '# N=', N, '  nrhs=', nrhs
    write(funit, '(a)') '# Solution for last load case k=nrhs'
    write(funit, '(a)') '# Columns: x  y  u_numerical  u_exact'
    do j = 1, N
      do i = 1, N
        if (i == 1 .or. i == N .or. j == 1 .or. j == N) then
          write(funit, '(4(1x,es14.6))') &
              real(i - 1, dp) * h, real(j - 1, dp) * h, 0.0_dp, &
              u_exact_k(real(i - 1, dp) * h, real(j - 1, dp) * h, nrhs)
        else
          k = (j - 2) * ndof_1d + (i - 2) + 1
          write(funit, '(4(1x,es14.6))') &
              real(i - 1, dp) * h, real(j - 1, dp) * h, u_all(k, nrhs), &
              u_exact_k(real(i - 1, dp) * h, real(j - 1, dp) * h, nrhs)
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
    s_solve    = real(t_solve,    dp) / real(clock_rate, dp)

    write(output_unit, '(/,a)') 'Timing summary (LU reuse advantage):'
    write(output_unit, '(a,g14.6,a)') '  Assembly (K + all b_k) : ', s_assemble, ' s'
    write(output_unit, '(a,g14.6,a)') '  LU factorization       : ', s_factor,   ' s  (done ONCE)'
    write(output_unit, '(a,g14.6,a,i0,a)') &
        '  All solves (', nrhs, ' calls)    : ', s_solve, ' s'
    if (s_solve > 0.0_dp) then
      write(output_unit, '(a,g14.6,a)') &
          '  Per-solve cost         : ', s_solve / real(nrhs, dp), ' s/call'
    end if
    write(output_unit, '(a,g14.6,a)') &
        '  Naive estimate         : ', s_factor * real(nrhs, dp), &
        ' s  (if factorized each time)'

  end subroutine run
end module multiple_loads_solver

! ==============================================================
! Main program: parse CLI, call run.
! ==============================================================
program multiple_loads
  use multiple_loads_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N, nrhs, ios, iarg, nargs, pos_count
  character(len=256) :: arg, outfile

  N       = 42
  nrhs    = 8
  outfile = 'multiple_loads.dat'

  nargs     = command_argument_count()
  iarg      = 1
  pos_count = 0

  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: multiple_loads [--help] [N] [nrhs] [output_file]'
      write(*, '(a)') '  N           total grid size (>= 3, default 42)'
      write(*, '(a)') '  nrhs        number of load cases (>= 1, default 8)'
      write(*, '(a)') '  output_file output file (default: multiple_loads.dat)'
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
        read(arg, *, iostat=ios) nrhs
        if (ios /= 0 .or. nrhs < 1) then
          write(error_unit, '(a)') 'ERROR: nrhs must be a positive integer'
          stop 1
        end if
      case (3)
        outfile = trim(arg)
      end select
    end if

    iarg = iarg + 1
  end do

  call run(N, nrhs, trim(outfile))

end program multiple_loads
