! SPDX-License-Identifier: GPL-2.0-only
!
! conv_diff_cds.f90 - Steady-state advection-diffusion on a unit square
!
! Solves the Convection-Diffusion equation using Central Differences:
!
!   u*phi_x + v*phi_y - alpha * (phi_xx + phi_yy) = f(x,y)   on (0,1) x (0,1)
!
! Flow field (Stagnation flow):
!   u = x,   v = -y
!
! Test Cases
! ----------
! CASE 1: Pure Dirichlet (Trigonometric MMS)
!   Exact: phi(x,y) = sin(pi * x) * cos(pi * y)
!   BCs:   phi = Exact on all boundaries (North, South, East, West)
!
! CASE 2: Mixed Dirichlet/Neumann (Trigonometric-Polynomial MMS)
!   Exact: phi(x,y) = cos^2(pi * x / 2) * (1 - y^2)
!   BCs:   North (y=1): phi = 0                      (Inlet)
!          West  (x=0): phi = 1 - y^2                (Source)
!          South (y=0): d(phi)/dy = 0                (Symmetry)
!          East  (x=1): d(phi)/dx = 0                (Zero Gradient)
!
! Discretisation
! --------------
! Total grid N x N. We solve for ALL nodes (interior + boundary).
!
! 5-point Central Differencing Scheme (2nd-Order, O(h^2)):
!   Advection: u_i * (phi_{i+1} - phi_{i-1}) / 2h
!   Diffusion: -alpha * (phi_{i+1} - 2phi_i + phi_{i-1}) / h^2
!
! For Neumann boundaries (South, East), a ghost-cell mirroring approach
! is used to fold the boundary coefficients back into the interior,
! preserving 2nd-order accuracy without breaking matrix symmetry.
!
! Usage
! -----
!   conv_diff_cds [--help] [--case 1|2] [N] [alpha] [output_file]
!

module conv_diff_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  ! === CASE 1: Original Trigonometric MMS (Pure Dirichlet) ===
  elemental function u_ex_1(x, y) result(phi)
    real(dp), intent(in) :: x, y
    real(dp) :: phi
    phi = sin(pi * x) * cos(pi * y)
  end function u_ex_1

  elemental function f_src_1(x, y, alpha) result(f)
    real(dp), intent(in) :: x, y, alpha
    real(dp) :: f
    f = pi * (x * cos(pi * x) * cos(pi * y) + y * sin(pi * x) * sin(pi * y)) &
        + 2.0_dp * alpha * pi**2 * sin(pi * x) * cos(pi * y)
  end function f_src_1

  ! ! === CASE 2: Ferziger/Peric Polynomial MMS (Mixed Boundaries) ===
  ! elemental function u_ex_2(x, y) result(phi)
  !   real(dp), intent(in) :: x, y
  !   real(dp) :: phi
  !   phi = (1.0_dp - x)**2 * (1.0_dp - y**2)
  ! end function u_ex_2

  ! elemental function f_src_2(x, y, alpha) result(f)
  !   real(dp), intent(in) :: x, y, alpha
  !   real(dp) :: f
  !   f = 2.0_dp * y**2 * (1.0_dp - x)**2 &
  !       - 2.0_dp * x * (1.0_dp - x) * (1.0_dp - y**2) &
  !       + 2.0_dp * alpha * ((1.0_dp - x)**2 - (1.0_dp - y**2))
  ! end function f_src_2

  ! === CASE 2: Mixed Boundaries (Trigonometric-Polynomial MMS) ===
  elemental function u_ex_2(x, y) result(phi)
    real(dp), intent(in) :: x, y
    real(dp) :: phi
    phi = (cos(pi * x / 2.0_dp)**2) * (1.0_dp - y**2)
  end function u_ex_2

  elemental function f_src_2(x, y, alpha) result(f)
    real(dp), intent(in) :: x, y, alpha
    real(dp) :: f
    real(dp) :: c_half_sq, one_minus_y2

    ! Precompute shared terms for efficiency
    c_half_sq = cos(pi * x / 2.0_dp)**2
    one_minus_y2 = 1.0_dp - y**2

    f = - (pi * x / 2.0_dp) * sin(pi * x) * one_minus_y2 &
        + 2.0_dp * y**2 * c_half_sq &
        + alpha * (pi**2 / 2.0_dp) * cos(pi * x) * one_minus_y2 &
        + 2.0_dp * alpha * c_half_sq
  end function f_src_2

  ! ---------------------------------------------------------------
  ! Root-mean-square of element-wise difference
  ! ---------------------------------------------------------------
  function rms_diff(a, b) result(rms)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: rms
    rms = norm2(a - b) / sqrt(real(size(a), dp))
  end function rms_diff

! ---------------------------------------------------------------
  ! Assemble the sparse matrix for all N x N nodes
  ! ---------------------------------------------------------------
  subroutine assemble_system(N, alpha, test_case, a, rnr, snr, b, nz)
    integer, intent(in) :: N, test_case
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: a(:), b(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    integer :: i, j, k_row
    real(dp) :: h, x, y, u, v
    real(dp) :: A_P, A_E, A_W, A_N, A_S
    logical :: is_dirichlet

    h = 1.0_dp / real(N - 1, dp)
    nz = 0
    b = 0.0_dp

    do j = 1, N
      y = real(j - 1, dp) * h
      v = -y
      do i = 1, N
        x = real(i - 1, dp) * h
        u = x
        k_row = (j - 1) * N + i

        is_dirichlet = .false.

        ! 1. Identify boundary nodes based on the test case
        if (test_case == 1) then
           if (i == 1 .or. i == N .or. j == 1 .or. j == N) is_dirichlet = .true.
        else if (test_case == 2) then
           if (i == 1 .or. j == N) is_dirichlet = .true.
        end if

        ! Calculate the base interior diagonal for scaling consistency
        A_P = 4.0_dp * alpha

        if (is_dirichlet) then
           ! Dirichlet Equation: A_P * phi = A_P * Exact
           call add(k_row, k_row, A_P)
           if (test_case == 1) b(k_row) = A_P * u_ex_1(x, y)
           if (test_case == 2) b(k_row) = A_P * u_ex_2(x, y)
        else
           ! 2. Interior & Neumann Equation (CDS)
           A_E = -alpha + u * h / 2.0_dp
           A_W = -alpha - u * h / 2.0_dp
           A_N = -alpha + v * h / 2.0_dp
           A_S = -alpha - v * h / 2.0_dp

           ! 3. Apply Ghost Cell mirroring for Case 2 Neumann boundaries
           if (test_case == 2) then
              if (i == N) then ! East boundary (Zero Gradient)
                 A_W = A_W + A_E
                 A_E = 0.0_dp
              end if
              if (j == 1) then ! South boundary (Symmetry)
                 A_N = A_N + A_S
                 A_S = 0.0_dp
              end if
           end if

           ! Insert into Sparse Matrix Structure
           call add(k_row, k_row, A_P)
           if (abs(A_W) > 0 .and. i > 1) call add(k_row, k_row - 1, A_W)
           if (abs(A_E) > 0 .and. i < N) call add(k_row, k_row + 1, A_E)
           if (abs(A_S) > 0 .and. j > 1) call add(k_row, k_row - N, A_S)
           if (abs(A_N) > 0 .and. j < N) call add(k_row, k_row + N, A_N)

           ! Set Right Hand Side Force
           if (test_case == 1) b(k_row) = (h**2) * f_src_1(x, y, alpha)
           if (test_case == 2) b(k_row) = (h**2) * f_src_2(x, y, alpha)
        end if
      end do
    end do

  contains
    subroutine add(r, c, val)
      integer, intent(in) :: r, c
      real(dp), intent(in) :: val
      nz = nz + 1
      rnr(nz) = r; snr(nz) = c; a(nz) = val
    end subroutine add
  end subroutine assemble_system

  ! ---------------------------------------------------------------
  ! Assemble, solve and report
  ! ---------------------------------------------------------------
  subroutine run(N, alpha, test_case, outfile)
    use y12m, only: y12ma
    integer, intent(in) :: N, test_case
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: outfile

    integer :: ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h
    real(dp), allocatable :: b(:), a(:), pivot(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail

    real(dp), allocatable :: u_full(:,:), u_ex_full(:,:)
    integer :: i, j, k
    real(dp) :: err_max, err_rms
    integer :: funit

    ndof = N * N
    h    = 1.0_dp / real(N - 1, dp)

    write(output_unit, '(a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, '  (', ndof, ' DOFs)'
    write(output_unit, '(a,i0)')   'Test Case           : ', test_case
    write(output_unit, '(a,f8.4)') 'Alpha (Diffusivity) : ', alpha

    ! Allocate generous workspace for LU fill-in (general non-symmetric)
    nz_max = 5 * ndof
    nn     = 40 * ndof
    nn1    = nn
    iha    = ndof

    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
    allocate(u_full(N, N), u_ex_full(N, N))

    aflag = 0.0_dp; iflag = 0; ifail = 0

    ! 1. Assemble
    call assemble_system(N, alpha, test_case, a, rnr, snr, b, nz)

    ! 2. Solve
    call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)

    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1
    end if

    ! 3. Map back and compute errors over the ENTIRE grid
    err_max = 0.0_dp
    do j = 1, N
      do i = 1, N
        k = (j - 1) * N + i
        u_full(i, j) = b(k)

        if (test_case == 1) u_ex_full(i, j) = u_ex_1(real(i - 1, dp) * h, real(j - 1, dp) * h)
        if (test_case == 2) u_ex_full(i, j) = u_ex_2(real(i - 1, dp) * h, real(j - 1, dp) * h)
      end do
    end do

    err_max = maxval(abs(u_full - u_ex_full))
    err_rms = rms_diff(u_full, u_ex_full)

    write(output_unit, '(a,es12.4)') 'Max pointwise error : ', err_max
    write(output_unit, '(a,es12.4)') 'RMS pointwise error : ', err_rms

    ! 4. Write Data
    open(newunit=funit, file=outfile, status='unknown', action='write')
    write(funit, '(a)') '# Columns: x   y   u_numerical   u_exact'
    do j = 1, N
      do i = 1, N
        write(funit, '(4(1x,es14.6))') real(i - 1, dp) * h, real(j - 1, dp) * h, &
            u_full(i, j), u_ex_full(i, j)
      end do
      write(funit, *)
    end do
    close(funit)
    write(output_unit, '(a)') 'Solution written to ' // trim(outfile)

  end subroutine run

end module conv_diff_solver

! ==============================================================
! Main program: parse CLI
! ==============================================================
program conv_diff_cds
  use conv_diff_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N, ios, iarg, nargs, pos_count
  integer :: test_case
  real(8) :: alpha
  character(len=256) :: arg, outfile

  ! Defaults
  N         = 42
  alpha     = 0.1d0
  test_case = 1
  outfile   = 'conv_diff.dat'

  nargs = command_argument_count()
  iarg = 1
  pos_count = 0

  ! Parse arguments (handles optional flags alongside positional inputs)
  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: conv_diff_cds [--help] [--case 1|2] [N] [alpha] [output_file]'
      write(*, '(a)') '  --case 1 : Pure Dirichlet (Trigonometric)'
      write(*, '(a)') '  --case 2 : Mixed Dirichlet/Neumann (Ferziger/Peric)'
      stop 0

    else if (trim(arg) == '--case') then
      iarg = iarg + 1
      if (iarg <= nargs) then
        call get_command_argument(iarg, arg)
        read(arg, *, iostat=ios) test_case
        if (ios /= 0 .or. (test_case /= 1 .and. test_case /= 2)) then
           write(error_unit,*) "ERROR: --case must be 1 or 2."
           stop 1
        end if
      end if

    else
      ! Positional arguments
      pos_count = pos_count + 1
      if (pos_count == 1) read(arg, *, iostat=ios) N
      if (pos_count == 2) read(arg, *, iostat=ios) alpha
      if (pos_count == 3) outfile = trim(arg)
    end if

    iarg = iarg + 1
  end do

  if (N < 3) then
    write(error_unit, '(a)') 'ERROR: N must be an integer >= 3'
    stop 1
  end if

  call run(N, alpha, test_case, trim(outfile))

end program conv_diff_cds