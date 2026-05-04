! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: Gemini
!
! darcy_mixed.f90 - Mixed Darcy Flow near a Stagnation Point
!
! Solves the saddle-point system:
!    q + \nabla p = 0
!    \nabla \cdot q = 0
!
! Discretization:
!    Staggered Marker-and-Cell (MAC) Grid (Equivalent to RT0-P0 Mixed FEM)
!    Domain: [0, 1] x [0, 1] | Origin (0,0) is a stagnation point.
!
! Solver Architecture:
!    - System: Indefinite saddle-point block system.
!    - Ordering: Cell-by-cell interleaved DOFs to minimize LU bandwidth/fill-in.
!    - Linear Algebra: y12ma (sparse direct LU decomposition).
!
! ==============================================================================
! Implementation Logic & Numerical Scaling
! ==============================================================================
! 1. Scaling & Symmetry:
!    The system is scaled for numerical balance rather than pure symmetry.
!    Momentum rows are O(h) and Continuity rows are O(1), resulting in a matrix
!    of the form [hI, B^T; B, 0]. While this improves the condition number for
!    direct solvers, it renders the global system non-symmetric.
!
! 2. Nullspace & Pressure Pinning:
!    The hydrostatic pressure nullspace is resolved by pinning cell (1,1).
!    The divergence equation for this cell is replaced by the identity
!    p(1,1) = p_exact. While this technically violates local mass conservation
!    at one point, the global error remains O(h^2).
!
! 3. Strong Boundary Enforcement:
!    Boundary pressure is enforced "strongly" within the momentum equations
!    using a half-step one-sided gradient: u + (p_cell - p_ext)/(h/2) = 0.
!    This avoids ghost-cell overhead but localizes an O(h^1) error at the
!    boundary faces, which is characteristic of MAC/RT0 on Dirichlet domains.
!
! ==============================================================================
! MAC Staggered Grid Architecture
! ==============================================================================
! To satisfy the LBB (inf-sup) stability condition and prevent checkerboarding,
! variables are staggered. This is equivalent to Raviart-Thomas (RT0-P0) FEM.
!
! Physical Shape of a Single Cell (i, j):
!
!                    v-velocity (TOP face): v(i, j+1)
!                             ------^------
!                             |           |
!      u-velocity (LEFT)  -> u(i,j)  p(i,j) -> u(i+1, j)  u-velocity (RIGHT)
!                             |  (center)  |
!                             ------^------
!                    v-velocity (BOTTOM face): v(i, j)
!
! Matrix Assembly (Interleaved):
! To minimize LU fill-in, DOFs are grouped by cell. For cell (i,j):
!   1. U-velocity (Left Face)
!   2. V-velocity (Bottom Face)
!   3. Pressure   (Cell Center)
!
! "Orphan" boundary faces (the absolute Right and Top edges) are appended to
! the end of the global array to preserve the triplet structure of the interior.
!
! ==============================================================================
! Exact Solution: Cubic Stagnation Flow
! ==============================================================================
! Simulates incompressible potential flow near a stagnation point (0,0).
!
! Equations:
!    p(x,y) = x^3 - 3xy^2
!    u(x,y) = -3x^2 + 3y^2
!    v(x,y) = 6xy
!
! Verification:
!    div(q) = du/dx + dv/dy = -6x + 6x = 0 (Exactly divergence-free).
!
! ==============================================================================
! Convergence Metrics
! ==============================================================================
! Errors are measured at primary DOF locations to verify RT0 theory:
!
! 1. Velocity (L_inf): Measured at face midpoints. Expected O(h^2) in the
!    interior, dropping to O(h^1) globally due to boundary layer pollution.
!
! 2. Pressure (L_inf): Measured at cell centers. Expected O(h^2).
! ==============================================================================

module darcy_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)

contains

  elemental function p_exact(x, y) result(p)
    real(dp), intent(in) :: x, y
    real(dp) :: p
    p = x**3 - 3.0_dp * x * y**2
  end function p_exact

  elemental function u_exact(x, y) result(u)
    real(dp), intent(in) :: x, y
    real(dp) :: u
    u = -3.0_dp * x**2 + 3.0_dp * y**2
  end function u_exact

  elemental function v_exact(x, y) result(v)
    real(dp), intent(in) :: x, y
    real(dp) :: v
    v = 6.0_dp * x * y
  end function v_exact

  ! ---------------------------------------------------------------
  ! Index Mapping: Interleaved (Cell-by-Cell) Ordering
  ! ---------------------------------------------------------------
  pure function idx_u(i, j, N) result(idx)
    integer, intent(in) :: i, j, N
    integer :: idx
    if (i <= N) then
      idx = (j - 1) * N * 3 + (i - 1) * 3 + 1
    else
      idx = 3 * N**2 + j
    end if
  end function idx_u

  pure function idx_v(i, j, N) result(idx)
    integer, intent(in) :: i, j, N
    integer :: idx
    if (j <= N) then
      idx = (j - 1) * N * 3 + (i - 1) * 3 + 2
    else
      idx = 3 * N**2 + N + i
    end if
  end function idx_v

  pure function idx_p(i, j, N) result(idx)
    integer, intent(in) :: i, j, N
    integer :: idx
    idx = (j - 1) * N * 3 + (i - 1) * 3 + 3
  end function idx_p

  ! ---------------------------------------------------------------
  ! Saddle-Point Matrix Assembly
  ! ---------------------------------------------------------------
  subroutine assemble_system(N, h, a, rnr, snr, nz, b, nz_max)
    integer, intent(in) :: N, nz_max
    real(dp), intent(in) :: h
    real(dp), intent(out) :: a(:), b(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    integer :: i, j, row
    real(dp) :: x, y

    nz = 0
    b = 0.0_dp

    ! 1. U-Velocity Faces (Momentum X: u + dp/dx = 0)
    do j = 1, N
      y = (real(j, dp) - 0.5_dp) * h
      do i = 1, N + 1
        row = idx_u(i, j, N)
        x = real(i - 1, dp) * h

        if (i == 1) then
          ! Left Boundary: u + (p_cell - p_ext)/(h/2) = 0 -> (h/2)u + p_cell = p_ext
          call add_nz_local(row, row, h / 2.0_dp)
          call add_nz_local(row, idx_p(i, j, N), 1.0_dp)
          b(row) = p_exact(x, y)
        else if (i == N + 1) then
          ! Right Boundary: u + (p_ext - p_prev)/(h/2) = 0 -> (h/2)u - p_prev = -p_ext
          call add_nz_local(row, row, h / 2.0_dp)
          call add_nz_local(row, idx_p(i - 1, j, N), -1.0_dp)
          b(row) = -p_exact(x, y)
        else
          ! Interior: u + (p_i - p_i-1)/h = 0 -> h*u + p_i - p_i-1 = 0
          call add_nz_local(row, row, h)
          call add_nz_local(row, idx_p(i - 1, j, N), -1.0_dp)
          call add_nz_local(row, idx_p(i, j, N), 1.0_dp)
        end if
      end do
    end do

    ! 2. V-Velocity Faces (Momentum Y: v + dp/dy = 0)
    do j = 1, N + 1
      y = real(j - 1, dp) * h
      do i = 1, N
        row = idx_v(i, j, N)
        x = (real(i, dp) - 0.5_dp) * h

        if (j == 1) then
          ! Bottom Boundary
          call add_nz_local(row, row, h / 2.0_dp)
          call add_nz_local(row, idx_p(i, j, N), 1.0_dp)
          b(row) = p_exact(x, y)
        else if (j == N + 1) then
          ! Top Boundary
          call add_nz_local(row, row, h / 2.0_dp)
          call add_nz_local(row, idx_p(i, j - 1, N), -1.0_dp)
          b(row) = -p_exact(x, y)
        else
          ! Interior
          call add_nz_local(row, row, h)
          call add_nz_local(row, idx_p(i, j - 1, N), -1.0_dp)
          call add_nz_local(row, idx_p(i, j, N), 1.0_dp)
        end if
      end do
    end do

    ! 3. Pressure Cells (Continuity: div q = 0)
    do j = 1, N
      do i = 1, N
        row = idx_p(i, j, N)

        ! Pin pressure at (1,1) to resolve nullspace
        if (i == 1 .and. j == 1) then
           call add_nz_local(row, row, 1.0_dp)
           b(row) = p_exact(0.5_dp*h, 0.5_dp*h)
        else
           ! Standard divergence: (u_east - u_west)/h + (v_north - v_south)/h = 0
           ! Multiplied by h to maintain symmetry with momentum blocks
           call add_nz_local(row, idx_u(i, j, N), -1.0_dp)
           call add_nz_local(row, idx_u(i + 1, j, N), 1.0_dp)
           call add_nz_local(row, idx_v(i, j, N), -1.0_dp)
           call add_nz_local(row, idx_v(i, j + 1, N), 1.0_dp)
           b(row) = 0.0_dp
        end if
      end do
    end do

  contains
    subroutine add_nz_local(r, c, val)
      integer, intent(in) :: r, c
      real(dp), intent(in) :: val
      if (nz >= nz_max) then
         write(error_unit,*) "FATAL: nz_max exceeded during assembly."
         stop 1
      end if
      nz = nz + 1
      rnr(nz) = r
      snr(nz) = c
      a(nz) = val
    end subroutine add_nz_local
  end subroutine assemble_system

  ! ---------------------------------------------------------------
  ! Main Driver
  ! ---------------------------------------------------------------
  subroutine run(N, outfile)
    use y12m, only: y12ma
    integer, intent(in) :: N
    character(len=*), intent(in) :: outfile

    integer :: ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h, err_p, err_u, err_v, u_val, v_val, p_val, x_val, y_val
    real(dp), allocatable :: b(:), a(:), pivot(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail, i, j, funit

    ndof = 3 * N**2 + 2 * N
    h    = 1.0_dp / real(N, dp)

    nz_max = 12 * N**2 + 8 * N
    nn     = 30 * nz_max ! Generous buffer for LU fill-in
    nn1    = nn
    iha    = ndof
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))

    ! 1. Assemble
    call assemble_system(N, h, a, rnr, snr, nz, b, nn)

    ! 2. Solve
    aflag = 0.0_dp; iflag = 0; ifail = 0
    call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)

    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1
    end if

    ! 3. Error Analysis
    err_p = 0.0_dp; err_u = 0.0_dp; err_v = 0.0_dp

    ! Validate Primary DOFs (Face Midpoints for Velocity)
    do j = 1, N
      y_val = (real(j, dp) - 0.5_dp) * h
      do i = 1, N + 1
        x_val = real(i - 1, dp) * h
        u_val = b(idx_u(i, j, N))
        err_u = max(err_u, abs(u_val - u_exact(x_val, y_val)))
      end do
    end do

    do j = 1, N + 1
      y_val = real(j - 1, dp) * h
      do i = 1, N
        x_val = (real(i, dp) - 0.5_dp) * h
        v_val = b(idx_v(i, j, N))
        err_v = max(err_v, abs(v_val - v_exact(x_val, y_val)))
      end do
    end do

    ! Pressure and Data Export
    open(newunit=funit, file=outfile, status='unknown', action='write')
    do j = 1, N
      y_val = (real(j, dp) - 0.5_dp) * h
      do i = 1, N
        x_val = (real(i, dp) - 0.5_dp) * h
        p_val = b(idx_p(i, j, N))
        err_p = max(err_p, abs(p_val - p_exact(x_val, y_val)))
        write(funit, '(3(1x,es14.6))') x_val, y_val, p_val
      end do
      write(funit, *)
    end do
    close(funit)

    write(output_unit, '(a,i0,a)') 'Results for N=', N, ':'
    write(output_unit, '(a,es12.4)') '  Max L_inf Error (P) : ', err_p
    write(output_unit, '(a,es12.4)') '  Max L_inf Error (U) : ', err_u
    write(output_unit, '(a,es12.4)') '  Max L_inf Error (V) : ', err_v

  end subroutine run
end module darcy_solver

program darcy_driver
  use darcy_solver, only: run
  implicit none
  integer :: N, ios
  character(len=256) :: arg, outfile

  N = 32; outfile = 'darcy.dat'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *, iostat=ios) N
  end if
  if (command_argument_count() >= 2) call get_command_argument(2, outfile)

  call run(N, trim(outfile))
end program darcy_driver
