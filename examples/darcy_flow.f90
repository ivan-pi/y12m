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
!    Staggered MAC Grid (Equivalent to RT0-P0 Mixed FEM)
!    Domain: [0, 1] x [0, 1]
!    Origin (0,0) is a stagnation point.
!
! Features:
!    - Cell-by-cell interleaved DOF ordering to minimize LU fill-in bandwidth.
!    - Explicit pressure anchoring to resolve the hydrostatic nullspace.
!    - Direct enforcement of Pressure-Flux boundary relations.
!
! Solver: y12ma (double precision sparse direct).
!
! ==============================================================================
! Mixed Darcy Flow: MAC Staggered Grid Architecture
! ==============================================================================
!
! To satisfy the LBB (inf-sup) stability condition and prevent unphysical
! pressure checkerboarding, this solver uses a staggered Marker-and-Cell (MAC)
! grid. This is mathematically equivalent to RT0-P0 mixed finite elements.
!
! Physical Shape of a Single Cell (i, j):
!
!         (v-velocity lives on TOP horizontal face)
!                         v(i, j+1)
!                      ------^------
!                      |           |
! (u-velocity on       |           |      (u-velocity on
!  LEFT vertical -> u(i,j)  p(i,j) -> u(i+1, j)  RIGHT vertical face)
!  face)               |  (center)  |
!                      |           |
!                      ------^------
!                         v(i, j)
!         (v-velocity lives on BOTTOM horizontal face)
!
! Degrees of Freedom (DOFs):
!    - Pressure (p) is a cell-centered scalar.
!    - Velocity (u, v) is a face-centered vector normal to the cell edges.
!
! Matrix Assembly (Interleaved):
!    To prevent massive LU fill-in, the DOFs are numbered cell-by-cell.
!    For a given cell (i,j), the storage order is:
!      1. Left face U-velocity
!      2. Bottom face V-velocity
!      3. Cell center Pressure
!    Right and Top boundaries of the global domain are appended at the end.
!
! ==============================================================================
! Boundary Treatment & Pressure Anchoring
! ==============================================================================
! 1. Boundary Conditions:
!    Pressure is imposed weakly at the boundaries through the momentum equations.
!    For a boundary face (e.g., x=0), the gradient is approximated using a
!    half-step: u + (p_cell - p_ext)/(h/2) = 0. This ensures a consistent flux
!    injection without the "ghost cell" overhead.
!
! 2. Nullspace Resolution:
!    A pure Darcy system with Neumann-like flux boundaries has a constant
!    pressure nullspace. To ensure a unique solution and prevent solver
!    instability, the first pressure DOF (1,1) is explicitly anchored to the
!    exact solution value, replacing its divergence equation.
!
! ==============================================================================
! Exact Solution: Cubic Stagnation Flow
! ==============================================================================
! This test case simulates incompressible potential flow near a stagnation point
! (the origin). It uses a cubic polynomial for pressure, which yields a perfectly
! smooth quadratic velocity field.
!
! Equations:
!    p(x,y) = x^3 - 3xy^2
!    u(x,y) = -dp/dx = -3x^2 + 3y^2
!    v(x,y) = -dp/dy = 6xy
!
! Verification of Incompressibility (Mass Conservation):
!    div(q) = du/dx + dv/dy = -6x + 6x = 0
!
! ==============================================================================
! Verification & Convergence (RT0-P0 / MAC Grid)
! ==============================================================================
! To validate the RT0-P0 theory, errors are measured directly at the
! native DOFs:
!
! 1. Velocity Error: Measured at face midpoints (primary DOFs).
!    On uniform grids, this avoids interpolation errors and typically yields
!    O(h^2) in the interior, dropping toward O(h^1) at the global boundaries
!    due to the half-step gradient approximation.
!
! 2. Pressure Error: Measured at cell centers.
!    Expected convergence is O(h^2) across the domain.
!
! ==============================================================================
! Solver Methodology: Saddle-Point vs. Schur Complement
! ==============================================================================
! NOTE: This driver solves the fully coupled, indefinite saddle-point system:
! [ M   B^T ] [ q ] = [ f ]
! [ B    0  ] [ p ] = [ g ]
!
! This is solved directly using y12ma. The system is scaled such that the
! momentum diagonal is O(h) and the divergence/gradient blocks are O(1),
! providing better numerical balance for the LU decomposition than h^2 scaling.
!
! ==============================================================================
! Indexing Strategy: The "Orphan" Boundary Faces
! ==============================================================================
! In a perfectly periodic domain (wrapping around like a torus), every cell
! would have exactly one left (U) face, one bottom (V) face, and one center (P).
! This would consume 3N^2 Degrees of Freedom (DOFs) in pure triplets.
!
! On a bounded square domain, we have "leftover" faces:
!    - The absolute Right wall has N extra U-faces.
!    - The absolute Top wall has N extra V-faces.
! Total DOFs = 3N^2 (Interior Triplets) + 2N (Boundary Orphans).
!
! How this is handled without destroying the LU bandwidth:
! 1. The first 3N^2 indices strictly interleave the interior/left/bottom
!    triplets: (U, V, P) for cell 1, (U, V, P) for cell 2, etc.
! 2. The 2N "orphan" boundary faces are appended to the very end. Because
!    these faces only couple with their immediate neighbor, they generate
!    negligible fill-in during the LU sweep.
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
