! SPDX-License-Identifier: GPL-2.0-only
!
! poisson_9pt.f90 - Steady-state heat diffusion on a unit square
!
! Solves the 2-D Laplace equation
!
!   Δu = 0   on  Ω = (0,1) × (0,1)
!
! with boundary conditions (a plate with a sinusoidal hot top edge):
!
!   u(x, 0) = 0           bottom  (cold)
!   u(x, 1) = sin(π x)    top     (sinusoidal hot edge)
!   u(0, y) = 0           left    (cold)
!   u(1, y) = 0           right   (cold)
!
! Exact solution:
!   u(x, y) = sin(π x) sinh(π y) / sinh(π)
!
! Discretisation
! --------------
! An N×N interior grid with uniform spacing h = 1/(N+1) is used.
! Interior grid points are indexed (i,j), i=1..N (x), j=1..N (y), with
! global DOF index k = (j-1)*N + i.
!
! The 9-point isotropic Laplacian (Mehrstellen) stencil is:
!
!   [1  4  1]         1
!   [4 -20  4]  ×  ———————
!   [1  4  1]        6h²
!
! so the assembled linear system is:
!   -20 u_{i,j}
!    +4 (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1})
!    +1 (u_{i+1,j+1} + u_{i-1,j+1} + u_{i+1,j-1} + u_{i-1,j-1}) = 0
!
! Boundary values are moved to the right-hand side.
!
! Solver
! ------
! y12ma (resolving to y12maf, the double-precision driver) performs
! Gaussian elimination without iterative refinement.
!
! Output
! ------
! Writes "poisson_9pt.dat" with columns:
!    x   y   u_numerical   u_exact
! for all grid points (including boundaries), with blank lines between
! scanlines so that gnuplot's "pm3d map" mode works directly.
!
! Usage after building:
!   ./poisson_9pt
!   gnuplot plot_poisson.gp      # produces poisson_9pt.png
!
program poisson_9pt
  use, intrinsic :: iso_fortran_env, only: dp => real64, &
      output_unit, error_unit
  use y12m, only: y12ma
  implicit none

  ! ---------------------------------------------------------------
  ! Grid parameters – change N to refine the grid
  ! ---------------------------------------------------------------
  integer, parameter :: N = 30              ! interior points per direction
  integer, parameter :: ndof = N * N        ! total degrees of freedom
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: h = 1.0_dp / real(N + 1, dp)

  ! ---------------------------------------------------------------
  ! Sparse storage
  !
  ! The 9-point stencil produces at most 9·ndof non-zeros.  The y12m
  ! solver maintains LU factors in both row-ordered (a, snr) and
  ! column-ordered (rnr) format, both of which grow as fill-in is
  ! generated.  For a 2-D Poisson matrix with natural ordering the
  ! half-bandwidth is N+1, so the LU fill-in is at most (N+1)·ndof
  ! entries per triangle.  Both nn (row list) and nn1 (column list)
  ! must accommodate the full LU factor.  The constant 10 below
  ! provides a generous safety margin.
  ! ---------------------------------------------------------------
  integer, parameter :: nnz_orig = 9 * ndof         ! non-zeros in A
  integer, parameter :: nn  = 10 * (N + 1) * ndof   ! a, snr workspace
  integer, parameter :: nn1 = 10 * (N + 1) * ndof   ! rnr workspace
  integer, parameter :: iha = ndof                   ! leading dim of ha

  ! ---------------------------------------------------------------
  ! Solver arrays (declared at program level → stored in data segment)
  ! ---------------------------------------------------------------
  real(dp) :: a(nn)
  real(dp) :: pivot(ndof)
  real(dp) :: aflag(8)
  real(dp) :: b(ndof)
  integer :: snr(nn)
  integer :: rnr(nn1)
  integer :: ha(iha, 11)
  integer :: iflag(10)
  integer :: nz, ifail

  ! ---------------------------------------------------------------
  ! Local variables
  ! ---------------------------------------------------------------
  character(len=*), parameter :: datafmt = '(4(1x,es14.6))'
  integer :: i, j, k, funit
  real(dp) :: xi, yj, u_exact, err
  real(dp) :: err_max, err_l2

  ! ==============================================================
  ! 1.  Assemble the sparse linear system  A u = b
  ! ==============================================================
  nz = 0
  b(1:ndof) = 0.0_dp

  do j = 1, N           ! y-direction (row) index
    do i = 1, N         ! x-direction (column) index
      k  = (j - 1) * N + i    ! global DOF index
      xi = real(i, dp) * h    ! x-coordinate
      yj = real(j, dp) * h    ! y-coordinate

      ! --- Centre: coefficient -20 ---
      nz = nz + 1
      rnr(nz) = k
      snr(nz) = k
      a(nz) = -20.0_dp

      ! --- East (i+1, j): coefficient +4 ---
      if (i < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + 1
        a(nz) = 4.0_dp
      end if
      ! Right boundary u=0: no contribution to b

      ! --- West (i-1, j): coefficient +4 ---
      if (i > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - 1
        a(nz) = 4.0_dp
      end if
      ! Left boundary u=0: no contribution to b

      ! --- North (i, j+1): coefficient +4 ---
      if (j < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + N
        a(nz) = 4.0_dp
      else
        ! Top boundary: u(x, 1) = sin(pi * xi)
        b(k) = b(k) - 4.0_dp * sin(pi * xi)
      end if

      ! --- South (i, j-1): coefficient +4 ---
      if (j > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - N
        a(nz) = 4.0_dp
      end if
      ! Bottom boundary u=0: no contribution to b

      ! --- North-East (i+1, j+1): coefficient +1 ---
      if (i < N .and. j < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + N + 1
        a(nz) = 1.0_dp
      else if (i < N .and. j == N) then
        ! Top boundary at (i+1, N+1): u = sin(pi*(xi+h))
        b(k) = b(k) - sin(pi * (xi + h))
      end if
      ! Right boundary or top-right corner: u=0

      ! --- North-West (i-1, j+1): coefficient +1 ---
      if (i > 1 .and. j < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + N - 1
        a(nz) = 1.0_dp
      else if (i > 1 .and. j == N) then
        ! Top boundary at (i-1, N+1): u = sin(pi*(xi-h))
        b(k) = b(k) - sin(pi * (xi - h))
      end if
      ! Left boundary (sin(0)=0) or top-left corner: u=0

      ! --- South-East (i+1, j-1): coefficient +1 ---
      if (i < N .and. j > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - N + 1
        a(nz) = 1.0_dp
      end if
      ! Right or bottom boundary: u=0

      ! --- South-West (i-1, j-1): coefficient +1 ---
      if (i > 1 .and. j > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - N - 1
        a(nz) = 1.0_dp
      end if
      ! Left or bottom boundary: u=0

    end do
  end do

  ! ==============================================================
  ! 2.  Solve with y12ma (double-precision driver, no iter. refinement)
  ! ==============================================================
  ifail = 0
  call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)

  if (ifail /= 0) then
    write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
    stop 1
  end if
  ! The solution vector u is now stored in b(1:ndof).

  ! ==============================================================
  ! 3.  Compute pointwise errors against the exact solution
  ! ==============================================================
  err_max = 0.0_dp
  err_l2  = 0.0_dp
  do j = 1, N
    do i = 1, N
      k  = (j - 1) * N + i
      xi = real(i, dp) * h
      yj = real(j, dp) * h
      u_exact = sin(pi * xi) * sinh(pi * yj) / sinh(pi)
      err = abs(b(k) - u_exact)
      err_max = max(err_max, err)
      err_l2  = err_l2 + err**2
    end do
  end do
  err_l2 = sqrt(err_l2 / real(ndof, dp))

  write(output_unit, '(a,i0,a,i0,a,i0,a)') &
      'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
  write(output_unit, '(a,es12.4)') 'Max pointwise error  : ', err_max
  write(output_unit, '(a,es12.4)') 'RMS error            : ', err_l2
  write(output_unit, '(a)')        'Exact solution: u = sin(pi*x)*sinh(pi*y)/sinh(pi)'

  ! ==============================================================
  ! 4.  Write gnuplot data file
  !
  ! Format: one block per y-scanline, blocks separated by blank lines.
  ! Each line: x  y  u_numerical  u_exact
  ! Includes boundary points so the full domain [0,1]×[0,1] is covered.
  ! ==============================================================
  open(newunit=funit, file='poisson_9pt.dat', status='replace', action='write')
  write(funit, '(a)') &
      '# 9-point isotropic stencil, Laplace equation on unit square'
  write(funit, '(a)') &
      '# BCs: u = sin(pi*x) on top edge; u = 0 on left/right/bottom'
  write(funit, '(a)') &
      '# Columns: x   y   u_numerical   u_exact'
  write(funit, '(a)') '#'

  do j = 0, N + 1
    yj = real(j, dp) * h
    do i = 0, N + 1
      xi = real(i, dp) * h
      u_exact = sin(pi * xi) * sinh(pi * yj) / sinh(pi)
      if (i == 0 .or. i == N + 1 .or. j == 0) then
        ! Left, right, or bottom boundary: u = 0
        write(funit, datafmt) xi, yj, 0.0_dp, u_exact
      else if (j == N + 1) then
        ! Top boundary: u = sin(pi*x)
        write(funit, datafmt) xi, yj, sin(pi * xi), u_exact
      else
        ! Interior point: solution from solver
        k = (j - 1) * N + i
        write(funit, datafmt) xi, yj, b(k), u_exact
      end if
    end do
    write(funit, *)   ! blank line separates scanlines for gnuplot pm3d
  end do
  close(funit)

  write(output_unit, '(a)') 'Solution written to poisson_9pt.dat'
  write(output_unit, '(a)') 'Run "gnuplot plot_poisson.gp" to produce poisson_9pt.png'

end program poisson_9pt
