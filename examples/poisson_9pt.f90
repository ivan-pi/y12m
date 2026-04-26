! SPDX-License-Identifier: GPL-2.0-only
!
! poisson_9pt.f90 - Steady-state heat diffusion on a unit square
!
! Solves the 2-D Laplace equation
!
!   Δu = 0   on  Ω = (0,1) × (0,1)
!
! with boundary conditions (a plate with a parabolic hot top edge):
!
!   u(x, 0) = 0                bottom  (cold)
!   u(x, 1) = 4 x (1 - x)     top     (parabolic hot edge, peak = 1)
!   u(0, y) = 0                left    (cold)
!   u(1, y) = 0                right   (cold)
!
! Exact solution (Fourier sine series, odd modes only):
!
!   u(x,y) = Σ_{n=1,3,5,...} [32 / (n π)³]  sin(n π x)  sinh(n π y) / sinh(n π)
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
! Boundary values are moved to the right-hand side.
!
! Solver
! ------
! y12ma (resolving to y12maf, double precision) performs Gaussian
! elimination without iterative refinement.
!
! Usage
! -----
!   poisson_9pt [N]           (integer N ≥ 2, default 20)
!   gnuplot plot_poisson.gp   (produces poisson_9pt.png)
!
program poisson_9pt
  use, intrinsic :: iso_fortran_env, only: dp => real64, &
      output_unit, error_unit
  use y12m, only: y12ma
  implicit none

  ! ---------------------------------------------------------------
  ! Grid parameters (N supplied at run time)
  ! ---------------------------------------------------------------
  integer :: N            ! interior grid points per direction
  integer :: ndof         ! total DOFs = N*N
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp) :: h           ! uniform grid spacing = 1/(N+1)

  ! ---------------------------------------------------------------
  ! Sparse-matrix workspace (dynamically allocated)
  !
  ! y12m requires nn >= 2*nz and nn1 >= nz (checked by y12mbe).
  ! During LU factorisation the semi-bandwidth grows to N+1, so
  ! both the row list (a, snr) and the column list (rnr) need room
  ! for fill-in:  nn = nn1 = max(3*nz, 2*(N+2)*ndof).
  ! ---------------------------------------------------------------
  integer :: nz_max           ! upper bound: 9*ndof for 9-point stencil
  integer :: nz               ! actual nnz after assembly
  integer :: nn, nn1, iha
  real(dp), allocatable :: a(:), pivot(:), b(:)
  integer, allocatable :: snr(:), rnr(:), ha(:,:)
  real(dp) :: aflag(8)
  integer :: iflag(10), ifail

  ! Saved originals for discrete residual (y12ma overwrites a, snr, rnr, b)
  real(dp), allocatable :: a0(:), b0(:), residual(:)
  integer, allocatable :: snr0(:), rnr0(:)

  ! ---------------------------------------------------------------
  ! Local scalars
  ! ---------------------------------------------------------------
  character(len=*), parameter :: datafmt = '(4(1x,es14.6))'
  integer :: i, j, k, p
  real(dp) :: xi, yj, err_max, err_l2, res_inf, res_2sq

  ! ==============================================================
  ! 1.  Read N from command line (default 20)
  ! ==============================================================
  call parse_n(N)
  ndof = N * N
  h = 1.0_dp / real(N + 1, dp)

  ! ==============================================================
  ! 2.  Allocate storage
  ! ==============================================================
  nz_max = 9 * ndof
  nn = max(3 * nz_max, 2 * (N + 2) * ndof)
  nn1 = nn
  iha = ndof
  allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
  allocate(a0(nz_max), b0(ndof), snr0(nz_max), rnr0(nz_max), residual(ndof))

  ! ==============================================================
  ! 3.  Assemble sparse system  A u = b  (row by row)
  ! ==============================================================
  b = 0.0_dp
  nz = 0
  do j = 1, N
    do i = 1, N
      k = (j - 1) * N + i
      xi = real(i, dp) * h

      ! Centre: coefficient -20
      nz = nz + 1
      rnr(nz) = k
      snr(nz) = k
      a(nz) = -20.0_dp

      ! West (i-1, j): +4, or left boundary (u = 0)
      if (i > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - 1
        a(nz) = 4.0_dp
      end if

      ! East (i+1, j): +4, or right boundary (u = 0)
      if (i < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + 1
        a(nz) = 4.0_dp
      end if

      ! South (i, j-1): +4, or bottom boundary (u = 0)
      if (j > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - N
        a(nz) = 4.0_dp
      end if

      ! North (i, j+1): +4, or top boundary u = g(xi)
      if (j < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + N
        a(nz) = 4.0_dp
      else
        b(k) = b(k) - 4.0_dp * bc_top(xi)
      end if

      ! SW (i-1, j-1): +1, or zero boundary
      if (i > 1 .and. j > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - N - 1
        a(nz) = 1.0_dp
      end if

      ! SE (i+1, j-1): +1, or zero boundary
      if (i < N .and. j > 1) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k - N + 1
        a(nz) = 1.0_dp
      end if

      ! NW (i-1, j+1): +1, or top boundary when j == N
      if (i > 1 .and. j < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + N - 1
        a(nz) = 1.0_dp
      else if (i > 1 .and. j == N) then
        b(k) = b(k) - bc_top(xi - h)
      end if

      ! NE (i+1, j+1): +1, or top boundary when j == N
      if (i < N .and. j < N) then
        nz = nz + 1
        rnr(nz) = k
        snr(nz) = k + N + 1
        a(nz) = 1.0_dp
      else if (i < N .and. j == N) then
        b(k) = b(k) - bc_top(xi + h)
      end if

    end do
  end do

  ! Save originals for residual check
  a0(1:nz) = a(1:nz)
  b0 = b
  snr0(1:nz) = snr(1:nz)
  rnr0(1:nz) = rnr(1:nz)

  ! ==============================================================
  ! 4.  Solve with y12ma (double-precision, no iterative refinement)
  ! ==============================================================
  call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)

  if (ifail /= 0) then
    write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
    stop 1
  end if
  ! Solution vector u_h is now in b(1:ndof).

  ! ==============================================================
  ! 5.  Diagnostics: pointwise error vs exact + discrete residual
  ! ==============================================================
  err_max = 0.0_dp
  err_l2 = 0.0_dp
  do j = 1, N
    do i = 1, N
      k  = (j - 1) * N + i
      xi = real(i, dp) * h
      yj = real(j, dp) * h
      err_max = max(err_max, abs(b(k) - u_exact(xi, yj)))
      err_l2  = err_l2 + (b(k) - u_exact(xi, yj))**2
    end do
  end do
  err_l2 = sqrt(err_l2 / real(ndof, dp))

  ! Discrete residual  r = b_orig - A_orig * u_h
  residual = b0
  do p = 1, nz
    residual(rnr0(p)) = residual(rnr0(p)) - a0(p) * b(snr0(p))
  end do
  res_inf  = maxval(abs(residual))
  res_2sq  = sum(residual**2)

  write(output_unit, '(a,i0,a,i0,a,i0,a)') &
      'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
  write(output_unit, '(a,es12.4)') 'Max pointwise error  : ', err_max
  write(output_unit, '(a,es12.4)') 'RMS pointwise error  : ', err_l2
  write(output_unit, '(a,es12.4)') 'Residual ||r||_inf   : ', res_inf
  write(output_unit, '(a,es12.4)') 'Residual ||r||_2^2   : ', res_2sq
  write(output_unit, '(a)') &
      'Exact: u = sum_{n odd} [32/(n pi)^3] sin(n pi x) sinh(n pi y)/sinh(n pi)'

  ! ==============================================================
  ! 6.  Write gnuplot data file
  ! ==============================================================
  call write_output(N, ndof, h, b)

  deallocate(a, pivot, b, snr, rnr, ha, a0, b0, snr0, rnr0, residual)

contains

  !> Boundary condition on the top edge y = 1: parabola 4x(1-x)
  elemental function bc_top(x) result(u)
    real(dp), intent(in) :: x
    real(dp) :: u
    u = 4.0_dp * x * (1.0_dp - x)
  end function bc_top

  !> Exact solution via Fourier sine series (odd modes only)
  !>
  !> u(x,y) = sum_{n odd} [32/(n pi)^3] sin(n pi x) sinh(n pi y)/sinh(n pi)
  function u_exact(x, y) result(u)
    real(dp), intent(in) :: x, y
    real(dp) :: u
    integer :: m
    real(dp) :: rn
    u = 0.0_dp
    do m = 1, 99, 2
      rn = real(m, dp)
      u = u + 32.0_dp / (rn * pi)**3 &
            * sin(rn * pi * x) * sinh(rn * pi * y) / sinh(rn * pi)
    end do
  end function u_exact

  !> Write solution to poisson_9pt.dat in gnuplot pm3d format.
  !>
  !> One scanline per y value; blank lines separate scanlines.
  !> Columns: x  y  u_numerical  u_exact
  subroutine write_output(n_grid, n_dof, dx, u)
    integer, intent(in) :: n_grid, n_dof
    real(dp), intent(in) :: dx, u(n_dof)
    integer :: funit, ii, jj, kk
    real(dp) :: xi2, yj2, u_num
    open(newunit=funit, file='poisson_9pt.dat', status='replace', action='write')
    write(funit, '(a)') &
        '# 9-point isotropic stencil, Laplace equation on unit square'
    write(funit, '(a)') &
        '# BC: u = 4*x*(1-x) on top edge (y=1); u = 0 elsewhere'
    write(funit, '(a)') &
        '# Columns: x   y   u_numerical   u_exact'
    write(funit, '(a)') '#'
    do jj = 0, n_grid + 1
      yj2 = real(jj, dp) * dx
      do ii = 0, n_grid + 1
        xi2 = real(ii, dp) * dx
        if (ii == 0 .or. ii == n_grid + 1 .or. jj == 0) then
          u_num = 0.0_dp
        else if (jj == n_grid + 1) then
          u_num = bc_top(xi2)
        else
          kk = (jj - 1) * n_grid + ii
          u_num = u(kk)
        end if
        write(funit, datafmt) xi2, yj2, u_num, u_exact(xi2, yj2)
      end do
      write(funit, *)
    end do
    close(funit)
    write(output_unit, '(a)') 'Solution written to poisson_9pt.dat'
    write(output_unit, '(a)') 'Run "gnuplot plot_poisson.gp" to produce poisson_9pt.png'
  end subroutine write_output

  !> Parse the grid size N from the first command-line argument.
  !> Falls back to DEFAULT_N = 20 if no argument is given.
  subroutine parse_n(n_out)
    integer, intent(out) :: n_out
    integer, parameter :: default_n = 20
    character(len=32) :: arg
    integer :: ios
    if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      read(arg, *, iostat=ios) n_out
      if (ios /= 0 .or. n_out < 2) then
        write(error_unit, '(a)') &
            'Usage: poisson_9pt [N]   (integer N >= 2, default 20)'
        stop 1
      end if
    else
      n_out = default_n
    end if
  end subroutine parse_n

end program poisson_9pt
