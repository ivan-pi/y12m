! SPDX-License-Identifier: GPL-2.0-only
!
! poisson_9pt.f90 - Steady-state heat diffusion on a unit square
!
! Solves the 2-D Laplace equation
!
!   Du = 0   on  (0,1) x (0,1)
!
! with boundary conditions (parabolic hot top edge):
!
!   u(x, 0) = 0                bottom  (cold)
!   u(x, 1) = 4 x (1 - x)     top     (parabolic, peak = 1 at x = 0.5)
!   u(0, y) = 0                left    (cold)
!   u(1, y) = 0                right   (cold)
!
! Exact solution (Fourier sine series, odd modes only):
!
!   u(x,y) = sum_{n=1,3,5,...} [32 / (n pi)^3]
!               x sin(n pi x) sinh(n pi y) / sinh(n pi)
!
! Discretisation
! --------------
! N x N interior grid, spacing h = 1/(N+1).
! DOF index k = (j-1)*N + i  for  i,j = 1..N (i: x, j: y).
!
! 9-point isotropic Laplacian (Mehrstellen) stencil:
!
!   [1  4  1]         1
!   [4 -20  4]  x  -------
!   [1  4  1]        6 h^2
!
! Boundary values are moved to the right-hand side.
!
! Solver: y12ma (double precision), no iterative refinement.
!
! Usage
! -----
!   poisson_9pt [N]           (integer N >= 2, default 20)
!   gnuplot plot_poisson.gp   (produces poisson_9pt.png)
!

! ==============================================================
! Module poisson_solver
!
! Provides 'run(N)': assembles, solves and reports the discrete
! Laplace problem.  The main program is intentionally trivial.
! ==============================================================
module poisson_solver
  use, intrinsic :: iso_fortran_env, only: dp => real64, &
      output_unit, error_unit
  use y12m, only: y12ma
  implicit none
  private
  public :: run

  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  !> Boundary condition on the top edge y = 1: parabola 4x(1-x)
  elemental function bc_top(x) result(u)
    real(dp), intent(in) :: x
    real(dp) :: u
    u = 4.0_dp * x * (1.0_dp - x)
  end function bc_top

  !> Exact solution via Fourier sine series (50 odd modes)
  !>
  !>   u = sum_{n odd} [32/(n pi)^3] sin(n pi x) sinh(n pi y) / sinh(n pi)
  elemental function u_exact(x, y) result(u)
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

  !> Write solution to poisson_9pt.dat in gnuplot pm3d block format.
  !>
  !> u_num: numerical solution at interior points, shape (n_grid, n_grid).
  !> u_ex:  exact solution at interior points, precomputed by the caller.
  !>
  !> One blank-line-separated scanline per y-value (boundaries included).
  !> Columns: x  y  u_numerical  u_exact
  subroutine write_output(n_grid, u_num, u_ex, dx)
    integer, intent(in) :: n_grid
    real(dp), intent(in) :: u_num(n_grid, n_grid)
    real(dp), intent(in) :: u_ex(n_grid, n_grid)
    real(dp), intent(in) :: dx
    character(len=*), parameter :: datafmt = '(4(1x,es14.6))'
    integer :: funit, ii, jj
    real(dp) :: xi2, yj2, u_n, u_e
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
        if (jj == 0 .or. ii == 0 .or. ii == n_grid + 1) then
          u_n = 0.0_dp
          u_e = u_exact(xi2, yj2)
        else if (jj == n_grid + 1) then
          u_n = bc_top(xi2)
          u_e = u_exact(xi2, yj2)
        else
          u_n = u_num(ii, jj)
          u_e = u_ex(ii, jj)
        end if
        write(funit, datafmt) xi2, yj2, u_n, u_e
      end do
      write(funit, *)
    end do
    close(funit)
  end subroutine write_output

  !> Assemble and solve the discrete Laplace problem, compute error
  !> diagnostics, and write the solution to poisson_9pt.dat.
  subroutine run(N)
    integer, intent(in) :: N

    ! -----------------------------------------------------------
    ! Workspace (allocatable; freed automatically on return)
    ! -----------------------------------------------------------
    integer :: ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h
    real(dp), allocatable, target :: b(:)
    real(dp), allocatable :: a(:), pivot(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail

    ! Saved originals for residual check (y12ma overwrites a, snr, rnr, b)
    real(dp), allocatable :: a0(:), b0(:), residual(:)
    integer, allocatable :: snr0(:), rnr0(:)

    ! Exact solution at interior points (saved once, reused in output)
    real(dp), allocatable :: u_ex(:,:)

    integer :: i, j, k, p
    real(dp) :: xi, err_max, err_rms, res_inf, res_2

    ndof = N * N
    h    = 1.0_dp / real(N + 1, dp)

    ! -----------------------------------------------------------
    ! Allocate workspace
    ! -----------------------------------------------------------
    nz_max = 9 * ndof
    nn     = max(3 * nz_max, 2 * (N + 2) * ndof)
    nn1    = nn
    iha    = ndof
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
    allocate(a0(nz_max), b0(ndof), snr0(nz_max), rnr0(nz_max), residual(ndof))
    allocate(u_ex(N, N))

    ! ==============================================================
    ! 1. Assemble sparse system  A u = b
    !
    ! Outer loop over x-index i, inner over y-index j.
    ! An if-else-if block handles the four boundary column/row cases;
    ! the else branch covers fully interior points.
    ! b is zero-initialised; only non-zero BC contributions are added.
    ! ==============================================================
    b  = 0.0_dp
    nz = 0
    do i = 1, N
      do j = 1, N
        k  = (j - 1) * N + i
        xi = real(i, dp) * h

        ! Centre: always -20
        nz      = nz + 1
        rnr(nz) = k
        snr(nz) = k
        a(nz)   = -20.0_dp

        if (i == 1) then
          ! ---- Left column ----------------------------------------
          ! No West / SW / NW entries (left boundary, u = 0)

          ! East
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + 1
          a(nz)   = 4.0_dp

          ! South (or bottom boundary, u = 0)
          if (j > 1) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k - N
            a(nz)   = 4.0_dp
          end if

          ! North (or top boundary)
          if (j < N) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k + N
            a(nz)   = 4.0_dp
          else
            b(k) = b(k) - 4.0_dp * bc_top(xi)
          end if

          ! SE (or bottom boundary, u = 0)
          if (j > 1) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k - N + 1
            a(nz)   = 1.0_dp
          end if

          ! NE (or top boundary)
          if (j < N) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k + N + 1
            a(nz)   = 1.0_dp
          else
            b(k) = b(k) - bc_top(xi + h)
          end if

        else if (i == N) then
          ! ---- Right column ---------------------------------------
          ! No East / NE / SE entries (right boundary, u = 0)

          ! West
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - 1
          a(nz)   = 4.0_dp

          ! South (or bottom boundary, u = 0)
          if (j > 1) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k - N
            a(nz)   = 4.0_dp
          end if

          ! North (or top boundary)
          if (j < N) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k + N
            a(nz)   = 4.0_dp
          else
            b(k) = b(k) - 4.0_dp * bc_top(xi)
          end if

          ! SW (or bottom boundary, u = 0)
          if (j > 1) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k - N - 1
            a(nz)   = 1.0_dp
          end if

          ! NW (or top boundary)
          if (j < N) then
            nz      = nz + 1
            rnr(nz) = k
            snr(nz) = k + N - 1
            a(nz)   = 1.0_dp
          else
            b(k) = b(k) - bc_top(xi - h)
          end if

        else if (j == 1) then
          ! ---- Bottom row (interior columns only) -----------------
          ! No South / SW / SE entries (bottom boundary, u = 0)

          ! West
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - 1
          a(nz)   = 4.0_dp

          ! East
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + 1
          a(nz)   = 4.0_dp

          ! North (j+1 = 2; this branch is only reachable when i > 1 and i < N, so N >= 3)
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + N
          a(nz)   = 4.0_dp

          ! NW
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + N - 1
          a(nz)   = 1.0_dp

          ! NE
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + N + 1
          a(nz)   = 1.0_dp

        else if (j == N) then
          ! ---- Top row (interior columns only) --------------------
          ! North / NW / NE lie on top boundary (u = bc_top)

          ! West
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - 1
          a(nz)   = 4.0_dp

          ! East
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + 1
          a(nz)   = 4.0_dp

          ! South (j-1 = N-1 >= 1)
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - N
          a(nz)   = 4.0_dp

          ! North: top boundary
          b(k) = b(k) - 4.0_dp * bc_top(xi)

          ! SW
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - N - 1
          a(nz)   = 1.0_dp

          ! SE
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - N + 1
          a(nz)   = 1.0_dp

          ! NW: top boundary
          b(k) = b(k) - bc_top(xi - h)

          ! NE: top boundary
          b(k) = b(k) - bc_top(xi + h)

        else
          ! ---- Interior point: all 8 neighbours are interior DOFs --

          ! West
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - 1
          a(nz)   = 4.0_dp

          ! East
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + 1
          a(nz)   = 4.0_dp

          ! South
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - N
          a(nz)   = 4.0_dp

          ! North
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + N
          a(nz)   = 4.0_dp

          ! SW
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - N - 1
          a(nz)   = 1.0_dp

          ! SE
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k - N + 1
          a(nz)   = 1.0_dp

          ! NW
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + N - 1
          a(nz)   = 1.0_dp

          ! NE
          nz      = nz + 1
          rnr(nz) = k
          snr(nz) = k + N + 1
          a(nz)   = 1.0_dp

        end if

      end do
    end do

    ! Save originals before y12ma overwrites a, snr, rnr, b
    a0(1:nz)   = a(1:nz)
    b0         = b
    snr0(1:nz) = snr(1:nz)
    rnr0(1:nz) = rnr(1:nz)

    ! ==============================================================
    ! 2. Solve with y12ma (double precision, no iterative refinement)
    ! ==============================================================
    ifail = 0
    call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)

    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if
    ! Solution u_h is now in b(1:ndof).

    ! ==============================================================
    ! 3. Error vs exact solution; save u_ex(i,j) for output later
    !
    ! A 2-D pointer u2d remaps b(1:ndof) to u2d(1:N,1:N) so that
    ! the solution is accessed as u2d(i,j), matching the stencil.
    ! ==============================================================
    err_max = 0.0_dp
    block
      real(dp), pointer :: u2d(:,:)
      u2d(1:N, 1:N) => b
      do j = 1, N
        do i = 1, N
          u_ex(i, j) = u_exact(real(i, dp) * h, real(j, dp) * h)
          err_max = max(err_max, abs(u2d(i, j) - u_ex(i, j)))
        end do
      end do
      err_rms = norm2(u2d(1:N, 1:N) - u_ex) / sqrt(real(ndof, dp))
    end block

    ! ==============================================================
    ! 4. Discrete residual  r = b0 - A0 * u_h
    ! ==============================================================
    residual = b0
    do p = 1, nz
      residual(rnr0(p)) = residual(rnr0(p)) - a0(p) * b(snr0(p))
    end do
    res_inf = maxval(abs(residual))
    res_2   = norm2(residual)

    write(output_unit, '(a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
    write(output_unit, '(a,es12.4)') 'Max pointwise error  : ', err_max
    write(output_unit, '(a,es12.4)') 'RMS pointwise error  : ', err_rms
    write(output_unit, '(a,es12.4)') 'Residual ||r||_inf   : ', res_inf
    write(output_unit, '(a,es12.4)') 'Residual ||r||_2     : ', res_2
    write(output_unit, '(a)') &
        'Exact: u = sum_{n odd} [32/(n pi)^3] sin(n pi x) sinh(n pi y)/sinh(n pi)'

    ! ==============================================================
    ! 5. Write gnuplot data file using 2-D view of b
    ! ==============================================================
    block
      real(dp), pointer :: u2d(:,:)
      u2d(1:N, 1:N) => b
      call write_output(N, u2d, u_ex, h)
    end block

    ! Report completion (logging here, not inside write_output)
    write(output_unit, '(a)') 'Solution written to poisson_9pt.dat'
    write(output_unit, '(a)') 'Run "gnuplot plot_poisson.gp" to produce poisson_9pt.png'

  end subroutine run

end module poisson_solver

! ==============================================================
! Main program: parse N from command line, call run(N).
! ==============================================================
program poisson_9pt
  use poisson_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N

  ! Parse grid size from command line, then delegate to run(N)
  call parse_n(N)
  call run(N)

contains

  !> Read grid size N from the first command-line argument (default 20).
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
