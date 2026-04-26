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
!               * sin(n pi x) sinh(n pi y) / sinh(n pi)
!
! Discretisation
! --------------
! Total grid N x N (including boundary nodes), spacing h = 1/(N-1).
! Interior nodes: i,j = 2..N-1 (1-indexed); ndof = (N-2)^2.
! DOF index: k(i,j) = (j-2)*(N-2) + (i-2) + 1.
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
!   poisson_9pt [--help] [N] [output_file]
!
!     N           total grid size (>= 3, default 22 -> 20x20 interior DOFs)
!     output_file output data file name (default: poisson_9pt.dat)
!

! ==============================================================
! Module poisson_solver
! ==============================================================
module poisson_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: pi = acos(-1.0_dp)

  ! 9-point Mehrstellen stencil weights, indexed as stencil_w(is,js)
  ! for offsets is (x/E-W) and js (y/N-S), each ranging from -1 to +1.
  real(dp), parameter :: stencil_w(-1:1, -1:1) = reshape( &
      [1.0_dp, 4.0_dp, 1.0_dp, &
       4.0_dp, -20.0_dp, 4.0_dp, &
       1.0_dp, 4.0_dp, 1.0_dp], [3, 3])

contains

  !> Boundary condition on the top edge y = 1: parabola 4x(1-x)
  elemental function bc_top(x) result(u)
    real(dp), intent(in) :: x
    real(dp) :: u
    u = 4.0_dp * x * (1.0_dp - x)
  end function bc_top

  !> Exact solution via Fourier sine series (50 odd modes).
  !>
  !> Derivation (separation of variables):
  !>
  !>   Write u(x,y) = X(x) Y(y).  The Laplace equation gives
  !>
  !>     X'' + lambda X = 0,   X(0) = X(1) = 0
  !>       =>  X_n(x) = sin(n pi x),   lambda_n = (n pi)^2
  !>
  !>     Y'' - lambda_n Y = 0,   Y_n(0) = 0
  !>       =>  Y_n(y) = sinh(n pi y)
  !>
  !>   Hence  u(x,y) = sum_n  a_n  sin(n pi x)  sinh(n pi y).
  !>
  !>   Matching the top BC u(x,1) = g(x) = 4x(1-x) requires
  !>
  !>     sum_n  a_n sinh(n pi)  sin(n pi x)  =  g(x)
  !>
  !>   The Fourier sine coefficients of g are
  !>
  !>     g_n = 2 integral_0^1  4x(1-x) sin(n pi x) dx
  !>
  !>   For odd n:  g_n = 32 / (n pi)^3;  for even n:  g_n = 0.
  !>   Therefore  a_n = g_n / sinh(n pi), giving
  !>
  !>     u(x,y) = sum_{n odd} [32/(n pi)^3] sin(n pi x) sinh(n pi y)/sinh(n pi)
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

  ! ---------------------------------------------------------------
  ! Sparse matrix-vector product: y <- y + alpha * A * x  (COO)
  ! ---------------------------------------------------------------
  subroutine dp_coo_gemv(m, nz, alpha, a, ia, ja, x, y)
    integer, intent(in) :: m, nz
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: a(nz)
    integer, intent(in) :: ia(nz), ja(nz)
    real(dp), intent(in) :: x(m)
    real(dp), intent(inout) :: y(m)
    integer :: p
    do p = 1, nz
      y(ia(p)) = y(ia(p)) + alpha * a(p) * x(ja(p))
    end do
  end subroutine dp_coo_gemv

  ! ---------------------------------------------------------------
  ! Root-mean-square of element-wise difference of two 2-D arrays
  ! ---------------------------------------------------------------
  function rms_diff(a, b) result(rms)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: rms
    rms = norm2(a - b) / sqrt(real(size(a), dp))
  end function rms_diff

  ! ---------------------------------------------------------------
  ! Assemble the stiffness matrix A in COO format.
  !
  ! Interior nodes i,j = 2..N-1 with DOF index k(i,j) = (j-2)*(N-2)+(i-2)+1.
  ! For each interior node and each stencil offset (is,js), the entry
  ! is added to A only when the neighbour is also interior.
  ! ---------------------------------------------------------------
  subroutine assemble_matrix(N, a, rnr, snr, nz)
    integer, intent(in) :: N
    real(dp), intent(out) :: a(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz
    integer :: i, j, is, js, ni, nj, k_row, k_col, ndof_1d
    ndof_1d = N - 2
    nz = 0
    do j = 2, N - 1
      do i = 2, N - 1
        k_row = (j - 2) * ndof_1d + (i - 2) + 1
        do js = -1, 1
          do is = -1, 1
            ni = i + is
            nj = j + js
            if (ni >= 2 .and. ni <= N - 1 .and. nj >= 2 .and. nj <= N - 1) then
              k_col   = (nj - 2) * ndof_1d + (ni - 2) + 1
              nz      = nz + 1
              rnr(nz) = k_row
              snr(nz) = k_col
              a(nz)   = stencil_w(is, js)
            end if
          end do
        end do
      end do
    end do
  end subroutine assemble_matrix

  ! ---------------------------------------------------------------
  ! Compute the right-hand side from Dirichlet boundary conditions.
  !
  ! Only the top edge (j == N, u = bc_top) is non-zero.  The nodes
  ! adjacent to it are the row j = N-1; their js=+1 stencil neighbours
  ! lie on the top boundary.  bc_top is zero at the two corners
  ! (x = 0 and x = 1), so no special branch is needed for them.
  ! ---------------------------------------------------------------
  subroutine compute_rhs(N, b)
    integer, intent(in) :: N
    real(dp), intent(inout) :: b(:)
    integer :: i, is, k_row, ndof_1d
    real(dp) :: h
    h       = 1.0_dp / real(N - 1, dp)
    ndof_1d = N - 2
    do i = 2, N - 1
      k_row = (N - 3) * ndof_1d + (i - 2) + 1
      do is = -1, 1
        b(k_row) = b(k_row) - stencil_w(is, 1) * bc_top(real(i + is - 1, dp) * h)
      end do
    end do
  end subroutine compute_rhs

  ! ---------------------------------------------------------------
  ! Write the solution to a gnuplot pm3d data file.
  !
  ! u_num(N,N): numerical solution on the full N x N grid.
  ! u_ex(N,N):  exact solution on the full N x N grid.
  ! filename:   output file path.
  ! header:     nhr comment lines (each up to 70 chars) written at the top.
  !
  ! Node (i,j) is at physical coordinates (x,y) = ((i-1)*h, (j-1)*h).
  ! Columns: x  y  u_numerical  u_exact
  ! ---------------------------------------------------------------
  subroutine write_output(N, u_num, u_ex, filename, header, nhr)
    integer, intent(in) :: N
    real(dp), intent(in) :: u_num(N, N)
    real(dp), intent(in) :: u_ex(N, N)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nhr
    character(len=70), intent(in) :: header(nhr)
    character(len=*), parameter :: datafmt = '(4(1x,es14.6))'
    real(dp) :: h
    integer :: funit, i, j, ih
    h = 1.0_dp / real(N - 1, dp)
    open(newunit=funit, file=filename, status='unknown', action='write')
    do ih = 1, nhr
      write(funit, '("# ",a)') trim(header(ih))
    end do
    write(funit, '(a)') '# Columns: x   y   u_numerical   u_exact'
    write(funit, '(a)') '#'
    do j = 1, N
      do i = 1, N
        write(funit, datafmt) real(i - 1, dp) * h, real(j - 1, dp) * h, &
            u_num(i, j), u_ex(i, j)
      end do
      write(funit, *)
    end do
    close(funit)
  end subroutine write_output

  ! ---------------------------------------------------------------
  ! Assemble, solve and report the discrete Laplace problem.
  !
  ! N       total grid size (including boundaries); interior = (N-2)^2 DOFs
  ! outfile output data file name
  ! ---------------------------------------------------------------
  subroutine run(N, outfile)
    use y12m, only: y12ma
    integer, intent(in) :: N
    character(len=*), intent(in) :: outfile

    ! -----------------------------------------------------------
    ! Workspace (all allocatable; freed automatically on return)
    ! -----------------------------------------------------------
    integer :: ndof_1d, ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h
    real(dp), allocatable :: b(:), a(:), pivot(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail

    ! Saved originals for residual check (y12ma overwrites a, snr, rnr, b)
    real(dp), allocatable :: a0(:), b0(:), residual(:)
    integer, allocatable :: snr0(:), rnr0(:)

    ! Full N x N solution and exact arrays (includes boundary nodes)
    real(dp), allocatable :: u_full(:,:), u_ex_full(:,:)

    integer :: i, j, k
    real(dp) :: err_max, err_rms, res_inf, res_2

    ! Timing
    integer :: t0, t1, clock_rate
    integer :: t_assemble, t_solve, t_error, t_output

    ndof_1d = N - 2
    ndof    = ndof_1d * ndof_1d
    h       = 1.0_dp / real(N - 1, dp)

    write(output_unit, '(a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'

    ! -----------------------------------------------------------
    ! Allocate workspace
    ! -----------------------------------------------------------
    nz_max = 9 * ndof
    nn     = max(3 * nz_max, 2 * (ndof_1d + 2) * ndof)
    nn1    = nn
    iha    = ndof
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
    allocate(a0(nz_max), b0(ndof), snr0(nz_max), rnr0(nz_max), residual(ndof))
    allocate(u_full(N, N), u_ex_full(N, N))

    call system_clock(count_rate=clock_rate)

    ! ==============================================================
    ! 1. Assemble A u = b
    ! ==============================================================
    call system_clock(t0)
    b  = 0.0_dp
    call assemble_matrix(N, a, rnr, snr, nz)
    call compute_rhs(N, b)

    ! Save originals before y12ma overwrites a, snr, rnr, b
    a0(1:nz) = a(1:nz)
    b0 = b
    snr0(1:nz) = snr(1:nz)
    rnr0(1:nz) = rnr(1:nz)
    call system_clock(t1)
    t_assemble = t1 - t0

    ! ==============================================================
    ! 2. Solve with y12ma (double precision, no iterative refinement)
    ! ==============================================================
    call system_clock(t0)
    ifail = 0
    call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if
    call system_clock(t1)
    t_solve = t1 - t0
    ! Solution u_h is now in b(1:ndof).

    ! ==============================================================
    ! 3. Error vs exact solution and discrete residual
    !
    ! Map b back into u_full(2:N-1,2:N-1) and fill boundary values.
    ! ==============================================================
    call system_clock(t0)

    ! Fill full solution array (boundary + interior)
    u_full = 0.0_dp
    do j = 2, N - 1
      do i = 2, N - 1
        k = (j - 2) * ndof_1d + (i - 2) + 1
        u_full(i, j) = b(k)
      end do
    end do
    do i = 1, N
      u_full(i, N) = bc_top(real(i - 1, dp) * h)  ! top edge
    end do

    ! Exact solution on the full grid
    do j = 1, N
      do i = 1, N
        u_ex_full(i, j) = u_exact(real(i - 1, dp) * h, real(j - 1, dp) * h)
      end do
    end do

    ! Pointwise errors at interior nodes only
    err_max = maxval(abs(u_full(2:N-1, 2:N-1) - u_ex_full(2:N-1, 2:N-1)))
    err_rms = rms_diff(u_full(2:N-1, 2:N-1), u_ex_full(2:N-1, 2:N-1))

    ! Discrete residual  r = b0 - A0 * u_h
    residual = b0
    call dp_coo_gemv(ndof, nz, alpha=-1.0_dp, a=a0, ia=rnr0, ja=snr0, x=b, y=residual)
    res_inf = maxval(abs(residual))
    res_2   = norm2(residual)

    call system_clock(t1)
    t_error = t1 - t0

    write(output_unit, '(a,es12.4)') 'Max pointwise error  : ', err_max
    write(output_unit, '(a,es12.4)') 'RMS pointwise error  : ', err_rms
    write(output_unit, '(a,es12.4)') 'Residual ||r||_inf   : ', res_inf
    write(output_unit, '(a,es12.4)') 'Residual ||r||_2     : ', res_2
    write(output_unit, '(a)') &
        'Exact: u = sum_{n odd} [32/(n pi)^3] sin(n pi x) sinh(n pi y)/sinh(n pi)'

    ! ==============================================================
    ! 4. Write gnuplot data file
    ! ==============================================================
    call system_clock(t0)
    block
      character(len=70) :: hdr(2)
      hdr(1) = '9-point isotropic stencil, Laplace equation on unit square'
      hdr(2) = 'BC: u = 4*x*(1-x) on top edge (y=1); u = 0 elsewhere'
      call write_output(N, u_full, u_ex_full, outfile, hdr, size(hdr))
    end block
    call system_clock(t1)
    t_output = t1 - t0

    write(output_unit, '(a)') 'Solution written to ' // trim(outfile)
    write(output_unit, '(a)') 'Run "gnuplot plot_poisson.gp" to produce poisson_9pt.png'

    ! ==============================================================
    ! Timing summary
    ! ==============================================================
    write(output_unit, '(/,a)') 'Timing summary:'
    write(output_unit, '(a,f10.6,a)') &
        '  Assembly     : ', real(t_assemble, dp) / real(clock_rate, dp), ' s'
    write(output_unit, '(a,f10.6,a)') &
        '  Solve        : ', real(t_solve, dp) / real(clock_rate, dp), ' s'
    write(output_unit, '(a,f10.6,a)') &
        '  Error+resid  : ', real(t_error, dp) / real(clock_rate, dp), ' s'
    write(output_unit, '(a,f10.6,a)') &
        '  Output       : ', real(t_output, dp) / real(clock_rate, dp), ' s'

  end subroutine run

end module poisson_solver

! ==============================================================
! Main program: parse CLI, call run(N, outfile).
! ==============================================================
program poisson_9pt
  use poisson_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N, ios, iarg, nargs
  character(len=256) :: arg, outfile
  integer, parameter :: default_n = 22
  character(len=*), parameter :: default_outfile = 'poisson_9pt.dat'

  N      = default_n
  outfile = default_outfile
  nargs   = command_argument_count()

  ! Check for --help flag anywhere in the argument list
  do iarg = 1, nargs
    call get_command_argument(iarg, arg)
    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: poisson_9pt [--help] [N] [output_file]'
      write(*, '(a)') '  N           total grid size (>= 3, default 22)'
      write(*, '(a)') '              22 -> 20x20 interior DOFs'
      write(*, '(a)') '  output_file output data file (default: poisson_9pt.dat)'
      stop 0
    end if
  end do

  ! Parse positional argument 1: N
  if (nargs >= 1) then
    call get_command_argument(1, arg)
    read(arg, *, iostat=ios) N
    if (ios /= 0 .or. N < 3) then
      write(error_unit, '(a)') 'ERROR: N must be an integer >= 3'
      write(error_unit, '(a)') 'Usage: poisson_9pt [--help] [N] [output_file]'
      stop 1
    end if
  end if

  ! Parse positional argument 2: output file name
  if (nargs >= 2) then
    call get_command_argument(2, outfile)
  end if

  call run(N, trim(outfile))

end program poisson_9pt
