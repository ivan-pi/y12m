! SPDX-License-Identifier: GPL-2.0-only
!
! poisson_3d.f90 - Steady-state diffusion on a unit cube
!
! Solves the 3-D Poisson equation
!
!   -nabla^2 u = f(x,y,z)   on  (0,1) x (0,1) x (0,1)
!
! with homogeneous Dirichlet boundary conditions  u = 0  on the boundary.
!
! Manufactured solution (parametrised by sine-mode triplet kx, ky, kz):
!
!   u_exact(x,y,z) = sin(kx pi x) sin(ky pi y) sin(kz pi z)
!
! The manufactured solution automatically satisfies the homogeneous Dirichlet
! BCs because sin vanishes at 0 and at every positive integer multiple of pi.
!
! Right-hand side (from -nabla^2 u_exact):
!
!   f(x,y,z) = (kx^2 + ky^2 + kz^2) pi^2 u_exact(x,y,z)
!
! Derivation: apply -nabla^2 to u_exact term by term:
!   -u_xx = kx^2 pi^2 sin(kx pi x) sin(ky pi y) sin(kz pi z)
!   -u_yy = ky^2 pi^2 sin(kx pi x) sin(ky pi y) sin(kz pi z)
!   -u_zz = kz^2 pi^2 sin(kx pi x) sin(ky pi y) sin(kz pi z)
!
! Discretisation
! --------------
! Total grid N x N x N (including boundary nodes), spacing h = 1/(N-1).
! Interior nodes: i=2..N-1, j=2..N-1, k=2..N-1.
! ndof = (N-2)^3.
!
! DOF index for interior node (i,j,k):
!   idx = (k-2)*ni*nj + (j-2)*ni + (i-2) + 1
!   where ni = nj = nk = N-2.
!
! i varies fastest (x-direction), then j (y), then k (z), so that:
!   x-neighbour (i+/-1) -> idx +/- 1
!   y-neighbour (j+/-1) -> idx +/- ni
!   z-neighbour (k+/-1) -> idx +/- ni*nj
!
! 7-point finite difference stencil for -nabla^2 u (unscaled form):
!   Diagonal (centre):    6
!   6 face neighbours:   -1
!   Right-hand side:      h^2 * f(x,y,z)
!
! Memory
! ------
! The LU fill-in bandwidth for a 3-D direct solve is b = (N-2)^2.
! Workspace is allocated as
!   nn = nn1 = max(3*nz_orig, 2*(b+2)*ndof)
! which is adequate for N up to roughly 20 (b=324, ndof=5832, nn~3.8 M).
! For larger problems consider reordering (e.g. nested dissection) or
! switching to an iterative/multilevel solver.
!
! Solver: y12ma (double precision), no iterative refinement.
!
! Usage
! -----
!   poisson_3d [--help] [N] [kx] [ky] [kz]
!
!     N           total grid size per direction (>= 3, default 12
!                 -> (N-2)^3 = 10^3 = 1000 interior DOFs)
!     kx,ky,kz   sine-mode indices (>= 1, default: 1 1 1)
!

! ==============================================================
! Module poisson_3d_solver
! ==============================================================
module poisson_3d_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  !> Manufactured exact solution.
  elemental function u_exact(x, y, z, kx, ky, kz) result(u)
    real(dp), intent(in) :: x, y, z
    integer, intent(in) :: kx, ky, kz
    real(dp) :: u
    u = sin(real(kx, dp) * pi * x) &
      * sin(real(ky, dp) * pi * y) &
      * sin(real(kz, dp) * pi * z)
  end function u_exact

  !> Right-hand side f = -nabla^2 u_exact.
  elemental function f_rhs(x, y, z, kx, ky, kz) result(f)
    real(dp), intent(in) :: x, y, z
    integer, intent(in) :: kx, ky, kz
    real(dp) :: f
    real(dp) :: lambda
    lambda = real(kx*kx + ky*ky + kz*kz, dp) * pi**2
    f = lambda * u_exact(x, y, z, kx, ky, kz)
  end function f_rhs

  ! ---------------------------------------------------------------
  ! Sparse matrix-vector product: y <- y + alpha * A * x  (COO format)
  ! ---------------------------------------------------------------
  subroutine coo_gemv(n, nz, alpha, a, ia, ja, x, y)
    integer, intent(in) :: n, nz
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: a(nz)
    integer, intent(in) :: ia(nz), ja(nz)
    real(dp), intent(in) :: x(n)
    real(dp), intent(inout) :: y(n)
    integer :: p
    do p = 1, nz
      y(ia(p)) = y(ia(p)) + alpha * a(p) * x(ja(p))
    end do
  end subroutine coo_gemv

  ! ---------------------------------------------------------------
  ! Assemble the stiffness matrix A in COO format.
  !
  ! 7-point stencil for -nabla^2 u on a uniform N x N x N grid.
  ! Unscaled form: diagonal = 6, off-diagonal = -1, RHS = h^2 * f.
  !
  ! Only interior nodes are included; boundary contributions are zero
  ! because u = 0 on the boundary (homogeneous Dirichlet).
  ! ---------------------------------------------------------------
  subroutine assemble_matrix(N, kx, ky, kz, a, rnr, snr, b, nz)
    integer, intent(in) :: N, kx, ky, kz
    real(dp), intent(out) :: a(:), b(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    integer :: i, j, k, ni, nj, idx
    real(dp) :: h, h2, x, y_c, z_c

    ni = N - 2
    nj = N - 2
    h  = 1.0_dp / real(N - 1, dp)
    h2 = h * h

    nz = 0

    do k = 2, N - 1
      z_c = real(k - 1, dp) * h
      do j = 2, N - 1
        y_c = real(j - 1, dp) * h
        do i = 2, N - 1
          x   = real(i - 1, dp) * h
          idx = (k - 2) * ni * nj + (j - 2) * ni + (i - 2) + 1

          ! Diagonal
          nz = nz + 1
          rnr(nz) = idx; snr(nz) = idx; a(nz) = 6.0_dp

          ! x-direction neighbours
          if (i > 2) then
            nz = nz + 1
            rnr(nz) = idx; snr(nz) = idx - 1; a(nz) = -1.0_dp
          end if
          if (i < N - 1) then
            nz = nz + 1
            rnr(nz) = idx; snr(nz) = idx + 1; a(nz) = -1.0_dp
          end if

          ! y-direction neighbours
          if (j > 2) then
            nz = nz + 1
            rnr(nz) = idx; snr(nz) = idx - ni; a(nz) = -1.0_dp
          end if
          if (j < N - 1) then
            nz = nz + 1
            rnr(nz) = idx; snr(nz) = idx + ni; a(nz) = -1.0_dp
          end if

          ! z-direction neighbours
          if (k > 2) then
            nz = nz + 1
            rnr(nz) = idx; snr(nz) = idx - ni * nj; a(nz) = -1.0_dp
          end if
          if (k < N - 1) then
            nz = nz + 1
            rnr(nz) = idx; snr(nz) = idx + ni * nj; a(nz) = -1.0_dp
          end if

          ! Right-hand side: unscaled form -> b = h^2 * f
          b(idx) = h2 * f_rhs(x, y_c, z_c, kx, ky, kz)

        end do
      end do
    end do
  end subroutine assemble_matrix

  ! ---------------------------------------------------------------
  ! Assemble, solve and report the 3-D Poisson problem.
  ! ---------------------------------------------------------------
  subroutine run(N, kx, ky, kz)
    use y12m, only: y12ma
    integer, intent(in) :: N, kx, ky, kz

    integer :: ni, nj, nk, ndof, nz_orig, nz, nn, nn1, iha, bandwidth
    real(dp) :: h, h2
    real(dp), allocatable :: a(:), b(:), pivot(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail

    ! Saved originals for residual check (y12ma overwrites a, snr, rnr, b)
    real(dp), allocatable :: a0(:), b0(:)
    integer, allocatable :: snr0(:), rnr0(:)
    real(dp), allocatable :: residual(:)

    integer :: i, j, k, idx
    real(dp) :: x, y_c, z_c, u_num, u_ex_val
    real(dp) :: err_max, err_rms, res_inf, res_2
    integer :: t0, t1, clock_rate
    integer :: t_assemble, t_solve, t_error

    ni = N - 2
    nj = N - 2
    nk = N - 2
    ndof = ni * nj * nk
    h  = 1.0_dp / real(N - 1, dp)
    h2 = h * h

    write(output_unit, '(a,i0,a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
    write(output_unit, '(a,i0,a,i0,a,i0)') &
        'Mode: kx=', kx, '  ky=', ky, '  kz=', kz

    ! Workspace allocation:
    ! LU fill-in bandwidth for this ordering = ni * nj = (N-2)^2.
    ! Formula mirrors the 2-D examples: max(3*nz_orig, 2*(bandwidth+2)*ndof).
    nz_orig   = 7 * ndof
    bandwidth = ni * nj
    nn        = max(3 * nz_orig, 2 * (bandwidth + 2) * ndof)
    nn1       = nn
    iha       = ndof

    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
    allocate(a0(nz_orig), b0(ndof), snr0(nz_orig), rnr0(nz_orig), residual(ndof))

    aflag = 0.0_dp
    iflag = 0
    ifail = 0

    call system_clock(count_rate=clock_rate)

    ! ==================================================================
    ! 1. Assemble A u = b
    ! ==================================================================
    call system_clock(t0)
    b = 0.0_dp
    call assemble_matrix(N, kx, ky, kz, a, rnr, snr, b, nz)

    ! Save originals before y12ma overwrites a, snr, rnr, b
    a0(1:nz)  = a(1:nz)
    b0        = b
    snr0(1:nz) = snr(1:nz)
    rnr0(1:nz) = rnr(1:nz)
    call system_clock(t1)
    t_assemble = t1 - t0

    ! ==================================================================
    ! 2. Solve with y12ma (double precision, no iterative refinement)
    ! ==================================================================
    call system_clock(t0)
    call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1
    end if
    call system_clock(t1)
    t_solve = t1 - t0
    ! Solution u_h is now in b(1:ndof).

    ! ==================================================================
    ! 3. Error vs exact solution and discrete residual
    ! ==================================================================
    call system_clock(t0)

    err_max = 0.0_dp
    err_rms = 0.0_dp
    do k = 2, N - 1
      z_c = real(k - 1, dp) * h
      do j = 2, N - 1
        y_c = real(j - 1, dp) * h
        do i = 2, N - 1
          x   = real(i - 1, dp) * h
          idx = (k - 2) * ni * nj + (j - 2) * ni + (i - 2) + 1
          u_num    = b(idx)
          u_ex_val = u_exact(x, y_c, z_c, kx, ky, kz)
          err_max  = max(err_max, abs(u_num - u_ex_val))
          err_rms  = err_rms + (u_num - u_ex_val)**2
        end do
      end do
    end do
    err_rms = sqrt(err_rms / real(ndof, dp))

    ! Discrete residual  r = b0 - A0 * u_h
    residual = b0
    call coo_gemv(ndof, nz, -1.0_dp, a0, rnr0, snr0, b, residual)
    res_inf = maxval(abs(residual))
    res_2   = norm2(residual)

    call system_clock(t1)
    t_error = t1 - t0

    write(output_unit, '(a,es12.4)') 'Max pointwise error  : ', err_max
    write(output_unit, '(a,es12.4)') 'RMS pointwise error  : ', err_rms
    write(output_unit, '(a,es12.4)') 'Residual ||r||_inf   : ', res_inf
    write(output_unit, '(a,es12.4)') 'Residual ||r||_2     : ', res_2

    ! ==================================================================
    ! Solver diagnostics (values set by y12ma)
    ! ==================================================================
    write(output_unit, '(a)') ''
    write(output_unit, '(a)') 'Solver diagnostics (set by y12ma):'
    write(output_unit, '(2x,a,es14.6)') 'Stability factor   AFLAG(1) = ', aflag(1)
    write(output_unit, '(2x,a,es14.6)') 'Drop tolerance     AFLAG(2) = ', aflag(2)
    write(output_unit, '(2x,a,es14.6)') 'Growth limit       AFLAG(3) = ', aflag(3)
    write(output_unit, '(2x,a,es14.6)') 'Pivot tolerance    AFLAG(4) = ', aflag(4)
    write(output_unit, '(2x,a,es14.6)') 'Growth factor      AFLAG(5) = ', aflag(5)
    write(output_unit, '(2x,a,es14.6)') 'Max element in A   AFLAG(6) = ', aflag(6)
    write(output_unit, '(2x,a,es14.6)') 'Max element in LU  AFLAG(7) = ', aflag(7)
    write(output_unit, '(2x,a,es14.6)') 'Min pivot          AFLAG(8) = ', aflag(8)
    write(output_unit, '(2x,a,i0)')     'Row garbage coll.  IFLAG(6) = ', iflag(6)
    write(output_unit, '(2x,a,i0)')     'Col garbage coll.  IFLAG(7) = ', iflag(7)
    write(output_unit, '(2x,a,i0)')     'Max nnz in A       IFLAG(8) = ', iflag(8)

    ! ==================================================================
    ! Timing summary
    ! ==================================================================
    write(output_unit, '(a)') ''
    write(output_unit, '(a)') 'Timing summary:'
    write(output_unit, '(a,f10.6,a)') &
        '  Assembly     : ', real(t_assemble, dp) / real(clock_rate, dp), ' s'
    write(output_unit, '(a,f10.6,a)') &
        '  Solve        : ', real(t_solve, dp) / real(clock_rate, dp), ' s'
    write(output_unit, '(a,f10.6,a)') &
        '  Error+resid  : ', real(t_error, dp) / real(clock_rate, dp), ' s'

  end subroutine run

end module poisson_3d_solver

! ==============================================================
! Main program: parse CLI, call run(N, kx, ky, kz).
! ==============================================================
program poisson_3d
  use poisson_3d_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N, kx, ky, kz, ios, iarg, nargs
  character(len=256) :: arg
  integer, parameter :: default_n  = 12
  integer, parameter :: default_kx = 1
  integer, parameter :: default_ky = 1
  integer, parameter :: default_kz = 1

  N  = default_n
  kx = default_kx
  ky = default_ky
  kz = default_kz

  nargs = command_argument_count()

  ! Check for --help flag anywhere in the argument list
  do iarg = 1, nargs
    call get_command_argument(iarg, arg)
    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: poisson_3d [--help] [N] [kx] [ky] [kz]'
      write(*, '(a)') '  N           total grid size per direction (>= 3, default 12)'
      write(*, '(a)') '              (N-2)^3 interior DOFs; e.g. N=12 -> 10^3 = 1000 DOFs'
      write(*, '(a)') '  kx,ky,kz   sine-mode indices (>= 1, default: 1 1 1)'
      stop 0
    end if
  end do

  ! Parse positional argument 1: N
  if (nargs >= 1) then
    call get_command_argument(1, arg)
    read(arg, *, iostat=ios) N
    if (ios /= 0 .or. N < 3) then
      write(error_unit, '(a)') 'ERROR: N must be an integer >= 3'
      stop 1
    end if
  end if

  ! Parse positional argument 2: kx
  if (nargs >= 2) then
    call get_command_argument(2, arg)
    read(arg, *, iostat=ios) kx
    if (ios /= 0 .or. kx < 1) then
      write(error_unit, '(a)') 'ERROR: kx must be an integer >= 1'
      stop 1
    end if
  end if

  ! Parse positional argument 3: ky
  if (nargs >= 3) then
    call get_command_argument(3, arg)
    read(arg, *, iostat=ios) ky
    if (ios /= 0 .or. ky < 1) then
      write(error_unit, '(a)') 'ERROR: ky must be an integer >= 1'
      stop 1
    end if
  end if

  ! Parse positional argument 4: kz
  if (nargs >= 4) then
    call get_command_argument(4, arg)
    read(arg, *, iostat=ios) kz
    if (ios /= 0 .or. kz < 1) then
      write(error_unit, '(a)') 'ERROR: kz must be an integer >= 1'
      stop 1
    end if
  end if

  call run(N, kx, ky, kz)

end program poisson_3d
