! SPDX-License-Identifier: GPL-2.0-only
!
! newton_bratu.f90 - Newton iteration for the 2-D Bratu problem
!
! Solves the nonlinear Bratu (Liouville-Gelfand) problem
!
!   -nabla^2 u = lambda exp(u)    on (0,1) x (0,1)
!
! with homogeneous Dirichlet boundary conditions  u = 0  on the boundary.
!
! The problem has two solutions for lambda < lambda_c ~= 6.808, no solution
! for lambda > lambda_c, and a unique solution for lambda = lambda_c.
!
! Newton iteration
! ----------------
! Linearising F(u) := K u - h^2 lambda exp(u) = 0  about u^k:
!
!   J(u^k) delta_u = -F(u^k)
!   u^{k+1} = u^k + delta_u
!
! where:
!   F(u)_k   = [K u]_k - h^2 lambda exp(u_k)        (residual)
!   J(u)     = K - h^2 lambda diag(exp(u))            (Jacobian)
!
! K is the 5-point Laplacian stiffness (diagonal = 4, off-diagonal = -1).
!
! Discretisation
! --------------
! Total grid N x N (including boundary nodes), spacing h = 1/(N-1).
! Interior nodes: ndof = (N-2)^2.
! DOF index: k(i,j) = (j-2)*(N-2) + (i-2) + 1.
!
! Solver strategy -- Case (iii): same sparsity, changing values
! -------------------------------------------------------------
! The Jacobian J(u^k) has the SAME non-zero positions as K at every
! Newton step (only the diagonal values change).
!
!   * First step: IFLAG(4) = 1.  y12mb computes and stores the Markowitz
!     column ordering in HA.  y12mc + y12md complete the first solve.
!   * Subsequent steps: IFLAG(4) = 2.  y12mb reuses the stored ordering
!     from HA; it skips the expensive ordering computation and only
!     inserts the new numerical values into the already-known positions.
!
! This avoids recomputing the Markowitz ordering at each Newton step,
! which is advantageous when the number of iterations is large or when
! the problem is solved repeatedly (parameter continuation).
!
! Usage
! -----
!   newton_bratu [--help] [N] [lambda] [max_iter] [output_file]
!
!     N           total grid size (>= 3, default 22)
!     lambda      reaction parameter (default 1.0, must be < 6.808)
!     max_iter    maximum Newton iterations (default 20)
!     output_file output data file name (default: newton_bratu.dat)
!

! ==============================================================
! Module bratu_solver
! ==============================================================
module bratu_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)

contains

  !> Compute the Bratu residual F(u) and its infinity norm.
  !>
  !> F_k = [K u]_k - h2lam * exp(u_k)
  !> where [K u]_k = 4 u_k - sum of neighbour values.
  subroutine compute_residual(N, h2lam, u, F, fnorm)
    integer, intent(in) :: N
    real(dp), intent(in) :: h2lam
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: F(:)
    real(dp), intent(out) :: fnorm

    integer :: ndof_1d, i, j, k
    real(dp) :: ku_k

    ndof_1d = N - 2
    fnorm = 0.0_dp

    do j = 2, N - 1
      do i = 2, N - 1
        k = (j - 2) * ndof_1d + (i - 2) + 1

        ! 5-point stencil action [K u]_k = 4 u_k - neighbours
        ku_k = 4.0_dp * u(k)
        if (i > 2)   ku_k = ku_k - u((j - 2) * ndof_1d + (i - 3) + 1)  ! west
        if (i < N - 1) ku_k = ku_k - u((j - 2) * ndof_1d + (i - 1) + 1)  ! east
        if (j > 2)   ku_k = ku_k - u((j - 3) * ndof_1d + (i - 2) + 1)  ! south
        if (j < N - 1) ku_k = ku_k - u((j - 1) * ndof_1d + (i - 2) + 1)  ! north

        F(k) = ku_k - h2lam * exp(u(k))
        fnorm = max(fnorm, abs(F(k)))
      end do
    end do
  end subroutine compute_residual

  !> Assemble the Jacobian J(u) = K - h2lam * diag(exp(u)) in COO format.
  !>
  !> Entries are emitted in a fixed row-major order (diagonal first in each
  !> row, then west/east/south/north) so that successive calls with the same
  !> N produce identical (rnr, snr) patterns -- a requirement for
  !> IFLAG(4) = 2 structural reuse.
  subroutine assemble_jacobian(N, h2lam, u, a, rnr, snr, nz)
    integer, intent(in) :: N
    real(dp), intent(in) :: h2lam
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: a(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    integer :: ndof_1d, i, j, k_row

    ndof_1d = N - 2
    nz = 0

    do j = 2, N - 1
      do i = 2, N - 1
        k_row = (j - 2) * ndof_1d + (i - 2) + 1

        ! Diagonal: 4 - h2lam * exp(u_k)
        nz = nz + 1
        rnr(nz) = k_row
        snr(nz) = k_row
        a(nz) = 4.0_dp - h2lam * exp(u(k_row))

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
  end subroutine assemble_jacobian

  ! ---------------------------------------------------------------
  ! Main driver: Newton iteration with structural LU reuse
  ! ---------------------------------------------------------------
  subroutine run(N, lam, max_iter, outfile)
    use y12m, only: y12mb, y12mc, y12md
    integer, intent(in) :: N, max_iter
    real(dp), intent(in) :: lam
    character(len=*), intent(in) :: outfile

    integer :: ndof_1d, ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h, h2lam, fnorm
    real(dp), allocatable :: a(:), pivot(:), u(:), F(:), delta_u(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10), ifail
    integer :: i, j, k, iter, funit
    integer :: t0, t1, clock_rate
    integer :: t_order_first, t_order_rest, t_factor, t_solve
    real(dp) :: s_order_first, s_order_rest, s_factor, s_solve
    real(dp) :: t0_iter, t1_iter
    logical :: converged
    real(dp), parameter :: tol = 1.0e-10_dp

    ndof_1d = N - 2
    ndof = ndof_1d * ndof_1d
    h = 1.0_dp / real(N - 1, dp)
    h2lam = h * h * lam

    write(output_unit, '(a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
    write(output_unit, '(a,es10.3)') 'lambda        : ', lam
    write(output_unit, '(a,i0)')     'Max Newton it : ', max_iter

    ! Workspace
    nz_max = 5 * ndof
    nn = max(3 * nz_max, 2 * (ndof_1d + 2) * ndof)
    nn1 = nn
    iha = ndof
    allocate(a(nn), pivot(ndof), u(ndof), F(ndof), delta_u(ndof))
    allocate(snr(nn), rnr(nn1), ha(ndof, 11))

    call system_clock(count_rate=clock_rate)

    ! Initial guess: u = 0
    u = 0.0_dp

    t_order_first = 0
    t_order_rest  = 0
    t_factor      = 0
    t_solve       = 0
    converged = .false.

    write(output_unit, '(/,a)') 'Newton iteration:'
    write(output_unit, '(a)') '  iter   ||F||_inf   t_y12mb(s)  t_y12mc(s)  t_y12md(s)'

    ! Initialise flags and tolerances ONCE outside the Newton loop.
    ! iflag(9) and iflag(10) are written by y12mc (iflag(4)=1) with the
    ! row/column fill-in counts and MUST be preserved across Newton steps
    ! so that y12mc (iflag(4)=2) can reuse them.
    iflag = 0
    iflag(2) = 3        ! Markowitz search width
    iflag(3) = 1        ! use column interchanges
    iflag(5) = 1        ! discard L (single RHS per Newton step)
    aflag(1) = 16.0_dp
    aflag(2) = 1.0e-12_dp
    aflag(3) = 1.0e+16_dp
    aflag(4) = 1.0e-12_dp
    aflag(5:8) = 0.0_dp

    do iter = 1, max_iter

      ! Compute residual F(u^k) and its norm
      call compute_residual(N, h2lam, u, F, fnorm)

      if (iter == 1) then
        write(output_unit, '(2x,i4,3x,es10.3)') 0, fnorm
      end if

      if (fnorm < tol) then
        converged = .true.
        exit
      end if

      ! RHS for the Newton system: delta_u = J^{-1} (-F)
      delta_u = -F

      ! Assemble Jacobian J(u^k) -- same (rnr, snr) pattern every iteration
      call assemble_jacobian(N, h2lam, u, a, rnr, snr, nz)

      ! ---------------------------------------------------------------
      ! Ordering step: IFLAG(4) = 1 on first iteration (compute and
      ! store Markowitz ordering), IFLAG(4) = 2 on subsequent iterations
      ! (reuse stored ordering from HA).
      ! ---------------------------------------------------------------
      iflag(4) = merge(1, 2, iter == 1)   ! 1: compute ordering; 2: reuse it

      call system_clock(t0)
      call y12mb(ndof, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
      call system_clock(t1)
      if (iter == 1) then
        t_order_first = t_order_first + (t1 - t0)
      else
        t_order_rest = t_order_rest + (t1 - t0)
      end if
      t0_iter = real(t1 - t0, dp)   ! store for per-iter output

      if (ifail /= 0) then
        write(error_unit, '(a,i0,a,i0)') &
            'ERROR: y12mb iter ', iter, ' returned IFAIL = ', ifail
        stop 1, quiet = .true.
      end if

      ! Factorize J(u^k); y12mc applies L^{-1} to delta_u = -F
      call system_clock(t0)
      call y12mc(ndof, nz, a, snr, nn, rnr, nn1, pivot, delta_u, ha, iha, &
          aflag, iflag, ifail)
      call system_clock(t1)
      t_factor = t_factor + (t1 - t0)
      t1_iter = real(t1 - t0, dp)   ! store for per-iter output

      if (ifail /= 0) then
        write(error_unit, '(a,i0,a,i0)') &
            'ERROR: y12mc iter ', iter, ' returned IFAIL = ', ifail
        stop 1, quiet = .true.
      end if

      ! Solve for delta_u (backward substitution with U)
      call system_clock(t0)
      call y12md(ndof, a, nn, delta_u, pivot, snr, ha, iha, iflag, ifail)
      call system_clock(t1)
      t_solve = t_solve + (t1 - t0)

      if (ifail /= 0) then
        write(error_unit, '(a,i0,a,i0)') &
            'ERROR: y12md iter ', iter, ' returned IFAIL = ', ifail
        stop 1, quiet = .true.
      end if

      ! Newton update
      u = u + delta_u

      ! Compute and print residual after update
      call compute_residual(N, h2lam, u, F, fnorm)
      write(output_unit, '(2x,i4,3x,es10.3,3(2x,es10.3))') iter, fnorm, &
          t0_iter / real(clock_rate, dp), &
          t1_iter / real(clock_rate, dp), &
          real(t1 - t0, dp) / real(clock_rate, dp)

    end do

    if (converged) then
      write(output_unit, '(a,i0,a,es10.3)') &
          'Converged in ', iter - 1, ' iterations,  ||F|| = ', fnorm
    else
      write(output_unit, '(a,i0,a)') &
          'WARNING: did NOT converge after ', max_iter, ' iterations'
    end if

    ! =====================================================================
    ! Write solution
    ! =====================================================================
    open(newunit=funit, file=outfile, status='unknown', action='write')
    write(funit, '(a)') '# Bratu problem -- Newton iteration with structural LU reuse'
    write(funit, '(a,i0,a,es10.3)') '# N=', N, '  lambda=', lam
    write(funit, '(a)') '# Columns: x  y  u_numerical'
    do j = 1, N
      do i = 1, N
        if (i == 1 .or. i == N .or. j == 1 .or. j == N) then
          write(funit, '(3(1x,es14.6))') &
              real(i - 1, dp) * h, real(j - 1, dp) * h, 0.0_dp
        else
          k = (j - 2) * ndof_1d + (i - 2) + 1
          write(funit, '(3(1x,es14.6))') &
              real(i - 1, dp) * h, real(j - 1, dp) * h, u(k)
        end if
      end do
      write(funit, *)
    end do
    close(funit)
    write(output_unit, '(a)') 'Solution written to ' // trim(outfile)

    ! =====================================================================
    ! Timing summary
    ! =====================================================================
    s_order_first = real(t_order_first, dp) / real(clock_rate, dp)
    s_order_rest  = real(t_order_rest,  dp) / real(clock_rate, dp)
    s_factor      = real(t_factor,      dp) / real(clock_rate, dp)
    s_solve       = real(t_solve,       dp) / real(clock_rate, dp)

    write(output_unit, '(/,a)') 'Timing summary (structural reuse advantage):'
    write(output_unit, '(a,g14.6,a)') '  y12mb iter 1 (compute ordering): ', s_order_first, ' s'
    if (iter > 2) then
      write(output_unit, '(a,i0,a,g14.6,a)') &
          '  y12mb iters 2..', iter - 1, ' (reuse ordering): ', s_order_rest, ' s total'
      if (s_order_first > 0.0_dp .and. iter > 2) then
        write(output_unit, '(a,g14.6,a)') &
            '  Average speedup per call    : ', &
            s_order_first / (s_order_rest / real(iter - 2, dp)), 'x'
      end if
    end if
    write(output_unit, '(a,g14.6,a)') '  y12mc total (factorizations)  : ', s_factor, ' s'
    write(output_unit, '(a,g14.6,a)') '  y12md total (back-sub solves) : ', s_solve,  ' s'

  end subroutine run
end module bratu_solver

! ==============================================================
! Main program: parse CLI, call run.
! ==============================================================
program newton_bratu
  use bratu_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N, max_iter, ios, iarg, nargs, pos_count
  real(kind(1.0d0)) :: lam
  character(len=256) :: arg, outfile

  N        = 22
  lam      = 1.0d0
  max_iter = 20
  outfile  = 'newton_bratu.dat'

  nargs     = command_argument_count()
  iarg      = 1
  pos_count = 0

  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: newton_bratu [--help] [N] [lambda] [max_iter] [output_file]'
      write(*, '(a)') '  N           total grid size (>= 3, default 22)'
      write(*, '(a)') '  lambda      reaction parameter (default 1.0, must be < 6.808)'
      write(*, '(a)') '  max_iter    max Newton iterations (default 20)'
      write(*, '(a)') '  output_file output file (default: newton_bratu.dat)'
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
        read(arg, *, iostat=ios) lam
        if (ios /= 0 .or. lam <= 0.0d0) then
          write(error_unit, '(a)') 'ERROR: lambda must be a positive real'
          stop 1
        end if
      case (3)
        read(arg, *, iostat=ios) max_iter
        if (ios /= 0 .or. max_iter < 1) then
          write(error_unit, '(a)') 'ERROR: max_iter must be a positive integer'
          stop 1
        end if
      case (4)
        outfile = trim(arg)
      end select
    end if

    iarg = iarg + 1
  end do

  call run(N, lam, max_iter, trim(outfile))

end program newton_bratu
