! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12ma_rand_dp.f90 - Random 9x9 matrix test for y12ma (double-precision).
!
! Generates a random 9x9 diagonally-dominant sparse matrix A and three random
! right-hand-side vectors using seed 20260417 (yyyymmdd format).
! Solves each system with the y12ma black-box solver and verifies correctness
! in two ways:
!
!   1. Residual check: max|b_i - A*x_i| < 1e-10 for each RHS column,
!      using a sparse triplet matrix-vector product.
!   2. Dense cross-check: solves the equivalent dense system with LAPACK
!      dgesv (using the built-in matmul to form b_dense = A_dense * x_dense)
!      and asserts max|x_sparse - x_dense| < 1e-10 for each column.
!
! Pass criteria: IFAIL=0, max residual < 1e-10, and max diff vs dgesv < 1e-10.
!
program test_y12ma_rand_dp
  use y12m
  implicit none

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      integer,          intent(in)    :: n, nrhs, lda, ldb
      double precision, intent(inout) :: a(lda, n)
      integer,          intent(out)   :: ipiv(n)
      double precision, intent(inout) :: b(ldb, nrhs)
      integer,          intent(out)   :: info
    end subroutine dgesv
  end interface

  integer, parameter :: N    = 9
  integer, parameter :: NRHS = 3
  integer, parameter :: NNP  = 200
  integer, parameter :: NN1P = 200
  integer, parameter :: NMAX = 12

  double precision :: a(NNP), a0(NNP), pivot(NMAX), aflag(8)
  double precision :: b(N), b0(N, NRHS), x_sparse(N, NRHS)
  double precision :: adense(N, N), alu(N, N), x_dense(N, NRHS)
  double precision :: ax(N), resid(N)
  integer :: snr(NNP), snr0(NNP), rnr(NN1P), rnr0(NN1P)
  integer :: ha(NMAX, 11), iflag(10), ipiv(N)
  integer :: z, nn, nn1, iha, ifail, nfail, irhs, info
  integer :: seed_size
  integer, allocatable :: seed_array(:)

  nfail = 0

  ! --- Seed the random number generator (date: 20260417) ---
  call random_seed(size=seed_size)
  allocate(seed_array(seed_size))
  seed_array = 20260417
  call random_seed(put=seed_array)
  deallocate(seed_array)

  ! --- Generate random 9x9 sparse diagonally-dominant matrix ---
  call build_rand_dp(N, NNP, NN1P, z, a0, snr0, rnr0, adense)

  ! --- Generate 3 random RHS vectors with values in [-1, 1] ---
  call random_number(b0)
  b0 = 2.0d0 * b0 - 1.0d0

  ! --- Sparse solve: one y12ma call per RHS column ---
  do irhs = 1, NRHS
    a(1:z)   = a0(1:z)
    snr(1:z) = snr0(1:z)
    rnr(1:z) = rnr0(1:z)
    b(1:N)   = b0(1:N, irhs)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(N, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0,a,i0)') 'FAIL y12ma_rand_dp rhs=', irhs, &
          ' ifail=', ifail
      nfail = nfail + 1
      cycle
    end if
    x_sparse(1:N, irhs) = b(1:N)
    ! Residual: b_orig - A * x_sparse
    call sparse_matmul(N, z, a0, snr0, rnr0, x_sparse(1:N, irhs), ax)
    resid = b0(1:N, irhs) - ax
    call check_dp('y12ma_rand_dp sparse resid rhs', irhs, N, resid, &
        1.0d-10, nfail)
  end do

  ! --- Dense solve: LAPACK dgesv ---
  alu     = adense
  x_dense = b0
  call dgesv(N, NRHS, alu, N, ipiv, x_dense, N, info)
  if (info /= 0) then
    write(*,'(a,i0)') 'FAIL y12ma_rand_dp dgesv info=', info
    nfail = nfail + 1
  else
    ! --- Cross-check: sparse solution vs dense solution ---
    do irhs = 1, NRHS
      call check_close('y12ma_rand_dp sparse vs dense rhs', irhs, N, &
          x_sparse(1:N, irhs), x_dense(1:N, irhs), 1.0d-10, nfail)
    end do
  end if

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_rand_dp tests PASSED'

contains

  ! Compute y = A*x where A is given in triplet format (a, rnr, snr).
  subroutine sparse_matmul(n, z, a, snr, rnr, x, y)
    integer,          intent(in)  :: n, z
    double precision, intent(in)  :: a(*), x(n)
    integer,          intent(in)  :: snr(*), rnr(*)
    double precision, intent(out) :: y(n)
    integer :: k
    y(1:n) = 0.0d0
    do k = 1, z
      y(rnr(k)) = y(rnr(k)) + a(k) * x(snr(k))
    end do
  end subroutine sparse_matmul

  ! Build a random n-by-n diagonally-dominant sparse matrix in triplet format.
  ! Each off-diagonal entry (i,j) is included with probability 1/3,
  ! with a value drawn uniformly from [-1, 1].
  ! The diagonal entry in row i is set to sum(|off-diag row i|) + 2,
  ! guaranteeing strict diagonal dominance and hence non-singularity.
  subroutine build_rand_dp(n, nnmax, nn1max, z, a, snr, rnr, adense)
    integer,          intent(in)  :: n, nnmax, nn1max
    integer,          intent(out) :: z
    double precision, intent(out) :: a(nnmax)
    integer,          intent(out) :: snr(nnmax), rnr(nn1max)
    double precision, intent(out) :: adense(n, n)
    double precision :: row_sum(n), rval
    integer :: i, j
    z       = 0
    row_sum = 0.0d0
    adense  = 0.0d0
    ! Off-diagonal entries: include with probability 1/3
    do i = 1, n
      do j = 1, n
        if (i == j) cycle
        call random_number(rval)
        if (rval < 1.0d0/3.0d0) then
          call random_number(rval)
          rval = 2.0d0 * rval - 1.0d0
          z = z + 1
          rnr(z) = i ; snr(z) = j ; a(z) = rval
          adense(i, j) = rval
          row_sum(i) = row_sum(i) + abs(rval)
        end if
      end do
    end do
    ! Diagonal: strictly diagonally dominant
    do i = 1, n
      z = z + 1
      rnr(z) = i ; snr(z) = i ; a(z) = row_sum(i) + 2.0d0
      adense(i, i) = row_sum(i) + 2.0d0
    end do
  end subroutine build_rand_dp

  ! Check max|resid| < tol; increment nfail and print PASS/FAIL.
  subroutine check_dp(label, irhs, n, resid, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: irhs, n
    double precision, intent(in)    :: resid(n), tol
    integer,          intent(inout) :: nfail
    double precision :: err
    err = maxval(abs(resid))
    if (err > tol) then
      write(*,'(3a,i0,a,es12.5)') 'FAIL ', label, '=', irhs, &
          ' max_resid=', err
      nfail = nfail + 1
    else
      write(*,'(3a,i0,a,es12.5)') 'PASS ', label, '=', irhs, &
          ' max_resid=', err
    end if
  end subroutine check_dp

  ! Check max|x1 - x2| < tol; increment nfail and print PASS/FAIL.
  subroutine check_close(label, irhs, n, x1, x2, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: irhs, n
    double precision, intent(in)    :: x1(n), x2(n), tol
    integer,          intent(inout) :: nfail
    double precision :: err
    err = maxval(abs(x1 - x2))
    if (err > tol) then
      write(*,'(3a,i0,a,es12.5)') 'FAIL ', label, '=', irhs, &
          ' max_diff=', err
      nfail = nfail + 1
    else
      write(*,'(3a,i0,a,es12.5)') 'PASS ', label, '=', irhs, &
          ' max_diff=', err
    end if
  end subroutine check_close

end program test_y12ma_rand_dp
