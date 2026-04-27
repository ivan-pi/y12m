! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4-5
!
! test_y12ma_rand.F90 - Random 9x9 matrix test for y12ma (single or double precision).
!
! Generates a random 9x9 diagonally-dominant sparse matrix A and three random
! right-hand-side vectors using seed 20260417 (yyyymmdd format).
! Solves each system with the y12ma black-box solver and verifies correctness
! in two ways:
!
!   1. Residual check: max|b_i - A*x_i| < tol for each RHS column,
!      using a sparse triplet matrix-vector product.
!   2. Dense cross-check: solves the equivalent dense system with LAPACK
!      sgesv/dgesv (via the lapack_gesv generic interface) and asserts
!      max|x_sparse - x_dense| < tol for each column.
!
! Compile with -DTEST_SINGLE_PRECISION for single-precision (sgesv); omit
! for double-precision (dgesv).
!
program test_y12ma_rand
  use y12m
  implicit none

#ifdef TEST_SINGLE_PRECISION
  integer, parameter :: wp = kind(1.0e0)
#else
  integer, parameter :: wp = kind(1.0d0)
#endif

#ifdef TEST_SINGLE_PRECISION
  real(wp), parameter :: tol = 1.0e-4_wp
#else
  real(wp), parameter :: tol = 1.0e-10_wp
#endif

  interface lapack_gesv
    subroutine sgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      implicit none
      integer,          intent(in)    :: n, nrhs, lda, ldb
      real(kind(1.0e0)), intent(inout) :: a(lda, n)
      integer,          intent(out)   :: ipiv(n)
      real(kind(1.0e0)), intent(inout) :: b(ldb, nrhs)
      integer,          intent(out)   :: info
    end subroutine sgesv
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      implicit none
      integer,          intent(in)    :: n, nrhs, lda, ldb
      real(kind(1.0d0)), intent(inout) :: a(lda, n)
      integer,          intent(out)   :: ipiv(n)
      real(kind(1.0d0)), intent(inout) :: b(ldb, nrhs)
      integer,          intent(out)   :: info
    end subroutine dgesv
  end interface lapack_gesv

  integer, parameter :: N    = 9
  integer, parameter :: NRHS = 3
  integer, parameter :: NNP  = 200
  integer, parameter :: NN1P = 200
  integer, parameter :: NMAX = 12

  real(wp) :: a(NNP), a0(NNP), pivot(NMAX), aflag(8)
  real(wp) :: b(N), b0(N, NRHS), x_sparse(N, NRHS)
  real(wp) :: adense(N, N), alu(N, N), x_dense(N, NRHS)
  real(wp) :: ax(N), resid(N)
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
  call build_rand(N, NNP, NN1P, z, a0, snr0, rnr0, adense)

  ! --- Generate 3 random RHS vectors with values in [-1, 1] ---
  call random_number(b0)
  b0 = 2.0_wp * b0 - 1.0_wp

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
      write(*,'(a,i0,a,i0)') 'FAIL y12ma_rand rhs=', irhs, &
          ' ifail=', ifail
      nfail = nfail + 1
      cycle
    end if
    x_sparse(1:N, irhs) = b(1:N)
    ! Residual: b_orig - A * x_sparse
    call sparse_matmul(N, z, a0, snr0, rnr0, x_sparse(1:N, irhs), ax)
    resid = b0(1:N, irhs) - ax
    call check_resid('y12ma_rand sparse resid rhs', irhs, N, resid, tol, nfail)
  end do

  ! --- Dense solve: LAPACK sgesv/dgesv via lapack_gesv ---
  alu     = adense
  x_dense = b0
  call lapack_gesv(N, NRHS, alu, N, ipiv, x_dense, N, info)
  if (info /= 0) then
    write(*,'(a,i0)') 'FAIL y12ma_rand lapack_gesv info=', info
    nfail = nfail + 1
  else
    ! --- Cross-check: sparse solution vs dense solution ---
    do irhs = 1, NRHS
      call check_close('y12ma_rand sparse vs dense rhs', irhs, N, &
          x_sparse(1:N, irhs), x_dense(1:N, irhs), tol, nfail)
    end do
  end if

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_rand tests PASSED'

contains

  ! Compute y = A*x where A is given in triplet format (a, rnr, snr).
  subroutine sparse_matmul(n, z, a, snr, rnr, x, y)
    integer,   intent(in)  :: n, z
    real(wp),  intent(in)  :: a(*), x(n)
    integer,   intent(in)  :: snr(*), rnr(*)
    real(wp),  intent(out) :: y(n)
    integer :: k
    y(1:n) = 0.0_wp
    do k = 1, z
      y(rnr(k)) = y(rnr(k)) + a(k) * x(snr(k))
    end do
  end subroutine sparse_matmul

  ! Build a random n-by-n diagonally-dominant sparse matrix in triplet format.
  ! Each off-diagonal entry (i,j) is included with probability 1/3,
  ! with a value drawn uniformly from [-1, 1].
  ! The diagonal entry in row i is set to sum(|off-diag row i|) + 2,
  ! guaranteeing strict diagonal dominance and hence non-singularity.
  subroutine build_rand(n, nnmax, nn1max, z, a, snr, rnr, adense)
    integer,  intent(in)  :: n, nnmax, nn1max
    integer,  intent(out) :: z
    real(wp), intent(out) :: a(nnmax)
    integer,  intent(out) :: snr(nnmax), rnr(nn1max)
    real(wp), intent(out) :: adense(n, n)
    real(wp) :: row_sum(n), rval
    integer :: i, j
    z       = 0
    row_sum = 0.0_wp
    adense  = 0.0_wp
    ! Off-diagonal entries: include with probability 1/3
    do i = 1, n
      do j = 1, n
        if (i == j) cycle
        call random_number(rval)
        if (rval < 1.0_wp/3.0_wp) then
          call random_number(rval)
          rval = 2.0_wp * rval - 1.0_wp
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
      rnr(z) = i ; snr(z) = i ; a(z) = row_sum(i) + 2.0_wp
      adense(i, i) = row_sum(i) + 2.0_wp
    end do
  end subroutine build_rand

  ! Check max|resid| < tol; increment nfail and print PASS/FAIL.
  subroutine check_resid(label, irhs, n, resid, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: irhs, n
    real(wp),         intent(in)    :: resid(n), tol
    integer,          intent(inout) :: nfail
    real(wp) :: err
    err = maxval(abs(resid))
    if (err > tol) then
      write(*,'(3a,i0,a,es12.5)') 'FAIL ', label, '=', irhs, &
          ' max_resid=', err
      nfail = nfail + 1
    else
      write(*,'(3a,i0,a,es12.5)') 'PASS ', label, '=', irhs, &
          ' max_resid=', err
    end if
  end subroutine check_resid

  ! Check max|x1 - x2| < tol; increment nfail and print PASS/FAIL.
  subroutine check_close(label, irhs, n, x1, x2, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: irhs, n
    real(wp),         intent(in)    :: x1(n), x2(n), tol
    integer,          intent(inout) :: nfail
    real(wp) :: err
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

end program test_y12ma_rand
