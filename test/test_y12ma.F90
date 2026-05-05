! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12ma.F90 - Tests for y12ma (black-box sparse solver).
!
! This program tests the y12ma generic interface (resolving to y12mae for
! single precision and y12maf for double precision) by solving sparse linear
! systems Ax=b of sizes n=1..10.
!
! Compile with -DTEST_SINGLE_PRECISION to use single precision (wp = kind(1.0e0)),
! otherwise double precision (wp = kind(1.0d0)) is used.
!
!   n=1:       expected IFAIL=12 -- the package requires n>=2; this tests the
!              error-diagnostic path of y12ma.
!   n=2,3,4:   tridiagonal (diag=3, off-diag=-1), solutions verified.
!   n=5:       5x5 arrow matrix (full first row/column plus diagonal rest).
!   n=6,7:     tridiagonal, solutions verified.
!   n=8:       8x8 arrow matrix.
!   n=9,10:    tridiagonal, solutions verified.
!
! Pass criteria: IFAIL=12 for n=1; IFAIL=0 and max|x_i-1| <= tol for n=2..10.
!   tol = 1e-4  (single precision)
!   tol = 1e-10 (double precision)
!
#ifdef TEST_SINGLE_PRECISION
program test_y12ma_sp
#else
program test_y12ma_dp
#endif
  use y12m, only: y12ma
  use y12m_test_helpers, only: build_tridiag, build_arrow, check_solution
  implicit none

#ifdef TEST_SINGLE_PRECISION
  integer, parameter :: wp = kind(1.0e0)
  real(wp), parameter :: tol_sol = 1.0e-4_wp
#else
  integer, parameter :: wp = kind(1.0d0)
  real(wp), parameter :: tol_sol = 1.0e-10_wp
#endif

  integer, parameter :: NMAX = 12
  integer, parameter :: NNP  = 400
  integer, parameter :: NN1P = 400

  real(wp) :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer :: n, z, nn, nn1, iha, ifail, nfail

  nfail = 0

  ! --- n=1: verify expected error IFAIL=12 (N<2 not supported by package) ---
  ! build_diag from y12m_test_helpers is sp-only; construct the 1x1 system
  ! manually so this case works for both precisions.
  n = 1 ; z = 1
  a(1) = 1.0_wp ; snr(1) = 1 ; rnr(1) = 1 ; b(1) = 1.0_wp
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
  if (ifail == 12) then
#ifdef TEST_SINGLE_PRECISION
    write(*,'(a)') 'PASS y12ma_sp n=1 error-diag: ifail=12 (N<2) as expected'
#else
    write(*,'(a)') 'PASS y12ma_dp n=1 error-diag: ifail=12 (N<2) as expected'
#endif
  else
#ifdef TEST_SINGLE_PRECISION
    write(*,'(a,i0)') 'FAIL y12ma_sp n=1 error-diag: expected ifail=12, got ', ifail
#else
    write(*,'(a,i0)') 'FAIL y12ma_dp n=1 error-diag: expected ifail=12, got ', ifail
#endif
    nfail = nfail + 1
  end if

  ! --- n=2..4: tridiagonal (diag=3, off-diag=-1) ---
  do n = 2, 4
    call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (2)
#ifdef TEST_SINGLE_PRECISION
      call check_solution('y12ma_sp n=2 tridiag', n, b, ifail, tol_sol, nfail)
#else
      call check_solution('y12ma_dp n=2 tridiag', n, b, ifail, tol_sol, nfail)
#endif
    case (3)
#ifdef TEST_SINGLE_PRECISION
      call check_solution('y12ma_sp n=3 tridiag', n, b, ifail, tol_sol, nfail)
#else
      call check_solution('y12ma_dp n=3 tridiag', n, b, ifail, tol_sol, nfail)
#endif
    case (4)
#ifdef TEST_SINGLE_PRECISION
      call check_solution('y12ma_sp n=4 tridiag', n, b, ifail, tol_sol, nfail)
#else
      call check_solution('y12ma_dp n=4 tridiag', n, b, ifail, tol_sol, nfail)
#endif
    end select
  end do

  ! --- n=5: arrow matrix ---
  n = 5
  call build_arrow(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
#ifdef TEST_SINGLE_PRECISION
  call check_solution('y12ma_sp n=5 arrow', n, b, ifail, tol_sol, nfail)
#else
  call check_solution('y12ma_dp n=5 arrow', n, b, ifail, tol_sol, nfail)
#endif

  ! --- n=6,7: tridiagonal ---
  do n = 6, 7
    call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (6)
#ifdef TEST_SINGLE_PRECISION
      call check_solution('y12ma_sp n=6 tridiag', n, b, ifail, tol_sol, nfail)
#else
      call check_solution('y12ma_dp n=6 tridiag', n, b, ifail, tol_sol, nfail)
#endif
    case (7)
#ifdef TEST_SINGLE_PRECISION
      call check_solution('y12ma_sp n=7 tridiag', n, b, ifail, tol_sol, nfail)
#else
      call check_solution('y12ma_dp n=7 tridiag', n, b, ifail, tol_sol, nfail)
#endif
    end select
  end do

  ! --- n=8: arrow matrix ---
  n = 8
  call build_arrow(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
#ifdef TEST_SINGLE_PRECISION
  call check_solution('y12ma_sp n=8 arrow', n, b, ifail, tol_sol, nfail)
#else
  call check_solution('y12ma_dp n=8 arrow', n, b, ifail, tol_sol, nfail)
#endif

  ! --- n=9,10: tridiagonal ---
  do n = 9, 10
    call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (9)
#ifdef TEST_SINGLE_PRECISION
      call check_solution('y12ma_sp n=9 tridiag', n, b, ifail, tol_sol, nfail)
#else
      call check_solution('y12ma_dp n=9 tridiag', n, b, ifail, tol_sol, nfail)
#endif
    case (10)
#ifdef TEST_SINGLE_PRECISION
      call check_solution('y12ma_sp n=10 tridiag', n, b, ifail, tol_sol, nfail)
#else
      call check_solution('y12ma_dp n=10 tridiag', n, b, ifail, tol_sol, nfail)
#endif
    end select
  end do

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
#ifdef TEST_SINGLE_PRECISION
  write(*,'(a)') 'All test_y12ma_sp tests PASSED'
#else
  write(*,'(a)') 'All test_y12ma_dp tests PASSED'
#endif

#ifdef TEST_SINGLE_PRECISION
end program test_y12ma_sp
#else
end program test_y12ma_dp
#endif
