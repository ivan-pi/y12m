! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12mb_mc_md_sp.f90 - Tests for the three-step single-precision API:
!   y12mb (prepare) -> y12mc (factorize) -> y12md (solve).
!
! This program tests the generic interfaces y12mb, y12mc, and y12md, each
! resolving to their single-precision variants (y12mbe, y12mce, y12mde).
! Five sparse systems of sizes n=1..5 exercise different code paths:
!
!   n=1: expected IFAIL=12 — y12mb requires n>=2; tests error diagnostics.
!   n=2: 2x2 full matrix [[4,1],[2,5]], b=[5,7], solution x=[1,1].
!   n=3: 3x3 tridiagonal (diag=3, off-diag=-1), solution verified.
!   n=4: 4x4 tridiagonal, solution verified.
!   n=5: 5x5 arrow matrix (full first row/column plus diagonal rest),
!        solution verified.
!
! Pass criteria: IFAIL=12 for n=1, IFAIL=0 and max|x_i-1|<1e-4 for n=2..5.
!
program test_y12mb_mc_md_sp
  use y12m
  use y12m_test_helpers, only: build_tridiag, build_arrow, build_diag, &
      build_full2, check_solution, solve3
  implicit none

  integer, parameter :: NMAX  = 10
  integer, parameter :: NNP   = 200
  integer, parameter :: NN1P  = 200

  real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer :: n, z, ifail, nfail

  nfail = 0

  ! --- n=1: verify expected error IFAIL=12 (N<2 is not supported) ---
  n = 1
  call build_diag(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  if (ifail == 12) then
    write(*,'(a)') 'PASS y12mb_mc_md_sp n=1 error-diag: ifail=12 as expected'
  else
    write(*,'(a,i0)') &
        'FAIL y12mb_mc_md_sp n=1 error-diag: expected ifail=12, got ', ifail
    nfail = nfail + 1
  end if

  ! --- n=2: 2x2 full matrix ---
  n = 2
  call build_full2(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_sp n=2 full', n, b, ifail, 1.0e-4, nfail)

  ! --- n=3: tridiagonal ---
  n = 3
  call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_sp n=3 tridiag', n, b, ifail, 1.0e-4, nfail)

  ! --- n=4: tridiagonal ---
  n = 4
  call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_sp n=4 tridiag', n, b, ifail, 1.0e-4, nfail)

  ! --- n=5: arrow matrix ---
  n = 5
  call build_arrow(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_sp n=5 arrow', n, b, ifail, 1.0e-4, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mb_mc_md_sp tests PASSED'

end program test_y12mb_mc_md_sp
