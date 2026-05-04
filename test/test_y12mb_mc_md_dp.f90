! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12mb_mc_md_dp.f90 - Tests for the three-step double-precision API:
!   y12mb (prepare) -> y12mc (factorize) -> y12md (solve).
!
! This program tests the generic interfaces y12mb, y12mc, and y12md, each
! resolving to their double-precision variants (y12mbf, y12mcf, y12mdf).
! Five sparse systems of sizes n=6..10 cover different sparsity patterns:
!
!   n=6:  6x6 tridiagonal (diag=3, off-diag=-1), solution verified.
!   n=7:  7x7 arrow matrix (full first row/column plus diagonal rest).
!   n=8:  8x8 tridiagonal, solution verified.
!   n=9:  9x9 tridiagonal, solution verified.
!   n=10: 10x10 pentadiagonal (bandwidth 2, diag=5, dist-1 off-diag=-1,
!         dist-2 off-diag=-0.5), solution verified.
!
! Pass criteria: IFAIL=0 and max|x_i-1|<1e-10 for all systems.
!
program test_y12mb_mc_md_dp
  use y12m
  use y12m_test_helpers, only: build_tridiag, build_arrow, build_pentadiag, &
      check_solution, solve3
  implicit none

  integer, parameter :: NMAX  = 12
  integer, parameter :: NNP   = 500
  integer, parameter :: NN1P  = 500

  double precision :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer          :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer          :: n, z, ifail, nfail

  nfail = 0

  ! --- n=6: tridiagonal ---
  n = 6
  call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_dp n=6 tridiag', n, b, ifail, 1.0d-10, nfail)

  ! --- n=7: arrow ---
  n = 7
  call build_arrow(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_dp n=7 arrow', n, b, ifail, 1.0d-10, nfail)

  ! --- n=8: tridiagonal ---
  n = 8
  call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_dp n=8 tridiag', n, b, ifail, 1.0d-10, nfail)

  ! --- n=9: tridiagonal ---
  n = 9
  call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_dp n=9 tridiag', n, b, ifail, 1.0d-10, nfail)

  ! --- n=10: pentadiagonal (bandwidth 2) ---
  n = 10
  call build_pentadiag(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_solution('y12mb_mc_md_dp n=10 pentadiag', n, b, ifail, 1.0d-10, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mb_mc_md_dp tests PASSED'

end program test_y12mb_mc_md_dp
