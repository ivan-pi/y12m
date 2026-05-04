! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12ma_sp.f90 - Tests for y12ma (single-precision black-box solver).
!
! This program tests the y12ma generic interface (resolving to y12mae)
! by solving five sparse linear systems Ax=b of sizes n=1..5.
! All solvable systems use the known solution x=[1,...,1], constructed by
! setting b = A*[1,...,1] so that verification is straightforward.
!
!   n=1: expected IFAIL=12 — the package requires n>=2; this tests the
!        error-diagnostic path of y12ma.
!   n=2: 2x2 tridiagonal (diag=3, off-diag=-1), solution verified.
!   n=3: 3x3 tridiagonal, solution verified.
!   n=4: 4x4 tridiagonal, solution verified.
!   n=5: 5x5 arrow matrix (full first row/column plus diagonal rest),
!        solution verified.
!
! Pass criteria: IFAIL=12 for n=1, IFAIL=0 and max|x_i-1|<1e-4 for n=2..5.
!
program test_y12ma_sp
  use y12m
  use y12m_test_helpers, only: build_tridiag, build_arrow, build_diag, check_solution
  implicit none

  integer, parameter :: NMAX  = 10
  integer, parameter :: NNP   = 200
  integer, parameter :: NN1P  = 200

  real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer :: n, z, nn, nn1, iha, ifail, nfail

  nfail = 0

  ! --- n=1: verify expected error IFAIL=12 (N<2 not supported by package) ---
  n = 1
  call build_diag(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
  if (ifail == 12) then
    write(*,'(a)') 'PASS y12ma_sp n=1 error-diag: ifail=12 (N<2) as expected'
  else
    write(*,'(a,i0)') 'FAIL y12ma_sp n=1 error-diag: expected ifail=12, got ', &
        ifail
    nfail = nfail + 1
  end if

  ! --- n=2..4: tridiagonal, diag=3, off-diag=-1 ---
  do n = 2, 4
    call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (2)
      call check_solution('y12ma_sp n=2 tridiag', n, b, ifail, 1.0e-4, nfail)
    case (3)
      call check_solution('y12ma_sp n=3 tridiag', n, b, ifail, 1.0e-4, nfail)
    case (4)
      call check_solution('y12ma_sp n=4 tridiag', n, b, ifail, 1.0e-4, nfail)
    end select
  end do

  ! --- n=5: 5x5 arrow matrix ---
  n = 5
  call build_arrow(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
  call check_solution('y12ma_sp n=5 arrow', n, b, ifail, 1.0e-4, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_sp tests PASSED'

end program test_y12ma_sp
