! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12ma_dp.f90 - Tests for y12ma (double-precision black-box solver).
!
! This program tests the y12ma generic interface (resolving to y12maf)
! by solving five sparse linear systems Ax=b of sizes n=6..10.
! All systems have the known solution x=[1,...,1], constructed by setting
! b = A*[1,...,1] so that verification is straightforward.
!
!   n=6,7:  6x6 and 7x7 tridiagonal (diag=3, off-diag=-1), solutions verified.
!   n=8:    8x8 arrow matrix (full first row/column plus diagonal rest).
!   n=9,10: 9x9 and 10x10 tridiagonal, solutions verified.
!
! Pass criteria: IFAIL=0 and max|x_i-1|<1e-10 for all systems.
!
program test_y12ma_dp
  use y12m
  use y12m_test_helpers, only: build_tridiag, build_arrow, check_solution
  implicit none

  integer, parameter :: NMAX  = 12
  integer, parameter :: NNP   = 400
  integer, parameter :: NN1P  = 400

  double precision :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer          :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer          :: n, z, nn, nn1, iha, ifail, nfail

  nfail = 0

  ! --- n=6,7: tridiagonal systems ---
  do n = 6, 7
    call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (6)
      call check_solution('y12ma_dp n=6 tridiag', n, b, ifail, 1.0d-10, nfail)
    case (7)
      call check_solution('y12ma_dp n=7 tridiag', n, b, ifail, 1.0d-10, nfail)
    end select
  end do

  ! --- n=8: arrow matrix ---
  n = 8
  call build_arrow(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
  call check_solution('y12ma_dp n=8 arrow', n, b, ifail, 1.0d-10, nfail)

  ! --- n=9,10: tridiagonal systems ---
  do n = 9, 10
    call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (9)
      call check_solution('y12ma_dp n=9 tridiag', n, b, ifail, 1.0d-10, nfail)
    case (10)
      call check_solution('y12ma_dp n=10 tridiag', n, b, ifail, 1.0d-10, nfail)
    end select
  end do

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_dp tests PASSED'

end program test_y12ma_dp

