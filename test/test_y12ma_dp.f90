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
    call build_tridiag_dp(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (6)
      call check_dp('y12ma_dp n=6 tridiag', n, b, ifail, 1.0d-10, nfail)
    case (7)
      call check_dp('y12ma_dp n=7 tridiag', n, b, ifail, 1.0d-10, nfail)
    end select
  end do

  ! --- n=8: arrow matrix ---
  n = 8
  call build_arrow_dp(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
  call check_dp('y12ma_dp n=8 arrow', n, b, ifail, 1.0d-10, nfail)

  ! --- n=9,10: tridiagonal systems ---
  do n = 9, 10
    call build_tridiag_dp(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (9)
      call check_dp('y12ma_dp n=9 tridiag', n, b, ifail, 1.0d-10, nfail)
    case (10)
      call check_dp('y12ma_dp n=10 tridiag', n, b, ifail, 1.0d-10, nfail)
    end select
  end do

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_dp tests PASSED'

contains

  ! Build an n-by-n tridiagonal (diag=3, off-diag=-1), double precision,
  ! solution x=[1,...,1].
  subroutine build_tridiag_dp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,          intent(in)  :: n, nnmax, nn1max
    integer,          intent(out) :: z
    double precision, intent(out) :: a(nnmax), b(n)
    integer,          intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) =  3.0d0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0d0
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0d0
    end do
    b(1:n) = 0.0d0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_tridiag_dp

  ! Build an n-by-n arrow matrix, double precision, solution x=[1,...,1].
  subroutine build_arrow_dp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,          intent(in)  :: n, nnmax, nn1max
    integer,          intent(out) :: z
    double precision, intent(out) :: a(nnmax), b(n)
    integer,          intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i, j
    z = 0
    do j = 1, n
      z = z + 1 ; rnr(z) = 1 ; snr(z) = j ; a(z) = 1.0d0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = 1   ; a(z) = 1.0d0
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) = dble(i + 1)
    end do
    b(1:n) = 0.0d0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_arrow_dp

  ! Check IFAIL=0 and max|b(1:n)-1|<tol; increment nfail on failure.
  subroutine check_dp(label, n, b, ifail, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, ifail
    double precision, intent(in)    :: b(n), tol
    integer,          intent(inout) :: nfail
    double precision :: err
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ' ifail=', ifail
      nfail = nfail + 1
      return
    end if
    err = maxval(abs(b(1:n) - 1.0d0))
    if (err > tol) then
      write(*,'(3a,es12.5)') 'FAIL ', label, ' max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es12.5)') 'PASS ', label, ' max_err=', err
    end if
  end subroutine check_dp

end program test_y12ma_dp
