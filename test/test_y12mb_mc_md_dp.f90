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
  call build_tridiag_dp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_dp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_dp('y12mb_mc_md_dp n=6 tridiag', n, b, ifail, 1.0d-10, nfail)

  ! --- n=7: arrow ---
  n = 7
  call build_arrow_dp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_dp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_dp('y12mb_mc_md_dp n=7 arrow', n, b, ifail, 1.0d-10, nfail)

  ! --- n=8: tridiagonal ---
  n = 8
  call build_tridiag_dp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_dp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_dp('y12mb_mc_md_dp n=8 tridiag', n, b, ifail, 1.0d-10, nfail)

  ! --- n=9: tridiagonal ---
  n = 9
  call build_tridiag_dp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_dp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_dp('y12mb_mc_md_dp n=9 tridiag', n, b, ifail, 1.0d-10, nfail)

  ! --- n=10: pentadiagonal (bandwidth 2) ---
  n = 10
  call build_pentadiag_dp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_dp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_dp('y12mb_mc_md_dp n=10 pentadiag', n, b, ifail, 1.0d-10, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mb_mc_md_dp tests PASSED'

contains

  ! Execute the three-step prepare->factorize->solve sequence (double precision).
  subroutine solve3_dp(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, &
      iha, aflag, iflag, ifail)
    integer,          intent(in)    :: n, z, nn, nn1, iha
    double precision, intent(inout) :: a(nn), b(n), pivot(n), aflag(8)
    integer,          intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    integer,          intent(out)   :: ifail

    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 1

    aflag(1) = 16.0d0
    aflag(2) = 1.0d-12
    aflag(3) = 1.0d+16
    aflag(4) = 1.0d-12

    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) return

    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
        aflag, iflag, ifail)
    if (ifail /= 0) return

    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  end subroutine solve3_dp

  ! n-by-n tridiagonal (diag=3, off-diag=-1), double precision, x=[1,...,1].
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

  ! n-by-n arrow matrix, double precision, solution x=[1,...,1].
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

  ! n-by-n pentadiagonal (diag=5, dist-1 off-diag=-1, dist-2 off-diag=-0.5),
  ! double precision, solution x=[1,...,1].
  subroutine build_pentadiag_dp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,          intent(in)  :: n, nnmax, nn1max
    integer,          intent(out) :: z
    double precision, intent(out) :: a(nnmax), b(n)
    integer,          intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i ; a(z) = 5.0d0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0d0
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0d0
    end do
    do i = 3, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-2 ; a(z) = -0.5d0
    end do
    do i = 1, n-2
      z = z + 1 ; rnr(z) = i ; snr(z) = i+2 ; a(z) = -0.5d0
    end do
    b(1:n) = 0.0d0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_pentadiag_dp

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

end program test_y12mb_mc_md_dp
