! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12ma_sp.f90 - Tests for y12mae: single-precision black-box solver.
!
! Five sparse systems of sizes n=1..5 with known solution x=[1,...,1] are
! constructed using tridiagonal (n=1..4) and arrow (n=5) matrices and solved
! via y12mae.  The test passes when IFAIL=0 and max|x_i-1| < tolerance.
!
program test_y12ma_sp
  implicit none

  ! Explicit interface for y12mae.  The actual legacy implementation uses
  ! implicit integer typing for IFAIL (variable name starts with 'i'),
  ! so IFAIL is declared integer here to match the implementation.
  interface
    subroutine y12mae(n, z, a, snr, nn, rnr, nn1, &
        pivot, ha, iha, aflag, iflag, b, ifail)
      implicit none
      integer, intent(in)    :: n, z, nn, nn1, iha
      real,    intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      real,    intent(inout) :: pivot(n)
      integer, intent(inout) :: ha(iha,11)
      real,    intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      real,    intent(inout) :: b(n)
      integer, intent(out)   :: ifail
    end subroutine y12mae
  end interface

  integer, parameter :: NMAX  = 10
  integer, parameter :: NNP   = 200
  integer, parameter :: NN1P  = 200

  real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer :: n, z, nn, nn1, iha, ifail, nfail

  nfail = 0

  ! --- n=1: verify expected error IFAIL=12 (N<2 not supported by package) ---
  ! This exercises the error-diagnostic path of y12mae.
  n = 1
  call build_diag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12mae(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
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
    call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12mae(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    select case (n)
    case (2)
      call check_sp('y12ma_sp n=2 tridiag', n, b, ifail, 1.0e-4, nfail)
    case (3)
      call check_sp('y12ma_sp n=3 tridiag', n, b, ifail, 1.0e-4, nfail)
    case (4)
      call check_sp('y12ma_sp n=4 tridiag', n, b, ifail, 1.0e-4, nfail)
    end select
  end do

  ! --- n=5: 5x5 arrow matrix ---
  n = 5
  call build_arrow_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  nn = NNP ; nn1 = NN1P ; iha = NMAX
  call y12mae(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
  call check_sp('y12ma_sp n=5 arrow', n, b, ifail, 1.0e-4, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_sp tests PASSED'

contains

  ! Build a 1x1 diagonal system: A = [diag], b = [diag], solution x = [1].
  subroutine build_diag_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1
      rnr(z) = i ; snr(z) = i ; a(z) = real(i + 2)
    end do
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_diag_sp

  ! Build an n-by-n tridiagonal system (diag=3, off-diag=-1) with
  ! solution x=[1,...,1].  Entries stored in arbitrary order.
  subroutine build_tridiag_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) =  3.0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0
    end do
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_tridiag_sp

  ! Build an n-by-n arrow matrix with solution x=[1,...,1].
  ! Structure: A(1,j)=1 for j=1..n (full first row),
  !            A(i,1)=1 for i=2..n (full first col),
  !            A(i,i)=i+1 for i=2..n (diagonal blocks).
  subroutine build_arrow_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i, j
    z = 0
    ! First row: A(1,j)=1 for j=1..n
    do j = 1, n
      z = z + 1 ; rnr(z) = 1 ; snr(z) = j ; a(z) = 1.0
    end do
    ! Rows 2..n: sub-diagonal A(i,1)=1 and diagonal A(i,i)=i+1
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = 1   ; a(z) = 1.0
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) = real(i + 1)
    end do
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_arrow_sp

  ! Check that ifail=0 and max|b(1:n)-1| < tol; update nfail accordingly.
  subroutine check_sp(label, n, b, ifail, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, ifail
    real,             intent(in)    :: b(n), tol
    integer,          intent(inout) :: nfail
    real :: err
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ' ifail=', ifail
      nfail = nfail + 1
      return
    end if
    err = maxval(abs(b(1:n) - 1.0))
    if (err > tol) then
      write(*,'(3a,es10.3)') 'FAIL ', label, ' max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es10.3)') 'PASS ', label, ' max_err=', err
    end if
  end subroutine check_sp

end program test_y12ma_sp
