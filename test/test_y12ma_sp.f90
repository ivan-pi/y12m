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
  call build_diag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
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
    call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
    nn = NNP ; nn1 = NN1P ; iha = NMAX
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
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
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b, ifail)
  call check_sp('y12ma_sp n=5 arrow', n, b, ifail, 1.0e-4, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_sp tests PASSED'

contains

  ! Build a 1x1 diagonal system: A=diag(i+2), b=row sums, solution x=[1].
  subroutine build_diag_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i ; a(z) = real(i + 2)
    end do
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_diag_sp

  ! Build an n-by-n tridiagonal (diag=3, off-diag=-1), solution x=[1,...,1].
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

  ! Build an n-by-n arrow matrix (full first row/col, diagonal rest),
  ! solution x=[1,...,1].
  subroutine build_arrow_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i, j
    z = 0
    do j = 1, n
      z = z + 1 ; rnr(z) = 1 ; snr(z) = j ; a(z) = 1.0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = 1   ; a(z) = 1.0
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) = real(i + 1)
    end do
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_arrow_sp

  ! Check IFAIL=0 and max|b(1:n)-1|<tol; increment nfail on failure.
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
