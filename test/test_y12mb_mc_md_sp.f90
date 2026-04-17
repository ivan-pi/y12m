! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12mb_mc_md_sp.f90 - Tests for the three-step single-precision API:
!   y12mbe (prepare) -> y12mce (factorize) -> y12mde (solve).
!
! Five sparse systems of sizes n=1..5 are solved, covering diagonal (n=1),
! tridiagonal (n=2,3,4), and arrow (n=5) matrix structures.
! Known solution x=[1,...,1] for all cases.
!
program test_y12mb_mc_md_sp
  implicit none

  interface
    subroutine y12mbe(n, z, a, snr, nn, rnr, nn1, &
        ha, iha, aflag, iflag, ifail)
      implicit none
      integer, intent(in)    :: n, z, nn, nn1, iha
      real,    intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      integer, intent(inout) :: ha(iha,11)
      real,    intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      integer, intent(out)   :: ifail
    end subroutine y12mbe

    subroutine y12mce(n, z, a, snr, nn, rnr, nn1, &
        pivot, b, ha, iha, aflag, iflag, ifail)
      implicit none
      integer, intent(in)    :: n, z, nn, nn1, iha
      real,    intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      real,    intent(inout) :: pivot(n)
      real,    intent(inout) :: b(n)
      integer, intent(inout) :: ha(iha,11)
      real,    intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      integer, intent(out)   :: ifail
    end subroutine y12mce

    subroutine y12mde(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      implicit none
      integer, intent(in)  :: n, nn, iha
      real,    intent(in)  :: a(nn)
      real,    intent(in)  :: pivot(n)
      integer, intent(in)  :: snr(nn)
      integer, intent(in)  :: ha(iha,11)
      integer, intent(in)  :: iflag(10)
      real,    intent(inout) :: b(n)
      integer, intent(out)   :: ifail
    end subroutine y12mde
  end interface

  integer, parameter :: NMAX  = 10
  integer, parameter :: NNP   = 200
  integer, parameter :: NN1P  = 200

  real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer :: n, z, ifail, nfail

  nfail = 0

  ! --- n=1: verify expected error IFAIL=12 (N<2 is not supported) ---
  ! This tests the error-diagnostic path of the three-step API.
  n = 1
  call build_diag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_sp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
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
  call build_full2_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_sp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_sp('y12mb_mc_md_sp n=2 full', n, b, ifail, 1.0e-4, nfail)

  ! --- n=3: tridiagonal ---
  n = 3
  call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_sp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_sp('y12mb_mc_md_sp n=3 tridiag', n, b, ifail, 1.0e-4, nfail)

  ! --- n=4: tridiagonal ---
  n = 4
  call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_sp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_sp('y12mb_mc_md_sp n=4 tridiag', n, b, ifail, 1.0e-4, nfail)

  ! --- n=5: arrow matrix ---
  n = 5
  call build_arrow_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  call solve3_sp(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, &
      NMAX, aflag, iflag, ifail)
  call check_sp('y12mb_mc_md_sp n=5 arrow', n, b, ifail, 1.0e-4, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mb_mc_md_sp tests PASSED'

contains

  ! Prepare-factorize-solve a single-precision sparse system.
  subroutine solve3_sp(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, &
      iha, aflag, iflag, ifail)
    integer, intent(in)    :: n, z, nn, nn1, iha
    real,    intent(inout) :: a(nn), b(n), pivot(n), aflag(8)
    integer, intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    integer, intent(out)   :: ifail

    ! IFLAG(1) must be >= 0 before first call of the package.
    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0  ! single system, no LU reuse
    iflag(5) = 1  ! discard L after factorization (saves memory)

    aflag(1) = 16.0     ! stability factor
    aflag(2) = 1.0e-12  ! drop tolerance
    aflag(3) = 1.0e+16  ! max growth factor
    aflag(4) = 1.0e-12  ! min pivot threshold

    call y12mbe(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) return

    call y12mce(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
        aflag, iflag, ifail)
    if (ifail /= 0) return

    call y12mde(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  end subroutine solve3_sp

  ! 1x1 diagonal system: A=[diag], b=[diag], solution x=[1].
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

  ! 2x2 full matrix: A = [[4,1],[2,5]], b = A*[1,1] = [5,7], x = [1,1].
  subroutine build_full2_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    ! Row 1: A(1,1)=4, A(1,2)=1
    z = z + 1 ; rnr(z) = 1 ; snr(z) = 1 ; a(z) = 4.0
    z = z + 1 ; rnr(z) = 1 ; snr(z) = 2 ; a(z) = 1.0
    ! Row 2: A(2,1)=2, A(2,2)=5
    z = z + 1 ; rnr(z) = 2 ; snr(z) = 1 ; a(z) = 2.0
    z = z + 1 ; rnr(z) = 2 ; snr(z) = 2 ; a(z) = 5.0
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_full2_sp

  ! n-by-n tridiagonal (diag=3, off-diag=-1), solution x=[1,...,1].
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

  ! n-by-n arrow matrix, solution x=[1,...,1].
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

  ! Check ifail=0 and max|b(1:n)-1| < tol; increment nfail on failure.
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

end program test_y12mb_mc_md_sp
