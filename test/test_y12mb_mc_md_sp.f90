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

  ! Execute the three-step prepare->factorize->solve sequence (single precision).
  subroutine solve3_sp(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, &
      iha, aflag, iflag, ifail)
    integer, intent(in)    :: n, z, nn, nn1, iha
    real,    intent(inout) :: a(nn), b(n), pivot(n), aflag(8)
    integer, intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    integer, intent(out)   :: ifail

    ! IFLAG(1) must be >= 0 before the first call of the package.
    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0  ! single system; no LU reuse
    iflag(5) = 1  ! discard L after factorization (saves memory)

    aflag(1) = 16.0     ! stability factor
    aflag(2) = 1.0e-12  ! drop tolerance
    aflag(3) = 1.0e+16  ! max growth factor
    aflag(4) = 1.0e-12  ! min pivot threshold

    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) return

    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
        aflag, iflag, ifail)
    if (ifail /= 0) return

    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  end subroutine solve3_sp

  ! 1x1 diagonal system: A=diag(i+2), b=row sums, solution x=[1].
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

  ! 2x2 full matrix: A=[[4,1],[2,5]], b=A*[1,1]=[5,7], solution x=[1,1].
  subroutine build_full2_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    z = z + 1 ; rnr(z) = 1 ; snr(z) = 1 ; a(z) = 4.0
    z = z + 1 ; rnr(z) = 1 ; snr(z) = 2 ; a(z) = 1.0
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

  ! n-by-n arrow matrix (full first row/col, diagonal rest), x=[1,...,1].
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

end program test_y12mb_mc_md_sp
