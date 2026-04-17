! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12mg_mh.f90 - Tests for y12mh (1-norm) and y12mg (condition number).
!
! This program tests three analysis routines:
!
!   y12mh (single precision, resolves to y12mhe):
!     Computes the 1-norm (max absolute column sum) of a 5x5 tridiagonal
!     matrix (diag=3, off-diag=-1).  Expected result: 5.0
!     (interior columns have absolute column sum |−1|+|3|+|−1| = 5).
!
!   y12mh (double precision, resolves to y12mhf):
!     Same test for a 10x10 tridiagonal matrix.  Expected result: 5.0.
!
!   y12mg (single precision, resolves to y12mge):
!     Estimates the reciprocal condition number (RCOND) of a 5x5 tridiagonal
!     matrix after LU factorization.  The test sequence is:
!       y12mh -> y12mb -> y12mc -> y12md -> y12mg
!     y12mg requires that IFLAG(5)/=1 (L factors must be preserved).
!     Pass criterion: RCOND in (0, 1].
!
program test_y12mg_mh
  use y12m
  implicit none

  integer, parameter :: NMAX  = 12
  integer, parameter :: NNP   = 400
  integer, parameter :: NN1P  = 400

  ! Single-precision working arrays
  real    :: a_sp(NNP), pivot_sp(NMAX), b_sp(NMAX), aflag_sp(8)
  real    :: work_sp(NMAX), anorm_sp, rcond_sp
  integer :: snr_sp(NNP), rnr_sp(NN1P), ha_sp(NMAX,11), iflag_sp(10)

  ! Double-precision working arrays
  double precision :: a_dp(NNP), work_dp(NMAX), anorm_dp
  integer          :: snr_dp(NNP)

  integer :: n, z, nn, nn1, iha, ifail, nfail

  nfail = 0

  ! =====================================================================
  ! Test y12mh (single): 1-norm of a 5x5 tridiagonal matrix.
  ! Absolute column sums: boundary cols = |3|+|-1| = 4,
  !                       interior cols = |-1|+|3|+|-1| = 5.
  ! 1-norm = max column sum = 5.
  ! =====================================================================
  n = 5
  call build_tridiag_sp(n, NNP, NN1P, z, a_sp, snr_sp, rnr_sp, b_sp)
  call y12mh(n, z, a_sp, snr_sp, work_sp, anorm_sp)
  if (abs(anorm_sp - 5.0) > 1.0e-5) then
    write(*,'(a,f10.5)') 'FAIL y12mhe n=5 expected anorm=5.0 got ', anorm_sp
    nfail = nfail + 1
  else
    write(*,'(a,f10.5)') 'PASS y12mhe n=5 anorm=', anorm_sp
  end if

  ! =====================================================================
  ! Test y12mh (double): 1-norm of a 10x10 tridiagonal matrix.
  ! Same structure; 1-norm = 5.
  ! =====================================================================
  n = 10
  call build_tridiag_dp(n, NNP, z, a_dp, snr_dp)
  call y12mh(n, z, a_dp, snr_dp, work_dp, anorm_dp)
  if (abs(anorm_dp - 5.0d0) > 1.0d-10) then
    write(*,'(a,f12.7)') 'FAIL y12mhf n=10 expected anorm=5.0 got ', anorm_dp
    nfail = nfail + 1
  else
    write(*,'(a,f12.7)') 'PASS y12mhf n=10 anorm=', anorm_dp
  end if

  ! =====================================================================
  ! Test y12mg: condition number estimate for a 5x5 tridiagonal matrix.
  ! Sequence: y12mh -> y12mb -> y12mc -> y12md -> y12mg.
  ! IFLAG(5)=2 is required to preserve L factors for y12mg.
  ! =====================================================================
  n = 5
  call build_tridiag_sp(n, NNP, NN1P, z, a_sp, snr_sp, rnr_sp, b_sp)
  nn = NNP ; nn1 = NN1P ; iha = NMAX

  ! Compute 1-norm before factorization modifies a_sp.
  call y12mh(n, z, a_sp, snr_sp, work_sp, anorm_sp)

  ! IFLAG(5)=2: keep L factors (required by y12mg).
  iflag_sp(1) = 0
  iflag_sp(2) = 3
  iflag_sp(3) = 1
  iflag_sp(4) = 0
  iflag_sp(5) = 2

  aflag_sp(1) = 16.0
  aflag_sp(2) = 1.0e-12
  aflag_sp(3) = 1.0e+16
  aflag_sp(4) = 1.0e-12

  call y12mb(n, z, a_sp, snr_sp, nn, rnr_sp, nn1, ha_sp, iha, &
      aflag_sp, iflag_sp, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL y12mg: y12mb ifail=', ifail
    nfail = nfail + 1
  else

    call y12mc(n, z, a_sp, snr_sp, nn, rnr_sp, nn1, pivot_sp, b_sp, &
        ha_sp, iha, aflag_sp, iflag_sp, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL y12mg: y12mc ifail=', ifail
      nfail = nfail + 1
    else

      call y12md(n, a_sp, nn, b_sp, pivot_sp, snr_sp, ha_sp, iha, &
          iflag_sp, ifail)
      if (ifail /= 0) then
        write(*,'(a,i0)') 'FAIL y12mg: y12md ifail=', ifail
        nfail = nfail + 1
      else

        ! Verify the triangular solve before estimating the condition number.
        call check_sp('y12mg_mh y12md n=5', n, b_sp, ifail, 1.0e-4, nfail)

        ! y12mg requires IFAIL=0 on entry.
        ifail = 0
        call y12mg(n, nn, a_sp, snr_sp, work_sp, pivot_sp, anorm_sp, &
            rcond_sp, iha, ha_sp, iflag_sp, ifail)
        if (ifail /= 0) then
          write(*,'(a,i0)') 'FAIL y12mg ifail=', ifail
          nfail = nfail + 1
        else if (rcond_sp <= 0.0 .or. rcond_sp > 1.0) then
          write(*,'(a,es10.3)') 'FAIL y12mg rcond out of (0,1]: ', rcond_sp
          nfail = nfail + 1
        else
          write(*,'(a,es10.3,a,es10.3)') &
              'PASS y12mg n=5 rcond=', rcond_sp, &
              ' cond=', 1.0 / rcond_sp
        end if

      end if
    end if
  end if

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mg_mh tests PASSED'

contains

  ! n-by-n tridiagonal (diag=3, off-diag=-1), single precision, x=[1,...,1].
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

  ! n-by-n tridiagonal (diag=3, off-diag=-1), double precision.
  ! Only a and snr are filled; row numbers are not needed for y12mh.
  subroutine build_tridiag_dp(n, nnmax, z, a, snr)
    integer,          intent(in)  :: n, nnmax
    integer,          intent(out) :: z
    double precision, intent(out) :: a(nnmax)
    integer,          intent(out) :: snr(nnmax)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; snr(z) = i   ; a(z) =  3.0d0
    end do
    do i = 2, n
      z = z + 1 ; snr(z) = i-1 ; a(z) = -1.0d0
    end do
    do i = 1, n-1
      z = z + 1 ; snr(z) = i+1 ; a(z) = -1.0d0
    end do
  end subroutine build_tridiag_dp

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

end program test_y12mg_mh
