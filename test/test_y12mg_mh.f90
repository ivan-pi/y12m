! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12mg_mh.f90 - Tests for the analysis routines:
!   y12mhe  - compute the 1-norm of a single-precision sparse matrix
!   y12mhf  - compute the 1-norm of a double-precision sparse matrix
!   y12mge  - estimate the reciprocal condition number (single precision)
!
! The test uses a 5x5 tridiagonal matrix (diag=3, off-diag=-1):
!   - 1-norm = 5 (interior columns have column sum 1+3+1=5)
!   - After LU factorization the condition number estimate must be positive
!     and RCOND must be in (0, 1].
!
! A 10x10 tridiagonal matrix is used for the double-precision norm test.
!
program test_y12mg_mh
  implicit none

  interface
    ! Single-precision 1-norm
    subroutine y12mhe(n, nz, a, snr, work, anorm)
      implicit none
      integer, intent(in)    :: n, nz
      real,    intent(in)    :: a(nz)
      integer, intent(in)    :: snr(nz)
      real,    intent(inout) :: work(n)
      real,    intent(out)   :: anorm
    end subroutine y12mhe

    ! Double-precision 1-norm
    subroutine y12mhf(n, nz, a, snr, work, anorm)
      implicit none
      integer,          intent(in)    :: n, nz
      double precision, intent(in)    :: a(nz)
      integer,          intent(in)    :: snr(nz)
      double precision, intent(inout) :: work(n)
      double precision, intent(out)   :: anorm
    end subroutine y12mhf

    ! Prepare (single precision)
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

    ! Factorize (single precision)
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

    ! Solve (single precision)
    subroutine y12mde(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      implicit none
      integer, intent(in)    :: n, nn, iha
      real,    intent(in)    :: a(nn)
      real,    intent(in)    :: pivot(n)
      integer, intent(in)    :: snr(nn)
      integer, intent(in)    :: ha(iha,11)
      integer, intent(in)    :: iflag(10)
      real,    intent(inout) :: b(n)
      integer, intent(out)   :: ifail
    end subroutine y12mde

    ! Condition number estimate (single precision)
    subroutine y12mge(n, nn, a, snr, w, pivot, anorm, rcond, &
        iha, ha, iflag, ifail)
      implicit none
      integer, intent(in)    :: n, nn, iha
      real,    intent(in)    :: a(nn)
      integer, intent(in)    :: snr(nn)
      real,    intent(inout) :: w(n)
      real,    intent(in)    :: pivot(n)
      real,    intent(in)    :: anorm
      real,    intent(out)   :: rcond
      integer, intent(in)    :: ha(iha,11)
      integer, intent(in)    :: iflag(5)
      integer, intent(inout) :: ifail
    end subroutine y12mge
  end interface

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
  ! Test y12mhe: 1-norm of a 5x5 tridiagonal matrix.
  ! For a tridiagonal with diag=3 and off-diag=-1, the absolute column
  ! sums are:  col 1 = |3|+|-1| = 4 (boundary),
  !            col j (1<j<n) = |-1|+|3|+|-1| = 5 (interior),
  !            col n = |-1|+|3| = 4 (boundary).
  ! So the 1-norm = 5.
  ! =====================================================================
  n = 5
  call build_tridiag_sp(n, NNP, NN1P, z, a_sp, snr_sp, rnr_sp, b_sp)
  ! y12mhe works on the raw (column-indexed) non-zero array and SNR.
  ! At this point a_sp and snr_sp hold the original matrix in triplet order.
  call y12mhe(n, z, a_sp, snr_sp, work_sp, anorm_sp)
  if (abs(anorm_sp - 5.0) > 1.0e-5) then
    write(*,'(a,f10.5)') 'FAIL y12mhe n=5 expected anorm=5.0 got ', anorm_sp
    nfail = nfail + 1
  else
    write(*,'(a,f10.5)') 'PASS y12mhe n=5 anorm=', anorm_sp
  end if

  ! =====================================================================
  ! Test y12mhf: 1-norm of a 10x10 tridiagonal matrix (double precision).
  ! Same structure; 1-norm = 5.
  ! =====================================================================
  n = 10
  call build_tridiag_dp(n, NNP, z, a_dp, snr_dp)
  call y12mhf(n, z, a_dp, snr_dp, work_dp, anorm_dp)
  if (abs(anorm_dp - 5.0d0) > 1.0d-10) then
    write(*,'(a,f12.7)') 'FAIL y12mhf n=10 expected anorm=5.0 got ', anorm_dp
    nfail = nfail + 1
  else
    write(*,'(a,f12.7)') 'PASS y12mhf n=10 anorm=', anorm_dp
  end if

  ! =====================================================================
  ! Test y12mge: condition number estimate for a 5x5 tridiagonal matrix.
  ! Sequence: y12mhe -> y12mbe -> y12mce -> y12mde -> y12mge.
  ! y12mge requires iflag(5) /= 1 (L factors must be preserved).
  ! =====================================================================
  n = 5
  call build_tridiag_sp(n, NNP, NN1P, z, a_sp, snr_sp, rnr_sp, b_sp)
  nn = NNP ; nn1 = NN1P ; iha = NMAX

  ! Compute 1-norm before factorization modifies a_sp.
  call y12mhe(n, z, a_sp, snr_sp, work_sp, anorm_sp)

  ! Setup flags: iflag(5)=2 keeps L factors so y12mge can be called.
  iflag_sp(1) = 0
  iflag_sp(2) = 3
  iflag_sp(3) = 1
  iflag_sp(4) = 0   ! single system
  iflag_sp(5) = 2   ! KEEP L factors (required for y12mge)

  aflag_sp(1) = 16.0
  aflag_sp(2) = 1.0e-12
  aflag_sp(3) = 1.0e+16
  aflag_sp(4) = 1.0e-12

  call y12mbe(n, z, a_sp, snr_sp, nn, rnr_sp, nn1, ha_sp, iha, &
      aflag_sp, iflag_sp, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL y12mge: y12mbe ifail=', ifail
    nfail = nfail + 1
  else

    call y12mce(n, z, a_sp, snr_sp, nn, rnr_sp, nn1, pivot_sp, b_sp, &
        ha_sp, iha, aflag_sp, iflag_sp, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL y12mge: y12mce ifail=', ifail
      nfail = nfail + 1
    else

      call y12mde(n, a_sp, nn, b_sp, pivot_sp, snr_sp, ha_sp, iha, &
          iflag_sp, ifail)
      if (ifail /= 0) then
        write(*,'(a,i0)') 'FAIL y12mge: y12mde ifail=', ifail
        nfail = nfail + 1
      else

        ! Verify solution before computing condition number
        call check_sp('y12mg_mh y12mde n=5', n, b_sp, ifail, 1.0e-4, nfail)

        ! Call y12mge with ifail=0 on entry (required by implementation)
        ifail = 0
        call y12mge(n, nn, a_sp, snr_sp, work_sp, pivot_sp, anorm_sp, &
            rcond_sp, iha, ha_sp, iflag_sp, ifail)
        if (ifail /= 0) then
          write(*,'(a,i0)') 'FAIL y12mge ifail=', ifail
          nfail = nfail + 1
        else if (rcond_sp <= 0.0 .or. rcond_sp > 1.0) then
          write(*,'(a,es10.3)') 'FAIL y12mge rcond out of (0,1]: ', rcond_sp
          nfail = nfail + 1
        else
          write(*,'(a,es10.3,a,es10.3)') &
              'PASS y12mge n=5 rcond=', rcond_sp, &
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
  ! Only a and snr are needed for y12mhf.
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

end program test_y12mg_mh
