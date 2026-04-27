! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4-5
!
! test_y12mb_errors_dp.f90 – Error-branch coverage for the double-precision
! y12mb (y12mbf) routine.
!
! Each sub-test supplies a deliberately invalid argument to force a specific
! positive IFAIL return value, then verifies the expected code is returned.
! A final sub-test calls y12mb with a fully valid 3×3 diagonal matrix and
! checks:
!   • IFAIL  =  0  on exit
!   • IFLAG(1) = -1 on exit  (state-flag set by Y12MB on success)
!   • AFLAG(6) = maxval(abs(A(1:Z)))  (maximum element set by Y12MB)
!
! Error codes exercised and their triggering conditions:
!
!   ifail=12  n < 2                               (n=1)
!   ifail=13  z ≤ 0                               (z=0)
!   ifail= 5  nn < 2*z                            (nn=5 with z=3)
!   ifail= 6  nn1 < z                             (nn1=2 with z=3)
!   ifail=14  n > z  (structurally singular)      (n=5, z=3)
!   ifail=15  iha < n                             (iha=2, n=3)
!   ifail=16  iflag(4) ∉ {0,1,2}                 (iflag(4)=3 and iflag(4)=-1)
!   ifail=24  column index out of [1,n]           (snr(1)=n+1)
!   ifail=25  row index out of [1,n]              (rnr(1)=n+1)
!   ifail=17  all-zero row found                  (row 2 has no entries)
!   ifail=18  all-zero column found               (column 2 has no entries)
!   ifail=11  duplicate (i,j) position            (two entries at (2,2))
!
! Note – iflag(1) precondition not validated by y12mbf:
!   The API documentation states that iflag(1) must be ≥ 0 on entry to Y12MB,
!   but the source code (y12mbf.f) never checks this.  Passing iflag(1) = -2
!   (the post-Y12MC state) produces ifail = 0 and iflag(1) = -1, just as a
!   fresh call would.  This is a discrepancy between the documented precondition
!   and the actual implementation.
!
program test_y12mb_errors_dp
  use y12m, only: y12mb
  implicit none

  ! Array dimensions — large enough for all test cases plus headroom for the
  ! do-20 copy (a(z+1:2*z)) that y12mbf performs before the early-exit guard.
  integer, parameter :: NMAX   = 12
  integer, parameter :: NNMAX  = 60
  integer, parameter :: NN1MAX = 60

  double precision :: a(NNMAX), aflag(8)
  integer          :: snr(NNMAX), rnr(NN1MAX), ha(NMAX,11), iflag(10)
  integer          :: n, z, nn, nn1, iha, ifail, nfail
  double precision :: amax_expected

  nfail = 0

  ! =========================================================================
  ! ifail = 12 : n < 2
  ! =========================================================================
  ! Use a minimal 1×1 setup; every other parameter is valid.
  n = 1
  z = 1
  nn  = 2
  nn1 = 1
  iha = NMAX
  a(1) = 1.0d0 ; snr(1) = 1 ; rnr(1) = 1
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=12 n<2', ifail, 12, nfail)

  ! =========================================================================
  ! ifail = 13 : z ≤ 0
  ! =========================================================================
  n   = 2
  z   = 0
  nn  = 0
  nn1 = 0
  iha = NMAX
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=13 z<=0', ifail, 13, nfail)

  ! =========================================================================
  ! ifail = 5 : nn < 2*z
  ! =========================================================================
  n   = 3
  z   = 3
  nn  = 5        ! 5 < 2*3 = 6  →  ifail=5
  nn1 = 3
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=2.0d0 ; snr(2)=2 ; rnr(2)=2
  a(3)=3.0d0 ; snr(3)=3 ; rnr(3)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=5 nn<2*z', ifail, 5, nfail)

  ! =========================================================================
  ! ifail = 6 : nn1 < z
  ! =========================================================================
  n   = 3
  z   = 3
  nn  = 6
  nn1 = 2        ! 2 < 3  →  ifail=6
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=2.0d0 ; snr(2)=2 ; rnr(2)=2
  a(3)=3.0d0 ; snr(3)=3 ; rnr(3)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=6 nn1<z', ifail, 6, nfail)

  ! =========================================================================
  ! ifail = 14 : n > z  (structurally singular: fewer nonzeros than equations)
  ! =========================================================================
  ! All early checks pass (n>=2, z>0, nn>=2*z, nn1>=z), but n=5 > z=3.
  n   = 5
  z   = 3
  nn  = 6
  nn1 = 3
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=2.0d0 ; snr(2)=2 ; rnr(2)=2
  a(3)=3.0d0 ; snr(3)=3 ; rnr(3)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=14 n>z', ifail, 14, nfail)

  ! =========================================================================
  ! ifail = 15 : iha < n
  ! =========================================================================
  n   = 3
  z   = 3
  nn  = 6
  nn1 = 3
  iha = 2        ! 2 < n=3  →  ifail=15
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=2.0d0 ; snr(2)=2 ; rnr(2)=2
  a(3)=3.0d0 ; snr(3)=3 ; rnr(3)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=15 iha<n', ifail, 15, nfail)

  ! =========================================================================
  ! ifail = 16 : iflag(4) ∉ {0,1,2}  —  tested with iflag(4)=3 and iflag(4)=-1
  ! =========================================================================
  n   = 3
  z   = 3
  nn  = 6
  nn1 = 3
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=2.0d0 ; snr(2)=2 ; rnr(2)=2
  a(3)=3.0d0 ; snr(3)=3 ; rnr(3)=3
  call set_flags(iflag)
  iflag(4) = 3   ! > 2  →  ifail=16
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=16 iflag(4)=3', ifail, 16, nfail)

  call set_flags(iflag)
  iflag(4) = -1  ! < 0  →  ifail=16
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=16 iflag(4)=-1', ifail, 16, nfail)

  ! =========================================================================
  ! ifail = 24 : column index outside [1,n]
  ! =========================================================================
  ! snr(1) = n+1 = 4 is out of range.  iha=NMAX so that the internal access
  ! ha(n+1, 6) stays within the declared array bounds.
  n   = 3
  z   = 3
  nn  = 6
  nn1 = 3
  iha = NMAX        ! NMAX >= n+1 to avoid out-of-bounds inside y12mbf
  a(1)=1.0d0 ; snr(1)=n+1 ; rnr(1)=1   ! bad column
  a(2)=2.0d0 ; snr(2)=2   ; rnr(2)=2
  a(3)=3.0d0 ; snr(3)=3   ; rnr(3)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=24 bad col index', ifail, 24, nfail)

  ! =========================================================================
  ! ifail = 25 : row index outside [1,n]
  ! =========================================================================
  ! rnr(1) = n+1 = 4 is out of range; column indices are all valid.
  n   = 3
  z   = 3
  nn  = 6
  nn1 = 3
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1   ; rnr(1)=n+1   ! bad row
  a(2)=2.0d0 ; snr(2)=2   ; rnr(2)=2
  a(3)=3.0d0 ; snr(3)=3   ; rnr(3)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=25 bad row index', ifail, 25, nfail)

  ! =========================================================================
  ! ifail = 17 : all-zero row
  ! =========================================================================
  ! 3×3 matrix with 4 entries; row 2 has no entries.
  !   (1,1)=1  (1,2)=1  (3,2)=1  (3,3)=1
  ! z=4 > n=3, so ifail=14 is NOT triggered (that check requires n > z).
  n   = 3
  z   = 4
  nn  = 8
  nn1 = 4
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=1.0d0 ; snr(2)=2 ; rnr(2)=1
  a(3)=1.0d0 ; snr(3)=2 ; rnr(3)=3
  a(4)=1.0d0 ; snr(4)=3 ; rnr(4)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=17 empty row', ifail, 17, nfail)

  ! =========================================================================
  ! ifail = 18 : all-zero column
  ! =========================================================================
  ! 3×3 matrix with 3 entries; column 2 has no entries.
  !   (1,1)=1  (2,3)=1  (3,3)=1
  ! All rows non-empty; column 2 check fires at i=2 in the row/col scan loop.
  n   = 3
  z   = 3
  nn  = 6
  nn1 = 3
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=1.0d0 ; snr(2)=3 ; rnr(2)=2
  a(3)=1.0d0 ; snr(3)=3 ; rnr(3)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=18 empty col', ifail, 18, nfail)

  ! =========================================================================
  ! ifail = 11 : duplicate (i,j) position
  ! =========================================================================
  ! 3×3 matrix with 4 entries; two entries share position (2,2).
  !   (1,1)=5  (2,2)=4  (2,2)=4[dup]  (3,3)=3
  ! All rows and columns non-empty; the duplicate triggers ifail=11 in the
  ! column-linked-list construction loop.
  n   = 3
  z   = 4
  nn  = 8
  nn1 = 4
  iha = NMAX
  a(1)=5.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=4.0d0 ; snr(2)=2 ; rnr(2)=2
  a(3)=4.0d0 ; snr(3)=2 ; rnr(3)=2   ! duplicate (2,2)
  a(4)=3.0d0 ; snr(4)=3 ; rnr(4)=3
  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  call check_ifail('ifail=11 duplicate entry', ifail, 11, nfail)

  ! =========================================================================
  ! Valid call: verify ifail=0, iflag(1)=-1, and aflag(6)=maxval(abs(A(1:Z)))
  ! =========================================================================
  ! 3×3 diagonal matrix: A = diag(1, 2, 4).
  !   maxval(abs(A(1:3))) = 4.0
  n   = 3
  z   = 3
  nn  = 6
  nn1 = 3
  iha = NMAX
  a(1)=1.0d0 ; snr(1)=1 ; rnr(1)=1
  a(2)=2.0d0 ; snr(2)=2 ; rnr(2)=2
  a(3)=4.0d0 ; snr(3)=3 ; rnr(3)=3
  amax_expected = maxval(abs(a(1:z)))   ! = 4.0

  call set_flags(iflag)
  ha = 0
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)

  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL valid: expected ifail=0, got ', ifail
    nfail = nfail + 1
  else
    write(*,'(a)') 'PASS valid: ifail=0'
  end if

  if (iflag(1) /= -1) then
    write(*,'(a,i0)') 'FAIL valid: expected iflag(1)=-1 on exit, got ', iflag(1)
    nfail = nfail + 1
  else
    write(*,'(a)') 'PASS valid: iflag(1)=-1 on exit'
  end if

  if (aflag(6) /= amax_expected) then
    write(*,'(a,es12.5,a,es12.5)') &
        'FAIL valid: aflag(6)=', aflag(6), '  expected ', amax_expected
    nfail = nfail + 1
  else
    write(*,'(a,es12.5)') 'PASS valid: aflag(6)=maxval(abs(A))=', aflag(6)
  end if

  ! -------------------------------------------------------------------------
  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mb_errors_dp tests PASSED'

contains

  ! Initialise IFLAG and AFLAG to the same defaults used by the existing tests.
  subroutine set_flags(iflag)
    integer, intent(out) :: iflag(10)
    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 1
    iflag(6:10) = 0
  end subroutine set_flags

  ! Verify ifail equals the expected error code; increment nfail on mismatch.
  subroutine check_ifail(label, ifail, expected, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: ifail, expected
    integer,          intent(inout) :: nfail
    if (ifail == expected) then
      write(*,'(3a,i0)') 'PASS ', label, ' ifail=', ifail
    else
      write(*,'(3a,i0,a,i0)') 'FAIL ', label, &
          ' expected ifail=', expected, ' got ', ifail
      nfail = nfail + 1
    end if
  end subroutine check_ifail

end program test_y12mb_errors_dp
