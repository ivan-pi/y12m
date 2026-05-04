! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4-5
!
! test_y12mb_errors.f90 - Error-branch coverage for y12mb.
!
! Each sub-test supplies a deliberately invalid argument to force a specific
! positive IFAIL return value, then verifies the expected code is returned.
! A final sub-test calls y12mb with a fully valid 3x3 diagonal matrix and
! checks:
!   * IFAIL    =  0  on exit
!   * IFLAG(1) = -1  on exit  (state-flag set by Y12MB on success)
!   * AFLAG(6) = maxval(abs(A(1:Z)))  (maximum element set by Y12MB)
!
! Error codes exercised and their triggering conditions:
!
!   ifail=12  n < 2                               (n=1)
!   ifail=13  z <= 0                              (z=0)
!   ifail= 5  nn < 2*z                            (nn=5 with z=3)
!   ifail= 6  nn1 < z                             (nn1=2 with z=3)
!   ifail=14  n > z  (structurally singular)      (n=5, z=3)
!   ifail=15  iha < n                             (iha=2, n=3)
!   ifail=16  iflag(4) not in {0,1,2}             (iflag(4)=3 and iflag(4)=-1)
!   ifail=24  column index out of [1,n]           (snr(1)=n+1)
!   ifail=25  row index out of [1,n]              (rnr(1)=n+1)
!   ifail=17  all-zero row found                  (row 2 has no entries)
!   ifail=18  all-zero column found               (column 2 has no entries)
!   ifail=11  duplicate (i,j) position            (two entries at (2,2))
!
! Note - iflag(1) precondition not validated by y12mbf:
!   The API documentation states that iflag(1) must be >= 0 on entry to Y12MB,
!   but the source code (y12mbf.f) never checks this.  Passing iflag(1) = -2
!   (the post-Y12MC state) produces ifail = 0 and iflag(1) = -1, just as a
!   fresh call would.  This is a discrepancy between the documented precondition
!   and the actual implementation.
!
program test_y12mb_errors
  use y12m, only: y12mb
  use y12m_test_helpers, only: set_flags, check_ifail
  implicit none

#ifdef TEST_SINGLE_PRECISION
  integer, parameter :: dp = kind(1.0e0)
#else
  integer, parameter :: dp = kind(1.0d0)
#endif

  ! Array capacities - y12mbf copies a(1:z) to a(z+1:2*z) and snr(1:z) to
  ! snr(z+1:2*z) while scanning column/row indices; nn >= 2*z is therefore
  ! required even for tests that expect an early-error return.
  integer, parameter :: NMAX = 12, NNMAX = 60, NN1MAX = 60

  real(dp) :: a(NNMAX), aflag(8)
  integer :: snr(NNMAX), rnr(NN1MAX), ha(NMAX,11), iflag(10)
  integer :: ifail, nfail

  nfail = 0

  ! =========================================================================
  ! ifail = 12 : n < 2
  ! =========================================================================
  block
    integer :: n = 1, z = 1, nn = 2, nn1 = 1, iha = NMAX
    a(1:z) = [11]
    snr(1:z) = [1]
    rnr(1:z) = [1]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=12 n<2', ifail, 12, nfail)
  end block

  ! =========================================================================
  ! ifail = 13 : z <= 0
  ! =========================================================================
  block
    integer :: n = 2, z = 0, nn = 0, nn1 = 0, iha = NMAX
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=13 z<=0', ifail, 13, nfail)
  end block

  ! =========================================================================
  ! ifail = 5 : nn < 2*z
  ! =========================================================================
  block
    integer :: n = 3, z = 3, nn = 5, nn1 = 3, iha = NMAX  ! 5 < 2*3=6
    a(1:z) = [11, 22, 33]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=5 nn<2*z', ifail, 5, nfail)
  end block

  ! =========================================================================
  ! ifail = 6 : nn1 < z
  ! =========================================================================
  block
    integer :: n = 3, z = 3, nn = 6, nn1 = 2, iha = NMAX  ! 2 < 3
    a(1:z) = [11, 22, 33]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=6 nn1<z', ifail, 6, nfail)
  end block

  ! =========================================================================
  ! ifail = 14 : n > z  (structurally singular: fewer nonzeros than equations)
  ! =========================================================================
  ! All early checks pass (n>=2, z>0, nn>=2*z, nn1>=z), but n=5 > z=3.
  block
    integer :: n = 5, z = 3, nn = 6, nn1 = 3, iha = NMAX
    a(1:z) = [11, 22, 33]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=14 n>z', ifail, 14, nfail)
  end block

  ! =========================================================================
  ! ifail = 15 : iha < n
  ! =========================================================================
  block
    integer :: n = 3, z = 3, nn = 6, nn1 = 3, iha = 2  ! 2 < n=3
    a(1:z) = [11, 22, 33]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=15 iha<n', ifail, 15, nfail)
  end block

  ! =========================================================================
  ! ifail = 16 : iflag(4) not in {0,1,2}  -  tested with iflag(4)=3 and iflag(4)=-1
  ! =========================================================================
  block
    integer :: n = 3, z = 3, nn = 6, nn1 = 3, iha = NMAX
    a(1:z) = [11, 22, 33]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    call set_flags(iflag)
    iflag(4) = 3   ! > 2  -> ifail=16
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=16 iflag(4)=3', ifail, 16, nfail)
    call set_flags(iflag)
    iflag(4) = -1  ! < 0  -> ifail=16
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=16 iflag(4)=-1', ifail, 16, nfail)
  end block

  ! =========================================================================
  ! ifail = 24 : column index outside [1,n]
  ! =========================================================================
  ! iha=NMAX ensures ha(n+1,6) inside y12mbf stays within declared bounds.
  block
    integer :: n = 3, z = 3, nn = 6, nn1 = 3, iha = NMAX
    a(1:z) = [11, 22, 33]
    snr(1:z) = [n+1, 2, 3]   ! snr(1) = 4 is out of [1,n]
    rnr(1:z) = [1, 2, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=24 bad col index', ifail, 24, nfail)
  end block

  ! =========================================================================
  ! ifail = 25 : row index outside [1,n]
  ! =========================================================================
  block
    integer :: n = 3, z = 3, nn = 6, nn1 = 3, iha = NMAX
    a(1:z) = [11, 22, 33]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [n+1, 2, 3]   ! rnr(1) = 4 is out of [1,n]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=25 bad row index', ifail, 25, nfail)
  end block

  ! =========================================================================
  ! ifail = 17 : all-zero row
  ! =========================================================================
  ! 3x3 matrix with 4 entries; row 2 has no entries.
  !   (1,1), (1,2), (3,2), (3,3) - row 2 is empty
  ! z=4 > n=3, so ifail=14 is NOT triggered (that check requires n > z).
  block
    integer :: n = 3, z = 4, nn = 8, nn1 = 4, iha = NMAX
    a(1:z) = [11, 11, 11, 11]
    snr(1:z) = [1, 2, 2, 3]
    rnr(1:z) = [1, 1, 3, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=17 empty row', ifail, 17, nfail)
  end block

  ! =========================================================================
  ! ifail = 18 : all-zero column
  ! =========================================================================
  ! 3x3 matrix with 3 entries; column 2 has no entries.
  !   (1,1), (2,3), (3,3) - column 2 is empty
  ! All rows non-empty; column 2 check fires in the row/col scan loop.
  block
    integer :: n = 3, z = 3, nn = 6, nn1 = 3, iha = NMAX
    a(1:z) = [11, 11, 11]
    snr(1:z) = [1, 3, 3]
    rnr(1:z) = [1, 2, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=18 empty col', ifail, 18, nfail)
  end block

  ! =========================================================================
  ! ifail = 11 : duplicate (i,j) position
  ! =========================================================================
  ! 3x3 matrix with 4 entries; two entries share position (2,2).
  !   (1,1), (2,2), (2,2)[dup], (3,3)
  ! All rows and columns non-empty; the duplicate triggers ifail=11 in the
  ! column-linked-list construction loop.
  block
    integer :: n = 3, z = 4, nn = 8, nn1 = 4, iha = NMAX
    a(1:z) = [55, 44, 44, 33]
    snr(1:z) = [1, 2, 2, 3]
    rnr(1:z) = [1, 2, 2, 3]
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=11 duplicate entry', ifail, 11, nfail)
  end block

  ! =========================================================================
  ! Valid call: verify ifail=0, iflag(1)=-1, and aflag(6)=max|A(1:z)|
  ! =========================================================================
  ! 3x3 diagonal matrix diag(11, 22, 44); maxval(abs(A(1:3))) = 44.
  block
    integer :: n = 3, z = 3, nn = 6, nn1 = 3, iha = NMAX
    real(dp) :: amax_expected
    a(1:z) = [11, 22, 44]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    amax_expected = maxval(abs(a(1:z)))
    call set_flags(iflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL valid: expected ifail=0, got ', ifail
      nfail = nfail + 1
    end if
    if (iflag(1) /= -1) then
      write(*,'(a,i0)') 'FAIL valid: expected iflag(1)=-1, got ', iflag(1)
      nfail = nfail + 1
    end if
    if (aflag(6) /= amax_expected .or. &
        maxval(abs(a(1:z))) /= amax_expected) then
      write(*,'(a,es12.5,a,es12.5)') &
          'FAIL valid: aflag(6)=', aflag(6), ' expected ', amax_expected
      nfail = nfail + 1
    end if
  end block

  ! =========================================================================
  ! TODO: Verify output arrays after a successful Y12MB call
  ! =========================================================================
  ! After a valid call, also check:
  !   snr(1:z)        - column indices, reordered by row
  !   rnr(1:z)        - row numbers ordered by column (column-linked list)
  !   ha(i,1), ha(i,3) - start/end positions of row i in a/snr
  !   ha(i,4), ha(i,6) - start/end positions of column i in rnr
  !   ha columns 7, 8, 11 - pivotal search metadata used by y12mc
  !
  ! These checks should be added once the expected output layout for a small
  ! reference matrix has been worked out.
  ! =========================================================================

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mb_errors tests PASSED'

end program test_y12mb_errors
