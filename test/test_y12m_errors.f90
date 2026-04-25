! SPDX-License-Identifier: GPL-2.0-only
!
! test_y12m_errors.f90 - test the error-diagnostic (IFAIL) paths of the
! y12m sparse-solver package.
!
! The documentation (doc, section 7) lists 25 distinct IFAIL values.
! Many of them are input-validation checks that can be triggered without
! performing any real arithmetic.  This test file covers:
!
!   y12mb error paths (checked in y12mbe/y12mbf before any elimination):
!     IFAIL=12  n < 2
!     IFAIL=13  z <= 0
!     IFAIL=5   nn < 2*z   (array too short)
!     IFAIL=6   nn1 < z    (array too short)
!     IFAIL=14  z < n      (structurally singular)
!     IFAIL=15  iha < n    (HA first dimension too small)
!     IFAIL=16  iflag(4) out of range (must be 0, 1, or 2)
!     IFAIL=17  empty row in the matrix
!     IFAIL=18  empty column in the matrix
!     IFAIL=24  column index out of range
!     IFAIL=25  row index out of range
!
!   y12mc error paths:
!     IFAIL=2   y12mb was not called first (iflag(1) != -1)
!     IFAIL=22  iflag(5) = 3 passed to y12mc (forbidden; only 1 or 2)
!
! Each test records a PASS or FAIL line and accumulates the failure count;
! the program exits with status 1 if any test failed.
!
program test_y12m_errors
  use y12m, only: y12mb, y12mc
  implicit none

  integer, parameter :: NNP  = 200
  integer, parameter :: NN1P = 200
  integer, parameter :: NMAX = 10

  integer :: nfail

  nfail = 0

  ! ---- y12mb / y12mbe error checks ----------------------------------------

  ! IFAIL=12: n < 2
  call test_mb_err('IFAIL=12 n<2', &
      n=1, z=1, nn=NNP, nn1=NN1P, iha=NMAX, &
      mode=0, want=12, nfail=nfail)

  ! IFAIL=13: z <= 0
  call test_mb_err('IFAIL=13 z<=0', &
      n=3, z=0, nn=NNP, nn1=NN1P, iha=NMAX, &
      mode=0, want=13, nfail=nfail)

  ! IFAIL=5: nn < 2*z  (use nn=5, z=4 so nn=5 < 2*4=8)
  call test_mb_err('IFAIL=5 nn<2z', &
      n=3, z=4, nn=5, nn1=NN1P, iha=NMAX, &
      mode=0, want=5, nfail=nfail)

  ! IFAIL=6: nn1 < z  (use nn1=3, z=4 so nn1=3 < z=4)
  call test_mb_err('IFAIL=6 nn1<z', &
      n=3, z=4, nn=NNP, nn1=3, iha=NMAX, &
      mode=0, want=6, nfail=nfail)

  ! IFAIL=14: z < n  (n=5, z=4 => z < n, with nn >= 2*z and nn1 >= z)
  call test_mb_err('IFAIL=14 z<n', &
      n=5, z=4, nn=NNP, nn1=NN1P, iha=NMAX, &
      mode=0, want=14, nfail=nfail)

  ! IFAIL=15: iha < n  (iha=3, n=5)
  call test_mb_err('IFAIL=15 iha<n', &
      n=5, z=5, nn=NNP, nn1=NN1P, iha=3, &
      mode=0, want=15, nfail=nfail)

  ! IFAIL=16: iflag(4) out of range (mode=3)
  call test_mb_err('IFAIL=16 mode=3', &
      n=3, z=3, nn=NNP, nn1=NN1P, iha=NMAX, &
      mode=3, want=16, nfail=nfail)

  ! IFAIL=17: empty row – n=4, z=4 elements in rows {1,2,4} only (row 3 missing)
  call test_mb_err17(nfail)

  ! IFAIL=18: empty column – n=4, z=4 elements cover all rows but not col 3
  call test_mb_err18(nfail)

  ! IFAIL=24: column index > n
  call test_mb_err24(nfail)

  ! IFAIL=25: row index > n
  call test_mb_err25(nfail)

  ! ---- y12mc / y12mce error checks ----------------------------------------

  ! IFAIL=2: y12mc called without preceding y12mb (iflag(1) not -1)
  call test_mc_no_mb(nfail)

  ! IFAIL=22: y12mc called with iflag(5)=3 (forbidden)
  call test_mc_iflag5_3(nfail)

  ! ---- Summary --------------------------------------------------------------
  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12m_errors tests PASSED'

contains

  ! -----------------------------------------------------------------------
  ! Generic y12mb error tester.
  ! Builds a valid (or deliberately invalid) matrix with the given sizes,
  ! calls y12mb with the specified parameters, and checks that IFAIL == want.
  ! The matrix elements themselves are irrelevant for these early-exit checks.
  ! -----------------------------------------------------------------------
  subroutine test_mb_err(label, n, z, nn, nn1, iha, mode, want, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, z, nn, nn1, iha, mode, want
    integer,          intent(inout) :: nfail

    real    :: a(NNP), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: i, nz

    ! Fill just enough of a/snr/rnr to not crash on VALID inputs.
    nz = max(0, min(z, NNP))
    do i = 1, nz
      a(i)   = real(i)
      snr(i) = mod(i-1, max(1, n)) + 1
      rnr(i) = mod(i-1, max(1, n)) + 1
    end do

    aflag(1) = 16.0  ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e16; aflag(4) = 1.0e-12

    iflag = 0
    iflag(4) = mode

    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)

    if (ifail == want) then
      write(*,'(3a)') 'PASS ', label, ': got expected IFAIL'
    else
      write(*,'(3a,i0,a,i0)') 'FAIL ', label, &
          ': expected IFAIL=', want, ' got ', ifail
      nfail = nfail + 1
    end if
  end subroutine test_mb_err

  ! -----------------------------------------------------------------------
  ! IFAIL=17: n=4, z=4, row 3 has no entries.
  !   Elements at (1,1), (2,2), (4,4), (1,2) => rows used: 1,1,2,4.
  !   Row 3 is missing => IFAIL=17.
  ! -----------------------------------------------------------------------
  subroutine test_mb_err17(nfail)
    integer, intent(inout) :: nfail

    real    :: a(NNP), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: n, z

    n = 4 ; z = 4
    a(1)=3.0 ; snr(1)=1 ; rnr(1)=1   ! (1,1)
    a(2)=3.0 ; snr(2)=2 ; rnr(2)=2   ! (2,2)
    a(3)=3.0 ; snr(3)=4 ; rnr(3)=4   ! (4,4) -- col 3 also empty; row 3 caught first
    a(4)=1.0 ; snr(4)=2 ; rnr(4)=1   ! (1,2)

    aflag(1) = 16.0  ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e16; aflag(4) = 1.0e-12
    iflag    = 0

    call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)

    if (ifail == 17) then
      write(*,'(a)') 'PASS IFAIL=17 empty row: got expected IFAIL'
    else
      write(*,'(a,i0,a,i0)') 'FAIL IFAIL=17 empty row: expected 17 got ', ifail
      nfail = nfail + 1
    end if
  end subroutine test_mb_err17

  ! -----------------------------------------------------------------------
  ! IFAIL=18: n=4, z=4, all rows have entries but column 3 has none.
  !   Elements at (1,1), (2,2), (3,4), (4,2) => cols used: 1,2,4,2.
  !   Column 3 is empty => IFAIL=18.
  ! -----------------------------------------------------------------------
  subroutine test_mb_err18(nfail)
    integer, intent(inout) :: nfail

    real    :: a(NNP), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: n, z

    n = 4 ; z = 4
    a(1)=3.0 ; snr(1)=1 ; rnr(1)=1   ! (1,1)
    a(2)=3.0 ; snr(2)=2 ; rnr(2)=2   ! (2,2)
    a(3)=1.0 ; snr(3)=4 ; rnr(3)=3   ! (3,4)
    a(4)=1.0 ; snr(4)=2 ; rnr(4)=4   ! (4,2)

    aflag(1) = 16.0  ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e16; aflag(4) = 1.0e-12
    iflag    = 0

    call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)

    if (ifail == 18) then
      write(*,'(a)') 'PASS IFAIL=18 empty column: got expected IFAIL'
    else
      write(*,'(a,i0)') 'FAIL IFAIL=18 empty column: expected 18 got ', ifail
      nfail = nfail + 1
    end if
  end subroutine test_mb_err18

  ! -----------------------------------------------------------------------
  ! IFAIL=24: column index out of range (snr(i) > n).
  !   n=3, z=3 diagonal matrix but snr(2)=5 > n=3.
  ! -----------------------------------------------------------------------
  subroutine test_mb_err24(nfail)
    integer, intent(inout) :: nfail

    real    :: a(NNP), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: n, z

    n = 3 ; z = 3
    a(1)=2.0 ; snr(1)=1 ; rnr(1)=1
    a(2)=2.0 ; snr(2)=5 ; rnr(2)=2   ! col 5 > n=3
    a(3)=2.0 ; snr(3)=3 ; rnr(3)=3

    aflag(1) = 16.0  ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e16; aflag(4) = 1.0e-12
    iflag    = 0

    call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)

    if (ifail == 24) then
      write(*,'(a)') 'PASS IFAIL=24 col>n: got expected IFAIL'
    else
      write(*,'(a,i0)') 'FAIL IFAIL=24 col>n: expected 24 got ', ifail
      nfail = nfail + 1
    end if
  end subroutine test_mb_err24

  ! -----------------------------------------------------------------------
  ! IFAIL=25: row index out of range (rnr(i) > n).
  !   n=3, z=3 diagonal matrix but rnr(2)=5 > n=3.
  ! -----------------------------------------------------------------------
  subroutine test_mb_err25(nfail)
    integer, intent(inout) :: nfail

    real    :: a(NNP), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: n, z

    n = 3 ; z = 3
    a(1)=2.0 ; snr(1)=1 ; rnr(1)=1
    a(2)=2.0 ; snr(2)=2 ; rnr(2)=5   ! row 5 > n=3
    a(3)=2.0 ; snr(3)=3 ; rnr(3)=3

    aflag(1) = 16.0  ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e16; aflag(4) = 1.0e-12
    iflag    = 0

    call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)

    if (ifail == 25) then
      write(*,'(a)') 'PASS IFAIL=25 row>n: got expected IFAIL'
    else
      write(*,'(a,i0)') 'FAIL IFAIL=25 row>n: expected 25 got ', ifail
      nfail = nfail + 1
    end if
  end subroutine test_mb_err25

  ! -----------------------------------------------------------------------
  ! IFAIL=2: y12mc called without first calling y12mb.
  ! iflag(1) is left at 0 (not set to -1 by y12mb), so y12mc detects the
  ! violation and returns IFAIL=2.
  ! -----------------------------------------------------------------------
  subroutine test_mc_no_mb(nfail)
    integer, intent(inout) :: nfail

    real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: n, z

    n = 3 ; z = 3
    a(1)=2.0 ; snr(1)=1 ; rnr(1)=1
    a(2)=2.0 ; snr(2)=2 ; rnr(2)=2
    a(3)=2.0 ; snr(3)=3 ; rnr(3)=3
    b(1:n) = 2.0

    aflag(1) = 16.0  ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e16; aflag(4) = 1.0e-12

    iflag    = 0   ! iflag(1)=0, NOT -1 (y12mb sets it to -1)
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 1

    call y12mc(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, NMAX, &
        aflag, iflag, ifail)

    if (ifail == 2) then
      write(*,'(a)') 'PASS IFAIL=2 mc without mb: got expected IFAIL'
    else
      write(*,'(a,i0)') 'FAIL IFAIL=2 mc without mb: expected 2 got ', ifail
      nfail = nfail + 1
    end if
  end subroutine test_mc_no_mb

  ! -----------------------------------------------------------------------
  ! IFAIL=22: y12mc called with iflag(5)=3.
  ! y12mc requires iflag(5) to be 1 or 2; value 3 is for y12md reuse.
  ! -----------------------------------------------------------------------
  subroutine test_mc_iflag5_3(nfail)
    integer, intent(inout) :: nfail

    real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: n, z

    n = 3 ; z = 3
    a(1)=2.0 ; snr(1)=1 ; rnr(1)=1
    a(2)=2.0 ; snr(2)=2 ; rnr(2)=2
    a(3)=2.0 ; snr(3)=3 ; rnr(3)=3
    b(1:n) = 2.0

    aflag(1) = 16.0  ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e16; aflag(4) = 1.0e-12

    ! First call y12mb so iflag(1) is set to -1.
    iflag    = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 1
    call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL IFAIL=22: y12mb unexpectedly failed with ', ifail
      nfail = nfail + 1
      return
    end if

    ! Now call y12mc with iflag(5)=3 (forbidden).
    iflag(5) = 3
    call y12mc(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, NMAX, &
        aflag, iflag, ifail)

    if (ifail == 22) then
      write(*,'(a)') 'PASS IFAIL=22 mc with iflag5=3: got expected IFAIL'
    else
      write(*,'(a,i0)') 'FAIL IFAIL=22 mc with iflag5=3: expected 22 got ', ifail
      nfail = nfail + 1
    end if
  end subroutine test_mc_iflag5_3

end program test_y12m_errors
