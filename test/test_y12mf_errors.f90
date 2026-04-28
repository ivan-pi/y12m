! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4-6
!
! test_y12mf_errors.f90 - Error-branch and output-parameter coverage for y12mf.
!
! Use test_y12md_errors.F90 as template.
!
! y12mf (implemented by y12mfe) is a black-box driver that wraps y12mb, y12mc,
! and y12md with iterative refinement.  Only a single-precision variant exists.
!
! DISCREPANCIES NOTED:
!
!   1. The problem statement gives b = (10, 11, 45, 33, -22, 31) with claimed
!      exact solution x = (1, 2, 3, 4, 5, 6).  However, for the 6x6 reference
!      matrix, A(4,:)*[1,2,3,4,5,6] = -2*1 + 20*4 - 2*5 = 68, not 33.
!      This test uses the corrected b = (10, 11, 45, 68, -22, 31).
!
!   2. The problem statement specifies A(2,2) = 13.  The reference matrix
!      used consistently throughout the test suite (test_y12mc_errors.F90,
!      test_y12md_errors.F90) has A(2,2) = 12; that value is used here.
!
! Error codes exercised and their triggering conditions:
!
!   ifail=10  iflag(5) = 1       (y12mf requires iflag(5) = 2 or 3)
!   ifail=23  iflag(11) = 1      (minimum iteration count is 2)
!   ifail=12  n = 1              (propagated from internal y12mb call)
!   ifail=13  nz = 0             (propagated from internal y12mb call)
!   ifail= 5  nn < 2*nz          (propagated from internal y12mb call)
!   ifail=15  iha < n            (propagated from internal y12mb call)
!
! Success postconditions verified (using the 6x6 reference matrix with
! AFLAG(1:4) = [128.0, 1e-3, 1e16, 1e-12], IFLAG(2:5) = [2,1,1,2],
! IFLAG(11) = 25):
!
!   * IFAIL     = 0
!   * IFLAG(12) in [1, IFLAG(11)]   actual refinement iterations
!   * AFLAG(9)  >= 0.0              max-norm of last correction vector d(p-1)
!   * AFLAG(10) >= 0.0              max-norm of last residual vector r(p-1)
!   * AFLAG(11) >= 0.0              max-norm of corrected solution x(p)
!   * AFLAG(5)  >= 1.0              growth factor (set by internal y12mc call)
!   * AFLAG(8)  > 0.0               smallest pivot magnitude
!   * Y(i) /= 0.0 for all i         pivot elements (diagonal of U)
!   * min|Y(1:n)| == AFLAG(8)       consistency between Y and AFLAG(8)
!   * B1 == original RHS            y12mf saves original b in B1 on exit
!   * X close to (1, 2, 3, 4, 5, 6) accurate solution
!
program test_y12mf_errors
  use y12m, only: y12mf
  implicit none

  integer, parameter :: NMAX   = 10
  integer, parameter :: NNMAX  = 200
  integer, parameter :: NN1MAX = 100
  integer, parameter :: NZMAX  = 50    ! capacity for a1/sn (>= nz of reference)
  integer, parameter :: N_REF  = 6
  integer, parameter :: Z_REF  = 15

  real    :: a(NNMAX), a1(NZMAX), aflag(11), b(NMAX), b1(NMAX), x(NMAX), y_piv(NMAX)
  integer :: snr(NNMAX), rnr(NN1MAX), sn(NZMAX), ha(NMAX,13), iflag(12)
  integer :: ifail, nfail

  nfail = 0

  ! =========================================================================
  ! ifail = 10 : iflag(5) = 1  (y12mf requires iflag(5) = 2 or 3)
  !
  ! y12mfe immediately returns ifail=10 when state = iflag(5) = 1 before any
  ! array modifications or delegation to y12mb/y12mc/y12md.
  ! =========================================================================
  blk10: block
    integer :: n = N_REF, nz = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    iflag(5) = 1   ! invalid for y12mf: must be 2 or 3
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=10 iflag(5)=1', ifail, 10, nfail)
  end block blk10

  ! =========================================================================
  ! ifail = 23 : iflag(11) < 2  (minimum iteration count is 2)
  !
  ! y12mfe checks it = iflag(11) < 2 immediately after the state check;
  ! returns ifail=23 without modifying any work array.
  ! =========================================================================
  blk23: block
    integer :: n = N_REF, nz = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    iflag(11) = 1   ! < 2 -> ifail=23
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=23 iflag(11)<2', ifail, 23, nfail)
  end block blk23

  ! =========================================================================
  ! ifail = 12 : n < 2  (propagated from internal y12mb call)
  ! =========================================================================
  blk12: block
    integer :: n = 1, nz = 1, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    a(1) = 5.0
    snr(1) = 1
    rnr(1) = 1
    b(1)   = 5.0
    call set_flags(iflag)
    call set_aflag(aflag)
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=12 n<2', ifail, 12, nfail)
  end block blk12

  ! =========================================================================
  ! ifail = 13 : nz <= 0  (propagated from internal y12mb call)
  ! =========================================================================
  blk13: block
    integer :: n = N_REF, nz = 0, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=13 nz<=0', ifail, 13, nfail)
  end block blk13

  ! =========================================================================
  ! ifail = 5 : nn < 2*nz  (propagated from internal y12mb call)
  ! =========================================================================
  blk5: block
    integer :: n = N_REF, nz = Z_REF, nn = 2*Z_REF - 1, nn1 = NN1MAX, iha = NMAX

    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=5 nn<2*nz', ifail, 5, nfail)
  end block blk5

  ! =========================================================================
  ! ifail = 15 : iha < n  (propagated from internal y12mb call)
  ! =========================================================================
  blk15: block
    integer :: n = N_REF, nz = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = N_REF - 1

    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=15 iha<n', ifail, 15, nfail)
  end block blk15

  ! =========================================================================
  ! Success: y12mf with 6x6 reference matrix
  !
  ! Settings from problem statement:
  !   AFLAG(1:4) = [128.0, 1.e-3, 1.e+16, 1.e-12]
  !   IFLAG(2:5) = [2, 1, 1, 2]
  !   IFLAG(11)  = 25
  !
  ! Expected solution: x = (1, 2, 3, 4, 5, 6).
  ! =========================================================================
  success: block
    integer :: n = N_REF, nz = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    integer :: i, max_it
    real    :: b_orig(NMAX), min_piv, err_sol, err_b1
    logical :: any_fail

    call load_ref6(a, snr, rnr, b)
    b_orig(1:N_REF) = b(1:N_REF)   ! save original RHS for postcondition check

    aflag(1)    = 128.0e0
    aflag(2)    = 1.0e-3
    aflag(3)    = 1.0e+16
    aflag(4)    = 1.0e-12
    aflag(5:11) = 0.0

    iflag(1)    = 0
    iflag(2)    = 2
    iflag(3)    = 1
    iflag(4)    = 1
    iflag(5)    = 2   ! keep L, required by y12mf
    iflag(6:12) = 0
    max_it      = 25
    iflag(11)   = max_it

    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)

    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL success: y12mf returned ifail=', ifail
      nfail = nfail + 1
      exit success
    end if

    any_fail = .false.

    ! IFLAG(12): actual iteration count must be in [1, max_it]
    if (iflag(12) < 1 .or. iflag(12) > max_it) then
      write(*,'(a,i0,a,i0,a)') 'FAIL success: IFLAG(12)=', iflag(12), &
          ' not in [1,', max_it, ']'
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! AFLAG(9) >= 0: max-norm of last correction vector d(p-1)
    if (aflag(9) < 0.0) then
      write(*,'(a,1pg12.5)') 'FAIL success: AFLAG(9) < 0, got ', aflag(9)
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! AFLAG(10) >= 0: max-norm of last residual vector r(p-1)
    if (aflag(10) < 0.0) then
      write(*,'(a,1pg12.5)') 'FAIL success: AFLAG(10) < 0, got ', aflag(10)
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! AFLAG(11) >= 0: max-norm of corrected solution x(p)
    if (aflag(11) < 0.0) then
      write(*,'(a,1pg12.5)') 'FAIL success: AFLAG(11) < 0, got ', aflag(11)
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! AFLAG(5) >= 1.0: growth factor (set by the internal y12mc call)
    if (aflag(5) < 1.0) then
      write(*,'(a,1pg12.5)') 'FAIL success: AFLAG(5) < 1.0, got ', aflag(5)
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! AFLAG(8) > 0: smallest pivot magnitude (set during factorization)
    if (aflag(8) <= 0.0) then
      write(*,'(a,1pg12.5)') 'FAIL success: AFLAG(8) <= 0, got ', aflag(8)
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! Y_PIV(i) /= 0 for all i: pivot elements must all be nonzero
    do i = 1, n
      if (y_piv(i) == 0.0) then
        write(*,'(a,i0)') 'FAIL success: Y_PIV(i) == 0.0 at i=', i
        nfail = nfail + 1
        any_fail = .true.
      end if
    end do

    ! min|Y_PIV(1:n)| must equal AFLAG(8) (pivot array is consistent with aflag)
    min_piv = minval(abs(y_piv(1:n)))
    if (abs(min_piv - aflag(8)) > 1.0e-4 * max(abs(aflag(8)), tiny(1.0))) then
      write(*,'(a,2(1pg12.5,1x))') &
          'FAIL success: min|Y_PIV| /= AFLAG(8), got/exp=', min_piv, aflag(8)
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! B1 must contain the original RHS on exit (saved by y12mf)
    err_b1 = maxval(abs(b1(1:n) - b_orig(1:n)))
    if (err_b1 > 1.0e-5) then
      write(*,'(a,1pg12.5)') &
          'FAIL success: B1 /= original RHS, max_err=', err_b1
      nfail = nfail + 1
      any_fail = .true.
    end if

    ! Solution X must be close to the exact values (1, 2, 3, 4, 5, 6)
    err_sol = maxval(abs(x(1:n) - [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
    if (err_sol > 1.0e-4) then
      write(*,'(a,1pg12.5)') 'FAIL success: solution max_err=', err_sol
      nfail = nfail + 1
      any_fail = .true.
    end if

    if (.not. any_fail) then
      write(*,'(a,i0,a,1pg10.3,a,1pg10.3)') &
          'PASS success: iter=', iflag(12), &
          ' sol_err=', err_sol, ' res_norm=', aflag(10)
    end if

  end block success

  ! =========================================================================
  ! Summary
  ! =========================================================================
  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mf_errors tests PASSED'

contains

  ! Initialise IFLAG to defaults for a fresh y12mf call.
  ! iflag(5) = 2 is required by y12mf (keeps L for residual computation).
  subroutine set_flags(iflag)
    integer, intent(out) :: iflag(12)
    iflag(1)    = 0   ! fresh state
    iflag(2)    = 2   ! Markowitz search width
    iflag(3)    = 1   ! general Markowitz pivoting
    iflag(4)    = 1   ! first system in a same-structure sequence
    iflag(5)    = 2   ! keep L (required by y12mf; iflag(5)=1 -> ifail=10)
    iflag(6:12) = 0
    iflag(11)   = 10  ! maximum refinement iterations
  end subroutine set_flags

  ! Initialise AFLAG to the settings from the problem statement.
  subroutine set_aflag(aflag)
    real, intent(out) :: aflag(11)
    aflag(1)    = 128.0e0
    aflag(2)    = 1.0e-3
    aflag(3)    = 1.0e+16
    aflag(4)    = 1.0e-12
    aflag(5:11) = 0.0
  end subroutine set_aflag

  ! Load the 6x6 reference matrix (COO format, Z_REF=15 non-zeros) and the
  ! corrected right-hand side into the working arrays.
  !
  ! The COO ordering is non-row-major to exercise y12mb's reordering logic
  ! (same ordering as in test_y12mc_errors.F90 and test_y12md_errors.F90).
  !
  ! Matrix (row-major view):
  !   A =
  !     10   0   0   0   0   0
  !      0  12  -3  -1   0   0
  !      0   0  15   0   0   0
  !     -2   0   0  20  -2   0
  !     -1   0   0  -5   1  -1
  !     -1  -2   0   0   0   6
  !
  ! b = (10, 11, 45, 68, -22, 31) -> exact solution x = (1, 2, 3, 4, 5, 6)
  ! Note: b(4) = 68, not 33 as in the problem statement (see discrepancy #1).
  subroutine load_ref6(a, snr, rnr, b)
    real,    intent(out) :: a(NNMAX), b(NMAX)
    integer, intent(out) :: snr(NNMAX), rnr(NN1MAX)
    rnr(1:Z_REF) = [1, 6, 6, 6, 2, 2, 2, 4, 5, 5, 5, 5, 4, 4, 3]
    snr(1:Z_REF) = [1, 6, 2, 1, 2, 3, 4, 1, 1, 6, 5, 4, 4, 5, 3]
    a(1:Z_REF)   = [ 10.0,  6.0, -2.0, -1.0, &
                     12.0, -3.0, -1.0, -2.0, &
                     -1.0, -1.0,  1.0, -5.0, &
                     20.0, -2.0, 15.0]
    b(1:N_REF)   = [10.0, 11.0, 45.0, 68.0, -22.0, 31.0]
  end subroutine load_ref6

  ! Verify ifail equals the expected value; increment nfail on mismatch.
  subroutine check_ifail(label, ifail, expected, nfail)
    character(len=*), intent(in) :: label
    integer, intent(in)          :: ifail, expected
    integer, intent(inout)       :: nfail
    if (ifail /= expected) then
      write(*,'(3a,i0,a,i0)') 'FAIL ', label, &
          ' expected ifail=', expected, ' got ', ifail
      nfail = nfail + 1
    end if
  end subroutine check_ifail

end program test_y12mf_errors
