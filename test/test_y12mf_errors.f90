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
! Only y12mf-specific logic is exercised here.  Error paths that are delegated
! to y12mb, y12mc, or y12md (e.g. ifail=12/13/5/15 etc.) are already covered
! in test_y12mb_errors.F90 and test_y12mc_errors.F90 and are not repeated.
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
!   3. The old `doc` states that B holds the last correction vector d(p-1) on
!      successful exit.  Inspection of y12mfe.f shows that after the refinement
!      loop B actually holds the last residual vector r(p) = b - A*x(p), whose
!      max-norm equals AFLAG(10).  This test checks max|B(1:n)| == AFLAG(10).
!
! Error codes exercised (y12mf-specific only):
!
!   ifail=10  iflag(5) = 1  (y12mf requires iflag(5) = 2 or 3; checked before
!                            any call to y12mb/y12mc/y12md)
!   ifail=23  iflag(11) = 0 (minimum iteration count is 2; IFLAG(11) > 1)
!   ifail=23  iflag(11) = 1 (same check, boundary case)
!
! y12mf-specific output invariants verified after a successful call:
!
!   * IFLAG(12) in [1, IFLAG(11)]   actual refinement iterations performed
!   * AFLAG(9)  >= 0.0              max-norm of last correction vector d(p-1)
!   * AFLAG(10) >= 0.0              max-norm of last residual vector r(p)
!   * AFLAG(11) >= 1.0              max-norm of corrected solution x(p)
!   * AFLAG(9)/AFLAG(11) small      relative-error estimate from the doc
!   * max|B(1:n)| == AFLAG(10)      B holds last residual (see discrepancy #3)
!   * B1 == original RHS            y12mf saves original b in B1 on exit
!   * X close to (1, 2, 3, 4, 5, 6) accurate solution for first RHS
!
! LU-reuse verified (IFLAG(5) = 3):
!
!   After the first solve the LU factorization stored in A/SNR/HA is reused
!   for a second RHS: b2 = (60, 45, 60, 44, -20, -10) -> x2 = (6,5,4,3,2,1).
!   The same y12mf-specific output invariants are checked for the second solve.
!
program test_y12mf_errors
  use y12m, only: y12mf
  use y12m_test_helpers, only: load_ref6, check_ifail
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
  ! y12mfe sets ifail=10 immediately when state = iflag(5) = 1, before any
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
  ! ifail = 23 : iflag(11) < 2  (restriction: IFLAG(11) > 1)
  !
  ! y12mfe checks it = iflag(11) < 2 immediately after the state check;
  ! returns ifail=23 without modifying any work array.
  ! Two values are tested: 0 and 1 (both below the minimum of 2).
  ! =========================================================================
  blk23a: block
    integer :: n = N_REF, nz = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    iflag(11) = 0   ! < 2 -> ifail=23
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=23 iflag(11)=0', ifail, 23, nfail)
  end block blk23a

  blk23b: block
    integer :: n = N_REF, nz = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    iflag(11) = 1   ! < 2 -> ifail=23 (boundary case)
    ha = 0
    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y_piv, aflag, iflag, ifail)
    call check_ifail('ifail=23 iflag(11)=1', ifail, 23, nfail)
  end block blk23b

  ! =========================================================================
  ! Success: first solve with 6x6 reference matrix
  !
  ! Settings from problem statement:
  !   AFLAG(1:4) = [128.0, 1.e-3, 1.e+16, 1.e-12]
  !   IFLAG(2:5) = [2, 1, 1, 2]
  !   IFLAG(11)  = 25
  !
  ! Expected solution: x = (1, 2, 3, 4, 5, 6).
  !
  ! After this block the factorization is reused in the LU-reuse block below.
  ! =========================================================================
  success: block
    integer :: n = N_REF, nz = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    integer :: max_it
    real    :: b_orig(NMAX)

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

    call check_mf_outputs('success', N_REF, max_it, &
        b, b1, x, b_orig, aflag, iflag, &
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], nfail)

  end block success

  ! =========================================================================
  ! LU-reuse: second solve with IFLAG(5) = 3
  !
  ! The LU factorization stored in A/SNR/HA by the first solve is reused.
  ! A1/SN hold the original matrix (row-ordered) and must not be altered.
  !
  ! b2 = (60, 45, 60, 44, -20, -10) -> exact solution x2 = (6, 5, 4, 3, 2, 1)
  ! =========================================================================
  reuse: block
    integer :: n = N_REF, max_it
    real    :: b_orig2(NMAX)

    b(1:N_REF)     = [60.0, 45.0, 60.0, 44.0, -20.0, -10.0]
    b_orig2(1:N_REF) = b(1:N_REF)

    ! AFLAG and most of IFLAG are unchanged from the first solve; only
    ! iflag(5) must be updated to signal that the factorization is available.
    iflag(5)  = 3   ! reuse existing LU
    iflag(12) = 0
    max_it    = iflag(11)

    call y12mf(n, a, snr, NNMAX, rnr, NN1MAX, a1, sn, Z_REF, ha, NMAX, &
        b, b1, x, y_piv, aflag, iflag, ifail)

    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL reuse: y12mf returned ifail=', ifail
      nfail = nfail + 1
      exit reuse
    end if

    call check_mf_outputs('reuse', N_REF, max_it, &
        b, b1, x, b_orig2, aflag, iflag, &
        [6.0, 5.0, 4.0, 3.0, 2.0, 1.0], nfail)

  end block reuse

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
    iflag(11)   = 25  ! maximum refinement iterations
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

  ! Check all y12mf-specific output invariants after a successful call.
  !
  ! Checked:
  !   IFLAG(12) in [1, max_it]
  !   AFLAG(9)  >= 0              max-norm of last correction vector d(p-1)
  !   AFLAG(10) >= 0              max-norm of last residual vector r(p)
  !   AFLAG(11) >= 1              max-norm of corrected solution (must be >= 1
  !                               since exact solution has at least one component >= 1)
  !   AFLAG(9)/AFLAG(11) < tol   relative-error estimate (doc Section 9 / AFLAG(11))
  !   max|B(1:n)| == AFLAG(10)   B holds last residual, its norm equals AFLAG(10)
  !                               (see discrepancy #3 in file header)
  !   B1 == b_orig                y12mf stores the original RHS in B1 on exit
  !   X close to x_exact          solution accuracy
  subroutine check_mf_outputs(label, n, max_it, b, b1, x, b_orig, aflag, iflag, &
      x_exact, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, max_it
    real,             intent(in)    :: b(NMAX), b1(NMAX), x(NMAX)
    real,             intent(in)    :: b_orig(NMAX), x_exact(n)
    real,             intent(in)    :: aflag(11)
    integer,          intent(in)    :: iflag(12)
    integer,          intent(inout) :: nfail

    real :: err_sol, err_b1, res_norm, rel_err

    ! IFLAG(12): actual iteration count must be in [1, max_it]
    if (iflag(12) < 1 .or. iflag(12) > max_it) then
      write(*,'(3a,i0,a,i0,a)') 'FAIL ', label, ': IFLAG(12)=', iflag(12), &
          ' not in [1,', max_it, ']'
      nfail = nfail + 1
    end if

    ! AFLAG(9) >= 0: max-norm of last correction vector d(p-1)
    if (aflag(9) < 0.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(9) < 0, got ', aflag(9)
      nfail = nfail + 1
    end if

    ! AFLAG(10) >= 0: max-norm of last residual vector r(p)
    if (aflag(10) < 0.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(10) < 0, got ', aflag(10)
      nfail = nfail + 1
    end if

    ! AFLAG(11) >= 1.0: max-norm of corrected solution must be at least 1
    ! (all test solutions have components with absolute value >= 1).
    if (aflag(11) < 1.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(11) < 1.0, got ', aflag(11)
      nfail = nfail + 1
    end if

    ! AFLAG(9)/AFLAG(11): relative-error estimate documented in `doc` Section 9.
    ! Should be well below single-precision machine epsilon for a converged solve.
    if (aflag(11) > 0.0) then
      rel_err = aflag(9) / aflag(11)
      if (rel_err > 1.0e-3) then
        write(*,'(3a,1pg12.5)') 'FAIL ', label, &
            ': AFLAG(9)/AFLAG(11) too large, got ', rel_err
        nfail = nfail + 1
      end if
    end if

    ! max|B(1:n)| must equal AFLAG(10): after the refinement loop y12mfe stores
    ! the last residual r(p) in B; its max-norm is AFLAG(10).
    res_norm = maxval(abs(b(1:n)))
    if (abs(res_norm - aflag(10)) > 1.0e-5 * max(aflag(10), tiny(1.0))) then
      write(*,'(3a,2(1pg12.5,1x))') 'FAIL ', label, &
          ': max|B| /= AFLAG(10), got/exp=', res_norm, aflag(10)
      nfail = nfail + 1
    end if

    ! B1 must contain the original RHS on exit (y12mf copies b -> b1 on entry)
    err_b1 = maxval(abs(b1(1:n) - b_orig(1:n)))
    if (err_b1 > 1.0e-5) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, &
          ': B1 /= original RHS, max_err=', err_b1
      nfail = nfail + 1
    end if

    ! Solution X must be close to the exact values
    err_sol = maxval(abs(x(1:n) - x_exact(1:n)))
    if (err_sol > 1.0e-4) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': solution max_err=', err_sol
      nfail = nfail + 1
    end if

    write(*,'(3a,i0,a,1pg10.3,a,1pg10.3,a,1pg10.3)') &
        'PASS ', label, ': iter=', iflag(12), &
        ' sol_err=', err_sol, ' res_norm=', aflag(10), &
        ' rel_err=', aflag(9)/max(aflag(11), tiny(1.0))
  end subroutine check_mf_outputs

end program test_y12mf_errors
