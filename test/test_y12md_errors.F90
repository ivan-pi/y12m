! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4-5
!
! test_y12md_errors.F90 - Error-branch and success coverage for y12md.
!
! Use test_y12mc_errors.F90 as template.
!
! DISCREPANCIES NOTED:
!
!   1. Problem-statement b1(4): The problem statement gives
!      b1 = (10, 11, 45, 33, -22, 31) with claimed exact solution
!      x = (1, 2, 3, 4, 5, 6).  However, for the 6x6 reference matrix,
!      A(4,:)*[1,2,3,4,5,6] = -2*1 + 20*4 - 2*5 = 68, not 33.
!      This test therefore uses the corrected b1 = (10, 11, 45, 68, -22, 31)
!      so that the exact solution x = (1, 2, 3, 4, 5, 6) can be verified.
!      The second right-hand side b2 = (60, 45, 60, 44, -20, -10) with
!      exact solution x = (6, 5, 4, 3, 2, 1) is consistent and unchanged.
!
!   2. y12mde (single-precision variant) does not initialise ifail = 0 before
!      the iflag(1) precondition check; y12mdf (double-precision) does set
!      ifail = 0 on entry.  On the successful path through y12mde, ifail is
!      never written.  Tests below pre-initialise ifail = 0 before each y12md
!      call so the postcondition ifail = 0 is verifiable with both variants.
!
!   3. y12md performs no input-validation checks beyond iflag(1) == -2.
!      In particular it does NOT verify:
!        NN >= 2*Z  (when iflag(5) = 1; note: the docs say NN > 2*Z --
!                    the exact boundary condition >= vs > is ambiguous)
!        NN >= 3*Z  (when iflag(5) = 2, L is stored alongside U)
!        IHA >= N
!      Commented-out test blocks at the end of this file document what those
!      checks would look like once added to the implementation.
!      See GitHub issue: "y12md: add input-validation checks for NN and IHA".
!
! Error codes exercised and their triggering conditions:
!
!   ifail=1  iflag(1) /= -2   (y12md called without a preceding y12mc call)
!
!   Note: beyond the iflag(1) check the y12md implementation contains no
!   further input validation.  All other invalid-input scenarios are
!   documented only as commented-out blocks below.
!
! Success postconditions verified:
!   * IFAIL    = 0     on exit from y12md
!   * IFLAG(1) = -2    unchanged on exit from y12md
!   * Solution x = (1,2,3,4,5,6) for b1  (iflag(5) = 1 and iflag(5) = 2)
!   * Multiple-RHS: x = (1,2,3,4,5,6) for b1, x = (6,5,4,3,2,1) for b2
!
program test_y12md_errors
  use y12m, only: y12mb, y12mc, y12md
  implicit none

#ifdef TEST_SINGLE_PRECISION
  integer, parameter :: dp = kind(1.0e0)
#else
  integer, parameter :: dp = kind(1.0d0)
#endif

  integer, parameter :: NMAX = 10, NNMAX = 200, NN1MAX = 100
  integer, parameter :: N_REF = 6, Z_REF = 15   ! reference matrix: 6x6, 15 non-zeros

  real(dp) :: a(NNMAX), aflag(8), pivot(NMAX), b(NMAX), atol
  integer  :: snr(NNMAX), rnr(NN1MAX), ha(NMAX,11), iflag(10)
  integer  :: ifail, nfail

#ifdef TEST_SINGLE_PRECISION
  atol = 1.0e-3_dp
#else
  atol = 1.0e-10_dp
#endif
  nfail = 0

  ! =========================================================================
  ! ifail = 1 : iflag(1) /= -2  (y12md called without a prior y12mc)
  !
  ! Any value of iflag(1) other than -2 triggers an immediate ifail=1 return.
  ! Tested with iflag(1) = 0 (fresh) and iflag(1) = -1 (after y12mb only).
  ! =========================================================================

  blk1a: block
    integer :: n = 3, nn = NNMAX, iha = NMAX
    a(1:3)   = [11.0_dp, 22.0_dp, 33.0_dp]
    snr(1:3) = [1, 2, 3]
    pivot(1:3) = [11.0_dp, 22.0_dp, 33.0_dp]
    b(1:n)   = [1.0_dp, 1.0_dp, 1.0_dp]
    ha = 0
    call set_flags(iflag)
    iflag(1) = 0     ! NOT -2  ->  y12md must return ifail=1
    ifail = 99       ! sentinel: must be overwritten to 1
    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
    call check_ifail('ifail=1 iflag(1)=0', ifail, 1, nfail)
  end block blk1a

  blk1b: block
    integer :: n = 3, nn = NNMAX, iha = NMAX
    a(1:3)   = [11.0_dp, 22.0_dp, 33.0_dp]
    snr(1:3) = [1, 2, 3]
    pivot(1:3) = [11.0_dp, 22.0_dp, 33.0_dp]
    b(1:n)   = [1.0_dp, 1.0_dp, 1.0_dp]
    ha = 0
    call set_flags(iflag)
    iflag(1) = -1    ! set by y12mb, but y12mc not yet called  ->  ifail=1
    ifail = 99
    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
    call check_ifail('ifail=1 iflag(1)=-1', ifail, 1, nfail)
  end block blk1b

  ! =========================================================================
  ! Success: y12mb -> y12mc(iflag(5)=1) -> y12md
  !
  ! iflag(5)=1: L is discarded after factorization.  y12mc has already
  ! transformed b to c = L^{-1}Pb; y12md performs only the U-solve followed
  ! by the column permutation.
  !
  ! Expected solution: x = (1, 2, 3, 4, 5, 6).
  ! =========================================================================
  success_no_l: block
    integer  :: n = N_REF, z = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    call setup_ref6(n, z, nn, nn1, iha, a, snr, rnr, b, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit success_no_l

    ! iflag(5) = 1 is already set by set_flags inside setup_ref6.
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0 .or. iflag(1) /= -2) then
      write(*,'(a,i0,a,i0)') 'FAIL success_no_l y12mc: ifail=', ifail, &
          ' iflag(1)=', iflag(1)
      nfail = nfail + 1
      exit success_no_l
    end if

    ifail = 0   ! pre-initialise (needed for y12mde, see discrepancy #2)
    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
    if (ifail /= 0 .or. iflag(1) /= -2) then
      write(*,'(a,i0,a,i0)') 'FAIL success_no_l y12md: ifail=', ifail, &
          ' iflag(1)=', iflag(1)
      nfail = nfail + 1
      exit success_no_l
    end if

    call check_solution('success_no_l', n, b(1:n), &
        [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], atol, nfail)

  end block success_no_l

  ! =========================================================================
  ! Success: y12mb -> y12mc(iflag(5)=2) -> y12md
  !
  ! iflag(5)=2: L is retained in the working arrays for potential
  ! multiple-RHS reuse.  The solve itself follows the same code path as
  ! the iflag(5)=1 case (only the U-solve is performed in y12md).
  !
  ! Expected solution: x = (1, 2, 3, 4, 5, 6).
  ! =========================================================================
  success_with_l: block
    integer  :: n = N_REF, z = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX

    call setup_ref6(n, z, nn, nn1, iha, a, snr, rnr, b, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit success_with_l

    iflag(5) = 2   ! keep L for multiple-RHS reuse
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0 .or. iflag(1) /= -2) then
      write(*,'(a,i0,a,i0)') 'FAIL success_with_l y12mc: ifail=', ifail, &
          ' iflag(1)=', iflag(1)
      nfail = nfail + 1
      exit success_with_l
    end if

    ifail = 0   ! pre-initialise (needed for y12mde, see discrepancy #2)
    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
    if (ifail /= 0 .or. iflag(1) /= -2) then
      write(*,'(a,i0,a,i0)') 'FAIL success_with_l y12md: ifail=', ifail, &
          ' iflag(1)=', iflag(1)
      nfail = nfail + 1
      exit success_with_l
    end if

    call check_solution('success_with_l', n, b(1:n), &
        [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], atol, nfail)

  end block success_with_l

  ! =========================================================================
  ! Multiple RHS: y12mb -> y12mc(iflag(5)=2) -> y12md(b1) ->
  !               [iflag(5)=3] -> y12md(b2)
  !
  ! After factorizing once with iflag(5)=2 (which keeps L), the factorization
  ! can be reused for additional right-hand sides by setting iflag(5)=3 before
  ! each subsequent y12md call.  y12md then applies the row permutation,
  ! performs the L-solve, performs the U-solve, and applies the column
  ! permutation.
  !
  ! b1 = (10, 11, 45, 68, -22, 31) -> x1 = (1, 2, 3, 4, 5, 6)
  ! b2 = (60, 45, 60, 44, -20, -10) -> x2 = (6, 5, 4, 3, 2, 1)
  ! =========================================================================
  multi_rhs: block
    integer  :: n = N_REF, z = Z_REF, nn = NNMAX, nn1 = NN1MAX, iha = NMAX, k
    real(dp) :: rhs(NMAX,2), expected(N_REF,2)
    character(len=12), parameter :: labels(2) = &
        [character(len=12) :: 'multi_rhs x1', 'multi_rhs x2']

    ! b1 = (10,11,45,68,-22,31) -> x1 = (1..6)
    ! b2 = (60,45,60,44,-20,-10) -> x2 = (6..1)
    rhs(1:N_REF, 2)     = [60.0_dp, 45.0_dp, 60.0_dp, 44.0_dp, -20.0_dp, -10.0_dp]
    expected(1:N_REF,1) = [1.0_dp,  2.0_dp,  3.0_dp,  4.0_dp,  5.0_dp,  6.0_dp]
    expected(1:N_REF,2) = [6.0_dp,  5.0_dp,  4.0_dp,  3.0_dp,  2.0_dp,  1.0_dp]

    call setup_ref6(n, z, nn, nn1, iha, a, snr, rnr, rhs(:,1), ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit multi_rhs

    ! Keep L: y12mc transforms rhs(:,1) to c = L^{-1}Pb1; rhs(:,2) is untouched.
    iflag(5) = 2
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, rhs(:,1), ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL multi_rhs y12mc: ifail=', ifail
      nfail = nfail + 1
      exit multi_rhs
    end if

    ! Loop over right-hand sides.
    ! k=1: iflag(5)=2, rhs(:,1) already holds c=L^{-1}Pb1 (U-solve + col-perm only).
    ! k=2: iflag(5)=3, rhs(:,2) holds raw b2 (full row-perm + L-solve + U-solve + col-perm).
    do k = 1, 2
      call y12md(n, a, nn, rhs(:,k), pivot, snr, ha, iha, iflag, ifail)
      if (ifail /= 0) then
        write(*,'(a,i0,a,i0)') 'FAIL multi_rhs y12md k=', k, ': ifail=', ifail
        nfail = nfail + 1
        exit multi_rhs
      end if
      call check_solution(labels(k), n, rhs(1:n,k), expected(1:n,k), atol, nfail)

      ! Prepare to reuse L for the next right-hand side.
      iflag(5) = 3
    end do

  end block multi_rhs

  ! =========================================================================
  ! TODO: Commented-out tests for checks not present in y12md
  !
  ! The following blocks would test input-validation paths that y12md does not
  ! currently implement.  Enabling them requires adding the corresponding
  ! checks to the y12md source (see discrepancy #3 above).
  !
  ! The expected error code for NN violations is suggested to be ifail=5
  ! (consistent with y12mb/y12mc), and ifail=15 for IHA < N (consistent with
  ! y12mb).  These are tentative and should be confirmed when the checks are
  ! added.
  ! =========================================================================

  ! --- NN < 2*Z when iflag(5) = 1 (L not stored) ---
  ! After a successful y12mb + y12mc(iflag(5)=1), call y12md with
  ! nn = 2*Z_REF - 1 = 29.  The implementation should return ifail=5.
  ! Note: it is unclear from the docs whether the limit is NN > 2*Z or
  ! NN >= 2*Z; once a check is implemented this test will clarify the boundary.
  !
  ! blk_nn_no_l: block
  !   integer :: n = N_REF, z = Z_REF, nn_bad = 2*Z_REF - 1, &
  !              nn = NNMAX, nn1 = NN1MAX, iha = NMAX
  !   call setup_ref6(n, z, nn, nn1, iha, a, snr, rnr, b, ha, aflag, iflag, &
  !                   ifail, nfail)
  !   if (ifail /= 0) exit blk_nn_no_l
  !   ! iflag(5) = 1 already set.
  !   call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, &
  !              iflag, ifail)
  !   if (ifail /= 0) exit blk_nn_no_l
  !   ifail = 0
  !   call y12md(n, a, nn_bad, b, pivot, snr, ha, iha, iflag, ifail)
  !   call check_ifail('nn<2z iflag(5)=1', ifail, 5, nfail)
  ! end block blk_nn_no_l

  ! --- NN < 3*Z when iflag(5) = 2 (L stored) ---
  ! After a successful y12mb + y12mc(iflag(5)=2), call y12md with
  ! nn = 3*Z_REF - 1 = 44.  The implementation should return ifail=5.
  !
  ! blk_nn_with_l: block
  !   integer :: n = N_REF, z = Z_REF, nn_bad = 3*Z_REF - 1, &
  !              nn = NNMAX, nn1 = NN1MAX, iha = NMAX
  !   call setup_ref6(n, z, nn, nn1, iha, a, snr, rnr, b, ha, aflag, iflag, &
  !                   ifail, nfail)
  !   if (ifail /= 0) exit blk_nn_with_l
  !   iflag(5) = 2
  !   call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, &
  !              iflag, ifail)
  !   if (ifail /= 0) exit blk_nn_with_l
  !   ifail = 0
  !   call y12md(n, a, nn_bad, b, pivot, snr, ha, iha, iflag, ifail)
  !   call check_ifail('nn<3z iflag(5)=2', ifail, 5, nfail)
  ! end block blk_nn_with_l

  ! --- IHA < N ---
  ! After a successful y12mb + y12mc, call y12md with iha_bad = N_REF - 1 = 5.
  ! The implementation should return ifail=15 (consistent with y12mb).
  !
  ! blk_iha: block
  !   integer :: n = N_REF, z = Z_REF, nn = NNMAX, nn1 = NN1MAX, &
  !              iha = NMAX, iha_bad = N_REF - 1
  !   call setup_ref6(n, z, nn, nn1, iha, a, snr, rnr, b, ha, aflag, iflag, &
  !                   ifail, nfail)
  !   if (ifail /= 0) exit blk_iha
  !   call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, &
  !              iflag, ifail)
  !   if (ifail /= 0) exit blk_iha
  !   ifail = 0
  !   call y12md(n, a, nn, b, pivot, snr, ha, iha_bad, iflag, ifail)
  !   call check_ifail('iha<n', ifail, 15, nfail)
  ! end block blk_iha

  ! =========================================================================
  ! Summary
  ! =========================================================================
  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12md_errors tests PASSED'

contains

  ! Initialise IFLAG to defaults for a fresh single-system solve.
  subroutine set_flags(iflag)
    integer, intent(out) :: iflag(10)
    iflag(1)    = 0   ! state flag: y12mb will set to -1, y12mc to -2
    iflag(2)    = 3   ! Markowitz search width
    iflag(3)    = 1   ! general Markowitz pivoting
    iflag(4)    = 0   ! no structure reuse
    iflag(5)    = 1   ! discard L after factorization
    iflag(6:10) = 0
  end subroutine set_flags

  ! Initialise AFLAG to sensible defaults (mirrors y12ma defaults).
  subroutine set_aflag(aflag)
    real(dp), intent(out) :: aflag(8)
    aflag(1) = 16.0_dp       ! stability factor
    aflag(2) = 1.0e-12_dp    ! drop tolerance
    aflag(3) = 1.0e+16_dp    ! growth threshold
    aflag(4) = 1.0e-12_dp    ! singularity threshold
    aflag(5:8) = 0.0_dp
  end subroutine set_aflag

  ! Load the 6x6 reference matrix (COO format, Z_REF=15) and b1 into the
  ! working arrays.  The COO ordering is non-row-major to exercise y12mb's
  ! reordering logic (same ordering as in test_y12mc_errors.F90).
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
  ! b1 = (10, 11, 45, 68, -22, 31)  ->  exact solution x = (1,2,3,4,5,6)
  ! (b1(4) = 68, not 33 as in the problem statement; see discrepancy #1)
  subroutine load_ref6(a, snr, rnr, b)
    real(dp), intent(out) :: a(NNMAX), b(NMAX)
    integer,  intent(out) :: snr(NNMAX), rnr(NN1MAX)
    rnr(1:Z_REF) = [1, 6, 6, 6, 2, 2, 2, 4, 5, 5, 5, 5, 4, 4, 3]
    snr(1:Z_REF) = [1, 6, 2, 1, 2, 3, 4, 1, 1, 6, 5, 4, 4, 5, 3]
    a(1:Z_REF)   = [10.0_dp,  6.0_dp, -2.0_dp, -1.0_dp, &
                    12.0_dp, -3.0_dp, -1.0_dp, -2.0_dp, &
                    -1.0_dp, -1.0_dp,  1.0_dp, -5.0_dp, &
                    20.0_dp, -2.0_dp, 15.0_dp]
    b(1:N_REF) = [10.0_dp, 11.0_dp, 45.0_dp, 68.0_dp, -22.0_dp, 31.0_dp]
  end subroutine load_ref6

  ! Load the 6x6 reference matrix and call y12mb.  Records a failure in nfail
  ! if y12mb returns non-zero ifail.  On return, the arrays are in the state
  ! required for a y12mc call; iflag(5) = 1 (caller may override before y12mc).
  subroutine setup_ref6(n, z, nn, nn1, iha, a, snr, rnr, b, ha, aflag, iflag, ifail, nfail)
    integer,  intent(in)    :: n, z, nn, nn1, iha
    real(dp), intent(inout) :: a(NNMAX), aflag(8), b(NMAX)
    integer,  intent(inout) :: snr(NNMAX), rnr(NN1MAX), ha(NMAX,11), iflag(10)
    integer,  intent(out)   :: ifail
    integer,  intent(inout) :: nfail
    call load_ref6(a, snr, rnr, b)
    call set_flags(iflag)
    call set_aflag(aflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL setup_ref6 y12mb: unexpected ifail=', ifail
      nfail = nfail + 1
    end if
  end subroutine setup_ref6

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

  ! Verify that max|x - expected| <= tol; increment nfail on failure.
  subroutine check_solution(label, n, x, expected, tol, nfail)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: n
    real(dp),         intent(in) :: x(n), expected(n), tol
    integer,          intent(inout) :: nfail
    real(dp) :: err
    err = maxval(abs(x(1:n) - expected(1:n)))
    if (err > tol) then
      write(*,'(3a,es12.5)') 'FAIL ', label, ': max|x-expected|=', err
      nfail = nfail + 1
    end if
  end subroutine check_solution

end program test_y12md_errors
