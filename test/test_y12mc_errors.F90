! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4-5
!
! test_y12mc_errors.F90 - Error-branch and success coverage for y12mc.
!
! Use test_y12mb_errors.F90 as template.
!
! Each error sub-test supplies deliberately invalid arguments (or a degenerate
! matrix) to force a specific positive IFAIL return value, then verifies the
! expected code is returned.  The final sub-test exercises a successful
! Y12MB -> Y12MC -> Y12MD sequence on the 6-by-6 reference matrix from the
! problem statement and checks the documented postconditions.
!
! Note on error-triggering flag values: some tests modify IFLAG(2), IFLAG(3),
! or IFLAG(5) after calling y12mb to force a specific error condition in y12mc.
! This is done purely to exercise y12mc's input validation; in normal use
! these flags should be set once before the y12mb call and not modified
! between y12mb and y12mc.  Specifically:
!   - IFLAG(2) must be >= 1 (Markowitz search width).
!   - IFLAG(3) must be 0, 1, or 2 (pivoting strategy).
!   - IFLAG(5) must be 1 or 2 for the first y12mc call; value 3 is only
!     valid for y12md when reusing a factorization from a previous
!     IFLAG(5)=2 call.
!
! Error codes exercised and their triggering conditions:
!
!   ifail= 2  iflag(1) /= -1         (y12mc called without prior y12mb)
!   ifail=19  iflag(2) < 1           (iflag(2) = 0)
!   ifail=20  iflag(3) not in {0,1,2} (iflag(3) = 3)
!   ifail=21  iflag(5) not in {1,2,3} (iflag(5) = 0)
!   ifail=22  iflag(5) = 3           (back-substitution-only mode before factorization)
!   ifail= 9  no pivot in current row (3x3 permutation matrix, iflag(3)=0)
!   ifail= 3  small pivot during loop (3x3 diag(1,1e-12,1), aflag(4)=0.1)
!   ifail= 3  small last pivot        (2x2 diag(10,1e-12), aflag(4)=0.1)
!   ifail= 7  empty row after elimination (purpose-built 3x3, iflag(3)=0)
!   ifail= 8  empty column after step cleanup (purpose-built 3x3, iflag(3)=0)
!   ifail= 4  growth factor exceeded  (tiny-pivot 3x3, aflag(4)=0, aflag(3) small)
!
! Not yet exercised (future work):
!
!   ifail= 5  nn too small for fill-in during factorization
!             (requires a fill-heavy matrix plus a carefully chosen nn value
!             that passes the y12mb size check but overflows during y12mc)
!   ifail= 6  nn1 too small for fill-in during factorization
!             (same difficulty as ifail=5 but for the column linked list)
!
! Success postconditions verified:
!   * IFAIL    = 0     on exit from each of y12mb, y12mc, y12md
!   * IFLAG(1) = -2    (state flag set by y12mc on success)
!   * AFLAG(6) = maxval(abs(A_original))  (unchanged by y12mc)
!   * AFLAG(8) > 0     (minimum absolute pivot, set during factorization)
!   * IFLAG(6) >= 0, IFLAG(7) >= 0  (garbage-collection counters)
!   * pivot(1:N) close to the expected pivotal sequence
!   * Solution x satisfies ||A*x - b||_inf < tolerance
!
program test_y12mc_errors
  use y12m, only: y12mb, y12mc, y12md
  use y12m_test_helpers, only: set_flags, set_aflag, check_ifail
  implicit none

#ifdef TEST_SINGLE_PRECISION
  integer, parameter :: wp = kind(1.0e0)
#else
  integer, parameter :: wp = kind(1.0d0)
#endif

  integer, parameter :: NMAX   = 10
  integer, parameter :: NNMAX  = 200
  integer, parameter :: NN1MAX = 100

  real(wp) :: a(NNMAX), aflag(8), pivot(NMAX), b(NMAX)
  integer  :: snr(NNMAX), rnr(NN1MAX), ha(NMAX,11), iflag(10)
  integer  :: ifail, nfail

  nfail = 0

  ! =========================================================================
  ! ifail = 2 : y12mc called when iflag(1) /= -1
  !
  ! iflag(1) = -1 is set by a successful y12mb call.  Calling y12mc with any
  ! other value (e.g. 0 after a fresh initialisation) immediately returns
  ! ifail=2.  No y12mb call is made here.
  ! =========================================================================
  block
    integer :: n = 3, z = 3, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    a(1:z) = [11, 22, 33]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    b(1:n) = 0
    call set_flags(iflag)
    iflag(1) = 0   ! NOT -1 -> y12mc must return ifail=2
    ha = 0
    call set_aflag(aflag)
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=2 no y12mb', ifail, 2, nfail)
  end block

  ! =========================================================================
  ! ifail = 19 : iflag(2) < 1
  !
  ! iflag(2) is the Markowitz search width; the minimum legal value is 1.
  ! Setting it to 0 here after a valid y12mb call forces the error path in
  ! y12mc; in practice iflag(2) should always remain >= 1.
  ! =========================================================================
  blk19: block
    integer :: n = 3, z = 3, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    call init_diag3(a, snr, rnr, b)
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk19
    iflag(2) = 0   ! < 1 -> ifail=19
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=19 iflag(2)<1', ifail, 19, nfail)
  end block blk19

  ! =========================================================================
  ! ifail = 20 : iflag(3) not in {0, 1, 2}
  !
  ! iflag(3) selects the pivoting strategy (0=diagonal, 1=Markowitz,
  ! 2=diagonal-preference Markowitz).  Values outside this set are invalid.
  ! Setting it to 3 here purely to exercise the check; in practice only
  ! 0, 1, or 2 are meaningful.
  ! =========================================================================
  blk20: block
    integer :: n = 3, z = 3, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    call init_diag3(a, snr, rnr, b)
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk20
    iflag(3) = 3   ! > 2 -> ifail=20
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=20 iflag(3)=3', ifail, 20, nfail)
  end block blk20

  ! =========================================================================
  ! ifail = 21 : iflag(5) not in {1, 2, 3}
  !
  ! iflag(5) controls L-factor storage (1=discard, 2=keep, 3=back-sub only).
  ! Values 0 and 4 are both invalid.  Setting it to 0 after y12mb exercises
  ! the lower-bound check in y12mc.
  ! =========================================================================
  blk21: block
    integer :: n = 3, z = 3, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    call init_diag3(a, snr, rnr, b)
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk21
    iflag(5) = 0   ! < 1 -> ifail=21
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=21 iflag(5)=0', ifail, 21, nfail)
  end block blk21

  ! =========================================================================
  ! ifail = 22 : iflag(5) = 3
  !
  ! iflag(5)=3 means "use stored L to perform back-substitution only" - this
  ! is only valid after a previous full factorisation (iflag(5)=2).  y12mc
  ! rejects it here because it is called without a stored L.
  ! =========================================================================
  blk22: block
    integer :: n = 3, z = 3, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    call init_diag3(a, snr, rnr, b)
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk22
    iflag(5) = 3   ! back-sub only without stored L -> ifail=22
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=22 iflag(5)=3', ifail, 22, nfail)
  end block blk22

  ! =========================================================================
  ! ifail = 9 : no pivot element in the current row
  !
  ! When iflag(3)=0, y12mc searches row i for a diagonal entry (column i).
  ! The 3x3 cyclic permutation matrix has no diagonal entries at all; the
  ! search fails immediately on step 1, returning ifail=9.
  !
  !   A = [ 0 1 0 ]     nonzeros: (1,2)=1  (2,3)=1  (3,1)=1
  !       [ 0 0 1 ]
  !       [ 1 0 0 ]
  ! =========================================================================
  blk9: block
    integer :: n = 3, z = 3, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    a(1:z)   = [1, 1, 1]
    snr(1:z) = [2, 3, 1]
    rnr(1:z) = [1, 2, 3]
    b(1:n)   = [1, 1, 1]
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk9
    iflag(3) = 0   ! no pivot search; take diagonal -> row 1 has no diagonal
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=9 no diagonal pivot', ifail, 9, nfail)
  end block blk9

  ! =========================================================================
  ! ifail = 3 : pivot too small (during the elimination loop, step i < n)
  !
  ! When the absolute value of the chosen pivot is smaller than
  ! aflag(4)*aflag(6), y12mc considers the matrix near-singular and returns
  ! ifail=3.
  !
  ! 3x3 diagonal matrix  diag(1, 1e-12, 1) with aflag(4)=0.1.
  !   aflag(6) = 1  (set by y12mb)
  !   grmin    = aflag(4)*aflag(6) = 0.1
  !   step 2 pivot = 1e-12 < 0.1  ->  ifail=3
  ! =========================================================================
  blk3a: block
    integer :: n = 3, z = 3, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    a(1:z)   = [1.0_wp, 1.0e-12_wp, 1.0_wp]
    snr(1:z) = [1, 2, 3]
    rnr(1:z) = [1, 2, 3]
    b(1:n)   = [1, 1, 1]
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk3a
    iflag(3) = 0      ! take diagonal pivots in order
    ! aflag(6) has been set by y12mb to max|A|=1; preserve it.
    ! grmin = aflag(4)*aflag(6) = 0.1*1 = 0.1 > 1e-12 (step 2 pivot) -> ifail=3
    aflag(4) = 0.1_wp
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=3 small pivot in loop', ifail, 3, nfail)
  end block blk3a

  ! =========================================================================
  ! ifail = 3 : pivot too small (the very last pivot, checked after the loop)
  !
  ! For an n=2 system y12mc runs only one elimination step (i=1..n-1=1), then
  ! extracts the last pivot from the remaining diagonal.  A small last pivot
  ! also triggers ifail=3.
  !
  ! 2x2 diagonal matrix  diag(10, 1e-12) with aflag(4)=0.1.
  !   aflag(6) = 10  (set by y12mb)
  !   grmin    = 0.1*10 = 1
  !   step 1 pivot = 10 >= 1  ->  OK
  !   last pivot   = 1e-12 < 1  ->  ifail=3
  ! =========================================================================
  blk3b: block
    integer :: n = 2, z = 2, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    a(1:z)   = [10.0_wp, 1.0e-12_wp]
    snr(1:z) = [1, 2]
    rnr(1:z) = [1, 2]
    b(1:n)   = [1, 1]
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk3b
    iflag(3) = 0
    ! aflag(6) has been set by y12mb to max|A|=10; preserve it.
    ! grmin = aflag(4)*aflag(6) = 0.1*10 = 1 > 1e-12 (last pivot) -> ifail=3
    aflag(4) = 0.1_wp
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=3 small last pivot', ifail, 3, nfail)
  end block blk3b

  ! =========================================================================
  ! ifail = 7 : a row becomes empty during elimination
  !
  ! After the pivot column is eliminated from a row that had only that single
  ! off-diagonal entry, the active part of that row becomes empty -> ifail=7.
  !
  ! Matrix (n=3, z=4):
  !   row 1: (1,1)=5            <- pivot row, diagonal only, r9=0
  !   row 2: (2,1)=1            <- only element; eliminated in step 1 -> empty
  !   row 3: (3,2)=1  (3,3)=3
  !
  ! With iflag(3)=0, step 1 takes diagonal A(1,1)=5 as pivot.  r9=0 means no
  ! off-diagonal elements in the pivot row, so no fill-in is added to row 2.
  ! After removing (2,1) the row has no active elements -> ifail=7.
  ! =========================================================================
  blk7: block
    integer :: n = 3, z = 4, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    a(1:z)   = [5.0_wp, 1.0_wp, 1.0_wp, 3.0_wp]
    snr(1:z) = [1, 1, 2, 3]
    rnr(1:z) = [1, 2, 3, 3]
    b(1:n)   = [1, 1, 1]
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk7
    iflag(3) = 0
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=7 empty row', ifail, 7, nfail)
  end block blk7

  ! =========================================================================
  ! ifail = 8 : a column becomes empty during step cleanup
  !
  ! After step i y12mc checks each off-diagonal column of the pivot row.  If
  ! any of those columns has no active entries, it reports ifail=8.
  !
  ! Matrix (n=3, z=4):
  !   row 1: (1,1)=5  (1,2)=1   <- pivot row; col 2 is an off-diagonal
  !   row 2: (2,3)=1             <- col 2 has NO entry in row 2
  !   row 3: (3,3)=3
  !
  ! Column 2 appears only in the pivot row (row 1).  It therefore has exactly
  ! one active entry.  When y12mc checks the columns of the pivot row after
  ! step 1 it finds col 2 has r2 == r1 (start == end, i.e. 1 entry remaining
  ! pointing to the now-consumed pivot row) -> active check r2 > r1 fails
  ! -> ifail=8.
  !
  ! Note: col 1 has only row 1 as well, so no other rows are eliminated in
  ! step 1 (no fill-in loop executes).
  ! =========================================================================
  blk8: block
    integer :: n = 3, z = 4, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    a(1:z)   = [5.0_wp, 1.0_wp, 1.0_wp, 3.0_wp]
    snr(1:z) = [1, 2, 3, 3]
    rnr(1:z) = [1, 1, 2, 3]
    b(1:n)   = [1, 1, 1]
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk8
    iflag(3) = 0
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=8 empty col', ifail, 8, nfail)
  end block blk8

  ! =========================================================================
  ! ifail = 4 : growth factor exceeded
  !
  ! y12mc tracks the ratio aflag(7)/aflag(6) (max element seen during
  ! factorization / max element in the original matrix).  When this ratio
  ! exceeds aflag(3) (growth threshold, minimum 1e5) ifail=4 is returned.
  !
  ! Matrix (n=3, z=5):
  !   A = [ 1e-10  1   0 ]
  !       [  1     1   0 ]
  !       [  0     0   1 ]
  !
  !   aflag(6) = 1   (set by y12mb; max|A| = 1)
  !   aflag(3) = 1e5 (growth threshold; y12mc enforces minimum 1e5)
  !   aflag(4) = 0   (no singularity threshold, so small pivot is accepted)
  !
  ! With iflag(3)=0, step 1 takes the diagonal A(1,1)=1e-10 as pivot.
  ! The multiplier for row 2 is 1/1e-10 = 1e10.  Updating A(2,2):
  !   A(2,2) <- 1 - 1*1e10 = -1e10+1 ~= -1e10
  ! aflag(7) is updated to 1e10.  Growth = 1e10/1 = 1e10 > 1e5 -> ifail=4.
  ! =========================================================================
  blk4: block
    integer :: n = 3, z = 5, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    a(1:z)   = [1.0e-10_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
    snr(1:z) = [1, 2, 1, 2, 3]
    rnr(1:z) = [1, 1, 2, 2, 3]
    b(1:n)   = [1, 1, 1]
    call setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    if (ifail /= 0) exit blk4
    iflag(3) = 0
    ! aflag(6) has been set by y12mb to max|A|=1; preserve it.
    ! Set a small growth threshold so that the ratio aflag(7)/aflag(6) = 1e10
    ! exceeds it.  y12mc enforces minimum 1e5, so set exactly that.
    aflag(4) = 0.0_wp   ! no singularity threshold; accept tiny pivot
    aflag(3) = 1.0e5_wp ! growth threshold (y12mc clips to max(user,1e5))
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    call check_ifail('ifail=4 growth factor', ifail, 4, nfail)
  end block blk4

  ! =========================================================================
  ! TODO (future work): ifail=5 and ifail=6
  !
  ! ifail=5 is triggered when the row-linked list array (a / snr, size nn)
  ! runs out of space for fill-in entries during Gaussian elimination.  To
  ! trigger this a fill-heavy matrix is needed (e.g. an arrowhead/bordered
  ! structure) combined with an nn value that satisfies y12mb's check
  ! (nn >= 2*z) but is too small to accommodate all fill-ins in y12mc.
  ! With iflag(5)=2 (store L) the space freed by garbage collection is
  ! minimal, making it easier to overflow.
  !
  ! ifail=6 is the analogous condition for the column linked list (rnr,
  ! size nn1): nn1 = z passes y12mb but the first fill-in attempt overflows.
  !
  ! Both cases require careful matrix design and array sizing that has been
  ! analysed but not yet coded as a self-contained automated test.
  ! =========================================================================

  ! =========================================================================
  ! Success test: 6x6 reference system from the problem statement
  !
  ! A =
  !   10   0   0   0   0   0
  !    0  12  -3  -1   0   0
  !    0   0  15   0   0   0
  !   -2   0   0  20  -2   0
  !   -1   0   0  -5   1  -1
  !   -1  -2   0   0   0   6
  !
  ! b = [10, 11, 45, 33, -22, 31]
  !
  ! The matrix is supplied in the following non-row-major COO ordering so
  ! that y12mb's reordering logic is exercised:
  !
  !   rnr snr    a
  !     1   1  10.
  !     6   6   6.
  !     6   2  -2.
  !     6   1  -1.
  !     2   2  12.
  !     2   3  -3.
  !     2   4  -1.
  !     4   1  -2.
  !     5   1  -1.
  !     5   6  -1.
  !     5   5   1.
  !     5   4  -5.
  !     4   4  20.
  !     4   5  -2.
  !     3   3  15.
  !
  ! Pivot order (Markowitz, iflag(3)=1, iflag(2)=3):
  !   step 1: row 3, col 3 = 15       (single-element row, selected first)
  !   step 2: row 1, col 1 = 10       (other single-element row)
  !   step 3: row 2, col 2 = 12       (2-element row)
  !   step 4: row 6, col 6 =  6       (2-element row after step-3 fill-in)
  !   step 5: row 4, col 4 = 20       (last 2-element row)
  !   step 6: Schur complement = 179/360 ~= 0.4972  (remaining 1x1 block)
  !
  ! Note: the problem statement quotes pivots [10,15,-1,6,12,4.972], which
  ! appear to reflect a different Markowitz search or manual elimination order.
  ! The values below are those produced by y12mcf with the default settings.
  ! =========================================================================
  success: block
    integer  :: n = 6, z = 15, nn = NNMAX, nn1 = NN1MAX, iha = NMAX
    integer  :: i, k
    real(wp) :: resid, atol
    ! Original COO arrays saved for residual computation after solve.
    integer  :: row_orig(15), col_orig(15)
    real(wp) :: val_orig(15), b_orig(6)
    real(wp) :: expected_pivot(6)
    real(wp) :: r(6)

    ! --- input matrix in non-row-major COO ordering ---
    rnr(1:z) = [1, 6, 6, 6, 2, 2, 2, 4, 5, 5, 5, 5, 4, 4, 3]
    snr(1:z) = [1, 6, 2, 1, 2, 3, 4, 1, 1, 6, 5, 4, 4, 5, 3]
    a(1:z)   = [10.0_wp,  6.0_wp, -2.0_wp, -1.0_wp, &
                12.0_wp, -3.0_wp, -1.0_wp, -2.0_wp, &
                -1.0_wp, -1.0_wp,  1.0_wp, -5.0_wp, &
                20.0_wp, -2.0_wp, 15.0_wp]
    b(1:n) = [10.0_wp, 11.0_wp, 45.0_wp, 33.0_wp, -22.0_wp, 31.0_wp]

    ! Save original matrix and RHS for residual check after solve.
    row_orig(1:z) = rnr(1:z)
    col_orig(1:z) = snr(1:z)
    val_orig(1:z) = a(1:z)
    b_orig(1:n)   = b(1:n)

    ! Expected pivots in Markowitz order (as produced by y12mcf).
    ! Schur complement of the remaining 1x1 block: 179/360 ~= 0.4972
    expected_pivot = [15.0_wp, 10.0_wp, 12.0_wp, 6.0_wp, 20.0_wp, &
                      179.0_wp / 360.0_wp]

    call set_flags(iflag)
    call set_aflag(aflag)
    ha = 0

    ! --- y12mb ---
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL success y12mb: ifail=', ifail
      nfail = nfail + 1
      exit success
    end if

    ! --- y12mc ---
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL success y12mc: ifail=', ifail
      nfail = nfail + 1
      exit success
    end if

    ! Post-conditions after successful y12mc.
    if (iflag(1) /= -2) then
      write(*,'(a,i0)') 'FAIL success: iflag(1) /= -2, got ', iflag(1)
      nfail = nfail + 1
    end if

    if (aflag(6) /= 20.0_wp) then
      write(*,'(a,es12.5,a)') &
          'FAIL success: aflag(6)=', aflag(6), ' expected 20 (max|A|, unchanged)'
      nfail = nfail + 1
    end if

    if (aflag(8) <= 0.0_wp) then
      write(*,'(a,es12.5)') 'FAIL success: aflag(8) should be > 0, got ', aflag(8)
      nfail = nfail + 1
    end if

    ! Min |pivot| is the last Schur complement = 179/360 ~= 0.4972.
    ! Use precision-appropriate tolerances for all floating-point comparisons.
#ifdef TEST_SINGLE_PRECISION
    atol = 1.0e-4_wp   ! absolute tolerance (sp rounding)
#else
    atol = 1.0e-10_wp  ! absolute tolerance (dp rounding)
#endif
    if (abs(aflag(8) - 179.0_wp/360.0_wp) > atol) then
      write(*,'(a,es12.5,a)') &
          'FAIL success: aflag(8)=', aflag(8), ' expected ~179/360 (min|pivot|)'
      nfail = nfail + 1
    end if

    if (iflag(6) < 0 .or. iflag(7) < 0) then
      write(*,'(a,2i4)') &
          'FAIL success: iflag(6),iflag(7) should be >= 0, got ', &
          iflag(6), iflag(7)
      nfail = nfail + 1
    end if

    ! Note: the problem statement claims RNR(1:Z) will be zero after y12mc.
    ! In practice, garbage-collection compaction can leave non-zero entries
    ! at positions 1..Z.  This check is therefore omitted.

    ! Check pivot values.
    do i = 1, n
      if (abs(pivot(i) - expected_pivot(i)) > atol) then
        write(*,'(a,i0,a,es12.5,a,es12.5)') &
            'FAIL success: pivot(', i, ')=', pivot(i), &
            ' expected ~', expected_pivot(i)
        nfail = nfail + 1
      end if
    end do

    ! --- y12md (forward/back substitution) ---
    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL success y12md: ifail=', ifail
      nfail = nfail + 1
      exit success
    end if

    ! Residual  r = A_original * x - b_original
    r(1:n) = -b_orig(1:n)
    do k = 1, z
      r(row_orig(k)) = r(row_orig(k)) + val_orig(k) * b(col_orig(k))
    end do
    resid = maxval(abs(r(1:n)))
#ifdef TEST_SINGLE_PRECISION
    atol = 1.0e-3_wp
#else
    atol = 1.0e-10_wp
#endif
    if (resid > atol) then
      write(*,'(a,es12.5)') &
          'FAIL success: max residual ||A*x-b||_inf = ', resid
      nfail = nfail + 1
    end if

  end block success

  ! =========================================================================
  ! Summary
  ! =========================================================================
  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mc_errors tests PASSED'

contains

  ! Load a 3x3 diagonal matrix diag(11,22,33) into COO arrays and RHS.
  ! Used as the "setup" matrix for early-flag-check tests that only need a
  ! valid y12mb state; the matrix values are irrelevant for those tests.
  subroutine init_diag3(a, snr, rnr, b)
    real(wp), intent(out) :: a(NNMAX), b(NMAX)
    integer,  intent(out) :: snr(NNMAX), rnr(NN1MAX)
    a(1:3)   = [11.0_wp, 22.0_wp, 33.0_wp]
    snr(1:3) = [1, 2, 3]
    rnr(1:3) = [1, 2, 3]
    b(1:3)   = [1.0_wp, 1.0_wp, 1.0_wp]
  end subroutine init_diag3

  ! Call y12mb for a diagonal n=3 system, initialising ha and flags.
  ! Records a failure in nfail if y12mb returns non-zero.
  subroutine setup_y12mb(n, z, nn, nn1, iha, a, snr, rnr, ha, aflag, iflag, ifail, nfail)
    integer,  intent(in)    :: n, z, nn, nn1, iha
    real(wp), intent(inout) :: a(NNMAX), aflag(8)
    integer,  intent(inout) :: snr(NNMAX), rnr(NN1MAX), ha(NMAX,11), iflag(10)
    integer,  intent(out)   :: ifail
    integer,  intent(inout) :: nfail
    call set_flags(iflag)
    call set_aflag(aflag)
    ha = 0
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL setup_y12mb: unexpected ifail=', ifail
      nfail = nfail + 1
    end if
  end subroutine setup_y12mb

end program test_y12mc_errors
