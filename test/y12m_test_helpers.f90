! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! y12m_test_helpers.f90 - Shared helper procedures for the y12m test suite.
!
! Provides generic interfaces for matrix builders, solution checkers, flag
! initialisers, reference-matrix loaders, three-step solve wrappers, and
! sparse matrix-vector products.  Each generic interface dispatches to sp
! (single-precision) and dp (double-precision) specialisations based on the
! kind of the real array arguments.
!
module y12m_test_helpers
  use y12m, only: y12mb, y12mc, y12md
  implicit none
  private

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! ---- Generic interfaces (public) ----------------------------------------
  public :: build_tridiag
  public :: build_arrow
  public :: build_pentadiag
  public :: build_diag
  public :: build_full2
  public :: check_solution
  public :: check_ifail
  public :: set_flags
  public :: set_aflag
  public :: load_ref6
  public :: setup_ref6
  public :: solve3
  public :: sparse_gemv
  public :: forward_error
  public :: backward_error

  !> n-by-n tridiagonal (diag=3, off-diag=-1), solution x=[1,...,1].
  interface build_tridiag
    module procedure build_tridiag_sp
    module procedure build_tridiag_dp
  end interface build_tridiag

  !> n-by-n arrow (full first row/col + diagonal rest), solution x=[1,...,1].
  interface build_arrow
    module procedure build_arrow_sp
    module procedure build_arrow_dp
  end interface build_arrow

  !> n-by-n pentadiagonal (bandwidth 2), solution x=[1,...,1]. dp only.
  interface build_pentadiag
    module procedure build_pentadiag_dp
  end interface build_pentadiag

  !> n-by-n diagonal, solution x=[1,...,1]. sp only.
  interface build_diag
    module procedure build_diag_sp
  end interface build_diag

  !> Fixed 2x2 full matrix [[4,1],[2,5]]. sp only.
  interface build_full2
    module procedure build_full2_sp
  end interface build_full2

  !> Check solution accuracy.
  !>
  !> Two distinct overloads (resolved at compile time by type of argument 4):
  !>
  !>   check_solution(label, n, x, ifail, tol, nfail)
  !>     Verifies IFAIL=0 and max|x(1:n)-1| <= tol.  Arg 4 is integer.
  !>
  !>   check_solution(label, n, x, expected, tol, nfail)
  !>     Verifies max|x(1:n)-expected(1:n)| <= tol.  Arg 4 is a real array.
  interface check_solution
    module procedure check_solution_ones_sp
    module procedure check_solution_ones_dp
    module procedure check_solution_exp_sp
    module procedure check_solution_exp_dp
  end interface check_solution

  !> Initialise AFLAG to standard solver defaults.
  interface set_aflag
    module procedure set_aflag_sp
    module procedure set_aflag_dp
  end interface set_aflag

  !> Load the 6x6 reference matrix (15 non-zeros) and corrected RHS b1.
  interface load_ref6
    module procedure load_ref6_sp
    module procedure load_ref6_dp
  end interface load_ref6

  !> Load reference matrix, initialise flags, and call y12mb.
  interface setup_ref6
    module procedure setup_ref6_sp
    module procedure setup_ref6_dp
  end interface setup_ref6

  !> Three-step y12mb -> y12mc -> y12md wrapper.
  interface solve3
    module procedure solve3_sp
    module procedure solve3_dp
  end interface solve3

  !> Sparse COO matrix-vector product: y := y + alpha * A * x.
  interface sparse_gemv
    module procedure sparse_gemv_sp
    module procedure sparse_gemv_dp
  end interface sparse_gemv

  !> Forward error: maximum absolute difference between x and expected.
  interface forward_error
    module procedure forward_error_sp
    module procedure forward_error_dp
  end interface forward_error

  !> Backward error: residual norm divided by matrix/vector problem scale.
  interface backward_error
    module procedure backward_error_sp
    module procedure backward_error_dp
  end interface backward_error

contains

  ! ===========================================================================
  ! Matrix builders
  ! ===========================================================================

  subroutine build_tridiag_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,    intent(in)  :: n, nnmax, nn1max
    integer,    intent(out) :: z
    real(sp),   intent(out) :: a(nnmax), b(n)
    integer,    intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) =  3.0_sp
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0_sp
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0_sp
    end do
    b(1:n) = 0.0_sp
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_tridiag_sp

  subroutine build_tridiag_dp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,    intent(in)  :: n, nnmax, nn1max
    integer,    intent(out) :: z
    real(dp),   intent(out) :: a(nnmax), b(n)
    integer,    intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) =  3.0_dp
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0_dp
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0_dp
    end do
    b(1:n) = 0.0_dp
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_tridiag_dp

  subroutine build_arrow_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,    intent(in)  :: n, nnmax, nn1max
    integer,    intent(out) :: z
    real(sp),   intent(out) :: a(nnmax), b(n)
    integer,    intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i, j
    z = 0
    do j = 1, n
      z = z + 1 ; rnr(z) = 1 ; snr(z) = j ; a(z) = 1.0_sp
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = 1   ; a(z) = 1.0_sp
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) = real(i + 1, kind=sp)
    end do
    b(1:n) = 0.0_sp
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_arrow_sp

  subroutine build_arrow_dp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,    intent(in)  :: n, nnmax, nn1max
    integer,    intent(out) :: z
    real(dp),   intent(out) :: a(nnmax), b(n)
    integer,    intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i, j
    z = 0
    do j = 1, n
      z = z + 1 ; rnr(z) = 1 ; snr(z) = j ; a(z) = 1.0_dp
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = 1   ; a(z) = 1.0_dp
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) = real(i + 1, kind=dp)
    end do
    b(1:n) = 0.0_dp
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_arrow_dp

  !> n-by-n pentadiagonal (diag=5, dist-1 off-diag=-1, dist-2 off-diag=-0.5).
  subroutine build_pentadiag_dp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,    intent(in)  :: n, nnmax, nn1max
    integer,    intent(out) :: z
    real(dp),   intent(out) :: a(nnmax), b(n)
    integer,    intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i ; a(z) = 5.0_dp
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0_dp
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0_dp
    end do
    do i = 3, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-2 ; a(z) = -0.5_dp
    end do
    do i = 1, n-2
      z = z + 1 ; rnr(z) = i ; snr(z) = i+2 ; a(z) = -0.5_dp
    end do
    b(1:n) = 0.0_dp
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_pentadiag_dp

  !> n-by-n diagonal (diag(i) = i+2, b = row sums), solution x=[1,...,1].
  subroutine build_diag_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,    intent(in)  :: n, nnmax, nn1max
    integer,    intent(out) :: z
    real(sp),   intent(out) :: a(nnmax), b(n)
    integer,    intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i ; a(z) = real(i + 2, kind=sp)
    end do
    b(1:n) = 0.0_sp
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_diag_sp

  !> 2x2 full matrix A=[[4,1],[2,5]], b=A*[1,1]=[5,7], solution x=[1,1].
  subroutine build_full2_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,    intent(in)  :: n, nnmax, nn1max
    integer,    intent(out) :: z
    real(sp),   intent(out) :: a(nnmax), b(n)
    integer,    intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    z = z + 1 ; rnr(z) = 1 ; snr(z) = 1 ; a(z) = 4.0_sp
    z = z + 1 ; rnr(z) = 1 ; snr(z) = 2 ; a(z) = 1.0_sp
    z = z + 1 ; rnr(z) = 2 ; snr(z) = 1 ; a(z) = 2.0_sp
    z = z + 1 ; rnr(z) = 2 ; snr(z) = 2 ; a(z) = 5.0_sp
    b(1:n) = 0.0_sp
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_full2_sp

  ! ===========================================================================
  ! Solution checkers
  ! ===========================================================================

  !> Verify IFAIL=0 and max|x(1:n)-1| <= tol; print PASS/FAIL; bump nfail.
  subroutine check_solution_ones_sp(label, n, x, ifail, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, ifail
    real(sp),         intent(in)    :: x(n), tol
    integer,          intent(inout) :: nfail
    real(sp) :: err
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ' ifail=', ifail
      nfail = nfail + 1
      return
    end if
    err = maxval(abs(x(1:n) - 1.0_sp))
    if (err > tol) then
      write(*,'(3a,es10.3)') 'FAIL ', label, ' max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es10.3)') 'PASS ', label, ' max_err=', err
    end if
  end subroutine check_solution_ones_sp

  !> Verify IFAIL=0 and max|x(1:n)-1| <= tol; print PASS/FAIL; bump nfail.
  subroutine check_solution_ones_dp(label, n, x, ifail, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, ifail
    real(dp),         intent(in)    :: x(n), tol
    integer,          intent(inout) :: nfail
    real(dp) :: err
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ' ifail=', ifail
      nfail = nfail + 1
      return
    end if
    err = maxval(abs(x(1:n) - 1.0_dp))
    if (err > tol) then
      write(*,'(3a,es12.5)') 'FAIL ', label, ' max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es12.5)') 'PASS ', label, ' max_err=', err
    end if
  end subroutine check_solution_ones_dp

  !> Verify max|x(1:n)-expected(1:n)| <= tol; print PASS/FAIL; bump nfail.
  subroutine check_solution_exp_sp(label, n, x, expected, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    real(sp),         intent(in)    :: x(n), expected(n), tol
    integer,          intent(inout) :: nfail
    real(sp) :: err
    err = maxval(abs(x - expected))
    if (err > tol) then
      write(*,'(3a,es10.3)') 'FAIL ', label, ': max|x-expected|=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es10.3)') 'PASS ', label, ': max|x-expected|=', err
    end if
  end subroutine check_solution_exp_sp

  !> Verify max|x(1:n)-expected(1:n)| <= tol; print PASS/FAIL; bump nfail.
  subroutine check_solution_exp_dp(label, n, x, expected, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    real(dp),         intent(in)    :: x(n), expected(n), tol
    integer,          intent(inout) :: nfail
    real(dp) :: err
    err = maxval(abs(x - expected))
    if (err > tol) then
      write(*,'(3a,es12.5)') 'FAIL ', label, ': max|x-expected|=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es12.5)') 'PASS ', label, ': max|x-expected|=', err
    end if
  end subroutine check_solution_exp_dp

  !> Verify ifail == expected; print FAIL if not; bump nfail on mismatch.
  subroutine check_ifail(label, ifail, expected, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: ifail, expected
    integer,          intent(inout) :: nfail
    if (ifail /= expected) then
      write(*,'(3a,i0,a,i0)') 'FAIL ', label, &
          ' expected ifail=', expected, ' got ', ifail
      nfail = nfail + 1
    end if
  end subroutine check_ifail

  ! ===========================================================================
  ! Flag initialisers
  ! ===========================================================================

  !> Initialise IFLAG(1:10) to default values for a fresh single-system solve.
  !>
  !> AFLAG is NOT touched: y12mb ignores all of AFLAG on input and only writes
  !> aflag(6) (max|A|) on exit.
  subroutine set_flags(iflag)
    integer, intent(out) :: iflag(10)
    iflag(1)    = 0   ! state flag: y12mb will set to -1, y12mc to -2
    iflag(2)    = 3   ! Markowitz search width
    iflag(3)    = 1   ! general Markowitz pivoting
    iflag(4)    = 0   ! no structure reuse (fresh factorization)
    iflag(5)    = 1   ! discard L after factorization
    iflag(6:10) = 0
  end subroutine set_flags

  !> Initialise AFLAG(1:8) to standard solver defaults.
  subroutine set_aflag_sp(aflag)
    real(sp), intent(out) :: aflag(8)
    aflag(1) = 16.0_sp      ! stability factor u
    aflag(2) = 1.0e-12_sp   ! drop tolerance
    aflag(3) = 1.0e+16_sp   ! growth threshold
    aflag(4) = 1.0e-12_sp   ! singularity threshold
    aflag(5:8) = 0.0_sp
  end subroutine set_aflag_sp

  !> Initialise AFLAG(1:8) to standard solver defaults.
  subroutine set_aflag_dp(aflag)
    real(dp), intent(out) :: aflag(8)
    aflag(1) = 16.0_dp      ! stability factor u
    aflag(2) = 1.0e-12_dp   ! drop tolerance
    aflag(3) = 1.0e+16_dp   ! growth threshold
    aflag(4) = 1.0e-12_dp   ! singularity threshold
    aflag(5:8) = 0.0_dp
  end subroutine set_aflag_dp

  ! ===========================================================================
  ! Reference 6x6 matrix helpers
  ! ===========================================================================
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
  ! b1 = (10, 11, 45, 68, -22, 31)  ->  exact solution x = (1, 2, 3, 4, 5, 6)
  ! Note: b1(4) = 68, not 33 as in the original problem statement.
  !
  ! The COO ordering is non-row-major to exercise y12mb's reordering logic.

  subroutine load_ref6_sp(a, snr, rnr, b)
    real(sp), intent(inout) :: a(*), b(*)
    integer,  intent(inout) :: snr(*), rnr(*)
    rnr(1:15) = [1, 6, 6, 6, 2, 2, 2, 4, 5, 5, 5, 5, 4, 4, 3]
    snr(1:15) = [1, 6, 2, 1, 2, 3, 4, 1, 1, 6, 5, 4, 4, 5, 3]
    a(1:15)   = [10.0_sp,  6.0_sp, -2.0_sp, -1.0_sp, &
                 12.0_sp, -3.0_sp, -1.0_sp, -2.0_sp, &
                 -1.0_sp, -1.0_sp,  1.0_sp, -5.0_sp, &
                 20.0_sp, -2.0_sp, 15.0_sp]
    b(1:6)    = [10.0_sp, 11.0_sp, 45.0_sp, 68.0_sp, -22.0_sp, 31.0_sp]
  end subroutine load_ref6_sp

  subroutine load_ref6_dp(a, snr, rnr, b)
    real(dp), intent(inout) :: a(*), b(*)
    integer,  intent(inout) :: snr(*), rnr(*)
    rnr(1:15) = [1, 6, 6, 6, 2, 2, 2, 4, 5, 5, 5, 5, 4, 4, 3]
    snr(1:15) = [1, 6, 2, 1, 2, 3, 4, 1, 1, 6, 5, 4, 4, 5, 3]
    a(1:15)   = [10.0_dp,  6.0_dp, -2.0_dp, -1.0_dp, &
                 12.0_dp, -3.0_dp, -1.0_dp, -2.0_dp, &
                 -1.0_dp, -1.0_dp,  1.0_dp, -5.0_dp, &
                 20.0_dp, -2.0_dp, 15.0_dp]
    b(1:6)    = [10.0_dp, 11.0_dp, 45.0_dp, 68.0_dp, -22.0_dp, 31.0_dp]
  end subroutine load_ref6_dp

  !> Load the 6x6 reference matrix and call y12mb.
  !>
  !> On exit the arrays are ready for a y12mc call.  iflag(5) is set to 1 by
  !> set_flags; the caller may override it before calling y12mc.  Records a
  !> failure in nfail when y12mb returns a non-zero ifail.
  subroutine setup_ref6_sp(n, z, nn, nn1, iha, a, snr, rnr, b, ha, &
      aflag, iflag, ifail, nfail)
    integer,  intent(in)    :: n, z, nn, nn1, iha
    real(sp), intent(inout) :: a(nn), aflag(8)
    real(sp), intent(inout) :: b(*)
    integer,  intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
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
  end subroutine setup_ref6_sp

  !> Load the 6x6 reference matrix and call y12mb.
  !>
  !> On exit the arrays are ready for a y12mc call.  iflag(5) is set to 1 by
  !> set_flags; the caller may override it before calling y12mc.  Records a
  !> failure in nfail when y12mb returns a non-zero ifail.
  subroutine setup_ref6_dp(n, z, nn, nn1, iha, a, snr, rnr, b, ha, &
      aflag, iflag, ifail, nfail)
    integer,  intent(in)    :: n, z, nn, nn1, iha
    real(dp), intent(inout) :: a(nn), aflag(8)
    real(dp), intent(inout) :: b(*)
    integer,  intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
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
  end subroutine setup_ref6_dp

  ! ===========================================================================
  ! Three-step solve wrappers
  ! ===========================================================================

  !> Execute the three-step y12mb -> y12mc -> y12md sequence (single precision).
  subroutine solve3_sp(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, &
      iha, aflag, iflag, ifail)
    integer,  intent(in)    :: n, z, nn, nn1, iha
    real(sp), intent(inout) :: a(nn), b(n), pivot(n), aflag(8)
    integer,  intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    integer,  intent(out)   :: ifail
    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0   ! single system; no LU reuse
    iflag(5) = 1   ! discard L after factorization
    aflag(1) = 16.0_sp    ! stability factor
    aflag(2) = 1.0e-12_sp ! drop tolerance
    aflag(3) = 1.0e+16_sp ! max growth factor
    aflag(4) = 1.0e-12_sp ! min pivot threshold
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) return
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
        aflag, iflag, ifail)
    if (ifail /= 0) return
    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  end subroutine solve3_sp

  !> Execute the three-step y12mb -> y12mc -> y12md sequence (double precision).
  subroutine solve3_dp(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, &
      iha, aflag, iflag, ifail)
    integer,  intent(in)    :: n, z, nn, nn1, iha
    real(dp), intent(inout) :: a(nn), b(n), pivot(n), aflag(8)
    integer,  intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    integer,  intent(out)   :: ifail
    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 1
    aflag(1) = 16.0_dp
    aflag(2) = 1.0e-12_dp
    aflag(3) = 1.0e+16_dp
    aflag(4) = 1.0e-12_dp
    call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) return
    call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
        aflag, iflag, ifail)
    if (ifail /= 0) return
    call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  end subroutine solve3_dp

  ! ===========================================================================
  ! Sparse matrix-vector product
  ! ===========================================================================

  !> Accumulate y := y + alpha * A * x for COO matrix (single precision).
  !>
  !> A is stored in coordinate format with nz entries; no ordering assumed.
  subroutine sparse_gemv_sp(m, nz, alpha, val, row, col, x, y)
    integer,  intent(in)    :: m, nz
    real(sp), intent(in)    :: alpha, val(nz)
    integer,  intent(in)    :: row(nz), col(nz)
    real(sp), intent(in)    :: x(*)
    real(sp), intent(inout) :: y(m)
    integer :: k
    do k = 1, nz
      y(row(k)) = y(row(k)) + alpha * val(k) * x(col(k))
    end do
  end subroutine sparse_gemv_sp

  !> Accumulate y := y + alpha * A * x for COO matrix (double precision).
  !>
  !> A is stored in coordinate format with nz entries; no ordering assumed.
  subroutine sparse_gemv_dp(m, nz, alpha, val, row, col, x, y)
    integer,  intent(in)    :: m, nz
    real(dp), intent(in)    :: alpha, val(nz)
    integer,  intent(in)    :: row(nz), col(nz)
    real(dp), intent(in)    :: x(*)
    real(dp), intent(inout) :: y(m)
    integer :: k
    do k = 1, nz
      y(row(k)) = y(row(k)) + alpha * val(k) * x(col(k))
    end do
  end subroutine sparse_gemv_dp

  ! ===========================================================================
  ! Error metrics
  ! ===========================================================================

  subroutine forward_error_sp(n, x, expected, err)
    integer,  intent(in)  :: n
    real(sp), intent(in)  :: x(n), expected(n)
    real(sp), intent(out) :: err
    err = maxval(abs(x(1:n) - expected(1:n)))
  end subroutine forward_error_sp

  subroutine forward_error_dp(n, x, expected, err)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x(n), expected(n)
    real(dp), intent(out) :: err
    err = maxval(abs(x(1:n) - expected(1:n)))
  end subroutine forward_error_dp

  subroutine backward_error_sp(n, nz, val, row, col, x, b, err)
    integer,  intent(in)  :: n, nz
    real(sp), intent(in)  :: val(nz), x(n), b(n)
    integer,  intent(in)  :: row(nz), col(nz)
    real(sp), intent(out) :: err

    integer :: k
    real(sp) :: a_inf, x_inf, b_inf, r_inf, denom
    real(sp) :: row_sum(n), residual(n), ax(n)

    row_sum = 0.0_sp
    do k = 1, nz
      row_sum(row(k)) = row_sum(row(k)) + abs(val(k))
    end do
    a_inf = maxval(row_sum)

    ax = 0.0_sp
    call sparse_gemv(n, nz, 1.0_sp, val, row, col, x, ax)
    residual = b - ax
    r_inf = maxval(abs(residual))

    x_inf = maxval(abs(x))
    b_inf = maxval(abs(b))
    ! Guard against exact-zero denominator.
    denom = max(a_inf * x_inf + b_inf, tiny(1.0_sp))
    err = r_inf / denom
  end subroutine backward_error_sp

  subroutine backward_error_dp(n, nz, val, row, col, x, b, err)
    integer,  intent(in)  :: n, nz
    real(dp), intent(in)  :: val(nz), x(n), b(n)
    integer,  intent(in)  :: row(nz), col(nz)
    real(dp), intent(out) :: err

    integer :: k
    real(dp) :: a_inf, x_inf, b_inf, r_inf, denom
    real(dp) :: row_sum(n), residual(n), ax(n)

    row_sum = 0.0_dp
    do k = 1, nz
      row_sum(row(k)) = row_sum(row(k)) + abs(val(k))
    end do
    a_inf = maxval(row_sum)

    ax = 0.0_dp
    call sparse_gemv(n, nz, 1.0_dp, val, row, col, x, ax)
    residual = b - ax
    r_inf = maxval(abs(residual))

    x_inf = maxval(abs(x))
    b_inf = maxval(abs(b))
    ! Guard against exact-zero denominator.
    denom = max(a_inf * x_inf + b_inf, tiny(1.0_dp))
    err = r_inf / denom
  end subroutine backward_error_dp

end module y12m_test_helpers
