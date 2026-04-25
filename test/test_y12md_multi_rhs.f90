! SPDX-License-Identifier: GPL-2.0-only
!
! test_y12md_multi_rhs.f90 - verify that the LU factorization produced by
! y12mb+y12mc can be reused to solve a second right-hand-side by calling
! y12md with IFLAG(5)=3 (no new y12mb/y12mc call required).
!
! The documentation (doc, section Y12MD Remark 2) says:
!   "If the LU decomposition of matrix A is available (i.e. a system with
!    matrix A has already been solved) then Y12MB and Y12MC should not be
!    called.  The user should only assign the new right hand side vector to
!    array B, set IFLAG(5) = 3 and call Y12MD."
!
! For this to work, the preceding y12mc call must have been made with
! IFLAG(5)=2 so that the L factors are retained in array A.
!
! Test matrix: 4-by-4 tridiagonal, diag=3, off-diag=-1.
!
!   RHS 1 (b1) = A * [1,1,1,1]^T = [2, 1, 1, 2]  => solution x1 = [1,1,1,1]
!   RHS 2 (b2) = A * [2,1,2,1]^T = [5,-1, 4, 1]  => solution x2 = [2,1,2,1]
!
! Both single-precision (y12mbe/y12mce/y12mde) and double-precision
! (y12mbf/y12mcf/y12mdf) variants are exercised.
!
program test_y12md_multi_rhs
  use y12m, only: y12mb, y12mc, y12md
  implicit none

  integer, parameter :: NNP  = 200
  integer, parameter :: NN1P = 200
  integer, parameter :: NMAX = 10

  integer :: nfail

  nfail = 0
  call test_sp('multi_rhs SP n=4 tridiag', 4, nfail)
  call test_dp('multi_rhs DP n=4 tridiag', 4, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12md_multi_rhs tests PASSED'

contains

  ! -----------------------------------------------------------------------
  ! Single-precision variant
  ! -----------------------------------------------------------------------
  subroutine test_sp(label, n, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    integer,          intent(inout) :: nfail

    real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer :: z, i
    real    :: b1(NMAX), b2(NMAX), x_ref2(NMAX), err

    ! ---- Build tridiagonal matrix ----
    call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr)

    ! ---- RHS 1: b1 = A*[1..1]^T ----
    b1(1) = 2.0 ; b1(2) = 1.0 ; b1(3) = 1.0 ; b1(4) = 2.0

    ! ---- RHS 2: b2 = A*[2,1,2,1]^T ----
    b2(1) = 5.0 ; b2(2) = -1.0 ; b2(3) = 4.0 ; b2(4) = 1.0

    ! Reference solution for RHS 2
    x_ref2(1) = 2.0 ; x_ref2(2) = 1.0 ; x_ref2(3) = 2.0 ; x_ref2(4) = 1.0

    ! ---- Initialise flags for single-system solve, keeping L factors ----
    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 2   ! keep L factors so the second solve can reuse them

    aflag(1) = 16.0     ; aflag(2) = 1.0e-12
    aflag(3) = 1.0e+16  ; aflag(4) = 1.0e-12

    ! ---- Prepare (y12mb) ----
    call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12mb failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    ! ---- Factorise (y12mc) with RHS b1 ----
    b(1:n) = b1(1:n)
    call y12mc(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, NMAX, &
        aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12mc failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    ! ---- First solve (y12md) with IFLAG(5)=2 ----
    call y12md(n, a, NNP, b, pivot, snr, ha, NMAX, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12md (first) failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    ! Check solution x1 ≈ [1,1,1,1]
    err = maxval(abs(b(1:n) - 1.0))
    if (err > 1.0e-4) then
      write(*,'(3a,es10.3)') 'FAIL ', label, ': first solve max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es10.3)') 'PASS ', label, ': first solve max_err=', err
    end if

    ! ---- Second solve: assign new RHS and set IFLAG(5)=3 ----
    b(1:n) = b2(1:n)
    iflag(5) = 3   ! reuse existing LU; y12mb and y12mc must NOT be called

    call y12md(n, a, NNP, b, pivot, snr, ha, NMAX, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12md (second) failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    ! Check solution x2 ≈ [2,1,2,1]
    err = maxval(abs(b(1:n) - x_ref2(1:n)))
    if (err > 1.0e-4) then
      write(*,'(3a,es10.3)') 'FAIL ', label, ': second solve max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es10.3)') 'PASS ', label, ': second solve max_err=', err
    end if
  end subroutine test_sp

  ! -----------------------------------------------------------------------
  ! Double-precision variant
  ! -----------------------------------------------------------------------
  subroutine test_dp(label, n, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    integer,          intent(inout) :: nfail

    double precision :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
    integer          :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10), ifail
    integer          :: z, i
    double precision :: b1(NMAX), b2(NMAX), x_ref2(NMAX), err

    call build_tridiag_dp(n, NNP, NN1P, z, a, snr, rnr)

    b1(1) = 2.0d0 ; b1(2) = 1.0d0 ; b1(3) = 1.0d0 ; b1(4) = 2.0d0
    b2(1) = 5.0d0 ; b2(2) = -1.0d0; b2(3) = 4.0d0 ; b2(4) = 1.0d0
    x_ref2(1) = 2.0d0 ; x_ref2(2) = 1.0d0
    x_ref2(3) = 2.0d0 ; x_ref2(4) = 1.0d0

    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 2

    aflag(1) = 16.0d0    ; aflag(2) = 1.0d-12
    aflag(3) = 1.0d+16   ; aflag(4) = 1.0d-12

    call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12mb failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    b(1:n) = b1(1:n)
    call y12mc(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, NMAX, &
        aflag, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12mc failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    call y12md(n, a, NNP, b, pivot, snr, ha, NMAX, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12md (first) failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    err = maxval(abs(b(1:n) - 1.0d0))
    if (err > 1.0d-10) then
      write(*,'(3a,es12.5)') 'FAIL ', label, ': first solve max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es12.5)') 'PASS ', label, ': first solve max_err=', err
    end if

    b(1:n) = b2(1:n)
    iflag(5) = 3

    call y12md(n, a, NNP, b, pivot, snr, ha, NMAX, iflag, ifail)
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': y12md (second) failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    err = maxval(abs(b(1:n) - x_ref2(1:n)))
    if (err > 1.0d-10) then
      write(*,'(3a,es12.5)') 'FAIL ', label, ': second solve max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es12.5)') 'PASS ', label, ': second solve max_err=', err
    end if
  end subroutine test_dp

  ! -----------------------------------------------------------------------
  ! Matrix builders: 4-by-4 tridiagonal (diag=3, off-diag=-1).
  ! Only fills a/snr/rnr; the caller provides b separately.
  ! -----------------------------------------------------------------------
  subroutine build_tridiag_sp(n, nnmax, nn1max, z, a, snr, rnr)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax)
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
  end subroutine build_tridiag_sp

  subroutine build_tridiag_dp(n, nnmax, nn1max, z, a, snr, rnr)
    integer,          intent(in)  :: n, nnmax, nn1max
    integer,          intent(out) :: z
    double precision, intent(out) :: a(nnmax)
    integer,          intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) =  3.0d0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0d0
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0d0
    end do
  end subroutine build_tridiag_dp

end program test_y12md_multi_rhs
