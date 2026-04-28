! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12ma_y12mf_bench_sp.f90
!
! Compare single-precision Y12MA (direct solve) and Y12MF (iterative
! refinement) on four benchmark matrices taken from the literature.
! All matrices are stored in COO (triplet) format used by y12m:
!   rnr(k) = row index of entry k
!   snr(k) = column index of entry k
!   a(k)   = value of entry k
!
! Test problems:
!   1. UMFPACK 5x5 — UMFPACK User Guide, section 5.2
!      n=5, nnz=12, x_expected = [1, 2, 3, 4, 5]
!      COO in column-sorted order (as given in the User Guide).
!   2. Harwell-Boeing 4x4 — HSL MA50 test problem
!      n=4, nnz=10, x_expected = [1, 1, 1, 1]
!   3. Pardiso 8x8 unsymmetric — Pardiso User Guide (pardiso_unsym.f)
!      CSR converted to COO; b = A * [1,...,1] (row sums).
!   4. Templates 6x6 — "Templates for the Solution of Linear Systems"
!      CSR converted to COO; b = A * [1,...,1] (row sums).
!
! Pass criteria: IFAIL=0 and max|x_i - x_expected_i| < 1e-4 for both
!               Y12MA and Y12MF on every test problem.
!
program test_y12ma_y12mf_bench_sp
  use y12m, only: y12ma, y12mf
  implicit none

  integer :: nfail
  nfail = 0

  ! -----------------------------------------------------------------------
  ! Test 1: UMFPACK 5x5 (UMFPACK User Guide, section 5.2)
  !   COO: column-sorted order (rnr=row, snr=column).
  !   Known solution: x = [1, 2, 3, 4, 5].
  ! -----------------------------------------------------------------------
  call run_test('umfpack_5x5', &
      n    = 5, &
      a0   = [2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.], &
      rnr0 = [1, 2, 1, 3, 5, 2, 3, 4, 5, 3, 2, 5], &
      snr0 = [1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5], &
      b0   = [8., 45., -3., 3., 19.], &
      xexp = [1., 2., 3., 4., 5.], &
      nfail = nfail)

  ! -----------------------------------------------------------------------
  ! Test 2: Harwell-Boeing 4x4 (HSL MA50 test problem)
  !   Known solution: x = [1, 1, 1, 1].
  ! -----------------------------------------------------------------------
  call run_test('harwell_boeing_4x4', &
      n    = 4, &
      a0   = [1., 4., 6., 7., 8., 9., 11., 12., 14., 16.], &
      rnr0 = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4], &
      snr0 = [1, 4, 2, 3, 4, 1, 3, 4, 2, 4], &
      b0   = [5., 21., 32., 30.], &
      xexp = [1., 1., 1., 1.], &
      nfail = nfail)

  ! -----------------------------------------------------------------------
  ! Test 3: Pardiso 8x8 unsymmetric (Pardiso User Guide, pardiso_unsym.f)
  !   Stored in CSR in the source; converted here to COO (row-sorted).
  !   b = A * [1,...,1] = row sums, so x_expected = [1,...,1].
  ! -----------------------------------------------------------------------
  call run_test('pardiso_8x8', &
      n    = 8, &
      a0   = [ 7.,  1.,  2.,  7., &
              -4.,  8.,  2., &
               1.,  5., &
               7.,  9., &
              -4., &
               7.,  3.,  8., &
               1., 11., &
              -3.,  2.,  5.], &
      rnr0 = [1, 1, 1, 1,  2, 2, 2,  3, 3,  4, 4,  5,  6, 6, 6,  7, 7,  8, 8, 8], &
      snr0 = [1, 3, 6, 7,  2, 3, 5,  3, 8,  4, 7,  2,  3, 6, 8,  2, 7,  3, 7, 8], &
      b0   = [17., 6., 6., 16., -4., 18., 12., 4.], &
      xexp = [1., 1., 1., 1., 1., 1., 1., 1.], &
      nfail = nfail)

  ! -----------------------------------------------------------------------
  ! Test 4: Templates 6x6 (Templates for the Solution of Linear Systems,
  !         https://netlib.org/linalg/html_templates/node91.html)
  !   Stored in CSR in the reference; converted here to COO (row-sorted).
  !   b = A * [1,...,1] = row sums, so x_expected = [1,...,1].
  ! -----------------------------------------------------------------------
  call run_test('templates_6x6', &
      n    = 6, &
      a0   = [10., -2., &
               3.,  9.,  3., &
               7.,  8.,  7., &
               3.,  8.,  7.,  5., &
               8.,  9.,  9., 13., &
               4.,  2., -1.], &
      rnr0 = [1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4, 4,  5, 5, 5, 5,  6, 6, 6], &
      snr0 = [1, 5,  1, 2, 6,  2, 3, 4,  1, 3, 4, 5,  2, 4, 5, 6,  2, 5, 6], &
      b0   = [8., 15., 22., 23., 39., 5.], &
      xexp = [1., 1., 1., 1., 1., 1.], &
      nfail = nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_y12mf_bench_sp tests PASSED'

contains

  ! Run both y12ma (direct) and y12mf (iterative refinement) on the matrix
  ! described by the COO triplet (a0, rnr0, snr0) with right-hand side b0.
  ! Solutions are compared against xexp; failures increment nfail.
  subroutine run_test(label, n, a0, rnr0, snr0, b0, xexp, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    real,             intent(in)    :: a0(:), b0(n), xexp(n)
    integer,          intent(in)    :: rnr0(:), snr0(:)
    integer,          intent(inout) :: nfail

    integer, parameter :: NNP  = 200
    integer, parameter :: NN1P = 200

    integer :: nz, i, ifail

    ! y12ma workspace (11-column ha, 8-element aflag, 10-element iflag)
    real    :: a_ma(NNP), b_ma(n), pivot_ma(n), aflag_ma(8)
    integer :: snr_ma(NNP), rnr_ma(NN1P), ha_ma(n, 11), iflag_ma(10)

    ! y12mf workspace (13-column ha, 11-element aflag, 12-element iflag).
    ! a1_mf and sn_mf preserve the original matrix values and column indices
    ! so that y12mf can compute residuals during iterative refinement.
    real    :: a_mf(NNP), b_mf(n), b1_mf(n), x_mf(n), y_mf(n), aflag_mf(11)
    real    :: a1_mf(size(a0))
    integer :: snr_mf(NNP), rnr_mf(NN1P), ha_mf(n, 13), iflag_mf(12)
    integer :: sn_mf(size(a0))

    nz = size(a0)

    ! ---- y12ma: direct solve ----
    ! y12ma sets AFLAG(1-4) and IFLAG(2-5) internally; only IFLAG(1) matters.
    a_ma(1:nz)   = a0
    snr_ma(1:nz) = snr0
    rnr_ma(1:nz) = rnr0
    b_ma(1:n)    = b0
    iflag_ma(1)  = 0
    call y12ma(n, nz, a_ma, snr_ma, NNP, rnr_ma, NN1P, pivot_ma, ha_ma, n, &
        aflag_ma, iflag_ma, b_ma, ifail)
    call check_result(label//' y12ma', n, b_ma, xexp, ifail, 1.0e-4, nfail)

    ! ---- y12mf: iterative refinement ----
    ! a1_mf and sn_mf must be initialised to zero; y12mf fills them with the
    ! original matrix data during the first call (IFLAG(5)=2).
    a_mf(1:nz)   = a0
    snr_mf(1:nz) = snr0
    rnr_mf(1:nz) = rnr0
    b_mf(1:n)    = b0
    b1_mf(1:n)   = 0.0
    x_mf(1:n)    = 0.0
    y_mf(1:n)    = 0.0
    do i = 1, nz
      a1_mf(i) = 0.0
      sn_mf(i) = 0
    end do

    aflag_mf(1)  = 16.0
    aflag_mf(2)  = 1.0e-12
    aflag_mf(3)  = 1.0e+16
    aflag_mf(4)  = 1.0e-12
    iflag_mf(1)  = 0
    iflag_mf(2)  = 3
    iflag_mf(3)  = 1
    iflag_mf(4)  = 1
    iflag_mf(5)  = 2   ! keep L factors (required by y12mf)
    iflag_mf(11) = 10  ! maximum refinement iterations

    call y12mf(n, a_mf, snr_mf, NNP, rnr_mf, NN1P, a1_mf, sn_mf, nz, &
        ha_mf, n, b_mf, b1_mf, x_mf, y_mf, aflag_mf, iflag_mf, ifail)
    call check_result(label//' y12mf', n, x_mf, xexp, ifail, 1.0e-4, nfail)

  end subroutine run_test

  ! Check IFAIL=0 and max|sol-xexp|<tol; increment nfail and print on failure.
  subroutine check_result(label, n, sol, xexp, ifail, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, ifail
    real,             intent(in)    :: sol(n), xexp(n), tol
    integer,          intent(inout) :: nfail
    real :: err
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ' ifail=', ifail
      nfail = nfail + 1
      return
    end if
    err = maxval(abs(sol(1:n) - xexp(1:n)))
    if (err > tol) then
      write(*,'(3a,es10.3)') 'FAIL ', label, ' max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es10.3)') 'PASS ', label, ' max_err=', err
    end if
  end subroutine check_result

end program test_y12ma_y12mf_bench_sp
