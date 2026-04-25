! SPDX-License-Identifier: GPL-2.0-only
!
! test_y12mf_output_params.f90 - verify the output-parameter invariants
! specific to y12mf (y12mfe), the iterative-refinement solver.
!
! Only a single-precision version of y12mf exists (y12mfe).
!
! The invariants checked after a successful y12mfe call are:
!
!   IFLAG(12) in [1, IFLAG(11)]   iterations actually performed
!   AFLAG(9)  >= 0                max-norm of last correction vector d(p-1)
!   AFLAG(10) >= 0                max-norm of last residual  vector r(p-1)
!   AFLAG(11) >= 0                max-norm of corrected solution vector x
!   AFLAG(5)  >= 1.0              growth factor (set by internal y12mc call)
!   AFLAG(8)  > 0                 smallest pivot magnitude
!   Y(i) all nonzero              Y holds the pivots (diagonal of U)
!   min|Y| == AFLAG(8)            consistency between Y and AFLAG(8)
!   B1 == original RHS            y12mf saves the RHS in B1 on exit
!   X == solution (accurate)      max|x - x_true| < tolerance
!
! Test matrices:
!   n=4 tridiagonal (diag=3, off-diag=-1), b = A*[1,1,1,1]^T => x=[1,1,1,1]
!   n=5 arrow matrix,                       b = A*[1,...,1]^T => x=[1,...,1]
!
program test_y12mf_output_params
  use y12m, only: y12mf
  implicit none

  integer, parameter :: NNP  = 300
  integer, parameter :: NN1P = 300
  integer, parameter :: NZMAX = 50   ! upper bound on nnz for allocating a1/sn
  integer, parameter :: NMAX = 10

  integer :: nfail

  nfail = 0
  call test_mf('y12mf SP n=4 tridiag', 4, 'tridiag', nfail)
  call test_mf('y12mf SP n=5 arrow',   5, 'arrow',   nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mf_output_params tests PASSED'

contains

  ! -----------------------------------------------------------------------
  ! Run y12mfe on a test matrix and check documented output invariants.
  ! -----------------------------------------------------------------------
  subroutine test_mf(label, n, mtype, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    character(len=*), intent(in)    :: mtype   ! 'tridiag' or 'arrow'
    integer,          intent(inout) :: nfail

    real    :: a(NNP), a1(NZMAX), b(NMAX), b1(NMAX), x(NMAX), y_piv(NMAX)
    real    :: aflag(11), b_orig(NMAX)
    integer :: snr(NNP), rnr(NN1P), sn(NZMAX), ha(NMAX,13), iflag(12)
    integer :: z, nz, ifail, i, max_it
    real    :: min_piv, err_sol, err_b1

    ! ---- Build matrix and RHS ----
    if (mtype == 'tridiag') then
      call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
    else
      call build_arrow_sp(n, NNP, NN1P, z, a, snr, rnr, b)
    end if
    nz = z
    b_orig(1:n) = b(1:n)   ! save original RHS for later comparison

    ! ---- Initialise flags and parameters ----
    max_it = 10

    aflag(1) = 16.0     ! stability factor
    aflag(2) = 1.0e-12  ! drop tolerance
    aflag(3) = 1.0e+16  ! growth limit
    aflag(4) = 1.0e-12  ! pivot threshold
    ! aflag(5..11) are output parameters; no need to set them

    iflag(1)  = 0
    iflag(2)  = 3
    iflag(3)  = 1
    iflag(4)  = 0
    iflag(5)  = 2   ! y12mf requires iflag(5)=2 (not 1)
    iflag(11) = max_it   ! max iterations

    ! ---- Call y12mfe ----
    call y12mf(n, a, snr, NNP, rnr, NN1P, a1, sn, nz, ha, NMAX, &
        b, b1, x, y_piv, aflag, iflag, ifail)

    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': solver failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    ! ---- IFLAG(12): iteration count must be in [1, max_it] ----
    if (iflag(12) < 1 .or. iflag(12) > max_it) then
      write(*,'(3a,i0,a,i0)') 'FAIL ', label, &
          ': IFLAG(12)=', iflag(12), ' not in [1,', max_it
      nfail = nfail + 1
    end if

    ! ---- AFLAG(9) >= 0: max-norm of last correction vector ----
    if (aflag(9) < 0.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(9) < 0, got ', aflag(9)
      nfail = nfail + 1
    end if

    ! ---- AFLAG(10) >= 0: max-norm of last residual vector ----
    if (aflag(10) < 0.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(10) < 0, got ', aflag(10)
      nfail = nfail + 1
    end if

    ! ---- AFLAG(11) >= 0: max-norm of corrected solution ----
    if (aflag(11) < 0.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(11) < 0, got ', aflag(11)
      nfail = nfail + 1
    end if

    ! ---- AFLAG(5) >= 1.0: growth factor ----
    if (aflag(5) < 1.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(5) < 1.0, got ', aflag(5)
      nfail = nfail + 1
    end if

    ! ---- AFLAG(8) > 0: smallest pivot magnitude ----
    if (aflag(8) <= 0.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(8) <= 0, got ', aflag(8)
      nfail = nfail + 1
    end if

    ! ---- Y_PIV contains pivot elements; min|Y| == AFLAG(8) ----
    do i = 1, n
      if (y_piv(i) == 0.0) then
        write(*,'(3a,i0)') 'FAIL ', label, ': Y_PIV(i)==0 at i=', i
        nfail = nfail + 1
      end if
    end do
    min_piv = minval(abs(y_piv(1:n)))
    if (abs(min_piv - aflag(8)) > 1.0e-4 * max(abs(aflag(8)), tiny(1.0))) then
      write(*,'(3a,2(1pg12.5,1x))') 'FAIL ', label, &
          ': min|Y|/=AFLAG(8), got/exp=', min_piv, aflag(8)
      nfail = nfail + 1
    end if

    ! ---- B1 must contain the original RHS on exit ----
    err_b1 = maxval(abs(b1(1:n) - b_orig(1:n)))
    if (err_b1 > 1.0e-5) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, &
          ': B1 /= original RHS, max_err=', err_b1
      nfail = nfail + 1
    end if

    ! ---- X must be close to the true solution [1,...,1] ----
    err_sol = maxval(abs(x(1:n) - 1.0))
    if (err_sol > 1.0e-4) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': solution max_err=', err_sol
      nfail = nfail + 1
    end if

    write(*,'(3a,i0,a,1pg10.3,a,1pg10.3)') 'PASS ', label, &
        ': iter=', iflag(12), ' sol_err=', err_sol, ' res_norm=', aflag(10)
  end subroutine test_mf

  ! -----------------------------------------------------------------------
  ! n-by-n tridiagonal (diag=3, off-diag=-1), RHS = A*[1..1]^T
  ! -----------------------------------------------------------------------
  subroutine build_tridiag_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
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
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_tridiag_sp

  ! -----------------------------------------------------------------------
  ! n-by-n arrow matrix (full first row/col, diagonal rest), x=[1..1]
  ! -----------------------------------------------------------------------
  subroutine build_arrow_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i, j
    z = 0
    do j = 1, n
      z = z + 1 ; rnr(z) = 1 ; snr(z) = j ; a(z) = 1.0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = 1   ; a(z) = 1.0
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) = real(i + 1)
    end do
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_arrow_sp

end program test_y12mf_output_params
