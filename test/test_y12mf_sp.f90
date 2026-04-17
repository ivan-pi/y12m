! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12mf_sp.f90 - Tests for y12mf (single-precision iterative refinement).
!
! This program tests y12mf (resolving to y12mfe) which combines Gaussian
! elimination with iterative refinement to improve solution accuracy.  Unlike
! y12ma/y12mb+y12mc+y12md, y12mf preserves the original matrix (in a1/sn)
! and uses it to compute residuals and refine x until convergence.
!
! Four tridiagonal (diag=3, off-diag=-1) systems are solved:
!   n=3, 5, 7, 10  — all with the known solution x=[1,...,1].
!
! The factored form is stored with IFLAG(5)=2 (keep L factors) as required
! by y12mf.  The refined solution is returned in array x (not b).
!
! Pass criteria: IFAIL=0 and max|x_i-1|<1e-4 for all systems.
!
program test_y12mf_sp
  use y12m
  implicit none

  integer, parameter :: NMAX  = 12
  integer, parameter :: NNP   = 400
  integer, parameter :: NN1P  = 400

  real    :: a(NNP), a1(NNP), b(NMAX), b1(NMAX), x(NMAX), y(NMAX)
  real    :: aflag(11)
  integer :: snr(NNP), sn(NNP), rnr(NN1P), ha(NMAX,13), iflag(12)
  integer :: n, z, nn, nn1, iha, ifail, nfail

  nfail = 0

  ! Test sizes: n = 3, 5, 7, 10
  call run_test(3)
  call run_test(5)
  call run_test(7)
  call run_test(10)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mf_sp tests PASSED'

contains

  subroutine run_test(n_in)
    integer, intent(in) :: n_in
    character(len=30) :: label
    integer :: i, z_loc

    n = n_in
    write(label,'(a,i0,a)') 'y12mf_sp n=', n, ' tridiag'

    call build_tridiag_sp(n, NNP, NN1P, z_loc, a, snr, rnr, b)
    z = z_loc
    nn = NNP ; nn1 = NN1P ; iha = NMAX

    ! Settings per the reference example (test2.f).
    ! IFLAG(5)=2 keeps L factors as required by y12mf.
    aflag(1)  = 16.0    ! stability factor
    aflag(2)  = 1.0e-12 ! drop tolerance
    aflag(3)  = 1.0e+16 ! max growth factor before abort
    aflag(4)  = 1.0e-12 ! min pivot threshold

    iflag(1)  = 0
    iflag(2)  = 3
    iflag(3)  = 1
    iflag(4)  = 1       ! preserve structure info for iterative refinement
    iflag(5)  = 2       ! keep L factors (required by y12mf)
    iflag(11) = 10      ! maximum number of refinement iterations

    x(1:n)  = 0.0
    y(1:n)  = 0.0
    b1(1:n) = 0.0
    do i = 1, z
      a1(i) = 0.0
      sn(i) = 0
    end do

    call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, z, ha, iha, &
        b, b1, x, y, aflag, iflag, ifail)

    call check_sp(trim(label), n, x, ifail, 1.0e-4, nfail)
  end subroutine run_test

  ! n-by-n tridiagonal (diag=3, off-diag=-1), solution x=[1,...,1].
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

  ! Check IFAIL=0 and max|x(1:n)-1|<tol; increment nfail on failure.
  subroutine check_sp(label, n, x, ifail, tol, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n, ifail
    real,             intent(in)    :: x(n), tol
    integer,          intent(inout) :: nfail
    real :: err
    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ' ifail=', ifail
      nfail = nfail + 1
      return
    end if
    err = maxval(abs(x(1:n) - 1.0))
    if (err > tol) then
      write(*,'(3a,es10.3)') 'FAIL ', label, ' max_err=', err
      nfail = nfail + 1
    else
      write(*,'(3a,es10.3)') 'PASS ', label, ' max_err=', err
    end if
  end subroutine check_sp

end program test_y12mf_sp
