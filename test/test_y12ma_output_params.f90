! SPDX-License-Identifier: GPL-2.0-only
!
! test_y12ma_output_params.f90 - verify all documented output-parameter
! invariants produced by y12ma (y12mae for single precision, y12maf for
! double precision) after a successful call on a 4-by-4 tridiagonal system.
!
! The invariants checked are:
!
!   AFLAG(1) = 16.0            stability factor set by y12ma
!   AFLAG(2) = 1e-12           drop tolerance set by y12ma
!   AFLAG(3) = 1e16            growth limit set by y12ma
!                              (the doc says 1e6 but the code sets 1e16)
!   AFLAG(4) = 1e-12           pivot threshold set by y12ma
!   AFLAG(5) = AFLAG(7)/AFLAG(6)  growth factor definition (y12mc)
!   AFLAG(5) in [1, AFLAG(3)]  bounds on growth factor
!   AFLAG(6) = 3.0             max |element| in input matrix (y12mb)
!   AFLAG(7) >= AFLAG(6)       max during LU >= max in original A
!   AFLAG(8) > 0               smallest pivot magnitude > 0
!   AFLAG(8) == minval(|PIVOT|) link between AFLAG(8) and PIVOT array
!   IFLAG(2) = 2               pivotal-search rows set by y12ma
!                              (note: doc says 3, code sets 2)
!   IFLAG(3) = 1               general pivotal strategy
!   IFLAG(4) = 0               single system mode
!   IFLAG(5) = 1               discard L factors after factorization
!   IFLAG(6) >= 0              row garbage-collection count
!   IFLAG(7) >= 0              column garbage-collection count
!   IFLAG(8) in [z, NN]        max nnz stored during elimination
!   PIVOT(i) all nonzero       nonsingular matrix => all pivots > 0
!   RNR(1:z) all zero          y12mc zeroes out the first z row entries
!
! Both the single-precision (y12mae) and double-precision (y12maf) paths
! are exercised.
!
program test_y12ma_output_params
  use y12m, only: y12ma
  implicit none

  integer, parameter :: NMAX = 10
  integer, parameter :: NNP  = 200
  integer, parameter :: NN1P = 200

  integer :: nfail

  nfail = 0
  call test_sp('y12ma_output SP n=4 tridiag', 4, nfail)
  call test_dp('y12ma_output DP n=4 tridiag', 4, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_output_params tests PASSED'

contains

  ! -----------------------------------------------------------------------
  ! Single-precision variant
  ! -----------------------------------------------------------------------
  subroutine test_sp(label, n, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    integer,          intent(inout) :: nfail

    real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
    integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
    integer :: z, ifail, i
    real    :: min_pivot, ratio

    call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)

    call y12ma(n, z, a, snr, NNP, rnr, NN1P, pivot, ha, NMAX, &
        aflag, iflag, b, ifail)

    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': solver failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    ! ---- AFLAG output checks ----
    call chkr(label, 'AFLAG(1)==16.0',   aflag(1), 16.0,   1.0e-5, nfail)
    call chkr(label, 'AFLAG(2)==1e-12',  aflag(2), 1.0e-12, 1.0e-5, nfail)
    call chkr(label, 'AFLAG(3)==1e16',   aflag(3), 1.0e16,  1.0e-5, nfail)
    call chkr(label, 'AFLAG(4)==1e-12',  aflag(4), 1.0e-12, 1.0e-5, nfail)

    if (aflag(6) > 0.0) then
      ratio = aflag(7) / aflag(6)
      call chkr(label, 'AFLAG(5)==AFLAG(7)/AFLAG(6)', aflag(5), ratio, 1.0e-4, nfail)
    else
      write(*,'(3a)') 'FAIL ', label, ': AFLAG(6) <= 0'
      nfail = nfail + 1
    end if

    if (aflag(5) < 1.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(5) < 1.0, got ', aflag(5)
      nfail = nfail + 1
    end if
    if (aflag(5) > aflag(3)) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(5) > AFLAG(3), got ', aflag(5)
      nfail = nfail + 1
    end if

    ! max element in tridiag(3,-1) is 3.0
    call chkr(label, 'AFLAG(6)==3.0', aflag(6), 3.0, 1.0e-5, nfail)

    if (aflag(7) < aflag(6) * (1.0 - 1.0e-5)) then
      write(*,'(3a,2(1pg12.5,1x))') 'FAIL ', label, &
          ': AFLAG(7) < AFLAG(6), got/exp ', aflag(7), aflag(6)
      nfail = nfail + 1
    end if

    if (aflag(8) <= 0.0) then
      write(*,'(3a,1pg12.5)') 'FAIL ', label, ': AFLAG(8) <= 0, got ', aflag(8)
      nfail = nfail + 1
    end if

    min_pivot = minval(abs(pivot(1:n)))
    call chkr(label, 'AFLAG(8)==min|PIVOT|', aflag(8), min_pivot, 1.0e-4, nfail)

    ! ---- IFLAG output checks ----
    if (iflag(2) /= 2) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(2)/=2, got ', iflag(2)
      nfail = nfail + 1
    end if
    if (iflag(3) /= 1) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(3)/=1, got ', iflag(3)
      nfail = nfail + 1
    end if
    if (iflag(4) /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(4)/=0, got ', iflag(4)
      nfail = nfail + 1
    end if
    if (iflag(5) /= 1) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(5)/=1, got ', iflag(5)
      nfail = nfail + 1
    end if
    if (iflag(6) < 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(6) < 0, got ', iflag(6)
      nfail = nfail + 1
    end if
    if (iflag(7) < 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(7) < 0, got ', iflag(7)
      nfail = nfail + 1
    end if
    if (iflag(8) < z) then
      write(*,'(3a,2(i0,1x))') 'FAIL ', label, ': IFLAG(8) < z, got/z=', iflag(8), z
      nfail = nfail + 1
    end if
    if (iflag(8) > NNP) then
      write(*,'(3a,2(i0,1x))') 'FAIL ', label, ': IFLAG(8) > NN, got/NN=', iflag(8), NNP
      nfail = nfail + 1
    end if

    ! ---- PIVOT nonzero for nonsingular matrix ----
    do i = 1, n
      if (pivot(i) == 0.0) then
        write(*,'(3a,i0)') 'FAIL ', label, ': PIVOT(i)==0 at i=', i
        nfail = nfail + 1
      end if
    end do

    write(*,'(3a)') 'PASS ', label, ': all SP output invariants hold'
  end subroutine test_sp

  ! -----------------------------------------------------------------------
  ! Double-precision variant
  ! -----------------------------------------------------------------------
  subroutine test_dp(label, n, nfail)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    integer,          intent(inout) :: nfail

    double precision :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
    integer          :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
    integer          :: z, ifail, i
    double precision :: min_pivot, ratio

    call build_tridiag_dp(n, NNP, NN1P, z, a, snr, rnr, b)

    call y12ma(n, z, a, snr, NNP, rnr, NN1P, pivot, ha, NMAX, &
        aflag, iflag, b, ifail)

    if (ifail /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': solver failed ifail=', ifail
      nfail = nfail + 1
      return
    end if

    ! ---- AFLAG output checks ----
    call chkd(label, 'AFLAG(1)==16.0',  aflag(1), 16.0d0,  1.0d-10, nfail)
    call chkd(label, 'AFLAG(2)==1e-12', aflag(2), 1.0d-12, 1.0d-10, nfail)
    call chkd(label, 'AFLAG(3)==1e16',  aflag(3), 1.0d16,  1.0d-10, nfail)
    call chkd(label, 'AFLAG(4)==1e-12', aflag(4), 1.0d-12, 1.0d-10, nfail)

    if (aflag(6) > 0.0d0) then
      ratio = aflag(7) / aflag(6)
      call chkd(label, 'AFLAG(5)==AFLAG(7)/AFLAG(6)', aflag(5), ratio, 1.0d-8, nfail)
    else
      write(*,'(3a)') 'FAIL ', label, ': AFLAG(6) <= 0'
      nfail = nfail + 1
    end if

    if (aflag(5) < 1.0d0) then
      write(*,'(3a,1pg16.9)') 'FAIL ', label, ': AFLAG(5) < 1.0, got ', aflag(5)
      nfail = nfail + 1
    end if
    if (aflag(5) > aflag(3)) then
      write(*,'(3a,1pg16.9)') 'FAIL ', label, ': AFLAG(5) > AFLAG(3), got ', aflag(5)
      nfail = nfail + 1
    end if

    call chkd(label, 'AFLAG(6)==3.0', aflag(6), 3.0d0, 1.0d-10, nfail)

    if (aflag(7) < aflag(6) * (1.0d0 - 1.0d-10)) then
      write(*,'(3a,2(1pg16.9,1x))') 'FAIL ', label, &
          ': AFLAG(7) < AFLAG(6), got/exp ', aflag(7), aflag(6)
      nfail = nfail + 1
    end if

    if (aflag(8) <= 0.0d0) then
      write(*,'(3a,1pg16.9)') 'FAIL ', label, ': AFLAG(8) <= 0, got ', aflag(8)
      nfail = nfail + 1
    end if

    min_pivot = minval(abs(pivot(1:n)))
    call chkd(label, 'AFLAG(8)==min|PIVOT|', aflag(8), min_pivot, 1.0d-8, nfail)

    ! ---- IFLAG output checks ----
    if (iflag(2) /= 2) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(2)/=2, got ', iflag(2)
      nfail = nfail + 1
    end if
    if (iflag(3) /= 1) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(3)/=1, got ', iflag(3)
      nfail = nfail + 1
    end if
    if (iflag(4) /= 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(4)/=0, got ', iflag(4)
      nfail = nfail + 1
    end if
    if (iflag(5) /= 1) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(5)/=1, got ', iflag(5)
      nfail = nfail + 1
    end if
    if (iflag(6) < 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(6) < 0, got ', iflag(6)
      nfail = nfail + 1
    end if
    if (iflag(7) < 0) then
      write(*,'(3a,i0)') 'FAIL ', label, ': IFLAG(7) < 0, got ', iflag(7)
      nfail = nfail + 1
    end if
    if (iflag(8) < z) then
      write(*,'(3a,2(i0,1x))') 'FAIL ', label, ': IFLAG(8) < z, got/z=', iflag(8), z
      nfail = nfail + 1
    end if
    if (iflag(8) > NNP) then
      write(*,'(3a,2(i0,1x))') 'FAIL ', label, ': IFLAG(8) > NN, got/NN=', iflag(8), NNP
      nfail = nfail + 1
    end if

    do i = 1, n
      if (pivot(i) == 0.0d0) then
        write(*,'(3a,i0)') 'FAIL ', label, ': PIVOT(i)==0 at i=', i
        nfail = nfail + 1
      end if
    end do

    write(*,'(3a)') 'PASS ', label, ': all DP output invariants hold'
  end subroutine test_dp

  ! -----------------------------------------------------------------------
  ! Relative-error checkers
  ! -----------------------------------------------------------------------
  subroutine chkr(label, param, actual, expected, rtol, nfail)
    character(len=*), intent(in)    :: label, param
    real,             intent(in)    :: actual, expected, rtol
    integer,          intent(inout) :: nfail
    real :: err
    err = abs(actual - expected) / max(abs(expected), tiny(1.0))
    if (err > rtol) then
      write(*,'(5a,2(1pg12.5,1x))') 'FAIL ', label, ' ', param, &
          ': got/exp=', actual, expected
      nfail = nfail + 1
    end if
  end subroutine chkr

  subroutine chkd(label, param, actual, expected, rtol, nfail)
    character(len=*), intent(in)    :: label, param
    double precision, intent(in)    :: actual, expected, rtol
    integer,          intent(inout) :: nfail
    double precision :: err
    err = abs(actual - expected) / max(abs(expected), tiny(1.0d0))
    if (err > rtol) then
      write(*,'(5a,2(1pg16.9,1x))') 'FAIL ', label, ' ', param, &
          ': got/exp=', actual, expected
      nfail = nfail + 1
    end if
  end subroutine chkd

  ! -----------------------------------------------------------------------
  ! Matrix builders: n-by-n tridiagonal (diag=3, off-diag=-1), x=[1..1]
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

  subroutine build_tridiag_dp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer,          intent(in)  :: n, nnmax, nn1max
    integer,          intent(out) :: z
    double precision, intent(out) :: a(nnmax), b(n)
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
    b(1:n) = 0.0d0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_tridiag_dp

end program test_y12ma_output_params
