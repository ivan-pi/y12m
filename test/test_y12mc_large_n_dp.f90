! SPDX-License-Identifier: GPL-2.0-only
!
! test_y12mc_large_n_dp.f90 - Regression test for the 32-bit integer overflow
!   bug in y12mcf (and y12mce).
!
! Bug: when n >= 46341, computing nr = n*n overflows a signed 32-bit integer,
! producing a negative result.  In the Markowitz pivot-selection loop every
! candidate Markowitz cost r3 >= 0 satisfies r3 > nr (negative), so the
! test "if(r3.gt.r)go to 150" always branches away and rcoll is never
! assigned.  The subsequent access ha(rcoll,10) at line 135 of y12mcf.f then
! reads an uninitialised (or garbage) index, causing an out-of-bounds crash.
!
! Fix: initialise nr with huge(nr) (= 2^31-1) instead of n*n.
!
! This test solves a tridiagonal n=50000 system (well above the overflow
! threshold of 46341) via the three-step API y12mb -> y12mc -> y12md using
! double-precision arithmetic and the Markowitz pivot search (iflag(3)=1).
! Without the fix the program aborts with an array-bounds error.
!
! Matrix: tridiagonal with diag=2, off-diag=-1 (1-D Laplacian on n nodes).
! Exact solution: x = [1, 1, ..., 1]^T.
! RHS:           b = A*x = [1, 0, 0, ..., 0, 1]^T.
!
program test_y12mc_large_n_dp
  use y12m
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! n is chosen so that n*n overflows INT32 (overflow threshold is 46341).
  integer, parameter :: n = 50000

  ! Storage: tridiagonal has z0 = 3n-2 entries.  Allocate nn = 5*z0 to
  ! leave ample room for fill-in (none is expected for this ordering, but
  ! the API requires the buffer).
  integer, parameter :: z0  = 3 * n - 2
  integer, parameter :: nn  = 5 * z0
  integer, parameter :: nn1 = nn
  integer, parameter :: iha = n

  real(dp), allocatable :: a(:), pivot(:), b(:)
  integer,  allocatable :: snr(:), rnr(:), ha(:,:)

  real(dp) :: aflag(8)
  integer  :: iflag(10)
  integer  :: z, ifail, i
  real(dp) :: err

  allocate(a(nn), pivot(n), b(n), snr(nn), rnr(nn1), ha(n,11))

  ! ---------------------------------------------------------------
  ! Build the tridiagonal system in COO (row-index, col-index) form.
  ! rnr(k) = row index,  snr(k) = column index,  a(k) = value.
  ! Entries need not be sorted; y12mb will reorder them.
  !
  ! Diagonal entries (value = 2).
  ! ---------------------------------------------------------------
  z = 0
  do i = 1, n
    z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) =  2.0_dp
  end do
  ! Sub-diagonal entries (value = -1).
  do i = 2, n
    z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0_dp
  end do
  ! Super-diagonal entries (value = -1).
  do i = 1, n-1
    z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0_dp
  end do

  ! RHS: b = A * [1,...,1]^T = [1, 0, 0, ..., 0, 1]^T
  b = 0.0_dp
  b(1) = 1.0_dp
  b(n) = 1.0_dp

  ! ---------------------------------------------------------------
  ! Initialise flags (same defaults as y12mae/y12maf).
  ! iflag(3)=1 selects the Markowitz pivot search, which is the code
  ! path that triggered the overflow bug.
  ! ---------------------------------------------------------------
  iflag(1) = 0
  iflag(2) = 3
  iflag(3) = 1
  iflag(4) = 0
  iflag(5) = 1

  aflag(1) = 16.0_dp
  aflag(2) = 1.0d-12
  aflag(3) = 1.0d+16
  aflag(4) = 1.0d-12

  ! ---------------------------------------------------------------
  ! Solve: y12mb (analyse) -> y12mc (factorise) -> y12md (back-sub)
  ! ---------------------------------------------------------------
  call y12mb(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12mc_large_n_dp: y12mb ifail=', ifail
    stop 1
  end if

  call y12mc(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
      aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12mc_large_n_dp: y12mc ifail=', ifail
    stop 1
  end if

  call y12md(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12mc_large_n_dp: y12md ifail=', ifail
    stop 1
  end if

  ! ---------------------------------------------------------------
  ! Verify: solution should be x = [1,...,1]^T.
  ! For n=50000 the tridiagonal 1-D Laplacian is well-conditioned
  ! (condition number ~ n^2 ~ 2.5e9), so errors up to O(n^2*eps)
  ! ~ 5e-7 are expected; a tolerance of 1e-5 is generous.
  ! ---------------------------------------------------------------
  err = maxval(abs(b(1:n) - 1.0_dp))
  if (err > 1.0d-5) then
    write(*,'(a,es12.5)') &
        'FAIL test_y12mc_large_n_dp: max_err=', err
    stop 1
  end if
  write(*,'(a,es12.5)') &
      'PASS test_y12mc_large_n_dp: max_err=', err

end program test_y12mc_large_n_dp
