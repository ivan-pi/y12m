! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot
!
! test_y12mc_large_n_dp.f90 - Regression test for the 32-bit Markowitz
! overflow bug in y12mcf/y12mce.
!
! Previously, both y12mcf and y12mce initialised the Markowitz cost sentinel
! as  nr = n*n  (32-bit integer multiplication).  For n >= 46341 (floor of
! sqrt(2^31-1)) this overflows to a negative value, causing the Markowitz
! search to incorrectly favour every candidate pivot and producing a wrong or
! segfaulting result.
!
! The fix replaces  nr = n*n  with  nr = huge(nr)  which always yields the
! largest representable 32-bit integer regardless of n.
!
! This test constructs a tridiagonal system of size n = 50000 (well above the
! overflow threshold of 46341), solves it with the three-step y12mb + y12mc +
! y12md double-precision API, and verifies:
!   * IFAIL = 0 on every call
!   * max|x(i) - 1| < 1e-10  (all-ones solution)
!
! The system is  A x = b  where
!   A  = tridiag(n; diag=3, off-diag=-1)
!   b(i) = sum of row i = row-sum RHS so that x=[1,...,1] is exact.
!
! Allocatable arrays are used to avoid stack overflow for large n.
!
program test_y12mc_large_n_dp
  use y12m, only: y12mb, y12mc, y12md
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! Problem dimensions
  integer, parameter :: N   = 50000
  integer, parameter :: NZ  = 3*N - 2    ! non-zeros of tridiagonal N×N matrix
  integer, parameter :: NN  = 2*NZ       ! room for LU (iflag(5)=1 -> NN >= 2*NZ)
  integer, parameter :: NN1 = 2*NZ       ! column linked list

  ! Working arrays (allocatable to avoid stack overflow)
  real(dp), allocatable :: a(:), pivot(:), b(:), aflag(:)
  integer,  allocatable :: snr(:), rnr(:), ha(:,:), iflag(:)

  integer :: i, z, ifail

  allocate(a(NN), pivot(N), b(N), aflag(8))
  allocate(snr(NN), rnr(NN1), ha(N,11), iflag(10))

  ! Build tridiagonal n×n (diag=3, off-diag=-1), solution x=[1,...,1].
  z = 0
  do i = 1, N
    z = z + 1 ; snr(z) = i ; rnr(z) = i ; a(z) =  3.0_dp
  end do
  do i = 2, N
    z = z + 1 ; snr(z) = i-1 ; rnr(z) = i ; a(z) = -1.0_dp
  end do
  do i = 1, N-1
    z = z + 1 ; snr(z) = i+1 ; rnr(z) = i ; a(z) = -1.0_dp
  end do
  ! b = row sums: b(1)=b(N)=2, b(i)=1 for interior i
  b(1:N) = 0.0_dp
  do i = 1, z
    b(rnr(i)) = b(rnr(i)) + a(i)
  end do

  ! Flags
  iflag(1)    = 0
  iflag(2)    = 3   ! Markowitz search width; triggers the overflow if unpatched
  iflag(3)    = 1
  iflag(4)    = 0
  iflag(5)    = 1   ! discard L: requires NN >= 2*NZ
  iflag(6:10) = 0

  aflag(1) = 16.0_dp
  aflag(2) = 1.0d-12
  aflag(3) = 1.0d+16
  aflag(4) = 1.0d-12
  aflag(5:8) = 0.0_dp

  ha = 0

  call y12mb(N, z, a, snr, NN, rnr, NN1, ha, N, aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12mc_large_n_dp: y12mb ifail=', ifail
    stop 1
  end if

  call y12mc(N, z, a, snr, NN, rnr, NN1, pivot, b, ha, N, aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12mc_large_n_dp: y12mc ifail=', ifail
    stop 1
  end if

  call y12md(N, a, NN, b, pivot, snr, ha, N, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12mc_large_n_dp: y12md ifail=', ifail
    stop 1
  end if

  ! Verify solution: max|x(i)-1| < 1e-10
  block
    real(dp) :: err
    err = maxval(abs(b(1:N) - 1.0_dp))
    if (err > 1.0d-10) then
      write(*,'(a,es12.5)') &
          'FAIL test_y12mc_large_n_dp: solution max_err=', err
      stop 1
    end if
    write(*,'(a,es12.5)') &
        'PASS test_y12mc_large_n_dp n=50000 tridiag max_err=', err
  end block

  deallocate(a, pivot, b, aflag, snr, rnr, ha, iflag)

end program test_y12mc_large_n_dp
