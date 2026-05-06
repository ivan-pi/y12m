! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot
!
! test_y12ma_hilbert.f90 - SP vs DP accuracy comparison on an ill-conditioned
! 6×6 Hilbert matrix.
!
! The Hilbert matrix H has entries H(i,j) = 1/(i+j-1).  For n=6 its
! condition number is approximately 1.5×10^7.  This means:
!
!   * A double-precision solve (machine epsilon ~2e-16) is expected to produce a
!     solution with a relative error around 1e-9 — comfortably accurate.
!
!   * A single-precision solve (machine epsilon ~1.2e-7) is expected to produce a
!     solution with relative error around 1.0 — essentially inaccurate.
!
! This test verifies:
!   * Both solves return IFAIL = 0.
!   * DP max|x - 1| < 1e-5  (tight tolerance; easily satisfied at double precision).
!   * SP max|x - 1| > 1e-2  (lower bound showing single precision is inadequate).
!   * DP error < SP error    (DP is strictly more accurate for this problem).
!
! Exact solution: x = [1, 1, 1, 1, 1, 1]^T.
! Right-hand side: b(i) = sum_{j=1}^{6} 1/(i+j-1)  (row sums of H).
!
program test_y12ma_hilbert
  use y12m, only: y12ma
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: N    = 6
  integer, parameter :: NZ   = N * N          ! 6x6 dense matrix: 36 non-zeros
  integer, parameter :: NMAX = N + 4          ! leading dimension for pivot/ha (>= N)
  integer, parameter :: NN   = 4 * NZ         ! workspace: y12ma uses iflag(5)=1 -> NN >= 2*NZ
  integer, parameter :: NN1  = 4 * NZ

  ! Double-precision working arrays
  double precision :: a_dp(NN), pivot_dp(NMAX), b_dp(NMAX), aflag_dp(8)
  integer          :: snr_dp(NN), rnr_dp(NN1), ha_dp(NMAX,11), iflag_dp(10)

  ! Single-precision working arrays
  real    :: a_sp(NN), pivot_sp(NMAX), b_sp(NMAX), aflag_sp(8)
  integer :: snr_sp(NN), rnr_sp(NN1), ha_sp(NMAX,11), iflag_sp(10)

  double precision :: err_dp
  real             :: err_sp

  integer :: i, j, z, ifail_dp, ifail_sp, nfail

  nfail = 0

  ! =========================================================================
  ! Build the 6x6 Hilbert matrix and RHS in double precision
  ! =========================================================================
  z = 0
  do i = 1, N
    do j = 1, N
      z = z + 1
      rnr_dp(z) = i
      snr_dp(z) = j
      a_dp(z)   = 1.0d0 / dble(i + j - 1)
    end do
  end do
  b_dp(1:N) = 0.0d0
  do z = 1, NZ
    b_dp(rnr_dp(z)) = b_dp(rnr_dp(z)) + a_dp(z)
  end do

  ! =========================================================================
  ! Double-precision solve
  ! =========================================================================
  aflag_dp = 0.0d0
  ha_dp    = 0
  call y12ma(N, NZ, a_dp, snr_dp, NN, rnr_dp, NN1, pivot_dp, ha_dp, NMAX, &
      aflag_dp, iflag_dp, b_dp, ifail_dp)

  if (ifail_dp /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12ma_hilbert dp: y12ma ifail=', ifail_dp
    nfail = nfail + 1
    err_dp = huge(1.0d0)
  else
    err_dp = maxval(abs(b_dp(1:N) - 1.0d0))
    if (err_dp > 1.0d-5) then
      write(*,'(a,es12.5)') &
          'FAIL test_y12ma_hilbert dp: max_err > 1e-5, got ', err_dp
      nfail = nfail + 1
    else
      write(*,'(a,es12.5)') &
          'PASS test_y12ma_hilbert dp: max|x-1|=', err_dp
    end if
  end if

  ! =========================================================================
  ! Build the same matrix and RHS in single precision
  ! =========================================================================
  z = 0
  do i = 1, N
    do j = 1, N
      z = z + 1
      rnr_sp(z) = i
      snr_sp(z) = j
      a_sp(z)   = 1.0e0 / real(i + j - 1)
    end do
  end do
  b_sp(1:N) = 0.0e0
  do z = 1, NZ
    b_sp(rnr_sp(z)) = b_sp(rnr_sp(z)) + a_sp(z)
  end do

  ! =========================================================================
  ! Single-precision solve
  ! =========================================================================
  aflag_sp = 0.0e0
  ha_sp    = 0
  call y12ma(N, NZ, a_sp, snr_sp, NN, rnr_sp, NN1, pivot_sp, ha_sp, NMAX, &
      aflag_sp, iflag_sp, b_sp, ifail_sp)

  if (ifail_sp /= 0) then
    write(*,'(a,i0)') 'FAIL test_y12ma_hilbert sp: y12ma ifail=', ifail_sp
    nfail = nfail + 1
    err_sp = huge(1.0e0)
  else
    err_sp = maxval(abs(b_sp(1:N) - 1.0e0))
    write(*,'(a,es12.5)') &
        'INFO test_y12ma_hilbert sp: max|x-1|=', err_sp
  end if

  ! =========================================================================
  ! Compare SP vs DP accuracy
  ! =========================================================================
  if (ifail_dp == 0 .and. ifail_sp == 0) then

    ! SP solution should be demonstrably inaccurate for this ill-conditioned system
    if (err_sp <= 1.0e-2) then
      write(*,'(a,es12.5)') &
          'FAIL test_y12ma_hilbert sp_vs_dp: SP error unexpectedly small, got ', err_sp
      nfail = nfail + 1
    else
      write(*,'(a,es12.5)') &
          'PASS test_y12ma_hilbert sp_vs_dp: SP error large as expected, got ', err_sp
    end if

    ! DP solution should be strictly more accurate than SP
    if (dble(err_sp) > err_dp) then
      write(*,'(a,es12.5,a,es12.5)') &
          'PASS test_y12ma_hilbert precision_comparison: dp_err=', err_dp, &
          ' < sp_err=', err_sp
    else
      write(*,'(a,es12.5,a,es12.5)') &
          'FAIL test_y12ma_hilbert precision_comparison: dp_err=', err_dp, &
          ' not less than sp_err=', err_sp
      nfail = nfail + 1
    end if

  end if

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_hilbert tests PASSED'

end program test_y12ma_hilbert
