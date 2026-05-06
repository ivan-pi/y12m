! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot
!
! test_y12ma_hilbert.f90 - Precision comparison on Hilbert systems.
!
! This test avoids fixed absolute error cutoffs.  Instead it compares
! forward and backward errors between single and double precision on the
! same ill-conditioned problem.
!
! It also verifies that larger Hilbert systems trigger y12ma's singularity
! detection path in single precision.
!
program test_y12ma_hilbert
  use y12m, only: y12ma, y12mb, y12mc
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: nfail

  nfail = 0

  ! Compare SP and DP on a moderate ill-conditioned problem.
  call compare_precision(8, nfail)

  ! Exercise singularity detection on larger Hilbert matrices in SP.
  call check_sp_singularity_detection(12, 40, nfail)

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12ma_hilbert tests PASSED'

contains

  subroutine compare_precision(n, nfail)
    integer, intent(in) :: n
    integer, intent(inout) :: nfail

    real(sp) :: ferr_sp, berr_sp
    real(dp) :: ferr_dp, berr_dp
    integer :: ifail_sp, ifail_dp

    call solve_hilbert_sp(n, ifail_sp, ferr_sp, berr_sp)
    call solve_hilbert_dp(n, ifail_dp, ferr_dp, berr_dp)

    if (ifail_sp /= 0) then
      write(*,'(a,i0,a,i0)') 'FAIL compare_precision n=', n, ': SP ifail=', ifail_sp
      nfail = nfail + 1
      return
    end if
    if (ifail_dp /= 0) then
      write(*,'(a,i0,a,i0)') 'FAIL compare_precision n=', n, ': DP ifail=', ifail_dp
      nfail = nfail + 1
      return
    end if

    write(*,'(a,i0,a,es12.5,a,es12.5)') 'INFO n=', n, ' SP ferr=', ferr_sp, ' berr=', berr_sp
    write(*,'(a,i0,a,es12.5,a,es12.5)') 'INFO n=', n, ' DP ferr=', ferr_dp, ' berr=', berr_dp

    if (.not. (ferr_dp < real(ferr_sp, dp))) then
      write(*,'(a,i0)') 'FAIL compare_precision: DP forward error not smaller at n=', n
      nfail = nfail + 1
    else
      write(*,'(a,i0)') 'PASS compare_precision: DP forward error smaller at n=', n
    end if

    if (.not. (berr_dp < real(berr_sp, dp))) then
      write(*,'(a,i0)') 'FAIL compare_precision: DP backward error not smaller at n=', n
      nfail = nfail + 1
    else
      write(*,'(a,i0)') 'PASS compare_precision: DP backward error smaller at n=', n
    end if
  end subroutine compare_precision

  subroutine check_sp_singularity_detection(n_begin, n_end, nfail)
    integer, intent(in) :: n_begin, n_end
    integer, intent(inout) :: nfail

    integer :: n, ifail
    logical :: found

    found = .false.
    do n = n_begin, n_end
      call factorize_hilbert_sp(n, ifail)
      if (ifail /= 0) then
        if (ifail == 3 .or. ifail == 7 .or. ifail == 8) then
          write(*,'(a,i0,a,i0)') 'PASS singularity detection in SP at n=', n, ' ifail=', ifail
          found = .true.
          exit
        end if
        write(*,'(a,i0,a,i0)') 'FAIL unexpected SP ifail at n=', n, ': ifail=', ifail
        nfail = nfail + 1
        return
      end if
    end do

    if (.not. found) then
      write(*,'(a,i0,a,i0,a)') 'FAIL no SP singularity detection observed for Hilbert n=', &
          n_begin, '..', n_end, ' (expected ifail 3/7/8 at some size)'
      nfail = nfail + 1
    end if
  end subroutine check_sp_singularity_detection

  subroutine factorize_hilbert_sp(n, ifail)
    integer, intent(in) :: n
    integer, intent(out) :: ifail

    integer :: nz, nn, nn1, iha, i, j, k
    real(sp), allocatable :: a(:), b(:), pivot(:), aflag(:)
    integer, allocatable :: snr(:), rnr(:), iflag(:), ha(:,:)

    nz = n*n
    nn = 6*nz
    nn1 = 6*nz
    iha = n

    allocate(a(nn), b(n), pivot(n), aflag(8))
    allocate(snr(nn), rnr(nn1), iflag(10), ha(iha,11))

    k = 0
    do i = 1, n
      do j = 1, n
        k = k + 1
        rnr(k) = i
        snr(k) = j
        a(k) = 1.0_sp / real(i + j - 1, sp)
      end do
    end do
    b = 0.0_sp

    iflag(1) = 0
    iflag(2) = 3
    iflag(3) = 1
    iflag(4) = 0
    iflag(5) = 1
    iflag(6:10) = 0

    aflag(1) = 16.0_sp
    aflag(2) = 1.0e-12_sp
    aflag(3) = 1.0e16_sp
    aflag(4) = 1.0e-4_sp
    aflag(5:8) = 0.0_sp
    ha = 0

    call y12mb(n, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail == 0) then
      call y12mc(n, nz, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
    end if

    deallocate(a, b, pivot, aflag, snr, rnr, iflag, ha)
  end subroutine factorize_hilbert_sp

  subroutine solve_hilbert_sp(n, ifail, fwd_err, bwd_err)
    integer, intent(in) :: n
    integer, intent(out) :: ifail
    real(sp), intent(out) :: fwd_err, bwd_err

    integer :: nz, nn, nn1, iha, i, j, k
    real(sp), allocatable :: a(:), pivot(:), b(:), b0(:), aflag(:)
    integer, allocatable :: snr(:), rnr(:), iflag(:), ha(:,:)
    real(sp) :: ax, r_inf, a_inf, x_inf, b_inf, denom

    nz = n*n
    nn = 6*nz
    nn1 = 6*nz
    iha = n

    allocate(a(nn), pivot(n), b(n), b0(n), aflag(8))
    allocate(snr(nn), rnr(nn1), iflag(10), ha(iha,11))

    k = 0
    do i = 1, n
      do j = 1, n
        k = k + 1
        rnr(k) = i
        snr(k) = j
        a(k) = 1.0_sp / real(i + j - 1, sp)
      end do
    end do

    b0 = 0.0_sp
    do k = 1, nz
      b0(rnr(k)) = b0(rnr(k)) + a(k)
    end do
    b = b0

    aflag = 0.0_sp
    ha = 0
    call y12ma(n, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)

    if (ifail /= 0) then
      fwd_err = huge(1.0_sp)
      bwd_err = huge(1.0_sp)
      deallocate(a, pivot, b, b0, aflag, snr, rnr, iflag, ha)
      return
    end if

    fwd_err = maxval(abs(b - 1.0_sp))

    r_inf = 0.0_sp
    a_inf = 0.0_sp
    b_inf = maxval(abs(b0))
    x_inf = maxval(abs(b))
    do i = 1, n
      ax = 0.0_sp
      do j = 1, n
        ax = ax + (1.0_sp / real(i + j - 1, sp)) * b(j)
      end do
      r_inf = max(r_inf, abs(b0(i) - ax))
    end do
    do i = 1, n
      ax = 0.0_sp
      do j = 1, n
        ax = ax + abs(1.0_sp / real(i + j - 1, sp))
      end do
      a_inf = max(a_inf, ax)
    end do

    denom = max(a_inf * x_inf + b_inf, tiny(1.0_sp))
    bwd_err = r_inf / denom

    deallocate(a, pivot, b, b0, aflag, snr, rnr, iflag, ha)
  end subroutine solve_hilbert_sp

  subroutine solve_hilbert_dp(n, ifail, fwd_err, bwd_err)
    integer, intent(in) :: n
    integer, intent(out) :: ifail
    real(dp), intent(out) :: fwd_err, bwd_err

    integer :: nz, nn, nn1, iha, i, j, k
    real(dp), allocatable :: a(:), pivot(:), b(:), b0(:), aflag(:)
    integer, allocatable :: snr(:), rnr(:), iflag(:), ha(:,:)
    real(dp) :: ax, r_inf, a_inf, x_inf, b_inf, denom

    nz = n*n
    nn = 6*nz
    nn1 = 6*nz
    iha = n

    allocate(a(nn), pivot(n), b(n), b0(n), aflag(8))
    allocate(snr(nn), rnr(nn1), iflag(10), ha(iha,11))

    k = 0
    do i = 1, n
      do j = 1, n
        k = k + 1
        rnr(k) = i
        snr(k) = j
        a(k) = 1.0_dp / real(i + j - 1, dp)
      end do
    end do

    b0 = 0.0_dp
    do k = 1, nz
      b0(rnr(k)) = b0(rnr(k)) + a(k)
    end do
    b = b0

    aflag = 0.0_dp
    ha = 0
    call y12ma(n, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)

    if (ifail /= 0) then
      fwd_err = huge(1.0_dp)
      bwd_err = huge(1.0_dp)
      deallocate(a, pivot, b, b0, aflag, snr, rnr, iflag, ha)
      return
    end if

    fwd_err = maxval(abs(b - 1.0_dp))

    r_inf = 0.0_dp
    a_inf = 0.0_dp
    b_inf = maxval(abs(b0))
    x_inf = maxval(abs(b))
    do i = 1, n
      ax = 0.0_dp
      do j = 1, n
        ax = ax + (1.0_dp / real(i + j - 1, dp)) * b(j)
      end do
      r_inf = max(r_inf, abs(b0(i) - ax))
    end do
    do i = 1, n
      ax = 0.0_dp
      do j = 1, n
        ax = ax + abs(1.0_dp / real(i + j - 1, dp))
      end do
      a_inf = max(a_inf, ax)
    end do

    denom = max(a_inf * x_inf + b_inf, tiny(1.0_dp))
    bwd_err = r_inf / denom

    deallocate(a, pivot, b, b0, aflag, snr, rnr, iflag, ha)
  end subroutine solve_hilbert_dp

end program test_y12ma_hilbert
