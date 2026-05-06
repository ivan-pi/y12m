! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot
!
! test_y12ma_hilbert.f90 - Precision comparison on Hilbert systems.
!
! This test avoids fixed absolute error cutoffs.  Instead it compares
! forward and backward errors between single and double precision on the
! same ill-conditioned problem.
!
! Hilbert matrix H has entries H(i,j)=1/(i+j-1).  It is a classic dense
! matrix that quickly becomes ill-conditioned as n grows.  We intentionally
! feed this dense matrix to y12m's sparse-direct path so we can stress
! precision behavior and singularity detection in a controlled way.
!
! It also verifies that larger Hilbert systems trigger y12ma's singularity
! detection path in single precision.
!
program test_y12ma_hilbert
  use y12m_test_helpers, only: forward_error, backward_error
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

    iflag(1) = 0  ! fresh start
    iflag(2) = 3  ! Markowitz search width
    iflag(3) = 1  ! general Markowitz pivoting
    iflag(4) = 0  ! no structure reuse
    iflag(5) = 1  ! discard L factors after factorization
    iflag(6:10) = 0

    aflag(1) = 16.0_sp     ! stability factor
    aflag(2) = 1.0e-12_sp  ! drop tolerance
    aflag(3) = 1.0e16_sp   ! growth threshold
    aflag(4) = 1.0e-4_sp   ! looser singularity threshold to trigger detection on larger Hilbert n
    aflag(5:8) = 0.0_sp
    ha = 0

    call y12mb(n, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
    if (ifail /= 0) return
    call y12mc(n, nz, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
  end subroutine factorize_hilbert_sp

  subroutine solve_hilbert_sp(n, ifail, fwd_err, bwd_err)
    integer, intent(in) :: n
    integer, intent(out) :: ifail
    real(sp), intent(out) :: fwd_err, bwd_err

    integer :: nz, nn, nn1, iha, i, j, k
    real(sp), allocatable :: a(:), a0(:), pivot(:), b(:), b0(:), aflag(:)
    integer, allocatable :: snr(:), rnr(:), snr0(:), rnr0(:), iflag(:), ha(:,:)
    real(sp) :: x_true(n)

    nz = n*n
    nn = 6*nz
    nn1 = 6*nz
    iha = n

    ! Keep copies of the original COO matrix: y12ma overwrites a/snr/rnr with LU data.
    allocate(a(nn), a0(nz), pivot(n), b(n), b0(n), aflag(8))
    allocate(snr(nn), rnr(nn1), snr0(nz), rnr0(nz), iflag(10), ha(iha,11))

    k = 0
    do i = 1, n
      do j = 1, n
        k = k + 1
        rnr(k) = i
        snr(k) = j
        a(k) = 1.0_sp / real(i + j - 1, sp)
        rnr0(k) = rnr(k)
        snr0(k) = snr(k)
        a0(k) = a(k)
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
      return
    end if

    ! Forward error: max componentwise distance from exact x=[1,...,1].
    x_true = 1.0_sp
    call forward_error(n, b, x_true, fwd_err)

    ! Backward error: smallest relative perturbation that explains residual.
    call backward_error(n, nz, a0, rnr0, snr0, b, b0, bwd_err)
  end subroutine solve_hilbert_sp

  subroutine solve_hilbert_dp(n, ifail, fwd_err, bwd_err)
    integer, intent(in) :: n
    integer, intent(out) :: ifail
    real(dp), intent(out) :: fwd_err, bwd_err

    integer :: nz, nn, nn1, iha, i, j, k
    real(dp), allocatable :: a(:), a0(:), pivot(:), b(:), b0(:), aflag(:)
    integer, allocatable :: snr(:), rnr(:), snr0(:), rnr0(:), iflag(:), ha(:,:)
    real(dp) :: x_true(n)

    nz = n*n
    nn = 6*nz
    nn1 = 6*nz
    iha = n

    ! Keep copies of the original COO matrix: y12ma overwrites a/snr/rnr with LU data.
    allocate(a(nn), a0(nz), pivot(n), b(n), b0(n), aflag(8))
    allocate(snr(nn), rnr(nn1), snr0(nz), rnr0(nz), iflag(10), ha(iha,11))

    k = 0
    do i = 1, n
      do j = 1, n
        k = k + 1
        rnr(k) = i
        snr(k) = j
        a(k) = 1.0_dp / real(i + j - 1, dp)
        rnr0(k) = rnr(k)
        snr0(k) = snr(k)
        a0(k) = a(k)
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
      return
    end if

    ! Forward error: max componentwise distance from exact x=[1,...,1].
    x_true = 1.0_dp
    call forward_error(n, b, x_true, fwd_err)

    ! Backward error: smallest relative perturbation that explains residual.
    call backward_error(n, nz, a0, rnr0, snr0, b, b0, bwd_err)
  end subroutine solve_hilbert_dp

end program test_y12ma_hilbert
