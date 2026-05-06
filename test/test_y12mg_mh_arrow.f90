! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot
!
! test_y12mg_mh_arrow.f90 - y12mh (1-norm) and y12mg (condition number)
! with non-tridiagonal sparse matrices.
!
! This test extends the coverage of test_y12mg_mh.f90 (which only uses
! tridiagonal matrices) by applying y12mh and the full
! y12mh -> y12mb -> y12mc -> y12md -> y12mg sequence to:
!
!   (a) A 6x6 arrow matrix (full first row and column, diagonal remainder).
!       The arrow matrix used here is:
!         row 1: all-ones [1, 1, 1, 1, 1, 1]
!         row i (i>=2): col 1 = 1, col i = i+1
!
!       Exact 1-norm (max absolute column sum):
!         col 1: sum = 1 (row 1) + 5*1 (rows 2-6) = 6
!         col j (j>=2): sum = 1 (row 1) + (j+1) (row j, diagonal) = j+2
!         For n=6: col 6 sum = 6+2 = 8
!         1-norm = max(6, 4, 5, 6, 7, 8) = 8
!
!   (b) A 8x8 arrow matrix (double precision).
!       1-norm for n=8: col 8 sum = 8+2 = 10  =>  1-norm = 10.
!
!   (c) A 10x10 pentadiagonal matrix (double precision).
!       Pentadiagonal: diag=5, dist-1 off-diag=-1, dist-2 off-diag=-0.5.
!       Column j sum = |a(j,j)| + |a(j-1,j)| + |a(j+1,j)|
!                      + |a(j-2,j)| + |a(j+2,j)|
!       Interior columns (3 <= j <= 8) contribute all five bands:
!         5 + 1 + 1 + 0.5 + 0.5 = 8
!       1-norm = 8
!
! The RCOND condition-number estimate is checked to be in (0, 1].
!
program test_y12mg_mh_arrow
  use y12m
  use y12m_test_helpers, only: build_arrow, build_pentadiag, check_solution
  implicit none

  integer, parameter :: NMAX = 12
  integer, parameter :: NNP  = 600
  integer, parameter :: NN1P = 600

  ! Single-precision working arrays
  real    :: a_sp(NNP), pivot_sp(NMAX), b_sp(NMAX), aflag_sp(8)
  real    :: work_sp(NMAX), anorm_sp, rcond_sp
  integer :: snr_sp(NNP), rnr_sp(NN1P), ha_sp(NMAX,11), iflag_sp(10)

  ! Double-precision working arrays
  double precision :: a_dp(NNP), pivot_dp(NMAX), b_dp(NMAX), aflag_dp(8)
  double precision :: work_dp(NMAX), anorm_dp, rcond_dp
  integer          :: snr_dp(NNP), rnr_dp(NN1P), ha_dp(NMAX,11), iflag_dp(10)

  integer :: n, z, nn, nn1, iha, ifail, nfail
  real    :: expected_norm_sp

  nfail = 0

  ! =====================================================================
  ! (a) y12mh on a 6x6 arrow matrix (single precision)
  !
  ! col 1 sum = n = 6
  ! col j (j>=2) sum = 1 + (j+1) = j+2
  ! For n=6: max(6, 4, 5, 6, 7, 8) = 8  ->  1-norm = 8
  ! =====================================================================
  n = 6
  expected_norm_sp = real(n + 2)   ! col n has the largest sum: 1 + (n+1) = n+2
  call build_arrow(n, NNP, NN1P, z, a_sp, snr_sp, rnr_sp, b_sp)
  call y12mh(n, z, a_sp, snr_sp, work_sp, anorm_sp)
  if (abs(anorm_sp - expected_norm_sp) > 1.0e-5) then
    write(*,'(a,f8.3,a,f8.3)') &
        'FAIL y12mh n=6 arrow: expected anorm=', expected_norm_sp, ' got ', anorm_sp
    nfail = nfail + 1
  else
    write(*,'(a,f8.3)') 'PASS y12mh n=6 arrow: anorm=', anorm_sp
  end if

  ! =====================================================================
  ! (b) y12mh on an 8x8 arrow matrix (double precision)
  !
  ! col 1 sum = n = 8
  ! col j (j>=2) sum = j+2
  ! For n=8: max(8, 4, 5, 6, 7, 8, 9, 10) = 10  ->  1-norm = 10
  ! =====================================================================
  n = 8
  call build_arrow(n, NNP, NN1P, z, a_dp, snr_dp, rnr_dp, b_dp)
  call y12mh(n, z, a_dp, snr_dp, work_dp, anorm_dp)
  if (abs(anorm_dp - 10.0d0) > 1.0d-10) then
    write(*,'(a,f10.5)') &
        'FAIL y12mh n=8 arrow: expected anorm=10.0 got ', anorm_dp
    nfail = nfail + 1
  else
    write(*,'(a,f10.5)') 'PASS y12mh n=8 arrow: anorm=', anorm_dp
  end if

  ! =====================================================================
  ! (c) y12mh on a 10x10 pentadiagonal matrix (double precision)
  !
  ! Band structure: main=5, dist1=-1, dist2=-0.5
  ! Interior col j (3<=j<=n-2) absolute column sum:
  !   |-0.5| + |-1| + |5| + |-1| + |-0.5| = 8
  ! Corner cols (1,2,n-1,n) have fewer terms but are <= 8.
  ! 1-norm = 8
  ! =====================================================================
  n = 10
  call build_pentadiag(n, NNP, NN1P, z, a_dp, snr_dp, rnr_dp, b_dp)
  call y12mh(n, z, a_dp, snr_dp, work_dp, anorm_dp)
  if (abs(anorm_dp - 8.0d0) > 1.0d-10) then
    write(*,'(a,f10.5)') &
        'FAIL y12mh n=10 pentadiag: expected anorm=8.0 got ', anorm_dp
    nfail = nfail + 1
  else
    write(*,'(a,f10.5)') 'PASS y12mh n=10 pentadiag: anorm=', anorm_dp
  end if

  ! =====================================================================
  ! y12mg condition-number estimate on a 6x6 arrow matrix (single precision)
  !
  ! Sequence: y12mh -> y12mb -> y12mc(iflag(5)=2) -> y12md -> y12mg
  ! IFLAG(5)=2 is required by y12mg (L factors must be preserved).
  ! Pass criterion: RCOND in (0, 1].
  ! =====================================================================
  n = 6
  call build_arrow(n, NNP, NN1P, z, a_sp, snr_sp, rnr_sp, b_sp)
  nn = NNP ; nn1 = NN1P ; iha = NMAX

  call y12mh(n, z, a_sp, snr_sp, work_sp, anorm_sp)

  iflag_sp(1) = 0 ; iflag_sp(2) = 3 ; iflag_sp(3) = 1
  iflag_sp(4) = 0 ; iflag_sp(5) = 2   ! keep L for y12mg
  aflag_sp(1) = 16.0  ; aflag_sp(2) = 1.0e-12
  aflag_sp(3) = 1.0e+16 ; aflag_sp(4) = 1.0e-12

  call y12mb(n, z, a_sp, snr_sp, nn, rnr_sp, nn1, ha_sp, iha, &
      aflag_sp, iflag_sp, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL y12mg n=6 arrow (sp): y12mb ifail=', ifail
    nfail = nfail + 1
  else
    call y12mc(n, z, a_sp, snr_sp, nn, rnr_sp, nn1, pivot_sp, b_sp, &
        ha_sp, iha, aflag_sp, iflag_sp, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL y12mg n=6 arrow (sp): y12mc ifail=', ifail
      nfail = nfail + 1
    else
      call y12md(n, a_sp, nn, b_sp, pivot_sp, snr_sp, ha_sp, iha, &
          iflag_sp, ifail)
      if (ifail /= 0) then
        write(*,'(a,i0)') 'FAIL y12mg n=6 arrow (sp): y12md ifail=', ifail
        nfail = nfail + 1
      else
        call check_solution('y12mg_mh_arrow n=6 arrow sp', n, b_sp, ifail, &
            1.0e-4, nfail)
        ifail = 0
        call y12mg(n, nn, a_sp, snr_sp, work_sp, pivot_sp, anorm_sp, &
            rcond_sp, iha, ha_sp, iflag_sp, ifail)
        if (ifail /= 0) then
          write(*,'(a,i0)') 'FAIL y12mg n=6 arrow (sp): y12mg ifail=', ifail
          nfail = nfail + 1
        else if (rcond_sp <= 0.0 .or. rcond_sp > 1.0) then
          write(*,'(a,es10.3)') &
              'FAIL y12mg n=6 arrow (sp): rcond out of (0,1]: ', rcond_sp
          nfail = nfail + 1
        else
          write(*,'(a,es10.3,a,es10.3)') &
              'PASS y12mg n=6 arrow (sp): rcond=', rcond_sp, &
              ' cond=', 1.0 / rcond_sp
        end if
      end if
    end if
  end if

  ! =====================================================================
  ! y12mg condition-number estimate on an 8x8 arrow matrix (double precision)
  ! =====================================================================
  n = 8
  call build_arrow(n, NNP, NN1P, z, a_dp, snr_dp, rnr_dp, b_dp)
  nn = NNP ; nn1 = NN1P ; iha = NMAX

  call y12mh(n, z, a_dp, snr_dp, work_dp, anorm_dp)

  iflag_dp(1) = 0 ; iflag_dp(2) = 3 ; iflag_dp(3) = 1
  iflag_dp(4) = 0 ; iflag_dp(5) = 2
  aflag_dp(1) = 16.0d0 ; aflag_dp(2) = 1.0d-12
  aflag_dp(3) = 1.0d+16 ; aflag_dp(4) = 1.0d-12

  call y12mb(n, z, a_dp, snr_dp, nn, rnr_dp, nn1, ha_dp, iha, &
      aflag_dp, iflag_dp, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL y12mg n=8 arrow (dp): y12mb ifail=', ifail
    nfail = nfail + 1
  else
    call y12mc(n, z, a_dp, snr_dp, nn, rnr_dp, nn1, pivot_dp, b_dp, &
        ha_dp, iha, aflag_dp, iflag_dp, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL y12mg n=8 arrow (dp): y12mc ifail=', ifail
      nfail = nfail + 1
    else
      call y12md(n, a_dp, nn, b_dp, pivot_dp, snr_dp, ha_dp, iha, &
          iflag_dp, ifail)
      if (ifail /= 0) then
        write(*,'(a,i0)') 'FAIL y12mg n=8 arrow (dp): y12md ifail=', ifail
        nfail = nfail + 1
      else
        call check_solution('y12mg_mh_arrow n=8 arrow dp', n, b_dp, ifail, &
            1.0d-10, nfail)
        ifail = 0
        call y12mg(n, nn, a_dp, snr_dp, work_dp, pivot_dp, anorm_dp, &
            rcond_dp, iha, ha_dp, iflag_dp, ifail)
        if (ifail /= 0) then
          write(*,'(a,i0)') 'FAIL y12mg n=8 arrow (dp): y12mg ifail=', ifail
          nfail = nfail + 1
        else if (rcond_dp <= 0.0d0 .or. rcond_dp > 1.0d0) then
          write(*,'(a,es12.5)') &
              'FAIL y12mg n=8 arrow (dp): rcond out of (0,1]: ', rcond_dp
          nfail = nfail + 1
        else
          write(*,'(a,es12.5,a,es12.5)') &
              'PASS y12mg n=8 arrow (dp): rcond=', rcond_dp, &
              ' cond=', 1.0d0 / rcond_dp
        end if
      end if
    end if
  end if

  ! =====================================================================
  ! y12mg condition-number estimate on a 10x10 pentadiagonal (double precision)
  ! =====================================================================
  n = 10
  call build_pentadiag(n, NNP, NN1P, z, a_dp, snr_dp, rnr_dp, b_dp)
  nn = NNP ; nn1 = NN1P ; iha = NMAX

  call y12mh(n, z, a_dp, snr_dp, work_dp, anorm_dp)

  iflag_dp(1) = 0 ; iflag_dp(2) = 3 ; iflag_dp(3) = 1
  iflag_dp(4) = 0 ; iflag_dp(5) = 2
  aflag_dp(1) = 16.0d0 ; aflag_dp(2) = 1.0d-12
  aflag_dp(3) = 1.0d+16 ; aflag_dp(4) = 1.0d-12

  call y12mb(n, z, a_dp, snr_dp, nn, rnr_dp, nn1, ha_dp, iha, &
      aflag_dp, iflag_dp, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL y12mg n=10 pentadiag (dp): y12mb ifail=', ifail
    nfail = nfail + 1
  else
    call y12mc(n, z, a_dp, snr_dp, nn, rnr_dp, nn1, pivot_dp, b_dp, &
        ha_dp, iha, aflag_dp, iflag_dp, ifail)
    if (ifail /= 0) then
      write(*,'(a,i0)') 'FAIL y12mg n=10 pentadiag (dp): y12mc ifail=', ifail
      nfail = nfail + 1
    else
      call y12md(n, a_dp, nn, b_dp, pivot_dp, snr_dp, ha_dp, iha, &
          iflag_dp, ifail)
      if (ifail /= 0) then
        write(*,'(a,i0)') 'FAIL y12mg n=10 pentadiag (dp): y12md ifail=', ifail
        nfail = nfail + 1
      else
        call check_solution('y12mg_mh_arrow n=10 pentadiag dp', n, b_dp, &
            ifail, 1.0d-10, nfail)
        ifail = 0
        call y12mg(n, nn, a_dp, snr_dp, work_dp, pivot_dp, anorm_dp, &
            rcond_dp, iha, ha_dp, iflag_dp, ifail)
        if (ifail /= 0) then
          write(*,'(a,i0)') 'FAIL y12mg n=10 pentadiag (dp): y12mg ifail=', ifail
          nfail = nfail + 1
        else if (rcond_dp <= 0.0d0 .or. rcond_dp > 1.0d0) then
          write(*,'(a,es12.5)') &
              'FAIL y12mg n=10 pentadiag (dp): rcond out of (0,1]: ', rcond_dp
          nfail = nfail + 1
        else
          write(*,'(a,es12.5,a,es12.5)') &
              'PASS y12mg n=10 pentadiag (dp): rcond=', rcond_dp, &
              ' cond=', 1.0d0 / rcond_dp
        end if
      end if
    end if
  end if

  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mg_mh_arrow tests PASSED'

end program test_y12mg_mh_arrow
