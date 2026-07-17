! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: Claude Code:claude-fable-5
!
! test_y12mf_dp.f90 - Tests for y12mf (double-precision iterative refinement).
!
! This program tests y12mf (resolving to y12mff) which combines Gaussian
! elimination with iterative refinement.  Residuals are accumulated in
! higher-than-double precision: native quad when available, otherwise a
! double-double compensated accumulation (see src/y12mff.F90).  The
! accumulator actually compiled in is reported via the y12mff_has_quad and
! y12mff_accum_kind module constants.
!
! Systems solved (all with known solution x = [1,...,1]):
!   tridiagonal (diag=3, off-diag=-1): n = 3, 5, 7, 10
!   arrow (full first row/col + diagonal): n = 6, 10
!
! The factored form is stored with IFLAG(5)=2 (keep L factors) as required
! by y12mf.  The refined solution is returned in array x (not b).
!
! Pass criteria: consistent accumulator constants, and IFAIL=0 with
! max|x_i-1| < 1e-13 for all systems.
!
program test_y12mf_dp
   use y12m, only: y12mf, y12mff_has_quad, y12mff_accum_kind
   use y12m_test_helpers, only: build_tridiag, build_arrow, check_solution
   implicit none

   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: NMAX = 12
   integer, parameter :: NNP = 400
   integer, parameter :: NN1P = 400
   real(dp), parameter :: TOL = 1.0e-13_dp

   real(dp) :: a(NNP), a1(NNP), b(NMAX), b1(NMAX), x(NMAX), y(NMAX)
   real(dp) :: aflag(11)
   integer :: snr(NNP), sn(NNP), rnr(NN1P), ha(NMAX,13), iflag(12)
   integer :: n, z, nn, nn1, iha, ifail, nfail

   nfail = 0

   ! Report and sanity-check the accumulator selected at compile time.
   if (y12mff_has_quad) then
      write(*,'(a,i0,a)') 'y12mff accumulator: quad precision (kind=', &
          y12mff_accum_kind, ')'
      if (y12mff_accum_kind /= selected_real_kind(33)) then
         write(*,'(a)') 'FAIL y12mff_accum_kind inconsistent with quad selection'
         nfail = nfail + 1
      end if
   else
      write(*,'(a,i0,a)') 'y12mff accumulator: double-double (kind=', &
          y12mff_accum_kind, ')'
      if (y12mff_accum_kind /= dp) then
         write(*,'(a)') 'FAIL y12mff_accum_kind inconsistent with double-double selection'
         nfail = nfail + 1
      end if
   end if

   ! Tridiagonal systems: n = 3, 5, 7, 10
   call run_tridiag(3)
   call run_tridiag(5)
   call run_tridiag(7)
   call run_tridiag(10)

   ! Arrow systems: n = 6, 10
   call run_arrow(6)
   call run_arrow(10)

   if (nfail /= 0) then
      write(*,'(i0,a)') nfail, ' test(s) FAILED'
      stop 1
   end if
   write(*,'(a)') 'All test_y12mf_dp tests PASSED'

contains

   subroutine run_tridiag(n_in)
      integer, intent(in) :: n_in
      character(len=30) :: label
      n = n_in
      write(label,'(a,i0,a)') 'y12mf_dp n=', n, ' tridiag'
      call build_tridiag(n, NNP, NN1P, z, a, snr, rnr, b)
      call solve_and_check(trim(label))
   end subroutine run_tridiag

   subroutine run_arrow(n_in)
      integer, intent(in) :: n_in
      character(len=30) :: label
      n = n_in
      write(label,'(a,i0,a)') 'y12mf_dp n=', n, ' arrow'
      call build_arrow(n, NNP, NN1P, z, a, snr, rnr, b)
      call solve_and_check(trim(label))
   end subroutine run_arrow

   !> Call y12mf on the COO system already loaded in a/snr/rnr/b and verify
   !> the refined solution against the known solution x = [1,...,1].
   subroutine solve_and_check(label)
      character(len=*), intent(in) :: label
      integer :: i

      nn = NNP
      nn1 = NN1P
      iha = NMAX

      ! Settings per the reference example (test2.f).
      ! IFLAG(5)=2 keeps L factors as required by y12mf.
      aflag(1) = 16.0_dp     ! stability factor
      aflag(2) = 1.0e-12_dp  ! drop tolerance
      aflag(3) = 1.0e+16_dp  ! max growth factor before abort
      aflag(4) = 1.0e-12_dp  ! min pivot threshold

      iflag(1) = 0
      iflag(2) = 3
      iflag(3) = 1
      iflag(4) = 1   ! preserve structure info for iterative refinement
      iflag(5) = 2   ! keep L factors (required by y12mf)
      iflag(11) = 10 ! maximum number of refinement iterations

      x(1:n) = 0.0_dp
      y(1:n) = 0.0_dp
      b1(1:n) = 0.0_dp
      do i = 1, z
         a1(i) = 0.0_dp
         sn(i) = 0
      end do

      call y12mf(n, a, snr, nn, rnr, nn1, a1, sn, z, ha, iha, &
          b, b1, x, y, aflag, iflag, ifail)

      call check_solution(label, n, x, ifail, TOL, nfail)
   end subroutine solve_and_check

end program test_y12mf_dp
