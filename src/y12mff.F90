! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: Claude Code:claude-fable-5
!
! y12mff.F90 - Double precision version of y12mfe.
!
! Solves a sparse system Ax = b by Gaussian elimination with subsequent
! iterative refinement.  The algorithm and the control flow follow the
! single precision routine y12mfe (src/legacy/y12mfe.f) exactly; only the
! working precision and the residual accumulation differ.
!
! y12mfe accumulates the residuals r = b - A*x in double precision and
! rounds back to single.  For a double precision working type the
! accumulator must be wider than double.  Two strategies are provided:
!
!   * quad precision, selected_real_kind(33), when the compiler offers it;
!   * a double-double compensated accumulation (TwoProd via IEEE_FMA plus
!     TwoSum, error terms carried in a low-order word) which yields a
!     residual as accurate as if it were computed in twice double
!     precision (Ogita/Rump/Oishi "Dot2").
!
! The strategy is fixed at compile time:
!
!   Y12M_ACCUM_QUAD  defined: compile the quad path only.
!   Y12M_ACCUM_DD    defined: compile the double-double path only.
!   neither          defined: compile both paths and pick quad when
!                             selected_real_kind(33) > 0 (the unused path
!                             is dead code behind a constant condition).
!                             This mode requires IEEE_FMA (Fortran 2018).
!
! The CMake build always defines exactly one of the two macros, driven by
! a configure-time capability check and the Y12M_FORCE_DOUBLE_DOUBLE
! option; the self-detecting mode serves build systems that do not
! preprocess capability checks (e.g. fpm).
!
! WARNING: the double-double path is destroyed by value-unsafe floating
! point optimisation (-ffast-math, -Ofast, ifx's default -fp-model=fast,
! nvfortran without -Kieee).  The compensation terms are algebraically
! zero, so a compiler licensed to reassociate eliminates them and the
! accumulation silently degrades to plain double precision.  The CMake
! build pins value-safe flags for this file; keep that in mind when
! compiling it by other means.
!
#if defined(Y12M_ACCUM_QUAD) && defined(Y12M_ACCUM_DD)
#error "Y12M_ACCUM_QUAD and Y12M_ACCUM_DD are mutually exclusive"
#endif
subroutine y12mff(n, a, snr, nn, rnr, nn1, a1, sn, nz, &
      ha, iha, b, b1, x, y, aflag, iflag, ifail)
#ifndef Y12M_ACCUM_QUAD
   use, intrinsic :: ieee_arithmetic, only: ieee_fma
#endif
   implicit none

   integer, parameter :: dp = kind(1.0d0)
#ifdef Y12M_ACCUM_DD
   logical, parameter :: use_quad = .false.
#else
   integer, parameter :: qp = selected_real_kind(33)
#ifdef Y12M_ACCUM_QUAD
   logical, parameter :: use_quad = .true.
#else
   logical, parameter :: use_quad = qp > 0
#endif
   ! When quad is unavailable hp falls back to dp; the quad path is then
   ! dead code behind the constant use_quad and is never executed.
   integer, parameter :: hp = merge(qp, dp, use_quad)
#endif

   integer, intent(in) :: n, nn, nn1, nz, iha
   real(dp), intent(inout) :: a(nn)
   integer, intent(inout) :: snr(nn)
   integer, intent(inout) :: rnr(nn1)
   real(dp), intent(inout) :: a1(nz)
   integer, intent(inout) :: sn(nz)
   integer, intent(inout) :: ha(iha,13)
   real(dp), intent(inout) :: b(n)
   real(dp), intent(inout) :: b1(n)
   real(dp), intent(inout) :: x(n)
   real(dp), intent(inout) :: y(n)
   real(dp), intent(inout) :: aflag(11)
   integer, intent(inout) :: iflag(12)
   integer, intent(out) :: ifail

   external y12mbf, y12mcf, y12mdf

   real(dp) :: d, dd, dres, gt1, gt2, xx, xm, dcor, er
   integer :: nres, state, kit, it
   integer :: i, j, l1, l2

   ! Store the non-zero elements, their column numbers, information about
   ! row starts, information about row ends, and the right-hand side.
   ifail = 0
   nres = 0
   dres = 0.0_dp
   state = iflag(5)
   kit = 1
   it = iflag(11)
   if (state == 1) ifail = 10
   if (it < 2) ifail = 23
   if (ifail /= 0) go to 160

   do i = 1, n
      b1(i) = b(i)
   end do

   if (state == 3) go to 70

   call y12mbf(n, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
   if (ifail /= 0) go to 160

   do i = 1, nz
      sn(i) = snr(i)
      a1(i) = a(i)
   end do

   do i = 1, n
      ha(i,12) = ha(i,1)
      ha(i,13) = ha(i,3)
   end do

   if (aflag(2) >= 0.0_dp) go to 60

   gt1 = aflag(6)
   do i = 1, n
      l1 = ha(i,1)
      l2 = ha(i,3)
      gt2 = 0.0_dp
      do j = l1, l2
         d = abs(a(j))
         if (gt2 < d) gt2 = d
      end do
      if (gt2 < gt1) gt1 = gt2
   end do
   aflag(2) = -gt1*aflag(2)

   ! Find the first solution.
60 call y12mcf(n, nz, a, snr, nn, rnr, nn1, y, b, ha, iha, aflag, iflag, ifail)
   if (ifail /= 0) go to 160

70 call y12mdf(n, a, nn, b, y, snr, ha, iha, iflag, ifail)
   if (ifail /= 0) go to 160

   ! Prepare the data in order to begin the iterations.
   dd = 0.0_dp
   do i = 1, n
      x(i) = b(i)
      xx = abs(b(i))
      if (dd < xx) dd = xx
   end do
   xm = dd
   if (dd == 0.0_dp) go to 160

   ! Begin to iterate.
90 d = dd
   dres = 0.0_dp
   do i = 1, n
      l1 = ha(i,12)
      l2 = ha(i,13)
#if defined(Y12M_ACCUM_DD)
      er = residual_dd(b1(i), l1, l2)
#elif defined(Y12M_ACCUM_QUAD)
      er = residual_quad(b1(i), l1, l2)
#else
      if (use_quad) then
         er = residual_quad(b1(i), l1, l2)
      else
         er = residual_dd(b1(i), l1, l2)
      end if
#endif
      ! Store residuals rounded to working (double) precision.
      b(i) = er
      xx = abs(er)
      if (dres < xx) dres = xx
   end do
   if (dres == 0.0_dp) go to 160
   if (nres == 1) go to 150
   if (dres > 1.0e4_dp*xm) go to 150

   kit = kit + 1
   iflag(5) = 3
   call y12mdf(n, a, nn, b, y, snr, ha, iha, iflag, ifail)
   if (ifail /= 0) go to 160

   ! Compute the uniform norm of the current correction vector.
   dd = 0.0_dp
   do i = 1, n
      xx = abs(b(i))
      if (dd < xx) dd = xx
   end do
   if (dd == 0.0_dp) go to 160

   ! Check the convergence criterion.
   if (dd > d .and. kit > 2) go to 160

   ! Calculate an improved solution.
   dcor = dd
   xm = 0.0_dp
   do i = 1, n
      x(i) = x(i) + b(i)
      xx = abs(x(i))
      if (xx > xm) xm = xx
   end do

   ! Check the stopping criteria.
   if (10.0_dp + dd/xm == 10.0_dp) go to 140
   if (kit < it) go to 90

   ! End of the iterations.
140 nres = 1
   go to 90

150 dd = abs(dd)

160 iflag(5) = state
   iflag(12) = kit
   aflag(9) = dd
   aflag(10) = dres
   aflag(11) = xm
   return

contains

#ifndef Y12M_ACCUM_DD
   !> Row residual b1i - sum_j a1(j)*x(sn(j)), j = l1..l2, accumulated in
   !> quad precision and rounded once to double.
   pure function residual_quad(b1i, l1, l2) result(er)
      real(dp), intent(in) :: b1i
      integer, intent(in) :: l1, l2
      real(dp) :: er
      real(hp) :: er_hp
      integer :: j, l3
      er_hp = real(b1i, hp)
      do j = l1, l2
         l3 = sn(j)
         er_hp = er_hp - real(a1(j), hp)*real(x(l3), hp)
      end do
      er = real(er_hp, dp)
   end function residual_quad
#endif

#ifndef Y12M_ACCUM_QUAD
   !> Row residual b1i - sum_j a1(j)*x(sn(j)), j = l1..l2, accumulated in
   !> double-double arithmetic (compensated dot product, Ogita/Rump/Oishi
   !> "Dot2").  Each product is split exactly into p + pe with IEEE_FMA
   !> (TwoProd); the subtraction from the high word er_hi is made exact
   !> with a TwoSum, and both error terms are collected in the low word
   !> er_lo.  The result is as accurate as a residual computed in twice
   !> double precision, using only *, +, - and one fma per element.
   !>
   !> Correctness relies on IEEE-compliant evaluation: the compensation
   !> terms are algebraically zero and must not be simplified away (see
   !> the file header note on -ffast-math).
   pure function residual_dd(b1i, l1, l2) result(er)
      real(dp), intent(in) :: b1i
      integer, intent(in) :: l1, l2
      real(dp) :: er
      real(dp) :: er_hi, er_lo, p, pe, s, z, e
      integer :: j, l3
      er_hi = b1i
      er_lo = 0.0_dp
      do j = l1, l2
         l3 = sn(j)
         ! TwoProd: p + pe = a1(j)*x(l3) exactly
         p = a1(j)*x(l3)
         pe = ieee_fma(a1(j), x(l3), -p)
         ! TwoSum: s + e = er_hi - p exactly
         s = er_hi - p
         z = s - er_hi
         e = (er_hi - (s - z)) - (p + z)
         er_hi = s
         er_lo = er_lo + (e - pe)
      end do
      er = er_hi + er_lo
   end function residual_dd
#endif

end subroutine y12mff
