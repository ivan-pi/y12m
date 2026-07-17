! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: Claude Code:claude-fable-5
!
! y12mff.F90 - Double precision version of y12mfe.
!
! Solves a sparse system Ax = b by Gaussian elimination with subsequent
! iterative refinement.  The algorithm is that of the single precision
! routine y12mfe (src/legacy/y12mfe.f); the goto-based control flow of the
! original is expressed here with a named block, a bounded refinement loop,
! and early exits, with identical semantics (same stopping tests, same
! IFLAG/AFLAG outputs).
!
! y12mfe accumulates the residuals r = b - A*x in double precision and
! rounds back to single.  For a double precision working type the
! accumulator must be wider than double.  Two strategies are provided:
!
!   * quad precision, selected_real_kind(33), when the compiler offers it;
!   * a double-double compensated accumulation (TwoProd + TwoSum, error
!     terms carried in a low-order word) which yields a residual as
!     accurate as if it were computed in twice double precision
!     (Ogita/Rump/Oishi "Dot2").
!
! The strategy is fixed at compile time:
!
!   Y12M_ACCUM_QUAD  defined: compile the quad path only.
!   Y12M_ACCUM_DD    defined: compile the double-double path only.
!   neither          defined: compile both paths and pick quad when
!                             selected_real_kind(33) > 0 (the unused path
!                             is dead code behind a constant condition).
!
! The exact product splitting (TwoProd) in the double-double path uses
! Dekker's algorithm by default: portable Fortran, no library calls.
! Define Y12M_USE_IEEE_FMA to use the Fortran 2018 IEEE_FMA intrinsic
! instead (one fma per element instead of Dekker's 17 flops; worthwhile
! only when the compiler inlines it to a hardware instruction rather than
! dispatching to a runtime call).
!
! The CMake build selects the macros from a configure-time capability
! check and the Y12M_FORCE_DOUBLE_DOUBLE / Y12M_USE_IEEE_FMA options; the
! self-detecting mode serves build systems that do not preprocess
! capability checks (e.g. fpm).
!
! WARNING: the double-double path relies on IEEE-compliant evaluation.
!
!   * Value-unsafe optimisation (-ffast-math, -Ofast, ifx's default
!     -fp-model=fast, nvfortran without -Kieee) eliminates the
!     compensation terms, which are algebraically zero, and the
!     accumulation silently degrades to plain double precision.
!   * FMA contraction (gfortran/flang default -ffp-contract=fast on
!     hardware with fma instructions) breaks Dekker's splitting: fusing
!     t - aj with t = splitter*aj skips exactly the rounding the split
!     depends on.  The IEEE_FMA variant is immune to contraction.
!
! The CMake build pins the required flags on this file (e.g. for gfortran
! -fno-unsafe-math-optimizations, plus -ffp-contract=off for the Dekker
! variant); when compiling the double-double path by other means, apply
! the equivalent flags yourself.
!
#if defined(Y12M_ACCUM_QUAD) && defined(Y12M_ACCUM_DD)
#error "Y12M_ACCUM_QUAD and Y12M_ACCUM_DD are mutually exclusive"
#endif
subroutine y12mff(n, a, snr, nn, rnr, nn1, a1, sn, nz, &
      ha, iha, b, b1, x, y, aflag, iflag, ifail)
#if !defined(Y12M_ACCUM_QUAD) && defined(Y12M_USE_IEEE_FMA)
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

   real(dp) :: d, dd, dres, gt1, gt2, xx, xm, er
   integer :: nres, state, kit, it
   integer :: i, j, l1, l2, pass

   ifail = 0
   nres = 0
   dd = 0.0_dp
   dres = 0.0_dp
   xm = 0.0_dp
   state = iflag(5)
   kit = 1
   it = iflag(11)

   solve: block

      if (state == 1) ifail = 10
      if (it < 2) ifail = 23
      if (ifail /= 0) exit solve

      ! Save the right-hand side; b is overwritten by the solves.
      b1(1:n) = b(1:n)

      ! state == 3 reuses an existing factorization: skip straight to the
      ! back-substitution for the new right-hand side.
      if (state /= 3) then

         call y12mbf(n, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
         if (ifail /= 0) exit solve

         ! Keep a row-ordered copy of the original matrix (and its row
         ! start/end positions) for the residual computation: y12mcf
         ! overwrites a/snr with the factors.
         a1(1:nz) = a(1:nz)
         sn(1:nz) = snr(1:nz)
         do i = 1, n
            ha(i,12) = ha(i,1)
            ha(i,13) = ha(i,3)
         end do

         ! A negative aflag(2) requests a relative drop tolerance: scale it
         ! by the smallest row maximum (but never above aflag(6) = max|A|).
         if (aflag(2) < 0.0_dp) then
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
         end if

         ! Factorize and find the first solution.
         call y12mcf(n, nz, a, snr, nn, rnr, nn1, y, b, ha, iha, &
             aflag, iflag, ifail)
         if (ifail /= 0) exit solve

      end if

      call y12mdf(n, a, nn, b, y, snr, ha, iha, iflag, ifail)
      if (ifail /= 0) exit solve

      ! Prepare the data in order to begin the iterations.
      x(1:n) = b(1:n)
      dd = maxval(abs(b(1:n)))
      xm = dd
      if (dd == 0.0_dp) exit solve

      ! Refinement.  Pass k computes the residual r = b1 - A*x(k-1); unless
      ! a stopping test fires it then performs the kit-th solve for the
      ! correction d(k) and updates x.  When nres is set (correction
      ! negligible, or iteration budget reached) one final residual-only
      ! pass runs so that the reported dres and the residual left in b
      ! belong to the final x; "it" therefore bounds the number of passes.
      ! kit counts the solves (initial solution included) and cannot serve
      ! as the loop index: the final pass does not advance it, and it is
      ! reported in iflag(12).
      refine: do pass = 1, it

         d = dd
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
         if (dres == 0.0_dp) exit refine

         ! Final residual-only pass done, or residual blow-up.
         if (nres == 1 .or. dres > 1.0e4_dp*xm) then
            dd = abs(dd)
            exit refine
         end if

         ! Solve for the correction (in b) and measure it.
         kit = kit + 1
         iflag(5) = 3
         call y12mdf(n, a, nn, b, y, snr, ha, iha, iflag, ifail)
         if (ifail /= 0) exit refine
         dd = maxval(abs(b(1:n)))
         if (dd == 0.0_dp) exit refine

         ! Divergence: the corrections grow after the second solve.
         if (dd > d .and. kit > 2) exit refine

         ! Calculate an improved solution.
         x(1:n) = x(1:n) + b(1:n)
         xm = maxval(abs(x(1:n)))

         ! Stopping criteria: negligible correction, or no iterations left.
         ! Either way one final residual-only pass follows.
         if (10.0_dp + dd/xm == 10.0_dp .or. kit >= it) nres = 1

      end do refine

   end block solve

   ! Restore the caller's flags and report the diagnostics.
   iflag(5) = state
   iflag(12) = kit
   aflag(9) = dd
   aflag(10) = dres
   aflag(11) = xm

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
   !> "Dot2").  Each product is split exactly into p + pe (TwoProd); the
   !> subtraction from the high word er_hi is made exact with a TwoSum,
   !> and both error terms are collected in the low word er_lo.  The
   !> result is as accurate as a residual computed in twice double
   !> precision.
   !>
   !> Correctness relies on IEEE-compliant evaluation; see the file
   !> header note on -ffast-math and -ffp-contract.
   pure function residual_dd(b1i, l1, l2) result(er)
      real(dp), intent(in) :: b1i
      integer, intent(in) :: l1, l2
      real(dp) :: er
#ifndef Y12M_USE_IEEE_FMA
      real(dp), parameter :: splitter = 134217729.0_dp  ! 2**27 + 1
      real(dp) :: t, a_hi, a_lo, x_hi, x_lo
#endif
      real(dp) :: er_hi, er_lo, aj, xj, p, pe, s, z, e
      integer :: j, l3
      er_hi = b1i
      er_lo = 0.0_dp
      do j = l1, l2
         l3 = sn(j)
         aj = a1(j)
         xj = x(l3)
         ! TwoProd: p + pe = aj*xj exactly
         p = aj*xj
#ifdef Y12M_USE_IEEE_FMA
         pe = ieee_fma(aj, xj, -p)
#else
         ! Dekker splitting: both factors are cut into 26-bit halves whose
         ! four partial products are exact, recovering the rounding error
         ! of p.  Overflows for |aj| or |xj| >= 2**996; residuals of such
         ! magnitude are far outside the solver's working range.
         t = splitter*aj
         a_hi = t - (t - aj)
         a_lo = aj - a_hi
         t = splitter*xj
         x_hi = t - (t - xj)
         x_lo = xj - x_hi
         pe = ((a_hi*x_hi - p) + a_hi*x_lo + a_lo*x_hi) + a_lo*x_lo
#endif
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
