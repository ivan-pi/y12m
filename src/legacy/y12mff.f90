! SPDX-License-Identifier: GPL-2.0-only
!
! y12mff.f90 - Double precision version of y12mfe.
!
! Solves large and sparse systems of linear equations by Gaussian
! elimination with subsequent iterative refinement.  This is the double
! precision counterpart to y12mfe: all working arrays are double precision.
!
! The residuals r_i = b1(i) - sum_j(a1(j)*x(sn(j))) are accumulated in
! higher-than-working precision so that the iterative correction is
! meaningful.  The implementation checks at compile time whether a real
! kind with at least 33 decimal digits (quad precision) is available.  If
! so, quad arithmetic is used for the accumulator.  Otherwise a compact
! double-double accumulation based on TwoSum and TwoProduct (via
! ieee_arithmetic's ieee_fma) gives approximately 106 binary digits.
!
subroutine y12mff(n, a, snr, nn, rnr, nn1, a1, sn, nz, &
    ha, iha, b, b1, x, y, aflag, iflag, ifail)
  use, intrinsic :: ieee_arithmetic, only: ieee_fma
  implicit none

  integer, intent(in) :: n, nn, nn1, nz, iha
  real(kind(1.0d0)), intent(inout) :: a(nn)
  integer, intent(inout) :: snr(nn)
  integer, intent(inout) :: rnr(nn1)
  real(kind(1.0d0)), intent(inout) :: a1(nz)
  integer, intent(inout) :: sn(nz)
  integer, intent(inout) :: ha(iha,13)
  real(kind(1.0d0)), intent(inout) :: b(n)
  real(kind(1.0d0)), intent(inout) :: b1(n)
  real(kind(1.0d0)), intent(inout) :: x(n)
  real(kind(1.0d0)), intent(inout) :: y(n)
  real(kind(1.0d0)), intent(inout) :: aflag(11)
  integer, intent(inout) :: iflag(12)
  integer, intent(out) :: ifail

  ! ---------------------------------------------------------------------------
  ! Precision selection
  ! ---------------------------------------------------------------------------
  ! Working precision (double precision).
  integer, parameter :: dp = kind(1.0d0)
  ! Try for true quad precision (at least 33 decimal digits).
  ! We deliberately skip 80-bit extended precision (which covers 18 digits
  ! but not 33) so that the fallback double-double path is used on x86-64
  ! when only the 80-bit type is available.
  integer, parameter :: qp = selected_real_kind(33)
  logical, parameter :: has_qp = (qp > 0)
  ! hp is the higher precision kind.  When quad is not available hp equals
  ! dp, but in that branch the double-double accumulator is used instead.
  integer, parameter :: hp = merge(qp, dp, has_qp)

  ! ---------------------------------------------------------------------------
  ! Local scalars
  ! ---------------------------------------------------------------------------
  real(dp) :: d, dd, dres, gt1, gt2, xx, xm, dcor
  integer :: nres, state, kit, it
  integer :: i, j, l1, l2, l3

  ! Higher-precision residual accumulator.
  ! er_hp    - quad accumulator (used when has_qp is true)
  ! er_hi/lo - double-double accumulator (used when has_qp is false)
  ! er       - the residual rounded back to working precision
  real(hp) :: er_hp
  real(dp) :: er_hi, er_lo, er

  ifail = 0
  nres  = 0
  dres  = 0.0_dp
  state = iflag(5)
  kit   = 1
  it    = iflag(11)
  if (state == 1) ifail = 10
  if (it < 2)     ifail = 23
  if (ifail /= 0) go to 160

  do i = 1, n
    b1(i) = b(i)
  end do

  if (state == 3) go to 70

  call y12mbf(n, nz, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  if (ifail /= 0) go to 160

  do i = 1, nz
    sn(i)  = snr(i)
    a1(i)  = a(i)
  end do

  do i = 1, n
    ha(i,12) = ha(i,1)
    ha(i,13) = ha(i,3)
  end do

  if (aflag(2) >= 0.0_dp) go to 60

  gt1 = aflag(6)
  do i = 1, n
    l1   = ha(i,1)
    l2   = ha(i,3)
    gt2  = 0.0_dp
    do j = l1, l2
      d = abs(a(j))
      if (gt2 < d) gt2 = d
    end do
    if (gt2 < gt1) gt1 = gt2
  end do
  aflag(2) = -gt1 * aflag(2)

60 call y12mcf(n, nz, a, snr, nn, rnr, nn1, y, b, ha, iha, aflag, iflag, ifail)
  if (ifail /= 0) go to 160

70 call y12mdf(n, a, nn, b, y, snr, ha, iha, iflag, ifail)
  if (ifail /= 0) go to 160

  ! Prepare data for iteration.
  dd = 0.0_dp
  do i = 1, n
    x(i) = b(i)
    xx   = abs(b(i))
    if (dd < xx) dd = xx
  end do
  xm = dd
  if (dd == 0.0_dp) go to 160

  ! Begin iterative refinement.
90 d    = dd
  dres = 0.0_dp

  do i = 1, n
    l1 = ha(i,12)
    l2 = ha(i,13)

    ! Compute residual r_i = b1(i) - sum_j a1(j)*x(sn(j))
    ! in higher-than-working precision.
    if (has_qp) then
      ! Quad precision path.
      er_hp = real(b1(i), hp)
      do j = l1, l2
        l3    = sn(j)
        er_hp = er_hp - real(a1(j), hp) * real(x(l3), hp)
      end do
      er = real(er_hp, dp)
    else
      ! Double-double path.
      ! The accumulator (er_hi, er_lo) represents the sum er_hi + er_lo with
      ! approximately twice the precision of double precision.
      er_hi = b1(i)
      er_lo = 0.0_dp
      do j = l1, l2
        l3 = sn(j)
        call dd_sub_product(er_hi, er_lo, a1(j), x(l3))
      end do
      er = er_hi + er_lo
    end if

    ! Store residual rounded to working precision.
    b(i) = er
    xx   = abs(er)
    if (dres < xx) dres = xx
  end do

  if (dres == 0.0_dp) go to 160
  if (nres == 1)       go to 150
  if (dres > 1.0e4_dp * xm) go to 150

  kit      = kit + 1
  iflag(5) = 3
  call y12mdf(n, a, nn, b, y, snr, ha, iha, iflag, ifail)
  if (ifail /= 0) go to 160

  ! Compute the uniform norm of the current correction.
  dd = 0.0_dp
  do i = 1, n
    xx = abs(b(i))
    if (dd < xx) dd = xx
  end do
  if (dd == 0.0_dp) go to 160

  ! Check convergence.
  if (dd > d .and. kit > 2) go to 160

  ! Update solution.
  dcor = dd
  xm   = 0.0_dp
  do i = 1, n
    x(i) = x(i) + b(i)
    xx   = abs(x(i))
    if (xx > xm) xm = xx
  end do

  ! Check stopping criterion.
  if (10.0_dp + dd/xm == 10.0_dp) go to 140
  if (kit < it) go to 90

140 nres = 1
  go to 90

150 dd = abs(dd)

160 iflag(5)  = state
  iflag(12) = kit
  aflag(9)  = dd
  aflag(10) = dres
  aflag(11) = xm
  return

contains

  ! Update the double-double accumulator (s_hi, s_lo) by subtracting a*b.
  !
  ! Uses TwoProduct (via ieee_fma) to compute the exact product a*b as a
  ! double-double pair (p, p_err), then uses TwoSum to subtract p from s_hi
  ! while collecting the rounding error into s_lo.
  !
  subroutine dd_sub_product(s_hi, s_lo, a, b)
    real(dp), intent(inout) :: s_hi, s_lo
    real(dp), intent(in) :: a, b
    real(dp) :: p, p_err, t, z

    ! TwoProduct: p + p_err = a * b  (exact, one rounding only)
    p     = a * b
    p_err = ieee_fma(a, b, -p)

    ! TwoSum: (s_hi + (-p)) = t + z  (exact)
    t = s_hi - p
    z = t - s_hi
    ! Rounding error of the subtraction s_hi - p:
    !   e = (s_hi - (t - z)) + ((-p) - z)
    s_lo = s_lo - p_err + ((s_hi - (t - z)) + (-p - z))
    s_hi = t
  end subroutine dd_sub_product

end subroutine y12mff
