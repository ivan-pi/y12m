! SPDX-License-Identifier: GPL-2.0-only
!
! y12m_c.f90 -- Thin bind(c) wrapper module for the y12m sparse solver.
!
! These routines provide C-callable entry points for the legacy y12m
! Fortran subroutines.  Calling conventions:
!
!   * Integer/float/double intent(in) scalar arguments are passed by value.
!   * Array arguments are always passed as C pointers (no value attribute).
!   * IFAIL is returned as the function return value.
!   * Routines without IFAIL (y12mhe, y12mhf) are exposed as void C functions.
!
! Array ha is the integer workspace used by y12m.  Its Fortran declaration is
! ha(iha, K) where K is 11 for y12mb/y12mc/y12md/y12ma/y12mg, 13 for y12mf
! and 3 for the y12mg condition-number subroutine.  From C, ha is passed as a
! flat int array of iha*K elements stored in Fortran column-major order.
!
module y12m_c
  use iso_c_binding, only: c_int, c_float, c_double
  implicit none
  private

contains

  ! ===========================================================================
  ! y12mae -- single-precision black-box solve (prepare + factorize + solve)
  ! ===========================================================================
  function y12mae_c(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b) result(ifail) bind(c, name="y12mae_c")
    integer(c_int), value, intent(in) :: n, z, nn, nn1, iha
    real(c_float),  intent(inout)     :: a(*), pivot(*), aflag(*), b(*)
    integer(c_int), intent(inout)     :: snr(*), rnr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mae
    call y12mae(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)
  end function y12mae_c

  ! ===========================================================================
  ! y12maf -- double-precision black-box solve (prepare + factorize + solve)
  ! ===========================================================================
  function y12maf_c(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
      aflag, iflag, b) result(ifail) bind(c, name="y12maf_c")
    integer(c_int), value, intent(in) :: n, z, nn, nn1, iha
    real(c_double), intent(inout)     :: a(*), pivot(*), aflag(*), b(*)
    integer(c_int), intent(inout)     :: snr(*), rnr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12maf
    call y12maf(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)
  end function y12maf_c

  ! ===========================================================================
  ! y12mbe -- single-precision prepare step
  ! ===========================================================================
  function y12mbe_c(n, z, a, snr, nn, rnr, nn1, ha, iha, &
      aflag, iflag) result(ifail) bind(c, name="y12mbe_c")
    integer(c_int), value, intent(in) :: n, z, nn, nn1, iha
    real(c_float),  intent(inout)     :: a(*), aflag(*)
    integer(c_int), intent(inout)     :: snr(*), rnr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mbe
    call y12mbe(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  end function y12mbe_c

  ! ===========================================================================
  ! y12mbf -- double-precision prepare step
  ! ===========================================================================
  function y12mbf_c(n, z, a, snr, nn, rnr, nn1, ha, iha, &
      aflag, iflag) result(ifail) bind(c, name="y12mbf_c")
    integer(c_int), value, intent(in) :: n, z, nn, nn1, iha
    real(c_double), intent(inout)     :: a(*), aflag(*)
    integer(c_int), intent(inout)     :: snr(*), rnr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mbf
    call y12mbf(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  end function y12mbf_c

  ! ===========================================================================
  ! y12mce -- single-precision factorize step (LU decomposition)
  ! ===========================================================================
  function y12mce_c(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
      aflag, iflag) result(ifail) bind(c, name="y12mce_c")
    integer(c_int), value, intent(in) :: n, z, nn, nn1, iha
    real(c_float),  intent(inout)     :: a(*), pivot(*), b(*), aflag(*)
    integer(c_int), intent(inout)     :: snr(*), rnr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mce
    call y12mce(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
  end function y12mce_c

  ! ===========================================================================
  ! y12mcf -- double-precision factorize step (LU decomposition)
  ! ===========================================================================
  function y12mcf_c(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, &
      aflag, iflag) result(ifail) bind(c, name="y12mcf_c")
    integer(c_int), value, intent(in) :: n, z, nn, nn1, iha
    real(c_double), intent(inout)     :: a(*), pivot(*), b(*), aflag(*)
    integer(c_int), intent(inout)     :: snr(*), rnr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mcf
    call y12mcf(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
  end function y12mcf_c

  ! ===========================================================================
  ! y12mde -- single-precision back-substitution (solve step)
  ! ===========================================================================
  function y12mde_c(n, a, nn, b, pivot, snr, ha, iha, iflag) result(ifail) &
      bind(c, name="y12mde_c")
    integer(c_int), value, intent(in) :: n, nn, iha
    real(c_float),  intent(in)        :: a(*), pivot(*)
    real(c_float),  intent(inout)     :: b(*)
    integer(c_int), intent(in)        :: snr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mde
    call y12mde(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  end function y12mde_c

  ! ===========================================================================
  ! y12mdf -- double-precision back-substitution (solve step)
  ! ===========================================================================
  function y12mdf_c(n, a, nn, b, pivot, snr, ha, iha, iflag) result(ifail) &
      bind(c, name="y12mdf_c")
    integer(c_int), value, intent(in) :: n, nn, iha
    real(c_double), intent(in)        :: a(*), pivot(*)
    real(c_double), intent(inout)     :: b(*)
    integer(c_int), intent(in)        :: snr(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mdf
    call y12mdf(n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  end function y12mdf_c

  ! ===========================================================================
  ! y12mfe -- single-precision solve with iterative refinement
  ! ===========================================================================
  function y12mfe_c(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
      b, b1, x, y, aflag, iflag) result(ifail) bind(c, name="y12mfe_c")
    integer(c_int), value, intent(in) :: n, nn, nn1, nz, iha
    real(c_float),  intent(inout)     :: a(*), a1(*), b(*), b1(*), x(*), y(*)
    real(c_float),  intent(inout)     :: aflag(*)
    integer(c_int), intent(inout)     :: snr(*), rnr(*), sn(*), ha(*), iflag(*)
    integer(c_int)                    :: ifail
    external :: y12mfe
    call y12mfe(n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, &
        b, b1, x, y, aflag, iflag, ifail)
  end function y12mfe_c

  ! ===========================================================================
  ! y12mge -- single-precision reciprocal condition number estimate
  !
  ! ifail_in: pass the IFAIL value returned by y12mce (or 0).  If non-zero
  !           on entry, the subroutine sets rcond = -1 and returns immediately
  !           with the same non-zero code.
  ! ===========================================================================
  function y12mge_c(n, nn, a, snr, w, pivot, anorm, rcond, iha, ha, &
      iflag, ifail_in) result(ifail) bind(c, name="y12mge_c")
    integer(c_int), value, intent(in) :: n, nn, iha, ifail_in
    real(c_float),  value, intent(in) :: anorm
    real(c_float),  intent(in)        :: a(*), pivot(*)
    integer(c_int), intent(in)        :: snr(*), ha(*), iflag(*)
    real(c_float),  intent(inout)     :: w(*)
    real(c_float),  intent(out)       :: rcond
    integer(c_int)                    :: ifail
    external :: y12mge
    ifail = ifail_in
    call y12mge(n, nn, a, snr, w, pivot, anorm, rcond, iha, ha, iflag, ifail)
  end function y12mge_c

  ! ===========================================================================
  ! y12mhe -- single-precision 1-norm of a sparse matrix (no IFAIL)
  ! ===========================================================================
  subroutine y12mhe_c(n, nz, a, snr, work, anorm) bind(c, name="y12mhe_c")
    integer(c_int), value, intent(in) :: n, nz
    real(c_float),  intent(in)        :: a(*)
    integer(c_int), intent(in)        :: snr(*)
    real(c_float),  intent(inout)     :: work(*)
    real(c_float),  intent(out)       :: anorm
    external :: y12mhe
    call y12mhe(n, nz, a, snr, work, anorm)
  end subroutine y12mhe_c

  ! ===========================================================================
  ! y12mhf -- double-precision 1-norm of a sparse matrix (no IFAIL)
  ! ===========================================================================
  subroutine y12mhf_c(n, nz, a, snr, work, anorm) bind(c, name="y12mhf_c")
    integer(c_int), value, intent(in) :: n, nz
    real(c_double), intent(in)        :: a(*)
    integer(c_int), intent(in)        :: snr(*)
    real(c_double), intent(inout)     :: work(*)
    real(c_double), intent(out)       :: anorm
    external :: y12mhf
    call y12mhf(n, nz, a, snr, work, anorm)
  end subroutine y12mhf_c

end module y12m_c
