module y12m

  implicit none

  private

  ! Factorize and solve
  public :: y12ma
  ! Low level routines: prepare, factorize, and solve
  public :: y12mb
  public :: y12mc
  public :: y12md
  ! Solve with Gaussian elimination and iterative refinement
  public :: y12mf

  public :: y12mg
  public :: y12mh

  ! -------------------------------------------------------------------------
  ! Precision constants for y12mff (double-precision iterative refinement).
  !
  ! y12mff accumulates residuals in higher-than-working precision.  The exact
  ! strategy is selected at compile time using the same logic as in y12mff.f90:
  !
  !   * If the processor supports a real kind with at least 33 significant
  !     decimal digits (i.e. selected_real_kind(33) > 0, quad precision),
  !     that kind is used directly.
  !   * Otherwise a double-double emulation based on ieee_fma is used; in
  !     that case y12mff_accum_kind equals kind(1.0d0).
  !
  ! 80-bit extended-precision (x86) is deliberately ignored: it has only 18
  ! decimal digits, not 33, so selected_real_kind(33) returns -1 on platforms
  ! where only the 80-bit extended type is available.
  ! -------------------------------------------------------------------------

  !> True when y12mff uses native quad precision for residual accumulation;
  !> false when it uses the double-double emulation via ieee_fma.
  logical, parameter, public :: y12mff_has_quad = (selected_real_kind(33) > 0)

  !> Kind parameter of the accumulator used by y12mff.
  !> When y12mff_has_quad is true this equals selected_real_kind(33);
  !> when it is false it equals kind(1.0d0) (double-double emulation).
  integer, parameter, public :: y12mff_accum_kind = &
      merge(selected_real_kind(33), kind(1.0d0), y12mff_has_quad)


  !> Solve sparse systems of linear algebraic equations by
  !> Gaussian elimination
  !>
  interface y12ma
    subroutine y12mae(n, z, a, snr, nn, rnr, nn1, &
        pivot, ha, iha, aflag, iflag, b, ifail)
      implicit none
      integer, intent(in) :: n, z, nn, nn1, iha
      real, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      real, intent(inout) :: pivot(n)
      integer, intent(inout) :: ha(iha,11)
      real, intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      real, intent(inout) :: b(n)
      integer, intent(out) :: ifail
    end subroutine
    subroutine y12maf(n, z, a, snr, nn, rnr, nn1, &
        pivot, ha, iha, aflag, iflag, b, ifail)
      implicit none
      integer, intent(in) :: n, z, nn, nn1, iha
      double precision, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      double precision, intent(inout) :: pivot(n)
      integer, intent(inout) :: ha(iha,11)
      double precision, intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      double precision, intent(inout) :: b(n)
      integer, intent(out) :: ifail
    end subroutine
  end interface

  !> Prepare a system of linear algebraic equations to be factorized
  !>
  interface y12mb
    subroutine y12mbe(n, z, a, snr, nn, rnr, nn1, &
        ha, iha, aflag, iflag, ifail)
      implicit none
      integer n, z, nn, nn1, iha, ifail
      real a(nn), aflag(8)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    end subroutine
    subroutine y12mbf(n, z, a, snr, nn, rnr, nn1, &
        ha, iha, aflag, iflag, ifail)
      implicit none
      integer n, z, nn, nn1, iha, ifail
      double precision a(nn), aflag(8)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    end subroutine
  end interface

  !> Factorize a matrix A into two triangular matrices L and U.
  !>
  interface y12mc
    subroutine y12mce(n, z, a, snr, nn, rnr, nn1, &
        pivot, b, ha, iha, aflag, iflag, ifail)
      implicit none
      integer, intent(in) :: n, z
      integer :: nn, nn1, iha, ifail
      real :: a(nn), pivot(n), b(n), aflag(8)
      integer :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    end subroutine
    subroutine y12mcf(n, z, a, snr, nn, rnr, nn1, &
        pivot, b, ha, iha, aflag, iflag, ifail)
      implicit none
      integer, intent(in) :: n, z
      integer :: nn, nn1, iha, ifail
      double precision :: a(nn), pivot(n), b(n), aflag(8)
      integer :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
    end subroutine
  end interface

  !> Solve sparse systems of linear equations
  !>
  interface y12md
    subroutine y12mde(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      implicit none
      integer, intent(in) :: n, nn, iha
      real, intent(in) :: a(nn)
      real, intent(inout) :: b(n)
      real, intent(in) :: pivot(n)
      integer, intent(in) :: snr(nn)
      integer, intent(in) :: ha(iha,11)
      integer, intent(in) :: iflag(10)
      integer, intent(out) :: ifail
    end subroutine
    subroutine y12mdf(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      implicit none
      integer, intent(in) :: n, nn, iha
      double precision, intent(in) :: a(nn)
      double precision, intent(inout) :: b(n)
      double precision, intent(in) :: pivot(n)
      integer, intent(in) :: snr(nn)
      integer, intent(in) :: ha(iha,11)
      integer, intent(in) :: iflag(10)
      integer, intent(out) :: ifail
    end subroutine
  end interface

  !> Solve large and sparse systems of linear algebraic equations by
  !> the use of Gaussian elimination and sparse matrix technique.
  !>
  !> Iterative refinement is applied in order to improve accuracy.
  !>
  interface y12mf
    subroutine y12mfe(n, a, snr, nn, rnr, nn1, a1,sn, nz, &
        ha, iha, b, b1, x, y, aflag,iflag, ifail)
      implicit none
      integer, intent(in) :: n, nn, nn1, nz, iha
      real, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      real, intent(inout) :: a1(nz)
      integer, intent(inout) :: sn(nz)
      integer, intent(inout) :: ha(iha,13)
      real, intent(inout) :: b(n)
      real, intent(inout) :: b1(n)
      real, intent(inout) :: x(n)
      real, intent(inout) :: y(n)
      real, intent(inout) :: aflag(11)
      integer, intent(inout) :: iflag(12)
      integer, intent(out) :: ifail
    end subroutine
    subroutine y12mff(n, a, snr, nn, rnr, nn1, a1, sn, nz, &
        ha, iha, b, b1, x, y, aflag, iflag, ifail)
      implicit none
      integer, intent(in) :: n, nn, nn1, nz, iha
      double precision, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      double precision, intent(inout) :: a1(nz)
      integer, intent(inout) :: sn(nz)
      integer, intent(inout) :: ha(iha,13)
      double precision, intent(inout) :: b(n)
      double precision, intent(inout) :: b1(n)
      double precision, intent(inout) :: x(n)
      double precision, intent(inout) :: y(n)
      double precision, intent(inout) :: aflag(11)
      integer, intent(inout) :: iflag(12)
      integer, intent(out) :: ifail
    end subroutine
  end interface

  !> Calculate the reciprocal condition number of matrix A
  !>
  interface y12mg
    subroutine y12mge(n,nn,a,snr,w,pivot,anorm,rcond,iha,ha,iflag,ifail)
      implicit none
      integer, intent(in) :: n, nn, iha
      real, intent(in) :: a(nn)
      integer, intent(in) :: snr(nn)
      real, intent(inout) :: w(n)
      real, intent(in) :: pivot(n)
      real, intent(in) :: anorm
      real, intent(out) :: rcond
      integer, intent(in) :: ha(iha,3)
      integer, intent(in) :: iflag(5)
      integer, intent(inout) :: ifail
    end subroutine
    subroutine y12mgf(n,nn,a,snr,w,pivot,anorm,rcond,iha,ha,iflag,ifail)
      implicit none
      integer, intent(in) :: n, nn, iha
      double precision, intent(in) :: a(nn)
      integer, intent(in) :: snr(nn)
      double precision, intent(inout) :: w(n)
      double precision, intent(in) :: pivot(n)
      double precision, intent(in) :: anorm
      double precision, intent(out) :: rcond
      integer, intent(in) :: ha(iha,3)
      integer, intent(in) :: iflag(5)
      integer, intent(inout) :: ifail
    end subroutine
  end interface

  !> Compute the one-norm of a sparse matrix A
  !>
  interface y12mh
    subroutine y12mhe(n,nz,a,snr,work,anorm)
      implicit none
      integer, intent(in) :: n, nz
      real, intent(in) :: a(nz)
      integer, intent(in) :: snr(nz)
      real, intent(inout) :: work(n)
      real, intent(out) :: anorm
    end subroutine
    subroutine y12mhf(n,nz,a,snr,work,anorm)
      implicit none
      integer, intent(in) :: n, nz
      double precision, intent(in) :: a(nz)
      integer, intent(in) :: snr(nz)
      double precision, intent(inout) :: work(n)
      double precision, intent(out) :: anorm
    end subroutine
  end interface

end module