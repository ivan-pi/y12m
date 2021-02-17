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


  !> Solve sparse systems of linear algebraic equations by
  !> Gaussian elimination
  !>
  interface y12ma
    subroutine y12mae(n, z, a, snr, nn, rnr, nn1, &
        pivot, ha, iha, aflag, iflag, b, ifail)
      integer, intent(in) :: n
      integer, intent(in) :: z
      real, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(in) :: nn
      integer, intent(inout) :: rnr(nn1)
      integer, intent(in) :: nn1
      real, intent(inout) :: pivot(n)
      integer, intent(inout) :: ha(iha,11)
      integer, intent(in) :: iha
      real, intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      real, intent(inout) :: b(n)
      real, intent(out) :: ifail
    end subroutine
    subroutine y12maf(n, z, a, snr, nn, rnr, nn1, &
        pivot, ha, iha, aflag, iflag, b, ifail)
      integer, intent(in) :: n
      integer, intent(in) :: z
      double precision, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(in) :: nn
      integer, intent(inout) :: rnr(nn1)
      integer, intent(in) :: nn1
      double precision, intent(inout) :: pivot(n)
      integer, intent(inout) :: ha(iha,11)
      integer, intent(in) :: iha
      double precision, intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      double precision, intent(inout) :: b(n)
      double precision, intent(out) :: ifail
    end subroutine
  end interface

  !> Prepare a system of linear algebraic equations to be factorized
  !>
  interface y12mb
    subroutine y12mbe(n, z, a, snr, nn, rnr, nn1, &
        ha, iha, aflag, iflag, ifail)
      real a(nn), aflag(8)
      integer snr(nn), rnr(nn1), ha(n,11), iflag(10)
      integer n, z, iha, ifail
    end subroutine
    subroutine y12mbf(n, z, a, snr, nn, rnr, nn1, &
        ha, iha, aflag, iflag, ifail)
      double precision a(nn), aflag(8)
      integer snr(nn), rnr(nn1), ha(n,11), iflag(10)
      integer n, z, iha, ifail
    end subroutine
  end interface

  !> Factorize a matrix A into two triangular matrices L and U.
  !>
  interface y12mc
    subroutine y12mce(n, z, a, snr, nn, rnr, nn1, &
        pivot, b, ha, iha, aflag, iflag, ifail)
      real :: a(nn), pivot(n), b(n), aflag(8)
      integer :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      integer :: n, z, nn, nn1, iha, ifail
    end subroutine
    subroutine y12mcf(n, z, a, snr, nn, rnr, nn1, &
        pivot, b, ha, iha, aflag, iflag, ifail)
      double precision :: a(nn), pivot(n), b(n), aflag(8)
      integer :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      integer :: n, z, nn, nn1, iha, ifail
    end subroutine
  end interface

  !> Solve sparse systems of linear equations
  !>
  interface y12md
    subroutine y12mde(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      integer, intent(in) :: n
      real, intent(in) :: a(nn)
      integer, intent(in) :: nn
      real, intent(inout) :: b(n)
      real, intent(in) :: pivot(n)
      integer, intent(in) :: snr(nn)
      integer, intent(in) :: ha(iha,11)
      integer, intent(in) :: iha
      integer, intent(in) :: iflag(10)
      integer, intent(out) :: ifail
    end subroutine
    subroutine y12mdf(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      integer, intent(in) :: n
      double precision, intent(in) :: a(nn)
      integer, intent(in) :: nn
      double precision, intent(inout) :: b(n)
      double precision, intent(in) :: pivot(n)
      integer, intent(in) :: snr(nn)
      integer, intent(in) :: ha(iha,11)
      integer, intent(in) :: iha
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
      integer, intent(in) :: n
      integer, intent(in) :: nn
      real, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(inout) :: rnr(nn1)
      integer, intent(in) :: nn1
      real, intent(inout) :: a1(nz)
      integer, intent(inout) :: sn(nz)
      integer, intent(in) :: nz
      integer, intent(inout) :: ha(iha,13)
      integer, intent(in) :: iha
      real, intent(inout) :: b(n)
      real, intent(inout) :: b1(n)
      real, intent(inout) :: x(n)
      real, intent(inout) :: y(n)
      real, intent(inout) :: aflag(11)
      integer, intent(inout) :: iflag(12)
      integer, intent(out) :: ifail
    end subroutine
  end interface

end module