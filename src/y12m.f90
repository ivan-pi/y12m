module y12m

  implicit none

  private

  public :: y12ma
  public :: y12mb
  public :: y12mc
  public :: y12md

  public :: y12mfe

  interface y12ma
    subroutine y12mae(n, z, a, snr, nn, rnr, nn1, &
      pivot, ha, iha, aflag, iflag, b, ifail)
      real a(nn), pivot(n), b(n), aflag(8)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      integer n, z, nn, nn1, iha, ifail
    end subroutine
    subroutine y12maf(n, z, a, snr, nn, rnr, nn1, &
        pivot, ha, iha, aflag, iflag, b, ifail)
      double precision a(nn), pivot(n), b(n), aflag(8)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      integer n, z, nn, nn1, iha, ifail
    end subroutine
  end interface


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

  interface y12md
    subroutine y12mde(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      real :: a(nn), pivot(n), b(n)
      integer :: snr(nn), ha(iha,11), iflag(10)
      integer :: n, nn, iha, ifail
    end subroutine
    subroutine y12mdf(n, a, nn, b, pivot, snr, &
        ha, iha, iflag, ifail)
      double precision :: a(nn), pivot(n), b(n)
      integer :: snr(nn), ha(iha,11), iflag(10)
      integer :: n, nn, iha, ifail
    end subroutine
  end interface

  interface
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