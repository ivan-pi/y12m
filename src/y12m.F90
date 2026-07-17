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

  ! Estimate reciprocal condition number
  public :: y12mg

  ! Compute the matrix one-norm
  public :: y12mh

  ! Error string
  public :: y12m_get_error_msg

  ! -------------------------------------------------------------------------
  ! Residual accumulation strategy of y12mff (double precision iterative
  ! refinement).  y12mff accumulates the residuals b - A*x in higher than
  ! working precision: native quad when available, otherwise a double-double
  ! (compensated) accumulation.  The selection is fixed at compile time and
  ! must match src/y12mff.F90, so the same Y12M_ACCUM_QUAD / Y12M_ACCUM_DD
  ! macros drive both files (the CMake build defines exactly one of them).
  ! -------------------------------------------------------------------------

  !> True when y12mff accumulates residuals in native quad precision;
  !> false when it uses the double-double emulation.  Note "uses", not
  !> "has": with the double-double path forced at build time this is false
  !> even on a compiler that supports quad precision.
#if defined(Y12M_ACCUM_DD)
  logical, parameter, public :: y12mff_uses_quad = .false.
#elif defined(Y12M_ACCUM_QUAD)
  logical, parameter, public :: y12mff_uses_quad = .true.
#else
  logical, parameter, public :: y12mff_uses_quad = selected_real_kind(33) > 0
#endif

  !> Kind parameter of the y12mff residual accumulator.  Equal to
  !> selected_real_kind(33) when y12mff_uses_quad is true; otherwise equal
  !> to kind(1.0d0), the working kind on which the double-double arithmetic
  !> is built.
  integer, parameter, public :: y12mff_accumulator_kind = &
      merge(selected_real_kind(33), kind(1.0d0), y12mff_uses_quad)


  !> Solve a sparse system of linear equations Ax=b by Gaussian elimination.
  !>
  !> All-in-one driver: calls y12mb, y12mc, and y12md in sequence.
  !> AFLAG(1-4) and IFLAG(2-5) are set to fixed internal defaults.
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

  !> Prepare the non-zero elements of the sparse matrix A for factorization.
  !>
  !> Reorders the input non-zeros into row-ordered and column-ordered lists
  !> and builds auxiliary data structures in HA for use by y12mc.
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

  !> Factorize the sparse matrix A into lower and upper triangular factors
  !> L and U using sparse Gaussian elimination.
  !>
  !> Pivot selection uses stability criterion AFLAG(1) and drop tolerance
  !> AFLAG(2). Set IFLAG(4) = 2 to reuse fill-in from a previous call.
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

  !> Solve the system Ax=b using the LU-factorization computed by y12mc.
  !>
  !> Performs forward substitution with L and back substitution with U.
  !> Applies row and column permutations if interchanges were used.
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

  !> Solve a large sparse system Ax=b by Gaussian elimination with
  !> iterative refinement to improve accuracy.
  !>
  !> Calls y12mb, y12mc, and y12md to compute the initial solution, then
  !> refines iteratively until convergence or IFLAG(11) iterations.
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

  !> Estimate the reciprocal condition number RCOND of the sparse matrix A.
  !>
  !> Must be called immediately after y12mc. Supply ANORM (the one-norm of A)
  !> obtained from y12mh. Returns RCOND = -1 if an input error is detected.
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

  !> Compute the one-norm (maximum absolute column sum) of the sparse matrix A.
  !>
  !> Call before y12mc to obtain ANORM for a subsequent call to y12mg.
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

contains

  subroutine y12m_get_error_msg(ifail, message, length)
      implicit none

      ! Arguments
      integer, intent(in) :: ifail
      character(len=*), intent(out) :: message
      integer, intent(out), optional :: length

      ! Local Table of Error Messages
      character(len=256), dimension(25), parameter :: ERROR_TABLE = [ &
          character(len=256) :: &
          "The coefficient matrix A is not factorized. Call Y12MC before Y12MD, or ensure IFLAG(1) is 0 or greater.", &
          "The coefficient matrix A is not ordered. Ensure Y12MB is called successfully before calling Y12MC.", &
          "A numerically singular pivot element has been selected. Matrix is likely singular; try increasing AFLAG(4).", &
          "Large growth factor detected. Elements grew too quickly; try using a smaller stability factor in AFLAG(1).", &
          "Length NN of arrays A and SNR is insufficient. Increase the size of NN and potentially NN1.", &
          "Length NN1 of array RNR is insufficient. Increase the size of NN1 and potentially NN.", &
          "A row with zeroes only in its active part found. If drop-tolerance is small, try decreasing AFLAG(2).", &
          "A column with zeroes only in its active part found. The matrix is likely singular. Check AFLAG(8) and AFLAG(5).", &
          "A pivotal element is missing. Decrease drop-tolerance AFLAG(2) or check if the strategy matches the matrix.", &
          "Invalid request to remove non-zeros of matrix L. Set IFLAG(5) to 2 instead of 1 before calling Y12MF.", &
          "Duplicate elements found at the same position. Check input data for overlapping row and column indices.", &
          "The number of equations is too small. N must be 2 or larger.", &
          "The number of non-zero elements is zero or negative. Check the NZ or Z parameter in the input data.", &
          "Number of non-zeros is smaller than the number of equations. The matrix is structurally singular.", &
          "Array HA dimension IHA is too small. Ensure IHA is greater than or equal to N.", &
          "Invalid IFLAG(4) value. Set this parameter to 0, 1, or 2.", &
          "Row with zeroes only found before elimination. The matrix is structurally singular and cannot be solved.", &
          "Column with zeroes only found before elimination. The matrix is structurally singular and cannot be solved.", &
          "IFLAG(2) is too small. Set IFLAG(2) to a positive integer; 3 is the recommended value.", &
          "IFLAG(3) is out of range. Valid options are 0, 1, or 2.", &
          "IFLAG(5) is out of range. Valid options are 1, 2, or 3.", &
          "Subroutines Y12MB and Y12MC require IFLAG(5) to be 1 or 2. Check setting before calling.", &
          "Allowed iterations IFLAG(11) is too low. Set IFLAG(11) to 2 or higher.", &
          "Column index error. Found an element with a column index outside the range 1 to N.", &
          "Row index error. Found an element with a row index outside the range 1 to N." &
      ]

      select case(ifail)
      case(0)
        associate(errmsg => "Success: no error detected.")
          if (present(length)) length = len(errmsg)
          message = errmsg
        end associate
      case(1:25)
        if (present(length)) length = len_trim(ERROR_TABLE(ifail))
        message = ERROR_TABLE(ifail)
      case default
        associate(errmsg => "Unknown error code.")
          if (present(length)) length = len(errmsg)
          message = errmsg
        end associate
      end select

  end subroutine

end module y12m