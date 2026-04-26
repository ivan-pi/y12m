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


  !> Solve a sparse system of linear equations Ax=b by Gaussian elimination.
  !>
  !> This is the all-in-one driver routine. It sets AFLAG(1-4) and IFLAG(2-5)
  !> to fixed internal defaults (overriding any values the caller may have
  !> placed there), then calls y12mb, y12mc, and y12md in sequence.
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

  !> Prepare the non-zero elements of a sparse matrix A in order to
  !> solve the system Ax=b using sparse matrix techniques.
  !>
  !> The routine counts non-zeros per row and column, copies the input
  !> non-zeros to the working portion of A and SNR, records row and column
  !> start positions in HA, reorders the elements into row-ordered and
  !> column-ordered lists, and stores auxiliary ordering information used
  !> by y12mc during Gaussian elimination.
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
  !> L and U using sparse matrix techniques and Gaussian elimination.
  !>
  !> The routine stores information needed to begin the elimination, then
  !> performs the elimination step by step: for each step it selects a
  !> pivot (subject to the stability criterion in AFLAG(1) and the drop
  !> tolerance in AFLAG(2)), applies row and column interchanges, removes
  !> the pivot element, updates the remaining active submatrix, handles
  !> fill-in by compacting the row-ordered and column-ordered lists when
  !> necessary, and updates the row-length ordering list used for pivot
  !> selection. Fill-in information from a previous factorization can be
  !> reused when IFLAG(4) = 2.
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
  !> When IFLAG(5) = 3 (solve with L and U), the routine first performs
  !> forward substitution with the lower triangular factor L, applying
  !> row-interchange permutations if interchanges were used during
  !> factorization. It then performs back substitution with the upper
  !> triangular factor U. Finally, if column interchanges were applied
  !> during elimination, the solution vector is reordered accordingly.
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

  !> Solve a large sparse system of linear equations Ax=b by Gaussian
  !> elimination with iterative refinement to improve accuracy.
  !>
  !> The routine saves the original non-zero values and their column
  !> indices for use in the refinement phase. It calls y12mb and y12mc
  !> to factorize A (unless the factorization is already available), then
  !> calls y12md to compute the first solution. Iterative refinement then
  !> proceeds by computing residuals (in extended precision), solving the
  !> correction equation using the existing factorization, accumulating
  !> corrections, and checking a stopping criterion based on the ratio of
  !> the correction norm to the solution norm. Iterations continue until
  !> convergence or until IFLAG(11) iterations have been performed.
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
  end interface

  !> Compute an estimate of the reciprocal condition number (RCOND) of the
  !> sparse matrix A, as defined by Dongarra, Bunch, Moler and Stewart
  !> (LINPACK User's Guide, SIAM, Philadelphia, 1979).
  !>
  !> Must be called immediately after y12mc while the LU-factorization of
  !> A is still intact. The one-norm of A (ANORM) must be supplied by the
  !> caller; it is most conveniently obtained by calling y12mh before
  !> calling y12mc.
  !>
  !> The algorithm solves U^T * w = e (where the components of e are +1 or
  !> -1 chosen to maximize the norm), then L^T * y = w, and finally
  !> (LU) * z = y. RCOND is estimated as (||y|| / ANORM) / ||z||.
  !> On successful exit RCOND contains an approximation to the reciprocal
  !> condition number. If an error is detected on entry, RCOND is set to -1.
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

  !> Compute the one-norm of a sparse matrix A stored in coordinate form.
  !>
  !> All parameters except ANORM have the same meaning as in the other
  !> routines of the y12m package. On exit, ANORM contains the one-norm
  !> (the maximum absolute column sum) of A. This routine should be called
  !> before y12mc so that ANORM is available for a subsequent call to y12mg.
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