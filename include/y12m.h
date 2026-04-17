/* SPDX-License-Identifier: GPL-2.0-only */
/*
 * y12m.h -- C interface for the y12m sparse solver library.
 *
 * Each routine follows Zlatev's original interface with the following
 * adjustments made for idiomatic C use:
 *
 *   - Integer and float/double intent(in) scalar arguments are passed by
 *     value instead of by pointer.
 *   - IFAIL is returned as the function's integer return value.
 *     Routines without IFAIL (y12mhe_c, y12mhf_c) return void.
 *   - Array arguments that are read-only (intent(in) in Fortran) are
 *     decorated with const.
 *
 * The function names carry a "_c" suffix to distinguish them from the
 * underlying Fortran legacy routines (y12mae, y12maf, ...) which are
 * accessible through the Fortran module interface.
 *
 * Two-dimensional Fortran arrays (ha) are stored in column-major order.
 * The flat index for Fortran element ha(i,j) is (j-1)*iha + (i-1).
 * From C, ha is simply an int array; the library manages all internal
 * indexing.  Allocate iha*11 ints for y12ma/y12mb/y12mc/y12md,
 * iha*13 ints for y12mf, and iha*3 ints are sufficient for y12mge.
 *
 * Index conventions
 * -----------------
 * All indices stored in SNR and RNR are 1-based (Fortran convention).
 * Array positions are also 1-based: the first non-zero value is a[0]
 * paired with snr[0] and rnr[0].
 *
 * Error codes (IFAIL)
 * -------------------
 * A return value of 0 indicates success.  Non-zero values signal errors;
 * see the y12m documentation for the full list of error codes.
 *
 * Usage with extern "C"
 * ---------------------
 * This header is safe to include from C++ translation units; the
 * declarations are wrapped in an extern "C" block automatically.
 */

#ifndef Y12M_H
#define Y12M_H

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * y12mae_c -- single-precision black-box sparse solve
 *
 * Solves the sparse n x n linear system A x = b by Gaussian elimination
 * in a single call (prepare + LU factorize + back-substitute).
 *
 * Parameters
 * ----------
 * n        Number of equations; must be >= 2.
 * z        Number of non-zero elements of A stored in a[0..z-1]; must
 *          satisfy 0 < z <= n*n.
 * a        [inout]  Array of length nn.  On entry a[0..z-1] contains the
 *          non-zero values of A.  The array is overwritten on exit.
 * snr      [inout]  Array of length nn.  On entry snr[0..z-1] contains the
 *          1-based column indices of the non-zero values.  Modified on exit.
 * nn       Length of arrays a and snr; must satisfy nn >= 2*z.
 * rnr      [inout]  Array of length nn1.  On entry rnr[0..z-1] contains
 *          the 1-based row indices of the non-zero values.  Modified on exit.
 * nn1      Length of array rnr; must satisfy nn1 >= z.
 * pivot    [inout]  Workspace array of length n.
 * ha       [inout]  Integer workspace of iha*11 elements stored in
 *          Fortran column-major order (see header notes).  iha must be >= n.
 * iha      First (leading) dimension of ha; must be >= n.
 * aflag    [inout]  Real array of length 8 carrying solver tolerances.
 *          y12mae_c overwrites aflag[0..3] with its own defaults; user-
 *          supplied values for these four entries are ignored.
 * iflag    [inout]  Integer array of length 10 carrying solver control flags.
 *          y12mae_c overwrites iflag[1..4] with its own defaults; user-
 *          supplied values for these four entries are ignored.
 * b        [inout]  Array of length n.  On entry the right-hand side vector
 *          b; on successful exit the solution vector x.
 *
 * Returns 0 on success; a non-zero error code on failure.
 * ========================================================================= */
int y12mae_c(int n, int z,
             float *a, int *snr, int nn,
             int *rnr, int nn1,
             float *pivot, int *ha, int iha,
             float *aflag, int *iflag,
             float *b);

/* =========================================================================
 * y12maf_c -- double-precision black-box sparse solve
 *
 * Identical to y12mae_c but uses double precision for all real arrays.
 * ========================================================================= */
int y12maf_c(int n, int z,
             double *a, int *snr, int nn,
             int *rnr, int nn1,
             double *pivot, int *ha, int iha,
             double *aflag, int *iflag,
             double *b);

/* =========================================================================
 * y12mbe_c -- single-precision prepare step
 *
 * Validates and reorders the sparse matrix data stored in a/snr/rnr in
 * preparation for LU factorization.  Must be called before y12mce_c.
 * On successful return iflag[0] == -1.
 *
 * Parameters: same as y12mae_c except no pivot and no b.  The aflag and
 * iflag arrays must be initialised by the caller before this call; see
 * the y12m documentation for recommended values.
 *
 * Returns 0 on success; a non-zero error code on failure.
 * ========================================================================= */
int y12mbe_c(int n, int z,
             float *a, int *snr, int nn,
             int *rnr, int nn1,
             int *ha, int iha,
             float *aflag, int *iflag);

/* =========================================================================
 * y12mbf_c -- double-precision prepare step
 *
 * Identical to y12mbe_c but uses double precision for all real arrays.
 * ========================================================================= */
int y12mbf_c(int n, int z,
             double *a, int *snr, int nn,
             int *rnr, int nn1,
             int *ha, int iha,
             double *aflag, int *iflag);

/* =========================================================================
 * y12mce_c -- single-precision LU factorize step
 *
 * Performs the LU factorization of the matrix prepared by y12mbe_c.
 * iflag[0] must equal -1 on entry (set by y12mbe_c).
 * On successful return iflag[0] == -2.
 *
 * Parameters: same as y12mae_c.
 *
 * Returns 0 on success; a non-zero error code on failure.
 * ========================================================================= */
int y12mce_c(int n, int z,
             float *a, int *snr, int nn,
             int *rnr, int nn1,
             float *pivot, float *b,
             int *ha, int iha,
             float *aflag, int *iflag);

/* =========================================================================
 * y12mcf_c -- double-precision LU factorize step
 *
 * Identical to y12mce_c but uses double precision for all real arrays.
 * ========================================================================= */
int y12mcf_c(int n, int z,
             double *a, int *snr, int nn,
             int *rnr, int nn1,
             double *pivot, double *b,
             int *ha, int iha,
             double *aflag, int *iflag);

/* =========================================================================
 * y12mde_c -- single-precision back-substitution step
 *
 * Solves the system using the LU factors produced by y12mce_c.
 * iflag[0] must equal -2 on entry (set by y12mce_c).
 *
 * Parameters
 * ----------
 * n        Number of equations.
 * a        [in]    LU factors as stored by y12mce_c in array a (length nn).
 * nn       Length of arrays a and snr.
 * b        [inout] On entry the right-hand side vector; on exit the solution.
 * pivot    [in]    Pivot values as stored by y12mce_c in array pivot (length n).
 * snr      [in]    Column-index array as left by y12mce_c (length nn).
 * ha       [in]    Integer workspace as left by y12mce_c (iha*11 elements).
 * iha      First dimension of ha.
 * iflag    [in]    Control flags as left by y12mce_c; iflag[0] must be -2.
 *
 * Returns 0 on success; a non-zero error code on failure.
 * ========================================================================= */
int y12mde_c(int n,
             const float *a, int nn,
             float *b,
             const float *pivot,
             const int *snr,
             const int *ha, int iha,
             const int *iflag);

/* =========================================================================
 * y12mdf_c -- double-precision back-substitution step
 *
 * Identical to y12mde_c but uses double precision for all real arrays.
 * ========================================================================= */
int y12mdf_c(int n,
             const double *a, int nn,
             double *b,
             const double *pivot,
             const int *snr,
             const int *ha, int iha,
             const int *iflag);

/* =========================================================================
 * y12mfe_c -- single-precision solve with iterative refinement
 *
 * Solves the sparse system Ax = b using Gaussian elimination with
 * iterative refinement to improve accuracy.  Requires additional storage
 * arrays a1, sn, b1, x, y and an enlarged ha workspace (iha*13 elements).
 *
 * Parameters
 * ----------
 * n        Number of equations; must be >= 2.
 * a        [inout]  Non-zero values of A (length nn).
 * snr      [inout]  1-based column indices (length nn).
 * nn       Length of arrays a and snr; must satisfy nn >= 2*z.
 * rnr      [inout]  1-based row indices (length nn1).
 * nn1      Length of array rnr; must satisfy nn1 >= z.
 * a1       [inout]  Auxiliary copy of non-zero values (length nz).
 * sn       [inout]  Auxiliary copy of column indices (length nz).
 * nz       Number of non-zero elements stored in a1/sn (original count).
 * ha       [inout]  Integer workspace of iha*13 elements (iha >= n).
 * iha      First (leading) dimension of ha; must be >= n.
 * b        [inout]  On entry the right-hand side; on exit the solution.
 * b1       [inout]  Workspace array of length n.
 * x        [inout]  Workspace / iterative-refinement output array (length n).
 * y        [inout]  Workspace array of length n.
 * aflag    [inout]  Real array of length 11 with solver flags/tolerances.
 * iflag    [inout]  Integer array of length 12 with solver control flags.
 *
 * Returns 0 on success; a non-zero error code on failure.
 * ========================================================================= */
int y12mfe_c(int n,
             float *a, int *snr, int nn,
             int *rnr, int nn1,
             float *a1, int *sn, int nz,
             int *ha, int iha,
             float *b, float *b1, float *x, float *y,
             float *aflag, int *iflag);

/* =========================================================================
 * y12mge_c -- single-precision reciprocal condition number estimate
 *
 * Estimates 1/cond(A) using the LU factorization produced by y12mce_c.
 * Call y12mhe_c first to obtain anorm, then call y12mge_c.
 *
 * Parameters
 * ----------
 * n        Number of equations.
 * nn       Length of arrays a and snr.
 * a        [in]    LU factors as stored by y12mce_c (length nn).
 * snr      [in]    Column indices as stored by y12mce_c (length nn).
 * w        [inout] Workspace array of length n.
 * pivot    [in]    Pivot values as stored by y12mce_c (length n).
 * anorm    [in]    1-norm of the original matrix A, as returned by y12mhe_c.
 * rcond    [out]   On exit an approximation to 1/cond(A).  Set to -1.0
 *          if ifail_in is non-zero on entry.
 * iha      First dimension of ha; must be >= n.
 * ha       [in]    Integer workspace as left by y12mce_c (iha*11 elements;
 *          only the first 3 columns are used by y12mge_c).
 * iflag    [in]    Control flags as left by y12mce_c (first 5 entries used).
 * ifail_in IFAIL value from the preceding call (pass 0 if no prior error).
 *          If non-zero, the subroutine sets *rcond = -1.0 and returns
 *          immediately with the same non-zero code.
 *
 * Returns 0 on success; a non-zero error code on failure.
 * ========================================================================= */
int y12mge_c(int n, int nn,
             const float *a, const int *snr,
             float *w,
             const float *pivot,
             float anorm, float *rcond,
             int iha, const int *ha,
             const int *iflag,
             int ifail_in);

/* =========================================================================
 * y12mhe_c -- single-precision 1-norm of a sparse matrix
 *
 * Computes the 1-norm (maximum absolute column sum) of sparse matrix A.
 * Typically called before y12mce_c to obtain anorm for a subsequent
 * y12mge_c condition-number estimate.
 *
 * Parameters
 * ----------
 * n        Number of columns of the matrix (must be >= 1).
 * nz       Number of non-zero elements (must be >= 1).
 * a        [in]    Non-zero values (length nz).
 * snr      [in]    1-based column indices of the non-zero values (length nz).
 * work     [inout] Workspace array of length n (overwritten on exit).
 * anorm    [out]   On exit the 1-norm of the matrix.
 * ========================================================================= */
void y12mhe_c(int n, int nz,
              const float *a, const int *snr,
              float *work, float *anorm);

/* =========================================================================
 * y12mhf_c -- double-precision 1-norm of a sparse matrix
 *
 * Identical to y12mhe_c but uses double precision for all real arrays.
 * ========================================================================= */
void y12mhf_c(int n, int nz,
              const double *a, const int *snr,
              double *work, double *anorm);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* Y12M_H */
