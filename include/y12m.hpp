/* SPDX-License-Identifier: GPL-2.0-only */
/*
 * y12m.hpp -- C++ interface for the y12m sparse solver library.
 *
 * Provides a thin C++ layer on top of the C interface declared in y12m.h:
 *
 *   - A generic (template) interface that dispatches to the float or double
 *     variant of each routine based on the template argument.
 *   - Precondition assertions (active unless NDEBUG is defined) that enforce
 *     the constraints documented in the y12m API.
 *
 * All wrappers live in the `y12m` namespace.
 *
 * Enforced preconditions
 * ----------------------
 * The following constraints are checked with assert() before every call:
 *   n >= 2        system must have at least two equations
 *   z > 0         matrix must have at least one non-zero element
 *   z <= n*n      cannot have more non-zeros than a dense matrix
 *   nn >= 2*z     A/SNR arrays must accommodate fill-in
 *   nn1 >= z      RNR array must hold all non-zero row indices
 *   iha >= n      HA leading dimension must be at least n
 *
 * For the factorize step the additional prerequisite
 *   iflag[0] == -1  (y12mbe/y12mbf must have been called first)
 * is also checked.
 *
 * For the backsolve step:
 *   iflag[0] == -2  (y12mce/y12mcf must have been called first)
 *
 * Example (double precision)
 * --------------------------
 *   #include "y12m.hpp"
 *
 *   // ... fill a, snr, rnr, b with the sparse system data ...
 *   int ifail = y12m::solve<double>(n, z, a, snr, nn, rnr, nn1,
 *                                   pivot, ha, iha, aflag, iflag, b);
 *   if (ifail != 0) { ... handle error ... }
 *   // b now contains the solution vector x
 */

#ifndef Y12M_HPP
#define Y12M_HPP

#include <cassert>
#include "y12m.h"

namespace y12m {

/* =========================================================================
 * solve<T>  --  black-box sparse solve (prepare + factorize + back-solve)
 *
 * Template wrapper around y12mae (T=float) / y12maf (T=double).
 * See y12m.h::y12mae for a full parameter description.
 * ========================================================================= */
template<typename T>
int solve(int n, int z, T *a, int *snr, int nn, int *rnr, int nn1,
          T *pivot, int *ha, int iha, T *aflag, int *iflag, T *b);

template<>
inline int solve<float>(int n, int z, float *a, int *snr, int nn, int *rnr,
                        int nn1, float *pivot, int *ha, int iha,
                        float *aflag, int *iflag, float *b)
{
    assert(n >= 2    && "y12m::solve<float>: n must be >= 2");
    assert(z > 0     && "y12m::solve<float>: z must be > 0");
    assert(z <= n*n  && "y12m::solve<float>: z must be <= n*n");
    assert(nn >= 2*z && "y12m::solve<float>: nn must be >= 2*z");
    assert(nn1 >= z  && "y12m::solve<float>: nn1 must be >= z");
    assert(iha >= n  && "y12m::solve<float>: iha must be >= n");
    return ::y12mae_c(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b);
}

template<>
inline int solve<double>(int n, int z, double *a, int *snr, int nn, int *rnr,
                         int nn1, double *pivot, int *ha, int iha,
                         double *aflag, int *iflag, double *b)
{
    assert(n >= 2    && "y12m::solve<double>: n must be >= 2");
    assert(z > 0     && "y12m::solve<double>: z must be > 0");
    assert(z <= n*n  && "y12m::solve<double>: z must be <= n*n");
    assert(nn >= 2*z && "y12m::solve<double>: nn must be >= 2*z");
    assert(nn1 >= z  && "y12m::solve<double>: nn1 must be >= z");
    assert(iha >= n  && "y12m::solve<double>: iha must be >= n");
    return ::y12maf_c(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b);
}

/* =========================================================================
 * prepare<T>  --  prepare step (y12mbe / y12mbf)
 *
 * Validates and reorders the sparse matrix data before LU factorization.
 * See y12m.h::y12mbe for a full parameter description.
 * ========================================================================= */
template<typename T>
int prepare(int n, int z, T *a, int *snr, int nn, int *rnr, int nn1,
            int *ha, int iha, T *aflag, int *iflag);

template<>
inline int prepare<float>(int n, int z, float *a, int *snr, int nn, int *rnr,
                          int nn1, int *ha, int iha, float *aflag, int *iflag)
{
    assert(n >= 2    && "y12m::prepare<float>: n must be >= 2");
    assert(z > 0     && "y12m::prepare<float>: z must be > 0");
    assert(z <= n*n  && "y12m::prepare<float>: z must be <= n*n");
    assert(nn >= 2*z && "y12m::prepare<float>: nn must be >= 2*z");
    assert(nn1 >= z  && "y12m::prepare<float>: nn1 must be >= z");
    assert(iha >= n  && "y12m::prepare<float>: iha must be >= n");
    return ::y12mbe_c(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag);
}

template<>
inline int prepare<double>(int n, int z, double *a, int *snr, int nn,
                           int *rnr, int nn1, int *ha, int iha,
                           double *aflag, int *iflag)
{
    assert(n >= 2    && "y12m::prepare<double>: n must be >= 2");
    assert(z > 0     && "y12m::prepare<double>: z must be > 0");
    assert(z <= n*n  && "y12m::prepare<double>: z must be <= n*n");
    assert(nn >= 2*z && "y12m::prepare<double>: nn must be >= 2*z");
    assert(nn1 >= z  && "y12m::prepare<double>: nn1 must be >= z");
    assert(iha >= n  && "y12m::prepare<double>: iha must be >= n");
    return ::y12mbf_c(n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag);
}

/* =========================================================================
 * factorize<T>  --  LU factorize step (y12mce / y12mcf)
 *
 * Performs the LU factorization; must be called after prepare().
 * iflag[0] must equal -1 on entry.
 * See y12m.h::y12mce for a full parameter description.
 * ========================================================================= */
template<typename T>
int factorize(int n, int z, T *a, int *snr, int nn, int *rnr, int nn1,
              T *pivot, T *b, int *ha, int iha, T *aflag, int *iflag);

template<>
inline int factorize<float>(int n, int z, float *a, int *snr, int nn,
                            int *rnr, int nn1, float *pivot, float *b,
                            int *ha, int iha, float *aflag, int *iflag)
{
    assert(n >= 2         && "y12m::factorize<float>: n must be >= 2");
    assert(z > 0          && "y12m::factorize<float>: z must be > 0");
    assert(nn >= 2*z      && "y12m::factorize<float>: nn must be >= 2*z");
    assert(nn1 >= z       && "y12m::factorize<float>: nn1 must be >= z");
    assert(iha >= n       && "y12m::factorize<float>: iha must be >= n");
    assert(iflag[0] == -1 && "y12m::factorize<float>: iflag[0] must be -1 "
                             "(call prepare first)");
    return ::y12mce_c(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag);
}

template<>
inline int factorize<double>(int n, int z, double *a, int *snr, int nn,
                             int *rnr, int nn1, double *pivot, double *b,
                             int *ha, int iha, double *aflag, int *iflag)
{
    assert(n >= 2         && "y12m::factorize<double>: n must be >= 2");
    assert(z > 0          && "y12m::factorize<double>: z must be > 0");
    assert(nn >= 2*z      && "y12m::factorize<double>: nn must be >= 2*z");
    assert(nn1 >= z       && "y12m::factorize<double>: nn1 must be >= z");
    assert(iha >= n       && "y12m::factorize<double>: iha must be >= n");
    assert(iflag[0] == -1 && "y12m::factorize<double>: iflag[0] must be -1 "
                             "(call prepare first)");
    return ::y12mcf_c(n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag);
}

/* =========================================================================
 * backsolve<T>  --  back-substitution step (y12mde / y12mdf)
 *
 * Applies the LU factors to solve for x; must be called after factorize().
 * iflag[0] must equal -2 on entry.
 * See y12m.h::y12mde for a full parameter description.
 * ========================================================================= */
template<typename T>
int backsolve(int n, const T *a, int nn, T *b, const T *pivot,
              const int *snr, const int *ha, int iha, const int *iflag);

template<>
inline int backsolve<float>(int n, const float *a, int nn, float *b,
                            const float *pivot, const int *snr,
                            const int *ha, int iha, const int *iflag)
{
    assert(n >= 2         && "y12m::backsolve<float>: n must be >= 2");
    assert(iha >= n       && "y12m::backsolve<float>: iha must be >= n");
    assert(iflag[0] == -2 && "y12m::backsolve<float>: iflag[0] must be -2 "
                             "(call factorize first)");
    return ::y12mde_c(n, a, nn, b, pivot, snr, ha, iha, iflag);
}

template<>
inline int backsolve<double>(int n, const double *a, int nn, double *b,
                             const double *pivot, const int *snr,
                             const int *ha, int iha, const int *iflag)
{
    assert(n >= 2         && "y12m::backsolve<double>: n must be >= 2");
    assert(iha >= n       && "y12m::backsolve<double>: iha must be >= n");
    assert(iflag[0] == -2 && "y12m::backsolve<double>: iflag[0] must be -2 "
                             "(call factorize first)");
    return ::y12mdf_c(n, a, nn, b, pivot, snr, ha, iha, iflag);
}

/* =========================================================================
 * norm1<T>  --  1-norm of a sparse matrix (y12mhe / y12mhf)
 *
 * Computes the maximum absolute column sum of sparse matrix A.
 * See y12m.h::y12mhe for a full parameter description.
 * ========================================================================= */
template<typename T>
void norm1(int n, int nz, const T *a, const int *snr, T *work, T *anorm);

template<>
inline void norm1<float>(int n, int nz, const float *a, const int *snr,
                         float *work, float *anorm)
{
    assert(n  >= 1 && "y12m::norm1<float>: n must be >= 1");
    assert(nz >= 1 && "y12m::norm1<float>: nz must be >= 1");
    ::y12mhe_c(n, nz, a, snr, work, anorm);
}

template<>
inline void norm1<double>(int n, int nz, const double *a, const int *snr,
                          double *work, double *anorm)
{
    assert(n  >= 1 && "y12m::norm1<double>: n must be >= 1");
    assert(nz >= 1 && "y12m::norm1<double>: nz must be >= 1");
    ::y12mhf_c(n, nz, a, snr, work, anorm);
}

/* =========================================================================
 * rcond  --  reciprocal condition number estimate (y12mge, float only)
 *
 * Estimates 1/cond(A) using the LU factorization from factorize().
 * Call norm1<float>() first to obtain anorm, then call this function.
 *
 * Parameters: see y12m.h::y12mge for a full description.
 * ifail_in defaults to 0 (no prior error).
 * ========================================================================= */
inline int rcond(int n, int nn,
                 const float *a, const int *snr,
                 float *w, const float *pivot,
                 float anorm, float *rcond_out,
                 int iha, const int *ha,
                 const int *iflag, int ifail_in = 0)
{
    assert(n   >= 2 && "y12m::rcond: n must be >= 2");
    assert(nn  >= n && "y12m::rcond: nn must be >= n");
    assert(iha >= n && "y12m::rcond: iha must be >= n");
    return ::y12mge_c(n, nn, a, snr, w, pivot, anorm, rcond_out,
                    iha, ha, iflag, ifail_in);
}

} /* namespace y12m */

#endif /* Y12M_HPP */
