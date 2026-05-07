/* SPDX-License-Identifier: GPL-2.0-only */
/*
 * y12m.h -- C interface for the y12m sparse solver library.
 *
 * Thin static-inline wrappers over the Fortran subroutines that follow the
 * de facto standard Fortran name-mangling convention (lower-case identifier
 * with a trailing underscore, e.g. y12mae becomes y12mae_).
 *
 * Interface adjustments relative to the original Fortran calling sequence:
 *
 *   - INTEGER and REAL/DOUBLE PRECISION intent(in) scalar arguments are
 *     passed by value.
 *   - IFAIL is returned as the function's integer return value instead of
 *     being an argument.  Routines without IFAIL (y12mhe, y12mhf) return void.
 *   - Array arguments that are read-only (Fortran intent(in)) are declared
 *     const.
 *
 * For parameter descriptions and calling conventions see docs/API.md.
 *
 * Two-dimensional Fortran array ha(iha, K) is stored in column-major order.
 * From C, declare it as a flat int array of iha*K elements:
 *
 *   iha*11 ints for y12ma / y12mb / y12mc / y12md / y12mg
 *   iha*13 ints for y12mf
 *
 * All indices in snr and rnr are 1-based (Fortran convention).
 *
 * This header is safe to include from C++ translation units: the declarations
 * are automatically wrapped in an extern "C" block.
 */

#ifndef Y12M_H
#define Y12M_H

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------------
 * Raw Fortran symbols -- de facto gfortran ABI (lower-case + trailing _).
 * All arguments are passed by pointer; Fortran has no const qualifier, so
 * intent(in) arrays that appear const in the public wrappers below are
 * declared non-const here.  The const casts in the wrappers are deliberate
 * and safe: the Fortran routines do not modify intent(in) arguments.
 * These are an implementation detail; call the wrappers below instead.
 * ------------------------------------------------------------------------- */

extern void y12mae_(int *n, int *z,
                    float *a, int *snr, int *nn,
                    int *rnr, int *nn1,
                    float *pivot, int *ha, int *iha,
                    float *aflag, int *iflag,
                    float *b, int *ifail);

extern void y12maf_(int *n, int *z,
                    double *a, int *snr, int *nn,
                    int *rnr, int *nn1,
                    double *pivot, int *ha, int *iha,
                    double *aflag, int *iflag,
                    double *b, int *ifail);

extern void y12mbe_(int *n, int *z,
                    float *a, int *snr, int *nn,
                    int *rnr, int *nn1,
                    int *ha, int *iha,
                    float *aflag, int *iflag, int *ifail);

extern void y12mbf_(int *n, int *z,
                    double *a, int *snr, int *nn,
                    int *rnr, int *nn1,
                    int *ha, int *iha,
                    double *aflag, int *iflag, int *ifail);

extern void y12mce_(int *n, int *z,
                    float *a, int *snr, int *nn,
                    int *rnr, int *nn1,
                    float *pivot, float *b,
                    int *ha, int *iha,
                    float *aflag, int *iflag, int *ifail);

extern void y12mcf_(int *n, int *z,
                    double *a, int *snr, int *nn,
                    int *rnr, int *nn1,
                    double *pivot, double *b,
                    int *ha, int *iha,
                    double *aflag, int *iflag, int *ifail);

extern void y12mde_(int *n,
                    float *a, int *nn,
                    float *b,
                    float *pivot, int *snr,
                    int *ha, int *iha,
                    int *iflag, int *ifail);

extern void y12mdf_(int *n,
                    double *a, int *nn,
                    double *b,
                    double *pivot, int *snr,
                    int *ha, int *iha,
                    int *iflag, int *ifail);

extern void y12mfe_(int *n,
                    float *a, int *snr, int *nn,
                    int *rnr, int *nn1,
                    float *a1, int *sn, int *nz,
                    int *ha, int *iha,
                    float *b, float *b1, float *x, float *y,
                    float *aflag, int *iflag, int *ifail);

extern void y12mge_(int *n, int *nn,
                    float *a, int *snr,
                    float *w, float *pivot,
                    float *anorm, float *rcond,
                    int *iha, int *ha,
                    int *iflag, int *ifail);

extern void y12mhe_(int *n, int *nz,
                    float *a, int *snr,
                    float *work, float *anorm);

extern void y12mhf_(int *n, int *nz,
                    double *a, int *snr,
                    double *work, double *anorm);

/* -------------------------------------------------------------------------
 * Public C interface.
 *
 * Scalar intent(in) arguments are passed by value; IFAIL is the return value.
 * For parameter descriptions see docs/API.md.
 * ------------------------------------------------------------------------- */

/* y12mae -- single-precision black-box solve (prepare + factorize + solve) */
static inline int y12mae(int n, int z,
                          float a[], int snr[], int nn,
                          int rnr[], int nn1,
                          float pivot[], int ha[], int iha,
                          float aflag[], int iflag[],
                          float b[])
{
    int ifail;
    y12mae_(&n, &z, a, snr, &nn, rnr, &nn1, pivot, ha, &iha,
            aflag, iflag, b, &ifail);
    return ifail;
}

/* y12maf -- double-precision black-box solve (prepare + factorize + solve) */
static inline int y12maf(int n, int z,
                          double a[], int snr[], int nn,
                          int rnr[], int nn1,
                          double pivot[], int ha[], int iha,
                          double aflag[], int iflag[],
                          double b[])
{
    int ifail;
    y12maf_(&n, &z, a, snr, &nn, rnr, &nn1, pivot, ha, &iha,
            aflag, iflag, b, &ifail);
    return ifail;
}

/* y12mbe -- single-precision prepare step */
static inline int y12mbe(int n, int z,
                          float a[], int snr[], int nn,
                          int rnr[], int nn1,
                          int ha[], int iha,
                          float aflag[], int iflag[])
{
    int ifail;
    y12mbe_(&n, &z, a, snr, &nn, rnr, &nn1, ha, &iha, aflag, iflag, &ifail);
    return ifail;
}

/* y12mbf -- double-precision prepare step */
static inline int y12mbf(int n, int z,
                          double a[], int snr[], int nn,
                          int rnr[], int nn1,
                          int ha[], int iha,
                          double aflag[], int iflag[])
{
    int ifail;
    y12mbf_(&n, &z, a, snr, &nn, rnr, &nn1, ha, &iha, aflag, iflag, &ifail);
    return ifail;
}

/* y12mce -- single-precision LU factorize step */
static inline int y12mce(int n, int z,
                          float a[], int snr[], int nn,
                          int rnr[], int nn1,
                          float pivot[], float b[],
                          int ha[], int iha,
                          float aflag[], int iflag[])
{
    int ifail;
    y12mce_(&n, &z, a, snr, &nn, rnr, &nn1, pivot, b, ha, &iha,
            aflag, iflag, &ifail);
    return ifail;
}

/* y12mcf -- double-precision LU factorize step */
static inline int y12mcf(int n, int z,
                          double a[], int snr[], int nn,
                          int rnr[], int nn1,
                          double pivot[], double b[],
                          int ha[], int iha,
                          double aflag[], int iflag[])
{
    int ifail;
    y12mcf_(&n, &z, a, snr, &nn, rnr, &nn1, pivot, b, ha, &iha,
            aflag, iflag, &ifail);
    return ifail;
}

/* y12mde -- single-precision back-substitution step */
static inline int y12mde(int n,
                          const float a[], int nn,
                          float b[],
                          const float pivot[], const int snr[],
                          const int ha[], int iha,
                          const int iflag[])
{
    int ifail;
    y12mde_(&n, (float *)a, &nn, b, (float *)pivot, (int *)snr,
            (int *)ha, &iha, (int *)iflag, &ifail);
    return ifail;
}

/* y12mdf -- double-precision back-substitution step */
static inline int y12mdf(int n,
                          const double a[], int nn,
                          double b[],
                          const double pivot[], const int snr[],
                          const int ha[], int iha,
                          const int iflag[])
{
    int ifail;
    y12mdf_(&n, (double *)a, &nn, b, (double *)pivot, (int *)snr,
            (int *)ha, &iha, (int *)iflag, &ifail);
    return ifail;
}

/* y12mfe -- single-precision solve with iterative refinement */
static inline int y12mfe(int n,
                          float a[], int snr[], int nn,
                          int rnr[], int nn1,
                          float a1[], int sn[], int nz,
                          int ha[], int iha,
                          float b[], float b1[], float x[], float y[],
                          float aflag[], int iflag[])
{
    int ifail;
    y12mfe_(&n, a, snr, &nn, rnr, &nn1, a1, sn, &nz, ha, &iha,
            b, b1, x, y, aflag, iflag, &ifail);
    return ifail;
}

/*
 * y12mge -- single-precision reciprocal condition number estimate.
 * ifail_in: IFAIL returned by the preceding y12mce call (pass 0 if none).
 * If ifail_in != 0 on entry, y12mge sets *rcond = -1 and returns immediately.
 */
static inline int y12mge(int n, int nn,
                          const float a[], const int snr[],
                          float w[], const float pivot[],
                          float anorm, float *rcond,
                          int iha, const int ha[],
                          const int iflag[], int ifail_in)
{
    int ifail = ifail_in;
    y12mge_(&n, &nn, (float *)a, (int *)snr, w, (float *)pivot,
            &anorm, rcond, &iha, (int *)ha, (int *)iflag, &ifail);
    return ifail;
}

/* y12mhe -- single-precision 1-norm of a sparse matrix */
static inline void y12mhe(int n, int nz,
                            const float a[], const int snr[],
                            float work[], float *anorm)
{
    y12mhe_(&n, &nz, (float *)a, (int *)snr, work, anorm);
}

/* y12mhf -- double-precision 1-norm of a sparse matrix */
static inline void y12mhf(int n, int nz,
                            const double a[], const int snr[],
                            double work[], double *anorm)
{
    y12mhf_(&n, &nz, (double *)a, (int *)snr, work, anorm);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* Y12M_H */
