/* SPDX-License-Identifier: GPL-2.0-only */
/*
 * y12m.h -- C interface for the y12m sparse solver library.
 *
 * Thin static-inline wrappers over the Fortran subroutines.
 * By default FNAME(name) maps to the common gfortran ABI naming
 * convention (lower-case identifier with trailing underscore).
 *
 * Interface adjustments relative to the original Fortran calling sequence:
 *
 *   - INTEGER and REAL/DOUBLE PRECISION intent(in) scalar arguments are
 *     passed by value.
 *   - IFAIL is returned as the function's integer return value instead of
 *     being an argument, except for Y12MG where IFAIL remains an in/out
 *     argument and the function result is RCOND.
 *   - Y12MH returns ANORM as the function result.
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

#ifndef FNAME
#define FNAME(name) name ## _
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void FNAME(y12mae)(int *n, int *z,
                          float *a, int *snr, int *nn,
                          int *rnr, int *nn1,
                          float *pivot, int *ha, int *iha,
                          float *aflag, int *iflag,
                          float *b, int *ifail);

extern void FNAME(y12maf)(int *n, int *z,
                          double *a, int *snr, int *nn,
                          int *rnr, int *nn1,
                          double *pivot, int *ha, int *iha,
                          double *aflag, int *iflag,
                          double *b, int *ifail);

extern void FNAME(y12mbe)(int *n, int *z,
                          float *a, int *snr, int *nn,
                          int *rnr, int *nn1,
                          int *ha, int *iha,
                          float *aflag, int *iflag, int *ifail);

extern void FNAME(y12mbf)(int *n, int *z,
                          double *a, int *snr, int *nn,
                          int *rnr, int *nn1,
                          int *ha, int *iha,
                          double *aflag, int *iflag, int *ifail);

extern void FNAME(y12mce)(int *n, int *z,
                          float *a, int *snr, int *nn,
                          int *rnr, int *nn1,
                          float *pivot, float *b,
                          int *ha, int *iha,
                          float *aflag, int *iflag, int *ifail);

extern void FNAME(y12mcf)(int *n, int *z,
                          double *a, int *snr, int *nn,
                          int *rnr, int *nn1,
                          double *pivot, double *b,
                          int *ha, int *iha,
                          double *aflag, int *iflag, int *ifail);

extern void FNAME(y12mde)(int *n,
                          float *a, int *nn,
                          float *b,
                          float *pivot, int *snr,
                          int *ha, int *iha,
                          int *iflag, int *ifail);

extern void FNAME(y12mdf)(int *n,
                          double *a, int *nn,
                          double *b,
                          double *pivot, int *snr,
                          int *ha, int *iha,
                          int *iflag, int *ifail);

extern void FNAME(y12mfe)(int *n,
                          float *a, int *snr, int *nn,
                          int *rnr, int *nn1,
                          float *a1, int *sn, int *nz,
                          int *ha, int *iha,
                          float *b, float *b1, float *x, float *y,
                          float *aflag, int *iflag, int *ifail);

extern void FNAME(y12mge)(int *n, int *nn,
                          float *a, int *snr,
                          float *w, float *pivot,
                          float *anorm, float *rcond,
                          int *iha, int *ha,
                          int *iflag, int *ifail);

extern void FNAME(y12mgf)(int *n, int *nn,
                          double *a, int *snr,
                          double *w, double *pivot,
                          double *anorm, double *rcond,
                          int *iha, int *ha,
                          int *iflag, int *ifail);

extern void FNAME(y12mhe)(int *n, int *nz,
                          float *a, int *snr,
                          float *work, float *anorm);

extern void FNAME(y12mhf)(int *n, int *nz,
                          double *a, int *snr,
                          double *work, double *anorm);

static inline int y12mae(int n, int z,
                         float a[], int snr[], int nn,
                         int rnr[], int nn1,
                         float pivot[], int ha[], int iha,
                         float aflag[], int iflag[],
                         float b[])
{
   int ifail;
   FNAME(y12mae)(&n, &z, a, snr, &nn, rnr, &nn1, pivot, ha, &iha,
                 aflag, iflag, b, &ifail);
   return ifail;
}

static inline int y12maf(int n, int z,
                         double a[], int snr[], int nn,
                         int rnr[], int nn1,
                         double pivot[], int ha[], int iha,
                         double aflag[], int iflag[],
                         double b[])
{
   int ifail;
   FNAME(y12maf)(&n, &z, a, snr, &nn, rnr, &nn1, pivot, ha, &iha,
                 aflag, iflag, b, &ifail);
   return ifail;
}

static inline int y12mbe(int n, int z,
                         float a[], int snr[], int nn,
                         int rnr[], int nn1,
                         int ha[], int iha,
                         float aflag[], int iflag[])
{
   int ifail;
   FNAME(y12mbe)(&n, &z, a, snr, &nn, rnr, &nn1, ha, &iha, aflag, iflag, &ifail);
   return ifail;
}

static inline int y12mbf(int n, int z,
                         double a[], int snr[], int nn,
                         int rnr[], int nn1,
                         int ha[], int iha,
                         double aflag[], int iflag[])
{
   int ifail;
   FNAME(y12mbf)(&n, &z, a, snr, &nn, rnr, &nn1, ha, &iha, aflag, iflag, &ifail);
   return ifail;
}

static inline int y12mce(int n, int z,
                         float a[], int snr[], int nn,
                         int rnr[], int nn1,
                         float pivot[], float b[],
                         int ha[], int iha,
                         float aflag[], int iflag[])
{
   int ifail;
   FNAME(y12mce)(&n, &z, a, snr, &nn, rnr, &nn1, pivot, b, ha, &iha,
                 aflag, iflag, &ifail);
   return ifail;
}

static inline int y12mcf(int n, int z,
                         double a[], int snr[], int nn,
                         int rnr[], int nn1,
                         double pivot[], double b[],
                         int ha[], int iha,
                         double aflag[], int iflag[])
{
   int ifail;
   FNAME(y12mcf)(&n, &z, a, snr, &nn, rnr, &nn1, pivot, b, ha, &iha,
                 aflag, iflag, &ifail);
   return ifail;
}

static inline int y12mde(int n,
                         const float a[], int nn,
                         float b[],
                         const float pivot[], const int snr[],
                         const int ha[], int iha,
                         const int iflag[])
{
   int ifail;
   FNAME(y12mde)(&n, (float *)a, &nn, b, (float *)pivot, (int *)snr,
                 (int *)ha, &iha, (int *)iflag, &ifail);
   return ifail;
}

static inline int y12mdf(int n,
                         const double a[], int nn,
                         double b[],
                         const double pivot[], const int snr[],
                         const int ha[], int iha,
                         const int iflag[])
{
   int ifail;
   FNAME(y12mdf)(&n, (double *)a, &nn, b, (double *)pivot, (int *)snr,
                 (int *)ha, &iha, (int *)iflag, &ifail);
   return ifail;
}

static inline int y12mfe(int n,
                         float a[], int snr[], int nn,
                         int rnr[], int nn1,
                         float a1[], int sn[], int nz,
                         int ha[], int iha,
                         float b[], float b1[], float x[], float y[],
                         float aflag[], int iflag[])
{
   int ifail;
   FNAME(y12mfe)(&n, a, snr, &nn, rnr, &nn1, a1, sn, &nz, ha, &iha,
                 b, b1, x, y, aflag, iflag, &ifail);
   return ifail;
}

static inline float y12mge(int n, int nn,
                           const float a[], const int snr[],
                           float w[], const float pivot[],
                           float anorm,
                           const int ha[], int iha,
                           const int iflag[], int *ifail)
{
   float rcond;
   FNAME(y12mge)(&n, &nn, (float *)a, (int *)snr, w, (float *)pivot,
                 &anorm, &rcond, &iha, (int *)ha, (int *)iflag, ifail);
   return rcond;
}

static inline double y12mgf(int n, int nn,
                            const double a[], const int snr[],
                            double w[], const double pivot[],
                            double anorm,
                            const int ha[], int iha,
                            const int iflag[], int *ifail)
{
   double rcond;
   FNAME(y12mgf)(&n, &nn, (double *)a, (int *)snr, w, (double *)pivot,
                 &anorm, &rcond, &iha, (int *)ha, (int *)iflag, ifail);
   return rcond;
}

static inline float y12mhe(int n, int nz,
                           const float a[], const int snr[],
                           float work[])
{
   float anorm;
   FNAME(y12mhe)(&n, &nz, (float *)a, (int *)snr, work, &anorm);
   return anorm;
}

static inline double y12mhf(int n, int nz,
                            const double a[], const int snr[],
                            double work[])
{
   double anorm;
   FNAME(y12mhf)(&n, &nz, (double *)a, (int *)snr, work, &anorm);
   return anorm;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* Y12M_H */
