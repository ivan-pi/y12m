/* SPDX-License-Identifier: GPL-2.0-only */
/*
 * test_y12m_c.c -- C interface smoke tests for the y12m sparse solver.
 *
 * Tests the C binding declared in y12m.h by solving small tridiagonal
 * systems with the known solution x = [1, ..., 1]:
 *
 *   1) Single-precision black-box solve via y12mae_c.
 *   2) Double-precision black-box solve via y12maf_c.
 *   3) Double-precision three-step API: y12mbf_c -> y12mcf_c -> y12mdf_c.
 *
 * Pass criteria: IFAIL == 0 and max|x_i - 1| < 1e-5 (float) or 1e-10 (double).
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "y12m.h"

/* ----- sizing constants -------------------------------------------------- */
#define N    6          /* system dimension                */
#define NNP  50         /* length of a/snr arrays          */
#define NN1P 50         /* length of rnr array             */
#define IHA  12         /* first dimension of ha           */

/* ----- helpers ----------------------------------------------------------- */

static int nfail = 0;

static void check_sp(const char *label, int n, const float *b,
                     int ifail, float tol)
{
    float err = 0.0f;
    int i;
    if (ifail != 0) {
        printf("FAIL %s  ifail=%d\n", label, ifail);
        ++nfail;
        return;
    }
    for (i = 0; i < n; i++) {
        float d = (float)fabs((double)(b[i] - 1.0f));
        if (d > err) err = d;
    }
    if (err > tol) {
        printf("FAIL %s  max_err=%g\n", label, (double)err);
        ++nfail;
    } else {
        printf("PASS %s  max_err=%g\n", label, (double)err);
    }
}

static void check_dp(const char *label, int n, const double *b,
                     int ifail, double tol)
{
    double err = 0.0;
    int i;
    if (ifail != 0) {
        printf("FAIL %s  ifail=%d\n", label, ifail);
        ++nfail;
        return;
    }
    for (i = 0; i < n; i++) {
        double d = fabs(b[i] - 1.0);
        if (d > err) err = d;
    }
    if (err > tol) {
        printf("FAIL %s  max_err=%g\n", label, err);
        ++nfail;
    } else {
        printf("PASS %s  max_err=%g\n", label, err);
    }
}

/*
 * Build an n×n tridiagonal system with diag=3, off-diag=-1.
 * The exact solution is x = [1,...,1], so b = A*x = row-sums.
 * Triplets are stored in COO format with 1-based indices.
 */
static void build_tridiag_sp(int n, int *z_out,
                              float *a, int *snr, int *rnr, float *b)
{
    int z = 0, i;
    for (i = 1; i <= n; i++) {
        a[z] =  3.0f; snr[z] = i;   rnr[z] = i; ++z;
    }
    for (i = 2; i <= n; i++) {
        a[z] = -1.0f; snr[z] = i-1; rnr[z] = i; ++z;
    }
    for (i = 1; i <= n-1; i++) {
        a[z] = -1.0f; snr[z] = i+1; rnr[z] = i; ++z;
    }
    for (i = 0; i < n; i++) b[i] = 0.0f;
    for (i = 0; i < z; i++) b[rnr[i]-1] += a[i];
    *z_out = z;
}

static void build_tridiag_dp(int n, int *z_out,
                              double *a, int *snr, int *rnr, double *b)
{
    int z = 0, i;
    for (i = 1; i <= n; i++) {
        a[z] =  3.0; snr[z] = i;   rnr[z] = i; ++z;
    }
    for (i = 2; i <= n; i++) {
        a[z] = -1.0; snr[z] = i-1; rnr[z] = i; ++z;
    }
    for (i = 1; i <= n-1; i++) {
        a[z] = -1.0; snr[z] = i+1; rnr[z] = i; ++z;
    }
    for (i = 0; i < n; i++) b[i] = 0.0;
    for (i = 0; i < z; i++) b[rnr[i]-1] += a[i];
    *z_out = z;
}

/* ----- tests ------------------------------------------------------------- */

/* 1. Single-precision black-box solve */
static void test_y12mae(void)
{
    float  a[NNP], pivot[IHA], b[N], aflag[8];
    int    snr[NNP], rnr[NN1P], ha[IHA * 11], iflag[10];
    int    z, ifail;

    build_tridiag_sp(N, &z, a, snr, rnr, b);
    ifail = y12mae_c(N, z, a, snr, NNP, rnr, NN1P,
                   pivot, ha, IHA, aflag, iflag, b);
    check_sp("y12mae_c n=6 tridiag", N, b, ifail, 1e-5f);
}

/* 2. Double-precision black-box solve */
static void test_y12maf(void)
{
    double a[NNP], pivot[IHA], b[N], aflag[8];
    int    snr[NNP], rnr[NN1P], ha[IHA * 11], iflag[10];
    int    z, ifail;

    build_tridiag_dp(N, &z, a, snr, rnr, b);
    ifail = y12maf_c(N, z, a, snr, NNP, rnr, NN1P,
                   pivot, ha, IHA, aflag, iflag, b);
    check_dp("y12maf_c n=6 tridiag", N, b, ifail, 1e-10);
}

/* 3. Double-precision three-step API: y12mbf_c -> y12mcf_c -> y12mdf_c */
static void test_three_step_dp(void)
{
    double a[NNP], pivot[IHA], b[N], aflag[8];
    int    snr[NNP], rnr[NN1P], ha[IHA * 11], iflag[10];
    int    z, ifail;

    build_tridiag_dp(N, &z, a, snr, rnr, b);

    /* Initialise aflag and iflag to the same defaults y12maf_c would set. */
    memset(aflag, 0, sizeof(aflag));
    memset(iflag, 0, sizeof(iflag));
    aflag[0] = 16.0;  aflag[1] = 1e-12; aflag[2] = 1e16; aflag[3] = 1e-12;
    iflag[1] = 2;     iflag[2] = 1;     iflag[3] = 0;    iflag[4] = 1;

    ifail = y12mbf_c(N, z, a, snr, NNP, rnr, NN1P, ha, IHA, aflag, iflag);
    if (ifail == 0)
        ifail = y12mcf_c(N, z, a, snr, NNP, rnr, NN1P,
                       pivot, b, ha, IHA, aflag, iflag);
    if (ifail == 0)
        ifail = y12mdf_c(N, a, NNP, b, pivot, snr, ha, IHA, iflag);

    check_dp("y12mbf_c+y12mcf_c+y12mdf_c n=6 tridiag", N, b, ifail, 1e-10);
}

/* ----- main -------------------------------------------------------------- */

int main(void)
{
    test_y12mae();
    test_y12maf();
    test_three_step_dp();

    if (nfail != 0) {
        printf("%d test(s) FAILED\n", nfail);
        return 1;
    }
    printf("All test_y12m_c tests PASSED\n");
    return 0;
}
