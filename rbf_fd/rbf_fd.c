#include "rbf_fd.h"
#include <assert.h>
#include <math.h>

// ---------------------------------------------------------------------------
// LAPACK backend
// Compile with -DRBF_USE_ACCELERATE to use Apple Accelerate (no LAPACKE).
// Default: LAPACKE interface (OpenBLAS, MKL, LAPACKE standalone, ...).
//
// rbf__dgetrf: LU factorize an n*n matrix (overwrites A, fills ipiv).
// rbf__dgetrs: solve A*X = B using a prior LU factorization (overwrites B).
// rbf__dgecon: estimate rcond = 1/(||A||_1 * ||A^-1||_1) from LU factors.
// All return 0 on success, > 0 on singular matrix, < 0 on argument error.
// ---------------------------------------------------------------------------
#ifdef RBF_USE_ACCELERATE
#  define ACCELERATE_NEW_LAPACK
#  include <Accelerate/Accelerate.h>
static inline int rbf__dgetrf(int n, double *A, int lda, int *ipiv) {
    int info;
    dgetrf_(&n, &n, A, &lda, ipiv, &info);
    return info;
}
static inline int rbf__dgetrs(int n, int nrhs, const double *A, int lda,
                               const int *ipiv, double *B, int ldb) {
    char trans = 'N';
    int info;
    // Fortran interface has no const; casts are safe (routine does not modify A or ipiv).
    dgetrs_(&trans, &n, &nrhs, (double *)A, &lda, (int *)ipiv, B, &ldb, &info);
    return info;
}
static inline int rbf__dgecon(int n, const double *A, int lda,
                               double anorm, double *rcond) {
    char norm = '1';
    int info;
    double work[4*n];   // VLA; n is typically small (< 200)
    int    iwork[n];    // VLA
    dgecon_(&norm, &n, (double *)A, &lda, &anorm, rcond, work, iwork, &info);
    return info;
}
static inline double rbf__dlange(int n, const double *A, int lda) {
    char norm = '1';
    double work;  // not referenced for 1-norm
    return dlange_(&norm, &n, &n, (double *)A, &lda, &work);
}
#else
#  include <lapacke.h>
static inline int rbf__dgetrf(int n, double *A, int lda, int *ipiv) {
    return LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, lda, ipiv);
}
static inline int rbf__dgetrs(int n, int nrhs, const double *A, int lda,
                               const int *ipiv, double *B, int ldb) {
    return LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, nrhs, A, lda, ipiv, B, ldb);
}
static inline int rbf__dgecon(int n, const double *A, int lda,
                               double anorm, double *rcond) {
    return LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', n, A, lda, anorm, rcond);
}
static inline double rbf__dlange(int n, const double *A, int lda) {
    return LAPACKE_dlange(LAPACK_COL_MAJOR, '1', n, n, A, lda);
}
#endif

#define MAX_P 10  // Maximum supported polynomial degree

// Integer power by repeated squaring: x^n, n >= 0.
static inline double ipow(double x, int n) {
    double result = 1.0;
    while (n > 0) {
        if (n & 1) result *= x;
        x *= x;
        n >>= 1;
    }
    return result;
}

// Evaluates the PHS kernel phi(r) = r^q, given r^2 and odd q >= 1.
//
// Decomposition: r^q = (r^2)^((q-1)/2) * r
//
// Cost vs. pow(sqrt(rsqr), q):
//   q=3:  1 sqrt, 1 mul
//   q=5:  1 sqrt, 2 muls
//   q=7:  1 sqrt, 3 muls
//   q=9:  1 sqrt, 3 muls  (repeated squaring: (9-1)/2 = 4 = 100b)
//
// Handles rsqr == 0.0 correctly: returns 0.0 for all valid q.
static inline double phs_kernel(double rsqr, int q) {
    return ipow(rsqr, (q - 1) / 2) * sqrt(rsqr);
}

static void rbf_eval_basis(rbf_poly_t rbf, double xc, double yc,
    int n, const double x[], const double y[], double b[])
{
    int k = 0;
    for (int i = 0; i < n; i++) {
        double rsqr = (xc - x[i])*(xc - x[i]) + (yc - y[i])*(yc - y[i]);
        b[k++] = phs_kernel(rsqr, rbf.q);
    }

    double px[MAX_P + 1], py[MAX_P + 1];
    px[0] = 1.0; py[0] = 1.0;
    for (int d = 1; d <= rbf.p; d++) {
        px[d] = px[d-1] * xc;
        py[d] = py[d-1] * yc;
    }

    for (int d = 0; d <= rbf.p; d++) {
        for (int i = 0; i <= d; i++) {
            b[k++] = px[i] * py[d - i];
        }
    }
}

// Evaluates the derivative D^der of the augmented basis at the origin (0,0).
//
// x[i], y[i] are stencil point positions in the local frame (xc = 0), so the
// displacement from stencil point i to the center is (0 - x[i], 0 - y[i]),
// which is why dr/dx = (xc - x[i])/r = -x[i]/r.
//
// Chain rule for a radial function phi(r), r = ||(xc,yc) - (x[i],y[i])||:
//   dphi/dx_k   = phi' * dr/dx_k
//   d^2phi/dx_k dx_l = phi'' * (dr/dx_k)*(dr/dx_l) + phi' * d^2r/dx_k dx_l
// where phi' = q*r^(q-1), phi'' = q*(q-1)*r^(q-2),
// dr/dx = -xi/r, d^2r/dx^2 = (r^2-xi^2)/r^3, d^2r/dxdy = -xi*yi/r^3.
//
// Supported derivative orders: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2).
// Higher-order entries are not implemented and write 0.
static void rbf_eval_der(rbf_poly_t rbf, rbf_deriv_t der,
    int n, const double x[], const double y[], double b[])
{
    int k = 0;
    const int ord = der.x + der.y;
    assert(ord <= 2);  // caller must validate; higher orders are not implemented
    const double dq = (double)rbf.q;

    for (int i = 0; i < n; i++) {
        // xi, yi: displacement from stencil point to center (= -stencil position)
        double xi = x[i], yi = y[i];
        double rsqr = xi*xi + yi*yi;
        double val;

        if (rsqr == 0.0) {
            val = 0.0;
        } else {
            double r    = sqrt(rsqr);
            double phi  = ipow(rsqr, (rbf.q - 1) / 2) * r;  // r^q
            double phi1 = dq * phi / r;                       // phi' = q*r^(q-1)

            if (ord == 0) {
                val = phi;
            } else if (ord == 1) {
                // dr/dx = -xi/r,  dr/dy = -yi/r
                val = phi1 * (der.x ? -xi : -yi) / r;
            } else {
                double phi2 = (dq - 1.0) * phi1 / r;  // phi'' = q*(q-1)*r^(q-2)
                double rx = -xi / r, ry = -yi / r;
                if (der.x == 2) {
                    // d^2r/dx^2 = (r^2 - xi^2) / r^3
                    val = phi2*rx*rx + phi1*(rsqr - xi*xi)/(rsqr*r);
                } else if (der.y == 2) {
                    // d^2r/dy^2 = (r^2 - yi^2) / r^3
                    val = phi2*ry*ry + phi1*(rsqr - yi*yi)/(rsqr*r);
                } else {
                    // d^2r/dxdy = -xi*yi / r^3
                    val = phi2*rx*ry - phi1*xi*yi/(rsqr*r);
                }
            }
        }
        b[k++] = val;
    }

    // Polynomial part: D^der (x^a * y^b) at (0,0) is non-zero only when (a,b)==(der.x,der.y).
    // Monomial x^a * y^b with a+b = d occupies flat index d*(d+1)/2 + a.
    const int m = rbf_poly_terms(rbf.p);
    for (int j = 0; j < m; j++) b[k + j] = 0.0;

    const int d = der.x + der.y;
    if (d <= rbf.p) {
        double fac = 1.0;
        for (int i = 1; i <= der.x; i++) fac *= i;
        for (int i = 1; i <= der.y; i++) fac *= i;
        b[k + d*(d+1)/2 + der.x] = fac;
    }
}

static void rbf_build_matrix(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[])
{
    const int m = rbf_poly_terms(rbf.p);

#define IDX(i, j) ((i) + (j) * ldm)

    // Fill columns 0..n-1: each column i holds the basis centered at x[i].
    // Rows 0..n-1   -> A block (RBF kernel matrix, symmetric)
    // Rows n..n+m-1 -> P block (polynomial Vandermonde)
    for (int i = 0; i < n; i++) {
        rbf_eval_basis(rbf, x[i], y[i], n, x, y, &M[IDX(0, i)]);
    }

    // Transpose P into the top-right P^T block (columns n..n+m-1, rows 0..n-1).
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            M[IDX(i, n + j)] = M[IDX(n + j, i)];
        }
    }

    // Zero the bottom-right m*m block.
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            M[IDX(n + i, n + j)] = 0.0;
        }
    }

#undef IDX
}

int rbf_factorize(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[])
{
    const int m  = rbf_poly_terms(rbf.p);
    const int nt = n + m;

    if (rbf.q < 1 || rbf.q % 2 == 0) return -1;
    if (rbf.p < 0 || rbf.p > MAX_P)  return -1;
    if (n < m)                        return -3;
    if (ldm < nt)                     return -5;

    rbf_build_matrix(rbf, n, x, y, ldm, M);

    int info = rbf__dgetrf(nt, M, ldm, ipiv);
    assert(info >= 0);  // info < 0 is a LAPACK argument error
    return info;
}

int rbf_factorize_rcond(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[],
    double *rcond, double *anorm)
{
    const int m  = rbf_poly_terms(rbf.p);
    const int nt = n + m;

    if (rbf.q < 1 || rbf.q % 2 == 0) return -1;
    if (rbf.p < 0 || rbf.p > MAX_P)  return -1;
    if (n < m)                        return -3;
    if (ldm < nt)                     return -5;

    rbf_build_matrix(rbf, n, x, y, ldm, M);

    // Compute 1-norm before dgetrf overwrites M.
    double norm = rbf__dlange(nt, M, ldm);
    if (anorm) *anorm = norm;

    int info = rbf__dgetrf(nt, M, ldm, ipiv);
    assert(info >= 0);
    if (info != 0) { *rcond = 0.0; return info; }

    info = rbf__dgecon(nt, M, ldm, norm, rcond);
    return info;
}

int rbf_stencil_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int num_ops, const rbf_deriv_t ops[],
    int ldw, double weights[])
{
    const int nt = rbf_system_size(rbf, n);

    if (ldw < nt) return -10;
    for (int q = 0; q < num_ops; q++) {
        if (ops[q].x + ops[q].y > 2) return -9;
        rbf_eval_der(rbf, ops[q], n, x, y, &weights[q * ldw]);
    }

    int info = rbf__dgetrs(nt, num_ops, M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int rbf_operator_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int num_terms, const rbf_deriv_t ders[], const double *coeffs,
    double weights[])
{
    const int nt = rbf_system_size(rbf, n);

    for (int i = 0; i < nt; i++) weights[i] = 0.0;

    double wrk[nt];  // VLA; nt is typically small (< 200)
    for (int t = 0; t < num_terms; t++) {
        if (ders[t].x + ders[t].y > 2) return -9;
        rbf_eval_der(rbf, ders[t], n, x, y, wrk);
        double c = coeffs ? coeffs[t] : 1.0;
        for (int i = 0; i < nt; i++) {
            weights[i] += c * wrk[i];
        }
    }

    int info = rbf__dgetrs(nt, 1, M, ldm, ipiv, weights, nt);
    assert(info == 0);
    return info;
}

int rbf_interpolation_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int np, const double xp[], const double yp[],
    int ldw, double weights[])
{
    const int nt = rbf_system_size(rbf, n);

    if (ldw < nt) return -11;

    for (int q = 0; q < np; q++) {
        rbf_eval_basis(rbf, xp[q], yp[q], n, x, y, &weights[q * ldw]);
    }

    int info = rbf__dgetrs(nt, np, M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int rbf_diff12_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int ldw, double weights[])
{
    static const rbf_deriv_t ops[5] = {{1,0}, {0,1}, {2,0}, {1,1}, {0,2}};
    return rbf_stencil_weights(rbf, n, x, y, ldm, M, ipiv, 5, ops, ldw, weights);
}
