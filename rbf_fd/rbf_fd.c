#include "rbf_fd.h"
#include <assert.h>
#include <lapacke.h>
#include <math.h>

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

static void rbf_basis(rbfpoly_t rbf, double xc, double yc,
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
// The stencil coordinates (x, y) must be in the local frame with the
// evaluation center at the origin.
//
// Supported derivative orders: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2).
// Higher-order derivatives write 0 (not implemented).
//
// Chain rule for a radial function phi(r):
//   dphi/dx   = phi' * dr/dx
//   d^2phi/dx^2 = phi'' * (dr/dx)^2 + phi' * d^2r/dx^2
// where phi' = q*r^(q-1), phi'' = q*(q-1)*r^(q-2),
// dr/dx = -xi/r, d^2r/dx^2 = (r^2 - xi^2)/r^3, d^2r/dxdy = -xi*yi/r^3.
static void rbf_eval_der(rbfpoly_t rbf, derivative_t der,
    int n, const double x[], const double y[], double b[])
{
    int k = 0;
    const int ord = der.qx + der.qy;
    const double dq = (double)rbf.q;

    for (int i = 0; i < n; i++) {
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
                val = phi1 * (der.qx ? -xi : -yi) / r;
            } else {
                double phi2 = (dq - 1.0) * phi1 / r;  // phi'' = q*(q-1)*r^(q-2)
                double rx = -xi / r, ry = -yi / r;
                if (der.qx == 2) {
                    // d^2r/dx^2 = (r^2 - xi^2) / r^3
                    val = phi2*rx*rx + phi1*(rsqr - xi*xi)/(rsqr*r);
                } else if (der.qy == 2) {
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

    // Polynomial part: D^der (x^a * y^b) at (0,0).
    // Non-zero only for the monomial (a,b) == (qx,qy); value is qx! * qy!.
    double fac = 1.0;
    for (int i = 1; i <= der.qx; i++) fac *= i;
    for (int i = 1; i <= der.qy; i++) fac *= i;

    for (int d = 0; d <= rbf.p; d++) {
        for (int i = 0; i <= d; i++) {
            b[k++] = (i == der.qx && (d - i) == der.qy) ? fac : 0.0;
        }
    }
}

int rbf_factorize(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[])
{
    const int m  = rbf_poly_terms(rbf.p);
    const int nt = n + m;

    if (rbf.q < 1 || rbf.q % 2 == 0) return -1;
    if (rbf.p < 0 || rbf.p > MAX_P)  return -1;
    if (n < m)                        return -3;
    if (ldm < nt)                     return -5;

#define IDX(i, j) ((i) + (j) * ldm)

    // Fill columns 0..n-1: each column i holds the basis centered at x[i].
    // Rows 0..n-1   -> A block (RBF kernel matrix, symmetric)
    // Rows n..n+m-1 -> P block (polynomial Vandermonde)
    for (int i = 0; i < n; i++) {
        rbf_basis(rbf, x[i], y[i], n, x, y, &M[IDX(0, i)]);
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

    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, nt, nt, M, ldm, ipiv);
    assert(info >= 0);  // info < 0 is a LAPACK argument error
    return info;
}

int rbf_derivative_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_ops, const derivative_t ops[],
    int ldw, double weights[])
{
    const int nt = rbf_system_size(rbf, n);

    if (ldw < nt) return -10;

    for (int q = 0; q < num_ops; q++) {
        rbf_eval_der(rbf, ops[q], n, x, y, &weights[q * ldw]);
    }

    int info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nt, num_ops,
                              M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int rbf_operator_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_terms, const derivative_t ders[], const double *coeffs,
    double weights[])
{
    const int nt = rbf_system_size(rbf, n);

    for (int i = 0; i < nt; i++) weights[i] = 0.0;

    double wrk[nt];  // VLA; nt is typically small (< 200)
    for (int t = 0; t < num_terms; t++) {
        rbf_eval_der(rbf, ders[t], n, x, y, wrk);
        double c = coeffs ? coeffs[t] : 1.0;
        for (int i = 0; i < nt; i++) {
            weights[i] += c * wrk[i];
        }
    }

    int info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nt, 1, M, ldm, ipiv, weights, nt);
    assert(info == 0);
    return info;
}

int rbf_interpolation_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int np, const double xp[], const double yp[],
    int ldw, double weights[])
{
    const int nt = rbf_system_size(rbf, n);

    if (ldw < nt) return -11;

    for (int q = 0; q < np; q++) {
        rbf_basis(rbf, xp[q], yp[q], n, x, y, &weights[q * ldw]);
    }

    int info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nt, np,
                              M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int rbf_standard_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int ldw, double weights[])
{
    static const derivative_t ops[5] = {{1,0}, {0,1}, {2,0}, {1,1}, {0,2}};
    return rbf_derivative_weights(rbf, n, x, y, ldm, M, ipiv, 5, ops, ldw, weights);
}
