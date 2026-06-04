#include "rbf_fd.h"
#include <assert.h>

// Evaluates the augmented RBF-FD basis (PHS + monomials) at (xc, yc),
// using stencil points (x[i], y[i]).  Output b[] has n + rbf_poly_terms(p)
// entries.
static void rbf_basis(const rbfpoly_t *rbf, double xc, double yc,
    int n, const double x[], const double y[], double b[])
{
    int k = 0;

    for (int i = 0; i < n; i++) {
        double rsqr = (xc - x[i])*(xc - x[i]) + (yc - y[i])*(yc - y[i]);
        b[k++] = phs_kernel(rsqr, rbf->q);
    }

    // Precompute monomial factors iteratively to avoid pow() calls.
    double px[MAX_P + 1], py[MAX_P + 1];
    px[0] = 1.0;  py[0] = 1.0;
    for (int d = 1; d <= rbf->p; d++) {
        px[d] = px[d-1] * xc;
        py[d] = py[d-1] * yc;
    }

    for (int d = 0; d <= rbf->p; d++) {
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
// RBF part uses the identity  r^(q-2) = r^q / r^2  to avoid negative
// exponents in phs_kernel for any valid odd q.
static void rbf_eval_der(const rbfpoly_t *rbf, const derivative_t *der,
    int n, const double x[], const double y[], double b[])
{
    int k = 0;
    double dq = (double)rbf->q;

    for (int i = 0; i < n; i++) {
        double xi = x[i], yi = y[i];
        double rsqr = xi*xi + yi*yi;
        double val;

        if (rsqr == 0.0) {
            val = 0.0;
        } else {
            double phi  = phs_kernel(rsqr, rbf->q); // r^q
            double phi2 = phi / rsqr;               // r^(q-2)

            if      (der->qx == 0 && der->qy == 0) {
                val = phi;
            } else if (der->qx == 1 && der->qy == 0) {
                // d/dx r^q at origin = -q * xi * r^(q-2)
                val = -dq * xi * phi2;
            } else if (der->qx == 0 && der->qy == 1) {
                val = -dq * yi * phi2;
            } else if (der->qx == 2 && der->qy == 0) {
                // d²/dx² r^q = q*r^(q-2) * (1 + (q-2)*xi²/r²)
                val = dq * phi2 * (1.0 + (dq - 2.0) * xi*xi / rsqr);
            } else if (der->qx == 0 && der->qy == 2) {
                val = dq * phi2 * (1.0 + (dq - 2.0) * yi*yi / rsqr);
            } else if (der->qx == 1 && der->qy == 1) {
                // d²/dxdy r^q = q*(q-2)*xi*yi * r^(q-4)
                //             = q*(q-2)*xi*yi * phi2/rsqr
                val = dq * (dq - 2.0) * xi * yi * phi2 / rsqr;
            } else {
                val = 0.0;  // higher order not implemented
            }
        }
        b[k++] = val;
    }

    // Polynomial part: d^qx d^qy (x^a y^b) at (0,0).
    // Non-zero only for the monomial (a,b) == (qx,qy); value is qx! * qy!.
    double fac = 1.0;
    for (int i = 1; i <= der->qx; i++) fac *= i;
    for (int i = 1; i <= der->qy; i++) fac *= i;

    for (int d = 0; d <= rbf->p; d++) {
        for (int i = 0; i <= d; i++) {
            b[k++] = (i == der->qx && (d - i) == der->qy) ? fac : 0.0;
        }
    }
}

int factorized_rbf_matrix(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[])
{
    const int m  = rbf_poly_terms(rbf->p);
    const int nt = n + m;

    if (rbf->q < 1 || rbf->q % 2 == 0) return -1;
    if (rbf->p < 0 || rbf->p > MAX_P)  return -1;
    if (n < m)                          return -3;
    if (ldm < nt)                       return -5;

#define IDX(i, j) ((i) + (j) * ldm)

    // Fill columns 0..n-1: each column i holds the basis centered at x[i].
    // Rows 0..n-1    → A block (RBF kernel matrix, symmetric)
    // Rows n..n+m-1  → P block (polynomial Vandermonde)
    for (int i = 0; i < n; i++) {
        rbf_basis(rbf, x[i], y[i], n, x, y, &M[IDX(0, i)]);
    }

    // Transpose P into the top-right P^T block (columns n..n+m-1, rows 0..n-1).
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            M[IDX(i, n + j)] = M[IDX(n + j, i)];
        }
    }

    // Zero the bottom-right m×m block.
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

int compute_derivative_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_ops, const derivative_t ops[],
    int ldw, double weights[])
{
    const int nt = n + rbf_poly_terms(rbf->p);

    for (int q = 0; q < num_ops; q++) {
        rbf_eval_der(rbf, &ops[q], n, x, y, &weights[q * ldw]);
    }

    int info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nt, num_ops,
                              M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int compute_operator_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_terms, const derivative_t ders[], const double coeffs[],
    double weights[])
{
    const int nt = n + rbf_poly_terms(rbf->p);

    for (int i = 0; i < nt; i++) weights[i] = 0.0;

    double wrk[nt];  // VLA; nt is typically small (< 200)
    for (int t = 0; t < num_terms; t++) {
        rbf_eval_der(rbf, &ders[t], n, x, y, wrk);
        double c = coeffs ? coeffs[t] : 1.0;
        for (int i = 0; i < nt; i++) {
            weights[i] += c * wrk[i];
        }
    }

    int info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nt, 1, M, ldm, ipiv, weights, nt);
    assert(info == 0);
    return info;
}

int compute_interpolation_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int npoints, const double xp[], const double yp[],
    int ldw, double weights[])
{
    const int nt = n + rbf_poly_terms(rbf->p);

    for (int q = 0; q < npoints; q++) {
        rbf_basis(rbf, xp[q], yp[q], n, x, y, &weights[q * ldw]);
    }

    int info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nt, npoints,
                              M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int compute_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int ldw, double weights[])
{
    static const derivative_t ops[5] = {{1,0}, {0,1}, {2,0}, {1,1}, {0,2}};
    return compute_derivative_weights(rbf, n, x, y, ldm, M, ipiv, 5, ops, ldw, weights);
}
