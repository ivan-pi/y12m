#ifndef RBF_FD_H
#define RBF_FD_H

#include <math.h>
#include <lapacke.h>

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

// Returns the number of 2D polynomial terms up to degree p: (p+1)(p+2)/2.
static inline int rbf_poly_terms(int p) {
    return ((p + 1) * (p + 2)) / 2;
}

typedef struct {
    int q;  // PHS exponent, must be odd and >= 1
    int p;  // Maximum polynomial degree, 0 <= p <= MAX_P
} rbfpoly_t;

// A derivative operator: d^qx/dx^qx * d^qy/dy^qy
typedef struct {
    int qx;
    int qy;
} derivative_t;

// Builds and factorizes the augmented RBF-FD weight matrix.
//
// Inputs: stencil coordinates (x, y) of length n, already in local frame
//         with the evaluation center at the origin.
//
// Returns 0 on success, < 0 on argument error, > 0 if the matrix is singular.
int factorized_rbf_matrix(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[]);

// Computes RBF-FD weights for num_ops derivative operators simultaneously.
//
// weights must have space for num_ops * ldw doubles (column-major, ldw >= nt).
// After the call, weights[q*ldw + i] (i < n) are the stencil weights for ops[q].
int compute_derivative_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_ops, const derivative_t ops[],
    int ldw, double weights[]);

// Computes RBF-FD weights for a general linear differential operator
// L = sum_{t} coeffs[t] * D^{ders[t]}.
//
// If coeffs is NULL every term is weighted 1.0.
// weights must have space for at least (n + rbf_poly_terms(p)) doubles.
int compute_operator_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_terms, const derivative_t ders[], const double coeffs[],
    double weights[]);

// Computes interpolation weights at npoints evaluation points (xp, yp).
//
// weights must have space for npoints * ldw doubles (column-major, ldw >= nt).
// The stencil coordinates (x, y) must match those used to build M.
int compute_interpolation_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int npoints, const double xp[], const double yp[],
    int ldw, double weights[]);

// Convenience wrapper: computes weights for the standard five operators
// dx, dy, dxx, dxy, dyy (in that column order).
//
// weights must have space for 5 * ldw doubles.
int compute_weights(
    const rbfpoly_t *rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int ldw, double weights[]);

#endif // RBF_FD_H
