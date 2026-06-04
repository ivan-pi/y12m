#ifndef RBF_FD_H
#define RBF_FD_H

// Returns the number of 2D polynomial basis terms up to degree p: (p+1)(p+2)/2.
static inline int rbf_poly_terms(int p) {
    return ((p + 1) * (p + 2)) / 2;
}

typedef struct {
    int q;  // PHS exponent, must be odd and >= 1
    int p;  // Maximum polynomial degree
} rbf_poly_t;

// Derivative operator d^(x+y) / dx^x dy^y.
typedef struct {
    int x;
    int y;
} rbf_deriv_t;

// Returns the total system size n + rbf_poly_terms(rbf.p).
static inline int rbf_system_size(rbf_poly_t rbf, int n) {
    return n + rbf_poly_terms(rbf.p);
}

// Builds and factorizes the augmented RBF-FD system matrix in-place.
//
// The stencil coordinates (x, y) must be in the local frame with the
// evaluation center at the origin.
//
// If rcond is non-NULL the reciprocal 1-norm condition number is written
// there on success (0.0 if the matrix is singular).
//
// Returns 0 on success, < 0 on argument error, > 0 if the matrix is singular.
int rbf_factorize(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[],
    double *rcond);

// Computes RBF-FD weights for num_ops derivative operators simultaneously.
//
// weights is column-major with leading dimension ldw >= rbf_system_size(rbf, n).
// After the call, weights[q*ldw + i] for i < n are the stencil weights for ops[q].
// Returns 0 on success, -9 if any operator has order x+y > 2, -10 if ldw < nt.
int rbf_stencil_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int num_ops, const rbf_deriv_t ops[],
    int ldw, double weights[]);

// Computes RBF-FD weights for a general linear differential operator
// L = sum_t coeffs[t] * D^{ders[t]}.
//
// If coeffs is NULL every coefficient defaults to 1.0.
// weights must hold at least rbf_system_size(rbf, n) doubles.
// Returns 0 on success, -9 if any derivative has order x+y > 2.
int rbf_operator_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int num_terms, const rbf_deriv_t ders[], const double *coeffs,
    double weights[]);

// Computes interpolation weights at np evaluation points (xp, yp).
//
// weights is column-major with leading dimension ldw >= rbf_system_size(rbf, n).
// The stencil coordinates (x, y) must match those used to build M.
// Returns 0 on success, -11 if ldw < nt.
int rbf_interpolation_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int np, const double xp[], const double yp[],
    int ldw, double weights[]);

// Convenience wrapper: weights for dx, dy, dxx, dxy, dyy (in that column order).
//
// weights must hold at least 5 * ldw doubles.
// Returns 0 on success, -10 if ldw < nt.
int rbf_diff12_weights(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int ldw, double weights[]);

#endif // RBF_FD_H
