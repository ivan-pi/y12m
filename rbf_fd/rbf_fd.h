#ifndef RBF_FD_H
#define RBF_FD_H

// Returns the number of 2D polynomial basis terms up to degree p: (p+1)(p+2)/2.
static inline int rbf_poly_terms(int p) {
    return ((p + 1) * (p + 2)) / 2;
}

typedef struct {
    int q;  // PHS exponent, must be odd and >= 1
    int p;  // Maximum polynomial degree
} rbfpoly_t;

typedef struct {
    int qx;
    int qy;
} derivative_t;

// Returns the total system size n + rbf_poly_terms(rbf.p).
static inline int rbf_system_size(rbfpoly_t rbf, int n) {
    return n + rbf_poly_terms(rbf.p);
}

// Builds and factorizes the augmented RBF-FD system matrix in-place.
//
// The stencil coordinates (x, y) must be in the local frame with the
// evaluation center at the origin.
//
// Returns 0 on success, < 0 on argument error, > 0 if the matrix is singular.
int rbf_factorize(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[]);

// Computes RBF-FD weights for num_ops derivative operators simultaneously.
//
// weights is column-major with leading dimension ldw >= rbf_system_size(rbf, n).
// After the call, weights[q*ldw + i] for i < n are the stencil weights for ops[q].
int rbf_derivative_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_ops, const derivative_t ops[],
    int ldw, double weights[]);

// Computes RBF-FD weights for a general linear differential operator
// L = sum_t coeffs[t] * D^{ders[t]}.
//
// If coeffs is NULL every coefficient defaults to 1.0.
// weights must hold at least rbf_system_size(rbf, n) doubles.
int rbf_operator_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int num_terms, const derivative_t ders[], const double *coeffs,
    double weights[]);

// Computes interpolation weights at np evaluation points (xp, yp).
//
// weights is column-major with leading dimension ldw >= rbf_system_size(rbf, n).
// The stencil coordinates (x, y) must match those used to build M.
int rbf_interpolation_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int np, const double xp[], const double yp[],
    int ldw, double weights[]);

// Convenience wrapper: weights for the standard five operators
// dx, dy, dxx, dxy, dyy (in that column order).
//
// weights must hold at least 5 * ldw doubles.
int rbf_standard_weights(
    rbfpoly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int *ipiv,
    int ldw, double weights[]);

#endif // RBF_FD_H
