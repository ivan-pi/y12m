#ifndef RBF_FD_H
#define RBF_FD_H

// Maximum supported polynomial augmentation degree.
// Can be raised at compile time (-DRBF_MAX_P=12), but larger values
// increase the VLA sizes in the implementation.
#ifndef RBF_MAX_P
#  define RBF_MAX_P 10
#endif

// Returns RBF_MAX_P as a runtime value.  Useful for FFI callers (Fortran,
// Python ctypes, Julia ccall) that cannot access C preprocessor macros.
inline int rbf_max_p(void) { return RBF_MAX_P; }

// Returns the number of 2D polynomial basis terms up to degree p: (p+1)(p+2)/2.
static inline int rbf_poly_terms(int p) {
    return ((p + 1) * (p + 2)) / 2;
}

typedef struct {
    int q;  // PHS exponent, must be odd and >= 1
    int p;  // Maximum polynomial degree, 0 <= p <= RBF_MAX_POLY_DEGREE
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
// Returns 0 on success, < 0 on argument error, > 0 if the matrix is singular.
int rbf_factorize(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[]);

// Builds and factorizes the matrix, then estimates the reciprocal 1-norm
// condition number via dgecon.
//
// rcond is required; anorm receives the 1-norm of the unmodified matrix
// (pass NULL to discard).  On a singular matrix *rcond is set to 0.0.
//
// Returns 0 on success, < 0 on argument error, > 0 if the matrix is singular.
int rbf_factorize_rcond(
    rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[],
    double *rcond, double *anorm);

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

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

// Leading-dimension alignment in doubles: ldm is always a multiple of this,
// keeping each matrix column on a SIMD-friendly boundary.
// Override at compile time: -DRBF_LD_ALIGN=16 for AVX-512, =1 to disable.
#ifndef RBF_LD_ALIGN
#  define RBF_LD_ALIGN 8
#endif

// Fills x[n] and y[n] with the stencil coordinates for stencil i in its
// local frame (evaluation centre at the origin).
typedef void (*rbf_gather_fn)(int i, int n, double *x, double *y, void *data);

// Stores the computed weights W (column-major, nops columns of length ldw)
// for stencil i into user storage.  W must be treated as read-only.
typedef void (*rbf_pack_fn)(int i, int n, int nops,
                            const double *W, int ldw, void *data);

// Computes nops columns of weights into W[ldw * nops] given a factored system.
// params carries any extra data (operator specs, eval points, etc.).
typedef void (*rbf_weights_fn)(rbf_poly_t rbf, int n,
                               const double *x, const double *y,
                               int ldm, const double *M, const int *ipiv,
                               int ldw, double *W, void *params);

// Generic assembly loop: factorizes each stencil, calls weights_fn to fill W,
// then calls pack to store results.  nops is the number of weight columns
// (determines the width of the W buffer passed to weights_fn and pack).
//
// Uses one thread per OpenMP thread (or one thread if OpenMP is absent).
// If elapsed is non-NULL, *elapsed receives the wall-clock time in seconds.
// Returns 0 on success.
// Returns -1 if rbf.q is not a positive odd integer or rbf.p is out of range.
// Returns -3 if n < rbf_poly_terms(rbf.p) (too few stencil points for the
//   polynomial degree; the system would be underdetermined).
// Aborts via assert if a stencil matrix is singular (ill-conditioned points).
int rbf_assemble(
    int nstencils, rbf_poly_t rbf, int n,
    rbf_gather_fn gather, void *gather_data,
    int nops, rbf_weights_fn weights_fn, void *weights_params,
    rbf_pack_fn pack, void *pack_data,
    double *elapsed);

// Convenience wrappers around rbf_assemble for the three standard weight modes.

// Weights for dx, dy, dxx, dxy, dyy (nops = 5).
int rbf_assemble_diff12(
    int nstencils, rbf_poly_t rbf, int n,
    rbf_gather_fn gather, void *gather_data,
    rbf_pack_fn pack, void *pack_data,
    double *elapsed);

// Weights for num_ops arbitrary derivative operators (nops = num_ops).
int rbf_assemble_stencil(
    int nstencils, rbf_poly_t rbf, int n,
    rbf_gather_fn gather, void *gather_data,
    int num_ops, const rbf_deriv_t ops[],
    rbf_pack_fn pack, void *pack_data,
    double *elapsed);

// Weights for a single linear operator L = sum_t coeffs[t] * D^{ders[t]} (nops = 1).
// If coeffs is NULL every coefficient defaults to 1.0.
int rbf_assemble_operator(
    int nstencils, rbf_poly_t rbf, int n,
    rbf_gather_fn gather, void *gather_data,
    int num_terms, const rbf_deriv_t ders[], const double *coeffs,
    rbf_pack_fn pack, void *pack_data,
    double *elapsed);

#endif // RBF_FD_H
