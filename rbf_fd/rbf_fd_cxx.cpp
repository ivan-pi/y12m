// rbf_fd_cxx.cpp — C++20 implementation of the rbf_fd C interface.
// Requires C++20 (std::span). Compile with -std=c++20.
// Link with the same LAPACK flags as rbf_fd.c; do not link both object files.
#include "rbf_fd.h"
#include <array>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

// ---------------------------------------------------------------------------
// LAPACK backend (same compile-time switches as rbf_fd.c)
// ---------------------------------------------------------------------------
#ifdef RBF_USE_ACCELERATE
#  define ACCELERATE_NEW_LAPACK
#  include <Accelerate/Accelerate.h>

static int rbf__dgetrf(int n, double* A, int lda, int* ipiv) noexcept {
    int info;
    dgetrf_(&n, &n, A, &lda, ipiv, &info);
    return info;
}
static int rbf__dgetrs(int n, int nrhs, const double* A, int lda,
                       const int* ipiv, double* B, int ldb) noexcept {
    char trans = 'N';
    int info;
    // Fortran interface has no const; casts are safe.
    dgetrs_(&trans, &n, &nrhs, const_cast<double*>(A), &lda,
            const_cast<int*>(ipiv), B, &ldb, &info);
    return info;
}
static int rbf__dgecon(int n, const double* A, int lda,
                       double anorm, double* rcond) {
    char norm = '1';
    int info;
    std::vector<double> work(4 * n);
    std::vector<int>    iwork(n);
    dgecon_(&norm, &n, const_cast<double*>(A), &lda,
            &anorm, rcond, work.data(), iwork.data(), &info);
    return info;
}
static double rbf__dlange(int n, const double* A, int lda) noexcept {
    char norm = '1';
    double work{};  // not referenced for 1-norm
    return dlange_(&norm, &n, &n, const_cast<double*>(A), &lda, &work);
}

#else
#  include <lapacke.h>

static int rbf__dgetrf(int n, double* A, int lda, int* ipiv) noexcept {
    return LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, lda, ipiv);
}
static int rbf__dgetrs(int n, int nrhs, const double* A, int lda,
                       const int* ipiv, double* B, int ldb) noexcept {
    return LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, nrhs, A, lda, ipiv, B, ldb);
}
static int rbf__dgecon(int n, const double* A, int lda,
                       double anorm, double* rcond) noexcept {
    return LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', n, A, lda, anorm, rcond);
}
static double rbf__dlange(int n, const double* A, int lda) noexcept {
    return LAPACKE_dlange(LAPACK_COL_MAJOR, '1', n, n, A, lda);
}
#endif

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
constexpr int max_poly_degree = RBF_MAX_P;

// Non-owning column-major matrix view.
struct ColMajorView {
    double* data;
    int ld;
    double& operator()(int i, int j) noexcept       { return data[i + j * ld]; }
    double  operator()(int i, int j) const noexcept { return data[i + j * ld]; }
};

[[nodiscard]] static double ipow(double x, int n) noexcept {
    double r = 1.0;
    for (; n > 0; n >>= 1, x *= x)
        if (n & 1) r *= x;
    return r;
}

[[nodiscard]] static double phs_kernel(double rsqr, int q) noexcept {
    return ipow(rsqr, (q - 1) / 2) * std::sqrt(rsqr);
}

static void rbf_eval_basis(rbf_poly_t rbf, double xc, double yc,
    std::span<const double> x, std::span<const double> y,
    std::span<double> b) noexcept
{
    const int n = static_cast<int>(x.size());
    int k = 0;
    for (int i = 0; i < n; ++i) {
        const double dx = xc - x[i], dy = yc - y[i];
        b[k++] = phs_kernel(dx*dx + dy*dy, rbf.q);
    }

    std::array<double, max_poly_degree + 1> px, py;
    px[0] = py[0] = 1.0;
    for (int d = 1; d <= rbf.p; ++d) {
        px[d] = px[d-1] * xc;
        py[d] = py[d-1] * yc;
    }
    for (int d = 0; d <= rbf.p; ++d)
        for (int i = 0; i <= d; ++i)
            b[k++] = px[i] * py[d - i];
}

static void rbf_eval_der(rbf_poly_t rbf, rbf_deriv_t der,
    std::span<const double> x, std::span<const double> y,
    std::span<double> b) noexcept
{
    const int    n   = static_cast<int>(x.size());
    const int    ord = der.x + der.y;
    const double dq  = static_cast<double>(rbf.q);
    assert(ord <= 2);

    for (int i = 0; i < n; ++i) {
        const double xi = x[i], yi = y[i];
        const double rsqr = xi*xi + yi*yi;
        double val;
        if (rsqr == 0.0) {
            val = 0.0;
        } else {
            const double r    = std::sqrt(rsqr);
            const double phi  = ipow(rsqr, (rbf.q - 1) / 2) * r;
            const double phi1 = dq * phi / r;
            if (ord == 0) {
                val = phi;
            } else if (ord == 1) {
                val = phi1 * (der.x ? -xi : -yi) / r;
            } else {
                const double phi2 = (dq - 1.0) * phi1 / r;
                const double rx = -xi / r, ry = -yi / r;
                if (der.x == 2)
                    val = phi2*rx*rx + phi1*(rsqr - xi*xi)/(rsqr*r);
                else if (der.y == 2)
                    val = phi2*ry*ry + phi1*(rsqr - yi*yi)/(rsqr*r);
                else
                    val = phi2*rx*ry - phi1*xi*yi/(rsqr*r);
            }
        }
        b[i] = val;
    }

    const int m = rbf_poly_terms(rbf.p);
    std::fill(b.begin() + n, b.begin() + n + m, 0.0);
    const int d = der.x + der.y;
    if (d <= rbf.p) {
        double fac = 1.0;
        for (int i = 1; i <= der.x; ++i) fac *= i;
        for (int i = 1; i <= der.y; ++i) fac *= i;
        b[n + d*(d+1)/2 + der.x] = fac;
    }
}

static void rbf_build_matrix(rbf_poly_t rbf,
    std::span<const double> x, std::span<const double> y,
    ColMajorView M) noexcept
{
    const int n  = static_cast<int>(x.size());
    const int m  = rbf_poly_terms(rbf.p);
    const int nt = n + m;

    for (int i = 0; i < n; ++i)
        rbf_eval_basis(rbf, x[i], y[i], x, y,
                       {&M(0, i), static_cast<std::size_t>(nt)});

    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            M(i, n + j) = M(n + j, i);

    for (int j = 0; j < m; ++j)
        std::fill_n(&M(n, n + j), m, 0.0);
}

// ---------------------------------------------------------------------------
// Public C interface
// ---------------------------------------------------------------------------
extern "C" {

int rbf_factorize(rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[])
{
    const int m  = rbf_poly_terms(rbf.p);
    const int nt = n + m;
    if (rbf.q < 1 || rbf.q % 2 == 0)      return -1;
    if (rbf.p < 0 || rbf.p > max_poly_degree) return -1;
    if (n < m)                              return -3;
    if (ldm < nt)                           return -5;

    rbf_build_matrix(rbf, {x, static_cast<std::size_t>(n)},
                          {y, static_cast<std::size_t>(n)},
                     ColMajorView{M, ldm});
    const int info = rbf__dgetrf(nt, M, ldm, ipiv);
    assert(info >= 0);
    return info;
}

int rbf_factorize_rcond(rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, double M[], int ipiv[],
    double* rcond, double* anorm)
{
    const int m  = rbf_poly_terms(rbf.p);
    const int nt = n + m;
    if (rbf.q < 1 || rbf.q % 2 == 0)      return -1;
    if (rbf.p < 0 || rbf.p > max_poly_degree) return -1;
    if (n < m)                              return -3;
    if (ldm < nt)                           return -5;

    rbf_build_matrix(rbf, {x, static_cast<std::size_t>(n)},
                          {y, static_cast<std::size_t>(n)},
                     ColMajorView{M, ldm});
    const double norm = rbf__dlange(nt, M, ldm);
    if (anorm) *anorm = norm;

    int info = rbf__dgetrf(nt, M, ldm, ipiv);
    assert(info >= 0);
    if (info != 0) { *rcond = 0.0; return info; }

    info = rbf__dgecon(nt, M, ldm, norm, rcond);
    return info;
}

int rbf_stencil_weights(rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int num_ops, const rbf_deriv_t ops[],
    int ldw, double weights[])
{
    const int nt = rbf_system_size(rbf, n);
    if (ldw < nt) return -10;

    const std::span<const double> xs{x, static_cast<std::size_t>(n)};
    const std::span<const double> ys{y, static_cast<std::size_t>(n)};
    for (int q = 0; q < num_ops; ++q) {
        if (ops[q].x + ops[q].y > 2) return -9;
        rbf_eval_der(rbf, ops[q], xs, ys,
                     {weights + q * ldw, static_cast<std::size_t>(nt)});
    }
    const int info = rbf__dgetrs(nt, num_ops, M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int rbf_operator_weights(rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int num_terms, const rbf_deriv_t ders[], const double* coeffs,
    double weights[])
{
    const int nt = rbf_system_size(rbf, n);
    std::fill_n(weights, nt, 0.0);

    const std::span<const double> xs{x, static_cast<std::size_t>(n)};
    const std::span<const double> ys{y, static_cast<std::size_t>(n)};
    std::vector<double> wrk(nt);
    for (int t = 0; t < num_terms; ++t) {
        if (ders[t].x + ders[t].y > 2) return -9;
        rbf_eval_der(rbf, ders[t], xs, ys, wrk);
        const double c = coeffs ? coeffs[t] : 1.0;
        for (int i = 0; i < nt; ++i)
            weights[i] += c * wrk[i];
    }
    const int info = rbf__dgetrs(nt, 1, M, ldm, ipiv, weights, nt);
    assert(info == 0);
    return info;
}

int rbf_interpolation_weights(rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int np, const double xp[], const double yp[],
    int ldw, double weights[])
{
    const int nt = rbf_system_size(rbf, n);
    if (ldw < nt) return -11;

    const std::span<const double> xs{x, static_cast<std::size_t>(n)};
    const std::span<const double> ys{y, static_cast<std::size_t>(n)};
    for (int q = 0; q < np; ++q)
        rbf_eval_basis(rbf, xp[q], yp[q], xs, ys,
                       {weights + q * ldw, static_cast<std::size_t>(nt)});
    const int info = rbf__dgetrs(nt, np, M, ldm, ipiv, weights, ldw);
    assert(info == 0);
    return info;
}

int rbf_diff12_weights(rbf_poly_t rbf,
    int n, const double x[], const double y[],
    int ldm, const double M[], const int ipiv[],
    int ldw, double weights[])
{
    constexpr std::array<rbf_deriv_t, 5> ops{{{1,0},{0,1},{2,0},{1,1},{0,2}}};
    return rbf_stencil_weights(rbf, n, x, y, ldm, M, ipiv,
                               static_cast<int>(ops.size()), ops.data(), ldw, weights);
}

} // extern "C"
