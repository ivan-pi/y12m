#include "rbf_fd.h"
#include <assert.h>

#if _OPENMP
#  include <omp.h>
#else
#  include <time.h>
#endif

static double walltime(void) {
#if _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
}

// ---------------------------------------------------------------------------
// Generic assembly loop
// ---------------------------------------------------------------------------

int rbf_assemble(
    int nstencils, const rbf_kernel_t *kernel, int n,
    rbf_gather_fn gather, void *gather_data,
    int nops, rbf_weights_fn weights_fn, void *weights_params,
    rbf_pack_fn pack, void *pack_data,
    double *elapsed)
{
    if (kernel->type != RBF_PHS)                       return -2;
    if (kernel->phs.q < 1 || kernel->phs.q % 2 == 0)  return -1;
    if (kernel->p < 0 || kernel->p > RBF_MAX_P)        return -1;

    rbf_poly_t rbf = { .q = kernel->phs.q, .p = kernel->p };

    if (n < rbf_poly_terms(rbf.p))                     return -3;

    const int nt  = rbf_system_size(rbf, n);
    const int ldm = ((nt + RBF_LD_ALIGN - 1) / RBF_LD_ALIGN) * RBF_LD_ALIGN;
    const int ldw = ldm;

    const double t0 = walltime();

#pragma omp parallel
    {
        double M[ldm * nt], W[ldw * nops];
        int    ipiv[nt];
        double x[n], y[n];
        double eps;   // per-stencil shape parameter; ignored for PHS

#pragma omp for schedule(static)
        for (int i = 0; i < nstencils; i++) {
            gather(i, n, x, y, &eps, gather_data);
            int ierr = rbf_factorize(rbf, n, x, y, ldm, M, ipiv);
            assert(ierr == 0);
            weights_fn(rbf, n, x, y, ldm, M, ipiv, ldw, W, weights_params);
            pack(i, n, nops, W, ldw, pack_data);
        }
    }

    if (elapsed) *elapsed = walltime() - t0;
    return 0;
}

// ---------------------------------------------------------------------------
// Weight callbacks and convenience wrappers
// ---------------------------------------------------------------------------

static void cb_diff12(rbf_poly_t rbf, int n,
    const double *x, const double *y,
    int ldm, const double *M, const int *ipiv,
    int ldw, double *W, void *params)
{
    (void)params;
    rbf_diff12_weights(rbf, n, x, y, ldm, M, ipiv, ldw, W);
}

int rbf_assemble_diff12(
    int nstencils, const rbf_kernel_t *kernel, int n,
    rbf_gather_fn gather, void *gather_data,
    rbf_pack_fn pack, void *pack_data,
    double *elapsed)
{
    return rbf_assemble(nstencils, kernel, n,
        gather, gather_data,
        5, cb_diff12, NULL,
        pack, pack_data, elapsed);
}

// ---- stencil weights -------------------------------------------------------

typedef struct { int num_ops; const rbf_deriv_t *ops; } stencil_params_t;

static void cb_stencil(rbf_poly_t rbf, int n,
    const double *x, const double *y,
    int ldm, const double *M, const int *ipiv,
    int ldw, double *W, void *params)
{
    stencil_params_t *p = params;
    rbf_stencil_weights(rbf, n, x, y, ldm, M, ipiv, p->num_ops, p->ops, ldw, W);
}

int rbf_assemble_stencil(
    int nstencils, const rbf_kernel_t *kernel, int n,
    rbf_gather_fn gather, void *gather_data,
    int num_ops, const rbf_deriv_t ops[],
    rbf_pack_fn pack, void *pack_data,
    double *elapsed)
{
    stencil_params_t p = {num_ops, ops};
    return rbf_assemble(nstencils, kernel, n,
        gather, gather_data,
        num_ops, cb_stencil, &p,
        pack, pack_data, elapsed);
}

// ---- operator weights ------------------------------------------------------

typedef struct {
    int num_terms;
    const rbf_deriv_t *ders;
    const double *coeffs;
} operator_params_t;

static void cb_operator(rbf_poly_t rbf, int n,
    const double *x, const double *y,
    int ldm, const double *M, const int *ipiv,
    int ldw, double *W, void *params)
{
    (void)ldw;
    operator_params_t *p = params;
    rbf_operator_weights(rbf, n, x, y, ldm, M, ipiv,
                         p->num_terms, p->ders, p->coeffs, W);
}

int rbf_assemble_operator(
    int nstencils, const rbf_kernel_t *kernel, int n,
    rbf_gather_fn gather, void *gather_data,
    int num_terms, const rbf_deriv_t ders[], const double *coeffs,
    rbf_pack_fn pack, void *pack_data,
    double *elapsed)
{
    operator_params_t p = {num_terms, ders, coeffs};
    return rbf_assemble(nstencils, kernel, n,
        gather, gather_data,
        1, cb_operator, &p,
        pack, pack_data, elapsed);
}

// ---- interpolation weights -------------------------------------------------

typedef struct {
    int np;
    const double *xp;
    const double *yp;
} interp_params_t;

static void cb_interp(rbf_poly_t rbf, int n,
    const double *x, const double *y,
    int ldm, const double *M, const int *ipiv,
    int ldw, double *W, void *params)
{
    interp_params_t *p = params;
    rbf_interpolation_weights(rbf, n, x, y, ldm, M, ipiv,
                              p->np, p->xp, p->yp, ldw, W);
}

int rbf_assemble_interp(
    int nstencils, const rbf_kernel_t *kernel, int n,
    rbf_gather_fn gather, void *gather_data,
    int np, const double xp[], const double yp[],
    rbf_pack_fn pack, void *pack_data,
    double *elapsed)
{
    interp_params_t p = {np, xp, yp};
    return rbf_assemble(nstencils, kernel, n,
        gather, gather_data,
        np, cb_interp, &p,
        pack, pack_data, elapsed);
}
