#include "rbf_fd.h"
#include <assert.h>

#if _OPENMP
#  include <omp.h>
#else
#  include <time.h>
#endif

static double wall_time(void) {
#if _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
}

int rbf_assemble_diff12(
    int nstencils, rbf_poly_t rbf, int n,
    rbf_gather_fn gather, void *gather_data,
    rbf_pack_fn pack, void *pack_data,
    double *elapsed)
{
    const int nt   = rbf_system_size(rbf, n);
    const int vlen = 8;
    const int ldm  = ((nt + vlen - 1) / vlen) * vlen;
    const int ldw  = ldm;

    const double t0 = wall_time();

#pragma omp parallel
    {
        // Declared inside the parallel region: each thread gets its own
        // stack allocation without any firstprivate copying.
        double M[ldm * nt], W[ldw * 5];
        int    ipiv[nt];
        double x[n], y[n];

#pragma omp for schedule(static)
        for (int i = 0; i < nstencils; i++) {
            gather(i, n, x, y, gather_data);
            int ierr = rbf_factorize(rbf, n, x, y, ldm, M, ipiv);
            assert(ierr == 0);
            rbf_diff12_weights(rbf, n, x, y, ldm, M, ipiv, ldw, W);
            pack(i, n, 5, W, ldw, pack_data);
        }
    }

    if (elapsed) *elapsed = wall_time() - t0;
    return 0;
}
