#include <stdio.h>
#include "rbf_fd.h"

int main(void)
{
    rbf_poly_t rbf = {.q = 3, .p = 2};

    // 3x3 uniform stencil in local coordinates (center at origin).
    // n=9, p=2 gives m=6 polynomial terms and nt=15.
    const int n = 9;
    const double h = 1.0;
    double x[] = {-h, 0.0,  h, -h, 0.0,  h, -h, 0.0,  h};
    double y[] = {-h,  -h, -h, 0.0, 0.0, 0.0,  h,   h,  h};

    // const int does not create a compile-time constant in C99 -- these are
    // VLAs.  Use #define or enum for compile-time-sized arrays; use malloc
    // for large stencils where stack space matters.
    const int nt = rbf_system_size(rbf, n);
    double M[nt * nt];
    int    ipiv[nt];

    int ierr = rbf_factorize(rbf, n, x, y, nt, M, ipiv);
    if (ierr != 0) {
        fprintf(stderr, "rbf_factorize failed: %d\n", ierr);
        return 1;
    }

    // --- Gradient weights: dx, dy ---
    rbf_deriv_t grad[2] = {{1, 0}, {0, 1}};
    double W[nt * 2];

    ierr = rbf_stencil_weights(rbf, n, x, y, nt, M, ipiv, 2, grad, nt, W);
    if (ierr != 0) { fprintf(stderr, "rbf_stencil_weights failed\n"); return 1; }

    printf("dx weights:");
    for (int i = 0; i < n; i++) printf(" %7.4f", W[i]);
    printf("\ndy weights:");
    for (int i = 0; i < n; i++) printf(" %7.4f", W[nt + i]);
    printf("\n\n");

    // --- Laplacian via rbf_operator_weights (coeffs = NULL => all 1.0) ---
    rbf_deriv_t lap_ders[2] = {{2, 0}, {0, 2}};
    double wlap[nt];

    ierr = rbf_operator_weights(rbf, n, x, y, nt, M, ipiv, 2, lap_ders, NULL, wlap);
    if (ierr != 0) { fprintf(stderr, "rbf_operator_weights failed\n"); return 1; }

    printf("Laplacian weights:");
    for (int i = 0; i < n; i++) printf(" %7.4f", wlap[i]);
    printf("\n\n");

    // --- Anisotropic operator: 2*dxx + 3*dyy ---
    double aniso_c[2] = {2.0, 3.0};
    double waniso[nt];

    ierr = rbf_operator_weights(rbf, n, x, y, nt, M, ipiv, 2, lap_ders, aniso_c, waniso);
    if (ierr != 0) { fprintf(stderr, "rbf_operator_weights (aniso) failed\n"); return 1; }

    printf("2*dxx + 3*dyy weights:");
    for (int i = 0; i < n; i++) printf(" %7.4f", waniso[i]);
    printf("\n");

    return 0;
}
