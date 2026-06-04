#include <stdio.h>
#include "rbf_fd.h"

int main(void)
{
    // 3×3 uniform stencil in local coordinates (center at origin).
    // n=9, p=2 gives m=6 polynomial terms and nt=15.
    rbfpoly_t rbf = {.q = 3, .p = 2};
    const int n  = 9;
    const double h = 1.0;

    double x[] = {-h, 0.0, h, -h, 0.0, h, -h, 0.0,  h};
    double y[] = {-h,  -h, -h, 0.0, 0.0, 0.0,  h,   h,  h};

    const int m  = rbf_poly_terms(rbf.p);  // 6
    const int nt = n + m;                  // 15

    double M[15 * 15];
    int    ipiv[15];

    int ierr = factorized_rbf_matrix(&rbf, n, x, y, nt, M, ipiv);
    if (ierr != 0) {
        fprintf(stderr, "factorized_rbf_matrix failed: %d\n", ierr);
        return 1;
    }

    // --- Gradient weights (dx, dy) via compute_derivative_weights ---
    derivative_t grad[2] = {{1, 0}, {0, 1}};
    double W[15 * 2];

    ierr = compute_derivative_weights(&rbf, n, x, y, nt, M, ipiv, 2, grad, nt, W);
    if (ierr != 0) { fprintf(stderr, "compute_derivative_weights failed\n"); return 1; }

    printf("dx weights (first %d stencil pts):", n);
    for (int i = 0; i < n; i++) printf(" %7.4f", W[i]);
    printf("\ndy weights (first %d stencil pts):", n);
    for (int i = 0; i < n; i++) printf(" %7.4f", W[nt + i]);
    printf("\n\n");

    // --- Laplacian via compute_operator_weights (coeffs = NULL → all 1.0) ---
    derivative_t lap_ders[2] = {{2, 0}, {0, 2}};
    double wlap[15];

    ierr = compute_operator_weights(&rbf, n, x, y, nt, M, ipiv, 2, lap_ders, NULL, wlap);
    if (ierr != 0) { fprintf(stderr, "compute_operator_weights failed\n"); return 1; }

    printf("Laplacian weights:");
    for (int i = 0; i < n; i++) printf(" %7.4f", wlap[i]);
    printf("\n\n");

    // --- Anisotropic operator: 2*dxx + 3*dyy with explicit coefficients ---
    derivative_t aniso_ders[2] = {{2, 0}, {0, 2}};
    double       aniso_c[2]    = {2.0, 3.0};
    double       waniso[15];

    ierr = compute_operator_weights(&rbf, n, x, y, nt, M, ipiv,
                                    2, aniso_ders, aniso_c, waniso);
    if (ierr != 0) { fprintf(stderr, "compute_operator_weights (aniso) failed\n"); return 1; }

    printf("2*dxx + 3*dyy weights:");
    for (int i = 0; i < n; i++) printf(" %7.4f", waniso[i]);
    printf("\n");

    return 0;
}
