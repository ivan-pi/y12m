/* SPDX-License-Identifier: GPL-2.0-only */
#include "y12m.h"

#include <math.h>
#include <stdio.h>

int main(void)
{
   {
      const int n = 3;
      const int nz = 7;
      const float a[] = {2.0f, -1.0f, -1.0f, 2.0f, -1.0f, -1.0f, 2.0f};
      const int snr[] = {1, 1, 2, 2, 2, 3, 3};
      float work[3] = {0.0f, 0.0f, 0.0f};
      float anorm = y12mhe(n, nz, a, snr, work);
      if (fabsf(anorm - 4.0f) > 1.0e-6f) {
         fprintf(stderr, "y12mhe returned %g (expected 4)\n", anorm);
         return 1;
      }
   }

   {
      const int n = 3;
      const int nz = 7;
      const double a[] = {2.0, -1.0, -1.0, 2.0, -1.0, -1.0, 2.0};
      const int snr[] = {1, 1, 2, 2, 2, 3, 3};
      double work[3] = {0.0, 0.0, 0.0};
      double anorm = y12mhf(n, nz, a, snr, work);
      if (fabs(anorm - 4.0) > 1.0e-12) {
         fprintf(stderr, "y12mhf returned %.17g (expected 4)\n", anorm);
         return 1;
      }
   }

   {
      const int n = 5;
      const int z = 12;
      const int nn = 50;
      const int nn1 = 50;
      const int iha = 5;
      float a[50] = {0.0f};
      int snr[50] = {0};
      int rnr[50] = {0};
      float b[5];
      float pivot[5];
      int ha[5 * 11];
      float aflag[8];
      int iflag[10];
      const float a0[12] = {
         2.0f, 3.0f, 3.0f, 4.0f, 6.0f, -1.0f, -3.0f, 2.0f, 1.0f, 4.0f, 2.0f, 1.0f
      };
      const int rnr0[12] = {1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5};
      const int snr0[12] = {1, 2, 1, 3, 5, 2, 3, 4, 3, 2, 3, 5};
      const float rhs[5] = {8.0f, 45.0f, -3.0f, 3.0f, 19.0f};
      float resid[5];
      float max_abs_resid;
      int ifail;
      int i;

      for (i = 0; i < z; ++i) {
         a[i] = a0[i];
         snr[i] = snr0[i];
         rnr[i] = rnr0[i];
      }
      for (i = 0; i < n; ++i) {
         b[i] = rhs[i];
         resid[i] = rhs[i];
      }

      ifail = y12mae(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b);
      if (ifail != 0) {
         fprintf(stderr, "y12mae returned ifail=%d\n", ifail);
         return 1;
      }

      for (i = 0; i < z; ++i) {
         const int row = rnr0[i] - 1;
         const int col = snr0[i] - 1;
         resid[row] = resid[row] - a0[i] * b[col];
      }
      max_abs_resid = 0.0f;
      for (i = 0; i < n; ++i) {
         const float abs_resid = fabsf(resid[i]);
         if (abs_resid > max_abs_resid) {
            max_abs_resid = abs_resid;
         }
      }
      if (max_abs_resid > 1.0e-4f) {
         fprintf(stderr, "y12mae max residual=%g\n", max_abs_resid);
         return 1;
      }
   }

   {
      const int n = 1;
      const int nn = 1;
      const int iha = 1;
      const float a[] = {1.0f};
      const int snr[] = {1};
      float w[] = {0.0f};
      const float pivot[] = {1.0f};
      const float anorm = 1.0f;
      const int ha[] = {0, 0, 0};
      const int iflag[] = {0, 0, 0, 0, 0};
      int ifail = 1;
      float rcond = y12mge(n, nn, a, snr, w, pivot, anorm, ha, iha, iflag, &ifail);
      if (fabsf(rcond + 1.0f) > 1.0e-6f) {
         fprintf(stderr, "y12mge returned rcond=%g (expected -1)\n", rcond);
         return 1;
      }
   }

   {
      const int n = 1;
      const int nn = 1;
      const int iha = 1;
      const double a[] = {1.0};
      const int snr[] = {1};
      double w[] = {0.0};
      const double pivot[] = {1.0};
      const double anorm = 1.0;
      const int ha[] = {0, 0, 0};
      const int iflag[] = {0, 0, 0, 0, 0};
      int ifail = 1;
      double rcond = y12mgf(n, nn, a, snr, w, pivot, anorm, ha, iha, iflag, &ifail);
      if (fabs(rcond + 1.0) > 1.0e-12) {
         fprintf(stderr, "y12mgf returned rcond=%.17g (expected -1)\n", rcond);
         return 1;
      }
   }

   return 0;
}
