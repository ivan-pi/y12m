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
