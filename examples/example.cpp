/**
 * DPLifting example — 5-variable bounded knapsack (same as README).
 *
 *   C = {x1,x2}, N0 = {x3,x4}, Nu = {x5}, lift order {x3,x4,x5}
 *   seed: 2x1 + x2 <= 4 on residual capacity b^2 = 18
 */
#include <DPLifting.h>
#include <cstdio>

static void print_cut(int n, const double* coef, double rhs) {
   std::printf("  ");
   for (int i = 0; i < n; i++) {
      if (coef[i] > EPS_DPL)
         std::printf("%.4f*x_%d + ", coef[i], i + 1);
   }
   std::printf("<= %.4f\n", rhs);
}

int main(void) {
   const int n = 5;
   double p[] = {2.0, 1.0, 0.0, 0.0, 0.0};
   double w[] = {8.0, 5.0, 4.0, 3.0, 5.0};
   double u[] = {2.0, 3.0, 6.0, 5.0, 1.0};
   int isuseub[] = {0, 0, 0, 0, 1};   /* Nu = {x5}: down-lift */
   int seed[] = {0, 1};               /* C = {x1, x2} */
   int order[] = {2, 3, 4};           /* lift x3, x4, x5 */
   double rhs = 0.0;

   DPLifting lift = {};
   if (!lifting(&lift, p, w, u, isuseub, 18.0, 1,
            seed, 2, order, 3, &rhs, 1, nullptr, n,
            /* threshold */ 0.0, /* duration */ 0.0, DPLIFTING_MODE_DPT)) {
      std::fprintf(stderr, "lifting failed\n");
      return 1;
   }

   std::printf("lifting() MODE_DPT (README 5-var example):\n");
   print_cut(n, p, rhs);
   std::printf("  (%.4f s)\n", lift.duration);

   /* Same instance via C ABI */
   double coef[] = {2.0, 1.0, 0.0, 0.0, 0.0};
   rhs = 0.0;
   if (dplifting_lift_cover(n, coef, w, u, isuseub, 18.0, 1,
            seed, 2, order, 3, &rhs, 1, 0.0, DPLIFTING_MODE_DPT, nullptr)
         != DPLIFTING_OK) {
      std::fprintf(stderr, "dplifting_lift_cover failed\n");
      return 1;
   }
   std::printf("dplifting_lift_cover:\n");
   print_cut(n, coef, rhs);
   return 0;
}
