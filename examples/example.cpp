/**
 * DLLifting usage examples:
 *   1. lifting()           — C++ API with forced or auto DL/DP mode
 *   2. dllifting_lift_cover — stable C ABI
 */
#include <DLLifting.h>
#include <cstdio>
#include <cstring>

static void print_cut(int n, const double* coef, double rhs) {
   std::printf("  ");
   for (int i = 0; i < n; i++) {
      if (coef[i] > EPS_DL)
         std::printf("%.4f*x_%d + ", coef[i], i + 1);
   }
   std::printf("<= %.4f\n", rhs);
}

static int demo_lifting_auto(void) {
   const int n = 3;
   double p[] = {1.0, 1.0, 0.0};
   double w[] = {2.0, 3.0, 4.0};
   double u[] = {1.0, 1.0, 1.0};
   int isuseub[] = {0, 0, 0};
   int seed[] = {0, 1};
   int order[] = {2};
   double rhs = 0.0;

   DLLifting lift = {};
   if (!lifting(&lift, p, w, u, isuseub, 4.0, 0, seed, 2, order, 1,
            &rhs, 1, nullptr, n, 10.0, 0.0, DLLIFTING_MODE_AUTO)) {
      std::fprintf(stderr, "lifting (AUTO) failed\n");
      return 1;
   }

   std::printf("lifting() with DLLIFTING_MODE_AUTO:\n");
   print_cut(n, p, rhs);
   std::printf("  (%.4f s)\n", lift.duration);
   return 0;
}

static int demo_force_mode(void) {
   const int n = 3;
   double p_dl[] = {1.0, 1.0, 0.0};
   double p_dp[] = {1.0, 1.0, 0.0};
   double w[] = {2.0, 3.0, 4.0};
   double u[] = {1.0, 1.0, 1.0};
   int isuseub[] = {0, 0, 0};
   int seed[] = {0, 1};
   int order[] = {2};
   double rhs_dl = 0.0, rhs_dp = 0.0;
   DLLifting lift = {};

   std::memset(&lift, 0, sizeof(lift));
   if (!lifting(&lift, p_dl, w, u, isuseub, 4.0, 0, seed, 2, order, 1,
            &rhs_dl, 1, nullptr, n, 200.0, 0.0, DLLIFTING_MODE_DL))
      return 1;

   std::memset(&lift, 0, sizeof(lift));
   if (!lifting(&lift, p_dp, w, u, isuseub, 4.0, 0, seed, 2, order, 1,
            &rhs_dp, 1, nullptr, n, 200.0, 0.0, DLLIFTING_MODE_DP))
      return 1;

   std::printf("forced DL vs DP (same instance):\n");
   std::printf("  DL: rhs=%.4f  coef=[%.4f, %.4f, %.4f]\n",
         rhs_dl, p_dl[0], p_dl[1], p_dl[2]);
   std::printf("  DP: rhs=%.4f  coef=[%.4f, %.4f, %.4f]\n",
         rhs_dp, p_dp[0], p_dp[1], p_dp[2]);
   return 0;
}

static int demo_c_api(void) {
   const int n = 3;
   double coef[] = {1.0, 1.0, 0.0};
   const double w[] = {2.0, 3.0, 4.0};
   const double u[] = {1.0, 1.0, 1.0};
   const int use_ub[] = {0, 0, 0};
   const int seed[] = {0, 1};
   const int order[] = {2};
   double rhs = 0.0;

   if (dllifting_lift_cover(n, coef, w, u, use_ub, 4.0, 0,
            seed, 2, order, 1, &rhs, 1, 10.0, DLLIFTING_MODE_AUTO, nullptr)
         != DLLIFTING_OK) {
      std::fprintf(stderr, "dllifting_lift_cover failed\n");
      return 1;
   }

   std::printf("dllifting_lift_cover (C ABI):\n");
   print_cut(n, coef, rhs);
   return 0;
}

int main(void) {
   if (demo_lifting_auto() != 0)
      return 1;
   if (demo_force_mode() != 0)
      return 1;
   if (demo_c_api() != 0)
      return 1;
   return 0;
}
