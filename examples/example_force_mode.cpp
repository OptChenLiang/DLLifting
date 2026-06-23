/**
 * Compare forced DL vs forced DP on the same cover instance (no threshold switching).
 */
#include <DLLifting.h>
#include <cstdio>
#include <cstring>

static int run_once(double* p, const double* w, const double* u,
      const int* isuseub, int isdl_mode, double* out_rhs, double* out_time)
{
   const int n = 3;
   int seed[] = {0, 1};
   int order[] = {2};
   double rhs = 0.0;
   DLLifting lift;
   std::memset(&lift, 0, sizeof(lift));
   int ok = lifting(&lift, p, const_cast<double*>(w), const_cast<double*>(u),
         const_cast<int*>(isuseub), 4.0, 0, seed, 2, order, 1, &rhs, 1,
         nullptr, n, 200.0, 0.0, isdl_mode);
   if(!ok)
      return 0;
   *out_rhs = rhs;
   *out_time = lift.duration;
   return 1;
}

int main()
{
   double p_dl[] = {1.0, 1.0, 0.0};
   double p_dp[] = {1.0, 1.0, 0.0};
   const double w[] = {2.0, 3.0, 4.0};
   const double u[] = {1.0, 1.0, 1.0};
   const int isuseub[] = {0, 0, 0};
   double rhs_dl = 0.0, rhs_dp = 0.0;
   double t_dl = 0.0, t_dp = 0.0;

   if(!run_once(p_dl, w, u, isuseub, DLLIFTING_MODE_DL, &rhs_dl, &t_dl)
         || !run_once(p_dp, w, u, isuseub, DLLIFTING_MODE_DP, &rhs_dp, &t_dp))
   {
      fprintf(stderr, "lifting failed\n");
      return 1;
   }

   printf("forced DL: rhs=%.4f  time=%.6f s  coef=[%.4f, %.4f, %.4f]\n",
         rhs_dl, t_dl, p_dl[0], p_dl[1], p_dl[2]);
   printf("forced DP: rhs=%.4f  time=%.6f s  coef=[%.4f, %.4f, %.4f]\n",
         rhs_dp, t_dp, p_dp[0], p_dp[1], p_dp[2]);

   if(rhs_dl != rhs_dp)
      fprintf(stderr, "warning: rhs differs (expected on some instances)\n");

   return 0;
}
