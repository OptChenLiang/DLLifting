/**
 * Minimal C++ example: lift a <= cover inequality and print the cut.
 */
#include "DLLifting.h"
#include <cstdio>

int main()
{
   const int n = 3;
   double p[] = {1.0, 1.0, 0.0};
   double w[] = {2.0, 3.0, 4.0};
   double u[] = {1.0, 1.0, 1.0};
   int isuseub[] = {0, 0, 0};
   int seed[] = {0, 1};
   int order[] = {2};
   double rhs = 0.0;
   const double cap = 4.0;

   DLLifting lift = {};
   if(!lifting(&lift, p, w, u, isuseub, cap, 0, seed, 2, order, 1,
            &rhs, 1, nullptr, n, 10.0, 0.0))
   {
      fprintf(stderr, "lifting failed\n");
      return 1;
   }

   printf("Lifted cut: ");
   for(int i = 0; i < n; i++)
      printf("%.4f*x_%d + ", p[i], i);
   printf("<= %.4f\n", rhs);
   return 0;
}
