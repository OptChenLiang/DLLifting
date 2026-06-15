#include "dllifting/dllifting_c.h"
#include <stdio.h>

int main(void)
{
   const int n = 3;
   double coef[] = {1.0, 1.0, 0.0};
   const double w[] = {2.0, 3.0, 4.0};
   const double u[] = {1.0, 1.0, 1.0};
   const int use_ub[] = {0, 0, 0};
   const int seed[] = {0, 1};
   const int order[] = {2};
   double rhs = 0.0;

   int rc = dllifting_lift_cover(
         n, coef, w, u, use_ub,
         4.0, 0,
         seed, 2, order, 1,
         &rhs, 1, 10.0, NULL);

   if(rc != DLLIFTING_OK)
   {
      fprintf(stderr, "dllifting_lift_cover failed: %d\n", rc);
      return 1;
   }

   printf("rhs = %.4f, coef = [%.4f, %.4f, %.4f]\n",
         rhs, coef[0], coef[1], coef[2]);
   return 0;
}
