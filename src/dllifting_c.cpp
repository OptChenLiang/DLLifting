#include <dllifting_c.h>
#include <DLLifting.h>

#include <cstring>

int dllifting_lift_cover(
    int n,
    double* coef,
    const double* weight,
    const double* ub,
    const int* use_ub,
    double cap,
    int is_subcap,
    const int* seed,
    int n_seed,
    const int* lifting_order,
    int n_order,
    double* rhs,
    int is_leq,
    double threshold,
    const double* x_frac)
{
   if(n <= 0 || coef == nullptr || weight == nullptr || ub == nullptr
         || use_ub == nullptr || rhs == nullptr)
      return DLLIFTING_ERR_ARGS;
   if(n_seed > 0 && seed == nullptr)
      return DLLIFTING_ERR_ARGS;
   if(n_order > 0 && lifting_order == nullptr)
      return DLLIFTING_ERR_ARGS;

   double* w = new double[n];
   double* u = new double[n];
   int* isuseub = new int[n];
   int* seed_mut = nullptr;
   int* order_mut = nullptr;
   if(w == nullptr || u == nullptr || isuseub == nullptr)
   {
      delete[] w;
      delete[] u;
      delete[] isuseub;
      return DLLIFTING_ERR_ALLOC;
   }

   for(int i = 0; i < n; i++)
   {
      w[i] = weight[i];
      u[i] = ub[i];
      isuseub[i] = use_ub[i];
   }

   if(n_seed > 0)
   {
      seed_mut = new int[n_seed];
      if(seed_mut == nullptr)
      {
         delete[] w;
         delete[] u;
         delete[] isuseub;
         return DLLIFTING_ERR_ALLOC;
      }
      std::memcpy(seed_mut, seed, (size_t)n_seed * sizeof(int));
   }

   if(n_order > 0)
   {
      order_mut = new int[n_order];
      if(order_mut == nullptr)
      {
         delete[] seed_mut;
         delete[] w;
         delete[] u;
         delete[] isuseub;
         return DLLIFTING_ERR_ALLOC;
      }
      std::memcpy(order_mut, lifting_order, (size_t)n_order * sizeof(int));
   }

   DLLifting lift;
   std::memset(&lift, 0, sizeof(lift));
   double rhs_val = 0.0;
   int ok = lifting(
         &lift,
         coef,
         w,
         u,
         isuseub,
         cap,
         is_subcap,
         seed_mut,
         n_seed,
         order_mut,
         n_order,
         &rhs_val,
         is_leq,
         const_cast<double*>(x_frac),
         n,
         threshold,
         0.0);

   delete[] order_mut;
   delete[] seed_mut;
   delete[] isuseub;
   delete[] u;
   delete[] w;

   if(!ok)
      return DLLIFTING_ERR_ALLOC;

   *rhs = rhs_val;
   return DLLIFTING_OK;
}
