/**
 * @file DLLifting.cpp
 * @brief Implementation of DL/DP hybrid knapsack lifting.
 **/

#include <DLLifting.h>
#include <cstring>

#define INITSIZE_LIFTING 5000000
#define CHECK 1
#ifndef REDUCTION
#  ifdef DLLIFTING_REDUCTION
#    define REDUCTION 1
#  endif
#endif

#ifdef DLLIFTING_DEBUG
#  define DLLIFTING_LOG(...) do { printf(__VA_ARGS__); } while (0)
#else
#  define DLLIFTING_LOG(...) ((void)0)
#endif

#ifdef DLTIME
extern double findtime;
extern double mergetime;

extern int niter;
extern double dptime;
#endif

int Lifting_Compress(DLLifting* lift, int begin);
int Lifting_Expand(DLLifting* lift);
static void Lifting_update(DLLifting* lift, int w, double p, int unbounded);

// Effective unboundedness: u is INF, or u >= 2 and a*(u+1) > b
static int Lifting_unbounded(const DLLifting* lift, DTwtype a, DTutype ub)
{
   if(ISZERO(a))
      return 0;
   if(ISINF(ub))
      return 1;
   if(ub < 1.5)
      return 0;
   return a * (ub + 1.0) > lift->maxcap;
}

void dllifting_policy_default(DLLiftingPolicy* pol)
{
   if(!pol)
      return;
   pol->rho_th = DLLIFTING_DEFAULT_RHO_TH;
   pol->beta_th = DLLIFTING_DEFAULT_BETA_TH;
   pol->u_bar_th = DLLIFTING_DEFAULT_UBAR_TH;
   pol->prefer_dl_fuzzy = 1;
}

void dllifting_compute_features(
      const DTwtype* w, const DTutype* u, int n, DTctype cap, DLLiftingFeatures* feat)
{
   if(!feat)
      return;
   feat->rho_w = 1.0;
   feat->beta = 0.0;
   feat->u_bar = 0.0;
   feat->w_mean = 0.0;
   feat->w_min = 0.0;
   feat->w_max = 0.0;
   if(!w || n <= 0)
      return;
   double wmin = w[0], wmax = w[0], wsum = 0.0, usum = 0.0;
   int i;
   for(i = 0; i < n; i++)
   {
      double wi = w[i];
      if(wi < wmin)
         wmin = wi;
      if(wi > wmax)
         wmax = wi;
      wsum += wi;
      if(u)
         usum += ISINF(u[i]) ? 1e6 : u[i];
   }
   feat->w_min = wmin;
   feat->w_max = wmax;
   feat->w_mean = wsum / (double)n;
   feat->u_bar = u ? (usum / (double)n) : 0.0;
   feat->rho_w = (wmin > EPS_DL) ? (wmax / wmin) : 1.0;
   feat->beta = (feat->w_mean > EPS_DL) ? (cap / feat->w_mean) : 0.0;
}

static void dllifting_policy_fill(DLLiftingPolicy* pol)
{
   if(pol->rho_th <= 0)
      pol->rho_th = DLLIFTING_DEFAULT_RHO_TH;
   if(pol->beta_th <= 0)
      pol->beta_th = DLLIFTING_DEFAULT_BETA_TH;
   if(pol->u_bar_th <= 0)
      pol->u_bar_th = DLLIFTING_DEFAULT_UBAR_TH;
}

int dllifting_select_backend(const DLLiftingFeatures* feat, const DLLiftingPolicy* pol)
{
   DLLiftingPolicy local;
   if(!pol)
      dllifting_policy_default(&local);
   else
   {
      local = *pol;
      dllifting_policy_fill(&local);
   }
   pol = &local;

   const double rho = feat ? feat->rho_w : 1.0;
   const double beta = feat ? feat->beta : 0.0;
   const double ubar = feat ? feat->u_bar : 0.0;

   if(ubar > 0 && ubar <= pol->u_bar_th + EPS_DL)
      return DLLIFTING_MODE_DL;
   if(rho >= pol->rho_th - EPS_DL)
      return DLLIFTING_MODE_DL;
   if(pol->prefer_dl_fuzzy)
   {
      const double band = 0.10;
      if(fabs(rho - pol->rho_th) <= band * pol->rho_th)
         return DLLIFTING_MODE_DL;
      if(fabs(beta - pol->beta_th) <= band * pol->beta_th)
         return DLLIFTING_MODE_DL;
   }
   if(beta < pol->beta_th - EPS_DL)
      return DLLIFTING_MODE_DL;
   return DLLIFTING_MODE_DP;
}

double Lifting_bar_b(const DLLifting* lift)
{
   if(!lift)
      return 0.0;
   double b = lift->subcap;
   if(b < 0)
      b = 0;
#ifdef REDUCTION
   if(lift->reduction_active && lift->tableleft > 0)
      return (b < lift->tableleft) ? b : lift->tableleft;
#endif
   return b;
}

/** τ > bar b → prefer DP (0); τ < bar b → prefer DL (1). */
static int Lifting_prefer_dl_from_barb(const DLLifting* lift)
{
   const double barb = Lifting_bar_b(lift);
   double tau = lift->switch_cap > 0 ? lift->switch_cap : lift->threshold;
   if(tau > barb + EPS_DL)
      return 0;
   if(tau + EPS_DL < barb)
      return 1;
   return lift->isDL ? 1 : 0;
}

/** MODE_THRESHOLD mid-lift: switch DL↔DP for one binary-split chunk. */
static void Lifting_add_chunk(DLLifting* lift, DTptype p, DTwtype ww, int geq_path)
{
   const int thr_mode = (lift->force_mode == DLLIFTING_MODE_THRESHOLD);
   if(lift->isDL)
   {
      if(thr_mode && !Lifting_prefer_dl_from_barb(lift))
      {
         Lifting_Expand(lift);
         Lifting_DPiter(lift, FLOOR_INT(ww), p);
         lift->isDL = 0;
      }
      else if(geq_path)
         Lifting_update(lift, FLOOR_INT(ww), p, 0);
      else
         Lifting_Mergesort(lift, p, ww);
   }
   else
   {
      if(thr_mode && Lifting_prefer_dl_from_barb(lift))
      {
         Lifting_Compress(lift, 0);
         Lifting_Mergesort(lift, p, ww);
         lift->isDL = 1;
      }
      else
         Lifting_DPiter(lift, FLOOR_INT(ww), p);
   }
}

static void Lifting_resolve_tau(DLLifting* lift, double threshold, double w_mean)
{
   if(threshold > 0)
      lift->switch_cap = threshold;
   else
      lift->switch_cap = lift->beta_th * (w_mean > EPS_DL ? w_mean : 1.0);
   lift->threshold = lift->switch_cap;
}

static void Lifting_apply_isdl_mode(DLLifting* lift, int isdl_mode,
      const DLLiftingFeatures* feat, const DLLiftingPolicy* pol, double cap)
{
   if(isdl_mode == DLLIFTING_MODE_DP)
   {
      lift->force_mode = DLLIFTING_MODE_DP;
      lift->isDL = 0;
   }
   else if(isdl_mode == DLLIFTING_MODE_DL)
   {
      lift->force_mode = DLLIFTING_MODE_DL;
      lift->isDL = 1;
   }
   else if(isdl_mode == DLLIFTING_MODE_THRESHOLD)
   {
      lift->force_mode = DLLIFTING_MODE_THRESHOLD;
      lift->isDL = (cap > lift->switch_cap + EPS_DL) ? 1 : 0;
   }
   else
   {
      lift->force_mode = DLLIFTING_MODE_AUTO;
      lift->isDL = (dllifting_select_backend(feat, pol) == DLLIFTING_MODE_DL) ? 1 : 0;
   }
}

#ifdef REDUCTION
// U^k = tableleft; b^k = subcap; bar b^k = min(subcap, U^k).  Slack m = max(0, subcap - U^k).
// When U^k > subcap the table is truncated to subcap (bar b^k = subcap); do not zero U^k.
static int Lifting_reduction(DLLifting* lift, double* m_out)
{
   if(!lift->reduction_active)
      return 0;
   if(lift->tableleft <= 0 || lift->subcap <= 0)
      return 0;
   double m = lift->subcap - lift->tableleft;
   if(m < 0)
      m = 0;
   *m_out = m;
   return 1;
}

static double Lifting_reduction_floor(const DLLifting* lift, double floor_w)
{
   if(!lift->reduction_active || lift->tableleft <= 0 || lift->subcap <= 0)
      return floor_w;
   double m = lift->subcap - lift->tableleft;
   return m > floor_w ? m : floor_w;
}

static int Lifting_reduction_start(const DLLifting* lift, int floor_j)
{
   if(!lift->reduction_active || lift->tableleft <= 0 || lift->subcap <= 0)
      return floor_j;
   double rm = lift->subcap - lift->tableleft;
   if(rm <= 0)
      return floor_j;
   int mb = CEIL_INT(rm);
   return mb > floor_j ? mb : floor_j;
}

static void Lifting_record_reduction_init(DLLifting* lift)
{
   lift->reduction_U_init = lift->tableleft;
   lift->reduction_b_init = lift->subcap;
   lift->reduction_usable = (lift->tableleft > 0 && lift->subcap > 0) ? 1 : 0;
}

static double Lifting_length(const DLLifting* lift, int i)
{
   if(Lifting_unbounded(lift, lift->w[i], lift->u[i]))
      return 0.0;
   return lift->u[i] * lift->w[i];
}
#endif


/* Min p-cost to cover residual demand `cap` on the current geq DP/DL table.
 * Residual queries must NOT be clamped to the original RHS (minweight): that
 * made Up-lifting see dp[b] for every b-ja < b and forced alpha=0. */
static DTptype Lifting_Geqfind(DLLifting* lift, DTctype cap)
{
   if(cap <= 0)
      return 0.0;
   if(lift->isDL)
      Lifting_Expand(lift);
   {
      int q = FLOOR_INT(cap);
      if(q < 0)
         return 0.0;
      if(q > (int)lift->cap)
         q = (int)lift->cap;
      return lift->dplist[q];
   }
}

/* Expand -> DPiter(Inf) -> Compress -> Expand. */
static void Lifting_update(DLLifting* lift, int w, double p, int unbounded)
{
   if(lift->isDL)
      Lifting_Expand(lift);
   if(unbounded)
      Lifting_DPiterInf(lift, w, p);
   else
      Lifting_DPiter(lift, w, p);
   if(lift->isDL)
   {
      /* Keep DL authoritative: MODE_DL, MODE_AUTO, or MODE_THRESHOLD while τ says DL. */
      if(lift->force_mode == DLLIFTING_MODE_DL ||
            lift->force_mode == DLLIFTING_MODE_AUTO ||
            (lift->force_mode == DLLIFTING_MODE_THRESHOLD && Lifting_prefer_dl_from_barb(lift)))
      {
         Lifting_Compress(lift, 0);
         Lifting_Expand(lift);
      }
      else if(lift->force_mode == DLLIFTING_MODE_THRESHOLD)
         lift->isDL = 0;
   }
}

// knapsack DP kernels
void Lifting_DPiterInf(DLLifting* lift, int w, double p)
{
   int c = FLOOR_INT(lift->cap);
   double* dp = lift->dplist;
#ifdef DLTIME
   double tmp = Lifting_GetTime();
#endif
   int j;
   if(lift->isleq)
   {
      for(j = w; j<= c; j++ )
      {
         if(dp[j] < dp[j-w] + p)
         {
            dp[j] = dp[j-w] + p;
         }
      }
   }
   else
   {
      for(j = 0; j <= w; j++ )
      {
         if(dp[j] > p) 
            dp[j] = p;
      }
      for(; j <= c; j++ )
      {
         if(dp[j] > dp[j-w] + p)
         {
            dp[j] = dp[j-w] + p;
         }
      }
   }
#ifdef DLTIME
   dptime += Lifting_GetTime() - tmp;
#endif
}

void Lifting_DPiter(DLLifting* lift, int w, double p)
{
   int c = FLOOR_INT(lift->cap);
   double* dp = lift->dplist;
#ifdef DLTIME
   double tmp = Lifting_GetTime();
#endif
   int j;
   if(lift->isleq)
   {
      double m = w;
#ifdef REDUCTION 
      m = Lifting_reduction_floor(lift, w);
#endif
      double g = dp[c];
      int n_soltable = 1;
      for(j = c; j>= m; j-- )
      {
         if(dp[j] < dp[j-w] + p)
         {
            dp[j] = dp[j-w] + p;
         }
         if(!ISEQ(dp[j], g))
         {
            g = dp[j];
            n_soltable ++;
         }
      }
      lift->n_soltable = n_soltable;
   }
   else
   {
      /* Bounded >= cover: one 0-1 item (weight w, cost p) into min-cost-to-cover-at-least-j.
       * Forward dp[j]=dp[j-w]+p would allow unbounded reuse of this copy — wrong.
       * Scratch = second half of dplist (size 2*maxsolsize from Alloc). */
      double* old = dp + lift->maxsolsize;
      int jstart = 0;
      for(j = 0; j <= c; j++)
         old[j] = dp[j];
#ifdef REDUCTION
      jstart = Lifting_reduction_start(lift, 0);
#endif
      for(j = jstart; j <= c; j++)
      {
         double cand;
         if(j <= w)
            cand = p; /* remaining demand 0 */
         else if(old[j - w] >= INF_DL / 2)
            continue;
         else
            cand = p + old[j - w];
         if(cand < dp[j])
            dp[j] = cand;
      }
   }
#ifdef DLTIME
   dptime += Lifting_GetTime() - tmp;
#endif
}

void Lifting_Print(DLLifting* lift)
{
   if(lift->isDL)
      Lifting_Printsum(lift);
   else
      Lifting_DPPrint(lift->dplist, lift->cap);
}

void Lifting_DPPrint(double* dp, int c)
{
#ifdef DLLIFTING_DEBUG
   double p = dp[0];
   int i = 0;
   printf( "w = 0, p = %.1f\n", p);
   for(i = 0; i<=c; i++)
   {
      if(!ISEQ(dp[i], p))
      {
         p = dp[i];
         printf( "w = %d, p = %.1f\n", i, p);
      }
   }
   printf("\n");
#else
   (void)dp;
   (void)c;
#endif
}

void Lifting_DPFree(double* dp)
{
   free(dp);
}

void Lifting_Printsum(DLLifting* lift)
{
#ifdef DLLIFTING_DEBUG
   int i = 0;
   printf("len: %d\n", lift->n_soltable);
   for(i = 0; i<lift->n_soltable; i++)
   {
      printf("wsum = %.2f, psum = %.2f\n", lift->wsum[i], lift->psum[i]);
   }
   printf("\n");
#else
   (void)lift;
#endif
}

void Lifting_Check(DLLifting* lift)
{
   int i = 0;
   for(i = 1; i<lift->n_soltable; i++)
   {
      if(!ISLT(lift->wsum[i-1], lift->wsum[i]) || !ISLT(lift->psum[i-1], lift->psum[i]))
         return;
   }
}


int Lifting_Alloc(DLLifting* lift, int len, int scale, double threshold)
{
   (void)scale;
   len = len+1;
   // Do not cap at INITSIZE_LIFTING: silent truncation caused buffer overruns when
   // b (or variable upper bounds) exceed 5e6 on the lambda benchmark grid.
   if(lift->dplist != nullptr || lift->psum1 != nullptr)
      Lifting_Free(lift);
   if(sizeof(DTptype) !=  sizeof(DTwtype))
   {
      lift->psum1 = (DTptype*) malloc ( 2 * len * sizeof(DTptype));
      if(lift->psum1 == nullptr)
         return 0;
      lift->wsum1 = (DTwtype*) malloc ( 2 * len * sizeof(DTwtype)); 
      if(lift->wsum1 == nullptr)
         return 0;
      lift->psum2 = lift->psum1 + len;
      lift->wsum2 = lift->wsum1 + len;
   }
   else
   {
      lift->psum1 = (DTptype*) malloc ( 4 * len * sizeof(DTptype));
      if(lift->psum1 == nullptr)
         return 0;
      lift->psum2 = lift->psum1 + len;
      lift->wsum1 = (DTwtype*)(lift->psum1 + 2*len);
      lift->wsum2 = lift->wsum1 + len;
   }
   
   lift->maxsolsize = len;
   lift->psum = lift->psum1;
   lift->wsum = lift->wsum1;

   lift->dplist = (double*) malloc ( (2* len)*sizeof(double)) ;
   if(lift->dplist == nullptr)
      return 0;
   /* Initial isDL overwritten by lifting() for AUTO/forced modes. */
   lift->isDL = 1;
   lift->threshold = threshold;
   lift->switch_cap = threshold;
   lift->maxsolsize = len; 
   return 1;
}

int Lifting_Realloc(DLLifting* lift, int len)
{
   DTptype* psum1;
   DTwtype* wsum1;
   DTptype* psum2;
   DTwtype* wsum2;
   if(sizeof(DTptype) !=  sizeof(DTwtype))
   {
      psum1 = (DTptype*) realloc (lift->psum1, 2 * len * sizeof(DTptype));
      if(psum1 == nullptr)
         return 0;
      psum2 = psum1 + len;
      wsum1 = (DTwtype*) realloc (lift->wsum1, 2 * len * sizeof(DTwtype)); 
      if(wsum1 == nullptr)
         return 0;
      wsum2 = wsum1 + len;
   }
   else
   {
      psum1 = (DTptype*) realloc (lift->psum1, 4 * len * sizeof(DTptype));
      if(psum1 == nullptr)
         return 0;
      psum2 = psum1 + len;
      wsum1 = (DTwtype*)(psum1 + 2*len);
      wsum2 = wsum1 + len;   
   }
   if(lift->psum == lift->psum1)
   {
      lift->psum = psum1;
      lift->wsum = wsum1;
   }
   else
   {
      lift->psum = psum2;
      lift->wsum = wsum2;
   }
   lift->psum1 = psum1;
   lift->psum2 = psum2;
   lift->wsum1 = wsum1;
   lift->wsum2 = wsum2;
   lift->maxsolsize = len;
   DLLIFTING_LOG("realloc len = %d\n", len);
   return 1;
}

// Reset dominated lists
int Lifting_Reset(DLLifting* lift, int len)
{
   (void)len;
   lift->n_soltable = 1;
   lift->wsum[0] = 0;
   lift->psum[0] = 0.0;


   int i = 0;
   if(lift->isleq)
   {
      for(i = 0; i<=lift->cap; i++)
      {
         lift->dplist[i] = 0;
      }
   }
   else
   {
      lift->dplist[0] = 0;
      for(i = 1; i<=lift->cap; i++)
         lift->dplist[i] = INF_DL;
   }
   return 1; 
}

int Lifting_Free(DLLifting* lift)
{
   Lifting_DPFree(lift->dplist);
   lift->dplist = nullptr;
   if(sizeof(DTptype) !=  sizeof(DTwtype))
   {
      if(lift->psum1 != nullptr)
         free(lift->psum1);
      if(lift->wsum1 != nullptr)
         free(lift->wsum1);
   }
   else
   {
      if(lift->psum1 != nullptr)
         free(lift->psum1);
   }
   lift->psum1 = lift->psum2 = nullptr;
   lift->wsum1 = lift->wsum2 = nullptr;
   lift->psum = lift->wsum = nullptr;
   return 1;
}

int Lifting_Calsubcap(DLLifting* lift)
{
   int i = 0;
   lift->subcap = lift->cap;
   for(i = 0; i<lift->n_liftingorder; i++)
   {
      if(lift->isuseub[lift->liftingorder[i]])
      {
         lift->subcap = lift->subcap - lift->u[lift->liftingorder[i]]*lift->w[lift->liftingorder[i]];
      }
#ifdef REDUCTION 
      else if(lift->reduction_active)
      {
         lift->tableleft += Lifting_length(lift, lift->liftingorder[i]); 
      }
#endif
   }
#ifdef REDUCTION 
   if(lift->reduction_active)
   {
      for(i = 0; i< lift->n_seed; i++)
      {
         lift->tableleft += Lifting_length(lift, lift->seed[i]); 
      }
   }
#endif
   return 1;
}


int Lifting_Calcap(DLLifting* lift)
{
   int i = 0;
   lift->cap = lift->subcap;
   for(i = 0; i<lift->n_liftingorder; i++)
   {
      if(lift->isuseub[lift->liftingorder[i]])
         lift->cap = lift->cap + lift->u[lift->liftingorder[i]]*lift->w[lift->liftingorder[i]];
#ifdef REDUCTION 
      else if(lift->reduction_active)
      {
         lift->tableleft += Lifting_length(lift, lift->liftingorder[i]); 
      }
#endif
   }
#ifdef REDUCTION 
   if(lift->reduction_active)
   {
      for(i = 0; i< lift->n_seed; i++)
      {
         lift->tableleft += Lifting_length(lift, lift->seed[i]); 
      }
   }
#endif
   return 1;
}

int Lifting_Wiszero(DLLifting* lift, DTptype p, DTwtype w, DTutype u)
{
   (void)w;
   if(ISZERO(p))
      return 1;
   if(lift->isleq)
   {
      if(ISINF(u))
      {
         lift->wsum[0] = 0;
         lift->psum[0] = INF_DL;
         lift->n_soltable = 1;
      }
      else
      {
         int i = 0;
         for(i = 0; i<lift->n_soltable;i++)
            lift->psum[i] = lift->psum[i] + p*u;
      }
   }
   return 1;
}

int Lifting_Piszero(DLLifting* lift, DTptype p, DTwtype w, DTutype u)
{
   if(ISZERO(w))
      return 1;
   if(lift->isleq == 0)
   {
      /* Zero-profit items still contribute weight. Always apply via the DP
       * 0-1 update (Lifting_update / DPiter): shifting DL wsum would drop the
       * zero-weight empty cover and disagree with the DP table. */
      DTutype k;
      DTutype uu = u;
      if(Lifting_unbounded(lift, w, u))
      {
         Lifting_update(lift, FLOOR_INT(w), 0.0, 1);
         return 1;
      }
      for(k = 1; uu != 0; k += k)
      {
         if(k > uu)
            k = uu;
         Lifting_update(lift, FLOOR_INT(w * k), 0.0, 0);
         uu -= k;
      }
      (void)p;
      return 1;
   }
   (void)p;
   return 1;
}

// Merge one bounded item into the DL table.
int Lifting_Mergesort(DLLifting* lift, DTptype p, DTwtype w)
{
   int i = 0, j = 0, k = 0;
   DTptype* oldpsum = nullptr;
   DTwtype* oldwsum = nullptr;
   DTptype* newpsum = nullptr;
   DTwtype* newwsum = nullptr;

   if(lift->psum == lift->psum1)
   {
      oldpsum = lift->psum1;
      oldwsum = lift->wsum1;
      newpsum = lift->psum2;
      newwsum = lift->wsum2;
   }
   else
   {
      oldpsum = lift->psum2;
      oldwsum = lift->wsum2;
      newpsum = lift->psum1;
      newwsum = lift->wsum1;
   }
   if(lift->isleq)
   {
#ifdef DLTIME
      double tmp = Lifting_GetTime();
#endif

      newwsum[k] = oldwsum[i];
      newpsum[k] = oldpsum[i];
#ifdef REDUCTION 
      double m = 0;
      int use_red = Lifting_reduction(lift, &m);
#endif
      double tmpwsum;      
      double tmppsum;      
      while( i < lift->n_soltable &&  j < lift->n_soltable && ISLE(newwsum[k], lift->cap))
      {
         tmpwsum = oldwsum[j] + w;
         tmppsum = oldpsum[j] + p;
         if(ISLE(oldwsum[i], tmpwsum))
         {
            if(ISGT(oldpsum[i], newpsum[k]))
            {
#ifdef REDUCTION 
               if(ISGT(oldwsum[i], newwsum[k]) && (!use_red || oldwsum[i] >= m - EPS_DL))
#else
                  if(ISGT(oldwsum[i], newwsum[k]))
#endif
                     k++;
               newwsum[k] = oldwsum[i];
               newpsum[k] = oldpsum[i];
            }
            i++;
         }
         else
         {
            if(ISGT(tmppsum, newpsum[k]))
            {
#ifdef REDUCTION 
               if(ISGT(tmpwsum, newwsum[k]) && (!use_red || tmpwsum >= m - EPS_DL))
#else
                  if(ISGT(tmpwsum, newwsum[k]))
#endif
                     k++;
               newwsum[k] = tmpwsum;
               newpsum[k] = tmppsum;
            }
            j++;
         }
      }
      while(j < lift->n_soltable && ISLE(newwsum[k], lift->cap))
      {
         tmpwsum = oldwsum[j] + w;
         tmppsum = oldpsum[j] + p;
         if(ISGT(tmppsum, newpsum[k]))
         {
#ifdef REDUCTION 
            if(ISGT(tmpwsum, newwsum[k]) && (!use_red || tmpwsum >= m - EPS_DL))
#else
               if(ISGT(tmpwsum, newwsum[k]))
#endif
                  k++;
            newwsum[k] = tmpwsum;
            newpsum[k] = tmppsum;
         }
         j++;
      }
      if(ISGT(newwsum[k], lift->cap))
      {
         k--;
      }
#ifdef DLTIME
      mergetime += Lifting_GetTime() - tmp;
#endif
   }
   else
   {
#ifdef DLTIME
      double tmp = Lifting_GetTime();
#endif
      newwsum[k] = oldwsum[i];
      newpsum[k] = oldpsum[i];
#ifdef REDUCTION 
      double m = 0;
      int use_red = Lifting_reduction(lift, &m);
#endif
      double tmpwsum;
      double tmppsum;
      while( i < lift->n_soltable &&  j < lift->n_soltable && ISLE(newwsum[k], lift->cap))
      {
         tmpwsum = oldwsum[j] + w;
         tmppsum = oldpsum[j] + p;
         if(ISLE(oldwsum[i], tmpwsum))
         {
            if(ISLT(oldpsum[i], newpsum[k]))
            {
#ifdef REDUCTION 
               if(ISGT(oldwsum[i], newwsum[k]) && (!use_red || oldwsum[i] >= m - EPS_DL))
#else
               if(ISGT(oldwsum[i], newwsum[k]))
#endif
                  k++;
               newwsum[k] = oldwsum[i];
               newpsum[k] = oldpsum[i];
            }
            i++;
         }
         else
         {
            if(ISGT(tmpwsum, newwsum[k]) || ISLT(tmppsum, newpsum[k]))
            {
#ifdef REDUCTION 
               if(ISGT(tmpwsum, newwsum[k]) && (!use_red || tmpwsum >= m - EPS_DL))
#else
               if(ISGT(tmpwsum, newwsum[k]))
#endif
                  k++;
               newwsum[k] = tmpwsum;
               newpsum[k] = tmppsum;
            }
            j++;
         }
      }
      while(j < lift->n_soltable && ISLE(newwsum[k], lift->cap))
      {
         tmpwsum = oldwsum[j] + w;
         tmppsum = oldpsum[j] + p;
         if(ISGT(tmpwsum, newwsum[k]) || ISLT(tmppsum, newpsum[k]))
         {
#ifdef REDUCTION 
            if(ISGT(tmpwsum, newwsum[k]) && (!use_red || tmpwsum >= m - EPS_DL))
#else
            if(ISGT(tmpwsum, newwsum[k]))
#endif
               k++;
            newwsum[k] = tmpwsum;
            newpsum[k] = tmppsum;
         }
         j++;
      }
      if(ISGT(newwsum[k], lift->cap))
         k--;
#ifdef DLTIME
      mergetime += Lifting_GetTime() - tmp;
#endif
   }

   lift->psum = newpsum;
   lift->wsum = newwsum;
   lift->n_soltable = k + 1;
   return 1;
}

int Lifting_Mergesortinf(DLLifting* lift, DTptype p, DTwtype w)
{
   int i = 0, j = 0, k = 0;
   DTptype* oldpsum = nullptr;
   DTwtype* oldwsum = nullptr;
   DTptype* newpsum = nullptr;
   DTwtype* newwsum = nullptr;

   int newsize = CEIL_INT(lift->maxcap)+1;
   if(newsize > lift->maxsolsize)
   {
      if( Lifting_Realloc(lift, newsize) != 1)
         return 0;
      lift->maxsolsize = newsize;
   }
   if(lift->psum == lift->psum1)
   {
      oldpsum = lift->psum1;
      oldwsum = lift->wsum1;
      newpsum = lift->psum2;
      newwsum = lift->wsum2;
   }
   else
   {
      oldpsum = lift->psum2;
      oldwsum = lift->wsum2;
      newpsum = lift->psum1;
      newwsum = lift->wsum1;
   }
   if( lift->isleq )
   {
      newwsum[k] = oldwsum[i];
      newpsum[k] = oldpsum[i];
      k++;
      while(i< lift->n_soltable)
      {
         if (ISGE(oldwsum[i], newwsum[k - 1]) && ISLE(oldpsum[i], newpsum[k - 1])) 
         {
            i++;
         }
         else if (ISGE(newwsum[j] + w, newwsum[k - 1]) && ISLE(newpsum[j] + p, newpsum[k - 1])) 
         {
            j++;
         }
         else if(ISLE(oldwsum[i], newwsum[j] + w) && ISGE(oldpsum[i], newpsum[j] + p )) 
         {
            newwsum[k] = oldwsum[i];
            newpsum[k] = oldpsum[i];
            i++;
            j++;
            k++;
         }
         else if(ISLE(newwsum[j] + w, oldwsum[i]) && ISGE(newpsum[j] + p, oldpsum[i]))  
         {
            newwsum[k] = newwsum[j] + w;
            newpsum[k] = newpsum[j] + p;
            i++;
            j++;
            k++;
         }
         else if(ISGE(newwsum[j] + w, oldwsum[i]) && ISGE(newpsum[j] + p, oldpsum[i]))   
         {
            newwsum[k] = oldwsum[i];
            newpsum[k] = oldpsum[i];
            i++;
            k++;
         }
         else if(ISLE(newwsum[j] + w, oldwsum[i]) && ISLE(newpsum[j] + p, oldpsum[i]))   
         {
            newwsum[k] = newwsum[j]+w;
            newpsum[k] = newpsum[j]+p;
            j++;
            k++;
         }
         else
         {
            break;
         }
         if (ISGE(newwsum[k - 1], lift->maxcap))
         {
            break;
         }
      }
      while(j < k && ISLT(newwsum[k-1], lift->maxcap))
      {
         if(ISGT(newwsum[j] + w, newwsum[k-1]) && ISGT(newpsum[j] + p, newpsum[k-1]))
         {
            newwsum[k] = newwsum[j] + w;
            newpsum[k] = newpsum[j] + p;
            j++;
            k++;
         }
         else
            j++;
      }
   }
   else
   {
      while( i< lift->n_soltable && ISLT(oldwsum[i], w) && ISLT(oldpsum[i], p))
      {
         newwsum[k] = oldwsum[i];
         newpsum[k] = oldpsum[i];
         k = k + 1;
         i = i + 1;
      }
      while(j< k)
      {
         newwsum[k] = newwsum[j] + w;
         newpsum[k] = newpsum[j] + p;
         if(i< lift->n_soltable && ISLT( oldwsum[i], newwsum[k]) )
         {
            newwsum[k] = oldwsum[i];
            newpsum[k] = oldpsum[i];
            i = i + 1; 
         }
         while(j< k)
         {
            if(i < lift->n_soltable && ISGE(oldwsum[i], newwsum[k]) && ISLE(oldpsum[i], newpsum[k]))
            {
               newwsum[k] = oldwsum[i];
               newpsum[k] = oldpsum[i];
               i++;
            }
            else if(i < lift->n_soltable && ISLE(oldwsum[i], newwsum[k]) && ISGE(oldpsum[i], newpsum[k]))
               i++;
            else if(ISGE(newwsum[j] + w, newwsum[k]) && ISLE(newpsum[j] + p, newpsum[k]))
            {
               newwsum[k] = newwsum[j] + w;
               newpsum[k] = newpsum[j] + p;
               j++;
            }
            else if(ISLE(newwsum[j] + w, newwsum[k]) && ISGE(newpsum[j] + p, newpsum[k]))
               j++;
            else
               break;
         }
         k++;
         if(ISGE(newwsum[k-1], lift->maxcap))
            break;
      }
   }
   lift->psum = newpsum;
   lift->wsum = newwsum;
   lift->n_soltable = k;
   return 1;
}

// Add item (p,w,u) with binary splitting; large-u items use unbounded DP.
int Lifting_Multiply(DLLifting* lift, DTptype p, DTwtype w, DTutype u)
{
   DTutype k; 
   if (lift->n_soltable == 0) 
      return 0;

   if ( ISZERO(p) ) 
   {
      Lifting_Piszero(lift, p, w, u);
      return 1;
   }
   if ( ISZERO(w) ) 
   {
      Lifting_Wiszero(lift, p, w, u);
      return 1;
   }

   if( lift->isleq )
   {
      if( Lifting_unbounded(lift, w, u) )
      {
         /* DP path maintains dplist; refresh psum before unbounded DL merge. */
         if(!lift->isDL)
            Lifting_Compress(lift, 0);
         Lifting_Mergesortinf(lift, p, w);
         if(!lift->isDL)
            Lifting_Expand(lift);
      }
      else
      {
         for(k = 1; u != 0; k += k) 
         {
            if (k > u) 
               k = u;
            Lifting_add_chunk(lift, p * k, w * k, 0);
            u -= k;
         }
      }
   }
   else
   {
      if( Lifting_unbounded(lift, w, u) )
      {
         Lifting_update(lift, FLOOR_INT(w), p, 1);
      }
      else
      {
         for (k = 1; u != 0; k += k) 
         {
            if (k >  u) 
               k = u;
            Lifting_add_chunk(lift, p * k, w * k, 1);
            u -= k;
         }
      }
   }
   lift->solvedsize++;
   return 1;
}

int Lifting_Findind(DLLifting* lift, DTctype cap, int begin, int end, int isleq)
{
   int i = begin, j = end, m = i;
   if( isleq )
   {
      if ( ISGT( lift->wsum[i], cap) ) 
               return -1;
            if ( ISLE( lift->wsum[j], cap) ) 
               return j;
      while( i < j-1 )
      {
         m = ( i + j )/2;
         if( ISLT( cap, lift->wsum[m] ) )
            j = m;
         else 
            i = m;
      }
      return i;
   }
   else
   {
      int best = -1;
      int t;
      for( t = begin; t <= end; t++ )
      {
         if( ISGE( lift->wsum[t], cap) )
         {
            if( best < 0 || ISLT( lift->psum[t], lift->psum[best] ) )
               best = t;
         }
      }
      return best;
   }
   return -1;
}

// Find the solution 
DTptype Lifting_Findsol(DLLifting* lift, DTctype cap, int begin, int end, int isleq)
{
   int ind = Lifting_Findind(lift, cap, begin, end, isleq);
   if( ind == -1)
      return INF_DL;
   return lift->psum[ind];
}

// Build seed table from seed items and return initial cover rhs.
DTptype Lifting_Calinitrhs(DLLifting* lift)
{
   int i = 0;

   for(i = 0; i< lift->n_seed; i++)
   {
      Lifting_Multiply(lift, lift->p[lift->seed[i]], lift->w[lift->seed[i]], lift->u[lift->seed[i]]); 
#ifdef REDUCTION 
      if(lift->reduction_active)
         lift->tableleft -= Lifting_length(lift, lift->seed[i]);
#endif
   }

   if(lift->isleq)
   {
      DTctype qcap = lift->subcap;
      if(lift->isDL)
      {
         i = Lifting_Findind(lift, qcap, 0, lift->n_soltable-1, 1);
         if( i == -1)
            return INF_DL;
         return lift->psum[i];
      }
      return lift->dplist[FLOOR_INT(qcap)];
   }
   return Lifting_Geqfind(lift, lift->subcap);
}

int Lifting_Init(
      Lifting* lift, 
      DTptype* p, DTwtype* w, DTutype* u, int* isuseub, 
      DTctype cap, int issubcap, 
      int * seed, int n_seed, 
      int* liftingorder, int n_liftingorder, 
      int isleq, double* x, DTctype maxcap, int n)
{
   lift->p = p;
   lift->w = w;
   lift->u = u;
   lift->isuseub = isuseub;
   lift->seed = seed;
   lift->liftingorder = liftingorder;
   lift->n_seed = n_seed;
   lift->n_liftingorder = n_liftingorder;
   lift->isleq = isleq;
   lift->x = x;
   lift->activity = 0;
   lift->maxcap = maxcap;
   lift->n = n;
#ifdef REDUCTION 
   lift->tableleft = 0;
   /* Fill U^k whenever +R might be used; finalize after Calsubcap. */
   {
      int any_unb = 0;
      int i;
      for(i = 0; i < n; i++)
      {
         if(Lifting_unbounded(lift, w[i], u[i]))
         {
            any_unb = 1;
            break;
         }
      }
      if(lift->reduction_request == DLLIFTING_RED_OFF || any_unb)
         lift->reduction_active = 0;
      else
         /* ON or AUTO: accumulate tableleft in Calsubcap; AUTO refined below. */
         lift->reduction_active = 1;
   }
#else
   lift->reduction_active = 0;
#endif

   if(issubcap == 1) 
   {
      lift->subcap = cap;
      lift->minweight = cap;
      Lifting_Calcap(lift);
   }
   else
   {
      lift->cap = cap;
      lift->minweight = cap;
      Lifting_Calsubcap(lift);
   }
#ifdef REDUCTION
   /* +R AUTO: same mid-lift rule — large residual (bar b > τ) enables reduction. */
   if(lift->reduction_request == DLLIFTING_RED_AUTO && lift->reduction_active)
   {
      double tau = lift->switch_cap;
      if(tau <= 0)
         tau = lift->threshold;
      double barb = Lifting_bar_b(lift);
      lift->reduction_active = (barb > tau + EPS_DL) ? 1 : 0;
   }
   Lifting_record_reduction_init(lift);
#else
   lift->reduction_usable = 0;
#endif
   {
      int i, v;
      for(i = 0; i < lift->n_seed; i++)
      {
         v = lift->seed[i];
         if(ISINF(lift->u[v]))
            continue;
         if(lift->w[v] * lift->u[v] > lift->cap)
            lift->cap = lift->w[v] * lift->u[v];
      }
      for(i = 0; i < lift->n_liftingorder; i++)
      {
         v = lift->liftingorder[i];
         if(ISINF(lift->u[v]))
            continue;
         if(lift->w[v] * lift->u[v] > lift->cap)
            lift->cap = lift->w[v] * lift->u[v];
      }
   }
   if(lift->cap > lift->maxcap)
      lift->maxcap = lift->cap;
   Lifting_Reset(lift, FLOOR_INT(lift->cap));
   return 1;
}

int Lifting_Iter(DLLifting* lift, DTptype p, DTwtype w, DTutype u)
{
   Lifting_Multiply(lift, p, w, u); 
   return 1;
}

// Up lifting 
int Lifting_Up(DLLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs)
{
   int solind = -1, j = 1, u0;
   double temp;
   if(lift->isleq)
   {
      *alpha = INF_DL; 
      solind = lift->n_soltable - 1;
      u0 = FLOOR_INT( MIN_DL(u, lift->subcap*1.0/a) ); 

      if(u0 == 0)
         *alpha = 0;
      else
      {
         for( j = 1; j <= u0; j++)
         {
            if(lift->isDL)
            {
#ifdef DLTIME
               double tmp = Lifting_GetTime();
#endif
               solind = Lifting_Findind(lift, lift->subcap - j*a, 0, solind, lift->isleq);
#ifdef DLTIME
               findtime += Lifting_GetTime() - tmp;
#endif
               if(solind < 0)
                  return 0;
               temp = (*rhs - lift->psum[solind])/j;

            }
            else
            {
#ifdef REDUCTION
               double rm;
               if(Lifting_reduction(lift, &rm) && lift->subcap - j * a < rm - EPS_DL)
                  continue;
#endif
               temp = (*rhs - lift->dplist[FLOOR_INT(lift->subcap - j*a)])/j;
            }
            if( temp < *alpha)
            {
               *alpha = temp;
            }
         }
      }

      if( ISZERO(*alpha) )
         *alpha = 0;
      if(!(*alpha >= 0 && *alpha < INF_DL/10))
         return 0;
      Lifting_Iter(lift, *alpha, a, u);
   }
   else
   {
      /* >= sequential up-lifting: alpha = max_j (rhs - z(b-ja))/j.
       * Always multiply into the table (even alpha=0): zero-coeff items still
       * contribute weight and must update z(·) for later variables. */
      *alpha = 0; 
      u0 = CEIL_INT( MIN_DL( u, lift->subcap*1.0/a) ); 
      if(u0 <= 0)
         return 0;
      for( j = 1; j<=u0; j++)
      {
         if(lift->subcap - j*a < 0)
            temp = (*rhs)/j;
         else
            temp = (*rhs - Lifting_Geqfind(lift, lift->subcap - j*a))/j;
         if(temp >= INF_DL / 10 || temp < -EPS_DL)
            continue;
         if( temp > *alpha)
            *alpha = temp;
      }
      if(ISZERO(*alpha))
         *alpha = 0;
      if(!(*alpha >= 0 && *alpha < INF_DL/10))
         return 0;
      Lifting_Iter(lift, *alpha, a, u); 
   }
   return 1;
}

// Down lifting
int Lifting_Down(DLLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs)
{
   int solind = -1, j = 1, u0;
   double temp;
   if(lift->isleq)
   {
      *alpha = -INF_DL; 
      u0 = CEIL_INT( MIN_DL( u, (lift->cap - lift->subcap)*1.0/a) ); 
      solind = 0;
      if(u0 <= 0)
         return 0;

      for( j = 1; j<=u0; j++)
      {
         if(lift->isDL)
         {
#ifdef DLTIME
            double tmp = Lifting_GetTime();
#endif
            solind = Lifting_Findind(lift, lift->subcap + j*a, solind, lift->n_soltable-1, lift->isleq);
#ifdef DLTIME
            findtime += Lifting_GetTime() - tmp;
#endif
            if(solind < 0)
               return 0;
            temp = (lift->psum[solind] - *rhs)/j;
         }
         else
         {
#ifdef REDUCTION
            double rm;
            if(Lifting_reduction(lift, &rm))
            {
               int capj = FLOOR_INT(lift->subcap + j * a);
               if(capj < Lifting_reduction_start(lift, 0))
                  continue;
            }
#endif
            temp = (lift->dplist[FLOOR_INT(lift->subcap + j*a)] - *rhs)/j;
         }
         if( temp > *alpha)
            *alpha = temp;
      }
      lift->subcap = lift->subcap + a*u;
      *rhs = *rhs + *alpha*u;

      if(!(*alpha >= 0 && *alpha < INF_DL/10))
         return 0;
      Lifting_Multiply(lift, *alpha, a, u); 
   }
   else
   {
      *alpha = INF_DL; 
      u0 = FLOOR_INT( MIN_DL( u, (lift->cap - lift->subcap)*1.0/a) ); 
      solind = 0;
      if(u0 == 0)
      {
         *alpha = *rhs;
      }
      else
      {
         for( j = 1; j<=u0; j++)
         {
            double dpv = Lifting_Geqfind(lift, lift->subcap + j*a);
            if(dpv >= INF_DL/2)
               temp = 0;  
            else
               temp = (dpv - *rhs)/j;
            if( temp < *alpha)
               *alpha = temp;
         }
      }
      lift->subcap = lift->subcap + a*u;
      *rhs = *rhs + *alpha*u;
      if(ISZERO(*alpha))
         *alpha = 0;
      if(!(*alpha >= 0 && *alpha < INF_DL/10))
         return 0;
      /* Always update table: alpha=0 items still change residual weight. */
      Lifting_Multiply(lift, *alpha, a, u); 
   }
   return 1;
}

// DP to DL
int Lifting_Compress(DLLifting* lift, int begin)
{
   int i = begin;
   double p = lift->dplist[i];
   int k = 1;
   lift->wsum[0] = begin;
   lift->psum[0] = p;
   for(; i<=lift->cap; i++)
   {
      if(!ISEQ(lift->dplist[i], p))
      {
         p = lift->dplist[i];
         lift->wsum[k] = i;
         lift->psum[k] = p;
         k++;
      }
   }
   lift->n_soltable = k;
   return 1;
}

// DL to DP
int Lifting_Expand(DLLifting* lift)
{
   int j;
   int k;
   if(lift->isleq == 0)
   {
      for(j = 0; j <= (int)lift->cap; j++)
         lift->dplist[j] = INF_DL;
      if(lift->n_soltable == 1 && ISZERO(lift->wsum[0]) && ISZERO(lift->psum[0]))
      {
         lift->dplist[0] = 0.0;
         return 1;
      }
   }
   else
   {
      j = (int)lift->wsum[0];
      k = 0;
      for(; k < lift->n_soltable - 1; k++)
      {
         for(; j < lift->wsum[k+1]; j++)
         {
            lift->dplist[j] = lift->psum[k];
         }
      }
      for(; j <= lift->cap; j++)
      {
         lift->dplist[j] = lift->psum[k];
      }
      return 1;
   }
   j = (int)lift->wsum[0];
   k = 0;
   for(; k < lift->n_soltable - 1; k++)
   {
      for(; j < lift->wsum[k+1]; j++)
      {
         lift->dplist[j] = lift->psum[k];
      }
   }
   for(; j <= lift->cap; j++)
   {
      lift->dplist[j] = lift->psum[k];
   }
   return 1;
}

// Sequential up- / down-lifting along liftingorder.
int Lifting_Lifting(DLLifting* lift, DTptype* rhs)
{
   int i = 0;
   for( i = 0; i<lift->n_liftingorder; i++)
   {
      if(lift->time_limit > 0.0
            && (Lifting_GetTime() - lift->t_start) > lift->time_limit)
         return 0;

      if(lift->isuseub[lift->liftingorder[i]] == 0)  
      {
         if(!Lifting_Up(lift, &lift->p[lift->liftingorder[i]],
               lift->w[lift->liftingorder[i]], lift->u[lift->liftingorder[i]], rhs))
            return 0;
#ifdef REDUCTION 
         if(lift->reduction_active)
            lift->tableleft -= Lifting_length(lift, lift->liftingorder[i]);
#endif
         if(!ISGE(lift->p[lift->liftingorder[i]], 0))
            return 0;
      }
      else
      {
         if(!Lifting_Down(lift, &lift->p[lift->liftingorder[i]],
               lift->w[lift->liftingorder[i]], lift->u[lift->liftingorder[i]], rhs))
            return 0;
         if(!ISGE(lift->p[lift->liftingorder[i]], 0))
            return 0;
      }
   }
   return 1;
}

// Public driver: allocate tables, initialise, lift seed + order, free. 
int lifting(
      DLLifting* lift, 
      DTptype* p, DTwtype* w, DTutype* u, int* isuseub, 
      DTctype cap, int issubcap, 
      int* seed, int n_seed, 
      int* liftingorder, int n_liftingorder, 
      double* rhs,
      int isleq, double* x, int n, double threshold, double duration, int isdl_mode)
{
   if(lift == nullptr)
      return 0;
   /* Preserve optional policy overrides across memset of workspace. */
   const double save_rho_th = lift->rho_th;
   const double save_beta_th = lift->beta_th;
   const double save_ubar_th = lift->u_bar_th;
   const int save_red_req = lift->reduction_request;

   std::memset(lift, 0, sizeof(*lift));
   lift->rho_th = save_rho_th;
   lift->beta_th = save_beta_th;
   lift->u_bar_th = save_ubar_th;
   lift->reduction_request = save_red_req;
   lift->time_limit = duration;
   lift->t_start = Lifting_GetTime();

   {
      DTctype allocCap = cap;
      int i, v;
      for(i = 0; i < n_seed; i++)
      {
         v = seed[i];
         DTctype bound = ISINF(u[v]) ? cap : w[v] * u[v];
         if(bound > allocCap)
            allocCap = bound;
      }
      for(i = 0; i < n_liftingorder; i++)
      {
         v = liftingorder[i];
         DTctype bound = ISINF(u[v]) ? cap : w[v] * u[v];
         if(bound > allocCap)
            allocCap = bound;
      }
      if(Lifting_Alloc(lift, FLOOR_INT(allocCap), 1, threshold) == 0)
         return 0;
   }

   {
      DLLiftingPolicy pol;
      dllifting_policy_default(&pol);
      if(lift->rho_th > 0)
         pol.rho_th = lift->rho_th;
      else
         lift->rho_th = pol.rho_th;
      if(lift->beta_th > 0)
         pol.beta_th = lift->beta_th;
      else
         lift->beta_th = pol.beta_th;
      if(lift->u_bar_th > 0)
         pol.u_bar_th = lift->u_bar_th;
      else
         lift->u_bar_th = pol.u_bar_th;

      DLLiftingFeatures feat;
      dllifting_compute_features(w, u, n, cap, &feat);
      lift->feat_rho = feat.rho_w;
      lift->feat_beta = feat.beta;
      lift->feat_ubar = feat.u_bar;

      Lifting_resolve_tau(lift, threshold, feat.w_mean);
      Lifting_apply_isdl_mode(lift, isdl_mode, &feat, &pol, cap);
   }

   Lifting_Init(lift, p, w, u, isuseub, cap, issubcap, seed, n_seed, liftingorder, n_liftingorder, isleq, x, cap, n);

   if(isdl_mode == DLLIFTING_MODE_THRESHOLD)
      lift->isDL = Lifting_prefer_dl_from_barb(lift) ? 1 : 0;

   lift->duration = 0;

   clock_t startTime = clock();
   lift->rhs = Lifting_Calinitrhs(lift);
   /* Seed set must cover subcap; otherwise init cut is undefined. */
   if(lift->rhs >= INF_DL / 2 || !(lift->rhs >= 0))
   {
      Lifting_Free(lift);
      return 0;
   }
   if(!Lifting_Lifting(lift, &lift->rhs))
   {
      Lifting_Free(lift);
      return 0;
   }
   clock_t endTime = clock();

   lift->duration  = (double) (endTime - startTime) / CLOCKS_PER_SEC;
   *rhs = lift->rhs;

   Lifting_Free(lift);

   return 1;
}

/* C ABI wrapper (formerly src/dllifting_c.cpp). */
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
      int isdl_mode,
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
         0.0,
         isdl_mode);

   delete[] order_mut;
   delete[] seed_mut;
   delete[] isuseub;
   delete[] u;
   delete[] w;

   if(!ok)
      return DLLIFTING_ERR_INTERNAL;

   *rhs = rhs_val;
   return DLLIFTING_OK;
}
