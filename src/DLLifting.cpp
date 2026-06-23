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

#ifdef DLTIME
extern double findtime;
extern double mergetime;

extern int niter;
extern double dptime;
#endif

int Lifting_Compress(DLLifting* lift, int begin);
int Lifting_Expand(DLLifting* lift);

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

#ifdef REDUCTION
// L = sum u_i*w_i; m = subcap - L.  Invalid when any variable is unbounded.
static int Lifting_reduction(const DLLifting* lift, double* m_out)
{
   if(!lift->reduction_active)
      return 0;
   if(lift->tableleft <= 0 || lift->tableleft > lift->subcap)
      return 0;
   *m_out = lift->subcap - lift->tableleft;
   return 1;
}

static double Lifting_length(const DLLifting* lift, int i)
{
   if(Lifting_unbounded(lift, lift->w[i], lift->u[i]))
      return 0.0;
   return lift->u[i] * lift->w[i];
}
#endif


static DTctype Lifting_Geqcap(const DLLifting* lift, DTctype cap)
{
   if(cap < lift->minweight)
      return lift->minweight;
   return cap;
}

static DTptype Lifting_Geqfind(DLLifting* lift, DTctype cap)
{
   DTctype qcap = Lifting_Geqcap(lift, cap);
   if(lift->isDL)
      Lifting_Expand(lift);
   return lift->dplist[FLOOR_INT(qcap)];
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
      if(lift->force_mode == DLLIFTING_MODE_DL || lift->threshold <= 100)
      {
         Lifting_Compress(lift, 0);
         Lifting_Expand(lift);
      }
      else if(lift->force_mode < 0)
         lift->isDL = false;
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
      {
         double rm;
         if(Lifting_reduction(lift, &rm) && rm > m)
            m = rm;
      }
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
      int jstart;
      for(j = 0; j <= w; j++ )
      {
         if(dp[j] > p) 
            dp[j] = p;
      }
#ifdef REDUCTION
      jstart = w + 1;
      {
         double rm;
         if(Lifting_reduction(lift, &rm))
         {
            int mb = CEIL_INT(rm);
            if(jstart < mb)
               jstart = mb;
         }
      }
      for(j = jstart; j <= c; j++ )
#else
      for(j = w + 1; j <= c; j++ )
#endif
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

void Lifting_Print(DLLifting* lift)
{
#ifdef REDUCTION 
   double m = lift->subcap - lift->tableleft;
#endif
   if(lift->isDL)
      Lifting_Printsum(lift);
   else
      Lifting_DPPrint(lift->dplist, lift->cap);
}

void Lifting_DPPrint(double* dp, int c)
{
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
}

void Lifting_DPFree(double* dp)
{
   free(dp);
}

void Lifting_Printsum(DLLifting* lift)
{
   int i = 0;
   printf("len: %d\n", lift->n_soltable);
   for(i = 0; i<lift->n_soltable; i++)
   {
      printf("wsum = %.2f, psum = %.2f\n", lift->wsum[i], lift->psum[i]);
   }
   printf("\n");
}

void Lifting_Check(DLLifting* lift)
{
   int i = 0;
   for(i = 1; i<lift->n_soltable; i++)
   {
      assert(ISLT(lift->wsum[i-1], lift->wsum[i]));
      assert(ISLT(lift->psum[i-1], lift->psum[i]));
   }
}


int Lifting_Alloc(DLLifting* lift, int len, int scale, double threshold)
{
   len = len+1;
   if(len > INITSIZE_LIFTING)
      len = INITSIZE_LIFTING;
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
   if(lift->dplist == NULL)
   {
      fprintf(stderr, "ERROR: malloc dplist\n");
      return 0;
   }
   lift->isDL = (threshold <= 100.0) ? 1 : 0;
   lift->threshold = threshold;
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
   printf("relloc len = %d\n", len);
   return 1;
}

// Reset dominated lists
int Lifting_Reset(DLLifting* lift, int len)
{
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
   int* ind = lift->liftingorder;
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
      if(lift->tableleft > lift->subcap)
         lift->tableleft = 0;
   }
#endif
   return 1;
}


int Lifting_Calcap(DLLifting* lift)
{
   int i = 0;
   int* ind = lift->liftingorder;
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
      if(lift->tableleft > lift->subcap)
         lift->tableleft = 0;
   }
#endif
   return 1;
}

int Lifting_Wiszero(DLLifting* lift, DTptype p, DTwtype w, DTutype u)
{
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
      if(ISINF(u))
      {
         lift->wsum[0] = INF_DL;
         lift->psum[0] = 0;
         lift->n_soltable = 1;
      }
      else
      {
         int i = 0;
         for(i = 0; i<lift->n_soltable;i++)
            lift->wsum[i] = lift->wsum[i] + w*u;
      }
   }
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

   int newsize = 2*lift->n_soltable+1;
   if(0 && newsize > lift->maxsolsize && newsize <= lift->maxcap)
   {
      newsize = 2*newsize;
      if( Lifting_Realloc(lift, newsize) != 1)
      {
         printf("ERROR: Lifting_Realloc\n ");
         exit(0);
      }
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
      {
         printf("ERROR: Lifting_Realloc\n ");
         exit(0);
      }
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
            if(lift->isDL)
            {
               if(lift->force_mode < 0 && lift->threshold > 100)
               {
                  Lifting_Expand(lift);
                  Lifting_DPiter(lift, FLOOR_INT(w*k), p*k);
                  lift->isDL = false;
               }
               else
                  Lifting_Mergesort(lift, p*k, w*k);
            }
            else
            {
               if(lift->force_mode < 0 && lift->threshold < 100)
               {
                  Lifting_Compress(lift);
                  Lifting_Mergesort(lift, p*k, w*k);
                  lift->isDL = true;
               }
               else
               {
                  Lifting_DPiter(lift, FLOOR_INT(w*k), p*k);
               }
            }
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

            if(lift->isDL)
               Lifting_update(lift, FLOOR_INT(w*k), p*k, 0);
            else
            {
               if(lift->force_mode < 0 && lift->threshold < 100)
               {
                  Lifting_Compress(lift);
                  Lifting_Mergesort(lift, p*k, w*k);
                  lift->isDL = true;
               }
               else
               {
                  Lifting_DPiter(lift, FLOOR_INT(w*k), p*k);
               }
            }
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
   int* ind = lift->seed;

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
         {
            printf("Find %.0f\n", qcap);
            assert(i >= 0);
            return INF_DL;
         }
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
   lift->reduction_active = 1;
   {
      int i;
      for(i = 0; i < n; i++)
      {
         if(Lifting_unbounded(lift, w[i], u[i]))
         {
            lift->reduction_active = 0;
            break;
         }
      }
   }
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
   bool iszero = 1;
   double w0;
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
               assert(solind >= 0);
               temp = (*rhs - lift->psum[solind])/j;

            }
            else
            {
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
      assert( *alpha >= 0 && *alpha < INF_DL/10);
      Lifting_Iter(lift, *alpha, a, u);
   }
   else
   {
      *alpha = 0; 
      u0 = CEIL_INT( MIN_DL( u, lift->subcap*1.0/a) ); 
      assert(u0 > 0);
      for( j = 1; j<=u0; j++)
      {
         temp = (*rhs - Lifting_Geqfind(lift, lift->subcap - j*a))/j;
         if( temp > *alpha)
            *alpha = temp;
      }
      assert(*alpha >= 0 && *alpha < INF_DL/10);
      if(!ISZERO(*alpha))
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
      assert(u0 > 0);

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
            assert(solind >= 0);
            temp = (lift->psum[solind] - *rhs)/j;
         }
         else
         {
            temp = (lift->dplist[FLOOR_INT(lift->subcap + j*a)] - *rhs)/j;
         }
         if( temp > *alpha)
            *alpha = temp;
      }
      lift->subcap = lift->subcap + a*u;
      *rhs = *rhs + *alpha*u;

      assert(*alpha >= 0 && *alpha < INF_DL/10);
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
      assert(*alpha >= 0 && *alpha < INF_DL/10);
      if(!ISZERO(*alpha))
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
   int has_nonzero_lifting = 0;
   for( i = 0; i<lift->n_liftingorder; i++)
   {

      int var_idx = lift->liftingorder[i];

       DTptype old_p = lift->p[var_idx]; 

      if(lift->isuseub[lift->liftingorder[i]] == 0)  
      {
         Lifting_Up(lift, &lift->p[lift->liftingorder[i]], lift->w[lift->liftingorder[i]], lift->u[lift->liftingorder[i]], rhs);
#ifdef REDUCTION 
         if(lift->reduction_active)
            lift->tableleft -= Lifting_length(lift, lift->liftingorder[i]);
#endif
         assert(ISGE(lift->p[lift->liftingorder[i]],0));
      }
      else
      {
         Lifting_Down(lift, &lift->p[lift->liftingorder[i]], lift->w[lift->liftingorder[i]], lift->u[lift->liftingorder[i]], rhs);
         assert(ISGE(lift->p[lift->liftingorder[i]],0));
      }
      DTptype new_p = lift->p[var_idx];
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
   /* Caller may pass an uninitialized struct; clear before Alloc/Free logic. */
   std::memset(lift, 0, sizeof(*lift));

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

   if(isdl_mode == DLLIFTING_MODE_DP)
   {
      lift->isDL = 0;
      lift->force_mode = DLLIFTING_MODE_DP;
   }
   else if(isdl_mode == DLLIFTING_MODE_DL)
   {
      lift->isDL = 1;
      lift->force_mode = DLLIFTING_MODE_DL;
   }
   else
   {
      lift->force_mode = DLLIFTING_MODE_AUTO;
   }

   Lifting_Init(lift, p, w, u, isuseub, cap, issubcap, seed, n_seed, liftingorder, n_liftingorder, isleq, x, cap, n); 


   lift->duration = 0;

   clock_t startTime = clock();
   lift->rhs = Lifting_Calinitrhs(lift);
   Lifting_Lifting(lift, &lift->rhs);
   clock_t endTime = clock();

   lift->duration  = (double) (endTime - startTime) / CLOCKS_PER_SEC;
   *rhs = lift->rhs;

   Lifting_Free(lift);

   return 1;
}
