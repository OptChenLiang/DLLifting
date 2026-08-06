/**
 * @file DLLifting.h
 * @brief DL / DP hybrid lifting for knapsack sets.
 *
 * Preferred public entry points:
 *   - lifting()              (C++)
 *   - dllifting_lift_cover() (C ABI)
 * Other Lifting_* symbols are internal helpers exposed for advanced use / tests.
 *
 * Version: see DLLIFTING_VERSION.
 */
#ifndef __DLLIFTING_H__
#define __DLLIFTING_H__

#define DLLIFTING_VERSION_MAJOR 1
#define DLLIFTING_VERSION_MINOR 2
#define DLLIFTING_VERSION_PATCH 1
#define DLLIFTING_VERSION "1.2.1"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#ifndef DLLIFTING_API
#  if defined(_WIN32) && defined(DLLIFTING_BUILD_SHARED)
#    ifdef DLLIFTING_EXPORTS
#      define DLLIFTING_API __declspec(dllexport)
#    else
#      define DLLIFTING_API __declspec(dllimport)
#    endif
#  else
#    define DLLIFTING_API
#  endif
#endif

#define EPS_DL 1e-6
#define INF_DL 1e+20
/** lifting(..., isdl_mode): auto DL/DP switch by threshold */
#define DLLIFTING_MODE_AUTO (-1)
/** Force DP table (isDL=0); no threshold switching */
#define DLLIFTING_MODE_DP    0
/** Force DL table (isDL=1); no threshold switching */
#define DLLIFTING_MODE_DL    1
#define MIN_DL(a,b) (a<=b? a:b) 
#define MAX_DL(a,b) (a>=b? a:b) 
#define FLOOR_DL(a) ( floor( a + EPS_DL ) )
#define CEIL_DL(a) ( ceil( a - EPS_DL ) )
#define FLOOR_INT(a) ( (long)floor( a + EPS_DL ) )
#define CEIL_INT(a) ( (long)ceil( a - EPS_DL ) )

inline double ABS_DL(double a)
{
   return fabs(a);
}

inline int ISLE(double a, double b)
{
   return a <= b + EPS_DL;
}

inline int ISLT(double a, double b)
{
   return a < b - EPS_DL;
}

inline int ISGE(double a, double b)
{
   return a + EPS_DL >= b;
}

inline int ISGT(double a, double b)
{
   return a - EPS_DL > b;
}

inline int ISZERO(double a)
{
   return ABS_DL(a) <= EPS_DL;
}

inline int ISINF(double a)
{
   return a >= INF_DL - 1;
}

inline int ISEQ(double a, double b)
{
   return ISZERO(a - b);
}

inline double Lifting_GetTime()
{
	return (double)clock()/(double)CLOCKS_PER_SEC;
}

typedef double DTptype;
typedef double DTrctype;
typedef double DTwtype;
typedef double DTutype;
typedef double DTctype;
   
typedef struct DLLifting
{
   int                isleq;          // 1: <= knapsack; 0: >= knapsack
   DTptype*           p;              // lifting coefficients 
   DTwtype*           w;              // knapsack weights 
   DTutype*           u;              // variable upper bounds 
   int*               seed;           // indices in the seed inequality
   int                n_seed;         // nembers of indices in the seed inequality
   int*               liftingorder;   // variables lifted sequence
   int                n_liftingorder; 
   int*               isuseub;        // 1: fix variable at upper bound; 0: fix variable at lower bound
   DTctype            subcap;         // residual knapsack capacity
   DTctype            cap;            // table size for DP / DL
   DTctype            maxcap;
   DTctype            minweight;      // original >= threshold (isleq = 0 only) 
   int                n;
   int                solvedsize;
   double*            x_old;
   double*            x;
   double             activity;
   double*            dplist;         // dense DP table when isDL = 0 

   DTptype*           psum;           // DL breakpoint profits 
   DTwtype*           wsum;           // DL breakpoint weights
   int                n_soltable;
   int                maxsolsize;
   double             rhs;
   double             tableleft;     // REDUCTION: U^k = sum u_i*w_i on remaining bounded vars
   int                reduction_active; // 0 if any variable is unbounded (disable REDUCTION)
   int                reduction_usable; // 1 if U^k>0 and subcap>0 after init (REDUCTION enabled)
   double             reduction_U_init; // U^k at init (after Calsubcap)
   double             reduction_b_init; // b^k = subcap at init

   DTptype*           psum1;
   DTwtype*           wsum1;
   DTptype*           psum2;
   DTwtype*           wsum2;

   int                force_mode;       // DLLIFTING_MODE_* ; >=0 disables threshold switching
   int                isDL;             // 1: DL table is authoritative; 0: dplist is */

   double             threshold;       // < 100: prefer DL; > 100: prefer DP */
   double             duration;
} Lifting;

/** Same type as @c Lifting (struct tag @c DLLifting). */
typedef Lifting DLLifting;

void Lifting_Printsum(DLLifting* lift);
void Lifting_Check(DLLifting* lift);

int Lifting_Alloc(DLLifting* lift, int len, int scale, double threshold);
int Lifting_Realloc(DLLifting* lift, int len);
int Lifting_Reset(DLLifting* lift, int len);
int Lifting_Free(DLLifting* lift);
int Lifting_Calsubcap(DLLifting* lift);
int Lifting_Calcap(DLLifting* lift);

// Initialise lifting state from knapsack data and lifting sequence
int Lifting_Init(
      DLLifting* lift, 
      DTptype* p, DTwtype* w, DTutype* u, int* isuseub, 
      DTctype cap, int issubcap, 
      int* basis, int n_basis, 
      int* liftingorder, int n_liftingorder, 
      int isleq, double* x, DTctype maxRhs, int n);

int Lifting_Wiszero(DLLifting* lift, DTptype p, DTwtype w, DTutype u);
int Lifting_Piszero(DLLifting* lift, DTptype p, DTwtype w, DTutype u);

// Merge a bounded item (p, w) into the DL table
int Lifting_Mergesort(DLLifting* lift, DTptype p, DTwtype w);

// Merge an item with effectively unbounded multiplicity
int Lifting_Mergesortinf(DLLifting* lift, DTptype p, DTwtype w);

// Add item (p, w) via binary splitting; respects force_mode / threshold
int Lifting_Multiply(DLLifting* lift, DTptype p, DTwtype w, DTutype u);

// Query DL table: last breakpoint with wsum <= cap or min profit with wsum >= cap
int Lifting_Findind(DLLifting* lift, DTctype cap, int begin, int end, int isleq);
DTptype Lifting_Findsol(DLLifting* lift, DTctype cap, int begin, int end, int isleq);

// Build seed table and return initial rhs of the cover inequality
DTptype Lifting_Calinitrhs(DLLifting* lift);

int Lifting_Iter(DLLifting* lift, DTptype p, DTwtype w, DTutype u);
int Lifting_Lifting(DLLifting* lift, DTptype* rhs);

/**
 * Run full coefficient lifting for a cover inequality
 *
 * @param cap      Knapsack rhs 
 * @param isSubCap If 1, cap is the subproblem capacity; else global rhs
 * @param seed     Variables fixed in the seed inequality
 * @param liftingorder  Remaining variables in lifting order
 * @param isLeq    1 for <= knapsack; 0 for >= knapsack
 * @param threshold  < 100 uses DL; > 100 uses DP (when isdl_mode = DLLIFTING_MODE_AUTO)
 * @param duration   Reserved; unused (timing is written to lift->duration)
 * @param isdl_mode  DLLIFTING_MODE_AUTO / _DP / _DL — force table without switching
 * @return 1 on success; writes lifted coefficients into p and rhs.
 */
DLLIFTING_API int lifting(
      DLLifting* lift,
      DTptype* p, DTwtype* w, DTutype* u, int* isuseub,
      DTctype cap, int isSubCap,
      int* seed, int n_seed,
      int* liftingorder, int n_liftingorder,
      double* rhs,
      int isLeq, double* x, int n, double threshold, double duration, int isdl_mode);

int lifting_lifting(DLLifting* lift, DTptype* alpha, DTwtype* a, DTutype* u, int* isuseub, DTptype *rhs, int n, int isleq);
void Lifting_Printsoltable(DTptype* psum, DTwtype* wsum, int n);

// Up-lifting
int Lifting_Up(DLLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs);

// Down lifting
int Lifting_Down(DLLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs);

int Lifting_Compress(DLLifting* lift, int begin);
int Lifting_Expand(DLLifting* lift);

//One unbounded knapsack DP
void Lifting_DPiterInf(DLLifting* lift, int w, double p);

// One bounded knapsack DP
void Lifting_DPiter(DLLifting* lift, int w, double p);
void Lifting_DPPrint(double* dp, int c);
void Lifting_DPFree(double* dp);

#ifdef __cplusplus
extern "C" {
#endif

/** Return codes for dllifting_lift_cover */
#define DLLIFTING_OK           0
#define DLLIFTING_ERR_ALLOC   -1
#define DLLIFTING_ERR_ARGS    -2
#define DLLIFTING_ERR_INTERNAL -3

/**
 * Lift a cover inequality for a single knapsack row (C ABI).
 *
 * Modifies @p coef in place and writes the lifted inequality rhs to @p rhs.
 *
 * @param isdl_mode  DLLIFTING_MODE_AUTO / _DP / _DL
 */
DLLIFTING_API int dllifting_lift_cover(
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
      const double* x_frac);

#ifdef __cplusplus
}
#endif

#endif
