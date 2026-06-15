/**
 * @file DLLifting.h
 * @brief DL / DP hybrid lifting for knapsack sets.
 *
 * Implements coefficient lifting for inequalities of the form
 *   sum_i p_i x_i <= rhs
 * subject to a knapsack constraint sum_i w_i x_i <= b (isleq = 1) or
 * sum_i w_i x_i >= b (isleq = 0).  
 * Subproblems are solved either by DL (threshold < 100) or by a DP (threshold > 100). 
 */
#ifndef __DLLIFTING_H__
#define __DLLIFTING_H__

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

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
   double             tableleft;     // REDUCTION: L=sum u_i*w_i; m=subcap-L
   int                reduction_active; // 0 if any variable is unbounded (disable REDUCTION)

   DTptype*           psum1;
   DTwtype*           wsum1;
   DTptype*           psum2;
   DTwtype*           wsum2;

   int                onlyDL;
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

// Add item (p, w) via binary splitting; switches DL/DP by threshold
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
 * @param threshold  < 100 uses DL; > 100 uses DP
 * @return 1 on success; writes lifted coefficients into p and rhs.
 */
DLLIFTING_API int lifting(
      DLLifting* lift,
      DTptype* p, DTwtype* w, DTutype* u, int* isuseub,
      DTctype cap, int isSubCap,
      int* seed, int n_seed,
      int* liftingorder, int n_liftingorder,
      double* rhs,
      int isLeq, double* x, int n, double threshold, double duration);

int lifting_lifting(DLLifting* lift, DTptype* alpha, DTwtype* a, DTutype* u, int* isuseub, DTptype *rhs, int n, int isleq);
void Lifting_Printsoltable(DTptype* psum, DTwtype* wsum, int n);

// Up-lifting
int Lifting_Up(DLLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs);

// Down lifting
int Lifting_Down(DLLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs);

int Lifting_Compress(DLLifting* lift, int begin = 0);
int Lifting_Expand(DLLifting* lift);

//One unbounded knapsack DP
void Lifting_DPiterInf(DLLifting* lift, int w, double p);

// One bounded knapsack DP
void Lifting_DPiter(DLLifting* lift, int w, double p);
void Lifting_DPPrint(double* dp, int c);
void Lifting_DPFree(double* dp);

#endif
