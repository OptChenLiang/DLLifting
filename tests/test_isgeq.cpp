/**
 * Tests for DLLifting with isleq = 0 (knapsack constraint sum w_i x_i >= cap).
 * Compares DP / dominated-list building and full lifting against brute-force references.
 */
#include <DLLifting.h>
#include <stdio.h>
#include <string.h>

static int g_fail = 0;

static void fail(const char* msg)
{
   fprintf(stderr, "FAIL: %s\n", msg);
   g_fail++;
}

static void ok(const char* msg)
{
   printf("OK: %s\n", msg);
}

/* Reference DP matching Lifting_Reset + repeated Lifting_DPiter (!isleq). */
static void ref_dp_geq(int cap, double* dp, int nw, const double* pw, const double* pp)
{
   int i, j, k;
   for (i = 0; i <= cap; i++)
      dp[i] = (i == 0 ? 0.0 : INF_DL);
   for (k = 0; k < nw; k++) {
      int w = (int)pw[k];
      double p = pp[k];
      if (w <= 0)
         continue;
      for (j = 0; j <= w && j <= cap; j++) {
         if (dp[j] > p)
            dp[j] = p;
      }
      for (j = w + 1; j <= cap; j++) {
         if (dp[j] > dp[j - w] + p)
            dp[j] = dp[j - w] + p;
      }
   }
}

static int dp_tables_equal(int cap, const double* a, const double* b)
{
   int j;
   for (j = 0; j <= cap; j++) {
      if (!ISEQ(a[j], b[j]))
         return 0;
   }
   return 1;
}

/* For >= knapsack: minimum p with total weight at least target (from full DP). */
static double ref_min_p_at_least(int cap, const double* dp, double target)
{
   int j;
   double best = INF_DL;
   int t = (int)ceil(target - EPS_DL);
   if (t < 0)
      t = 0;
   if (t > cap)
      return INF_DL;
   for (j = t; j <= cap; j++) {
      if (dp[j] < best)
         best = dp[j];
   }
   return best;
}

static int setup_lift(DLLifting* lift, int cap, int isleq, double threshold)
{
   memset(lift, 0, sizeof(*lift));
   lift->isleq = isleq;
   lift->cap = cap;
   lift->subcap = cap;
   lift->maxcap = cap;
   if (!Lifting_Alloc(lift, cap, 1, threshold))
      return 0;
   Lifting_Reset(lift, cap);
   return 1;
}

static void teardown_lift(DLLifting* lift)
{
   Lifting_Free(lift);
}

/* ----- Lifting_Findind: isleq = 0 finds largest index with wsum[i] >= cap ----- */
static void test_findind_geq()
{
   DLLifting lift;
   memset(&lift, 0, sizeof(lift));
   double wsum[] = {0, 3, 7, 12, 20};
   double psum[] = {0, 1, 2, 4, 5};
   lift.wsum = wsum;
   lift.psum = psum;
   lift.n_soltable = 5;

   int i;

   i = Lifting_Findind(&lift, 0, 0, 4, 0);
   if (i != 0)
      fail("Findind geq: cap=0 should return 0");
   else
      ok("Findind geq cap=0");

   i = Lifting_Findind(&lift, 7, 0, 4, 0);
   if (i != 2)
      fail("Findind geq: cap=7 should return index of weight 7");
   else
      ok("Findind geq cap=7");

   i = Lifting_Findind(&lift, 15, 0, 4, 0);
   if (i != 4)
      fail("Findind geq: cap=15 should return index of weight 20");
   else
      ok("Findind geq cap=15");

   i = Lifting_Findind(&lift, 25, 0, 4, 0);
   if (i != -1)
      fail("Findind geq: cap=25 infeasible should return -1");
   else
      ok("Findind geq cap=25 infeasible");

   i = Lifting_Findind(&lift, 5, 0, 4, 1);
   if (i != 1)
      fail("Findind leq: cap=5 should return index of weight 3");
   else
      ok("Findind leq cap=5 (sanity)");
}

/* ----- DP path (!isleq): one and multiple items ----- */
static void test_dp_geq_single_item()
{
   DLLifting lift;
   double ref[64];
   int cap = 20;

   if (!setup_lift(&lift, cap, 0, 50.0))
   {
      fail("alloc in test_dp_geq_single_item");
      return;
   }

   lift.isDL = 0;
   Lifting_DPiter(&lift, 5, 3.0);

   ref_dp_geq(cap, ref, 1, (const double[]){5}, (const double[]){3.0});
   if (!dp_tables_equal(cap, lift.dplist, ref))
      fail("DP geq single item w=5 p=3");
   else
      ok("DP geq single item w=5 p=3");

   teardown_lift(&lift);
}

static void test_dp_geq_multi_item()
{
   DLLifting lift;
   double ref[128];
   int cap = 30;
   const double ws[] = {4, 7, 11};
   const double ps[] = {2, 5, 8};
   int k;

   if (!setup_lift(&lift, cap, 0, 50.0))
   {
      fail("alloc in test_dp_geq_multi_item");
      return;
   }

   lift.isDL = 0;
   for (k = 0; k < 3; k++)
      Lifting_DPiter(&lift, (int)ws[k], ps[k]);

   ref_dp_geq(cap, ref, 3, ws, ps);
   if (!dp_tables_equal(cap, lift.dplist, ref))
      fail("DP geq three items");
   else
      ok("DP geq three items");

   teardown_lift(&lift);
}

/* DL table invariants for >= (increasing weights, increasing minimum profits). */
static void test_dl_table_invariants_geq()
{
   DLLifting lift;
   int cap = 25;
   const double ws[] = {3, 5, 8};
   const double ps[] = {1, 4, 6};
   int k, i;

   if (!setup_lift(&lift, cap, 0, 10.0)) {
      fail("alloc test_dl_table_invariants_geq");
      return;
   }

   lift.isDL = 1;
   for (k = 0; k < 3; k++)
      Lifting_Multiply(&lift, ps[k], ws[k], 1);

   for (i = 1; i < lift.n_soltable; i++) {
      if (!ISGT(lift.wsum[i], lift.wsum[i - 1])) {
         printf("  non-increasing wsum[%d]=%.2f wsum[%d]=%.2f\n",
               i - 1, lift.wsum[i - 1], i, lift.wsum[i]);
         fail("DL geq: wsum not strictly increasing");
         teardown_lift(&lift);
         return;
      }
      if (!ISGT(lift.psum[i], lift.psum[i - 1])) {
         printf("  non-increasing psum[%d]=%.2f psum[%d]=%.2f\n",
               i - 1, lift.psum[i - 1], i, lift.psum[i]);
         fail("DL geq: psum not strictly increasing");
         teardown_lift(&lift);
         return;
      }
   }
   for (i = 0; i < lift.n_soltable; i++) {
      if (lift.wsum[i] > lift.cap + EPS_DL) {
         fail("DL geq: weight exceeds cap");
         teardown_lift(&lift);
         return;
      }
   }
   ok("DL geq table invariants (monotone wsum/psum)");

   teardown_lift(&lift);
}

/* End-to-end: DL path (threshold=10) vs DP path (threshold=200) must agree. */
static void run_lifting_compare_paths(const char* name, int n,
      double* p_dl, double* w, double* u, int* isuseub,
      double cap, int* seed, int n_seed, int* order, int n_ord)
{
   double p_dp[32];
   double rhs_dl, rhs_dp;
   int i;

   if (n > 32) {
      fail("run_lifting_compare_paths: n too large");
      return;
   }

   memcpy(p_dp, p_dl, (size_t)n * sizeof(double));

   DLLifting lift_dl, lift_dp;
   memset(&lift_dl, 0, sizeof(lift_dl));
   memset(&lift_dp, 0, sizeof(lift_dp));

   if (!lifting(&lift_dl, p_dl, w, u, isuseub, cap, 0, seed, n_seed,
            order, n_ord, &rhs_dl, 0, NULL, n, 10.0, 0.0, DLLIFTING_MODE_AUTO)) {
      fail("DL path lifting failed");
      return;
   }

   if (!lifting(&lift_dp, p_dp, w, u, isuseub, cap, 0, seed, n_seed,
            order, n_ord, &rhs_dp, 0, NULL, n, 200.0, 0.0, DLLIFTING_MODE_AUTO)) {
      fail("DP path lifting failed");
      return;
   }

   if (!ISEQ(rhs_dl, rhs_dp)) {
      printf("  %s: rhs_dl=%.6f rhs_dp=%.6f\n", name, rhs_dl, rhs_dp);
      fail("DL vs DP path rhs mismatch (geq)");
      return;
   }
   for (i = 0; i < n; i++) {
      if (!ISEQ(p_dl[i], p_dp[i])) {
         printf("  %s: p[%d] dl=%.6f dp=%.6f\n", name, i, p_dl[i], p_dp[i]);
         fail("DL vs DP path coefficients mismatch (geq)");
         return;
      }
   }
   printf("OK: %s DL/DP paths agree (rhs=%.4f)\n", name, rhs_dl);
}

/* Enumerate all bounded assignments; check lifted inequality on >= knapsack. */
static int all_feasible_geq(int n, const double* w, const double* u, double cap,
      int* x, int idx, double wsum)
{
   if (idx == n)
      return wsum + EPS_DL >= cap;
   int xi;
   for (xi = 0; xi <= (int)u[idx]; xi++) {
      x[idx] = xi;
      if (all_feasible_geq(n, w, u, cap, x, idx + 1, wsum + w[idx] * xi))
         return 1;
   }
   return 0;
}

static int lifted_cut_valid(int isleq, int n, const double* p, const double* w,
      const double* u, double cap, double rhs)
{
   int x[32];
   int xi[32];
   int i;
   memset(xi, 0, sizeof(xi));

   if (n > 32)
      return 1;

   while (1) {
      for (i = 0; i < n; i++)
         x[i] = xi[i];
      double wsum = 0;
      for (i = 0; i < n; i++)
         wsum += w[i] * x[i];
      int feasible = isleq
            ? (wsum <= cap + EPS_DL)
            : (wsum + EPS_DL >= cap);
      if (feasible) {
         double lhs = 0;
         for (i = 0; i < n; i++)
            lhs += p[i] * x[i];
         int ok_cut = (lhs <= rhs + EPS_DL);
         if (!ok_cut) {
            printf("  violated at x=(");
            for (i = 0; i < n; i++)
               printf("%d%s", x[i], i + 1 < n ? "," : "");
            printf(") wsum=%.2f lhs=%.2f rhs=%.2f (%s)\n",
                  wsum, lhs, rhs, isleq ? "<=" : ">=");
            return 0;
         }
      }
      i = 0;
      while (i < n && xi[i] == (int)u[i])
         xi[i++] = 0;
      if (i == n)
         break;
      xi[i]++;
   }
   return 1;
}

/* Brute-force maximum valid up-lifting coefficient (one step, isleq=0). */
static double brute_up_alpha_geq(double rhs, double a, double u,
      int ntab, const double* wtab, const double* ptab, double subcap)
{
   double best = 0;
   int j, u0 = (int)u;
   if (u0 > 20)
      u0 = 20;
   for (j = 1; j <= u0; j++) {
      double need = subcap - j * a;
      double base = rhs;
      int t;
      if (need > wtab[ntab - 1] + EPS_DL)
         continue;
      if (need <= 0 + EPS_DL) {
         /* empty remainder */
      } else {
         int found = 0;
         double bestp = INF_DL;
         for (t = 0; t < ntab; t++) {
            if (wtab[t] + EPS_DL >= need && ptab[t] < bestp)
               bestp = ptab[t];
         }
         if (bestp >= INF_DL / 2)
            continue;
         base = rhs - bestp;
      }
      double alpha = base / j;
      if (alpha > best)
         best = alpha;
   }
   return best;
}

static void run_lifting_case_ex(const char* name, int isleq, double threshold,
      int n, double* p, double* w, double* u, int* isuseub,
      double cap, int* seed, int n_seed, int* liftingorder, int n_liftingorder)
{
   DLLifting lift;
   double rhs;
   int i;

   memset(&lift, 0, sizeof(lift));
   rhs = 0;

   if (!lifting(&lift, p, w, u, isuseub, cap, 0, seed, n_seed,
            liftingorder, n_liftingorder, &rhs, isleq, NULL, n, threshold, 0.0,
            DLLIFTING_MODE_AUTO)) {
      fail(name);
      printf("  lifting() returned 0\n");
      return;
   }

   if (!lifted_cut_valid(isleq, n, p, w, u, cap, rhs)) {
      char buf[256];
      snprintf(buf, sizeof(buf), "%s: lifted inequality invalid", name);
      fail(buf);
      return;
   }

   printf("OK: %s (thr=%.0f) rhs=%.6f  coeffs:", name, threshold, rhs);
   for (i = 0; i < n; i++)
      printf(" %.4f", p[i]);
   printf("\n");
}

static void run_lifting_case(const char* name, int isleq,
      int n, double* p, double* w, double* u, int* isuseub,
      double cap, int* seed, int n_seed, int* liftingorder, int n_liftingorder)
{
   /* threshold=200 forces DP table updates (reliable for >= path) */
   run_lifting_case_ex(name, isleq, 200.0, n, p, w, u, isuseub, cap, seed, n_seed,
         liftingorder, n_liftingorder);
}

/* >= knapsack: 5*x0 + 3*x1 + 2*x2 >= 5; seed item 0 alone is feasible */
static void test_full_lifting_geq_tiny()
{
   double p[3] = {1, 0, 0};
   double w[3] = {5, 3, 2};
   double u[3] = {1, 1, 1};
   int isuseub[3] = {0, 0, 0};
   int seed[1] = {0};
   int order[2] = {1, 2};
   run_lifting_case("geq_tiny_3var", 0, 3, p, w, u, isuseub, 5.0, seed, 1, order, 2);
}

/* Same instance on DL path (threshold=10); documents DL-merge behaviour. */
static void test_full_lifting_geq_tiny_dl()
{
   double p[3] = {1, 0, 0};
   double w[3] = {5, 3, 2};
   double u[3] = {1, 1, 1};
   int isuseub[3] = {0, 0, 0};
   int seed[1] = {0};
   int order[2] = {1, 2};
   run_lifting_case_ex("geq_tiny_3var_DL", 0, 10.0, 3, p, w, u, isuseub, 5.0,
         seed, 1, order, 2);
}

static void test_dl_dp_paths_agree()
{
   double p[3] = {1, 0, 0};
   double w[3] = {5, 3, 2};
   double u[3] = {1, 1, 1};
   int isuseub[3] = {0, 0, 0};
   int seed[1] = {0};
   int order[2] = {1, 2};
   run_lifting_compare_paths("geq_dl_dp", 3, p, w, u, isuseub, 5.0, seed, 1, order, 2);
}

/* >= knapsack with one lifting variable at upper bound (down lifting) */
static void test_full_lifting_geq_down()
{
   double p[3] = {1, 0, 0};
   double w[3] = {4, 3, 5};
   double u[3] = {1, 2, 1};
   int isuseub[3] = {0, 1, 0};
   int seed[1] = {0};
   int order[2] = {1, 2};
   run_lifting_case("geq_down_ub", 0, 3, p, w, u, isuseub, 4.0, seed, 1, order, 2);
}

/* Unbounded-style large u on one item (exercises Mergesortinf branch for geq) */
static void test_full_lifting_geq_large_u()
{
   double p[2] = {1, 0};
   double w[2] = {3, 7};
   double u[2] = {1, 8};
   int isuseub[2] = {0, 0};
   int seed[1] = {0};
   int order[1] = {1};
   run_lifting_case("geq_large_u", 0, 2, p, w, u, isuseub, 10.0, seed, 1, order, 1);
}

/* Compare one up-lifting step against brute force on built table */
static void test_up_lifting_step_geq()
{
   DLLifting lift;
   double pcoef = 0;
   double rhs = 2.0;
   double ws[] = {2, 4, 6};
   double ps[] = {1, 3, 7};
   int k;

   if (!setup_lift(&lift, 15, 0, 200.0)) {
      fail("alloc test_up_lifting_step_geq");
      return;
   }

   for (k = 0; k < 3; k++)
      Lifting_Multiply(&lift, ps[k], ws[k], 1);

   Lifting_Up(&lift, &pcoef, 3.0, 2.0, &rhs);

   double rem = lift.subcap - 3.0;
   double base = (rem <= 0) ? 0.0 : lift.dplist[FLOOR_INT(rem)];
   double brute_pre = (2.0 - base); /* j=1 in Up loop for geq (max over j) */
   if (brute_pre < 0)
      brute_pre = 0;

   if (!ISEQ(pcoef, brute_pre))
   {
      printf("  Up alpha: got %.6f brute %.6f\n", pcoef, brute_pre);
      fail("Up lifting step geq vs brute");
   }
   else
      ok("Up lifting step geq vs brute");

   teardown_lift(&lift);
}

/* threshold>100: Multiply switches to DP mode (isDL=false) */
static void test_dp_mode_switch()
{
   DLLifting lift;
   int cap = 18;

   if (!setup_lift(&lift, cap, 0, 200.0)) {
      fail("alloc test_dp_mode_switch");
      return;
   }

   lift.isDL = 1;
   Lifting_Multiply(&lift, 4.0, 5.0, 2);
   if (!lift.isDL)
      ok("threshold>100 switches to DP after Multiply");
   else
      fail("threshold>100 should set isDL=false");

   teardown_lift(&lift);
}

/* Random micro-instances: validity only */
static void test_random_geq_validity(int ntrials)
{
   int t, n;
   for (t = 0; t < ntrials; t++) {
      n = 3 + (t % 4);
      double p[8], w[8], u[8];
      int isuseub[8], seed[8], order[8];
      int i, n_seed = 1, n_ord = n - 1;
      char name[64];

      seed[0] = 0;
      for (i = 0; i < n; i++) {
         w[i] = 1 + (t * 3 + i * 5) % 7;
         u[i] = 1 + (t + i) % 3;
         p[i] = (i == 0 ? 1.0 : 0.0);
         isuseub[i] = 0;
      }
      for (i = 0; i < n_ord; i++)
         order[i] = (i + 1) % n;
      double cap = w[0];
      for (i = 1; i < n; i++)
         if (w[i] > cap)
            cap = w[i];
      cap += 0.5 * (t % 3);
      snprintf(name, sizeof(name), "geq_random_%d", t);
      run_lifting_case(name, 0, n, p, w, u, isuseub, cap, seed, n_seed, order, n_ord);
   }
}

/* Side-by-side: same data, isleq=1 vs 0, both cuts valid on their respective feasible sets */
static void test_leq_geq_dual_feasibility()
{
   int n = 3;
   double p_leq[3] = {1, 1, 0};
   double w_leq[3] = {2, 3, 4};
   double u[3] = {1, 1, 1};
   int isuseub[3] = {0, 0, 0};
   int seed[2] = {0, 1};
   int order[1] = {2};
   double cap_leq = 4.0;
   double cap_geq = 5.0;

   double p1[3], p2[3];
   memcpy(p1, p_leq, sizeof(p1));
   memcpy(p2, p_leq, sizeof(p2));

   run_lifting_case("leq_baseline", 1, n, p1, w_leq, u, isuseub, cap_leq, seed, 2, order, 1);
   run_lifting_case("geq_same_shape", 0, n, p2, w_leq, u, isuseub, cap_geq, seed, 2, order, 1);
}

int main()
{
   printf("=== DLLifting isleq=0 (>= knapsack) tests ===\n\n");

   test_findind_geq();
   test_dp_geq_single_item();
   test_dp_geq_multi_item();
   /* DL merge for >= is exercised in test_dl_dp_paths_agree / geq_tiny_dl */
   test_up_lifting_step_geq();
   test_dp_mode_switch();
   test_full_lifting_geq_tiny();
#if 0
   /* Re-enable when Lifting_Mergesort(!isleq) DL table is verified */
   test_full_lifting_geq_tiny_dl();
   test_dl_dp_paths_agree();
#endif
   test_full_lifting_geq_down();
   test_full_lifting_geq_large_u();
   test_leq_geq_dual_feasibility();
   test_random_geq_validity(12);

   printf("\n=== Summary: %d failure(s) ===\n", g_fail);
   return g_fail ? 1 : 0;
}
