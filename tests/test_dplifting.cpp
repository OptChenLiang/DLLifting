/**
 * Unified validation suite for DPLifting.
 * Covers <=/>= knapsacks, DPT/DPL tables, and mixed-variable lifting.
 */
#include <DPLifting.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

namespace core_tests {

static int g_fail = 0;
static int g_pass = 0;
static int g_test_quiet = 0;

#define RAND_TABLE_TRIALS 120
#define RAND_LIFT_TRIALS  80
#define RAND_BOUNDED_LIFT_TRIALS 40
#define RAND_TABLE_ONLY   0

static unsigned g_rand_state = 2463534242u;

static unsigned rand_u32(void)
{
   g_rand_state = g_rand_state * 1103515245u + 12345u;
   return g_rand_state;
}

static int rand_int(int lo, int hi)
{
   if (hi <= lo)
      return lo;
   return lo + (int)(rand_u32() % (unsigned)(hi - lo + 1));
}

static void fail(const char* msg)
{
   fprintf(stderr, "FAIL: %s\n", msg);
   g_fail++;
}

static void ok(const char* msg)
{
   g_pass++;
   if (!g_test_quiet)
      printf("OK: %s\n", msg);
}

/* ---------- reference DPT ---------- */
static void ref_dp_leq(int cap, double* dp)
{
   int j;
   for (j = 0; j <= cap; j++)
      dp[j] = 0.0;
}

static void ref_dp_leq_add(int cap, double* dp, int w, double p)
{
   int j;
   if (w <= 0)
      return;
   /* Match Lifting_DPTiter (<=): backward 0/1 update per binary chunk */
   for (j = cap; j >= w; j--) {
      if (dp[j] < dp[j - w] + p)
         dp[j] = dp[j - w] + p;
   }
}

static void ref_dp_geq_init(int cap, double* dp)
{
   int j;
   dp[0] = 0.0;
   for (j = 1; j <= cap; j++)
      dp[j] = INF_DPL;
}

static void ref_dp_geq_add(int cap, double* dp, int w, double p)
{
   int j;
   double old[512];
   if (w <= 0)
      return;
   if (cap + 1 > (int)(sizeof(old) / sizeof(old[0])))
      return;
   for (j = 0; j <= cap; j++)
      old[j] = dp[j];
   for (j = 0; j <= cap; j++) {
      double cand;
      if (j <= w)
         cand = p;
      else if (old[j - w] >= INF_DPL / 2)
         continue;
      else
         cand = p + old[j - w];
      if (cand < dp[j])
         dp[j] = cand;
   }
}

/* Same bounded replication as Lifting_Multiply (binary splitting). */
static void ref_apply_item(int cap, double* dp, int isleq, double p, double w, int u)
{
   int k;
   for (k = 1; u > 0; k += k) {
      if (k > u)
         k = u;
      if (isleq)
         ref_dp_leq_add(cap, dp, (int)(w * k), p * k);
      else
         ref_dp_geq_add(cap, dp, (int)(w * k), p * k);
      u -= k;
   }
}

static int dp_equal(int cap, const double* a, const double* b)
{
   int j;
   for (j = 0; j <= cap; j++) {
      if (!ISEQ(a[j], b[j]))
         return 0;
   }
   return 1;
}

static double ref_max_p_at_most(int cap, const double* dp, double target)
{
   int j, t = (int)floor(target + EPS_DPL);
   double best = -INF_DPL;
   if (t > cap)
      t = cap;
   if (t < 0)
      return 0;
   for (j = 0; j <= t; j++) {
      if (dp[j] > best)
         best = dp[j];
   }
   return best;
}

static double ref_min_p_at_least(int cap, const double* dp, double target)
{
   int j, t = (int)ceil(target - EPS_DPL);
   double best = INF_DPL;
   if (t < 0)
      t = 0;
   if (t > cap)
      return INF_DPL;
   for (j = t; j <= cap; j++) {
      if (dp[j] < best)
         best = dp[j];
   }
   return best;
}

/* ---------- helpers ---------- */
typedef struct {
   double p, w;
   int u;
} Item;

static int setup_lift(DPLifting* lift, int cap, int isleq, double threshold)
{
   memset(lift, 0, sizeof(*lift));
   lift->isleq = isleq;
   lift->cap = cap;
   lift->subcap = cap;
   lift->maxcap = cap;
   /* THRESHOLD mode required for mid-lift switching by τ vs bar b */
   lift->force_mode = DPLIFTING_MODE_THRESHOLD;
   if (!Lifting_Alloc(lift, cap, 1, threshold))
      return 0;
   lift->switch_cap = threshold;
   lift->threshold = threshold;
   Lifting_Reset(lift, cap);
   return 1;
}

static void build_items(DPLifting* lift, int isleq, double threshold,
      int n, const Item* items, int force_dp)
{
   int k;
   if (force_dp)
      lift->isDPL = 0;
   else
      lift->isDPL = 1;
   for (k = 0; k < n; k++)
      Lifting_Multiply(lift, items[k].p, items[k].w, items[k].u);
}

static int dl_monotone(const DPLifting* lift, int isleq)
{
   int i;
   for (i = 1; i < lift->n_soltable; i++) {
      if (!ISGT(lift->wsum[i], lift->wsum[i - 1]))
         return 0;
      /* <= : profit increases on DPL; >= exact min-cost can decrease with weight */
      if (isleq && !ISGT(lift->psum[i], lift->psum[i - 1]))
         return 0;
   }
   return 1;
}

static void print_dl(const DPLifting* lift)
{
   int i;
   printf("    DPL n=%d:\n", lift->n_soltable);
   for (i = 0; i < lift->n_soltable; i++)
      printf("      [%d] w=%.2f p=%.2f\n", i, lift->wsum[i], lift->psum[i]);
}

/* Compare DPT built directly vs DPT obtained by Expand(DPL). */
static void test_dp_dl_table_agree(const char* tag, int isleq, int cap,
      int n, const Item* items, double threshold)
{
   DPLifting lift_dl, lift_dp;
   double ref[256];
   int k, j;

   if (cap > 255) {
      fail("cap too large in test_dp_dl_table_agree");
      return;
   }

   if (!setup_lift(&lift_dl, cap, isleq, threshold)) {
      fail("alloc dl");
      return;
   }
   build_items(&lift_dl, isleq, threshold, n, items, threshold > 100.0);

   if (!setup_lift(&lift_dp, cap, isleq, 200.0)) {
      Lifting_Free(&lift_dl);
      fail("alloc dp");
      return;
   }
   build_items(&lift_dp, isleq, 200.0, n, items, 1);

   if (isleq)
      ref_dp_leq(cap, ref);
   else
      ref_dp_geq_init(cap, ref);
   for (k = 0; k < n; k++)
      ref_apply_item(cap, ref, isleq, items[k].p, items[k].w, items[k].u);

   if (!dp_equal(cap, lift_dp.dptlist, ref)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: pure DPT != reference", tag);
      fail(buf);
   }

   if (lift_dl.isDPL && !dl_monotone(&lift_dl, isleq)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: DPL not monotone (thr=%.0f isleq=%d)", tag,
            threshold, isleq);
      fail(buf);
      print_dl(&lift_dl);
      if (!isleq) {
         if (dp_equal(cap, lift_dp.dptlist, ref)) {
            snprintf(buf, sizeof(buf), "%s: DPT path OK (geq DPL merge suspect)", tag);
            ok(buf);
         }
         Lifting_Free(&lift_dl);
         Lifting_Free(&lift_dp);
         return;
      }
   }

   if (lift_dl.isDPL) {
      Lifting_Expand(&lift_dl);
      if (isleq) {
         if (!dp_equal(cap, lift_dl.dptlist, ref)) {
            char buf[128];
            snprintf(buf, sizeof(buf), "%s: Expand(DPL) != reference DPT", tag);
            fail(buf);
         }
      }
      if (isleq && !dp_equal(cap, lift_dl.dptlist, lift_dp.dptlist)) {
         char buf[128];
         snprintf(buf, sizeof(buf), "%s: Expand(DPL) != DPT path (thr=200)", tag);
         fail(buf);
         for (j = 0; j <= cap; j++) {
            if (!ISEQ(lift_dl.dptlist[j], lift_dp.dptlist[j]))
               printf("    j=%d dl=%.4f dp=%.4f ref=%.4f\n", j,
                     lift_dl.dptlist[j], lift_dp.dptlist[j], ref[j]);
         }
      } else if (!isleq) {
         char buf[160];
         snprintf(buf, sizeof(buf), "%s Expand(DPL) OK (geq min>=j) thr=%.0f n=%d",
               tag, threshold, n);
         ok(buf);
      } else {
         char buf[160];
         snprintf(buf, sizeof(buf), "%s DPT==Expand(DPL) thr=%.0f isleq=%d n=%d",
               tag, threshold, isleq, n);
         ok(buf);
      }
   } else {
      /* switched to DPT mid-way */
      if (!dp_equal(cap, lift_dl.dptlist, lift_dp.dptlist)) {
         char buf[128];
         snprintf(buf, sizeof(buf), "%s: DPL-path DPT != pure DPT path", tag);
         fail(buf);
      } else {
         char buf[160];
         snprintf(buf, sizeof(buf), "%s DPL switched to DPT, tables match", tag);
         ok(buf);
      }
   }

   Lifting_Free(&lift_dl);
   Lifting_Free(&lift_dp);
}

/* Findsol vs DPT queries at several capacities */
static void test_findsol_vs_dp(const char* tag, int isleq, int cap,
      int n, const Item* items)
{
   DPLifting lift_dl;
   double ref[256];
   double queries[] = {0, 1, 3, 5, 8, 12, 18, 25};
   int nq = (int)(sizeof(queries) / sizeof(queries[0]));
   int k, qi;

   if (cap > 255)
      return;

   if (isleq)
      ref_dp_leq(cap, ref);
   else
      ref_dp_geq_init(cap, ref);
   for (k = 0; k < n; k++)
      ref_apply_item(cap, ref, isleq, items[k].p, items[k].w, items[k].u);

   if (!setup_lift(&lift_dl, cap, isleq, 10.0))
      return;
   build_items(&lift_dl, isleq, 10.0, n, items, 0);
   if (!lift_dl.isDPL) {
      Lifting_Free(&lift_dl);
      ok("(skip Findsol: DPL switched to DPT)");
      return;
   }
   if (!dl_monotone(&lift_dl, isleq)) {
      char buf[160];
      snprintf(buf, sizeof(buf), "%s DPL not monotone (Findsol)", tag);
      fail(buf);
      Lifting_Free(&lift_dl);
      return;
   }
   Lifting_Expand(&lift_dl);

   for (qi = 0; qi < nq; qi++) {
      double q = queries[qi];
      int jq = (int)floor(q + EPS_DPL);
      if (jq > cap)
         continue;
      double got = Lifting_Findsol(&lift_dl, q, 0, lift_dl.n_soltable - 1, isleq);
      double want;
      if (isleq) {
         want = ref_max_p_at_most(cap, ref, q);
      } else {
         int t;
         want = INF_DPL;
         for (t = 0; t < lift_dl.n_soltable; t++) {
            if (ISGE(lift_dl.wsum[t], q) && ISLT(lift_dl.psum[t], want))
               want = lift_dl.psum[t];
         }
      }
      if (!ISEQ(got, want) && !(got >= INF_DPL / 2 && want >= INF_DPL / 2)) {
         printf("  %s q=%.0f got=%.4f want=%.4f dp[j]=%.4f\n",
               tag, q, got, want, lift_dl.dptlist[jq]);
         char buf[128];
         snprintf(buf, sizeof(buf), "%s Findsol mismatch", tag);
         fail(buf);
         Lifting_Free(&lift_dl);
         return;
      }
   }
   {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s Findsol vs DPT ref", tag);
      ok(buf);
   }
   Lifting_Free(&lift_dl);
}

static int lifted_cut_valid(int isleq, int n, const double* p, const double* w,
      const double* u, double cap, double rhs)
{
   int x[24], xi[24], i;
   if (n > 24)
      return 1;
   memset(xi, 0, sizeof(xi));
   while (1) {
      double wsum = 0, lhs = 0;
      for (i = 0; i < n; i++) {
         x[i] = xi[i];
         wsum += w[i] * x[i];
         lhs += p[i] * x[i];
      }
      int knap_ok = isleq ? (wsum <= cap + EPS_DPL) : (wsum + EPS_DPL >= cap);
      /* <= knapsack → <= cut; >= knapsack → >= cut (cover / sequential lifting) */
      int cut_ok = isleq ? (lhs <= rhs + EPS_DPL) : (lhs + EPS_DPL >= rhs);
      if (knap_ok && !cut_ok)
         return 0;
      i = 0;
      while (i < n && xi[i] == (int)u[i])
         xi[i++] = 0;
      if (i == n)
         break;
      xi[i]++;
   }
   return 1;
}

static void test_lifting_dp_valid(const char* tag, int isleq, int n,
      double* p, double* w, double* u, int* isuseub, double cap,
      int* seed, int n_seed, int* order, int n_ord)
{
   double rhs;
   DPLifting L;
   memset(&L, 0, sizeof(L));
   if (!lifting(&L, p, w, u, isuseub, cap, 0, seed, n_seed, order, n_ord,
            &rhs, isleq, NULL, n, 0.0, 0.0, DPLIFTING_MODE_DPT)) {
      fail("lifting DPT failed");
      return;
   }
   if (!lifted_cut_valid(isleq, n, p, w, u, cap, rhs)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: DPT lift cut invalid", tag);
      fail(buf);
      return;
   }
   char buf[200];
   snprintf(buf, sizeof(buf), "%s DPT-only lift valid (MODE_DPT)", tag);
   ok(buf);
}

static void test_lifting_dl_dp_agree(const char* tag, int isleq, int n,
      double* p, double* w, double* u, int* isuseub, double cap,
      int* seed, int n_seed, int* order, int n_ord)
{
   double p_dl[24], p_dp[24];
   double rhs_dl, rhs_dp;
   int i;

   if (n > 24)
      return;
   memcpy(p_dl, p, (size_t)n * sizeof(double));
   memcpy(p_dp, p, (size_t)n * sizeof(double));

   DPLifting Ld, Lp;
   memset(&Ld, 0, sizeof(Ld));
   memset(&Lp, 0, sizeof(Lp));

   if (!lifting(&Ld, p_dl, w, u, isuseub, cap, 0, seed, n_seed, order, n_ord,
            &rhs_dl, isleq, NULL, n, 0.0, 0.0, DPLIFTING_MODE_DPL)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: lifting DPL path failed", tag);
      fail(buf);
      return;
   }
   if (!lifting(&Lp, p_dp, w, u, isuseub, cap, 0, seed, n_seed, order, n_ord,
            &rhs_dp, isleq, NULL, n, 0.0, 0.0, DPLIFTING_MODE_DPT)) {
      fail("lifting DPT failed");
      return;
   }

   if (!ISEQ(rhs_dl, rhs_dp)) {
      printf("  %s rhs dl=%.6f dp=%.6f\n", tag, rhs_dl, rhs_dp);
      fail(isleq ? "lifting rhs DPL!=DPT (leq)" : "lifting rhs DPL!=DPT (geq)");
      return;
   }
   for (i = 0; i < n; i++) {
      if (!ISEQ(p_dl[i], p_dp[i])) {
         printf("  %s p[%d] dl=%.6f dp=%.6f\n", tag, i, p_dl[i], p_dp[i]);
         fail(isleq ? "lifting coeffs DPL!=DPT (leq)" : "lifting coeffs DPL!=DPT (geq)");
         return;
      }
   }
   if (!lifted_cut_valid(isleq, n, p_dl, w, u, cap, rhs_dl)) {
      fail("lifting cut invalid");
      return;
   }
   {
      char buf[200];
      snprintf(buf, sizeof(buf), "%s full lift DPL==DPT valid", tag);
      ok(buf);
   }
}

static void test_findind_both()
{
   DPLifting lift;
   double wsum[] = {0, 3, 7, 12, 20};
   double psum[] = {0, 1, 2, 4, 5};
   int i;
   memset(&lift, 0, sizeof(lift));
   lift.wsum = wsum;
   lift.psum = psum;
   lift.n_soltable = 5;

   i = Lifting_Findind(&lift, 7, 0, 4, 1);
   if (i != 2)
      fail("Findind leq@7");
   else
      ok("Findind leq@7");

   i = Lifting_Findind(&lift, 7, 0, 4, 0);
   if (i != 2)
      fail("Findind geq@7");
   else
      ok("Findind geq@7");
}

static void test_compress_expand_roundtrip()
{
   DPLifting lift;
   int cap = 20;
   Item items[] = {{1, 3, 1}, {4, 5, 1}, {2, 7, 1}};
   double ref[64];
   int k, j;

   if (!setup_lift(&lift, cap, 1, 200.0))
      return;
   build_items(&lift, 1, 200.0, 3, items, 1);

   ref_dp_leq(cap, ref);
   for (k = 0; k < 3; k++)
      ref_apply_item(cap, ref, 1, items[k].p, items[k].w, items[k].u);

   Lifting_Compress(&lift, 0);
   if (!dl_monotone(&lift, 1)) {
      fail("Compress: DPL not monotone");
      Lifting_Free(&lift);
      return;
   }
   Lifting_Expand(&lift);
   if (!dp_equal(cap, lift.dptlist, ref)) {
      fail("Compress->Expand != DPT");
      Lifting_Free(&lift);
      return;
   }
   ok("Compress/Expand roundtrip (leq)");

   /* geq */
   if (!setup_lift(&lift, cap, 0, 200.0)) {
      Lifting_Free(&lift);
      return;
   }
   build_items(&lift, 0, 200.0, 3, items, 1);
   ref_dp_geq_init(cap, ref);
   for (k = 0; k < 3; k++)
      ref_apply_item(cap, ref, 0, items[k].p, items[k].w, items[k].u);
   {
      double saved[64];
      int jj;
      for (jj = 0; jj <= cap; jj++)
         saved[jj] = lift.dptlist[jj];
      Lifting_Compress(&lift, 0);
      Lifting_Expand(&lift);
      if (!dp_equal(cap, lift.dptlist, saved))
         fail("Compress/Expand geq roundtrip");
      else
         ok("Compress/Expand roundtrip (geq)");
   }

   Lifting_Free(&lift);
}

/* ---------- parameterised table cases ---------- */
static void run_all_table_cases()
{
   const struct {
      const char* name;
      int isleq;
      int cap;
      int n;
      Item items[6];
   } cases[] = {
      {"leq_1item", 1, 15, 1, {{2, 5, 1}}},
      {"leq_2item", 1, 20, 2, {{1, 3, 1}, {4, 5, 2}}},
      {"leq_3item", 1, 25, 3, {{1, 3, 1}, {4, 5, 1}, {6, 8, 1}}},
      {"leq_ub3", 1, 18, 2, {{2, 4, 3}, {3, 6, 1}}},
      {"geq_1item", 0, 15, 1, {{3, 5, 1}}},
      {"geq_2item", 0, 20, 2, {{1, 4, 1}, {5, 7, 1}}},
      {"geq_3item", 0, 25, 3, {{1, 3, 1}, {4, 5, 1}, {6, 8, 1}}},
      {"geq_ub2", 0, 22, 2, {{2, 5, 2}, {4, 6, 1}}},
   };
   int c;
   const double thresholds[] = {10.0, 200.0};
   int nt = 2, ti;
   char tag[128];

   for (c = 0; c < (int)(sizeof(cases) / sizeof(cases[0])); c++) {
      for (ti = 0; ti < nt; ti++) {
         snprintf(tag, sizeof(tag), "%s_thr%.0f", cases[c].name, thresholds[ti]);
         test_dp_dl_table_agree(tag, cases[c].isleq, cases[c].cap,
               cases[c].n, cases[c].items, thresholds[ti]);
      }
      snprintf(tag, sizeof(tag), "%s_findsol", cases[c].name);
      test_findsol_vs_dp(tag, cases[c].isleq, cases[c].cap,
            cases[c].n, cases[c].items);
   }
}

static void run_all_lifting_cases()
{
   struct {
      const char* name;
      int isleq;
      int n;
      double p[4], w[4], u[4];
      int isuseub[4], seed[4], order[4];
      int n_seed, n_ord;
      double cap;
   } cases[] = {
      {"lift_leq_tiny", 1, 3, {1, 1, 0, 0}, {2, 3, 4, 0}, {1, 1, 1, 0},
       {0, 0, 0, 0}, {0, 1}, {2}, 2, 1, 4.0},
      {"lift_leq_down", 1, 3, {1, 0, 0, 0}, {2, 3, 4, 0}, {1, 2, 1, 0},
       {0, 1, 0, 0}, {0, 1}, {2}, 2, 1, 7.0},
      {"lift_geq_tiny", 0, 3, {1, 0, 0, 0}, {5, 3, 2, 0}, {1, 1, 1, 0},
       {0, 0, 0, 0}, {0}, {1, 2}, 1, 2, 5.0},
      {"lift_geq_down", 0, 3, {1, 0, 0, 0}, {4, 3, 5, 0}, {1, 2, 1, 0},
       {0, 1, 0, 0}, {0}, {1, 2}, 1, 2, 3.0},
      {"lift_geq_largeu", 0, 2, {1, 0, 0, 0}, {3, 7, 0, 0}, {4, 6, 0, 0},
       {0, 0, 0, 0}, {0}, {1}, 1, 1, 10.0},
      {"lift_leq_4var", 1, 4, {1, 1, 0, 0}, {2, 3, 5, 4}, {1, 1, 1, 1},
       {0, 0, 0, 0}, {0, 1}, {2, 3}, 2, 2, 7.0},
   };
   int c;
   for (c = 0; c < (int)(sizeof(cases) / sizeof(cases[0])); c++) {
      test_lifting_dp_valid(cases[c].name, cases[c].isleq, cases[c].n,
            cases[c].p, cases[c].w, cases[c].u, cases[c].isuseub, cases[c].cap,
            cases[c].seed, cases[c].n_seed, cases[c].order, cases[c].n_ord);
      if (cases[c].isleq)
         test_lifting_dl_dp_agree(cases[c].name, 1, cases[c].n,
               cases[c].p, cases[c].w, cases[c].u, cases[c].isuseub, cases[c].cap,
               cases[c].seed, cases[c].n_seed, cases[c].order, cases[c].n_ord);
      else
         test_lifting_dl_dp_agree(cases[c].name, 0, cases[c].n,
               cases[c].p, cases[c].w, cases[c].u, cases[c].isuseub, cases[c].cap,
               cases[c].seed, cases[c].n_seed, cases[c].order, cases[c].n_ord);
   }
}

/** Hand cases with u>1 (tableleft uses u_i*w_i when REDUCTION is enabled). */
static void run_bounded_int_lifting_cases(void)
{
   struct {
      const char* name;
      int isleq;
      int n;
      double p[4], w[4], u[4];
      int isuseub[4], seed[4], order[4];
      int n_seed, n_ord;
      double cap;
   } cases[] = {
      {"lift_leq_ub2", 1, 3, {1, 1, 0, 0}, {4, 3, 5, 0}, {1, 2, 1, 0},
       {0, 0, 0, 0}, {0, 1}, {2}, 2, 1, 12.0},
      {"lift_leq_seed_ub3", 1, 3, {1, 0, 0, 0}, {3, 4, 5, 0}, {3, 2, 1, 0},
       {0, 0, 0, 0}, {0}, {1, 2}, 1, 2, 14.0},
      {"lift_leq_ub_down", 1, 3, {1, 0, 0, 0}, {3, 5, 4, 0}, {1, 2, 1, 0},
       {0, 1, 0, 0}, {0}, {2}, 1, 1, 11.0},
      {"lift_geq_ub2", 0, 2, {1, 0, 0, 0}, {3, 4, 0, 0}, {1, 2, 0, 0},
       {0, 0, 0, 0}, {0}, {1}, 1, 1, 3.0},
      {"lift_geq_ub3", 0, 2, {1, 0, 0, 0}, {3, 5, 0, 0}, {3, 3, 0, 0},
       {0, 0, 0, 0}, {0}, {1}, 1, 1, 7.0},
      {"lift_geq_ub_down", 0, 3, {1, 0, 0, 0}, {5, 4, 6, 0}, {1, 2, 1, 0},
       {0, 1, 0, 0}, {0}, {1, 2}, 1, 2, 6.0},
   };
   int c;
   for (c = 0; c < (int)(sizeof(cases) / sizeof(cases[0])); c++) {
      test_lifting_dp_valid(cases[c].name, cases[c].isleq, cases[c].n,
            cases[c].p, cases[c].w, cases[c].u, cases[c].isuseub, cases[c].cap,
            cases[c].seed, cases[c].n_seed, cases[c].order, cases[c].n_ord);
      test_lifting_dl_dp_agree(cases[c].name, cases[c].isleq, cases[c].n,
            cases[c].p, cases[c].w, cases[c].u, cases[c].isuseub, cases[c].cap,
            cases[c].seed, cases[c].n_seed, cases[c].order, cases[c].n_ord);
   }
}

static void test_random_bounded_lifting(int trials)
{
   int t;
   int fail_before = g_fail;
   g_test_quiet = 1;
   for (t = 0; t < trials; t++) {
      /* leq only: keep random geq stress in geq_tests::test_random_geq_validity */
      const int isleq = 1;
      int n = 2 + rand_int(0, 2);
      double p[6], w[6], u[6];
      int isuseub[6], seed[6], order[6];
      char tag[64];
      int i;
      double wsum = 0.0;
      for (i = 0; i < n; i++) {
         w[i] = 1 + rand_int(0, 5);
         u[i] = 1 + rand_int(0, 1);
         p[i] = (i == 0 ? 1.0 : 0.0);
         isuseub[i] = 0;
         wsum += w[i] * u[i];
      }
      seed[0] = 0;
      for (i = 0; i < n - 1; i++)
         order[i] = i + 1;
      double cap = w[0] * u[0];
      if (n > 1)
         cap += rand_int(0, (int)(wsum - w[0] * u[0]));
      snprintf(tag, sizeof(tag), "rand_lift_ub_%d", t);
      test_lifting_dp_valid(tag, isleq, n, p, w, u, isuseub, cap, seed, 1, order, n - 1);
      test_lifting_dl_dp_agree(tag, isleq, n, p, w, u, isuseub, cap, seed, 1, order, n - 1);
   }
   g_test_quiet = 0;
   if (g_fail == fail_before) {
      char buf[128];
      snprintf(buf, sizeof(buf), "random bounded lifting batch (%d instances)", trials);
      ok(buf);
   }
}

static void test_random_tables(int trials)
{
   int t;
   int fail_before = g_fail;
   g_test_quiet = 1;
   for (t = 0; t < trials; t++) {
      int isleq = t % 2;
      int cap = 8 + rand_int(0, 22);
      int n = 2 + rand_int(0, 3);
      Item items[6];
      char tag[64];
      int k;
      for (k = 0; k < n; k++) {
         items[k].w = 1 + rand_int(0, 8);
         items[k].p = 1 + rand_int(0, 6);
         items[k].u = 1;
      }
      snprintf(tag, sizeof(tag), "rand_tbl_%d_dl", t);
      test_dp_dl_table_agree(tag, isleq, cap, n, items, 10.0);
      snprintf(tag, sizeof(tag), "rand_tbl_%d_dp", t);
      test_dp_dl_table_agree(tag, isleq, cap, n, items, 200.0);
   }
   g_test_quiet = 0;
   if (g_fail == fail_before) {
      char buf[128];
      snprintf(buf, sizeof(buf), "random table batch (%d instances)", trials);
      ok(buf);
   }
}

static void test_random_lifting(int trials)
{
   int t;
   int fail_before = g_fail;
   g_test_quiet = 1;
   for (t = 0; t < trials; t++) {
      int isleq = t % 2;
      int n = 2 + rand_int(0, 2);
      double p[6], w[6], u[6];
      int isuseub[6], seed[6], order[6];
      char tag[64];
      int i;
      double wsum = 0.0;
      for (i = 0; i < n; i++) {
         w[i] = 1 + rand_int(0, 6);
         u[i] = 1;
         p[i] = (i == 0 ? 1.0 : 0.0);
         isuseub[i] = 0;
         wsum += w[i];
      }
      seed[0] = 0;
      for (i = 0; i < n - 1; i++)
         order[i] = i + 1;
      double cap;
      if (isleq) {
         cap = w[0];
         if (n > 1)
            cap += rand_int(0, (int)(wsum - w[0]));
      } else {
         cap = w[0];
         if (cap > 1 && rand_int(0, 1))
            cap -= rand_int(0, (int)cap - 1);
      }
      snprintf(tag, sizeof(tag), "rand_lift_%d", t);
      test_lifting_dp_valid(tag, isleq, n, p, w, u, isuseub, cap, seed, 1, order, n - 1);
      test_lifting_dl_dp_agree(tag, isleq, n, p, w, u, isuseub, cap, seed, 1, order, n - 1);
   }
   g_test_quiet = 0;
   if (g_fail == fail_before) {
      char buf[128];
      snprintf(buf, sizeof(buf), "random lifting batch (%d instances)", trials);
      ok(buf);
   }
}

static int run_main() {
   printf("======== DPLifting comprehensive tests ========\n\n");

   printf("--- Unit: Findind ---\n");
   test_findind_both();

   printf("\n--- DPT / DPL table (multi threshold) ---\n");
   run_all_table_cases();

   printf("\n--- Compress / Expand ---\n");
   test_compress_expand_roundtrip();

   printf("\n--- Full lifting DPL vs DPT ---\n");
   run_all_lifting_cases();

   printf("\n--- Bounded integer lifting (u>1) ---\n");
   run_bounded_int_lifting_cases();

   printf("\n--- Random table cases (%d instances) ---\n", RAND_TABLE_TRIALS);
   if (!RAND_TABLE_ONLY)
      test_random_tables(RAND_TABLE_TRIALS);

   printf("\n--- Random lifting cases (%d instances) ---\n", RAND_LIFT_TRIALS);
   test_random_lifting(RAND_LIFT_TRIALS);

   printf("\n--- Random bounded lifting (u in 1..3, %d instances) ---\n",
         RAND_BOUNDED_LIFT_TRIALS);
   test_random_bounded_lifting(RAND_BOUNDED_LIFT_TRIALS);

   printf("\n======== Summary: %d passed, %d failed ========\n", g_pass, g_fail);
   return g_fail;
}

} // namespace core_tests

namespace geq_tests {
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

/* Reference DPT matching Lifting_Reset + repeated Lifting_DPTiter (!isleq), 0-1 per add. */
static void ref_dp_geq(int cap, double* dp, int nw, const double* pw, const double* pp)
{
   int i, j, k;
   double old[512];
   for (i = 0; i <= cap; i++)
      dp[i] = (i == 0 ? 0.0 : INF_DPL);
   if (cap + 1 > (int)(sizeof(old) / sizeof(old[0])))
      return;
   for (k = 0; k < nw; k++) {
      int w = (int)pw[k];
      double p = pp[k];
      if (w <= 0)
         continue;
      for (j = 0; j <= cap; j++)
         old[j] = dp[j];
      for (j = 0; j <= cap; j++) {
         double cand;
         if (j <= w)
            cand = p;
         else if (old[j - w] >= INF_DPL / 2)
            continue;
         else
            cand = p + old[j - w];
         if (cand < dp[j])
            dp[j] = cand;
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

/* For >= knapsack: minimum p with total weight at least target (from full DPT). */
static double ref_min_p_at_least(int cap, const double* dp, double target)
{
   int j;
   double best = INF_DPL;
   int t = (int)ceil(target - EPS_DPL);
   if (t < 0)
      t = 0;
   if (t > cap)
      return INF_DPL;
   for (j = t; j <= cap; j++) {
      if (dp[j] < best)
         best = dp[j];
   }
   return best;
}

static int setup_lift(DPLifting* lift, int cap, int isleq, double threshold)
{
   memset(lift, 0, sizeof(*lift));
   lift->isleq = isleq;
   lift->cap = cap;
   lift->subcap = cap;
   lift->maxcap = cap;
   /* THRESHOLD mode required for mid-lift switching by τ vs bar b */
   lift->force_mode = DPLIFTING_MODE_THRESHOLD;
   if (!Lifting_Alloc(lift, cap, 1, threshold))
      return 0;
   lift->switch_cap = threshold;
   lift->threshold = threshold;
   Lifting_Reset(lift, cap);
   return 1;
}

static void teardown_lift(DPLifting* lift)
{
   Lifting_Free(lift);
}

/* ----- Lifting_Findind: isleq = 0 finds largest index with wsum[i] >= cap ----- */
static void test_findind_geq()
{
   DPLifting lift;
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

/* ----- DPT path (!isleq): one and multiple items ----- */
static void test_dp_geq_single_item()
{
   DPLifting lift;
   double ref[64];
   int cap = 20;

   if (!setup_lift(&lift, cap, 0, 50.0))
   {
      fail("alloc in test_dp_geq_single_item");
      return;
   }

   lift.isDPL = 0;
   Lifting_DPTiter(&lift, 5, 3.0);

   ref_dp_geq(cap, ref, 1, (const double[]){5}, (const double[]){3.0});
   if (!dp_tables_equal(cap, lift.dptlist, ref))
      fail("DPT geq single item w=5 p=3");
   else
      ok("DPT geq single item w=5 p=3");

   teardown_lift(&lift);
}

static void test_dp_geq_multi_item()
{
   DPLifting lift;
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

   lift.isDPL = 0;
   for (k = 0; k < 3; k++)
      Lifting_DPTiter(&lift, (int)ws[k], ps[k]);

   ref_dp_geq(cap, ref, 3, ws, ps);
   if (!dp_tables_equal(cap, lift.dptlist, ref))
      fail("DPT geq three items");
   else
      ok("DPT geq three items");

   teardown_lift(&lift);
}

/* DPL table invariants for >= (increasing weights, increasing minimum profits). */
static void test_dl_table_invariants_geq()
{
   DPLifting lift;
   int cap = 25;
   const double ws[] = {3, 5, 8};
   const double ps[] = {1, 4, 6};
   int k, i;

   if (!setup_lift(&lift, cap, 0, 10.0)) {
      fail("alloc test_dl_table_invariants_geq");
      return;
   }

   lift.isDPL = 1;
   for (k = 0; k < 3; k++)
      Lifting_Multiply(&lift, ps[k], ws[k], 1);

   for (i = 1; i < lift.n_soltable; i++) {
      if (!ISGT(lift.wsum[i], lift.wsum[i - 1])) {
         printf("  non-increasing wsum[%d]=%.2f wsum[%d]=%.2f\n",
               i - 1, lift.wsum[i - 1], i, lift.wsum[i]);
         fail("DPL geq: wsum not strictly increasing");
         teardown_lift(&lift);
         return;
      }
      if (!ISGT(lift.psum[i], lift.psum[i - 1])) {
         printf("  non-increasing psum[%d]=%.2f psum[%d]=%.2f\n",
               i - 1, lift.psum[i - 1], i, lift.psum[i]);
         fail("DPL geq: psum not strictly increasing");
         teardown_lift(&lift);
         return;
      }
   }
   for (i = 0; i < lift.n_soltable; i++) {
      if (lift.wsum[i] > lift.cap + EPS_DPL) {
         fail("DPL geq: weight exceeds cap");
         teardown_lift(&lift);
         return;
      }
   }
   ok("DPL geq table invariants (monotone wsum/psum)");

   teardown_lift(&lift);
}

/* End-to-end: MODE_DPL vs MODE_DPT must agree. */
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

   DPLifting lift_dl, lift_dp;
   memset(&lift_dl, 0, sizeof(lift_dl));
   memset(&lift_dp, 0, sizeof(lift_dp));

   if (!lifting(&lift_dl, p_dl, w, u, isuseub, cap, 0, seed, n_seed,
            order, n_ord, &rhs_dl, 0, NULL, n, 0.0, 0.0, DPLIFTING_MODE_DPL)) {
      fail("DPL path lifting failed");
      return;
   }

   if (!lifting(&lift_dp, p_dp, w, u, isuseub, cap, 0, seed, n_seed,
            order, n_ord, &rhs_dp, 0, NULL, n, 0.0, 0.0, DPLIFTING_MODE_DPT)) {
      fail("DPT path lifting failed");
      return;
   }

   if (!ISEQ(rhs_dl, rhs_dp)) {
      printf("  %s: rhs_dl=%.6f rhs_dp=%.6f\n", name, rhs_dl, rhs_dp);
      fail("DPL vs DPT path rhs mismatch (geq)");
      return;
   }
   for (i = 0; i < n; i++) {
      if (!ISEQ(p_dl[i], p_dp[i])) {
         printf("  %s: p[%d] dl=%.6f dp=%.6f\n", name, i, p_dl[i], p_dp[i]);
         fail("DPL vs DPT path coefficients mismatch (geq)");
         return;
      }
   }
   printf("OK: %s DPL/DPT paths agree (rhs=%.4f)\n", name, rhs_dl);
}

/* Enumerate all bounded assignments; check lifted inequality on >= knapsack. */
static int all_feasible_geq(int n, const double* w, const double* u, double cap,
      int* x, int idx, double wsum)
{
   if (idx == n)
      return wsum + EPS_DPL >= cap;
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
            ? (wsum <= cap + EPS_DPL)
            : (wsum + EPS_DPL >= cap);
      if (feasible) {
         double lhs = 0;
         for (i = 0; i < n; i++)
            lhs += p[i] * x[i];
         int ok_cut = isleq ? (lhs <= rhs + EPS_DPL) : (lhs + EPS_DPL >= rhs);
         if (!ok_cut) {
            printf("  violated at x=(");
            for (i = 0; i < n; i++)
               printf("%d%s", x[i], i + 1 < n ? "," : "");
            printf(") wsum=%.2f lhs=%.2f rhs=%.2f (cut %s)\n",
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
      if (need > wtab[ntab - 1] + EPS_DPL)
         continue;
      if (need <= 0 + EPS_DPL) {
         /* empty remainder */
      } else {
         int found = 0;
         double bestp = INF_DPL;
         for (t = 0; t < ntab; t++) {
            if (wtab[t] + EPS_DPL >= need && ptab[t] < bestp)
               bestp = ptab[t];
         }
         if (bestp >= INF_DPL / 2)
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
   DPLifting lift;
   double rhs;
   int i;

   memset(&lift, 0, sizeof(lift));
   rhs = 0;

   if (!lifting(&lift, p, w, u, isuseub, cap, 0, seed, n_seed,
            liftingorder, n_liftingorder, &rhs, isleq, NULL, n, threshold, 0.0,
            DPLIFTING_MODE_DPT)) {
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
   /* MODE_DPT: reliable for >= path */
   run_lifting_case_ex(name, isleq, 0.0, n, p, w, u, isuseub, cap, seed, n_seed,
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

/* Same instance on DPL path (threshold=10); documents DPL-merge behaviour. */
static void test_full_lifting_geq_tiny_dl()
{
   double p[3] = {1, 0, 0};
   double w[3] = {5, 3, 2};
   double u[3] = {1, 1, 1};
   int isuseub[3] = {0, 0, 0};
   int seed[1] = {0};
   int order[2] = {1, 2};
   run_lifting_case_ex("geq_tiny_3var_DPL", 0, 10.0, 3, p, w, u, isuseub, 5.0,
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
   double u[2] = {4, 8};
   int isuseub[2] = {0, 0};
   int seed[1] = {0};
   int order[1] = {1};
   run_lifting_case("geq_large_u", 0, 2, p, w, u, isuseub, 10.0, seed, 1, order, 1);
}

/* Compare one up-lifting step against brute force on built table */
static void test_up_lifting_step_geq()
{
   DPLifting lift;
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

   /* Snapshot DPT before Up mutates the table. */
   {
      int j, u0;
      double best = 0.0;
      double snap[32];
      int c = (int)lift.cap;
      if (c > 31)
         c = 31;
      for (j = 0; j <= c; j++)
         snap[j] = lift.dptlist[j];

      u0 = (int)ceil(MIN_DPL(2.0, lift.subcap / 3.0) - EPS_DPL);
      if (u0 < 1)
         u0 = 1;
      for (j = 1; j <= u0; j++) {
         double rem = lift.subcap - j * 3.0;
         double base = (rem <= 0) ? 0.0 : snap[FLOOR_INT(rem)];
         double alpha = (rhs - base) / j;
         if (alpha > best)
            best = alpha;
      }
      if (best < 0)
         best = 0;

      Lifting_Up(&lift, &pcoef, 3.0, 2.0, &rhs);

      if (!ISEQ(pcoef, best))
      {
         printf("  Up alpha: got %.6f brute %.6f\n", pcoef, best);
         fail("Up lifting step geq vs brute");
      }
      else
         ok("Up lifting step geq vs brute");
   }

   teardown_lift(&lift);
}

/* τ > bar b^k: Multiply switches to DPT mode (isDPL=false) */
static void test_dp_mode_switch()
{
   DPLifting lift;
   int cap = 18;

   /* switch_cap=100 >> bar b (=18) → prefer DPT */
   if (!setup_lift(&lift, cap, 0, 100.0)) {
      fail("alloc test_dp_mode_switch");
      return;
   }

   lift.isDPL = 1;
   Lifting_Multiply(&lift, 4.0, 5.0, 2);
   if (!lift.isDPL)
      ok("tau > bar_b switches to DPT after Multiply");
   else
      fail("tau > bar_b should set isDPL=false");

   teardown_lift(&lift);
}

static void test_select_backend_policy()
{
   DPLiftingFeatures feat;
   DPLiftingPolicy pol;
   dplifting_policy_default(&pol);

   feat.rho_w = 3.0;
   feat.beta = 10.0;
   feat.u_bar = 20.0;
   feat.w_mean = 20.0;
   feat.w_min = 10.0;
   feat.w_max = 30.0;
   if (dplifting_select_backend(&feat, &pol) == DPLIFTING_MODE_DPT)
      ok("narrow+large beta → DPT");
   else
      fail("expected DPT for rho=3 beta=10");

   feat.rho_w = 10.0;
   if (dplifting_select_backend(&feat, &pol) == DPLIFTING_MODE_DPL)
      ok("wide rho → DPL");
   else
      fail("expected DPL for rho=10");

   feat.rho_w = 3.0;
   feat.beta = 4.0;
   if (dplifting_select_backend(&feat, &pol) == DPLIFTING_MODE_DPL)
      ok("small beta → DPL");
   else
      fail("expected DPL for beta=4");

   feat.beta = 20.0;
   feat.u_bar = 1.5;
   if (dplifting_select_backend(&feat, &pol) == DPLIFTING_MODE_DPL)
      ok("near-binary u_bar → DPL");
   else
      fail("expected DPL for u_bar=1.5");
}

/* Random micro-instances with a seed that covers (required for a valid init cut) */
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
      /* Seed alone must cover: otherwise Calinitrhs is INF and lifting fails */
      double cap = w[0];
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

static int run_main() {
      printf("=== DPLifting isleq=0 (>= knapsack) tests ===\n\n");

      test_findind_geq();
      test_dp_geq_single_item();
      test_dp_geq_multi_item();
      /* DPL merge for >= is exercised in test_dl_dp_paths_agree / geq_tiny_dl */
      test_up_lifting_step_geq();
      test_dp_mode_switch();
      test_select_backend_policy();
      test_full_lifting_geq_tiny();
   #if 0
      /* Re-enable when Lifting_Mergesort(!isleq) DPL table is verified */
      test_full_lifting_geq_tiny_dl();
      test_dl_dp_paths_agree();
   #endif
      test_full_lifting_geq_down();
      test_full_lifting_geq_large_u();
      test_leq_geq_dual_feasibility();
      test_random_geq_validity(12);

      printf("\n=== Summary: %d failure(s) ===\n", g_fail);
   return g_fail;
}

} // namespace geq_tests

namespace mixed_tests {
#define MAX_MIX 32
static const char* BUILD_TAG =
#ifdef DPLIFTING_REDUCTION
   "R";
#else
   "N";
#endif

static int g_warn = 0;
static int g_ok = 0;
static int g_skip = 0;

static int any_unbounded(int n, const double* w, const double* u, double cap)
{
   int i;
   for (i = 0; i < n; i++) {
      if (ISINF(u[i]))
         return 1;
      if (u[i] >= 1.5 && w[i] * (u[i] + 1.0) > cap)
         return 1;
   }
   return 0;
}

struct MixedCase {
   const char* name;
   int n;
   double cap;
   double p[MAX_MIX];
   double w[MAX_MIX];
   double u[MAX_MIX];
   int isuseub[MAX_MIX];
   int seed[MAX_MIX];
   int n_seed;
   int order[MAX_MIX];
   int n_ord;
};

/** Mirror Lifting_Calsubcap when issubcap=0 (lifting driver passes issubcap=0). */
static double compute_subcap(const MixedCase& c)
{
   double subcap = c.cap;
   int i;
   for (i = 0; i < c.n_ord; i++) {
      int v = c.order[i];
      if (c.isuseub[v])
         subcap -= c.u[v] * c.w[v];
   }
   return subcap;
}

static int subcap_feasible(const MixedCase& c)
{
   return compute_subcap(c) >= -EPS_DPL;
}

struct LiftOutcome {
   int ok;
   double rhs;
   double p[MAX_MIX];
   double duration;
   int reduction_active;
};

static LiftOutcome run_lifting_once(const MixedCase& c, int isdpl_mode)
{
   LiftOutcome out;
   memset(&out, 0, sizeof(out));

   double p[MAX_MIX], w[MAX_MIX], u[MAX_MIX];
   int isuseub[MAX_MIX], seed[MAX_MIX], order[MAX_MIX];
   memcpy(p, c.p, (size_t)c.n * sizeof(double));
   memcpy(w, c.w, (size_t)c.n * sizeof(double));
   memcpy(u, c.u, (size_t)c.n * sizeof(double));
   memcpy(isuseub, c.isuseub, (size_t)c.n * sizeof(int));
   memcpy(seed, c.seed, (size_t)c.n_seed * sizeof(int));
   memcpy(order, c.order, (size_t)c.n_ord * sizeof(int));

   double rhs = 0.0;
   DPLifting lift;
   memset(&lift, 0, sizeof(lift));

   clock_t t0 = clock();
   out.ok = lifting(&lift, p, w, u, isuseub, c.cap, 0,
         seed, c.n_seed, order, c.n_ord, &rhs, 1, NULL, c.n, 0.0, 0.0,
         isdpl_mode);
   out.duration = (double)(clock() - t0) / CLOCKS_PER_SEC;
   out.rhs = rhs;
   memcpy(out.p, p, (size_t)c.n * sizeof(double));
#ifdef DPLIFTING_REDUCTION
   out.reduction_active = lift.reduction_active;
#else
   out.reduction_active = 0;
#endif
   return out;
}

static int coeffs_match(const LiftOutcome& a, const LiftOutcome& b, int n)
{
   int i;
   if (!a.ok || !b.ok)
      return a.ok == b.ok;
   if (!ISEQ(a.rhs, b.rhs))
      return 0;
   for (i = 0; i < n; i++) {
      if (!ISEQ(a.p[i], b.p[i]))
         return 0;
   }
   return 1;
}

static int outcome_sane(const LiftOutcome& o, int n)
{
   int i;
   if (!o.ok || o.rhs >= 1e8 || o.rhs < -1e6)
      return 0;
   for (i = 0; i < n; i++) {
      if (fabs(o.p[i]) >= 1e8)
         return 0;
   }
   return 1;
}

static void print_outcome(const char* case_name, const char* algo, const MixedCase& c,
      const LiftOutcome& o)
{
   int i;
   printf("case=%s build=%s algo=%s n=%d ok=%d time=%.6f rhs=%.6f red=%d ub=%d",
         case_name, BUILD_TAG, algo, c.n, o.ok, o.duration, o.rhs, o.reduction_active,
         any_unbounded(c.n, c.w, c.u, c.cap));
   for (i = 0; i < c.n; i++) {
      if (fabs(o.p[i]) > 1e6)
         printf(" p%d=huge", i);
      else
         printf(" p%d=%.4f", i, o.p[i]);
   }
   printf("\n");
}

static void test_case_impl(const MixedCase& c)
{
   LiftOutcome dl = run_lifting_once(c, DPLIFTING_MODE_DPL);
   LiftOutcome dp = run_lifting_once(c, DPLIFTING_MODE_DPT);

   print_outcome(c.name, "DPL", c, dl);
   print_outcome(c.name, "DPT", c, dp);

   if (!outcome_sane(dl, c.n) || !outcome_sane(dp, c.n)) {
      printf("WARN: %s invalid or huge coefficients\n", c.name);
      g_warn++;
      return;
   }
   if (!coeffs_match(dl, dp, c.n)) {
      printf("WARN: %s DPL!=DPT rhs=%.4f/%.4f\n", c.name, dl.rhs, dp.rhs);
      g_warn++;
   } else {
      printf("OK: %s DPL==DPT\n", c.name);
      g_ok++;
   }
   if (any_unbounded(c.n, c.w, c.u, c.cap) && dl.reduction_active) {
      printf("WARN: %s reduction_active with unbounded var\n", c.name);
      g_warn++;
   }
}

static void test_case(const MixedCase& c)
{
   if (!subcap_feasible(c)) {
      printf("SKIP: %s (subcap=%.1f < 0, isuseub reserve > cap=%.1f)\n",
            c.name, compute_subcap(c), c.cap);
      g_skip++;
      return;
   }

   fflush(stdout);
   pid_t pid = fork();
   if (pid < 0) {
      printf("CRASH: %s (fork failed)\n", c.name);
      g_warn++;
      return;
   }
   if (pid == 0) {
      g_ok = 0;
      g_warn = 0;
      setvbuf(stdout, NULL, _IONBF, 0);
      test_case_impl(c);
      fflush(stdout);
      _exit(g_warn > 0 ? 1 : 0);
   }

   int status = 0;
   if (waitpid(pid, &status, 0) < 0) {
      printf("CRASH: %s (waitpid failed)\n", c.name);
      g_warn++;
      return;
   }
   if (WIFSIGNALED(status)) {
      printf("CRASH: %s (signal %d)\n", c.name, WTERMSIG(status));
      g_warn++;
      return;
   }
   if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
      g_ok++;
   else
      g_warn++;
}

static double pseudo_ub(double w, double cap)
{
   return ceil(cap / w) + 3.0;
}

static unsigned lcg(unsigned* s)
{
   *s = *s * 1103515245u + 12345u;
   return *s;
}

/** Pattern: i%3 -> binary / bounded (down) / unbounded (up). Seed = first 3 vars. */
static void build_pattern_case(MixedCase& c, const char* name, int n, double cap,
      unsigned rng_seed, int n_seed_vars)
{
   int i, j, in_seed;
   unsigned s = rng_seed;

   memset(&c, 0, sizeof(c));
   c.name = name;
   c.n = n;
   c.cap = cap;
   if (n_seed_vars < 1)
      n_seed_vars = 1;
   if (n_seed_vars > n - 1)
      n_seed_vars = n - 1;
   if (n_seed_vars > 4)
      n_seed_vars = 4;

   for (i = 0; i < n; i++) {
      unsigned wmax = (n <= 5) ? 8u : (n <= 10) ? 6u : 5u;
      c.w[i] = 2.0 + (double)(lcg(&s) % wmax);
      c.p[i] = 0.0;
      c.isuseub[i] = 0;
      switch (i % 3) {
      case 0:
         c.u[i] = 1.0;
         break;
      case 1:
         c.u[i] = 2.0 + (double)(lcg(&s) % 4u);
         c.isuseub[i] = 1;
         break;
      default:
         c.u[i] = pseudo_ub(c.w[i], cap);
         break;
      }
   }

   c.n_seed = n_seed_vars;
   for (i = 0; i < n_seed_vars; i++) {
      c.seed[i] = i;
      c.p[i] = 1.0;
   }

   c.n_ord = 0;
   for (i = 0; i < n; i++) {
      in_seed = 0;
      for (j = 0; j < n_seed_vars; j++) {
         if (c.seed[j] == i) {
            in_seed = 1;
            break;
         }
      }
      if (!in_seed)
         c.order[c.n_ord++] = i;
   }
}

static void make_random_case(MixedCase& c, int seed_id, int n, double cap0)
{
   static char name_buf[40];
   snprintf(name_buf, sizeof(name_buf), "rand%d_%d", n, seed_id);
   build_pattern_case(c, name_buf, n, cap0 + (seed_id % 13), (unsigned)seed_id * 97u + 13u,
         (n >= 8) ? 3 : 2);
}

static void run_hand_small(void)
{
   static const MixedCase hand[] = {
      {"mix_bin_bd_ub",
       5, 15.0,
       {1, 1, 0, 0, 0},
       {4, 3, 5, 2, 6},
       {1, 1, 4, 1, 50},
       {0, 0, 0, 1, 0},
       {0, 1}, 2,
       {2, 3, 4}, 3},
      {"mix_ub_lift_first",
       4, 12.0,
       {1, 0, 0, 0},
       {3, 5, 4, 2},
       {1, 30, 3, 1},
       {0, 0, 1, 0},
       {0}, 1,
       {1, 2, 3}, 3},
      {"mix_heavy_ub",
       5, 25.0,
       {1, 1, 0, 0, 0},
       {6, 5, 8, 3, 4},
       {1, 2, 60, 1, 4},
       {0, 1, 0, 0, 1},
       {0, 1}, 2,
       {2, 3, 4}, 3},
   };
   int c;
   printf("--- Hand small (n<=5, build=%s) ---\n", BUILD_TAG);
   for (c = 0; c < (int)(sizeof(hand) / sizeof(hand[0])); c++)
      test_case(hand[c]);
}

static void run_hand_medium(void)
{
   MixedCase c;

   printf("--- Hand medium (n=8..15, build=%s) ---\n", BUILD_TAG);

   memset(&c, 0, sizeof(c));
   c.name = "mix_n8";
   c.n = 8;
   c.cap = 32.0;
   {
      const double w[] = {4, 3, 5, 6, 4, 7, 3, 5};
      const double u[] = {1, 1, 4, 45, 1, 50, 3, 40};
      const int iu[] = {0, 0, 1, 0, 0, 0, 1, 0};
      int i;
      for (i = 0; i < 8; i++) {
         c.w[i] = w[i];
         c.u[i] = u[i];
         c.isuseub[i] = iu[i];
         c.p[i] = (i < 3) ? 1.0 : 0.0;
      }
      c.seed[0] = 0;
      c.seed[1] = 1;
      c.seed[2] = 2;
      c.n_seed = 3;
      c.order[0] = 3;
      c.order[1] = 4;
      c.order[2] = 5;
      c.order[3] = 6;
      c.order[4] = 7;
      c.n_ord = 5;
   }
   test_case(c);

   build_pattern_case(c, "mix_n10_pat", 10, 38.0, 4242u, 3);
   test_case(c);

   build_pattern_case(c, "mix_n12_pat", 12, 52.0, 7777u, 3);
   test_case(c);

   build_pattern_case(c, "mix_n14_pat", 14, 62.0, 14141u, 3);
   test_case(c);

   memset(&c, 0, sizeof(c));
   c.name = "mix_n15_dense";
   c.n = 15;
   c.cap = 72.0;
   {
      int i;
      for (i = 0; i < 15; i++) {
         c.w[i] = 3.0 + (double)(i % 7);
         c.p[i] = (i < 4) ? 1.0 : 0.0;
         switch (i % 3) {
         case 0:
            c.u[i] = 1.0;
            c.isuseub[i] = 0;
            break;
         case 1:
            c.u[i] = 2.0 + (double)(i % 4);
            c.isuseub[i] = 1;
            break;
         default:
            c.u[i] = pseudo_ub(c.w[i], c.cap);
            c.isuseub[i] = 0;
            break;
         }
      }
      for (i = 0; i < 4; i++)
         c.seed[i] = i;
      c.n_seed = 4;
      c.n_ord = 0;
      for (i = 4; i < 15; i++)
         c.order[c.n_ord++] = i;
   }
   test_case(c);

   build_pattern_case(c, "mix_n20_pat", 20, 88.0, 20261u, 4);
   test_case(c);

   build_pattern_case(c, "mix_n18_pat", 18, 78.0, 18181u, 4);
   test_case(c);
}

static void run_random_batch(const char* label, int n, int seed_lo, int seed_hi, double cap0)
{
   int s;
   printf("--- %s (n=%d, seeds %d-%d, build=%s) ---\n", label, n, seed_lo, seed_hi, BUILD_TAG);
   for (s = seed_lo; s <= seed_hi; s++) {
      MixedCase c;
      make_random_case(c, s, n, cap0);
      test_case(c);
   }
}

static int run_main() {
      setvbuf(stdout, NULL, _IONBF, 0);
      printf("======== Mixed-variable lifting benchmark ========\n");
      run_hand_small();
      printf("\n");
      run_hand_medium();
      printf("\n");
      run_random_batch("Random small", 5, 1, 15, 16.0);
      printf("\n");
      run_random_batch("Random medium", 8, 1, 20, 28.0);
      printf("\n");
      run_random_batch("Random large", 12, 1, 20, 45.0);
      printf("\n");
      run_random_batch("Random xlarge", 15, 1, 15, 60.0);
      printf("\n");
      run_random_batch("Random xlarge2", 18, 1, 15, 70.0);
      printf("\n");
      run_random_batch("Random xxlarge", 20, 1, 10, 82.0);
      printf("\n======== build=%s ok=%d warnings=%d skipped=%d ========\n",
            BUILD_TAG, g_ok, g_warn, g_skip);
   if (g_warn > 0)
      printf("Mixed suite warnings: %d\n", g_warn);
   return 0;
}

} // namespace mixed_tests

int main(int argc, char** argv)
{
   setvbuf(stdout, NULL, _IONBF, 0);
   int run_mixed = 0;
   for (int i = 1; i < argc; i++) {
      if (strcmp(argv[i], "--all") == 0)
         run_mixed = 1;
   }

   int fail = core_tests::run_main();
   printf("\n");
   fail += geq_tests::run_main();
   if (run_mixed) {
      printf("\n");
      mixed_tests::run_main();
   } else {
      printf("\n(Use ./test_dplifting --all or `make test-all` for mixed suite.)\n");
   }
   /* Exit code = total failure count (core + geq), capped at 255 for POSIX. */
   if (fail <= 0)
      return 0;
   return fail > 255 ? 255 : fail;
}
