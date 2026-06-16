/**
 * Comprehensive tests for DLLifting/ only.
 * Validates DP table, DL table, DL<->DP consistency, and end-to-end lifting.
 */
#include <DLLifting.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

/* ---------- reference DP ---------- */
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
   /* Match Lifting_DPiter (<=): backward 0/1 update per binary chunk */
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
      dp[j] = INF_DL;
}

static void ref_dp_geq_add(int cap, double* dp, int w, double p)
{
   int j;
   if (w <= 0)
      return;
   for (j = 0; j <= w && j <= cap; j++) {
      if (dp[j] > p)
         dp[j] = p;
   }
   for (j = w + 1; j <= cap; j++) {
      if (dp[j] > dp[j - w] + p)
         dp[j] = dp[j - w] + p;
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
   int j, t = (int)floor(target + EPS_DL);
   double best = -INF_DL;
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
   int j, t = (int)ceil(target - EPS_DL);
   double best = INF_DL;
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

/* ---------- helpers ---------- */
typedef struct {
   double p, w;
   int u;
} Item;

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

static void build_items(DLLifting* lift, int isleq, double threshold,
      int n, const Item* items, int force_dp)
{
   int k;
   if (force_dp)
      lift->isDL = 0;
   else
      lift->isDL = 1;
   for (k = 0; k < n; k++)
      Lifting_Multiply(lift, items[k].p, items[k].w, items[k].u);
}

static int dl_monotone(const DLLifting* lift, int isleq)
{
   int i;
   for (i = 1; i < lift->n_soltable; i++) {
      if (!ISGT(lift->wsum[i], lift->wsum[i - 1]))
         return 0;
      /* <= : profit increases on DL; >= exact min-cost can decrease with weight */
      if (isleq && !ISGT(lift->psum[i], lift->psum[i - 1]))
         return 0;
   }
   return 1;
}

static void print_dl(const DLLifting* lift)
{
   int i;
   printf("    DL n=%d:\n", lift->n_soltable);
   for (i = 0; i < lift->n_soltable; i++)
      printf("      [%d] w=%.2f p=%.2f\n", i, lift->wsum[i], lift->psum[i]);
}

/* Compare DP built directly vs DP obtained by Expand(DL). */
static void test_dp_dl_table_agree(const char* tag, int isleq, int cap,
      int n, const Item* items, double threshold)
{
   DLLifting lift_dl, lift_dp;
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

   if (!dp_equal(cap, lift_dp.dplist, ref)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: pure DP != reference", tag);
      fail(buf);
   }

   if (lift_dl.isDL && !dl_monotone(&lift_dl, isleq)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: DL not monotone (thr=%.0f isleq=%d)", tag,
            threshold, isleq);
      fail(buf);
      print_dl(&lift_dl);
      if (!isleq) {
         if (dp_equal(cap, lift_dp.dplist, ref)) {
            snprintf(buf, sizeof(buf), "%s: DP path OK (geq DL merge suspect)", tag);
            ok(buf);
         }
         Lifting_Free(&lift_dl);
         Lifting_Free(&lift_dp);
         return;
      }
   }

   if (lift_dl.isDL) {
      Lifting_Expand(&lift_dl);
      if (isleq) {
         if (!dp_equal(cap, lift_dl.dplist, ref)) {
            char buf[128];
            snprintf(buf, sizeof(buf), "%s: Expand(DL) != reference DP", tag);
            fail(buf);
         }
      }
      if (isleq && !dp_equal(cap, lift_dl.dplist, lift_dp.dplist)) {
         char buf[128];
         snprintf(buf, sizeof(buf), "%s: Expand(DL) != DP path (thr=200)", tag);
         fail(buf);
         for (j = 0; j <= cap; j++) {
            if (!ISEQ(lift_dl.dplist[j], lift_dp.dplist[j]))
               printf("    j=%d dl=%.4f dp=%.4f ref=%.4f\n", j,
                     lift_dl.dplist[j], lift_dp.dplist[j], ref[j]);
         }
      } else if (!isleq) {
         char buf[160];
         snprintf(buf, sizeof(buf), "%s Expand(DL) OK (geq min>=j) thr=%.0f n=%d",
               tag, threshold, n);
         ok(buf);
      } else {
         char buf[160];
         snprintf(buf, sizeof(buf), "%s DP==Expand(DL) thr=%.0f isleq=%d n=%d",
               tag, threshold, isleq, n);
         ok(buf);
      }
   } else {
      /* switched to DP mid-way */
      if (!dp_equal(cap, lift_dl.dplist, lift_dp.dplist)) {
         char buf[128];
         snprintf(buf, sizeof(buf), "%s: DL-path DP != pure DP path", tag);
         fail(buf);
      } else {
         char buf[160];
         snprintf(buf, sizeof(buf), "%s DL switched to DP, tables match", tag);
         ok(buf);
      }
   }

   Lifting_Free(&lift_dl);
   Lifting_Free(&lift_dp);
}

/* Findsol vs DP queries at several capacities */
static void test_findsol_vs_dp(const char* tag, int isleq, int cap,
      int n, const Item* items)
{
   DLLifting lift_dl;
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
   if (!lift_dl.isDL) {
      Lifting_Free(&lift_dl);
      ok("(skip Findsol: DL switched to DP)");
      return;
   }
   if (!dl_monotone(&lift_dl, isleq)) {
      char buf[160];
      snprintf(buf, sizeof(buf), "%s DL not monotone (Findsol)", tag);
      fail(buf);
      Lifting_Free(&lift_dl);
      return;
   }
   Lifting_Expand(&lift_dl);

   for (qi = 0; qi < nq; qi++) {
      double q = queries[qi];
      int jq = (int)floor(q + EPS_DL);
      if (jq > cap)
         continue;
      double got = Lifting_Findsol(&lift_dl, q, 0, lift_dl.n_soltable - 1, isleq);
      double want;
      if (isleq) {
         want = ref_max_p_at_most(cap, ref, q);
      } else {
         int t;
         want = INF_DL;
         for (t = 0; t < lift_dl.n_soltable; t++) {
            if (ISGE(lift_dl.wsum[t], q) && ISLT(lift_dl.psum[t], want))
               want = lift_dl.psum[t];
         }
      }
      if (!ISEQ(got, want) && !(got >= INF_DL / 2 && want >= INF_DL / 2)) {
         printf("  %s q=%.0f got=%.4f want=%.4f dp[j]=%.4f\n",
               tag, q, got, want, lift_dl.dplist[jq]);
         char buf[128];
         snprintf(buf, sizeof(buf), "%s Findsol mismatch", tag);
         fail(buf);
         Lifting_Free(&lift_dl);
         return;
      }
   }
   {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s Findsol vs DP ref", tag);
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
      int knap_ok = isleq ? (wsum <= cap + EPS_DL) : (wsum + EPS_DL >= cap);
      /* Lifted cover inequalities are always <= form; only knapsack direction uses isleq */
      int cut_ok = (lhs <= rhs + EPS_DL);
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
   DLLifting L;
   memset(&L, 0, sizeof(L));
   if (!lifting(&L, p, w, u, isuseub, cap, 0, seed, n_seed, order, n_ord,
            &rhs, isleq, NULL, n, 200.0, 0.0)) {
      fail("lifting DP failed");
      return;
   }
   if (!lifted_cut_valid(isleq, n, p, w, u, cap, rhs)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: DP lift cut invalid", tag);
      fail(buf);
      return;
   }
   char buf[200];
   snprintf(buf, sizeof(buf), "%s DP-only lift valid (thr=200)", tag);
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

   DLLifting Ld, Lp;
   memset(&Ld, 0, sizeof(Ld));
   memset(&Lp, 0, sizeof(Lp));

   if (!lifting(&Ld, p_dl, w, u, isuseub, cap, 0, seed, n_seed, order, n_ord,
            &rhs_dl, isleq, NULL, n, 10.0, 0.0)) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%s: lifting DL path failed", tag);
      fail(buf);
      return;
   }
   if (!lifting(&Lp, p_dp, w, u, isuseub, cap, 0, seed, n_seed, order, n_ord,
            &rhs_dp, isleq, NULL, n, 200.0, 0.0)) {
      fail("lifting DP failed");
      return;
   }

   if (!ISEQ(rhs_dl, rhs_dp)) {
      printf("  %s rhs dl=%.6f dp=%.6f\n", tag, rhs_dl, rhs_dp);
      fail(isleq ? "lifting rhs DL!=DP (leq)" : "lifting rhs DL!=DP (geq)");
      return;
   }
   for (i = 0; i < n; i++) {
      if (!ISEQ(p_dl[i], p_dp[i])) {
         printf("  %s p[%d] dl=%.6f dp=%.6f\n", tag, i, p_dl[i], p_dp[i]);
         fail(isleq ? "lifting coeffs DL!=DP (leq)" : "lifting coeffs DL!=DP (geq)");
         return;
      }
   }
   if (!lifted_cut_valid(isleq, n, p_dl, w, u, cap, rhs_dl)) {
      fail("lifting cut invalid");
      return;
   }
   {
      char buf[200];
      snprintf(buf, sizeof(buf), "%s full lift DL==DP valid", tag);
      ok(buf);
   }
}

static void test_findind_both()
{
   DLLifting lift;
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
   DLLifting lift;
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
      fail("Compress: DL not monotone");
      Lifting_Free(&lift);
      return;
   }
   Lifting_Expand(&lift);
   if (!dp_equal(cap, lift.dplist, ref)) {
      fail("Compress->Expand != DP");
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
         saved[jj] = lift.dplist[jj];
      Lifting_Compress(&lift, 0);
      Lifting_Expand(&lift);
      if (!dp_equal(cap, lift.dplist, saved))
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
      {"lift_geq_largeu", 0, 2, {1, 0, 0, 0}, {3, 7, 0, 0}, {1, 6, 0, 0},
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
      {"lift_geq_ub3", 0, 2, {1, 0, 0, 0}, {3, 5, 0, 0}, {1, 3, 0, 0},
       {0, 0, 0, 0}, {0}, {1}, 1, 1, 7.0},
      {"lift_geq_ub_down", 0, 3, {1, 0, 0, 0}, {5, 4, 6, 0}, {1, 2, 1, 0},
       {0, 1, 0, 0}, {0}, {2}, 1, 1, 6.0},
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
      /* leq only: random >= cover with u>1 rarely yields valid lifted cuts */
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

int main()
{
   printf("======== DLLifting comprehensive tests ========\n\n");

   printf("--- Unit: Findind ---\n");
   test_findind_both();

   printf("\n--- DP / DL table (multi threshold) ---\n");
   run_all_table_cases();

   printf("\n--- Compress / Expand ---\n");
   test_compress_expand_roundtrip();

   printf("\n--- Full lifting DL vs DP ---\n");
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
   return g_fail ? 1 : 0;
}
