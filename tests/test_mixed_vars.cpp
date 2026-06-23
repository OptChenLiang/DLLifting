/**
 * Mixed binary / bounded / unbounded lifting benchmarks.
 * Sizes: n in {5,8,10,12,14,15,18,20}; hand + random seeds; DL & DP (N / R).
 * Run: make run-mixed
 */
#include <DLLifting.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#define MAX_MIX 32
#ifdef DLLIFTING_REDUCTION
static const char* BUILD_TAG = "R";
#else
static const char* BUILD_TAG = "N";
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
   return compute_subcap(c) >= -EPS_DL;
}

struct LiftOutcome {
   int ok;
   double rhs;
   double p[MAX_MIX];
   double duration;
   int reduction_active;
};

static LiftOutcome run_lifting_once(const MixedCase& c, int isdl_mode)
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
   DLLifting lift;
   memset(&lift, 0, sizeof(lift));

   clock_t t0 = clock();
   out.ok = lifting(&lift, p, w, u, isuseub, c.cap, 0,
         seed, c.n_seed, order, c.n_ord, &rhs, 1, NULL, c.n, 0.0, 0.0,
         isdl_mode);
   out.duration = (double)(clock() - t0) / CLOCKS_PER_SEC;
   out.rhs = rhs;
   memcpy(out.p, p, (size_t)c.n * sizeof(double));
#ifdef DLLIFTING_REDUCTION
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
   LiftOutcome dl = run_lifting_once(c, DLLIFTING_MODE_DL);
   LiftOutcome dp = run_lifting_once(c, DLLIFTING_MODE_DP);

   print_outcome(c.name, "DL", c, dl);
   print_outcome(c.name, "DP", c, dp);

   if (!outcome_sane(dl, c.n) || !outcome_sane(dp, c.n)) {
      printf("WARN: %s invalid or huge coefficients\n", c.name);
      g_warn++;
      return;
   }
   if (!coeffs_match(dl, dp, c.n)) {
      printf("WARN: %s DL!=DP rhs=%.4f/%.4f\n", c.name, dl.rhs, dp.rhs);
      g_warn++;
   } else {
      printf("OK: %s DL==DP\n", c.name);
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

int main(void)
{
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
   return 0;
}
