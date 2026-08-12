# DLLifting

**DLLifting** (v1.3.0) is a standalone C/C++ library for **DL/DP hybrid coefficient lifting** on general knapsack sets. Optional capacity reduction (**+R**) is independent of the DL/DP backend. Prefer the C++ API `lifting()` or the C ABI `dllifting_lift_cover()`.

## Scope

- General upper bounds (bounded and unbounded variables)
- Sequential up/down lifting with **DL** or **DP** subproblem tables
- Three DL/DP controls: hard force, τ-threshold policy, or automatic feature map
- Optional **+R** (compile with `REDUCTION=1`; runtime ON/OFF/AUTO)
- No external solver dependency

**Maturity:** \(\sum w x \le b\) (`is_leq = 1`) is the production path.  
\(\sum w x \ge b\) (`is_leq = 0`) is supported and covered by tests.

## Requirements

- C++11 (GCC or Clang; MSVC untested)
- GNU Make

## Build

| Target | Command |
| ------ | ------- |
| Library | `make` → `libdllifting.so` |
| Tests | `make test` |
| Extended tests | `make test-all` |
| Example | `make example && ./example` |

```bash
make && make test
make clean && make REDUCTION=0 && make test REDUCTION=0   # build without +R
```

```bash
g++ -O2 my_app.cpp -Iinclude -L. -ldllifting -Wl,-rpath,'$ORIGIN' -lm -o my_app
make install   # default PREFIX=$HOME/.local
```

## Parameter model (v1.3)

Two **orthogonal** axes. Do not combine `MODE_AUTO` with `threshold` to mean “force DL/DP”.

### Axis 1 — DL / DP (`isdl_mode`)

| Mode | Meaning |
| ---- | ------- |
| `DLLIFTING_MODE_DL` / `_DP` | **Manual hard force.** No mid-lift switch. |
| `DLLIFTING_MODE_THRESHOLD` | **Manual τ policy.** `threshold` = capacity \(\tau\). Initial and mid-lift use \(\tau\) vs \(\bar b^k=\min(b^k,U^k)\): \(\tau>\bar b^k\) → DP, \(\tau<\bar b^k\) → DL. |
| `DLLIFTING_MODE_AUTO` | **Default automatic.** Select from row features \((\rho_w,\beta,\bar u)\). **Ignores `threshold` for mode.** No mid-lift. |

Defaults for AUTO (override via `lift->rho_th` / `beta_th` / `u_bar_th` before the call; `0` → library defaults):

| Symbol | Default | Rule of thumb |
| ------ | ------- | ------------- |
| \(\rho_w=w_{\max}/w_{\min}\) | \(\rho_{\mathrm{th}}=6\) | \(\rho_w\ge\rho_{\mathrm{th}}\) → DL |
| \(\beta=b/\mathrm{mean}(w)\) | \(\beta_{\mathrm{th}}=6\) | \(\beta<\beta_{\mathrm{th}}\) → DL; else if narrow \(w\) → DP |
| \(\bar u=\mathrm{mean}(u)\) | \(\bar u_{\mathrm{th}}=3\) | near-binary / small \(\bar u\) → DL |
| Fuzzy band | 10% of thresholds | favors DL |

Helpers: `dllifting_compute_features`, `dllifting_select_backend`, `dllifting_policy_default`.

### Axis 2 — Reduction +R (`lift->reduction_request`)

Requires compile-time `REDUCTION=1` (default in `make`). Always disabled if any variable is treated as unbounded.

| Value | Meaning |
| ----- | ------- |
| `DLLIFTING_RED_ON` / `_OFF` | **Manual** |
| `DLLIFTING_RED_AUTO` (0, default) | Enable iff \(\bar b^0 > \tau\) (large residual) |

### Role of `threshold`

| Context | Role of `threshold` |
| ------- | ------------------- |
| `MODE_THRESHOLD` | \(\tau\) for initial + mid-lift |
| `RED_AUTO` | \(\tau\) for the \(\bar b^0 > \tau\) test |
| `MODE_AUTO` / `_DL` / `_DP` | **Not** used to choose the backend |

If `threshold <= 0`, \(\tau=\beta_{\mathrm{th}}\cdot\mathrm{mean}(w)\).

### Quick recipes

```cpp
DLLifting lift = {};
// Default production: AUTO backend + AUTO +R
lifting(&lift, ..., /*threshold*/ 0.0, 0.0, DLLIFTING_MODE_AUTO);

// Hard DP / DL (benchmarks)
lifting(&lift, ..., 0.0, 0.0, DLLIFTING_MODE_DP);
lifting(&lift, ..., 0.0, 0.0, DLLIFTING_MODE_DL);

// Manual τ with mid-lift (e.g. τ = 200)
lifting(&lift, ..., 200.0, 0.0, DLLIFTING_MODE_THRESHOLD);

// Force +R off
lift.reduction_request = DLLIFTING_RED_OFF;
lifting(&lift, ..., 0.0, 0.0, DLLIFTING_MODE_AUTO);
```

## API

Header: `include/DLLifting.h` (shim: `dllifting_c.h`). Version: `DLLIFTING_VERSION` (`"1.3.0"`).

### `lifting` (C++)

```cpp
int lifting(
      DLLifting* lift,
      double* p, double* w, double* u, int* isuseub,
      double cap, int isSubCap,
      int* seed, int n_seed,
      int* liftingorder, int n_liftingorder,
      double* rhs,
      int isLeq, double* x, int n,
      double threshold, double duration, int isdl_mode);
```

| Parameter | Meaning |
| --------- | ------- |
| `lift` | Workspace. Optionally set `reduction_request`, `rho_th`, `beta_th`, `u_bar_th` before the call. On success `lift->duration` is CPU seconds (struct is cleared/rebuilt inside; policy fields are preserved). |
| `p` | Seed coefficients in → lifted coefficients out |
| `w`, `u` | Weights and upper bounds |
| `isuseub` | `1` = down-lift (fix at UB), `0` = up-lift |
| `cap` / `isSubCap` | Capacity \(b\), or residual if `isSubCap=1` |
| `seed` / `liftingorder` | Cover indices and remaining lift order |
| `rhs` | Lifted RHS (out) |
| `isLeq` | `1`: \(\le\) knapsack; `0`: \(\ge\) (experimental) |
| `x` | Optional fractional point; may be `NULL` |
| `threshold` | \(\tau\) for `MODE_THRESHOLD` and `RED_AUTO`; see table above |
| `duration` | Optional time limit (seconds); `<=0` = unlimited |
| `isdl_mode` | `MODE_AUTO` / `MODE_THRESHOLD` / `MODE_DL` / `MODE_DP` |

**Return:** `1` success, `0` failure.

### `dllifting_lift_cover` (C)

Same semantics; returns `DLLIFTING_OK` (0) or `DLLIFTING_ERR_*`.

## Example

Knapsack set \(\mathcal{X}=\{x\in\mathbb{Z}_+^5: 8x_1+5x_2+4x_3+3x_4+5x_5\le 23,\ \ldots\}\).  
Seed \(2x_1+x_2\le 4\) on a restricted face with residual capacity \(18\); lift order \(\{3,4,5\}\).

```cpp
#include <DLLifting.h>
#include <cstdio>

int main() {
   const int n = 5;
   double p[] = {2.0, 1.0, 0.0, 0.0, 0.0};
   double w[] = {8.0, 5.0, 4.0, 3.0, 5.0};
   double u[] = {2.0, 3.0, 6.0, 5.0, 1.0};
   int isuseub[] = {0, 0, 0, 0, 1};
   int seed[] = {0, 1};
   int order[] = {2, 3, 4};
   double rhs = 0.0;

   DLLifting lift = {};
   if (!lifting(&lift, p, w, u, isuseub, 18.0, 1,
            seed, 2, order, 3, &rhs, 1, nullptr, n,
            /* threshold */ 0.0, /* duration */ 0.0, DLLIFTING_MODE_DP))
      return 1;

   for (int i = 0; i < n; i++)
      if (p[i] > EPS_DL)
         std::printf("%.4f*x_%d + ", p[i], i + 1);
   std::printf("<= %.4f  (%.4f s)\n", rhs, lift.duration);
   return 0;
}
```

Expected: coef `[2, 1, 1, 0.5, 1.5]`, rhs `5.5`.

```bash
make && make example && ./example
```

## Layout

```
DLLifting/
├── include/DLLifting.h
├── include/dllifting_c.h
├── src/DLLifting.cpp
├── examples/example.cpp
├── tests/test_dllifting.cpp
├── Makefile
├── CHANGELOG.md
└── LICENSE
```

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Xintong Wang, Liang Chen, Yu-hong Dai.

## Citation

```
Xintong Wang et al. DLLifting: DL/DP hybrid lifting for general knapsack set.
https://159.226.92.34:8000/wangxintong/dllifting (version 1.3.0).
```
