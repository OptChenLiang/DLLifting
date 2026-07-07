# Integrating DLLifting with MIP solvers

DLLifting is a **library**, not a solver plugin. You call it from a **cut callback** (or a custom separator) when you have:

- a knapsack row `sum_i w_i x_i <= b` or `>= b`,
- a seed cover inequality `sum_{i in S} p_i x_i <= rhs` (always in `<=` form after lifting),
- a lifting order for the remaining variables.

The lifted output is new coefficients in `p[]` and an updated `rhs`.

## Data mapping

| DLLifting | Typical MIP meaning |
|-----------|---------------------|
| `p[i]` | Cut coefficient for variable `i` |
| `w[i]` | Knapsack weight (often `a_i` from row) |
| `u[i]` | Variable upper bound |
| `isuseub[i]` | 1 = down-lift (fix at UB), 0 = up-lift |
| `seed[]` | Variables in the initial cover |
| `liftingorder[]` | Order for sequential lifting |
| `isLeq` | 1 if knapsack is `<= b`, 0 if `>= b` |
| `threshold` | Used when `isdl_mode = DLLIFTING_MODE_AUTO` (e.g. 10 for DL, 200 for DP) |
| `isdl_mode` | `DLLIFTING_MODE_AUTO`, `_DL`, or `_DP` — see [API.md](API.md) |

## C++ (recommended)

```cpp
#include <DLLifting.h>

void add_lifted_cover_cut(/* solver handles */, int n, ...) {
  double rhs = 0;
  DLLifting ctx;
  int ok = lifting(&ctx, p, w, u, isuseub, cap, 0,
                   seed, n_seed, order, n_order,
                   &rhs, is_leq, x_frac, n, 10.0, 0.0,
                   DLLIFTING_MODE_AUTO);
  if (!ok) return;
  // solver->addCut(n, indices, p, '<=', rhs);
}
```

## C API

```c
#include <dllifting_c.h>

double rhs;
int rc = dllifting_lift_cover(n, coef, weight, ub, use_ub,
    cap, 0, seed, n_seed, order, n_order,
    &rhs, is_leq, 10.0, DLLIFTING_MODE_AUTO, x_frac);
```

## IBM CPLEX

Use a **lazy constraint** or **user cut** callback (`CPXsetlazyconstraintcallbackfunc` / `CPXsetusercutcallbackfunc`).

1. In the callback, read the fractional solution `x`.
2. Run your **knapsack separation** to find a violated cover (e.g. `ExactSeparation/knapsep`).
3. Build `seed`, `liftingorder`, `p`, `w`, `u`, `isuseub`.
4. Call `lifting()` or `dllifting_lift_cover()`.
5. Add the cut with `CPXcutcallbackadd` / `CPXXcallbackaddusercuts`.

Link: `libdllifting.so` (or static link of `src/*.cpp`) + solver libraries. No solver header is required inside DLLifting itself.

## Gurobi

Use `GRBsetcallback` with `GRB_CB_MIPSOL` or `GRB_CB_MIPNODE`.

1. `GRBgetcallbackinfo(..., GRB_CB_MIPSOL_SOL, ...)` for `x`.
2. Separate knapsack → lift with DLLifting.
3. `GRBcbcut` or `GRBXaddconstr` with `GRB_LESS_EQUAL` and lifted coefficients.

## FICO Xpress / others

Same pattern: fractional solution → separate → `lifting()` → add linear `<=` cut.

## Linking

Build and install the shared library:

```bash
cd DLLifting
make
make install PREFIX=$HOME/.local
```

Link your separator:

```bash
g++ -O2 my_separator.cpp -I$HOME/.local/include -L$HOME/.local/lib \
    -ldllifting -Wl,-rpath,'$HOME/.local/lib' -lm -o my_separator
```

## Notes

- Cuts are always emitted as **`sum p_i x_i <= rhs`**; for `>=` knapsacks the lifting still produces a valid `<=` cutting plane for the covering set reformulation used here.
- Set `u[i]` to a large value (e.g. `1e20`) for unbounded integers; items with `w*(u+1) > cap` are treated as unbounded internally.
- Enable reduction at build time: `make REDUCTION=1` (default).
