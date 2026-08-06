# DLLifting

**DLLifting** (v1.2.1) is a standalone C/C++ library for **DL/DP hybrid coefficient lifting** on knapsack cover inequalities. It supports optional capacity reduction (**DL-R** / **DP-R**) and can be embedded in MIP solvers via cut callbacks or custom separators.

## Scope

- Knapsack set with general upper bounds (bounded and unbounded items)
- Sequential up/down lifting with hybrid **DL** / **DP** subproblem solvers
- Forced modes **DL**, **DP**, or threshold-based **AUTO** switching
- Optional reduction variants **DL-R** / **DP-R** (compile-time `REDUCTION=1`)
- Preferred API: C++ `lifting()`, C ABI `dllifting_lift_cover()`
- No dependency on external optimization libraries

**Maturity:** `sum w x <= b` (`is_leq = 1`) is the supported / tested path (`make test`).  
`sum w x >= b` (`is_leq = 0`) is **experimental** — extended checks via `make test-all` may still report failures.

## Requirements

- C++11 compiler (GCC or Clang; MSVC untested)
- GNU Make

## Build

| Target | Command | What it does |
| ------ | ------- | ------------ |
| Library | `make` | Build `libdllifting.so` |
| Unit tests | `make test` | Core `<=` **and** geq `>=`; exit code = failure count (e.g. 5 geq fails → exit 5) |
| Extended tests | `make test-all` | Also mixed-variable suite |
| Example | `make example` | Build `./example` (then run it) |

```bash
make
make test
make example && ./example
```

Disable reduction:

```bash
make clean && make REDUCTION=0 && make test REDUCTION=0
```

Link your program:

```bash
g++ -O2 my_app.cpp -Iinclude -L. -ldllifting -Wl,-rpath,'$ORIGIN' -lm -o my_app
```

Install (default `PREFIX=$HOME/.local`):

```bash
make install
# or: make install PREFIX=/usr/local
```

## API

Public header: `include/DLLifting.h` (optional shim: `dllifting_c.h`).  
Version: `DLLIFTING_VERSION` (`"1.2.1"`).

### Lifting modes

| Constant | Effect |
| -------- | ------ |
| `DLLIFTING_MODE_AUTO` | Use `threshold` (`< 100` → DL, `> 100` → DP); may switch on bounded items |
| `DLLIFTING_MODE_DL` | Force DL; `threshold` ignored |
| `DLLIFTING_MODE_DP` | Force DP; `threshold` ignored |

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

| Parameter | Type | Direction | Meaning |
| --------- | ---- | --------- | ------- |
| `lift` | `DLLifting*` | in/out | Workspace; on success `lift->duration` holds CPU seconds |
| `p` | `double*` | in/out | Seed coefficients in, lifted coefficients out (length `n`) |
| `w` | `double*` | in | Knapsack weights \(a_j\) |
| `u` | `double*` | in | Variable upper bounds |
| `isuseub` | `int*` | in | Per variable: `1` = down-lift (fix at UB), `0` = up-lift |
| `cap` | `double` | in | Capacity \(b\), or residual subcapacity if `isSubCap=1` |
| `isSubCap` | `int` | in | `1` if `cap` is residual after fixing \(N_u\); else `0` |
| `seed` | `int*` | in | Indices of seed (cover) variables |
| `n_seed` | `int` | in | Length of `seed` |
| `liftingorder` | `int*` | in | Remaining variables in lift order |
| `n_liftingorder` | `int` | in | Length of `liftingorder` |
| `rhs` | `double*` | out | Lifted inequality right-hand side |
| `isLeq` | `int` | in | `1`: \(\sum w x \le b\); `0`: \(\sum w x \ge b\) (experimental) |
| `x` | `double*` | in | Optional fractional point; may be `NULL` |
| `n` | `int` | in | Number of variables |
| `threshold` | `double` | in | Used only if `isdl_mode == AUTO`: `&lt;100` → DL, `&gt;100` → DP |
| `duration` | `double` | in | Reserved / unused (timing goes to `lift->duration`) |
| `isdl_mode` | `int` | in | `DLLIFTING_MODE_AUTO` / `_DL` / `_DP` |

**Return:** `1` on success, `0` on failure (allocation or internal error).

### `dllifting_lift_cover` (C)

```c
int dllifting_lift_cover(
      int n, double* coef,
      const double* weight, const double* ub, const int* use_ub,
      double cap, int is_subcap,
      const int* seed, int n_seed,
      const int* lifting_order, int n_order,
      double* rhs,
      int is_leq, double threshold, int isdl_mode,
      const double* x_frac);    /* may be NULL */
```

**Return:** `DLLIFTING_OK` (0); `DLLIFTING_ERR_ARGS` / `DLLIFTING_ERR_ALLOC` / `DLLIFTING_ERR_INTERNAL` on failure.

### Example

Consider the knapsack set

\[
\mathcal{X}=\Big\{ \boldsymbol{x}\in\mathbb{Z}_+^5:\;
8x_1+5x_2+4x_3+3x_4+5x_5\le 23,\;
x_1\le 2,\; x_2\le 3,\; x_3\le 6,\; x_4\le 5,\; x_5\le 1 \Big\}.
\]

Let \(C=\{1,2\}\), \(N_0=\{3,4\}\), \(N_u=\{5\}\), and lift in the order \(\{3,4,5\}\).
The seed inequality \(2x_1+x_2\le 4\) is valid on the restricted set

\[
\mathcal{X}(N_0,N_u)=\Big\{ \boldsymbol{x}\in\mathbb{Z}_+^5:\;
8x_1+5x_2\le 18,\;
x_1\le 2,\; x_2\le 3,\; x_3=0,\; x_4=0,\; x_5=1 \Big\},
\]

with \(N_0^2=\{3,4\}\), \(N_u^2=\{5\}\), \(b^2=18\), \(\beta^2=4\).

**Input** (0-based indices: \(x_1\!\mapsto\!0,\ldots,x_5\!\mapsto\!4\))

| Field | Value |
| ----- | ----- |
| `n` | 5 |
| `weight` | `[8, 5, 4, 3, 5]` |
| `coef` (seed) | `[2, 1, 0, 0, 0]` |
| `ub` | `[2, 3, 6, 5, 1]` |
| `isuseub` | `[0, 0, 0, 0, 1]` — \(N_u=\{x_5\}\): down-lift |
| `capacity` | `18` (\(b^2\)) |
| `is_subcap` | `1` |
| `seed` | `[0, 1]` |
| `lifting_order` | `[2, 3, 4]` |
| `threshold` | `10.0` — under `AUTO`: `<100` → DL, `>100` → DP; ignored if mode is `DL`/`DP` |
| `isdl_mode` | `DLLIFTING_MODE_DP` |

**Output:** coef `[2, 1, 1, 0.5, 1.5]`, rhs `5.5`  
(i.e. \(2x_1+x_2+x_3+\tfrac12 x_4+\tfrac32 x_5\le \tfrac{11}{2}\)).

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
            /* threshold */ 10.0, /* duration */ 0.0, DLLIFTING_MODE_DP))
      return 1;

   for (int i = 0; i < n; i++) {
      if (p[i] > EPS_DL)
         std::printf("%.4f*x_%d + ", p[i], i + 1);
   }
   std::printf("<= %.4f  (%.4f s)\n", rhs, lift.duration);
   return 0;
}
```

```bash
make && make example && ./example
```

## Project layout

```
dllifting/
├── include/DLLifting.h      # public header (C++ + C ABI)
├── include/dllifting_c.h    # compatibility include
├── src/DLLifting.cpp        # single translation unit
├── examples/example.cpp
├── tests/test_dllifting.cpp
├── Makefile
├── CHANGELOG.md
├── .gitlab-ci.yml
└── LICENSE
```

## License

MIT License — see [LICENSE](LICENSE). Copyright (c) 2026 Xintong Wang and contributors.

## Citation

If you use DLLifting in research, please cite your related knapsack lifting / separation paper and this repository:

```
Xintong Wang et al. DLLifting: DL/DP hybrid lifting for knapsack cover inequalities.
https://159.226.92.34:8000/wangxintong/dllifting (version 1.2.1).
```
