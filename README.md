# DLLifting

**DLLifting** is a standalone C/C++ library for **DL/DP hybrid coefficient lifting** on knapsack cover inequalities. It supports both `<=` and `>=` knapsack rows, optional capacity reduction (**DP_R** / **DL_R**), and can be embedded in MIP solvers via cut callbacks or custom separators.

## Scope

- Knapsack set with general upper bounds (bounded and unbounded items)
- Sequential up/down lifting with hybrid **DL** / **DP** subproblem solvers
- Forced modes **DL**, **DP**, or threshold-based **AUTO** switching
- Optional reduction variants **DL_R** / **DP_R** (compile-time `REDUCTION=1`)
- Stable C ABI (`dllifting_lift_cover`) and full C++ API (`lifting`)
- No dependency on external optimization libraries

## Requirements

- C++11 compiler (GCC, Clang, or MSVC)
- GNU Make

## Build

| Target | Command | What it does |
| ------ | ------- | ------------ |
| Library only | `make` | Compile sources and produce the shared library `libdllifting.so`. |
| Unit tests | `make test` | Build and run the unified test suite (`tests/test_dllifting.cpp`). |
| Example program | `make example` | Build the `example` binary. Execute `./example` afterward. |

Compile-time reduction support is enabled by default (`REDUCTION=1`). To disable:

```bash
make REDUCTION=0
make test REDUCTION=0
```

Link your own program against the shared library:

```bash
g++ -O2 my_app.cpp -Iinclude -L. -ldllifting -Wl,-rpath,'$ORIGIN' -lm -o my_app
```

Install headers and the shared library (optional):

```bash
make install PREFIX=$HOME/.local
```

## API

All public symbols are declared in `include/DLLifting.h`.

### Lifting modes

| Constant | Effect |
| -------- | ------ |
| `DLLIFTING_MODE_AUTO` | Threshold-based DL↔DP switch on bounded items (default) |
| `DLLIFTING_MODE_DL` | Force DL table; no switching |
| `DLLIFTING_MODE_DP` | Force DP table; no switching |

### `lifting` (C++)

Full sequential lifting with prescribed seed and lifting order. See `include/DLLifting.h`.

### `dllifting_lift_cover` (C)

Stable C wrapper for solver callbacks; see `include/DLLifting.h`.

**Return value:** `DLLIFTING_OK` (0) on success; negative error code otherwise.

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

**Input** (0-based indices in code: \(x_1\!\mapsto\!0,\ldots,x_5\!\mapsto\!4\))


| Field | Value |
| ----- | ----- |
| `n` | 5 |
| `weight` | `[8, 5, 4, 3, 5]` |
| `coef` (seed) | `[2, 1, 0, 0, 0]` |
| `ub` | `[2, 3, 6, 5, 1]` |
| `isuseub` | `[0, 0, 0, 0, 1]` — \(N_u=\{x_5\}\): down-lift |
| `capacity` | `18` (\(b^2\)) |
| `is_subcap` | `1` |
| `seed` | `[0, 1]` (\(x_1,x_2\)) |
| `lifting_order` | `[2, 3, 4]` (\(x_3,x_4,x_5\)) |
| `isdl_mode` | `DLLIFTING_MODE_DP` (integer weights) |

**Output**


| Field | Value |
| ----- | ----- |
| `coef` (lifted) | `[2, 1, 1, 0.5, 1.5]` |
| `rhs` | `5.5` |

Sequential lifting on \(\{x_3,x_4,x_5\}\): \(\alpha_3=1\), \(\alpha_4=\tfrac12\), \(\alpha_5=\tfrac32\).

```cpp
#include <DLLifting.h>
#include <cstdio>

int main() {
   const int n = 5;
   double p[] = {2.0, 1.0, 0.0, 0.0, 0.0};
   double w[] = {8.0, 5.0, 4.0, 3.0, 5.0};
   double u[] = {2.0, 3.0, 6.0, 5.0, 1.0};
   int isuseub[] = {0, 0, 0, 0, 1};   /* N_u = {x5}: down-lift */
   int seed[] = {0, 1};               /* C = {x1, x2} */
   int order[] = {2, 3, 4};           /* lift x3, x4, x5 */
   double rhs = 0.0;

   DLLifting lift = {};
   if (!lifting(&lift, p, w, u, isuseub, 18.0, 1,
            seed, 2, order, 3, &rhs, 1, nullptr, n,
            10.0, 0.0, DLLIFTING_MODE_DP))
      return 1;

   for (int i = 0; i < n; i++) {
      if (p[i] > EPS_DL)
         std::printf("%.4f*x_%d + ", p[i], i + 1);
   }
   std::printf("<= %.4f  (%.4f s)\n", rhs, lift.duration);
   return 0;
}
```

The final cut is \(2x_1+x_2+x_3+\tfrac12 x_4+\tfrac32 x_5\le \tfrac{11}{2}\).

Run this example:

```bash
make
make example && ./example
```

Expected cut: \(2x_1+x_2+x_3+\tfrac12 x_4+\tfrac32 x_5\le \tfrac{11}{2}\).

## Project layout

```
dllifting/
├── include/DLLifting.h   # public header (C++ + C ABI)
├── src/DLLifting.cpp     # lifting + C wrapper (single translation unit)
├── libdllifting.so       # built by make (shared library)
├── examples/example.cpp  # README 5-var instance
├── tests/test_dllifting.cpp
├── Makefile
└── LICENSE
```

## License

MIT License — see [LICENSE](LICENSE). Copyright (c) 2026 Xintong Wang and contributors.

## Citation

If you use DLLifting in research, please cite the associated knapsack lifting / separation work from your publication and link to this repository.
