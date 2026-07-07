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
| Library only | `make` | Compile sources and produce the shared library `libdllifting.so`. Does not build or run tests or examples. |
| Unit tests | `make test` | Build `libdllifting.so` and run `test_dllifting` + `test_isgeq`. Fails if any test fails. |
| Extended tests | `make test-all` | Run `make test` plus mixed-variable regression (`test_mixed_vars`, `test_mixed_vars_r`). |
| Example program | `make example` | Build the `example` binary (see `examples/example.cpp`). Execute `./example` afterward. |

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

Full sequential lifting with prescribed seed and lifting order. See `include/DLLifting.h` for the complete signature and `Lifting` / `DLLifting` workspace type.

### `dllifting_lift_cover` (C)

Stable C wrapper for solver callbacks:

```c
int dllifting_lift_cover(
      int n, double* coef,
      const double* weight, const double* ub, const int* use_ub,
      double cap, int is_subcap,
      const int* seed, int n_seed,
      const int* lifting_order, int n_order,
      double* rhs,
      int is_leq, double threshold, int isdl_mode,
      const double* x_frac);
```

**Return value:** `DLLIFTING_OK` (0) on success; negative error code otherwise.

### Example

```cpp
#include <DLLifting.h>
#include <cstdio>

int main() {
   const int n = 3;
   double p[] = {1.0, 1.0, 0.0};
   double w[] = {2.0, 3.0, 4.0};
   double u[] = {1.0, 1.0, 1.0};
   int isuseub[] = {0, 0, 0};
   int seed[] = {0, 1};
   int order[] = {2};
   double rhs = 0.0;

   DLLifting lift = {};
   if (!lifting(&lift, p, w, u, isuseub, 4.0, 0, seed, 2, order, 1,
            &rhs, 1, nullptr, n, 10.0, 0.0, DLLIFTING_MODE_AUTO))
      return 1;

   for (int i = 0; i < n; i++)
      if (p[i] > EPS_DL)
         std::printf("%.4f*x_%d + ", p[i], i + 1);
   std::printf("<= %.4f  (%.4f s)\n", rhs, lift.duration);
   return 0;
}
```

More demos: `examples/example.cpp` (`make example && ./example`).

Extended API notes, migration guide, and solver integration: see `docs/`.

## Project layout

```
dllifting/
├── include/DLLifting.h   # public header (C++ + C ABI)
├── src/DLLifting.cpp     # lifting implementation
├── src/dllifting_c.cpp   # C wrapper
├── libdllifting.so       # built by make (shared library)
├── examples/             # usage examples
├── tests/                # unit and regression tests
├── docs/                 # API, migration, integration notes
├── Makefile
└── LICENSE
```

## License

MIT License — see [LICENSE](LICENSE). Based on work by Igor Vasilyev and contributors.

## Citation

If you use DLLifting in research, please cite the associated knapsack lifting / separation work from your publication and link to this repository.
