# DLLifting API reference

## Mode constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `DLLIFTING_MODE_AUTO` | `-1` | Legacy behaviour: initial table from `threshold`; bounded items may switch DL↔DP when `threshold` crosses 100. |
| `DLLIFTING_MODE_DP` | `0` | Force DP table (`isDL = 0`) for the whole lift; no switching. |
| `DLLIFTING_MODE_DL` | `1` | Force DL table (`isDL = 1`) for the whole lift; no switching. |

Defined in `include/DLLifting.h` (C++) and `include/dllifting_c.h` (C).

## C++ entry point

```cpp
int lifting(
    DLLifting* lift,
    DTptype* p, DTwtype* w, DTutype* u, int* isuseub,
    DTctype cap, int isSubCap,
    int* seed, int n_seed,
    int* liftingorder, int n_liftingorder,
    double* rhs,
    int isLeq, double* x, int n,
    double threshold,
    double duration,
    int isdl_mode);
```

### Parameters (selected)

| Parameter | Description |
|-----------|-------------|
| `p`, `w`, `u` | Cut coefficients, knapsack weights, variable upper bounds (length `n`). |
| `isuseub` | Per variable: `1` = down-lift (fix at UB), `0` = up-lift. |
| `cap` | Knapsack capacity `b`. |
| `isSubCap` | `1` if `cap` is residual subcapacity after fixing the seed; else global rhs. |
| `seed`, `n_seed` | Variables in the initial cover. |
| `liftingorder`, `n_liftingorder` | Remaining variables in sequential lifting order. |
| `rhs` | Output: lifted inequality right-hand side. |
| `isLeq` | `1` for `sum w_i x_i <= b`; `0` for `>= b`. |
| `x` | Optional fractional point (may be `nullptr`). |
| `threshold` | Used only when `isdl_mode == DLLIFTING_MODE_AUTO` (see below). |
| `duration` | Reserved; currently unused (pass `0.0`). |
| `isdl_mode` | Table mode; see constants above. |

### Return value

`1` on success (coefficients written into `p`, `rhs` updated); `0` on allocation or internal failure.

## C entry point

```c
int dllifting_lift_cover(
    int n,
    double* coef,
    const double* weight,
    const double* ub,
    const int* use_ub,
    double cap,
    int is_subcap,
    const int* seed,
    int n_seed,
    const int* lifting_order,
    int n_order,
    double* rhs,
    int is_leq,
    double threshold,
    int isdl_mode,
    const double* x_frac);
```

Return codes: `DLLIFTING_OK`, `DLLIFTING_ERR_ALLOC`, `DLLIFTING_ERR_ARGS`.

## Choosing a mode

### `DLLIFTING_MODE_AUTO` (default for integrators)

- Initial table: DL if `threshold <= 100`, else DP.
- During **bounded** binary-split updates only, the implementation may switch
  DL↔DP when `threshold` crosses 100.
- **Unbounded** items (`w * (u+1) > cap`) always follow the current `isDL` flag;
  with typical large thresholds they stay on DL, so DL and DP timings can coincide
  unless you force a mode.

### `DLLIFTING_MODE_DL` / `DLLIFTING_MODE_DP`

Use these for benchmarks or when you want a fixed algorithm regardless of `threshold`:

| Item type | Forced DL | Forced DP |
|-----------|-----------|-----------|
| Unbounded | `Lifting_Mergesortinf` only | `Compress` → `Mergesortinf` → `Expand` |
| Bounded | `Lifting_Mergesort` | `Lifting_DPiter` |
| Up / down lift | Query `psum` | Query `dplist` |

`threshold` is ignored for table choice when mode is `_DL` or `_DP`.

### Recommended values

| Use case | `isdl_mode` | `threshold` |
|----------|-------------|-------------|
| Production cut separator (legacy) | `DLLIFTING_MODE_AUTO` | `10` (DL) or `200` (DP) |
| Benchmark pure DL | `DLLIFTING_MODE_DL` | any (ignored) |
| Benchmark pure DP | `DLLIFTING_MODE_DP` | any (ignored) |
| `>=` knapsack (known DL issues) | `DLLIFTING_MODE_DP` | — |

## `DLLifting` context fields (read-only after lift)

| Field | Meaning |
|-------|---------|
| `duration` | CPU time of the last `lifting()` call (seconds). |
| `isDL` | Final table type (`1` = DL, `0` = DP). |
| `force_mode` | Mode used (`AUTO` / `_DP` / `_DL`). |
| `threshold` | Threshold stored at allocation. |

## Examples

- Combined demos (C++ + C ABI + forced mode): `examples/example.cpp` (`make example && ./example`)
