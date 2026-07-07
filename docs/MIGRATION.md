# Migrating to 1.1.0

Version **1.1.0** extends the public API with an explicit **lifting mode** parameter.
If you upgrade from **1.0.0**, update all call sites as below.

## C++: `lifting()`

**Before (1.0.0):**

```cpp
lifting(&lift, p, w, u, isuseub, cap, 0,
        seed, n_seed, order, n_order, &rhs,
        is_leq, x, n, threshold, duration);
```

**After (1.1.0):**

```cpp
lifting(&lift, p, w, u, isuseub, cap, 0,
        seed, n_seed, order, n_order, &rhs,
        is_leq, x, n, threshold, duration,
        DLLIFTING_MODE_AUTO);   // preserve old threshold switching
```

To **force** a table type (no DL↔DP switching on bounded variables):

```cpp
lifting(..., threshold, duration, DLLIFTING_MODE_DL);
lifting(..., threshold, duration, DLLIFTING_MODE_DP);
```

## C: `dllifting_lift_cover()`

**Before (1.0.0):**

```c
dllifting_lift_cover(n, coef, weight, ub, use_ub,
    cap, is_subcap, seed, n_seed, order, n_order,
    &rhs, is_leq, threshold, x_frac);
```

**After (1.1.0):**

```c
dllifting_lift_cover(n, coef, weight, ub, use_ub,
    cap, is_subcap, seed, n_seed, order, n_order,
    &rhs, is_leq, threshold, DLLIFTING_MODE_AUTO, x_frac);
```

## Build / install

After upgrading, rebuild and install:

```bash
make clean && make
make install PREFIX=$HOME/.local
```

Link with `-I$PREFIX/include -L$PREFIX/lib -ldllifting`.

## Behaviour changes (non-breaking if you pass `AUTO`)

- Forced modes affect **both** bounded and unbounded subproblems.
- `DLLIFTING_MODE_AUTO` matches 1.0.0 threshold switching on bounded items only;
  unbounded behaviour is unchanged from 1.0.0.

## Checklist

- [ ] Add `isdl_mode` to every `lifting()` call
- [ ] Add `isdl_mode` to every `dllifting_lift_cover()` call (before `x_frac`)
- [ ] Rebuild and run your regression tests
- [ ] For `>=` knapsacks, prefer `DLLIFTING_MODE_DP` (see [TESTING.md](TESTING.md))
