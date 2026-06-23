# DLLifting

[![CI](https://github.com/wangxintong216-Cinty/dllifting/actions/workflows/ci.yml/badge.svg)](https://github.com/wangxintong216-Cinty/dllifting/actions/workflows/ci.yml)

**DLLifting** is a standalone C/C++ library for **DL/DP hybrid coefficient lifting** on knapsack cover inequalities. It supports both `<=` and `>=` knapsack rows and can be embedded in MIP solvers via cut callbacks or custom separators.

## Features

- Hybrid **DL** (threshold &lt; 100) and **DP** (threshold &gt; 100) subproblem solvers
- **Forced mode** (`DLLIFTING_MODE_DL` / `_DP`) for bounded and unbounded items — no threshold switching (v1.1.0)
- Optional **capacity reduction** (`tableleft`) for faster lifting
- C++ API (`lifting()`) and stable **C ABI** (`dllifting_lift_cover()`)
- CMake package config for `find_package(dllifting)`
- No external solver dependencies

## Quick start

### Build and install

```bash
git clone https://github.com/wangxintong216-Cinty/dllifting.git
cd dllifting
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local \
         -DDLLIFTING_REDUCTION=ON \
         -DDLLIFTING_BUILD_EXAMPLES=ON \
         -DDLLIFTING_BUILD_TESTS=ON
cmake --build .
ctest          # optional
cmake --install .
```

### Legacy Makefile (tests only)

```bash
make run          # main test suite
make run-all      # test_dllifting + test_isgeq
make run-mixed    # mixed-variable benchmark
```

## Use in your project

### CMake

```cmake
find_package(dllifting CONFIG REQUIRED)
add_executable(my_cut_gen my_cut_gen.cpp)
target_link_libraries(my_cut_gen PRIVATE dllifting::dllifting)
```

```cpp
#include <DLLifting.h>
```

### Manual link

```bash
g++ -O2 my_app.cpp -I/path/to/dllifting/include -L/path/to/lib -ldllifting -lm -o my_app
```

### C API

```c
#include <dllifting_c.h>
```

See `examples/example_basic.cpp`, `examples/example_c_api.c`, and `examples/example_force_mode.cpp`.

## Public API

| Entry | Language | Description |
|-------|----------|-------------|
| `lifting()` | C++ | Full sequential lifting; see `include/DLLifting.h` |
| `dllifting_lift_cover()` | C | Wrapper in `include/dllifting_c.h` |

**Lifting mode** (last argument / before `x_frac`):

| Constant | Effect |
|----------|--------|
| `DLLIFTING_MODE_AUTO` | Threshold-based DL↔DP switch on bounded items (default for integrators) |
| `DLLIFTING_MODE_DL` | Force DL table; no switching |
| `DLLIFTING_MODE_DP` | Force DP table; no switching |

Full reference: [docs/API.md](docs/API.md). Upgrading from 1.0.0: [docs/MIGRATION.md](docs/MIGRATION.md). Release notes: [CHANGELOG.md](CHANGELOG.md).

Solver integration (CPLEX, Gurobi, …): [docs/INTEGRATION_SOLVERS.md](docs/INTEGRATION_SOLVERS.md).

## Project layout

```
dllifting/
├── include/             # Public headers
├── src/                 # Library implementation
├── examples/            # Minimal usage examples
├── tests/               # Unit and regression tests
├── docs/                # API, migration, integration, testing
├── cmake/               # CMake package config template
└── CHANGELOG.md         # Release history
```

## CMake options

| Option | Default | Meaning |
|--------|---------|---------|
| `DLLIFTING_BUILD_SHARED` | OFF | Build shared library |
| `DLLIFTING_REDUCTION` | ON | Capacity reduction (`tableleft`) |
| `DLLIFTING_BUILD_EXAMPLES` | ON | Example binaries |
| `DLLIFTING_BUILD_TESTS` | ON | `test_dllifting`, `test_isgeq`, … |

## Testing

See [docs/TESTING.md](docs/TESTING.md) for test coverage and known limitations on `>=` knapsacks.

## License

MIT License — see [LICENSE](LICENSE). Based on work by Igor Vasilyev and contributors.

## Citation

If you use this library in research, please cite the associated knapsack lifting / separation work from your publication and link to this repository.
