# DLLifting

DL / DP hybrid **coefficient lifting** for knapsack cover inequalities (`<=` and `>=` knapsacks).

## Build and install

```bash
cd DLLifting
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local \
         -DDLLIFTING_REDUCTION=ON \
         -DDLLIFTING_BUILD_EXAMPLES=ON \
         -DDLLIFTING_BUILD_TESTS=ON
cmake --build .
ctest          # optional
cmake --install .
```

Legacy Makefile (tests only, no install):

```bash
make run
make run-reduction   # if REDUCTION enabled via -DDLLIFTING_REDUCTION in CMake lib
```

## Use in your project

### CMake

```cmake
find_package(dllifting CONFIG REQUIRED PATHS "$ENV{HOME}/.local/lib/cmake/dllifting")
add_executable(my_cut_gen my_cut_gen.cpp)
target_link_libraries(my_cut_gen PRIVATE dllifting::dllifting)
```

Include: `#include <dllifting/DLLifting.h>` or add `-I.../include` and `#include "DLLifting.h"`.

### Manual link

```bash
g++ -O2 my_app.cpp -I/path/to/DLLifting -L/path/to/lib -ldllifting -lm -o my_app
```

## Public API

| Entry | Language | Description |
|-------|----------|-------------|
| `lifting()` | C++ | Full lift; see `DLLifting.h` |
| `dllifting_lift_cover()` | C | Wrapper in `dllifting/dllifting_c.h` |

Examples: `examples/example_basic.cpp`, `examples/example_c_api.c`.

Solver callbacks (CPLEX, Gurobi, …): see [docs/INTEGRATION_SOLVERS.md](docs/INTEGRATION_SOLVERS.md).

## Files

| File | Role |
|------|------|
| `DLLifting.h` | C++ API |
| `DLLifting.cpp` | Implementation |
| `include/dllifting/dllifting_c.h` | C API |
| `src/dllifting_c.cpp` | C wrapper |

## Options (CMake)

| Option | Default | Meaning |
|--------|---------|---------|
| `DLLIFTING_BUILD_SHARED` | OFF | Shared library |
| `DLLIFTING_REDUCTION` | ON | Capacity reduction (`tableleft`) |
| `DLLIFTING_BUILD_EXAMPLES` | ON | Example binaries |
| `DLLIFTING_BUILD_TESTS` | ON | `test_dllifting`, `test_isgeq` |

## License

Same terms as the parent `kplifting` project / original author notice in `ExactSeparation`.
