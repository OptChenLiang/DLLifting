# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-05-28

### Added

- **`isdl_mode` parameter** on `lifting()` and `dllifting_lift_cover()` to force DL or DP
  subproblem tables for the entire lift, without threshold-based switching.
- Mode constants: `DLLIFTING_MODE_AUTO`, `DLLIFTING_MODE_DP`, `DLLIFTING_MODE_DL`
  (exported in both `DLLifting.h` and `dllifting_c.h`).
- Internal field `DLLifting::force_mode`; when set to `_DL` or `_DP`, bounded and unbounded
  multiply/update paths stay on the chosen table type.
- Example `examples/example_force_mode.cpp` comparing forced DL vs forced DP on the same instance.
- Documentation: [docs/API.md](docs/API.md), [docs/MIGRATION.md](docs/MIGRATION.md).

### Changed

- **Breaking:** `lifting()` gains a trailing `int isdl_mode` argument (after `duration`).
- **Breaking:** `dllifting_lift_cover()` gains `int isdl_mode` before `x_frac`.
- `Lifting_Multiply` and `Lifting_update` respect `force_mode` for both bounded and unbounded items.
- README and solver integration guide updated for the new API.

### Migration

See [docs/MIGRATION.md](docs/MIGRATION.md). Existing callers should pass
`DLLIFTING_MODE_AUTO` as the last argument to preserve threshold switching.

## [1.0.0] - 2026-05

### Added

- Initial standalone release: hybrid DL/DP lifting, optional capacity reduction,
  C++ and C APIs, CMake package config, examples, and test suite.
