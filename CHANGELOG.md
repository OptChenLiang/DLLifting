# Changelog

## [1.2.0] - 2026-08-06

### Changed

- Public-release hygiene: library no longer calls `exit()` on realloc failure; debug
  `printf` is gated behind `DLLIFTING_DEBUG`; release builds use `-DNDEBUG`.
- `make test` runs the core (`<=`) suite only and fails the process on any core failure.
- `make test-all` runs extended `>=` and mixed suites (`./test_dllifting --all`).
- README documents `threshold` / `isdl_mode`, full `lifting` / `dllifting_lift_cover`
  signatures, and marks `>=` knapsacks as experimental.
- `make install` installs both `DLLifting.h` and `dllifting_c.h`; default `PREFIX=$(HOME)/.local`.
- Single translation unit: C ABI lives in `src/DLLifting.cpp`.
- Version macros: `DLLIFTING_VERSION` (`1.2.0`).

### Fixed

- Removed dead realloc/`exit(0)` path in `Lifting_Mergesort`.

## [1.1.0] - 2026-05-28

### Added

- `isdl_mode` on `lifting()` / `dllifting_lift_cover()` (`AUTO` / `DL` / `DP`).

## [1.0.0] - 2026-05

### Added

- Initial standalone release.
