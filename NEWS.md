# treeowen 0.2.1

## Bug fix (cpp=FALSE)

* Root cause identified and fixed: `.onLoad()` and `.try_load_cpp()`
  were calling `.Call("treeowen_ping", PACKAGE = pkgname)`, but
  `R_useDynamicSymbols(dll, FALSE)` is set in `R_init_treeowen()`,
  which means **only names explicitly listed in `CallEntries`** are
  reachable via `.Call()`.  The `CallEntries` table registers the
  symbol as `"_treeowen_treeowen_ping"` (with the `_treeowen_` prefix
  that Rcpp adds), not as the bare `"treeowen_ping"`.  The mismatch
  caused the ping call to fail silently, setting `cpp=FALSE` every
  time the package was loaded.

  Fix: both `.onLoad()` and `.try_load_cpp()` now call
  `.Call("_treeowen_treeowen_ping", PACKAGE = pkgname)`, which exactly
  matches the registered name.

# treeowen 0.2.0

## Bug fixes

* **C++ not loaded (`cpp=FALSE`)**: `RcppExports.R` and `RcppExports.cpp`
  were missing from the package, so the shared library was never registered
  with R. Added both files with correct function signatures and the
  `R_init_treeowen()` registration entry point. The package now compiles
  and links correctly on all platforms.

* **`data.table` `:=` operator fails (`cedta()` error)**: `data.table` was
  listed under `Suggests` (optional), so R's calling-environment check
  (`cedta()`) blocked `:=` at runtime even when the package was installed.
  Moved `data.table` to `Imports` and added `import(data.table)` to
  `NAMESPACE`.

* **`could not find function ".xgboost_unify_compat"`**: Functions whose
  names begin with `.` are hidden from the R search path even when exported.
  Renamed all three wrappers:
  - `.xgboost_unify_compat()`  →  `xgboost_unify_compat()`
  - `.lightgbm_unify_compat()` →  `lightgbm_unify_compat()`
  - `.ranger_unify_compat()`   →  `ranger_unify_compat()`

* **`Rcpp.h` not found on macOS arm64**: `src/Makevars` now calls
  `Rscript -e "cat(Rcpp:::CxxFlags())"` at build time to inject the correct
  `-I` path, which fixes compilation on `/opt/R/arm64` and similar
  non-standard R installations.

## Documentation

* All man pages: replaced LaTeX math (`\eqn{}`, `\sum`, `\subseteq`, etc.)
  with plain prose.
* Every `\examples{}` block now includes the required `library()` calls so
  examples are self-contained.
* `README.md` rewritten in plain prose with version-compatibility tables for
  XGBoost, LightGBM, and Ranger.

# treeowen 0.1.0

## Initial release

* Exact and Monte Carlo Owen value computation for XGBoost, LightGBM, and
  Ranger.
* Hierarchy-guided group aggregation with tree-aware dynamic programming.
* Flat beeswarm visualization (`treeowen_beeswarm()`).
* Hierarchical beeswarm visualization (`treeowen_hierarchical_beeswarm()`).
* Feature and group importance (`treeowen_importance()`).
* Reference enumerators (`treeowen_exact_enum()`,
  `treeowen_exact_tuvalues()`).
* Custom auxiliary tree builders (`build_hierarchy_tree_binary()`,
  `build_hierarchy_tree_from_layers()`).
* Bundled dataset: `immuno` — gut microbiome data from 219 cancer
  immunotherapy patients across five cohorts.
