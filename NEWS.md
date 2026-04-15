# treeowen 0.2.1

## Bug fix

* Fixed persistent `cpp=FALSE` error on macOS arm64:
  - Changed `R_useDynamicSymbols(dll, FALSE)` to `TRUE` in
    `src/RcppExports.cpp` so dynamic symbol lookup is allowed.
    This removes the strict requirement that every `.Call()` name
    must exactly match the `R_CallMethodDef` table entry, which was
    causing `treeowen_ping` lookups to fail silently on some platforms.
  - Updated `.onLoad()` and `.try_load_cpp()` to first call
    `treeowen_ping()` via the `RcppExports.R` wrapper (the most
    reliable path), then fall back to `_treeowen_treeowen_ping` and
    `treeowen_ping` direct `.Call()` variants.

# treeowen 0.2.1

## Bug fix

* Fixed `cpp=FALSE` / C++ XPtr missing error that persisted after v0.2.0.

  **Root cause**: After `R_registerRoutines()` + `R_useDynamicSymbols(dll, FALSE)`,
  all C symbols are registered under the `_treeowen_` prefix
  (e.g. `_treeowen_treeowen_ping`). The `.onLoad` hook and `.try_load_cpp()`
  were still calling `.Call("treeowen_ping")` (without the prefix), which
  failed silently and set `cpp=FALSE`.

  **Fix**:
  - `zzz.R` `.onLoad`: now calls `.Call("_treeowen_treeowen_ping")` first,
    with a fallback to the plain name.
  - `R/utils.R` `.try_load_cpp()`: same two-step probe.
  - `src/RcppExports.cpp` `CallEntries`: also registers `"treeowen_ping"`
    as an alias so both calling conventions work.

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
