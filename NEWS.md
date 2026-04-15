# treeowen 0.1.2

## Bug fix

* Fixed compilation failure on macOS arm64 (Apple clang 17, /opt/R/arm64
  layout) where `R CMD INSTALL` did not inject the Rcpp include path
  automatically despite `LinkingTo: Rcpp` being set in DESCRIPTION.
  `src/Makevars` now calls `Rscript -e "cat(Rcpp:::CxxFlags())"` at
  build time to retrieve the correct `-I` flag on every platform.

# treeowen 0.1.1

## New features

* Added `.lightgbm_unify_compat()` — robust replacement for
  `treeshap::lightgbm.unify()`. Handles all known failure modes across
  lightgbm 3.x / 4.x and treeshap 0.3–0.4:
  - Three-level fallback: treeshap → `lgb.model.dt.tree()` → `dump_model()` JSON.
  - Column name aliases: `split_gain`/`value`, `threshold`/`split_point`,
    `count`/`internal_count`/`leaf_count`.
  - Three child-ID encoding schemes (negative leaf index, positive node
    index, string `"N<k>"/"L<k>"`), all resolved to row indices.
  - Feature names from Booster private field → `model$feature_name()` →
    `colnames(data)` fallback chain.
  - Internal helper `.lgb_parse_dump_model()` for JSON tree dump parsing.

* Added `.ranger_unify_compat()` — robust replacement for
  `treeshap::ranger.unify()`. Handles all known failure modes across
  ranger 0.11–0.16+ and treeshap 0.3–0.4:
  - Three-level fallback: treeshap → patched `treeInfo()` +
    `treeshap:::ranger_unify.common()` → direct `model$forest` parsing.
  - `treeInfo()` leaf prediction column: `pred.1` (old ranger) and
    `prediction` (ranger >= 0.14) both handled.
  - `treeInfo()` column name aliases: `nodeID`/`node`, `splitvarID`/
    `splitVarID`, `splitval`/`splitVal`, `terminal`/`isTerminal`.
  - `child.nodeIDs` stored as list-of-2-vectors (old) or separate
    `leftChild`/`rightChild` columns (new); both formats parsed.
  - Leaf predictions as scalar, matrix, or list-of-vectors; probability for
    class 1 extracted in all cases.
  - Internal helpers `.ranger_treeinfo_norm()` and
    `.ranger_tree_from_forest()`.
  - Non-probability forest: warning (not error).
  - `write.forest = FALSE`: informative error.

* Improved `.xgboost_unify_compat()`:
  - `xgb.model.dt.tree()` argument name: tries both `model=` and
    `xgb_model=` (xgboost 2.x rename).
  - treeshap bypassed entirely; model parsed directly from tree dump.

## Tests

* Added comprehensive unit tests (`tests/testthat/`):
  - `test-xgboost-unify.R` — 9 tests covering all [XGB-1..7] failure modes.
  - `test-lightgbm-unify.R` — 12 tests covering all [LGB-1..7] failure modes.
  - `test-ranger-unify.R` — 12 tests covering all [RNG-1..8] failure modes.
  - `test-treeowen-core.R` — 9 integration tests (efficiency axiom for all
    three learners, importance, print, method dispatch).
  - `test-utils.R` — 14 unit tests for internal helpers (Shapley weights,
    bitset operations, hierarchy builders, cache management).

## Documentation

* README rewritten in plain prose — all LaTeX formulas replaced with
  intuitive descriptions. Added version compatibility tables for xgboost,
  lightgbm, and ranger.
* All man pages (`man/*.Rd`): `\eqn{}` and `\deqn{}` math expressions
  replaced with plain-text equivalents.
* `treeowen.Rd`: updated `unified_model` parameter documentation to
  reference the new compatibility wrappers.
* `treeowen-package.Rd`: description rewritten in plain prose.

# treeowen 0.1.0

## Initial release

* Exact and Monte Carlo Owen value computation for XGBoost, LightGBM, and Ranger.
* Hierarchy-guided group aggregation with tree-aware dynamic programming.
* Flat beeswarm visualization (`treeowen_beeswarm()`).
* Hierarchical beeswarm visualization (`treeowen_hierarchical_beeswarm()`).
* Feature and group importance (`treeowen_importance()`).
* Reference enumerators (`treeowen_exact_enum()`, `treeowen_exact_tuvalues()`).
* Custom auxiliary tree builders (`build_hierarchy_tree_binary()`,
  `build_hierarchy_tree_from_layers()`).
* Bundled dataset: `immuno` — gut microbiome data from 219 cancer
  immunotherapy patients across five cohorts.
