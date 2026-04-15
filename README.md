# treeowen

Title: Efficient Game-Theoretic Explanations for Tree-Based Ensembles via Owen Values

Version: 0.1.3

Date: 2026-04-15

Author: Hyunwook Koh

Maintainer: Hyunwook Koh (hyunwook.koh@stonybrook.edu)

Description: This R package (treeowen) provides exact and Monte Carlo algorithms for computing Owen values in tree-based ensemble models, including gradient boosting machines (XGBoost, LightGBM) and random forests (Ranger). The Owen value extends the Shapley value to settings where features are organized into *a priori* groups, allocating attribution in two stages: first across groups, then within each group. This two-stage structure yields coherent explanations at both the group level and the individual feature level simultaneously. The package also includes global importance measures and visualization tools (flat and hierarchical beeswarm plots).

Depends: R(>= 4.1.0), Rcpp, parallel, grid, grDevices, stats, utils

License: MIT

## Reference

* Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles via Owen Values. (_In Review_)

## Troubleshooting Tips

If you have any problems using this R package, please report them via Issues (https://github.com/hk1785/treeowen/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).

## Prerequisites

xgboost, lightgbm, ranger, treeshap, ggplot2, ggbeeswarm, patchwork, scales
```r
install.packages(c("xgboost", "lightgbm", "ranger", "treeshap",
                   "ggplot2", "ggbeeswarm", "patchwork", "scales"))
```

## Installation

```r
library(devtools)
install_github("hk1785/treeowen", force = TRUE)
```

---------------------------------------------------------------------------------------------------------------------------------------

## 📋 Table of Contents

### 1. Main Functions
* :mag: **`treeowen`**: Computes Owen values, allocating model attribution in two stages — first across feature groups, then within each group.
* :mag: **`treeowen_importance`**: Derives feature-level and group-level importance scores from Owen values obtained using `treeowen`.

### 2. Visualization Tools
* :mag: **`treeowen_beeswarm`**: Flat beeswarm plot visualizing the distribution of Owen values at the feature level, group level, or both side-by-side.
* :mag: **`treeowen_hierarchical_beeswarm`**: Hierarchical beeswarm plot in which each group header row is followed by its nested feature rows, assembled into a multi-column patchwork layout.

### 3. Example Dataset
* :mag: **`immuno`**: Gut microbiome relative abundance data from 219 cancer patients treated with immune checkpoint inhibitor (ICI) therapy, pooled from five independent cohorts.

---------------------------------------------------------------------------------------------------------------------------------------

## 1. Model unification wrappers

Before computing Owen values you must convert your trained model into a unified format that treeowen can read. Use the wrapper functions below **instead** of calling `treeshap::xgboost.unify()`, `treeshap::lightgbm.unify()`, or `treeshap::ranger.unify()` directly. The wrappers automatically handle all known version-compatibility problems across different releases of each library.

### xgboost_unify_compat(model, data)

Converts a trained `xgb.Booster` to the unified format.

```r
unified <- xgboost_unify_compat(xgb_model, X)
```

Handled failure modes:

- `xgb.model.dt.tree()` argument was renamed from `model` to `xgb_model` in some builds; both are tried.
- The `ID` column is absent in older xgboost (before 1.7) and is reconstructed automatically from the `Tree` and `Node` columns.
- The split-gain column may be named `Quality` instead of `Gain` in newer xgboost; renamed automatically.
- `Yes`, `No`, and `Missing` columns may be stored as node-ID strings rather than row indices; converted automatically.
- `model$feature_names` may be NULL for models trained from unnamed matrices; falls back to `colnames(data)`.
- treeshap::xgboost.unify() is bypassed entirely because it calls internal xgboost functions that have changed between versions.

### lightgbm_unify_compat(model, data)

Converts a trained `lgb.Booster` to the unified format.

```r
unified <- lightgbm_unify_compat(lgb_model, X)
```

Handled failure modes:

- Tries `treeshap::lightgbm.unify()` first. If that fails, falls back to `lgb.model.dt.tree()`. If that also fails, parses the raw JSON dump from `model$dump_model()`.
- Column names changed across lightgbm versions: `split_gain` vs `value`, `threshold` vs `split_point`, `count` vs `internal_count` or `leaf_count`; all are normalised automatically.
- Three different child-node ID encodings depending on lightgbm version: negative-indexed leaf IDs, positive node IDs, and string format `"N<k>"` / `"L<k>"`; all are decoded to row indices.
- Feature names extracted from the Booster's private field, then `model$feature_name()`, then `colnames(data)` as a last resort.
- The best-round slot from `lgb.cv()` may be named `best_iter` or `best_iteration`; both are tried.
- The `verbose` argument may need to be passed inside `params` (older lightgbm) or directly to `lgb.cv()` (newer); use `-1L` in both places to be safe.

### ranger_unify_compat(model, data)

Converts a trained `ranger` model to the unified format. The model must have been trained with `probability = TRUE`.

```r
unified <- ranger_unify_compat(ranger_model, X)
```

Handled failure modes:

- Tries `treeshap::ranger.unify()` first. If that fails, patches `treeInfo()` output and calls `treeshap:::ranger_unify.common()` directly. If that also fails, builds the unified object from `model$forest` internals.
- `ranger::treeInfo()` returns `pred.1` in older ranger and `prediction` in newer ranger (0.14+); both are handled.
- `treeInfo()` column names vary across ranger versions: `nodeID` vs `node`, `splitvarID` vs `splitVarID`, `splitval` vs `splitVal`, `terminal` vs `isTerminal`; all aliases are normalised.
- `child.nodeIDs` may be a list of length-2 vectors (older ranger) or separate `leftChild` / `rightChild` columns (newer ranger); both formats are handled.
- Leaf predictions in probability forests may be a scalar, a matrix, or a list of vectors; the probability for class 1 is extracted in all cases.
- Models trained with `write.forest = FALSE` raise an informative error.
- Models trained with `probability = FALSE` emit a warning but still attempt unification.

---

## 2. treeowen()

Computes Owen values for every observation.

### Syntax

```r
treeowen(
  unified_model,
  x,
  groups,
  method            = c("auto", "exact", "approx"),
  hierarchy         = NULL,
  n_inner_mc        = 64L,
  inner_antithetic  = TRUE,
  target_se_inner   = 1e-3,
  min_inner_mc      = 32L,
  max_inner_mc      = 512L,
  chunk_size_inner  = 131072L,
  max_bytes         = 512 * 1024^2,
  inner_bitmask_max = TREEOWEN_INNER_BITMASK_DEFAULT,
  auto_exact_max_m  = TREEOWEN_AUTO_EXACT_MAX_M,
  check_efficiency  = FALSE,
  efficiency_tol    = 1e-6,
  dp_progress       = TRUE,
  dp_print_every    = 1L,
  n_cores           = 1L,
  use_cpp           = TRUE,
  verbose           = TRUE
)
```

### Key arguments

`unified_model` — Output of a `*_unify_compat()` wrapper.

`x` — Data frame or matrix of observations. Must contain the same features as the model.

`groups` — Named list of character vectors. Each element names the features belonging to that group. Every feature in `x` must appear in exactly one group (a full partition is required).

`method` — How to compute Owen values. `"auto"` uses the exact algorithm for groups with fewer features than `auto_exact_max_m` and Monte Carlo for larger groups. `"exact"` always uses exact. `"approx"` always uses Monte Carlo.

`hierarchy` — Structure of the auxiliary binary tree over groups. `NULL` (default) builds a balanced binary tree. Accepts a nested list or a matrix for custom structures. Does not affect the correctness of results, only computation speed.

`n_inner_mc`, `min_inner_mc`, `max_inner_mc`, `target_se_inner` — Control the number of Monte Carlo samples when using the approximate method.

`n_cores` — Number of CPU cores. Values greater than 1 use `parallel::mclapply()` on Linux and macOS. On Windows, keep this at 1.

`use_cpp` — Must be `TRUE`. The C++ backend is required.

### Return value

A `treeowen_result` list with:

- `owens` — Data frame of Owen values (rows = observations, columns = features).
- `baseline` — Average model prediction across all training observations.
- `v_all` — Model prediction for each observation.
- `groups` — The group partition used.
- `observations` — The input data `x`.
- `stats` — Diagnostics: number of trees, method used, group sizes, timing.
- `note` — Summary string.

For every observation, the sum of its Owen values equals its prediction minus the baseline.

### Example

```r
library(treeowen)
 
# ── Synthetic data ─────────────────────────────────────────────────────────
set.seed(42)
feat_names <- paste0("F", rep(1:5, times = 5), "G", rep(1:5, each = 5))
X <- as.data.frame(matrix(rnorm(50 * 25), 50, 25,
                           dimnames = list(NULL, feat_names)))
Y <- as.integer(0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2 > 0)
 
# Feature partition: 5 groups of 5 features each
groups <- setNames(lapply(1:5, function(k) paste0("F", 1:5, "G", k)),
                   paste0("G", 1:5))
 
# ── XGBoost ────────────────────────────────────────────────────────────────
library(xgboost)
model_xgb <- xgboost(xgb.DMatrix(as.matrix(X), label = Y),
                     nrounds = 100, max_depth = 3, eta = 0.1,
                     objective = "binary:logistic", verbose = 0)
result_xgb <- treeowen(xgboost_unify_compat(model_xgb, X), X, groups)
print(result_xgb)
 
# ── LightGBM ───────────────────────────────────────────────────────────────
library(lightgbm)
model_lgb <- lgb.train(
  params  = list(objective = "binary", learning_rate = 0.1,
                 max_depth = 3L, num_leaves = 7L, verbose = -1L),
  data    = lgb.Dataset(as.matrix(X), label = Y),
  nrounds = 100L, verbose = -1L)
result_lgb <- treeowen(lightgbm_unify_compat(model_lgb, X), X, groups)
print(result_lgb)
 
# ── Ranger ─────────────────────────────────────────────────────────────────
library(ranger)
model_rng <- ranger(.y ~ ., num.trees = 100L, max.depth = 3L,
                    probability = TRUE, keep.inbag = TRUE, seed = 42L,
                    data = cbind(X, .y = factor(Y, c(0,1), c("neg","pos"))))
result_rng <- treeowen(ranger_unify_compat(model_rng, X), X, groups)
print(result_rng)
```

---

## 3. treeowen_importance()

Derives importance scores from a `treeowen_result` object. Importance is the mean absolute Owen value across observations.

### Syntax

```r
treeowen_importance(
  result,
  type      = c("both", "feature", "group"),
  group_agg = c("sum_abs", "mean_abs", "sum", "mean"),
  sort      = TRUE,
  normalize = FALSE
)
```

### Return value

A list with `feature` (data frame with columns `feature` and `importance`), `group` (data frame with columns `group` and `importance`), and `group_attr` (matrix of group-level Owen values used by visualization functions).

### Example

```r
# result_xgb is produced by treeowen() — see Section 2 example above
imp <- treeowen_importance(result_xgb, type = "both", sort = TRUE)
print(imp$group)
```

---

## 4. treeowen_beeswarm()

Produces a beeswarm plot of Owen values. Each point is one observation. Point colour encodes the feature value (or a statistic over the group's features), so you can see at a glance which feature values are associated with positive or negative attribution.

### Syntax

```r
treeowen_beeswarm(
  result,
  level            = c("feature", "group", "both"),
  top_n_feature    = 20L,
  top_n_group      = NULL,
  group_color_stat = c("sum", "mean"),
  color_low        = "#0052A5",
  color_high       = "#DC2626",
  color_alpha      = 0.65,
  point_size       = 1.5,
  point_alpha      = 0.5,
  qlims            = c(0.05, 0.95),
  ...
)
```

Returns a `ggplot` object when `level = "feature"` or `"group"`, or a `treeowen_plots` list (with a `combined` patchwork figure) when `level = "both"`.

### Example

```r
library(ggplot2)
library(ggbeeswarm)  # for geom_quasirandom
library(patchwork)   # required for level = "both"

# result_xgb is produced by treeowen() — see Section 2 example above
p   <- treeowen_beeswarm(result_xgb, level = "feature", top_n_feature = 15L)
out <- treeowen_beeswarm(result_xgb, level = "both", top_n_feature = 10L, top_n_group = 5L)
print(out$combined)
```

---

## 5. treeowen_hierarchical_beeswarm()

Produces a hierarchical beeswarm layout in which each group header row is followed by beeswarm rows for the individual features in that group. The figure is assembled as a multi-column patchwork and can be saved to PDF and PNG.

### Syntax

```r
treeowen_hierarchical_beeswarm(
  ow_result,
  imp,
  top_n_group    = NULL,
  top_n_feat     = NULL,
  n_col          = 2L,
  show_colorbar  = TRUE,
  width_in       = 10,
  max_h_per_page = 20,
  lname          = "model",
  save_path      = NULL,
  dpi            = 300L,
  ...
)
```

`imp` must be from `treeowen_importance(type = "both")`. Returns a list of `patchwork` objects (one per page). When `save_path` is given, PDF and PNG files are written there.

### Example

```r
library(ggplot2)
library(ggbeeswarm)  # required
library(patchwork)   # required

# result_xgb is produced by treeowen() — see Section 2 example above
imp   <- treeowen_importance(result_xgb, type = "both")
pages <- treeowen_hierarchical_beeswarm(result_xgb, imp, top_n_group = 10L, n_col = 2L,
                                         save_path = "output/", lname = "xgboost")
print(pages[[1]])
```

---

## 6. Example dataset: immuno

Gut microbiome relative abundances from 219 cancer immunotherapy patients across five cohorts. Features are microbial species. The outcome is treatment response (1 = responder, 0 = non-responder).

```r
data(immuno)
immuno$X      # 219 x 953 matrix of relative abundances
immuno$Y      # numeric vector of length 219
immuno$groups # named list mapping genus names to species feature names
```