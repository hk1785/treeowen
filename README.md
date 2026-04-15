# treeowen

Title: Efficient Game-Theoretic Explanations for Tree-Based Ensembles via Owen Values

Version: 0.1.0

Date: 2026-04-15

Author: Hyunwook Koh

Maintainer: Hyunwook Koh (hyunwook.koh@stonybrook.edu)

Description: This R package (treeowen) provides exact and Monte Carlo algorithms for computing Owen values in tree-based ensemble models, including gradient boosting machines (XGBoost, LightGBM) and random forests (Ranger). The Owen value extends the Shapley value to settings where features are organized into *a priori* groups, allocating attribution in two stages: first across groups, then within each group. This two-stage structure yields coherent explanations at both the group level and the individual feature level simultaneously. The package also includes global importance measures and visualization tools (flat and hierarchical beeswarm plots).

Depends: R(>= 4.1.0), Rcpp, parallel, grid, grDevices, stats, utils

License: MIT

## Reference

* Koh H. _Efficient Game-Theoretic Explanations for Tree-Based Ensembles via Owen Values._ (_In Review_)

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

## 1. Main Functions

## :mag: treeowen

### Description
Computes exact or Monte Carlo Owen values for every observation in `x` using a trained tree-based ensemble and a feature partition `groups`. The Owen value allocates model attribution in two stages: an outer Shapley game across groups and an inner Shapley game within each group.

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
  dp_newline_every  = 10L,
  eta_alpha         = 0.2,
  n_cores           = 1L,
  use_cpp           = TRUE,
  verbose           = TRUE
)
```

### Arguments
* _unified_model_ - A unified model object produced by `treeshap::xgboost.unify()`, `treeshap::lightgbm.unify()`, or `treeshap::ranger.unify()`.
* _x_ - A numeric data frame or matrix of observations (_n_ rows × _p_ columns). Must contain the same features as `unified_model`.
* _groups_ - A named list of character vectors. Each element names the features belonging to that group. Together they must form a complete, non-overlapping partition of all features in `unified_model`.
* _method_ - Algorithm selection: `"auto"` (default; exact for groups with |G_k| < `auto_exact_max_m`, Monte Carlo otherwise), `"exact"` (force exact for all groups; requires |G_k| ≤ 60), or `"approx"` (force Monte Carlo for all groups).
* _hierarchy_ - Controls the auxiliary binary tree over groups used to organise the outer Shapley summation. `NULL` (default) builds a balanced binary tree safe for any _K_. A nested list, numeric vector, or integer matrix are also accepted. The choice does not affect correctness, only runtime.
* _n_inner_mc_ - Integer. Initial Monte Carlo permutation pairs per context when `method = "approx"` or `"auto"`. Default `64L`.
* _inner_antithetic_ - Logical. Antithetic sampling to reduce Monte Carlo variance. Default `TRUE`.
* _target_se_inner_ - Numeric. Target standard error for adaptive Monte Carlo stopping. Default `1e-3`.
* _min_inner_mc_ - Integer. Minimum Monte Carlo pairs per context. Default `32L`.
* _max_inner_mc_ - Integer. Maximum Monte Carlo pairs per context. Default `512L`.
* _chunk_size_inner_ - Integer. Batch size for inner forest evaluation. Default `131072L`.
* _max_bytes_ - Numeric. Memory limit (bytes) for batched evaluation. Default `512 * 1024^2` (512 MB).
* _inner_bitmask_max_ - Integer. Maximum |G_k| for bitmask-based exact enumeration. Default `TREEOWEN_INNER_BITMASK_DEFAULT` (60).
* _auto_exact_max_m_ - Integer. Groups with |G_k| ≥ `auto_exact_max_m` use Monte Carlo when `method = "auto"`. Default `TREEOWEN_AUTO_EXACT_MAX_M` (30).
* _check_efficiency_ - Logical. Verify the efficiency axiom Σψ_i = f(x) − E[f(X)] after computation. Default `FALSE`.
* _efficiency_tol_ - Numeric. Tolerance for the efficiency check. Default `1e-6`.
* _dp_progress_ - Logical. Show a progress bar. Default `TRUE`.
* _dp_print_every_ - Integer. Progress bar update frequency (rows). Default `1L`.
* _dp_newline_every_ - Integer. Deprecated; kept for backwards compatibility. Default `10L`.
* _eta_alpha_ - Numeric. Exponential moving average smoothing factor for ETA estimation. Default `0.2`.
* _n_cores_ - Integer. Parallel cores. Default `1L`. Values > 1 use `parallel::mclapply()` on Linux/macOS. Windows requires `n_cores = 1L`.
* _use_cpp_ - Logical. Use compiled C++ routines (strongly recommended). Default `TRUE`.
* _verbose_ - Logical. Print diagnostic messages. Default `TRUE`.

### Values
An object of class `"treeowen_result"` (a named list) with the following components:
* _owens_ - A data frame (_n_ × _p_) of Owen values. Column names match the features in `x`.
* _observations_ - A data frame (_n_ × _p_) of the input observations.
* _groups_ - The `groups` argument as supplied.
* _baseline_ - Numeric vector of length _n_: the empty-set forest value v_x(∅) = E[f(X)] per observation.
* _v_all_ - Numeric vector of length _n_: the full-set value v_x(N) = f(x) per observation.
* _stats_ - Named list of diagnostics: _n_, _p_, _K_, tree counts, effective method, group methods, C++ usage, and parallelism settings.
* _note_ - Character string summarising the run configuration.

### Example
Import requisite R packages
```r
library(treeowen)
library(xgboost)
library(treeshap)
```
Generate synthetic data (5 groups × 5 features, _n_ = 50)
```r
set.seed(42)
feat_names <- paste0("F", rep(1:5, times = 5), "G", rep(1:5, each = 5))
X <- as.data.frame(matrix(rnorm(50 * 25), 50, 25,
                          dimnames = list(NULL, feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2
Y      <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:5, function(k) paste0("F", 1:5, "G", k)),
                   paste0("G", 1:5))
```
Fit XGBoost and unify
```r
xgb_mod <- xgboost(data      = xgb.DMatrix(as.matrix(X), label = Y),
                   nrounds   = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
unified <- xgboost.unify(xgb_mod, X)
```
Compute Owen values
```r
result <- treeowen(unified, X, groups, method = "auto",
                   dp_progress = FALSE, verbose = TRUE)
print(result)
dim(result$owens)         # 50 x 25
head(result$owens[, 1:5]) # first 5 Owen values per observation
```
Verify efficiency axiom (should be < 1e-6 for exact)
```r
eff_err <- abs(rowSums(result$owens) - (result$v_all - result$baseline))
cat(sprintf("Max efficiency error: %.2e\n", max(eff_err)))
```
More Details
```r
?treeowen
```

## :mag: treeowen_importance

### Description
Derives feature-level and group-level importance scores from a `treeowen` result. Feature importance is the mean absolute Owen value across observations; group importance is the mean absolute sum of within-group Owen values (group-level attribution).

### Syntax
```r
treeowen_importance(
  result,
  type      = c("both", "feature", "group"),
  group_agg = c("sum_abs", "sum_abs_feat", "both"),
  sort      = TRUE,
  normalize = FALSE
)
```

### Arguments
* _result_ - An object of class `"treeowen_result"` from `treeowen`.
* _type_ - Character. Which scores to compute: `"both"` (default; both feature and group importance), `"feature"` (feature importance only), or `"group"` (group importance only).
* _group_agg_ - Character. How to compute group importance: `"sum_abs"` (default; I_k = (1/n) Σ_j |Ψ_k^(j)|, where Ψ_k^(j) = Σ_{i∈G_k} ψ_i^(j)), `"sum_abs_feat"` (I_k = (1/|G_k|) Σ_{i∈G_k} I_i, average of within-group feature importances), or `"both"` (compute both; `importance` column uses `"sum_abs"` for sorting; an extra column `importance_sum_abs_feat` is appended).
* _sort_ - Logical. Sort by decreasing importance. Default `TRUE`.
* _normalize_ - Logical. Normalize scores to sum to 1. Default `FALSE`.

### Values
An object of class `"treeowen_importance"` (a named list) with the following components:
* _feature_ - Data frame with columns `feature` and `importance`, sorted by decreasing importance. Present when `type %in% c("both", "feature")`.
* _group_ - Data frame with columns `group` and `importance` (plus `importance_sum_abs_feat` when `group_agg = "both"`), sorted by decreasing importance. Present when `type %in% c("both", "group")`.
* _group_attr_ - Numeric matrix (_n_ × _K_) of group-level Owen attributions Ψ_k^(j). Always present; used internally by plotting functions.
* _group_pos_ - Named list: group name → column indices in `result$owens`.
* _group_names_ - Character vector of group names.
* _feat_names_ - Character vector of feature names.
* _n_, _K_, _p_ - Integer scalars: number of observations, groups, and features.
* _type_, _group_agg_ - The arguments as supplied.

### Example
Import requisite R packages
```r
library(treeowen)
library(xgboost)
library(treeshap)
```
(Generate synthetic data and fit XGBoost as shown in `?treeowen`)
```r
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=5), "G", rep(1:5, each=5))
X <- as.data.frame(matrix(rnorm(50*25), 50, 25, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2
Y      <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:5, function(k) paste0("F",1:5,"G",k)), paste0("G",1:5))
xgb_mod <- xgboost(data    = xgb.DMatrix(as.matrix(X), label=Y),
                   nrounds = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
result  <- treeowen(xgboost.unify(xgb_mod, X), X, groups,
                    method = "auto", dp_progress = FALSE, verbose = FALSE)
```
Compute feature and group importance
```r
imp <- treeowen_importance(result, type = "both", group_agg = "sum_abs")
print(imp)
head(imp$feature, 5L)  # top 5 features
imp$group              # group importance
```
Normalised importance (sums to 1)
```r
imp_norm <- treeowen_importance(result, type = "group", normalize = TRUE)
imp_norm$group
```
More Details
```r
?treeowen_importance
```

---------------------------------------------------------------------------------------------------------------------------------------

## 2. Visualization Tools

## :mag: treeowen_beeswarm

### Description
Produces SHAP-style beeswarm plots of Owen value distributions at the feature level, the group level, or both panels side by side. Points are coloured by the feature's raw value (feature panel) or by a group-level summary statistic (group panel).

### Syntax
```r
treeowen_beeswarm(
  result,
  level             = c("feature", "group", "both"),
  group_agg         = c("sum_abs", "sum_abs_feat"),
  normalize_imp     = FALSE,
  top_n_feature     = 10L,
  top_n_group       = 10L,
  group_color_stat  = c("sum", "mean", "custom"),
  group_color_fn    = NULL,
  color_low         = "#0052A5",
  color_high        = "#DC2626",
  color_alpha       = 0.6,
  qlims             = c(0.05, 0.95),
  point_size        = 1.6,
  point_alpha       = 0.55,
  quasirandom_width = 0.4,
  title_feature     = NULL,
  title_group       = NULL,
  xlab_feature      = "Owen Value",
  xlab_group        = "Group Owen Value",
  legend_label_feat = "Value",
  legend_label_grp  = NULL,
  plot_title_size   = 16,
  axis_title_x_size = 11,
  axis_text_x_size  = 10,
  axis_text_y_size  = 12,
  legend_text_size  = 9,
  legend_title_size = 9,
  margin_t          = 6,
  margin_r          = 6,
  margin_b          = 6,
  margin_l          = 8,
  legend_position   = "bottom",
  legend_direction  = "horizontal",
  legend_barwidth   = grid::unit(170, "pt"),
  legend_barheight  = grid::unit(10,  "pt"),
  verbose           = FALSE
)
```

### Arguments
* _result_ - An object of class `"treeowen_result"` from `treeowen`.
* _level_ - Character. Panel(s) to produce: `"feature"` (default; one beeswarm row per feature), `"group"` (one row per group), or `"both"` (group | feature side-by-side via `patchwork`; requires `patchwork`).
* _group_agg_ - Character. Aggregation for sort order: `"sum_abs"` (default; mean absolute group attribution) or `"sum_abs_feat"` (mean of within-group feature importances). `"both"` is not permitted (ambiguous sort order).
* _normalize_imp_ - Logical. Normalise importance before sorting. Default `FALSE`.
* _top_n_feature_ - Integer or `NULL`. Number of top features to show. Default `10L`.
* _top_n_group_ - Integer or `NULL`. Number of top groups to show. Default `10L`.
* _group_color_stat_ - Character. Colour statistic for group panel: `"sum"` (default; sum of raw feature values), `"mean"`, or `"custom"`.
* _group_color_fn_ - Function or `NULL`. Required when `group_color_stat = "custom"`: `function(x_mat, feat_idx) -> numeric(n)`.
* _color_low_ - Hex colour for the low end of the gradient. Default `"#0052A5"` (blue).
* _color_high_ - Hex colour for the high end. Default `"#DC2626"` (red).
* _color_alpha_ - Numeric in [0,1]. Gradient transparency. Default `0.6`.
* _qlims_ - Numeric vector of length 2. Quantile limits for colour clipping. Default `c(0.05, 0.95)`.
* _point_size_ - Numeric. Point size. Default `1.6`.
* _point_alpha_ - Numeric in [0,1]. Point transparency. Default `0.55`.
* _quasirandom_width_ - Numeric. Quasirandom jitter spread. Default `0.4`.
* _title_feature_ - Character or `NULL`. Feature panel title.
* _title_group_ - Character or `NULL`. Group panel title.
* _xlab_feature_ - Character. x-axis label for the feature panel. Default `"Owen Value"`.
* _xlab_group_ - Character. x-axis label for the group panel. Default `"Group Owen Value"`.
* _legend_label_feat_ - Character. Colour-bar label for the feature panel. Default `"Value"`.
* _legend_label_grp_ - Character or `NULL`. Colour-bar label for the group panel. `NULL` auto-generates from `group_color_stat`.
* _plot_title_size_ - Numeric (pt). Default `16`.
* _axis_title_x_size_ - Numeric (pt). Default `11`.
* _axis_text_x_size_ - Numeric (pt). Default `10`.
* _axis_text_y_size_ - Numeric (pt). Default `12`.
* _legend_text_size_ - Numeric (pt). Default `9`.
* _legend_title_size_ - Numeric (pt). Default `9`.
* _margin_t_, _margin_r_, _margin_b_, _margin_l_ - Numeric (pt). Plot margins. Defaults: `6`, `6`, `6`, `8`.
* _legend_position_ - Character. Default `"bottom"`.
* _legend_direction_ - Character. Default `"horizontal"`.
* _legend_barwidth_ - `grid::unit`. Default `grid::unit(170, "pt")`.
* _legend_barheight_ - `grid::unit`. Default `grid::unit(10, "pt")`.
* _verbose_ - Logical. Print diagnostic messages. Default `FALSE`.

### Values
* When `level = "feature"` or `"group"`: a `ggplot` object.
* When `level = "both"`: an object of class `"treeowen_plots"` (a named list) with components: _combined_ (patchwork: group | feature side-by-side), _feature_ (standalone feature ggplot), _group_ (standalone group ggplot), _feat_ordered_ (features in display order), _grp_ordered_ (groups in display order).

### Example
Import requisite R packages
```r
library(treeowen)
library(xgboost)
library(treeshap)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
```
(Generate synthetic data, fit XGBoost, and compute Owen values as shown in `?treeowen`)
```r
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=5), "G", rep(1:5, each=5))
X <- as.data.frame(matrix(rnorm(50*25), 50, 25, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2
Y      <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:5, function(k) paste0("F",1:5,"G",k)), paste0("G",1:5))
xgb_mod <- xgboost(data    = xgb.DMatrix(as.matrix(X), label=Y),
                   nrounds = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
result  <- treeowen(xgboost.unify(xgb_mod, X), X, groups,
                    method = "auto", dp_progress = FALSE, verbose = FALSE)
```
Draw a feature-level beeswarm plot
```r
p_feat <- treeowen_beeswarm(result, level = "feature", top_n_feature = 15L,
                             xlab_feature = "Owen Value (XGBoost)")
p_feat
```
Draw a group-level beeswarm plot
```r
p_grp <- treeowen_beeswarm(result, level = "group", top_n_group = 5L,
                            group_color_stat = "sum")
p_grp
```
Draw both panels side-by-side
```r
out <- treeowen_beeswarm(result, level = "both",
                          top_n_feature = 10L, top_n_group = 5L)
out$combined
```
More Details
```r
?treeowen_beeswarm
```

## :mag: treeowen_hierarchical_beeswarm

### Description
Produces a hierarchical beeswarm plot in which each group and its nested features are displayed as individual rows, assembled into a multi-column layout using `patchwork`. Groups are ordered by decreasing group importance; features within each group by decreasing feature importance. An optional standalone colour-bar strip is appended at the bottom of each page.

### Syntax
```r
treeowen_hierarchical_beeswarm(
  ow_result,
  imp,
  top_n_group      = NULL,
  top_n_feat       = NULL,
  n_col            = 2L,
  group_color_stat = "sum",
  color_low        = "#0052A5",
  color_high       = "#DC2626",
  color_alpha      = 0.65,
  qlims            = c(0.05, 0.95),
  point_size_grp   = 2.2,
  point_size_feat  = 1.2,
  point_alpha      = 0.50,
  row_h_grp        = 0.55,
  row_h_feat       = 0.22,
  width_in         = 14,
  max_h_per_page   = 25,
  feat_axis_size   = 7,
  grp_axis_size    = 8.5,
  xlab             = "Owen Value",
  show_colorbar    = TRUE,
  colorbar_h       = 0.35,
  lname            = "treeowen",
  save_path        = NULL,
  dpi              = 300L,
  verbose          = FALSE
)
```

### Arguments
* _ow_result_ - An object of class `"treeowen_result"` from `treeowen`.
* _imp_ - An object of class `"treeowen_importance"` from `treeowen_importance` computed with `type = "both"` (must contain both `$feature` and `$group`).
* _top_n_group_ - Integer or `NULL`. Number of top groups by group importance to display. `NULL` (default) shows all groups.
* _top_n_feat_ - Integer or `NULL`. Restrict features to those among the top `top_n_feat` by feature importance that also belong to the selected groups. `NULL` (default) shows all features in the selected groups.
* _n_col_ - Integer. Number of layout columns. Default `2L`.
* _group_color_stat_ - Character. Colour statistic for group rows: `"sum"` (default; sum of raw feature values in the group) or `"mean"`.
* _color_low_ - Hex colour for the low gradient end. Default `"#0052A5"` (blue).
* _color_high_ - Hex colour for the high gradient end. Default `"#DC2626"` (red).
* _color_alpha_ - Numeric in [0,1]. Gradient transparency. Default `0.65`.
* _qlims_ - Numeric vector of length 2. Quantile limits for colour clipping. Default `c(0.05, 0.95)`.
* _point_size_grp_ - Numeric. Group row point size. Default `2.2`.
* _point_size_feat_ - Numeric. Feature row point size. Default `1.2`.
* _point_alpha_ - Numeric in [0,1]. Point transparency. Default `0.50`.
* _row_h_grp_ - Numeric. Group row height (inches). Default `0.55`.
* _row_h_feat_ - Numeric. Feature row height (inches). Default `0.22`.
* _width_in_ - Numeric. Total plot width (inches). Default `14`.
* _max_h_per_page_ - Numeric. Maximum plot height (inches) per page. When exceeded, output is split across multiple pages. Default `25`.
* _feat_axis_size_ - Numeric (pt). Feature row y-axis text size. Default `7`.
* _grp_axis_size_ - Numeric (pt). Group row y-axis text size. Default `8.5`.
* _xlab_ - Character. x-axis label. Default `"Owen Value"`.
* _show_colorbar_ - Logical. Append a standalone horizontal colour-bar strip at the bottom of each page. Default `TRUE`.
* _colorbar_h_ - Numeric. Colour-bar height (inches). Default `0.35`.
* _lname_ - Character. File-name prefix (used only when `save_path` is not `NULL`). Default `"treeowen"`.
* _save_path_ - Character or `NULL`. Directory path for PDF and PNG output. `NULL` (default) returns the list visibly without saving.
* _dpi_ - Integer. PNG resolution. Default `300L`.
* _verbose_ - Logical. Print diagnostic messages. Default `FALSE`.

### Values
A list of `patchwork` objects, one element per page. Each element is the assembled plot for that page (beeswarm columns + optional colour-bar strip). Returned visibly when `save_path = NULL`; invisibly when `save_path` is set (files written as `fig_<lname>_beeswarm_hierarchical[_<i>].pdf` and `.png`).

### Example
Import requisite R packages
```r
library(treeowen)
library(xgboost)
library(treeshap)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
```
(Generate synthetic data, fit XGBoost, and compute Owen values as shown in `?treeowen`)
```r
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=5), "G", rep(1:5, each=5))
X <- as.data.frame(matrix(rnorm(50*25), 50, 25, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2
Y      <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:5, function(k) paste0("F",1:5,"G",k)), paste0("G",1:5))
xgb_mod <- xgboost(data    = xgb.DMatrix(as.matrix(X), label=Y),
                   nrounds = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
result  <- treeowen(xgboost.unify(xgb_mod, X), X, groups,
                    method = "auto", dp_progress = FALSE, verbose = FALSE)
```
Compute importance (type = "both" is required)
```r
imp <- treeowen_importance(result, type = "both")
```
Draw a hierarchical beeswarm plot (return patchwork list)
```r
pages <- treeowen_hierarchical_beeswarm(
  ow_result        = result,
  imp              = imp,
  top_n_group      = 5L,
  top_n_feat       = NULL,
  n_col            = 2L,
  group_color_stat = "sum",
  show_colorbar    = TRUE,
  width_in         = 10,
  max_h_per_page   = 15,
  verbose          = TRUE
)
cat("Pages:", length(pages), "\n")
pages[[1]]
```
Save to disk (PDF + PNG)
```r
treeowen_hierarchical_beeswarm(
  ow_result  = result,
  imp        = imp,
  top_n_group = 5L,
  n_col      = 2L,
  width_in   = 10,
  lname      = "xgb_synthetic",
  save_path  = tempdir(),
  dpi        = 150L,
  verbose    = TRUE
)
```
More Details
```r
?treeowen_hierarchical_beeswarm
```

---------------------------------------------------------------------------------------------------------------------------------------

## 3. Example Dataset

## :mag: immuno

### Description
Relative abundance data of gut microbial taxa from 219 cancer patients treated with immune checkpoint inhibitor (ICI) therapy, pooled from five independent cohorts. The dataset includes a pre-built group partition `groups` that maps each of the 18 genera to its member species, ready for direct use with `treeowen`.

### Usage
```r
data(immuno)
```

### Format
Loading `data(immuno)` places two objects in the workspace:
* _immuno_ - A numeric data frame with 219 rows (patients) and 190 columns (microbial features: species).
* _groups_ - A named list with 18 elements. Each element is a character vector of species column names belonging to that genus (e.g., `"Akkermansia"`, `"Bacteroides"`). Together, `groups` forms a complete, non-overlapping partition of all 190 columns of `immuno`, suitable for direct use as the `groups` argument of `treeowen`.

### References
* Koh H. _Efficient Game-Theoretic Explanations for Tree-Based Ensembles via Owen Values._ (_In Review_)

### Example
```r
data(immuno)

# Dimensions
dim(immuno)          # 219 x 190
length(groups)       # 18 genera
names(groups)[1:5]   # first five genus names

# Response labels from row names
resp_labels <- sub(".*_(R|NR)$", "\\1", rownames(immuno))
table(resp_labels)   # 65 responders (R) and 154 non-responders (NR).

# Group sizes (species per genus)
sort(lengths(groups), decreasing = TRUE)
```
More Details
```r
?immuno
```
