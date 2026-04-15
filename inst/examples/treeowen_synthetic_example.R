# =============================================================================
# treeowen: Synthetic Data Example
# =============================================================================
#
# Setup
#   - 5 groups:    G1 … G5
#   - 5 features per group: F1G1–F5G1, …, F1G5–F5G5  (25 features total)
#   - n = 50 observations, binary outcome Y
#   - Three learners: XGBoost, LightGBM, Ranger
#
# Workflow (one section per exported function / feature)
#    §1  Generate synthetic data
#    §2  Fit models + unify via treeshap
#    §3  treeowen()                    — Owen values
#    §4  treeowen_importance()         — feature & group importance
#    §5  treeowen_beeswarm()           — flat beeswarm plots
#    §6  treeowen_hierarchical_beeswarm()  — hierarchical beeswarm
#    §7  treeowen_exact_enum()         — model-agnostic reference enumerator
#    §8  treeowen_exact_tuvalues()     — tree-aware reference enumerator
#    §9  build_hierarchy_tree_*()      — custom auxiliary tree
#   §10  clear_inner_enum_cache()      — cache management
#   §11  Multicore (n_cores > 1)       — Linux / macOS only
#
# Required packages (install once):
#   install.packages(c("xgboost", "lightgbm", "ranger", "treeshap",
#                      "ggplot2", "ggbeeswarm", "patchwork", "scales"))
# =============================================================================

library(treeowen)

# ── 0. Dependency check ───────────────────────────────────────────────────────
needed       <- c("xgboost", "lightgbm", "ranger", "treeshap",
                  "ggplot2", "ggbeeswarm", "patchwork")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_pkgs))
  stop("Install missing packages first:\n  install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")

# =============================================================================
# §1. Synthetic data
# =============================================================================
set.seed(42)

n_obs    <- 50L
n_groups <- 5L
n_feat_g <- 5L              # features per group
p        <- n_groups * n_feat_g   # 25 features total

# Feature names: F{j}G{k}  (j = feature within group, k = group index)
# rep(1:5, times=5)  → 1 2 3 4 5 | 1 2 3 4 5 | ...  (sequence repeated)
# rep(1:5, each=5)   → 1 1 1 1 1 | 2 2 2 2 2 | ...  (element repeated)
# → F1G1 F2G1 F3G1 F4G1 F5G1 | F1G2 ... | ... | F1G5 ... F5G5
feat_names <- paste0("F", rep(seq_len(n_feat_g), times = n_groups),
                     "G", rep(seq_len(n_groups),  each  = n_feat_g))
stopifnot(length(feat_names) == p, !anyDuplicated(feat_names))

# Feature matrix (data frame; required by treeshap unifiers)
X <- as.data.frame(
  matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p,
         dimnames = list(NULL, feat_names))
)

# Binary outcome: mainly driven by G1 (F1G1, F2G1) and G3 (F1G3, F3G3)
log_odds <- 0.8 * X$F1G1 - 0.6 * X$F2G1 +
            0.5 * X$F1G3 + 0.4 * X$F3G3 +
            0.3 * X$F1G2 + 0.1 * rnorm(n_obs)
Y <- as.integer(log_odds > 0)

cat(sprintf("Data: n=%d  p=%d  K=%d  Pr(Y=1)=%.2f\n",
            n_obs, p, n_groups, mean(Y)))

# Group partition (each group gets exactly 5 features)
groups <- setNames(
  lapply(seq_len(n_groups), function(k)
    paste0("F", seq_len(n_feat_g), "G", k)),
  paste0("G", seq_len(n_groups))
)

# Verify partition is a complete, non-overlapping cover
all_grouped <- unlist(groups, use.names = FALSE)
stopifnot(
  setequal(all_grouped, feat_names),
  !anyDuplicated(all_grouped)
)
cat("Group partition verified: complete and non-overlapping.\n")

# =============================================================================
# §2. Fit models + unify via treeshap
# =============================================================================
# Notes on treeshap requirements:
#   XGBoost:  xgboost.unify(model, X)  — X = feature data frame (no Y)
#   LightGBM: lightgbm.unify(model, X) — X = feature data frame (no Y)
#   Ranger:   ranger.unify(model, X)   — X = feature data frame (no Y);
#             model must be trained with probability=TRUE and keep.inbag=TRUE

# ── §2a. XGBoost ─────────────────────────────────────────────────────────────
xgb_mat   <- xgboost::xgb.DMatrix(as.matrix(X), label = Y)
xgb_model <- xgboost::xgboost(
  data        = xgb_mat,
  nrounds     = 50L,
  max_depth   = 3L,
  eta         = 0.1,
  objective   = "binary:logistic",
  eval_metric = "logloss",
  verbose     = 0L
)
unified_xgb <- treeshap::xgboost.unify(xgb_model, X)
cat("XGBoost: trained and unified.\n")

# ── §2b. LightGBM ────────────────────────────────────────────────────────────
lgb_ds    <- lightgbm::lgb.Dataset(as.matrix(X), label = Y)
lgb_model <- lightgbm::lgb.train(
  params  = list(
    objective     = "binary",
    metric        = "binary_logloss",
    learning_rate = 0.1,
    max_depth     = 3L,
    num_leaves    = 7L,   # ≤ 2^max_depth = 8
    verbose       = -1L
  ),
  data    = lgb_ds,
  nrounds = 50L
)
unified_lgb <- treeshap::lightgbm.unify(lgb_model, X)
cat("LightGBM: trained and unified.\n")

# ── §2c. Ranger ──────────────────────────────────────────────────────────────
# Use ".Y" as outcome name to avoid any accidental clash with feature names
ranger_df    <- cbind(X, .Y = factor(Y))
ranger_model <- ranger::ranger(
  .Y ~ .,
  data         = ranger_df,
  num.trees    = 50L,
  max.depth    = 3L,
  probability  = TRUE,   # required by treeshap::ranger.unify
  keep.inbag   = TRUE,   # required by treeshap::ranger.unify
  write.forest = TRUE    # default; stated explicitly for clarity
)
unified_rng <- treeshap::ranger.unify(ranger_model, X)  # X only — no .Y column
cat("Ranger: trained and unified.\n")

# =============================================================================
# §3. treeowen() — compute Owen values
# =============================================================================
# method = "auto" (default): uses exact for |G_k| < TREEOWEN_AUTO_EXACT_MAX_M
#   (30 by default), Monte Carlo otherwise.
# All groups here have |G_k| = 5, so exact is used for every group.

cat("\n--- §3a treeowen() [XGBoost, method='auto'] ---\n")
result_xgb <- treeowen(
  unified_model = unified_xgb,
  x             = X,
  groups        = groups,
  method        = "auto",    # default; exact for all groups (|G_k|=5 < 30)
  dp_progress   = FALSE,
  verbose       = TRUE
)
print(result_xgb)

cat("\n--- §3b treeowen() [LightGBM, method='exact'] ---\n")
result_lgb <- treeowen(
  unified_model = unified_lgb,
  x             = X,
  groups        = groups,
  method        = "exact",
  dp_progress   = FALSE,
  verbose       = TRUE
)

cat("\n--- §3c treeowen() [Ranger, method='exact'] ---\n")
result_rng <- treeowen(
  unified_model = unified_rng,
  x             = X,
  groups        = groups,
  method        = "exact",
  dp_progress   = FALSE,
  verbose       = TRUE
)

# Efficiency axiom: sum_i psi_i(x^j) == f(x^j) - E[f(X)] for all j
ow_mat  <- as.matrix(result_xgb$owens)
phi_sum <- rowSums(ow_mat)
gap     <- result_xgb$v_all - result_xgb$baseline
eff_err <- abs(phi_sum - gap)
cat(sprintf("\nEfficiency check (XGBoost): max|err| = %.2e  (should be < 1e-6)\n",
            max(eff_err)))
if (max(eff_err) > 1e-3)
  warning("Efficiency error > 1e-3; consider increasing n_inner_mc.")

# =============================================================================
# §4. treeowen_importance() — feature & group importance
# =============================================================================
cat("\n--- §4 treeowen_importance() ---\n")

# type = "both": compute feature AND group importance in one call
# group_agg = "sum_abs": I_k = mean_j |Psi_k^(j)|  (paper default)
imp_xgb <- treeowen_importance(
  result    = result_xgb,
  type      = "both",
  group_agg = "sum_abs",
  sort      = TRUE,
  normalize = FALSE
)
print(imp_xgb)

cat("\nTop-10 features (mean |Owen value|):\n")
print(head(imp_xgb$feature, 10L))

cat("\nGroup importances:\n")
print(imp_xgb$group)

# Normalised: sum of importance = 1 within each level
imp_norm <- treeowen_importance(result_xgb, type = "group", normalize = TRUE)
cat("\nNormalised group importance (sums to 1):\n")
print(imp_norm$group)
stopifnot(abs(sum(imp_norm$group$importance) - 1) < 1e-10)

# =============================================================================
# §5. treeowen_beeswarm() — flat beeswarm plots
# =============================================================================
# Returns a ggplot object (level="feature" or "group") or a
# treeowen_plots list (level="both").
cat("\n--- §5 treeowen_beeswarm() ---\n")

# ── §5a. Feature-level beeswarm ───────────────────────────────────────────────
p_feat <- treeowen_beeswarm(
  result        = result_xgb,
  level         = "feature",   # one row per feature
  top_n_feature = 15L,         # show top 15 features by importance
  color_low     = "#0052A5",   # blue = low feature value
  color_high    = "#DC2626",   # red  = high feature value
  xlab_feature  = "Owen Value (XGBoost)",
  verbose       = FALSE
)
print(p_feat)   # ggplot object

# ── §5b. Group-level beeswarm ────────────────────────────────────────────────
p_grp <- treeowen_beeswarm(
  result           = result_xgb,
  level            = "group",   # one row per group
  top_n_group      = 5L,
  group_color_stat = "sum",     # colour = sum of feature values in group
  xlab_group       = "Group Owen Value (XGBoost)",
  verbose          = FALSE
)
print(p_grp)

# ── §5c. Both panels (group | connector | feature) ───────────────────────────
# Requires patchwork.  Returns a treeowen_plots list.
out_both <- treeowen_beeswarm(
  result          = result_xgb,
  level           = "both",
  top_n_feature   = 10L,
  top_n_group     = 5L,
  connector_width = 0.10,   # relative width of the arrow connector panel
  verbose         = FALSE
)
print(out_both$combined)               # patchwork figure
if (!is.null(out_both$feature)) print(out_both$feature)
if (!is.null(out_both$group))   print(out_both$group)

# =============================================================================
# §6. treeowen_hierarchical_beeswarm() — hierarchical beeswarm
# =============================================================================
# imp must be computed with type = "both" (needs both $feature and $group).
cat("\n--- §6 treeowen_hierarchical_beeswarm() ---\n")

imp_lgb <- treeowen_importance(result_lgb, type = "both")

# ── §6a. Return patchwork list (no files saved) ───────────────────────────────
pages <- treeowen_hierarchical_beeswarm(
  ow_result        = result_lgb,
  imp              = imp_lgb,
  top_n_group      = 5L,     # top 5 groups by importance
  top_n_feat       = NULL,   # all features within those groups (5 × 5 = 25 rows)
  n_col            = 2L,     # 2-column layout
  group_color_stat = "sum",  # colour = sum of feature values
  show_colorbar    = TRUE,   # append standalone colorbar strip at bottom
  colorbar_h       = 0.35,   # colorbar height (inches)
  row_h_grp        = 0.55,   # height per group-header row (inches)
  row_h_feat       = 0.22,   # height per feature row (inches)
  width_in         = 10,
  max_h_per_page   = 15,     # split into multiple pages if > 15 in tall
  verbose          = TRUE
)
cat(sprintf("Pages returned: %d\n", length(pages)))
stopifnot(length(pages) >= 1L)
print(pages[[1]])

# ── §6b. Save PDF + PNG ───────────────────────────────────────────────────────
# lname sets the file-name prefix; save_path sets the target directory.
tmp_dir   <- tempdir()
invisible(
  treeowen_hierarchical_beeswarm(
    ow_result   = result_lgb,
    imp         = imp_lgb,
    top_n_group = 5L,
    n_col       = 2L,
    width_in    = 10,
    lname       = "lgb_synthetic",   # → fig_lgb_synthetic_beeswarm_hierarchical.{pdf,png}
    save_path   = tmp_dir,
    dpi         = 150L,
    verbose     = TRUE
  )
)
saved_pdf <- file.path(tmp_dir, "fig_lgb_synthetic_beeswarm_hierarchical.pdf")
stopifnot(file.exists(saved_pdf))
cat(sprintf("Saved: %s\n", saved_pdf))

# =============================================================================
# §7. treeowen_exact_enum() — model-agnostic reference enumerator
#     WARNING: exponential in K and max|G_k|; feasible only for tiny problems.
# =============================================================================
cat("\n--- §7 treeowen_exact_enum() [K=2 groups, n=5 rows] ---\n")

# IMPORTANT: unified_xgb was trained on ALL 25 features.
# Subset ROWS only; keep all columns so the dp_model is consistent.
# Passing only G1+G2 column-subset would mis-specify the model.
small_x      <- X[seq_len(5L), ]        # 5 rows, all 25 columns
small_groups <- groups[c("G1", "G2")]   # only 2 groups requested

result_enum <- treeowen_exact_enum(
  unified_model = unified_xgb,
  x             = small_x,
  groups        = small_groups,
  dp_progress   = FALSE,
  verbose       = FALSE
)
print(result_enum)
# Only the G1/G2 columns will have non-trivial Owen values;
# other features (G3–G5) are not in `groups` so they are
# treated as known (conditioned on) in every coalition.

# =============================================================================
# §8. treeowen_exact_tuvalues() — tree-aware reference enumerator
# =============================================================================
cat("\n--- §8 treeowen_exact_tuvalues() [K=2 groups, n=5 rows] ---\n")

result_tu <- treeowen_exact_tuvalues(
  unified_model = unified_xgb,
  x             = small_x,
  groups        = small_groups,
  dp_progress   = FALSE,
  verbose       = FALSE
)
print(result_tu)

# Consistency check: treeowen(exact) vs treeowen_exact_enum()
result_main_small <- treeowen(
  unified_model = unified_xgb,
  x             = small_x,
  groups        = small_groups,
  method        = "exact",
  dp_progress   = FALSE,
  verbose       = FALSE
)

# Compare only the grouped features (G1 + G2)
g12_feats <- c(small_groups$G1, small_groups$G2)
v_main    <- as.numeric(as.matrix(result_main_small$owens)[, g12_feats])
v_enum    <- as.numeric(as.matrix(result_enum$owens)[, g12_feats])

# Guard against degenerate (near-zero) case before computing correlation
if (sd(v_main) > 1e-12 && sd(v_enum) > 1e-12) {
  r <- cor(v_main, v_enum)
  cat(sprintf("\nPearson r (treeowen vs exact_enum, G1+G2): %.8f  (should be > 0.9999)\n", r))
  stopifnot(r > 0.9999)
} else {
  message("Owen values near-zero (degenerate model); skipping correlation check.")
}

# =============================================================================
# §9. build_hierarchy_tree_binary() / build_hierarchy_tree_from_layers()
# =============================================================================
# The auxiliary tree B over groups is used for computational efficiency only;
# it does not affect the theoretical correctness of the Owen values.
# hierarchy = NULL (default in treeowen()) always produces a balanced tree.
cat("\n--- §9 build_hierarchy_tree_*() ---\n")

K <- length(groups)   # 5

# build_hierarchy_tree_from_layers(): main entry point called by treeowen()
# hierarchy = NULL → balanced binary tree (lowest max-depth for K groups)
tree_lay <- build_hierarchy_tree_from_layers(K, hierarchy = NULL)
cat("Balanced tree via build_hierarchy_tree_from_layers(K, NULL):\n")
str(tree_lay, max.level = 2)

# build_hierarchy_tree_binary(): low-level builder for nested-list trees.
# Accepts a NESTED LIST where each leaf is a single integer group index and
# each internal node is list(left = ..., right = ...).
# Example: ((G1, G2), (G3, (G4, G5)))
tree_custom_list <- build_hierarchy_tree_binary(
  list(
    left  = list(left = 1L, right = 2L),   # meta-group {G1, G2}
    right = list(
      left  = 3L,                            # G3
      right = list(left = 4L, right = 5L)   # meta-group {G4, G5}
    )
  )
)
cat("\nCustom nested-list tree via build_hierarchy_tree_binary():\n")
str(tree_custom_list, max.level = 2)

# treeowen() with default hierarchy (NULL → balanced)
result_bal <- treeowen(
  unified_model = unified_xgb,
  x             = X[seq_len(5L), ],
  groups        = groups,
  method        = "exact",
  hierarchy     = NULL,   # balanced binary tree (default)
  dp_progress   = FALSE,
  verbose       = FALSE
)
cat("treeowen() with default hierarchy: OK\n")

# treeowen() with custom nested-list hierarchy
result_custom <- treeowen(
  unified_model = unified_xgb,
  x             = X[seq_len(5L), ],
  groups        = groups,
  method        = "exact",
  hierarchy     = list(          # ((G1,G2), (G3,(G4,G5)))
    left  = list(left = 1L, right = 2L),
    right = list(left = 3L, right = list(left = 4L, right = 5L))
  ),
  dp_progress   = FALSE,
  verbose       = FALSE
)
cat("treeowen() with custom nested-list hierarchy: OK\n")

# Both hierarchies give the same Owen values (correctness is invariant to B)
max_diff_hier <- max(abs(as.matrix(result_bal$owens) -
                         as.matrix(result_custom$owens)))
cat(sprintf("Max |balanced - custom| Owen diff: %.2e  (should be ~0)\n",
            max_diff_hier))

# =============================================================================
# §10. clear_inner_enum_cache()
# =============================================================================
cat("\n--- §10 clear_inner_enum_cache() ---\n")

# Populate the cache by running treeowen()
invisible(treeowen(
  unified_model = unified_xgb,
  x             = X[seq_len(3L), ],
  groups        = groups,
  method        = "exact",
  dp_progress   = FALSE,
  verbose       = FALSE
))

# Access the cache environment via the package namespace (safer than :::)
pkg_ns    <- getNamespace("treeowen")
cache_env <- pkg_ns$.INNER_ENUM_CACHE
if (is.null(cache_env))
  stop("Cannot access .INNER_ENUM_CACHE — check package installation")

n_before <- length(ls(cache_env))
cat(sprintf("Cache entries before clearing: %d\n", n_before))

clear_inner_enum_cache()

n_after <- length(ls(cache_env))
cat(sprintf("Cache entries after  clearing: %d\n", n_after))
stopifnot(n_after == 0L)

# =============================================================================
# §11. Multicore (n_cores > 1) — Linux/macOS only
#      Uses parallel::mclapply(); Windows requires n_cores = 1L.
# =============================================================================
cat("\n--- §11 Multicore ---\n")

n_phys  <- max(1L, parallel::detectCores(logical = FALSE))
n_cores <- if (.Platform$OS.type == "unix" && n_phys >= 2L) 2L else 1L

if (n_cores > 1L) {
  cat(sprintf("Physical cores: %d — running with n_cores = %d\n",
              n_phys, n_cores))
  result_par <- treeowen(
    unified_model = unified_xgb,
    x             = X,
    groups        = groups,
    method        = "exact",
    n_cores       = n_cores,
    dp_progress   = FALSE,
    verbose       = TRUE
  )
  # Exact algorithms produce bit-for-bit identical results across cores;
  # floating-point may differ by up to ~1e-14 due to summation order.
  max_diff <- max(abs(as.matrix(result_par$owens) -
                      as.matrix(result_xgb$owens)))
  cat(sprintf("Max |parallel - serial| Owen diff: %.2e  (should be < 1e-8)\n",
              max_diff))
  stopifnot(max_diff < 1e-8)
  cat("Parallel result matches serial result. OK\n")
} else {
  cat(sprintf(
    "Skipping parallel test (OS: %s, physical cores: %d).\n",
    .Platform$OS.type, n_phys))
  cat("On Linux/macOS with >= 2 cores, set n_cores = 2L (or more) in treeowen().\n")
}

# =============================================================================
# Summary
# =============================================================================
cat("\n=============================================================\n")
cat("All treeowen example sections completed successfully.\n")
cat("=============================================================\n")
