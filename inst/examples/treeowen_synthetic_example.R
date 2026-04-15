###############################################################################
# treeowen — Complete Example
# Three learners: XGBoost, LightGBM, Ranger
###############################################################################

library(treeowen)
library(treeshap)

# ── Packages (install once if needed) ─────────────────────────────────────────
# install.packages(c("xgboost", "lightgbm", "ranger", "treeshap",
#                    "ggplot2", "ggbeeswarm", "patchwork"))

###############################################################################
# 1. Synthetic data
###############################################################################
set.seed(42)

# 100 features in 20 groups of 5
feat_names <- paste0("F", rep(1:5, times = 20), "G", rep(1:20, each = 5))
X <- as.data.frame(
  matrix(rnorm(100 * 100), 100, 100, dimnames = list(NULL, feat_names))
)

# Binary outcome driven mainly by G1 and G3
log_odds <- 0.8 * X$F1G1 - 0.6 * X$F2G1 + 0.5 * X$F1G3 + 0.3 * X$F1G2
Y <- as.integer(log_odds > 0)

# Feature partition: 20 groups of 5 features each
groups <- setNames(
  lapply(1:20, function(k) paste0("F", 1:5, "G", k)),
  paste0("G", 1:20)
)

cat(sprintf("Data: n=%d  p=%d  K=%d  Pr(Y=1)=%.2f\n",
            nrow(X), ncol(X), length(groups), mean(Y)))


###############################################################################
# 2. XGBoost
###############################################################################
library(xgboost)

cat("\n── XGBoost ──────────────────────────────────────────────────────\n")

dm        <- xgb.DMatrix(as.matrix(X), label = Y)
model_xgb <- xgboost(dm, nrounds = 100, max_depth = 3, eta = 0.1,
                     objective = "binary:logistic", verbose = 0)

# Unify (robust wrapper handles all xgboost version differences)
unified_xgb <- xgboost_unify_compat(model_xgb, X)

# Compute Owen values
result_xgb <- treeowen(
  unified_model = unified_xgb,
  x             = X,
  groups        = groups,
  method        = "auto",   # exact for |G_k| < 30, MC otherwise
  dp_progress   = FALSE,
  verbose       = TRUE
)
print(result_xgb)

# Efficiency check: sum(Owen values) == f(x) - E[f(X)] for every row
eff_err <- abs(rowSums(result_xgb$owens) -
                 (result_xgb$v_all - result_xgb$baseline))
cat(sprintf("XGBoost efficiency max error: %.2e\n", max(eff_err)))

# Importance
imp_xgb <- treeowen_importance(result_xgb, type = "both",
                               group_agg = "sum_abs", sort = TRUE)
cat("\nTop 5 features:\n"); print(head(imp_xgb$feature, 5))
cat("\nGroup importance:\n");  print(imp_xgb$group)

# Beeswarm plots
p_feat_xgb <- treeowen_beeswarm(result_xgb, level = "feature",
                                top_n_feature = 10L)
p_grp_xgb  <- treeowen_beeswarm(result_xgb, level = "group")
p_both_xgb <- treeowen_beeswarm(result_xgb, level = "both",
                                top_n_feature = 10L, top_n_group = 5L)
# print(p_feat_xgb)
# print(p_grp_xgb)
# print(p_both_xgb$combined)

# Hierarchical beeswarm
pages_xgb <- treeowen_hierarchical_beeswarm(
  ow_result     = result_xgb,
  imp           = imp_xgb,
  top_n_group   = 5L,
  top_n_feature = 100L,
  n_col         = 2L,
  verbose       = TRUE
)
# print(pages_xgb[[1]])


###############################################################################
# 3. LightGBM
###############################################################################
library(lightgbm)

cat("\n── LightGBM ─────────────────────────────────────────────────────\n")

ds        <- lgb.Dataset(as.matrix(X), label = Y)
model_lgb <- lgb.train(
  params  = list(objective     = "binary",
                 metric        = "binary_logloss",
                 learning_rate = 0.1,
                 max_depth     = 3L,
                 num_leaves    = 7L,
                 verbose       = -1L),
  data    = ds,
  nrounds = 100L,
  verbose = -1L
)

# Unify (robust wrapper: treeshap → lgb.model.dt.tree → dump_model fallback)
unified_lgb <- lightgbm_unify_compat(model_lgb, X)

# Compute Owen values
result_lgb <- treeowen(
  unified_model = unified_lgb,
  x             = X,
  groups        = groups,
  method        = "auto",
  dp_progress   = FALSE,
  verbose       = TRUE
)
print(result_lgb)

# Efficiency check
eff_err_lgb <- abs(rowSums(result_lgb$owens) -
                     (result_lgb$v_all - result_lgb$baseline))
cat(sprintf("LightGBM efficiency max error: %.2e\n", max(eff_err_lgb)))

# Importance
imp_lgb <- treeowen_importance(result_lgb, type = "both",
                               group_agg = "sum_abs", sort = TRUE)
cat("\nTop 5 features:\n"); print(head(imp_lgb$feature, 5))
cat("\nGroup importance:\n");  print(imp_lgb$group)

# Beeswarm plots
p_feat_lgb <- treeowen_beeswarm(result_lgb, level = "feature",
                                top_n_feature = 10L)
p_both_lgb <- treeowen_beeswarm(result_lgb, level = "both",
                                top_n_feature = 10L, top_n_group = 5L)
# print(p_feat_lgb)
# print(p_both_lgb$combined)

# Hierarchical beeswarm
pages_lgb <- treeowen_hierarchical_beeswarm(
  ow_result     = result_lgb,
  imp           = imp_lgb,
  top_n_group   = 5L,
  top_n_feature = 100L,
  n_col         = 2L,
  verbose       = TRUE
)
# print(pages_lgb[[1]])


###############################################################################
# 4. Ranger
###############################################################################
library(ranger)

cat("\n── Ranger ───────────────────────────────────────────────────────\n")

# Ranger requires probability = TRUE for numeric leaf predictions
df_train  <- cbind(X, .y = factor(Y, c(0, 1), c("neg", "pos")))
model_rng <- ranger(
  formula      = .y ~ .,
  data         = df_train,
  num.trees    = 100L,
  max.depth    = 3L,
  probability  = TRUE,    # required: produces probability estimates
  keep.inbag   = TRUE,    # required by some treeshap versions
  write.forest = TRUE,
  seed         = 42L
)

# Unify (robust wrapper: treeshap → treeInfo patch → forest fallback)
# Pass X without the outcome column
unified_rng <- ranger_unify_compat(model_rng, X)

# Compute Owen values
result_rng <- treeowen(
  unified_model = unified_rng,
  x             = X,
  groups        = groups,
  method        = "auto",
  dp_progress   = FALSE,
  verbose       = TRUE
)
print(result_rng)

# Efficiency check
eff_err_rng <- abs(rowSums(result_rng$owens) -
                     (result_rng$v_all - result_rng$baseline))
cat(sprintf("Ranger efficiency max error: %.2e\n", max(eff_err_rng)))

# Importance
imp_rng <- treeowen_importance(result_rng, type = "both",
                               group_agg = "sum_abs", sort = TRUE)
cat("\nTop 5 features:\n"); print(head(imp_rng$feature, 5))
cat("\nGroup importance:\n");  print(imp_rng$group)

# Beeswarm plots
p_feat_rng <- treeowen_beeswarm(result_rng, level = "feature",
                                top_n_feature = 10L)
p_both_rng <- treeowen_beeswarm(result_rng, level = "both",
                                top_n_feature = 10L, top_n_group = 5L)
# print(p_feat_rng)
# print(p_both_rng$combined)

# Hierarchical beeswarm
pages_rng <- treeowen_hierarchical_beeswarm(
  ow_result     = result_rng,
  imp           = imp_rng,
  top_n_group   = 5L,
  top_n_feature = 100L,
  n_col         = 2L,
  verbose       = TRUE
)
# print(pages_rng[[1]])


###############################################################################
# 5. Cross-learner comparison
###############################################################################
cat("\n── Cross-learner group importance comparison ─────────────────────\n")

imp_compare <- data.frame(
  group    = imp_xgb$group$group,
  xgboost  = imp_xgb$group$importance,
  lightgbm = imp_lgb$group$importance[match(imp_xgb$group$group,
                                            imp_lgb$group$group)],
  ranger   = imp_rng$group$importance[match(imp_xgb$group$group,
                                            imp_rng$group$group)]
)
imp_compare$mean <- rowMeans(imp_compare[, -1], na.rm = TRUE)
imp_compare      <- imp_compare[order(-imp_compare$mean), ]
print(imp_compare, digits = 4, row.names = FALSE)

cat("\nDone.\n")
