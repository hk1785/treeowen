# tests/testthat/test-treeowen-core.R
#
# Integration tests: treeowen() end-to-end with all three learners.
# Each test runs a minimal train → unify → treeowen() → check efficiency pipeline.
# Also tests group definition, importance, and result structure.

skip_if_no_cpp <- function() {
  if (!isTRUE(.try_load_cpp())) testthat::skip("C++ backend not available")
}
skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    testthat::skip(sprintf("Package '%s' not available", pkg))
}

# ── Shared synthetic data ─────────────────────────────────────────────────
.make_core_data <- function(n = 40L, p = 10L, seed = 42L) {
  set.seed(seed)
  X <- as.data.frame(matrix(rnorm(n * p), n, p,
                              dimnames = list(NULL, paste0("x", seq_len(p)))))
  Y <- as.integer(rowSums(X[, 1:3]) > 0)
  # 5 groups of 2 features each
  groups <- list(
    G1 = c("x1",  "x2"),
    G2 = c("x3",  "x4"),
    G3 = c("x5",  "x6"),
    G4 = c("x7",  "x8"),
    G5 = c("x9",  "x10")
  )
  list(X = X, Y = Y, groups = groups, p = p, n = n)
}

# ──────────────────────────────────────────────────────────────────────────────
# TEST 1: XGBoost end-to-end — efficiency axiom holds
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core-XGB] treeowen() efficiency axiom holds for XGBoost", {
  skip_if_no_cpp()
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  d  <- .make_core_data()
  dm <- xgboost::xgb.DMatrix(as.matrix(d$X), label = d$Y)
  set.seed(42)
  m  <- xgboost::xgboost(data = dm, nrounds = 20L, max_depth = 2L,
                          objective = "binary:logistic", verbose = 0)
  u  <- xgboost_unify_compat(m, d$X)

  res <- treeowen(u, d$X, d$groups, method = "exact",
                  dp_progress = FALSE, verbose = FALSE)

  expect_s3_class(res, "treeowen_result")
  expect_equal(nrow(res$owens), d$n)
  expect_equal(ncol(res$owens), d$p)

  phi_sum <- rowSums(as.matrix(res$owens))
  gap     <- res$v_all - res$baseline
  eff_err <- abs(phi_sum - gap)
  expect_true(max(eff_err, na.rm = TRUE) < 1e-4,
              info = sprintf("Max efficiency error: %.2e", max(eff_err, na.rm = TRUE)))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 2: LightGBM end-to-end — efficiency axiom holds
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core-LGB] treeowen() efficiency axiom holds for LightGBM", {
  skip_if_no_cpp()
  skip_if_not_installed("lightgbm")

  d  <- .make_core_data()
  ds <- lightgbm::lgb.Dataset(as.matrix(d$X), label = d$Y)
  params <- list(objective = "binary", num_leaves = 7L, verbose = -1L,
                 learning_rate = 0.1)
  set.seed(42)
  m  <- lightgbm::lgb.train(params = params, data = ds, nrounds = 20L,
                              verbose = -1L)
  u  <- lightgbm_unify_compat(m, d$X)

  res <- treeowen(u, d$X, d$groups, method = "exact",
                  dp_progress = FALSE, verbose = FALSE)

  expect_s3_class(res, "treeowen_result")
  phi_sum <- rowSums(as.matrix(res$owens))
  gap     <- res$v_all - res$baseline
  eff_err <- abs(phi_sum - gap)
  expect_true(max(eff_err, na.rm = TRUE) < 1e-4,
              info = sprintf("Max efficiency error: %.2e", max(eff_err, na.rm = TRUE)))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 3: Ranger end-to-end — efficiency axiom holds
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core-RNG] treeowen() efficiency axiom holds for Ranger", {
  skip_if_no_cpp()
  skip_if_not_installed("ranger")

  d  <- .make_core_data()
  df <- cbind(d$X, .y = factor(d$Y, c(0,1), c("neg","pos")))
  set.seed(42)
  m  <- ranger::ranger(.y ~ ., data = df, num.trees = 20L, max.depth = 3L,
                        probability = TRUE, keep.inbag = TRUE,
                        write.forest = TRUE, seed = 42L)
  u  <- ranger_unify_compat(m, d$X)

  res <- treeowen(u, d$X, d$groups, method = "exact",
                  dp_progress = FALSE, verbose = FALSE)

  expect_s3_class(res, "treeowen_result")
  phi_sum <- rowSums(as.matrix(res$owens))
  gap     <- res$v_all - res$baseline
  eff_err <- abs(phi_sum - gap)
  expect_true(max(eff_err, na.rm = TRUE) < 1e-4,
              info = sprintf("Max efficiency error: %.2e", max(eff_err, na.rm = TRUE)))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 4: result structure — all required fields present
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core] treeowen_result has all required fields", {
  skip_if_no_cpp()
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  d  <- .make_core_data()
  dm <- xgboost::xgb.DMatrix(as.matrix(d$X), label = d$Y)
  set.seed(1)
  m  <- xgboost::xgboost(data = dm, nrounds = 10L, max_depth = 2L,
                          objective = "binary:logistic", verbose = 0)
  u  <- xgboost_unify_compat(m, d$X)
  res <- treeowen(u, d$X, d$groups, method = "exact",
                  dp_progress = FALSE, verbose = FALSE)

  expect_true(is.data.frame(res$owens))
  expect_true(is.numeric(res$baseline))
  expect_true(is.numeric(res$v_all))
  expect_true(is.list(res$stats))
  expect_true(is.character(res$note))
  expect_equal(length(res$baseline), d$n)
  expect_equal(length(res$v_all),    d$n)

  req_stats <- c("n", "p", "K", "trees", "method_effective",
                 "n_exact_groups", "n_approx_groups")
  expect_true(all(req_stats %in% names(res$stats)),
              info = paste("Missing stats:", paste(setdiff(req_stats, names(res$stats)), collapse = ", ")))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 5: treeowen_importance() returns correct structure
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core] treeowen_importance() returns correct structure", {
  skip_if_no_cpp()
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  d  <- .make_core_data()
  dm <- xgboost::xgb.DMatrix(as.matrix(d$X), label = d$Y)
  set.seed(1)
  m  <- xgboost::xgboost(data = dm, nrounds = 10L, max_depth = 2L,
                          objective = "binary:logistic", verbose = 0)
  u   <- xgboost_unify_compat(m, d$X)
  res <- treeowen(u, d$X, d$groups, method = "exact",
                  dp_progress = FALSE, verbose = FALSE)
  imp <- treeowen_importance(res, type = "both", sort = TRUE, normalize = FALSE)

  expect_true(is.data.frame(imp$feature))
  expect_true(is.data.frame(imp$group))
  expect_true("feature"    %in% colnames(imp$feature))
  expect_true("importance" %in% colnames(imp$feature))
  expect_true("group"      %in% colnames(imp$group))
  expect_true("importance" %in% colnames(imp$group))
  expect_equal(nrow(imp$feature), d$p)
  expect_equal(nrow(imp$group),   length(d$groups))

  # Importances must be non-negative
  expect_true(all(imp$feature$importance >= 0))
  expect_true(all(imp$group$importance   >= 0))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 6: normalize=TRUE sums to 1
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core] treeowen_importance(normalize=TRUE) group importances sum to 1", {
  skip_if_no_cpp()
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  d  <- .make_core_data()
  dm <- xgboost::xgb.DMatrix(as.matrix(d$X), label = d$Y)
  set.seed(1)
  m  <- xgboost::xgboost(data = dm, nrounds = 10L, max_depth = 2L,
                          objective = "binary:logistic", verbose = 0)
  u   <- xgboost_unify_compat(m, d$X)
  res <- treeowen(u, d$X, d$groups, method = "exact",
                  dp_progress = FALSE, verbose = FALSE)
  imp <- treeowen_importance(res, type = "group", normalize = TRUE)

  expect_equal(sum(imp$group$importance), 1, tolerance = 1e-10)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 7: method="approx" runs without error
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core] treeowen() method='approx' completes without error", {
  skip_if_no_cpp()
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  d  <- .make_core_data(n = 20L)
  dm <- xgboost::xgb.DMatrix(as.matrix(d$X), label = d$Y)
  set.seed(1)
  m  <- xgboost::xgboost(data = dm, nrounds = 10L, max_depth = 2L,
                          objective = "binary:logistic", verbose = 0)
  u  <- xgboost_unify_compat(m, d$X)

  expect_no_error(
    treeowen(u, d$X, d$groups, method = "approx",
             n_inner_mc = 16L, min_inner_mc = 8L, max_inner_mc = 32L,
             dp_progress = FALSE, verbose = FALSE)
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 8: method="auto" uses exact for small groups, approx for large
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core] treeowen() method='auto' splits methods correctly", {
  skip_if_no_cpp()
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  d  <- .make_core_data()
  dm <- xgboost::xgb.DMatrix(as.matrix(d$X), label = d$Y)
  set.seed(1)
  m  <- xgboost::xgboost(data = dm, nrounds = 10L, max_depth = 2L,
                          objective = "binary:logistic", verbose = 0)
  u  <- xgboost_unify_compat(m, d$X)

  # All groups have size 2 < auto_exact_max_m=3 → all exact
  res <- treeowen(u, d$X, d$groups, method = "auto", auto_exact_max_m = 3L,
                  dp_progress = FALSE, verbose = FALSE)
  expect_equal(res$stats$method_effective, "exact")
  expect_equal(res$stats$n_approx_groups, 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 9: print.treeowen_result does not error
# ──────────────────────────────────────────────────────────────────────────────
test_that("[Core] print.treeowen_result runs without error", {
  skip_if_no_cpp()
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  d  <- .make_core_data(n = 15L)
  dm <- xgboost::xgb.DMatrix(as.matrix(d$X), label = d$Y)
  set.seed(1)
  m  <- xgboost::xgboost(data = dm, nrounds = 5L, max_depth = 2L,
                          objective = "binary:logistic", verbose = 0)
  u   <- xgboost_unify_compat(m, d$X)
  res <- treeowen(u, d$X, d$groups, method = "exact",
                  dp_progress = FALSE, verbose = FALSE)

  expect_output(print(res), regexp = "TreeOwen result")
})
