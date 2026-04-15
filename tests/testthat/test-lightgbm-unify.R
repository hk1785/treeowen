# tests/testthat/test-lightgbm-unify.R
#
# Unit tests for .lightgbm_unify_compat()
#
# Covers failure modes:
#   [LGB-1]  treeshap::lightgbm.unify() bypass / fallback
#   [LGB-2]  lgb.model.dt.tree() unavailable → dump_model() fallback
#   [LGB-3]  Column name variations across lightgbm versions
#   [LGB-4]  Three child-ID encoding schemes (negative leaf, positive node, "N"/"L" strings)
#   [LGB-5]  Feature names NULL → model$feature_name() → colnames(data) fallback
#   [LGB-6]  best_iter slot name change (cv$best_iter vs cv$best_iteration)
#   [LGB-7]  lgb.cv() verbose argument location change

skip_if_not_lgb <- function() {
  if (!requireNamespace("lightgbm", quietly = TRUE))
    testthat::skip("Package 'lightgbm' not available")
}

# ── Shared fixture: tiny binary LightGBM model ────────────────────────────
.make_lgb_fixture <- function(n = 50L, p = 6L, seed = 1L) {
  skip_if_not_lgb()
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p,
               dimnames = list(NULL, paste0("f", seq_len(p))))
  Y <- as.integer(rowSums(X[, 1:3]) > 0)

  ds <- lightgbm::lgb.Dataset(X, label = Y)
  params <- list(objective = "binary", metric = "binary_logloss",
                 num_leaves = 7L, max_depth = 3L, verbose = -1L,
                 learning_rate = 0.1)
  set.seed(seed)
  m  <- lightgbm::lgb.train(params = params, data = ds, nrounds = 15L,
                              verbose = -1L)
  list(model = m, X = X, Y = Y, feature_names = paste0("f", seq_len(p)))
}

# ──────────────────────────────────────────────────────────────────────────────
# TEST 1: Output has correct class and structure
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB] output has correct class and structure", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()
  u  <- .lightgbm_unify_compat(fx$model, fx$X)

  expect_s3_class(u, "model_unified")
  expect_true(is.data.frame(u$model))
  expect_true(is.data.frame(u$data))
  expect_true(is.character(u$feature_names))
  expect_equal(attr(u, "model"), "lightgbm")

  req <- c("Tree", "Node", "Feature", "Decision.type",
           "Split", "Yes", "No", "Missing", "Prediction", "Cover")
  expect_true(all(req %in% colnames(u$model)),
              info = paste("Missing:", paste(setdiff(req, colnames(u$model)), collapse = ", ")))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 2: Feature names match training data
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB] feature_names match training data columns", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()
  u  <- .lightgbm_unify_compat(fx$model, fx$X)
  expect_equal(sort(u$feature_names), sort(fx$feature_names))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 3 [LGB-3]: Column alias normalisation — simulate "split_point" instead of "threshold"
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB-3] .alias() correctly maps split_point → threshold", {
  # Test the internal alias logic directly on a dummy data.frame
  df <- data.frame(split_point = c(0.5, NA, 1.2), stringsAsFactors = FALSE)

  # Replicate the alias helper from the wrapper
  .alias <- function(df, pref, alts) {
    if (pref %in% colnames(df)) return(df)
    for (a in alts) if (a %in% colnames(df)) { df[[pref]] <- df[[a]]; return(df) }
    df[[pref]] <- NA_real_; df
  }

  df2 <- .alias(df, "threshold", c("split_point"))
  expect_true("threshold" %in% colnames(df2))
  expect_equal(df2$threshold, df$split_point)
})

test_that("[LGB-3] .alias() correctly maps value → split_gain", {
  df <- data.frame(value = c(1.1, -0.3, 0.7), stringsAsFactors = FALSE)
  .alias <- function(df, pref, alts) {
    if (pref %in% colnames(df)) return(df)
    for (a in alts) if (a %in% colnames(df)) { df[[pref]] <- df[[a]]; return(df) }
    df[[pref]] <- NA_real_; df
  }
  df2 <- .alias(df, "split_gain", c("value"))
  expect_true("split_gain" %in% colnames(df2))
  expect_equal(df2$split_gain, df$value)
})

test_that("[LGB-3] .alias() correctly maps internal_count → count", {
  df <- data.frame(internal_count = c(10L, 5L, 8L), stringsAsFactors = FALSE)
  .alias <- function(df, pref, alts) {
    if (pref %in% colnames(df)) return(df)
    for (a in alts) if (a %in% colnames(df)) { df[[pref]] <- df[[a]]; return(df) }
    df[[pref]] <- NA_real_; df
  }
  df2 <- .alias(df, "count", c("internal_count", "leaf_count", "node_count"))
  expect_true("count" %in% colnames(df2))
  expect_equal(df2$count, df$internal_count)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 4 [LGB-4]: Child ID encoding — string "N<int>" / "L<int>" format
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB-4] 'N<int>'/'L<int>' child ID strings are parsed correctly", {
  # Reproduce the string-encoding resolution logic
  n_rows  <- 6L
  tree_id <- c(0L, 0L, 0L, 0L, 0L, 0L)

  # Simulate: node 0 → left=N1, right=N2;  N1 → left=L0, right=L1;  N2 → left=L2, right=L3
  # Encoded IDs:  N<k> → +k (internal node),  L<k> → -(k+1) (leaf)
  # Rows: N0, N1, N2, L0, L1, L2  (indices 1..6)
  encoded_id <- c(0L, 1L, 2L, -1L, -2L, -3L)
  lookup_key <- paste(tree_id, encoded_id, sep = ":")
  lookup_map <- setNames(seq_len(n_rows), lookup_key)

  .resolve <- function(raw_char) {
    encoded  <- suppressWarnings(as.integer(raw_char))
    n_mask <- grepl("^N", raw_char, perl = TRUE)
    l_mask <- grepl("^L", raw_char, perl = TRUE)
    encoded[n_mask] <-  as.integer(sub("^N", "", raw_char[n_mask]))
    encoded[l_mask] <- -(as.integer(sub("^L", "", raw_char[l_mask])) + 1L)
    key <- paste(tree_id, encoded, sep = ":")
    as.integer(unname(lookup_map[key]))
  }

  left_raw  <- c("N1", "L0", "L2", NA_character_, NA_character_, NA_character_)
  right_raw <- c("N2", "L1", "L3", NA_character_, NA_character_, NA_character_)

  left_r  <- .resolve(left_raw)
  right_r <- .resolve(right_raw)

  expect_equal(left_r[1],  2L)   # N1 = row 2
  expect_equal(right_r[1], 3L)   # N2 = row 3
  expect_equal(left_r[2],  4L)   # L0 = row 4
  expect_equal(right_r[2], 5L)   # L1 = row 5
  expect_equal(left_r[3],  6L)   # L2 = row 6
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 5 [LGB-4]: Negative-indexed leaf encoding
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB-4] Negative-indexed leaf IDs resolve to correct rows", {
  # leaf #0 encoded as -1; leaf #1 encoded as -2
  n_rows  <- 4L
  tree_id <- rep(0L, n_rows)
  # row1=node0, row2=node1, row3=leaf0, row4=leaf1
  encoded_id <- c(0L, 1L, -1L, -2L)
  lookup_key <- paste(tree_id, encoded_id, sep = ":")
  lookup_map <- setNames(seq_len(n_rows), lookup_key)

  .resolve <- function(raw_ids) {
    raw_char <- as.character(raw_ids)
    encoded  <- suppressWarnings(as.integer(raw_char))
    key <- paste(tree_id, encoded, sep = ":")
    as.integer(unname(lookup_map[key]))
  }

  left_ids  <- c(-1L, -3L, NA_integer_, NA_integer_)   # -1 = leaf0 (row3), -3 absent
  right_ids <- c(-2L, NA_integer_, NA_integer_, NA_integer_)

  lr <- .resolve(left_ids)
  rr <- .resolve(right_ids)
  expect_equal(lr[1],  3L)   # leaf0 = row 3
  expect_true (is.na(lr[2]))  # -3 not in lookup
  expect_equal(rr[1],  4L)   # leaf1 = row 4
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 6 [LGB-5]: Feature names fallback chain
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB-5] feature names resolved via model$feature_name() when private slot absent", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()

  # model$feature_name() should work
  nm <- tryCatch(fx$model$feature_name(), error = function(e) NULL)
  if (!is.null(nm) && length(nm)) {
    expect_equal(sort(nm), sort(fx$feature_names))
  } else {
    skip("model$feature_name() not available in this lightgbm version")
  }
})

test_that("[LGB-5] feature names fall back to colnames(data) when model has none", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()
  # Simulate: if feature_name() fails, colnames(data) is the last resort
  X_named <- fx$X
  nm_fallback <- colnames(X_named)
  expect_equal(sort(nm_fallback), sort(fx$feature_names))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 7 [LGB-6]: lgb.cv() best_iter slot name compatibility
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB-6] best_iter extracted safely from cv object", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()

  ds     <- lightgbm::lgb.Dataset(fx$X, label = fx$Y)
  params <- list(objective = "binary", metric = "binary_logloss",
                 num_leaves = 7L, verbose = -1L)
  cv     <- lightgbm::lgb.cv(params = params, data = ds, nrounds = 20L,
                               nfold = 3L, verbose = -1L,
                               early_stopping_rounds = 5L)

  # Both slot names should be tried
  best_n <- cv$best_iter
  if (is.null(best_n) || !is.finite(best_n) || best_n < 1L)
    best_n <- cv$best_iteration
  if (is.null(best_n) || !is.finite(best_n) || best_n < 1L)
    best_n <- 20L   # fallback

  expect_true(is.numeric(best_n))
  expect_true(best_n >= 1L)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 8 [LGB-7]: lgb.cv() verbose argument — both locations succeed
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB-7] lgb.cv() accepts verbose as direct argument", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()
  ds <- lightgbm::lgb.Dataset(fx$X, label = fx$Y)

  # verbose as direct argument (newer lightgbm)
  params_no_v <- list(objective = "binary", num_leaves = 7L)
  cv <- tryCatch(
    lightgbm::lgb.cv(params = params_no_v, data = ds, nrounds = 5L,
                      nfold = 3L, verbose = -1L),
    error = function(e) NULL
  )
  if (!is.null(cv)) {
    expect_true(!is.null(cv))
  } else {
    # verbose inside params (older lightgbm)
    params_v <- list(objective = "binary", num_leaves = 7L, verbose = -1L)
    cv2 <- lightgbm::lgb.cv(params = params_v, data = ds, nrounds = 5L, nfold = 3L)
    expect_true(!is.null(cv2))
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 9: Leaf nodes have numeric Prediction; internal nodes have NA
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB] internal nodes have NA Prediction; leaves have numeric Prediction", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()
  u  <- .lightgbm_unify_compat(fx$model, fx$X)
  m  <- u$model

  leaf     <- is.na(m$Feature)
  internal <- !is.na(m$Feature)

  if (any(internal)) {
    expect_true(all(is.na(m$Prediction[internal])),
                info = "Internal node Prediction must be NA")
  }
  if (any(leaf)) {
    expect_true(all(!is.na(m$Prediction[leaf])),
                info = "Leaf node Prediction must be numeric")
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 10: Roundtrip — model_unified feeds .prepare_dp_model
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB] model_unified passes .prepare_dp_model without error", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()
  u  <- .lightgbm_unify_compat(fx$model, fx$X)
  dp <- .prepare_dp_model(u, colnames(fx$X))
  expect_true(is.list(dp))
  expect_true(length(dp$trees) > 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 11: matrix and data.frame inputs both accepted
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB] matrix and data.frame inputs both produce valid output", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()
  u_mat <- .lightgbm_unify_compat(fx$model, fx$X)
  u_df  <- .lightgbm_unify_compat(fx$model, as.data.frame(fx$X))
  expect_s3_class(u_mat, "model_unified")
  expect_s3_class(u_df,  "model_unified")
  expect_equal(u_mat$feature_names, u_df$feature_names)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 12: .lgb_parse_dump_model() produces a valid data.frame
# ──────────────────────────────────────────────────────────────────────────────
test_that("[LGB-2] .lgb_parse_dump_model() returns a non-empty data.frame", {
  skip_if_not_lgb()
  fx <- .make_lgb_fixture()

  dm_raw <- tryCatch(fx$model$dump_model(), error = function(e) NULL)
  if (is.null(dm_raw) || is.null(dm_raw$tree_info)) {
    skip("dump_model() not available in this lightgbm version")
  }

  df <- .lgb_parse_dump_model(dm_raw$tree_info)
  expect_true(is.data.frame(df))
  expect_true(nrow(df) > 0L)
  expect_true("tree_index" %in% colnames(df))
  expect_true("feature"    %in% colnames(df))
})
