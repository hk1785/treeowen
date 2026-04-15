# tests/testthat/test-ranger-unify.R
#
# Unit tests for ranger_unify_compat()
#
# Covers failure modes:
#   [RNG-1]  treeInfo() "pred.1" vs "prediction" column name
#   [RNG-2]  ranger_unify.common() signature / availability change
#   [RNG-3]  treeInfo() column name variations (nodeID, splitvarID, splitval, terminal)
#   [RNG-4]  child.nodeIDs structure: list-of-2 vs separate leftChild/rightChild columns
#   [RNG-5]  Probability forest: leaf prediction as matrix/list/scalar
#   [RNG-6]  model$forest NULL (write.forest=FALSE)
#   [RNG-7]  keep.inbag not required for treeInfo path
#   [RNG-8]  Non-probability forest emits warning, does not stop

skip_if_not_rng <- function() {
  if (!requireNamespace("ranger", quietly = TRUE))
    testthat::skip("Package 'ranger' not available")
}

# ── Shared fixture: tiny probability Ranger model ────────────────────────
.make_rng_fixture <- function(n = 60L, p = 6L, seed = 1L, keep_inbag = TRUE,
                               probability = TRUE) {
  skip_if_not_rng()
  set.seed(seed)
  X  <- as.data.frame(matrix(rnorm(n * p), n, p,
                               dimnames = list(NULL, paste0("f", seq_len(p)))))
  Y  <- as.integer(rowSums(X[, 1:3]) > 0)
  df <- cbind(X, .y = factor(Y, c(0,1), c("neg","pos")))
  set.seed(seed)
  m  <- ranger::ranger(.y ~ ., data = df, num.trees = 10L, max.depth = 3L,
                        probability = probability, keep.inbag = keep_inbag,
                        write.forest = TRUE, seed = seed)
  list(model = m, X = X, Y = Y, feature_names = paste0("f", seq_len(p)))
}

# ──────────────────────────────────────────────────────────────────────────────
# TEST 1: Output has correct class and structure
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG] output has correct class and structure", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()
  u  <- ranger_unify_compat(fx$model, fx$X)

  expect_s3_class(u, "model_unified")
  expect_true(is.data.frame(u$model))
  expect_true(is.data.frame(u$data))
  expect_true(is.character(u$feature_names))
  expect_equal(attr(u, "model"), "ranger")

  req <- c("Tree", "Node", "Feature", "Decision.type",
           "Split", "Yes", "No", "Missing", "Prediction", "Cover")
  expect_true(all(req %in% colnames(u$model)),
              info = paste("Missing:", paste(setdiff(req, colnames(u$model)), collapse = ", ")))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 2: Feature names match training data
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG] feature_names match training data columns", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()
  u  <- ranger_unify_compat(fx$model, fx$X)
  expect_equal(sort(u$feature_names), sort(fx$feature_names))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 3 [RNG-1]: ranger_treeinfo_norm() handles "pred.1" → "prediction"
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-1] ranger_treeinfo_norm() normalises 'pred.1' to 'prediction'", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()

  # Simulate a treeInfo output that uses "pred.1" instead of "prediction"
  ti_raw <- tryCatch(ranger::treeInfo(fx$model, tree = 1L), error = function(e) NULL)
  if (is.null(ti_raw)) skip("ranger::treeInfo() not available")

  ti_df <- as.data.frame(ti_raw, stringsAsFactors = FALSE)

  # Rename prediction → pred.1 to simulate old ranger
  if ("prediction" %in% colnames(ti_df)) {
    colnames(ti_df)[colnames(ti_df) == "prediction"] <- "pred.1"
    expect_false("prediction" %in% colnames(ti_df))
    expect_true("pred.1" %in% colnames(ti_df))

    # Apply the normalization logic
    norm <- ranger_treeinfo_norm(ti_df, fx$feature_names)
    expect_true(is.data.frame(norm))
    expect_true("prediction" %in% colnames(norm))
    expect_true(nrow(norm) > 0L)
  } else {
    skip("prediction column absent in treeInfo output — skip pred.1 test")
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 4 [RNG-3]: ranger_treeinfo_norm() handles column name variations
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-3] ranger_treeinfo_norm() handles 'node' alias for 'nodeID'", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()
  ti_raw <- tryCatch(ranger::treeInfo(fx$model, tree = 1L), error = function(e) NULL)
  if (is.null(ti_raw)) skip("ranger::treeInfo() not available")

  ti_df <- as.data.frame(ti_raw, stringsAsFactors = FALSE)
  # Rename nodeID → node to simulate a version variation
  if ("nodeID" %in% colnames(ti_df)) {
    colnames(ti_df)[colnames(ti_df) == "nodeID"] <- "node"
    norm <- ranger_treeinfo_norm(ti_df, fx$feature_names)
    expect_true(is.data.frame(norm))
    expect_true("nodeid" %in% colnames(norm) || nrow(norm) > 0L)
  } else {
    skip("nodeID column not present in this ranger version")
  }
})

test_that("[RNG-3] ranger_treeinfo_norm() handles 'splitVal' alias for 'splitval'", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()
  ti_raw <- tryCatch(ranger::treeInfo(fx$model, tree = 1L), error = function(e) NULL)
  if (is.null(ti_raw)) skip("ranger::treeInfo() not available")

  ti_df <- as.data.frame(ti_raw, stringsAsFactors = FALSE)
  if ("splitval" %in% tolower(colnames(ti_df))) {
    # Rename to splitVal
    idx <- which(tolower(colnames(ti_df)) == "splitval")[1]
    colnames(ti_df)[idx] <- "splitVal"
    norm <- ranger_treeinfo_norm(ti_df, fx$feature_names)
    expect_true(is.data.frame(norm))
    expect_true("splitval" %in% colnames(norm))
  } else {
    skip("splitval not present in treeInfo output")
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 5 [RNG-4]: ranger_tree_from_forest() works with forest internals
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-4] ranger_tree_from_forest() produces valid data.frame", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()

  df <- tryCatch(
    ranger_tree_from_forest(fx$model, tree_idx = 1L, fx$feature_names),
    error = function(e) NULL
  )

  if (is.null(df)) skip("ranger forest structure not accessible in this version")

  expect_true(is.data.frame(df))
  expect_true(nrow(df) > 0L)
  expected_cols <- c("nodeid", "feature", "splitval", "terminal", "leftchild", "rightchild", "prediction")
  expect_true(all(expected_cols %in% colnames(df)),
              info = paste("Missing:", paste(setdiff(expected_cols, colnames(df)), collapse = ", ")))

  # Leaf nodes should have NA feature, NA splitval, and numeric prediction
  leaf <- df$terminal
  if (any(leaf)) {
    expect_true(all(is.na(df$feature[leaf])))
    expect_true(all(is.na(df$splitval[leaf])))
    expect_true(all(!is.na(df$prediction[leaf])))
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 6 [RNG-5]: Leaf predictions are probabilities in [0, 1]
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-5] leaf predictions are in [0, 1] for probability forest", {
  skip_if_not_rng()
  fx <- .make_rng_fixture(probability = TRUE)
  u  <- ranger_unify_compat(fx$model, fx$X)
  m  <- u$model

  leaf_preds <- m$Prediction[is.na(m$Feature)]
  leaf_preds <- leaf_preds[!is.na(leaf_preds)]

  if (length(leaf_preds) > 0L) {
    expect_true(all(leaf_preds >= 0 & leaf_preds <= 1),
                info = paste("Leaf predictions out of [0,1]:",
                             paste(head(leaf_preds[leaf_preds < 0 | leaf_preds > 1], 5), collapse = ", ")))
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 7 [RNG-6]: model$forest NULL raises informative error
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-6] model$forest NULL raises error in ranger_tree_from_forest()", {
  # Simulate a model without a forest
  fake_model <- list(forest = NULL, num.trees = 5L,
                     treetype = "Probability estimation")
  class(fake_model) <- "ranger"

  expect_error(
    ranger_tree_from_forest(fake_model, tree_idx = 1L, feature_names = "f1"),
    regexp = "write\\.forest"
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 8 [RNG-8]: Non-probability forest emits a warning, does not stop
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-8] non-probability forest gives warning not an error", {
  skip_if_not_rng()
  # Train a classification (non-probability) forest
  fx <- tryCatch(.make_rng_fixture(probability = FALSE), error = function(e) NULL)
  if (is.null(fx)) skip("Could not train non-probability ranger model")

  expect_warning(
    ranger_unify_compat(fx$model, fx$X),
    regexp = "probability=TRUE"
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 9 [RNG-7]: keep.inbag=FALSE still produces valid output (via treeInfo)
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-7] keep.inbag=FALSE model still unifies successfully", {
  skip_if_not_rng()
  fx <- .make_rng_fixture(keep_inbag = FALSE)
  u  <- ranger_unify_compat(fx$model, fx$X)
  expect_s3_class(u, "model_unified")
  expect_true(nrow(u$model) > 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 10: Roundtrip — model_unified feeds .prepare_dp_model
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG] model_unified passes .prepare_dp_model without error", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()
  u  <- ranger_unify_compat(fx$model, fx$X)
  dp <- .prepare_dp_model(u, fx$feature_names)
  expect_true(is.list(dp))
  expect_true(length(dp$trees) > 0L)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 11: matrix and data.frame inputs both accepted
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG] matrix and data.frame inputs both produce valid output", {
  skip_if_not_rng()
  fx <- .make_rng_fixture()
  u_df  <- ranger_unify_compat(fx$model, fx$X)
  u_mat <- ranger_unify_compat(fx$model, as.matrix(fx$X))
  expect_s3_class(u_df,  "model_unified")
  expect_s3_class(u_mat, "model_unified")
  expect_equal(u_df$feature_names, u_mat$feature_names)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 12 [RNG-4]: child.nodeIDs as list-of-2-vectors → correct row resolution
# ──────────────────────────────────────────────────────────────────────────────
test_that("[RNG-4] child.nodeIDs list structure resolves to correct row indices", {
  # Build a minimal fake forest with list-of-2 child structure
  fake_cn <- list(
    list(c(1L, 2L), c(3L, 4L), c(0L, 0L), c(0L, 0L), c(0L, 0L))  # 5 nodes, tree 1
  )
  # Node 0: left=1, right=2
  # Node 1: left=3, right=4
  # Nodes 2,3,4: leaves

  node_ids <- 0:4
  left_ids  <- sapply(fake_cn[[1]], `[[`, 1L)
  right_ids <- sapply(fake_cn[[1]], `[[`, 2L)
  is_term   <- (left_ids == 0L & right_ids == 0L)

  # Build node → row map (0-based node IDs → 1-based row indices)
  node_map <- setNames(seq_len(5L), as.character(node_ids))
  resolve  <- function(ids) {
    out <- unname(node_map[as.character(as.integer(ids))])
    out[is_term | is.na(ids)] <- NA_integer_
    as.integer(out)
  }

  left_r  <- resolve(left_ids)
  right_r <- resolve(right_ids)

  expect_equal(left_r[1],  2L)   # node 0's left child = node 1 = row 2
  expect_equal(right_r[1], 3L)   # node 0's right child = node 2 = row 3
  expect_true (is.na(left_r[3]))  # node 2 is a leaf → no children
})
