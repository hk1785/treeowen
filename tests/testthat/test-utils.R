# tests/testthat/test-utils.R
#
# Unit tests for internal utility functions in utils.R
# Tests: .prepare_dp_model, .validate_inputs, group/bitset helpers,
#        .safe_shapley_w, .popcount helpers, and model_unified attribute checks.

# ──────────────────────────────────────────────────────────────────────────────
# Helper: build a minimal model_unified object by hand
# ──────────────────────────────────────────────────────────────────────────────
.make_minimal_unified <- function() {
  # Single-tree: root splits on feature "a"; left=leaf(0.8), right=leaf(0.2)
  model_df <- data.frame(
    Tree          = c(0L, 0L, 0L),
    Node          = c(0L, 1L, 2L),
    Feature       = c("a", NA, NA),
    Decision.type = factor(c("<=", NA, NA), levels = c("<=", "<")),
    Split         = c(0.0, NA, NA),
    Yes           = c(2L, NA, NA),
    No            = c(3L, NA, NA),
    Missing       = c(NA_integer_, NA, NA),
    Prediction    = c(NA, 0.8, 0.2),
    Cover         = c(10, 5, 5),
    stringsAsFactors = FALSE
  )
  ret <- list(
    model         = model_df,
    data          = data.frame(a = rnorm(5), b = rnorm(5)),
    feature_names = c("a", "b")
  )
  class(ret) <- "model_unified"
  attr(ret, "missing_support") <- FALSE
  attr(ret, "model")           <- "xgboost"
  ret
}

# ──────────────────────────────────────────────────────────────────────────────
# TEST 1: .prepare_dp_model() basic structure
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .prepare_dp_model() returns correct structure", {
  u  <- .make_minimal_unified()
  dp <- .prepare_dp_model(u, c("a", "b"))

  expect_true(is.list(dp))
  expect_true("trees"      %in% names(dp))
  expect_true("feat_names" %in% names(dp))
  expect_true("has_missing" %in% names(dp))
  expect_equal(length(dp$trees), 1L)
  expect_equal(dp$feat_names, c("a", "b"))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 2: .prepare_dp_model() tree node count
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .prepare_dp_model() tree has correct number of nodes", {
  u  <- .make_minimal_unified()
  dp <- .prepare_dp_model(u, c("a", "b"))
  t1 <- dp$trees[[1L]]
  expect_equal(length(t1$nodes_global), 3L)  # root + 2 leaves
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 3: .validate_inputs() passes valid inputs
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .validate_inputs() passes for valid inputs", {
  u  <- .make_minimal_unified()
  X  <- data.frame(a = 1:5, b = 1:5)
  g  <- list(G1 = "a", G2 = "b")
  expect_true(isTRUE(.validate_inputs(u, X, g)))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 4: .validate_inputs() rejects incomplete partition
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .validate_inputs() rejects incomplete feature partition", {
  u <- .make_minimal_unified()
  X <- data.frame(a = 1:5, b = 1:5)
  g <- list(G1 = "a")   # "b" not assigned
  expect_error(.validate_inputs(u, X, g), regexp = "FULL partition|unassigned|uncovered")
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 5: .validate_inputs() rejects duplicate feature assignment
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .validate_inputs() rejects features in multiple groups", {
  u <- .make_minimal_unified()
  X <- data.frame(a = 1:5, b = 1:5)
  g <- list(G1 = c("a", "b"), G2 = c("a"))  # "a" in two groups
  expect_error(.validate_inputs(u, X, g), regexp = "multiple groups|duplicat")
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 6: .safe_shapley_w() — known values
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .safe_shapley_w() gives correct values for small m", {
  # For m=2: w(s=0,m=2) = 0!*1!/2! = 0.5
  #          w(s=1,m=2) = 1!*0!/2! = 0.5
  expect_equal(.safe_shapley_w(0L, 2L), 0.5, tolerance = 1e-12)
  expect_equal(.safe_shapley_w(1L, 2L), 0.5, tolerance = 1e-12)

  # For m=3: w(0,3)=2/6=1/3, w(1,3)=1/6, w(2,3)=2/6=1/3
  expect_equal(.safe_shapley_w(0L, 3L), 2/6, tolerance = 1e-12)
  expect_equal(.safe_shapley_w(1L, 3L), 1/6, tolerance = 1e-12)
  expect_equal(.safe_shapley_w(2L, 3L), 2/6, tolerance = 1e-12)

  # Weights sum to 1 for any m
  for (m in 1:5) {
    ws <- vapply(0:(m-1), function(s) .safe_shapley_w(s, m), numeric(1))
    expect_equal(sum(ws), 1.0, tolerance = 1e-10,
                 info = paste("m =", m))
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 7: .popcount_int() — correct popcount table
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .popcount_int() returns correct bit counts", {
  pc <- .popcount_int(7L)   # values 0..7
  expect_equal(pc[1L], 0L)   # 0 = 000
  expect_equal(pc[2L], 1L)   # 1 = 001
  expect_equal(pc[3L], 1L)   # 2 = 010
  expect_equal(pc[4L], 2L)   # 3 = 011
  expect_equal(pc[5L], 1L)   # 4 = 100
  expect_equal(pc[6L], 2L)   # 5 = 101
  expect_equal(pc[7L], 2L)   # 6 = 110
  expect_equal(pc[8L], 3L)   # 7 = 111
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 8: bitset helpers — round-trip from IDs to mask and back
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] bitset round-trip: .bitset_from_ids / .bitset_has", {
  K    <- 10L
  ids  <- c(1L, 3L, 7L, 10L)
  mask <- .bitset_from_ids(ids, K)

  for (g in seq_len(K)) {
    expected <- g %in% ids
    expect_equal(.bitset_has(mask, g), expected,
                 info = paste("g =", g))
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 9: .bitset_count() correctness
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .bitset_count() returns correct number of set bits", {
  K    <- 8L
  ids  <- c(2L, 4L, 8L)
  mask <- .bitset_from_ids(ids, K)
  expect_equal(.bitset_count(mask, K), 3L)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 10: .bitset_or() merges two bitsets correctly
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] .bitset_or() correctly merges two bitsets", {
  K  <- 6L
  m1 <- .bitset_from_ids(c(1L, 3L), K)
  m2 <- .bitset_from_ids(c(3L, 5L), K)
  mo <- .bitset_or(m1, m2)
  expect_true (.bitset_has(mo, 1L))
  expect_false(.bitset_has(mo, 2L))
  expect_true (.bitset_has(mo, 3L))
  expect_false(.bitset_has(mo, 4L))
  expect_true (.bitset_has(mo, 5L))
  expect_false(.bitset_has(mo, 6L))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 11: build_hierarchy_tree_from_layers() — balanced tree for K=4
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] build_hierarchy_tree_from_layers(K=4) returns balanced tree", {
  K    <- 4L
  tree <- build_hierarchy_tree_from_layers(K, hierarchy = NULL)
  expect_true(is.list(tree))
  expect_equal(tree$type, "node")
  # All 4 leaf IDs present
  leaves <- sort(as.integer(tree$leaves))
  expect_equal(leaves, 1:4)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 12: build_hierarchy_tree_binary() — custom nested list
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] build_hierarchy_tree_binary() handles nested list correctly", {
  tree <- build_hierarchy_tree_binary(
    list(left = list(left = 1L, right = 2L),
         right = 3L)
  )
  expect_equal(tree$type, "node")
  leaves <- sort(as.integer(tree$leaves))
  expect_equal(leaves, 1:3)
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 13: model_unified attr checks
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] model_unified attributes are set correctly by compat wrappers", {
  skip_if_not_installed <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) testthat::skip(sprintf("Package '%s' not available", pkg))
  }
  skip_if_not_installed("xgboost")
  skip_if_not_installed("data.table")

  set.seed(1)
  X  <- matrix(rnorm(40 * 4), 40, 4, dimnames = list(NULL, paste0("f", 1:4)))
  Y  <- as.integer(rowSums(X) > 0)
  dm <- xgboost::xgb.DMatrix(X, label = Y)
  m  <- xgboost::xgboost(dm, nrounds = 5L, objective = "binary:logistic", verbose = 0)
  u  <- .xgboost_unify_compat(m, X)

  expect_equal(attr(u, "model"), "xgboost")
  expect_true(isTRUE(attr(u, "missing_support")))
})

# ──────────────────────────────────────────────────────────────────────────────
# TEST 14: clear_inner_enum_cache() empties the cache
# ──────────────────────────────────────────────────────────────────────────────
test_that("[utils] clear_inner_enum_cache() empties .INNER_ENUM_CACHE", {
  # Populate cache
  .get_inner_cache(2L)
  .get_inner_cache(3L)

  cache_env <- getNamespace("treeowen")$.INNER_ENUM_CACHE
  expect_true(length(ls(cache_env)) >= 1L)

  clear_inner_enum_cache()
  expect_equal(length(ls(cache_env)), 0L)
})
