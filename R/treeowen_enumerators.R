# treeowen/R/treeowen_enumerators.R
# Reference enumerator functions for testing and benchmarking.
# Source: treeowen_main.R v5.14, L1596–1784

#' Model-Agnostic Owen Value Enumerator (Reference Implementation)
#'
#' Computes exact Owen values by exhaustive combinatorial enumeration,
#' treating the ensemble as a black box. Intended as a correctness reference
#' for small problems; exponential in both \eqn{K} and \eqn{|G_k|}.
#' Use only for \eqn{K \leq 5} and \eqn{|G_k| \leq 5}.
#'
#' @param unified_model A unified model object (see \code{\link{treeowen}}).
#' @param x A data frame or matrix of observations.
#' @param groups A named list of character vectors partitioning all features.
#' @param chunk_size_inner Integer. Batch size for inner evaluation.
#'   Default \code{131072L}.
#' @param max_bytes Numeric. Memory limit for batched evaluation (bytes).
#'   Default \code{512 * 1024^2}.
#' @param dp_progress Logical. Show progress bar. Default \code{FALSE}.
#' @param use_cpp Logical. Use compiled C++ routines. Default \code{TRUE}.
#' @param verbose Logical. Print diagnostics. Default \code{FALSE}.
#' @return A \code{"treeowen_result"} object identical in structure to
#'   the output of \code{\link{treeowen}}.
#' @seealso \code{\link{treeowen}}, \code{\link{treeowen_exact_tuvalues}}
#' @export
treeowen_exact_enum <- function(unified_model, x, groups,
                                 chunk_size_inner = 131072L,
                                 max_bytes        = 512 * 1024^2,
                                 dp_progress      = FALSE,
                                 use_cpp          = TRUE,
                                 verbose          = FALSE) {
  stopifnot(inherits(unified_model, "model_unified"))
  if (isTRUE(use_cpp) && !.CPP_LOADED) .try_load_cpp()
  cpp_ok <- isTRUE(use_cpp) && .CPP_AVAILABLE

  x <- as.data.frame(x)
  if (!is.null(unified_model$feature_names))
    x <- x[, intersect(colnames(x), unified_model$feature_names), drop = FALSE]
  feat_names <- colnames(x); p <- ncol(x); n <- nrow(x)
  if (!p) stop("No usable columns in x.")
  groups  <- Filter(length, groups); g_feats <- unlist(groups, use.names = FALSE)
  if (length(setdiff(g_feats, feat_names)) > 0)  stop("groups contain unknown features.")
  if (length(g_feats) != length(unique(g_feats))) stop("groups must form a partition.")
  if (length(setdiff(feat_names, g_feats)) > 0)  stop("Each x column must belong to exactly one group.")
  K <- length(groups)
  if (K > TREEOWEN_EXACT_ENUM_K_MAX)
    stop(sprintf("Exact enum requires K <= %d (C++ limit: d_minus < 30). Got K=%d.",
                 TREEOWEN_EXACT_ENUM_K_MAX, K))
  if (verbose) message(sprintf("treeowen_exact_enum: n=%d, p=%d, K=%d | cpp=%s", n, p, K, cpp_ok))

  dp_model <- .prepare_dp_model(unified_model, feat_names)
  maxj     <- if (length(dp_model$used_feat_cols)) max(dp_model$used_feat_cols) else 0L
  fpos     <- setNames(seq_len(p), feat_names)
  group_pos <- lapply(groups, function(g) as.integer(fpos[g]))

  if (cpp_ok && exists("inner_shapley_enum_one_group_cpp", mode = "function")) {
    ow <- matrix(0, n, p); base <- numeric(n); vall <- numeric(n)
    others_list_all <- lapply(seq_len(K), function(g)
      lapply(setdiff(seq_len(K), g), function(og) as.integer(group_pos[[og]])))
    trees_xptr <- .prepare_trees_xptr_once(dp_model)
    for (irow in seq_len(n)) {
      xr <- as.numeric(x[irow, , drop = TRUE])
      base[irow] <- eval_v_knownset_cpp(trees_xptr, xr, integer(0), maxj)
      vall[irow] <- eval_v_knownset_cpp(trees_xptr, xr, seq_len(p), maxj)
      ow_row <- numeric(p)
      for (g in seq_len(K)) {
        phi_g <- inner_shapley_enum_one_group_cpp(
          trees_xptr         = trees_xptr,
          x_row              = xr,
          maxj               = as.integer(maxj),
          group_feats_g      = as.integer(group_pos[[g]]),
          group_feats_others = others_list_all[[g]],
          K                  = as.integer(K),
          p                  = as.integer(p)
        )
        ow_row <- ow_row + as.numeric(phi_g)
      }
      ow[irow, ] <- ow_row
    }
    colnames(ow) <- feat_names; rownames(ow) <- rownames(x)
    ts  <- .model_tree_stats(unified_model)
    res <- list(owens = as.data.frame(ow), baseline = base, v_all = vall,
                unified_model = unified_model, observations = x, groups = groups,
                note  = sprintf("treeowen_exact_enum | n=%d, p=%d, K=%d | cpp=%s", n, p, K, cpp_ok),
                stats = list(n = n, p = p, K = K,
                             trees = ts$ntrees, total_nodes = ts$total_nodes))
    class(res) <- c("treeowen_result", "list")
    return(res)
  }
  stop("[TreeOwen] treeowen_exact_enum: C++ required. Run .try_load_cpp() first.")
}


#' Tree-Aware Owen Value Enumerator (Reference Implementation)
#'
#' Computes exact Owen values by exhaustive enumeration using tree-aware
#' dynamic programming (Step 2 of \code{\link{treeowen}} only; no outer
#' hierarchy aggregation). Intended as a correctness reference; still
#' exponential in \eqn{K}.
#' Use only for \eqn{K \leq 5} and \eqn{|G_k| \leq 5}.
#'
#' @param unified_model A unified model object (see \code{\link{treeowen}}).
#' @param x A data frame or matrix of observations.
#' @param groups A named list of character vectors partitioning all features.
#' @param chunk_size_inner Integer. Batch size for inner evaluation.
#'   Default \code{131072L}.
#' @param max_bytes Numeric. Memory limit for batched evaluation (bytes).
#'   Default \code{512 * 1024^2}.
#' @param inner_bitmask_max Integer. Maximum group size for bitmask enumeration.
#'   Default \code{TREEOWEN_INNER_BITMASK_DEFAULT}.
#' @param dp_progress Logical. Show progress bar. Default \code{FALSE}.
#' @param use_cpp Logical. Use compiled C++ routines. Default \code{TRUE}.
#' @param verbose Logical. Print diagnostics. Default \code{FALSE}.
#' @return A \code{"treeowen_result"} object identical in structure to
#'   the output of \code{\link{treeowen}}.
#' @seealso \code{\link{treeowen}}, \code{\link{treeowen_exact_enum}}
#' @export
treeowen_exact_tuvalues <- function(unified_model, x, groups,
                                     chunk_size_inner  = 131072L,
                                     max_bytes         = 512 * 1024^2,
                                     inner_bitmask_max = TREEOWEN_INNER_BITMASK_DEFAULT,
                                     dp_progress       = FALSE,
                                     use_cpp           = TRUE,
                                     verbose           = FALSE) {
  if (isTRUE(use_cpp) && !.CPP_LOADED) .try_load_cpp()
  cpp_ok <- isTRUE(use_cpp) && .CPP_AVAILABLE

  x <- as.data.frame(x)
  if (!is.null(unified_model$feature_names))
    x <- x[, intersect(colnames(x), unified_model$feature_names), drop = FALSE]
  feat_names <- colnames(x); p <- ncol(x); n <- nrow(x)
  if (!p) stop("No usable columns in x.")
  groups <- Filter(length, groups); K <- length(groups)

  group_pos <- lapply(groups, function(g) {
    idx <- match(g, feat_names)
    if (anyNA(idx))
      stop(sprintf("treeowen_exact_tuvalues: group feature(s) not found in x: %s",
                   paste(g[is.na(idx)], collapse = ", ")))
    as.integer(idx)
  })

  dp_model <- .prepare_dp_model(unified_model, feat_names)
  maxj     <- if (length(dp_model$used_feat_cols)) max(dp_model$used_feat_cols) else 0L
  if (!cpp_ok)
    stop("[TreeOwen] treeowen_exact_tuvalues: C++ required. Run .try_load_cpp() first.")
  trees_xptr_cache <- .prepare_trees_xptr_once(dp_model)

  make_char_func <- function(xr) {
    function(coalition) {
      coalition <- as.integer(coalition)
      coalition <- coalition[!is.na(coalition) & coalition >= 1L & coalition <= p]
      eval_v_knownset_cpp(trees_xptr_cache, xr, coalition, as.integer(maxj))
    }
  }

  .one_row_ma <- function(xr) {
    char_func <- make_char_func(xr)
    phi <- numeric(p)
    for (k in seq_len(K)) {
      feats_k <- group_pos[[k]]; m_k <- length(feats_k)
      if (m_k == 0L) next
      others  <- setdiff(seq_len(K), k); n_other <- length(others)
      w_T_tab <- vapply(0L:(K - 1L),   function(t) .safe_shapley_w(t, K),   numeric(1L))
      w_S_tab <- vapply(0L:(m_k - 1L), function(s) .safe_shapley_w(s, m_k), numeric(1L))
      for (tmask in 0L:(bitwShiftL(1L, n_other) - 1L)) {
        t_bits <- as.logical(intToBits(tmask))[seq_len(n_other)]
        T_grp  <- others[which(t_bits)]
        wT     <- w_T_tab[length(T_grp) + 1L]
        T_feat <- if (length(T_grp))
          as.integer(unlist(group_pos[T_grp], use.names = FALSE)) else integer(0)
        for (ii in seq_len(m_k)) {
          feat_i  <- feats_k[ii]; Gminus <- feats_k[-ii]; n_inner <- length(Gminus)
          for (smask in 0L:(bitwShiftL(1L, n_inner) - 1L)) {
            s_bits <- as.logical(intToBits(smask))[seq_len(n_inner)]
            S_feat <- as.integer(Gminus[which(s_bits)])
            wS     <- w_S_tab[length(S_feat) + 1L]
            base_coal <- unique(c(T_feat, S_feat))
            v0 <- char_func(base_coal)
            v1 <- char_func(unique(c(base_coal, feat_i)))
            phi[feat_i] <- phi[feat_i] + wT * wS * (v1 - v0)
          }
        }
      }
    }
    phi
  }

  ow <- matrix(0, n, p); base <- numeric(n); vall <- numeric(n)
  for (irow in seq_len(n)) {
    xr            <- as.numeric(x[irow, , drop = TRUE])
    char_func_row <- make_char_func(xr)
    base[irow]    <- char_func_row(integer(0))
    vall[irow]    <- char_func_row(seq_len(p))
    ow[irow, ]    <- .one_row_ma(xr)
  }

  colnames(ow) <- feat_names; rownames(ow) <- rownames(x)
  ts  <- .model_tree_stats(unified_model)
  res <- list(owens = as.data.frame(ow), baseline = base, v_all = vall,
              unified_model = unified_model, observations = x, groups = groups,
              note  = sprintf("treeowen_exact_tuvalues | n=%d, p=%d, K=%d | cpp=%s", n, p, K, cpp_ok),
              stats = list(n = n, p = p, K = K,
                           trees = ts$ntrees, total_nodes = ts$total_nodes))
  class(res) <- c("treeowen_result", "list")
  res
}
