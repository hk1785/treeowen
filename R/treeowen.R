# treeowen/R/treeowen.R
# Main exported function: treeowen()

#' Compute Owen Values for Tree-Based Ensembles
#'
#' Computes exact or Monte Carlo Owen values for every observation in \code{x}
#' using a trained tree-based ensemble \code{unified_model} and a feature
#' partition \code{groups}.
#'
#' @param unified_model A unified model object produced by
#'   \code{xgboost_unify_compat()}, \code{lightgbm_unify_compat()}, or
#'   \code{ranger_unify_compat()}.  These wrappers are robust replacements
#'   for \code{treeshap::xgboost.unify()}, \code{treeshap::lightgbm.unify()},
#'   and \code{treeshap::ranger.unify()} respectively, handling version
#'   incompatibilities across xgboost >= 1.6, lightgbm >= 3.x / 4.x, and
#'   ranger >= 0.11 / 0.14+.
#' @param x A data frame or matrix of observations (rows = observations,
#'   columns = features). Must contain the same features as \code{unified_model}.
#' @param groups A named list of character vectors. Each element gives the
#'   feature names belonging to that group. Must partition all features in
#'   \code{unified_model}.
#' @param method Character. One of \code{"auto"} (default), \code{"exact"}, or
#'   \code{"approx"}. When \code{"auto"}, the exact algorithm is used for
#'   groups with \code{|G_k| < auto_exact_max_m} and Monte Carlo for larger
#'   groups.
#' @param hierarchy Controls the auxiliary binary tree \eqn{\mathcal{B}} over
#'   groups. Accepted values:
#'   \describe{
#'     \item{\code{NULL}}{(default) Balanced binary tree over the \eqn{K} groups.}
#'     \item{Named list}{Nested-list tree compatible with
#'       \code{build_hierarchy_tree_binary()}.}
#'     \item{Numeric vector or \eqn{K}-column integer matrix}{Layer-based merge
#'       specification; each row merges groups sharing the same integer label.}
#'   }
#'   The hierarchy does not affect theoretical correctness; it influences runtime.
#' @param n_inner_mc Integer. Initial number of Monte Carlo permutation pairs
#'   for the inner estimator. Default \code{64L}.
#' @param inner_antithetic Logical. Use antithetic sampling. Default \code{TRUE}.
#' @param target_se_inner Numeric. Target standard error for adaptive MC
#'   stopping. Default \code{1e-3}.
#' @param min_inner_mc Integer. Minimum number of MC pairs per group context.
#'   Default \code{32L}.
#' @param max_inner_mc Integer. Maximum number of MC pairs per group context.
#'   Default \code{512L}.
#' @param chunk_size_inner Integer. Chunk size for batched inner evaluation.
#'   Default \code{131072L}.
#' @param max_bytes Numeric. Memory limit (bytes) for inner evaluation batches.
#'   Default \code{512 * 1024^2} (512 MB).
#' @param inner_bitmask_max Integer. Maximum group size for bitmask-based exact
#'   enumeration. Default \code{TREEOWEN_INNER_BITMASK_DEFAULT} (60).
#' @param auto_exact_max_m Integer. Groups with \code{|G_k| >= auto_exact_max_m}
#'   use Monte Carlo when \code{method = "auto"}. Default
#'   \code{TREEOWEN_AUTO_EXACT_MAX_M} (30).
#' @param check_efficiency Logical. Check the Owen value efficiency axiom after
#'   computation. Default \code{FALSE}.
#' @param efficiency_tol Numeric. Tolerance for the efficiency check.
#'   Default \code{1e-6}.
#' @param dp_progress Logical. Display a progress bar. Default \code{TRUE}.
#' @param dp_print_every Integer. Progress bar update frequency.
#'   Default \code{1L}.
#' @param dp_newline_every Integer. Deprecated; retained for backwards
#'   compatibility only. Default \code{10L}.
#' @param eta_alpha Numeric. EMA smoothing for ETA estimation. Default
#'   \code{0.2}.
#' @param n_cores Integer. Number of parallel cores. Default \code{1L}.
#'   Values \code{> 1} use \code{parallel::mclapply()} on Linux/macOS;
#'   Windows requires \code{n_cores = 1L}.
#' @param use_cpp Logical. Use compiled C++ routines. Default \code{TRUE}.
#' @param verbose Logical. Print diagnostic messages. Default \code{TRUE}.
#'
#' @return An object of class \code{"treeowen_result"} (a list) with components:
#' \describe{
#'   \item{\code{owens}}{Data frame of Owen values (n × p).}
#'   \item{\code{observations}}{Data frame of input observations (n × p).}
#'   \item{\code{groups}}{The \code{groups} argument.}
#'   \item{\code{baseline}}{Numeric vector \eqn{v_x(\emptyset) = E[f(X)]}.}
#'   \item{\code{v_all}}{Numeric vector \eqn{v_x(N) = f(x)} for each observation.}
#'   \item{\code{stats}}{A list of diagnostic statistics.}
#'   \item{\code{note}}{Character summary of the run configuration.}
#' }
#'
#' @references
#' Owen, G. (1977). Values of games with a priori unions.
#' In: Mathematical Economics and Game Theory. Springer, Berlin.
#'
#' Lundberg, S. M. et al. (2020). From local explanations to global
#' understanding with explainable AI for trees. \emph{Nature Machine
#' Intelligence}, 2(1), 56–67.
#'
#' @seealso \code{\link{treeowen_importance}}, \code{\link{treeowen_beeswarm}},
#'   \code{\link{treeowen_hierarchical_beeswarm}}
#'
#' @examples
#' \dontrun{
#' unified <- xgboost_unify_compat(model, X)
#' # or: unified <- lightgbm_unify_compat(lgb_model, X)
#' # or: unified <- ranger_unify_compat(ranger_model, X)
#' groups  <- list(G1 = c("f1", "f2"), G2 = c("f3", "f4"))
#' result  <- treeowen(unified, X, groups)
#' print(result)
#' }
#'
#' @export
treeowen <- function(unified_model, x, groups,
                     method             = c("auto", "exact", "approx"),
                     hierarchy          = NULL,
                     n_inner_mc         = 64L,
                     inner_antithetic   = TRUE,
                     target_se_inner    = 1e-3,
                     min_inner_mc       = 32L,
                     max_inner_mc       = 512L,
                     chunk_size_inner   = 131072L,
                     max_bytes          = 512 * 1024^2,
                     inner_bitmask_max  = TREEOWEN_INNER_BITMASK_DEFAULT,
                     auto_exact_max_m   = TREEOWEN_AUTO_EXACT_MAX_M,
                     check_efficiency   = FALSE,
                     efficiency_tol     = 1e-6,
                     dp_progress        = TRUE,
                     dp_print_every     = 1L,
                     dp_newline_every   = 10L,
                     eta_alpha          = 0.2,
                     n_cores            = 1L,
                     use_cpp            = TRUE,
                     verbose            = TRUE) {

  method <- match.arg(method)
  if (isTRUE(use_cpp) && !.CPP_LOADED) .try_load_cpp()
  cpp_ok <- isTRUE(use_cpp) && .CPP_AVAILABLE

  cm <- .prepare_common(unified_model, x, groups, max_bytes, chunk_size_inner)
  p  <- cm$p; n <- cm$n; K <- cm$K
  group_pos  <- cm$group_pos; maxj <- cm$maxj; dp_model <- cm$dp_model

  hier_tree <- if (is.null(hierarchy)) {
    build_hierarchy_tree_from_layers(K, NULL)
  } else if (is.list(hierarchy) || is.matrix(hierarchy) || is.numeric(hierarchy)) {
    build_hierarchy_tree_from_layers(K, hierarchy)
  } else {
    stop("hierarchy must be NULL, a list, a numeric vector, or a matrix.")
  }
  ann     <- .annotate_tree(hier_tree, K)
  max_deg <- .max_degree_hier_tree(hier_tree)

  group_sizes <- vapply(group_pos, length, integer(1L))

  group_methods <- if (method == "auto") {
    ifelse(group_sizes >= as.integer(auto_exact_max_m), "approx", "exact")
  } else {
    rep(method, K)
  }

  effective_method <- if (all(group_methods == "exact"))       "exact"
                      else if (all(group_methods == "approx")) "approx"
                      else                                     "mixed"

  n_exact_groups  <- sum(group_methods == "exact")
  n_approx_groups <- sum(group_methods == "approx")

  if (isTRUE(verbose)) {
    message(sprintf(
      "[TreeOwen] n=%d p=%d K=%d | method=%s (effective=%s) | cpp=%s | n_cores=%d",
      n, p, K, method, effective_method, cpp_ok, as.integer(n_cores)))
    if (effective_method == "mixed")
      message(sprintf(
        "[TreeOwen] per-group methods: %d exact (|G_k| < %d), %d approx (|G_k| >= %d)",
        n_exact_groups, as.integer(auto_exact_max_m),
        n_approx_groups, as.integer(auto_exact_max_m)))
  }

  n_inner_mc   <- as.integer(n_inner_mc);   if (!is.finite(n_inner_mc)   || n_inner_mc   <= 0L) n_inner_mc   <- 64L
  min_inner_mc <- as.integer(min_inner_mc); if (!is.finite(min_inner_mc) || min_inner_mc <= 0L) min_inner_mc <- 32L
  max_inner_mc <- as.integer(max_inner_mc); if (!is.finite(max_inner_mc) || max_inner_mc < min_inner_mc) max_inner_mc <- max(min_inner_mc, 512L)

  trees_xptr_cache <- if (cpp_ok) .prepare_trees_xptr_once(dp_model) else NULL

  leaf_contexts_all <- .precompute_leaf_contexts(ann, K, group_pos, maxj)

  leaf_contexts <- Filter(function(ctx) abs(ctx$weight) > .Machine$double.eps,
                          leaf_contexts_all)
  if (isTRUE(verbose) && length(leaf_contexts) < length(leaf_contexts_all))
    message(sprintf("[TreeOwen] leaf_contexts: %d → %d (zero-weight removed)",
                    length(leaf_contexts_all), length(leaf_contexts)))

  x_mat        <- as.matrix(cm$x)
  all_cols_int <- as.integer(unlist(group_pos, use.names = FALSE))

  .make_exact_solver <- function() {
    function(ctx, x_row_dbl) {
      .inner_shapley_in_group_given_context_direct(
        dp_model = dp_model, x_row = x_row_dbl, g = ctx$g,
        group_pos = group_pos, outer_known_cols = ctx$outer_known_cols,
        maxj = maxj, p = p, chunk_size = cm$chunk_size_inner,
        max_bytes = cm$max_bytes, inner_bitmask_max = inner_bitmask_max,
        use_cpp = cpp_ok, trees_xptr = trees_xptr_cache)
    }
  }

  .make_approx_solver <- function() {
    function(ctx, x_row_dbl) {
      .inner_perm_mc_one_group_direct(
        dp_model = dp_model, x_row = x_row_dbl, g = ctx$g,
        group_pos = group_pos, outer_known_cols = ctx$outer_known_cols,
        maxj = maxj, p = p, n_inner_mc = n_inner_mc,
        inner_antithetic = inner_antithetic, target_se = target_se_inner,
        min_mc = min_inner_mc, max_mc = max_inner_mc,
        use_cpp = cpp_ok, trees_xptr = trees_xptr_cache)
    }
  }

  .group_solvers <- vector("list", K)
  for (.gi in seq_len(K)) {
    .group_solvers[[.gi]] <- if (group_methods[.gi] == "exact") .make_exact_solver()
                             else                               .make_approx_solver()
  }

  is_exact <- (effective_method == "exact")

  compute_one_row <- function(irow, use_cpp_worker = cpp_ok) {
    x_row_dbl <- as.numeric(x_mat[irow, ])
    if (!use_cpp_worker || is.null(trees_xptr_cache))
      stop("[TreeOwen] compute_one_row: C++ XPtr missing. Ensure C++ is loaded.")
    base <- eval_v_knownset_cpp(trees_xptr_cache, x_row_dbl, integer(0),   as.integer(maxj))
    vall <- eval_v_knownset_cpp(trees_xptr_cache, x_row_dbl, all_cols_int, as.integer(maxj))
    phi <- numeric(p)
    for (ctx in leaf_contexts) {
      phi <- phi + ctx$weight * .group_solvers[[ctx$g]](ctx, x_row_dbl)
    }
    if (isTRUE(check_efficiency)) {
      lhs <- sum(phi); rhs <- vall - base
      if (is.finite(lhs) && is.finite(rhs) && abs(lhs - rhs) > as.numeric(efficiency_tol))
        warning(sprintf("Efficiency check failed (row %d)%s: sum(phi)=%.12g vs %.12g",
                        irow,
                        if (!is_exact) " [MC approx or mixed; deviation expected]" else "",
                        lhs, rhs))
    }
    list(phi = phi, base = base, vall = vall)
  }

  # Chunked worker: processes a block of row indices. Used by
  # .run_row_loop_parallel() for the doParallel/foreach (Windows) path and
  # as the inner kernel for mclapply-based dispatch on Linux/macOS.
  make_chunk_worker <- function(allow_cpp) {
    local({
      xptr_w <- NULL
      function(irows) {
        if (is.null(xptr_w)) {
          if (!allow_cpp)
            stop("[TreeOwen] Windows foreach worker: XPtr serialization not possible. Use n_cores=1 or Linux/macOS.")
          xptr_w <<- .prepare_trees_xptr_once(dp_model)
        }
        call_exact_solver <- function(ctx, x_row_dbl)
          .inner_shapley_in_group_given_context_direct(
            dp_model = dp_model, x_row = x_row_dbl, g = ctx$g,
            group_pos = group_pos, outer_known_cols = ctx$outer_known_cols,
            maxj = maxj, p = p, chunk_size = cm$chunk_size_inner,
            max_bytes = cm$max_bytes, inner_bitmask_max = inner_bitmask_max,
            use_cpp = TRUE, trees_xptr = xptr_w)
        call_approx_solver <- function(ctx, x_row_dbl)
          .inner_perm_mc_one_group_direct(
            dp_model = dp_model, x_row = x_row_dbl, g = ctx$g,
            group_pos = group_pos, outer_known_cols = ctx$outer_known_cols,
            maxj = maxj, p = p, n_inner_mc = n_inner_mc,
            inner_antithetic = inner_antithetic, target_se = target_se_inner,
            min_mc = min_inner_mc, max_mc = max_inner_mc,
            use_cpp = TRUE, trees_xptr = xptr_w)

        out <- vector("list", length(irows))
        for (li in seq_along(irows)) {
          irow      <- irows[[li]]
          x_row_dbl <- as.numeric(x_mat[irow, ])
          base <- eval_v_knownset_cpp(xptr_w, x_row_dbl, integer(0),   as.integer(maxj))
          vall <- eval_v_knownset_cpp(xptr_w, x_row_dbl, all_cols_int, as.integer(maxj))
          phi  <- numeric(p)
          for (ctx in leaf_contexts) {
            solver_fn <- if (group_methods[ctx$g] == "exact") call_exact_solver else call_approx_solver
            phi <- phi + ctx$weight * solver_fn(ctx, x_row_dbl)
          }
          if (isTRUE(check_efficiency)) {
            lhs <- sum(phi); rhs <- vall - base
            if (is.finite(lhs) && is.finite(rhs) && abs(lhs - rhs) > as.numeric(efficiency_tol))
              warning(sprintf("Efficiency check failed (row %d)%s: sum(phi)=%.12g vs %.12g",
                              irow,
                              if (!is_exact) " [MC approx or mixed; deviation expected]" else "",
                              lhs, rhs))
          }
          out[[li]] <- list(phi = phi, base = base, vall = vall)
        }
        out
      }
    })
  }

  # Per-row worker (wraps a single index). Used by pbmcapply-based dispatch
  # on Linux/macOS where each row is its own mclapply task so the progress
  # bar increments for every row.
  make_row_worker <- function(allow_cpp) {
    chunk_w <- make_chunk_worker(allow_cpp)
    function(irow) chunk_w(irow)[[1L]]
  }

  compute_chunk_fork      <- make_chunk_worker(allow_cpp = cpp_ok)
  compute_row_fork        <- make_row_worker(allow_cpp = cpp_ok)
  compute_chunk_win       <- make_chunk_worker(allow_cpp = FALSE)

  phase_tag <- if (effective_method == "exact") "DP" else if (effective_method == "approx") "MC" else "DP/MC"
  lr <- .run_row_loop(n,
                      compute_one_row_serial  = compute_one_row,
                      compute_row_fork        = compute_row_fork,
                      compute_chunk_fork      = compute_chunk_fork,
                      compute_chunk_win       = compute_chunk_win,
                      n_cores = n_cores, verbose = verbose,
                      dp_progress = dp_progress, dp_print_every = dp_print_every,
                      dp_newline_every = dp_newline_every, eta_alpha = eta_alpha,
                      K = K, chunk_size_inner = cm$chunk_size_inner,
                      chunk_cap = cm$chunk_cap, phase_tag = phase_tag,
                      cpp_available = cpp_ok)

  .make_treeowen_result(
    lr$ow, lr$basev, lr$vallv, cm$feat_names,
    unified_model, cm$x, cm$groups,
    stats = list(p = p, K = K, n = n, trees = cm$ts$ntrees,
                 total_nodes = cm$ts$total_nodes,
                 chunk_size_inner = cm$chunk_size_inner, max_bytes = cm$max_bytes,
                 max_nn = cm$max_nn, method_requested = method,
                 method_effective = effective_method,
                 group_methods = group_methods,
                 n_exact_groups = n_exact_groups,
                 n_approx_groups = n_approx_groups,
                 inner_bitmask_max = as.integer(inner_bitmask_max),
                 n_cores = as.integer(n_cores),
                 cpp_used = cpp_ok, max_deg = max_deg, max_m = cm$max_m,
                 auto_exact_max_m = as.integer(auto_exact_max_m),
                 n_inner_mc = n_inner_mc,
                 n_leaf_contexts = length(leaf_contexts),
                 n_leaf_contexts_raw = length(leaf_contexts_all)),
    note = sprintf("treeowen | method=%s (effective=%s, exact=%d/approx=%d groups) | n_cores=%d | cpp=%s",
                   method, effective_method,
                   n_exact_groups, n_approx_groups,
                   as.integer(n_cores), cpp_ok)
  )
}


#' @export
print.treeowen_result <- function(x, ...) {
  s <- x$stats
  cat(sprintf("TreeOwen result | n=%d, p=%d, K=%d\n", s$n, s$p, s$K))
  cat(sprintf("  method   : %s\n", x$note))
  cat(sprintf("  baseline : mean=%.6g  range=[%.6g, %.6g]\n",
              mean(x$baseline), min(x$baseline), max(x$baseline)))
  cat(sprintf("  v_all    : mean=%.6g  range=[%.6g, %.6g]\n",
              mean(x$v_all), min(x$v_all), max(x$v_all)))
  cat(sprintf("  trees    : %d (%d total nodes)\n", s$trees, s$total_nodes))
  invisible(x)
}
