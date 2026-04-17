# treeowen/R/treeowen_importance.R

#' Compute Feature and Group Importance from Owen Values
#'
#' Derives feature-level and group-level importance scores from a
#' \code{\link{treeowen}} result object.
#'
#' @param result An object of class \code{"treeowen_result"} returned by
#'   \code{\link{treeowen}}.
#' @param type Character. Which importance scores to compute:
#'   \code{"both"} (default), \code{"feature"}, or \code{"group"}.
#' @param group_agg Character. How to aggregate within-group feature importances:
#'   \describe{
#'     \item{\code{"sum_abs"}}{(default) \eqn{I_k = (1/n)\sum_j |\Psi_k^{(j)}|}.}
#'     \item{\code{"sum_abs_feat"}}{\eqn{I_k = (1/|G_k|)\sum_{i \in G_k} I_i}.}
#'     \item{\code{"both"}}{ Compute both; \code{importance} column uses
#'       \code{"sum_abs"} for sorting.}
#'   }
#' @param sort Logical. Sort by decreasing importance. Default \code{TRUE}.
#' @param normalize Logical. Normalize importance scores to sum to 1.
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{"treeowen_importance"} (a list) with
#'   components:
#' \describe{
#'   \item{\code{feature}}{Data frame with columns \code{feature} and
#'     \code{importance} (when \code{type \%in\% c("both","feature")}).}
#'   \item{\code{group}}{Data frame with columns \code{group} and
#'     \code{importance} (when \code{type \%in\% c("both","group")}).}
#'   \item{\code{group_attr}}{n × K matrix of group-level Owen attributions.}
#'   \item{\code{group_pos}}{Named list mapping group names to column indices.}
#'   \item{\code{group_names}, \code{feat_names}, \code{n}, \code{K},
#'     \code{p}, \code{type}, \code{group_agg}}{Metadata.}
#' }
#'
#' @seealso \code{\link{treeowen}}, \code{\link{treeowen_beeswarm}},
#'   \code{\link{treeowen_hierarchical_beeswarm}}
#'
#' @examples
#' \dontrun{
#' result <- treeowen(unified, X, groups)
#' imp    <- treeowen_importance(result)
#' print(imp)
#' head(imp$feature, 10)
#' head(imp$group)
#' }
#'
#' @export
treeowen_importance <- function(
    result,
    type        = c("both", "feature", "group"),
    group_agg   = c("sum_abs", "sum_abs_feat", "both"),
    sort        = TRUE,
    normalize   = FALSE
) {
  if (!inherits(result, c("treeowen_result", "list")))
    stop("result must be the output of treeowen() (class 'treeowen_result' or list).")
  if (is.null(result$owens) || is.null(result$groups))
    stop("result must contain $owens and $groups.")

  type      <- match.arg(type)
  group_agg <- match.arg(group_agg)

  ow_mat  <- as.matrix(result$owens)
  groups  <- result$groups
  n       <- nrow(ow_mat)
  p       <- ncol(ow_mat)
  K       <- length(groups)
  feat_names  <- colnames(ow_mat)
  group_names <- if (!is.null(names(groups))) names(groups) else paste0("G", seq_len(K))

  out <- list()

  if (type %in% c("both", "feature")) {
    feat_imp <- colMeans(abs(ow_mat))
    if (isTRUE(normalize) && sum(feat_imp) > 0) feat_imp <- feat_imp / sum(feat_imp)
    df_feat <- data.frame(feature = feat_names, importance = as.numeric(feat_imp),
                          stringsAsFactors = FALSE)
    if (isTRUE(sort)) df_feat <- df_feat[order(df_feat$importance, decreasing = TRUE), ]
    rownames(df_feat) <- NULL
    out$feature <- df_feat
  }

  fpos <- setNames(seq_len(p), feat_names)
  group_pos <- lapply(groups, function(g) {
    idx <- fpos[as.character(g)]; as.integer(idx[!is.na(idx)])
  })

  group_attr <- matrix(0, n, K, dimnames = list(rownames(ow_mat), group_names))
  for (k in seq_len(K)) {
    gp <- group_pos[[k]]
    if (length(gp)) group_attr[, k] <- rowSums(ow_mat[, gp, drop = FALSE])
  }

  if (type %in% c("both", "group")) {
    grp_imp_net  <- colMeans(abs(group_attr))
    grp_imp_feat <- vapply(seq_len(K), function(k) {
      gp <- group_pos[[k]]
      if (!length(gp)) return(0)
      mean(colMeans(abs(ow_mat[, gp, drop = FALSE])))
    }, numeric(1))

    if (group_agg == "sum_abs") {
      grp_imp <- grp_imp_net
    } else if (group_agg == "sum_abs_feat") {
      grp_imp <- grp_imp_feat
    } else {
      grp_imp <- grp_imp_net
    }

    if (isTRUE(normalize) && sum(grp_imp) > 0) grp_imp <- grp_imp / sum(grp_imp)

    df_grp <- data.frame(group = group_names, importance = as.numeric(grp_imp),
                          stringsAsFactors = FALSE)
    if (group_agg == "both") {
      norm_feat <- if (isTRUE(normalize) && sum(grp_imp_feat) > 0)
                     grp_imp_feat / sum(grp_imp_feat)
                   else grp_imp_feat
      df_grp$importance_sum_abs_feat <- as.numeric(norm_feat)
    }
    if (isTRUE(sort)) df_grp <- df_grp[order(df_grp$importance, decreasing = TRUE), ]
    rownames(df_grp) <- NULL
    out$group      <- df_grp
    out$group_attr <- group_attr
  } else {
    out$group_attr <- group_attr
  }

  out$group_pos   <- group_pos
  out$group_names <- group_names
  out$feat_names  <- feat_names
  out$n           <- n
  out$K           <- K
  out$p           <- p
  out$type        <- type
  out$group_agg   <- group_agg
  class(out) <- c("treeowen_importance", "list")
  out
}


#' @export
print.treeowen_importance <- function(x, digits = 4, ...) {
  cat(sprintf("TreeOwen Importance | n=%d, p=%d, K=%d\n", x$n, x$p, x$K))
  if (!is.null(x$feature)) {
    cat("\n[Feature Importance] (top 10)\n")
    print(head(x$feature, 10), digits = digits, row.names = FALSE)
  }
  if (!is.null(x$group)) {
    cat("\n[Group Importance]\n")
    print(x$group, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
