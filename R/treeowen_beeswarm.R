# treeowen/R/treeowen_beeswarm.R
# Exported function: treeowen_beeswarm()
# Renamed from treeowen_plot() in treeowen_main.R v5.14, L1972–2441.
# All internal references to "treeowen_plot" updated to "treeowen_beeswarm".

#' Beeswarm Plots of Owen Values (Feature and Group Level)
#'
#' Produces beeswarm plots of Owen value distributions at the feature level,
#' the group level, or both (side-by-side). Points are
#' coloured by the feature's raw value (feature plot) or by a group-level
#' summary statistic (group plot).
#'
#' @param result An object of class \code{"treeowen_result"} returned by
#'   \code{\link{treeowen}}.
#' @param level Character. Plotting level: \code{"feature"} (default),
#'   \code{"group"}, or \code{"both"}.
#' @param group_agg Character. Aggregation method for group importance used to
#'   determine sort order: \code{"sum_abs"} (default) or
#'   \code{"sum_abs_feat"}. \emph{Note:} \code{"both"} is not allowed here.
#' @param normalize_imp Logical. Normalize importance before sorting.
#'   Default \code{FALSE}.
#' @param top_n_feature Integer or \code{NULL}. Top features to display.
#'   Default \code{10L}.
#' @param top_n_group Integer or \code{NULL}. Top groups to display.
#'   Default \code{10L}.
#' @param group_color_stat Character. Statistic for group-level point colour:
#'   \code{"sum"} (default), \code{"mean"}, or \code{"custom"}.
#' @param group_color_fn Function or \code{NULL}. Required when
#'   \code{group_color_stat = "custom"}; signature
#'   \code{function(x_mat, feat_idx) -> numeric(n)}.
#' @param color_low,color_high Character. Hex colours for low/high values.
#'   Defaults: \code{"#0052A5"} (blue) and \code{"#DC2626"} (red).
#' @param color_alpha Numeric in \eqn{[0,1]}. Colour gradient transparency.
#'   Default \code{0.6}.
#' @param qlims Numeric vector of length 2. Quantile limits for colour
#'   clipping. Default \code{c(0.05, 0.95)}.
#' @param point_size Numeric. Point size. Default \code{1.6}.
#' @param point_alpha Numeric in \eqn{[0,1]}. Point transparency.
#'   Default \code{0.55}.
#' @param quasirandom_width Numeric. Spread of the quasirandom jitter.
#'   Default \code{0.4}.
#' @param title_feature,title_group Character or \code{NULL}. Plot titles.
#' @param xlab_feature,xlab_group Character. x-axis labels.
#' @param legend_label_feat Character. Colour-bar label for the feature panel.
#'   Default \code{"Value"}.
#' @param legend_label_grp Character or \code{NULL}. Colour-bar label for the
#'   group panel. \code{NULL} generates an automatic label.
#' @param plot_title_size,axis_title_x_size,axis_text_x_size,axis_text_y_size,legend_text_size,legend_title_size
#'   Numeric. Text sizes (pt).
#' @param margin_t,margin_r,margin_b,margin_l Numeric. Plot margins (pt).
#' @param legend_position,legend_direction Legend position and direction.
#' @param legend_barwidth,legend_barheight \code{grid::unit} objects for
#'   colour-bar dimensions.
#' @param verbose Logical. Print diagnostic messages. Default \code{FALSE}.
#'
#' @return
#' When \code{level = "feature"} or \code{level = "group"}: a \code{ggplot}
#' object.
#'
#' When \code{level = "both"}: an object of class \code{"treeowen_plots"}
#' (a list) with components:
#' \describe{
#'   \item{\code{combined}}{Patchwork of group | feature panels.}
#'   \item{\code{feature}}{Standalone feature \code{ggplot}.}
#'   \item{\code{group}}{Standalone group \code{ggplot}.}
#'   \item{\code{feat_ordered}}{Displayed features in importance order.}
#'   \item{\code{grp_ordered}}{Displayed groups in importance order.}
#' }
#'
#' @details
#' Requires \pkg{ggplot2} and \pkg{grid}. \pkg{ggbeeswarm} is strongly
#' recommended (falls back to \code{geom_jitter} if unavailable).
#' \pkg{patchwork} is required for \code{level = "both"}.
#'
#' @seealso \code{\link{treeowen}}, \code{\link{treeowen_importance}},
#'   \code{\link{treeowen_hierarchical_beeswarm}}
#'
#' @examples
#' \dontrun{
#' result <- treeowen(unified, X, groups)
#' p <- treeowen_beeswarm(result, level = "feature", top_n_feature = 20L)
#' print(p)
#' out <- treeowen_beeswarm(result, level = "both")
#' print(out$combined)
#' }
#'
#' @export
treeowen_beeswarm <- function(
    result,
    level              = c("feature", "group", "both"),
    group_agg          = c("sum_abs", "sum_abs_feat"),
    normalize_imp      = FALSE,
    top_n_feature      = 10L,
    top_n_group        = 10L,
    group_color_stat   = c("sum", "mean", "custom"),
    group_color_fn     = NULL,
    color_low          = "#0052A5",
    color_high         = "#DC2626",
    color_alpha        = 0.6,
    qlims              = c(0.05, 0.95),
    point_size         = 1.6,
    point_alpha        = 0.55,
    quasirandom_width  = 0.4,
    title_feature      = NULL,
    title_group        = NULL,
    xlab_feature       = "Owen Value",
    xlab_group         = "Owen Value",
    legend_label_feat  = "Value",
    legend_label_grp   = NULL,
    plot_title_size    = 16,
    axis_title_x_size  = 11,
    axis_text_x_size   = 10,
    axis_text_y_size   = 12,
    legend_text_size   = 9,
    legend_title_size  = 9,
    margin_t           = 6,
    margin_r           = 6,
    margin_b           = 6,
    margin_l           = 8,
    legend_position    = "bottom",
    legend_direction   = "horizontal",
    legend_barwidth    = grid::unit(170, "pt"),
    legend_barheight   = grid::unit(10,  "pt"),
    verbose            = FALSE
) {
  if (!.pkg_ok("ggplot2"))
    stop("ggplot2 is required for treeowen_beeswarm().")
  if (!.pkg_ok("grid"))
    stop("grid is required for treeowen_beeswarm(). It ships with base R.")

  has_beeswarm  <- .pkg_ok("ggbeeswarm")
  has_scales    <- .pkg_ok("scales")
  has_patchwork <- .pkg_ok("patchwork")

  if (!has_beeswarm && isTRUE(verbose))
    message("[treeowen_beeswarm] ggbeeswarm not found; using geom_jitter fallback.")
  if (!has_scales && isTRUE(verbose))
    message("[treeowen_beeswarm] scales not found; using internal oob fallback.")

  level            <- match.arg(level)
  group_agg        <- match.arg(group_agg)
  group_color_stat <- match.arg(group_color_stat)

  if (level == "both" && !has_patchwork)
    stop("patchwork is required for level='both'. Install with install.packages('patchwork').")

  if (!inherits(result, c("treeowen_result", "list")) ||
      is.null(result$owens) || is.null(result$observations) || is.null(result$groups))
    stop("result must be the output of treeowen() containing $owens, $observations, $groups.")

  stopifnot(length(qlims) == 2, is.numeric(qlims), all(is.finite(qlims)),
            all(qlims >= 0), all(qlims <= 1), qlims[1] < qlims[2])
  if (!is.numeric(point_alpha) || point_alpha < 0 || point_alpha > 1)
    stop("point_alpha must be in [0, 1].")
  if (!is.numeric(color_alpha) || color_alpha < 0 || color_alpha > 1)
    stop("color_alpha must be in [0, 1].")
  if (group_color_stat == "custom" && !is.function(group_color_fn))
    stop("group_color_stat='custom' requires group_color_fn to be a function.")
  if (!is.null(top_n_feature) && (!is.numeric(top_n_feature) || top_n_feature < 1))
    stop("top_n_feature must be a positive integer or NULL.")
  if (!is.null(top_n_group) && (!is.numeric(top_n_group) || top_n_group < 1))
    stop("top_n_group must be a positive integer or NULL.")

  imp_type <- switch(level, feature = "feature", group = "group", both = "both")
  imp <- treeowen_importance(result, type = imp_type, group_agg = group_agg,
                              sort = TRUE, normalize = normalize_imp)

  ow_mat        <- as.matrix(result$owens)
  x_mat         <- as.matrix(result$observations)
  n             <- nrow(ow_mat)
  p             <- ncol(ow_mat)
  group_attr    <- imp$group_attr
  group_pos     <- imp$group_pos
  K             <- imp$K
  grp_names_all <- imp$group_names

  .oob_squish <- if (has_scales) {
    scales::squish
  } else {
    function(x, range = c(0, 1), only.finite = TRUE) pmin(pmax(x, range[1]), range[2])
  }

  .safe_lims <- function(lims) {
    if (!all(is.finite(lims)) || diff(lims) < .Machine$double.eps)
      c(lims[1] - 0.5, lims[1] + 0.5)
    else lims
  }

  col_low_a  <- grDevices::adjustcolor(color_low,  alpha.f = color_alpha)
  col_high_a <- grDevices::adjustcolor(color_high, alpha.f = color_alpha)

  .beeswarm_geom <- function() {
    if (has_beeswarm) {
      ggbeeswarm::geom_quasirandom(groupOnX = TRUE, width = quasirandom_width,
                                    size = point_size, alpha = point_alpha)
    } else {
      ggplot2::geom_jitter(height = quasirandom_width, width = 0,
                           size = point_size, alpha = point_alpha)
    }
  }

  .color_scale <- function(label, lims) {
    lims <- .safe_lims(lims)
    ggplot2::scale_color_gradient(
      low = col_low_a, high = col_high_a, limits = lims, oob = .oob_squish,
      name = label,
      guide = ggplot2::guide_colorbar(direction = legend_direction,
                                       barwidth = legend_barwidth,
                                       barheight = legend_barheight,
                                       ticks = FALSE))
  }

  .common_theme <- function(hide_left_y = FALSE, hide_right_y = FALSE) {
    th <- ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(size = plot_title_size, face = "bold"),
        axis.title.x     = ggplot2::element_text(size = axis_title_x_size),
        axis.text.x      = ggplot2::element_text(size = axis_text_x_size),
        axis.text.y      = ggplot2::element_text(size = axis_text_y_size, face = "italic"),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin      = ggplot2::margin(t = margin_t, r = margin_r, b = margin_b, l = margin_l),
        legend.position  = legend_position, legend.direction = legend_direction,
        legend.box       = legend_direction, legend.box.just = "center",
        legend.text      = ggplot2::element_text(size = legend_text_size),
        legend.title     = ggplot2::element_text(size = legend_title_size))
    if (hide_left_y)
      th <- th + ggplot2::theme(
        plot.margin = ggplot2::margin(t = margin_t, r = margin_r, b = margin_b, l = 2),
        axis.text.y.left = ggplot2::element_blank(), axis.ticks.y.left = ggplot2::element_blank())
    if (hide_right_y)
      th <- th + ggplot2::theme(
        plot.margin = ggplot2::margin(t = margin_t, r = 2, b = margin_b, l = margin_l),
        axis.text.y.right = ggplot2::element_blank(), axis.ticks.y.right = ggplot2::element_blank())
    th
  }

  .apply_top_n <- function(ordered_vec, top_n) {
    if (!is.null(top_n) && is.finite(as.numeric(top_n)) && as.numeric(top_n) > 0)
      head(ordered_vec, as.integer(top_n))
    else ordered_vec
  }

  if (is.null(legend_label_grp))
    legend_label_grp <- switch(group_color_stat,
      "sum"    = "Value",
      "mean"   = "Value",
      "custom" = "Group\nstat")

  # ── feature data ─────────────────────────────────────────────────────────────
  feat_ordered <- NULL; df_long <- NULL; feat_lims <- NULL
  if (level %in% c("feature", "both")) {
    feat_ordered <- .apply_top_n(imp$feature$feature, top_n_feature)
    rows <- vector("list", length(feat_ordered))
    for (ii in seq_along(feat_ordered)) {
      fn  <- feat_ordered[ii]; col <- match(fn, colnames(ow_mat)); if (is.na(col)) next
      xcol <- if (fn %in% colnames(x_mat)) as.numeric(x_mat[, fn]) else rep(NA_real_, n)
      rows[[ii]] <- data.frame(feature = fn, owen_value = as.numeric(ow_mat[, col]),
                                feat_val = xcol, stringsAsFactors = FALSE)
    }
    df_long <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(df_long) || nrow(df_long) == 0)
      stop("No valid features found for feature beeswarm. Check top_n_feature.")
    df_long$feature <- factor(df_long$feature, levels = rev(feat_ordered))
    feat_lims <- .safe_lims(stats::quantile(df_long$feat_val, probs = qlims,
                                             na.rm = TRUE, names = FALSE, type = 7))
  }

  # ── group colour stat ─────────────────────────────────────────────────────────
  color_stat_mat <- NULL; grp_ordered <- NULL; df_glong <- NULL; grp_lims <- NULL
  if (level %in% c("group", "both")) {
    color_stat_mat <- matrix(NA_real_, n, K, dimnames = list(NULL, grp_names_all))
    for (k in seq_len(K)) {
      gp <- group_pos[[k]]; if (!length(gp)) { color_stat_mat[, k] <- 0; next }
      gp <- gp[gp >= 1L & gp <= p]; xsub <- x_mat[, gp, drop = FALSE]
      color_stat_mat[, k] <- switch(group_color_stat,
        "sum"    = rowSums(xsub),
        "mean"   = rowMeans(xsub),
        "custom" = {
          res <- tryCatch(as.numeric(group_color_fn(x_mat, gp)),
                          error = function(e) { warning(sprintf("[treeowen_beeswarm] group_color_fn error: %s", conditionMessage(e))); rep(NA_real_, n) })
          if (length(res) != n) rep(NA_real_, n) else res
        })
    }
    grp_ordered <- .apply_top_n(imp$group$group, top_n_group)
    rows_g <- vector("list", length(grp_ordered))
    for (ii in seq_along(grp_ordered)) {
      gn <- grp_ordered[ii]; col <- match(gn, grp_names_all); if (is.na(col)) next
      rows_g[[ii]] <- data.frame(group = gn, owen_value = as.numeric(group_attr[, col]),
                                  color_stat = as.numeric(color_stat_mat[, col]),
                                  stringsAsFactors = FALSE)
    }
    df_glong <- do.call(rbind, Filter(Negate(is.null), rows_g))
    if (is.null(df_glong) || nrow(df_glong) == 0)
      stop("No valid groups found for group beeswarm. Check top_n_group.")
    df_glong$group <- factor(df_glong$group, levels = rev(grp_ordered))
    grp_lims <- .safe_lims(stats::quantile(df_glong$color_stat, probs = qlims,
                                            na.rm = TRUE, names = FALSE, type = 7))
  }

  # ── standalone feature beeswarm ───────────────────────────────────────────────
  p_feat <- NULL
  if (level %in% c("feature", "both")) {
    p_feat <- ggplot2::ggplot(df_long, ggplot2::aes(x = owen_value, y = feature, color = feat_val)) +
      .beeswarm_geom() + .color_scale(legend_label_feat, feat_lims) +
      ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
      ggplot2::labs(title = title_feature, x = xlab_feature, y = NULL, color = "") +
      .common_theme()
  }

  # ── standalone group beeswarm ────────────────────────────────────────────────
  p_grp <- NULL
  if (level %in% c("group", "both")) {
    p_grp <- ggplot2::ggplot(df_glong, ggplot2::aes(x = owen_value, y = group, color = color_stat)) +
      .beeswarm_geom() + .color_scale(legend_label_grp, grp_lims) +
      ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
      ggplot2::labs(title = title_group, x = xlab_group, y = NULL, color = "") +
      .common_theme()
  }

  if (level == "feature") return(p_feat)
  if (level == "group")   return(p_grp)

  # ── both mode ────────────────────────────────────────────────────────────────
  n_feat <- length(feat_ordered); n_grp <- length(grp_ordered)

  p_grp_b <- ggplot2::ggplot(df_glong, ggplot2::aes(x = owen_value, y = group, color = color_stat)) +
    .beeswarm_geom() + .color_scale(legend_label_grp, grp_lims) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
    ggplot2::labs(title = title_group, x = xlab_group, y = NULL, color = "") +
    .common_theme(hide_right_y = TRUE)

  p_feat_b <- ggplot2::ggplot(df_long, ggplot2::aes(x = owen_value, y = feature, color = feat_val)) +
    .beeswarm_geom() + .color_scale(legend_label_feat, feat_lims) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
    ggplot2::labs(title = title_feature, x = xlab_feature, y = NULL, color = "") +
    ggplot2::scale_y_discrete(position = "right") + .common_theme(hide_left_y = TRUE)

  # group | feature side-by-side (no connector panel)
  combined <- (p_grp_b + p_feat_b) +
    patchwork::plot_layout(widths = c(1, 1))

  out <- list(combined = combined, feature = p_feat, group = p_grp,
              feat_ordered = feat_ordered, grp_ordered = grp_ordered)
  class(out) <- c("treeowen_plots", "list")
  out
}


#' @export
print.treeowen_plots <- function(x, ...) {
  if (!is.null(x$combined)) print(x$combined)
  else {
    if (!is.null(x$feature)) print(x$feature)
    if (!is.null(x$group))   print(x$group)
  }
  invisible(x)
}
