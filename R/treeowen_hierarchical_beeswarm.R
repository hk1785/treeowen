# treeowen/R/treeowen_hierarchical_beeswarm.R
# Exported function: treeowen_hierarchical_beeswarm()
# General-purpose hierarchical beeswarm for any grouped feature set.
# Group and feature labels are displayed as-is (no domain-specific transforms).

#' Hierarchical Beeswarm Plot of Owen Values (Group × Feature)
#'
#' Produces a hierarchical beeswarm plot where each group and its nested
#' features are displayed as individual rows, assembled into a multi-column
#' layout using \pkg{patchwork}. Groups are ordered by decreasing group
#' importance; features within each group by decreasing feature importance.
#'
#' @param ow_result An object of class \code{"treeowen_result"} returned by
#'   \code{\link{treeowen}}.
#' @param imp An object of class \code{"treeowen_importance"} returned by
#'   \code{\link{treeowen_importance}} (must include both \code{$feature} and
#'   \code{$group}).
#' @param top_n_group Integer or \code{NULL}. Top groups to display by group
#'   importance. \code{NULL} (default) displays all groups.
#' @param top_n_feature Integer or \code{NULL}. Restrict features to those among
#'   the top \code{top_n_feature} by feature importance that also belong to the
#'   selected groups. \code{NULL} (default) displays all features.
#' @param n_col Integer. Number of layout columns. Default \code{2L}.
#' @param group_color_stat Character. How to compute the colour statistic for
#'   each group row: \code{"sum"} (default) or \code{"mean"} of the raw
#'   feature values within the group.
#' @param color_low,color_high Character. Hex colours for the low and high
#'   ends of the colour gradient. Defaults: blue \code{"#0052A5"} and red
#'   \code{"#DC2626"}.
#' @param color_alpha Numeric in \eqn{[0,1]}. Colour gradient transparency.
#'   Default \code{0.65}.
#' @param qlims Numeric vector of length 2. Quantile limits for colour
#'   clipping. Default \code{c(0.05, 0.95)}.
#' @param point_size_grp,point_size_feat Numeric. Point sizes for group and
#'   feature rows respectively. Defaults: \code{2.2} / \code{1.2}.
#' @param point_alpha Numeric in \eqn{[0,1]}. Point transparency.
#'   Default \code{0.50}.
#' @param row_h_grp,row_h_feat Numeric. Row heights (inches) for group and
#'   feature rows respectively. Defaults: \code{0.55} / \code{0.22}.
#' @param width_in Numeric. Total plot width in inches. Default \code{14}.
#' @param max_h_per_page Numeric. Maximum plot height (inches) per page.
#'   When the total height exceeds this, the output is split across multiple
#'   pages. Default \code{25}.
#' @param feat_axis_size,grp_axis_size Numeric. y-axis text sizes (pt) for
#'   feature and group rows respectively. Defaults: \code{7} / \code{8.5}.
#' @param xlab Character. x-axis label. Default \code{"Owen Value"}.
#' @param show_colorbar Logical. Append a standalone colour-bar strip at the
#'   bottom of each page. Default \code{TRUE}.
#' @param colorbar_h Numeric. Height (inches) of the colour-bar strip.
#'   Default \code{0.28}.
#' @param lname Character. Prefix for output file names. Only used when
#'   \code{save_path} is not \code{NULL}. Default \code{"treeowen"}.
#' @param save_path Character or \code{NULL}. Directory for PDF/PNG output.
#'   \code{NULL} (default) suppresses saving; patchwork objects are returned.
#' @param dpi Integer. PNG resolution. Default \code{300L}.
#' @param verbose Logical. Print diagnostic messages. Default \code{FALSE}.
#'
#' @return A list of \pkg{patchwork} objects, one per page. Returned
#'   invisibly when \code{save_path} is not \code{NULL}.
#'
#' @details
#' Requires \pkg{ggplot2}, \pkg{ggbeeswarm}, and \pkg{patchwork}.
#'
#' Group and feature labels are displayed exactly as they appear in
#' \code{ow_result$groups} — no string transformation is applied. When
#' \code{show_colorbar = TRUE} (default), a compact horizontal colour-bar
#' strip is appended below the beeswarm grid on each page, showing the
#' gradient scale used for all points.
#'
#' @seealso \code{\link{treeowen}}, \code{\link{treeowen_importance}},
#'   \code{\link{treeowen_beeswarm}}
#'
#' @examples
#' \dontrun{
#' result <- treeowen(unified, X, groups)
#' imp    <- treeowen_importance(result, type = "both")
#' pages  <- treeowen_hierarchical_beeswarm(
#'   ow_result   = result, imp = imp,
#'   top_n_group = 10L, top_n_feature = 100L, n_col = 2L
#' )
#' print(pages[[1]])
#' # Save to disk
#' treeowen_hierarchical_beeswarm(
#'   ow_result  = result, imp = imp,
#'   save_path  = "output/", lname = "mymodel"
#' )
#' }
#'
#' @export
treeowen_hierarchical_beeswarm <- function(
    ow_result,
    imp,
    top_n_group      = NULL,
    top_n_feature       = NULL,
    n_col            = 2L,
    group_color_stat = "sum",
    color_low        = "#0052A5",
    color_high       = "#DC2626",
    color_alpha      = 0.65,
    qlims            = c(0.05, 0.95),
    point_size_grp   = 2.2,
    point_size_feat  = 1.2,
    point_alpha      = 0.50,
    row_h_grp        = 0.55,
    row_h_feat       = 0.22,
    width_in         = 14,
    max_h_per_page   = 25,
    feat_axis_size   = 7,
    grp_axis_size    = 8.5,
    xlab             = "Owen Value",
    show_colorbar    = TRUE,
    colorbar_h       = 0.28,
    lname            = "treeowen",
    save_path        = NULL,
    dpi              = 300L,
    verbose          = FALSE
) {
  if (!requireNamespace("ggplot2",    quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) stop("ggbeeswarm is required.")
  if (!requireNamespace("patchwork",  quietly = TRUE)) stop("patchwork is required.")

  ow_mat    <- as.matrix(ow_result$owens)
  x_mat     <- as.matrix(ow_result$observations)
  groups    <- ow_result$groups
  n         <- nrow(ow_mat)
  feat_cols <- colnames(ow_mat)

  # ── 1. group order ─────────────────────────────────────────────────────────
  grp_imp   <- imp$group
  grp_order <- if (is.null(top_n_group) || !is.finite(top_n_group) || top_n_group < 1L)
    grp_imp$group
  else
    head(grp_imp$group, as.integer(top_n_group))
  K <- length(grp_order)
  if (isTRUE(verbose))
    message(sprintf("[hierarchical] Displaying %d / %d groups", K, nrow(grp_imp)))

  # ── 2. feature → group mapping ─────────────────────────────────────────────
  feat_imp <- imp$feature
  feat2grp <- setNames(character(length(feat_cols)), feat_cols)
  for (gn in names(groups))
    for (fn in groups[[gn]])
      if (fn %in% feat_cols) feat2grp[fn] <- gn

  # ── 3. features per group (feature importance order) ───────────────────────
  top_feats <- if (is.null(top_n_feature) || !is.finite(top_n_feature) || top_n_feature < 1L)
    NULL
  else
    head(feat_imp$feature, as.integer(top_n_feature))

  feat_by_grp <- lapply(grp_order, function(gn) {
    feats <- feat_imp$feature[feat2grp[feat_imp$feature] == gn]
    if (!is.null(top_feats)) feats <- intersect(feats, top_feats)
    feats
  })
  names(feat_by_grp) <- grp_order

  # ── 4. global colour range ─────────────────────────────────────────────────
  all_vis <- intersect(feat_cols, colnames(x_mat))
  x_vals  <- as.numeric(x_mat[, all_vis, drop = FALSE])
  clim_lo <- stats::quantile(x_vals, qlims[1], na.rm = TRUE)
  clim_hi <- stats::quantile(x_vals, qlims[2], na.rm = TRUE)
  if (abs(clim_hi - clim_lo) < .Machine$double.eps) {
    clim_lo <- clim_lo - 0.5; clim_hi <- clim_hi + 0.5
  }
  sq       <- function(x) pmin(pmax(x, clim_lo), clim_hi)
  col_lo_a <- grDevices::adjustcolor(color_low,  alpha.f = color_alpha)
  col_hi_a <- grDevices::adjustcolor(color_high, alpha.f = color_alpha)

  grp_attr      <- imp$group_attr
  grp_names_all <- colnames(grp_attr)

  # oob helper (scales optional)
  .oob <- if (requireNamespace("scales", quietly = TRUE)) {
    scales::squish
  } else {
    function(x, range = c(0, 1), ...) pmin(pmax(x, range[1]), range[2])
  }

  # ── 5. group colour statistic ──────────────────────────────────────────────
  .grp_color_vals <- function(gn) {
    gf <- intersect(groups[[gn]], colnames(x_mat))
    if (!length(gf)) return(rep(0, n))
    if (group_color_stat == "sum") {
      if (length(gf) == 1L) as.numeric(x_mat[, gf]) else rowSums(x_mat[, gf, drop = FALSE])
    } else {
      if (length(gf) == 1L) as.numeric(x_mat[, gf]) else rowMeans(x_mat[, gf, drop = FALSE])
    }
  }

  # ── 6. shared x range (visible data only) ─────────────────────────────────
  vis_feat_ids     <- unlist(feat_by_grp, use.names = FALSE)
  vis_feat_cols_ow <- intersect(vis_feat_ids, colnames(ow_mat))
  vis_grp_cols     <- intersect(grp_order, colnames(grp_attr))
  all_ow <- if (length(vis_feat_cols_ow))
    as.numeric(ow_mat[, vis_feat_cols_ow, drop = FALSE]) else numeric(0)
  grp_ow <- if (length(vis_grp_cols))
    as.numeric(grp_attr[, vis_grp_cols, drop = FALSE]) else numeric(0)
  xlim   <- range(c(all_ow, grp_ow), na.rm = TRUE)
  xlim   <- xlim + diff(xlim) * c(-0.05, 0.05)

  # ── 7. single-row beeswarm plot factory ───────────────────────────────────
  # Labels are passed in as-is (already the final display strings).
  # is_group controls visual styling: group rows have bold text and grey background.
  .make_row_plot <- function(ow_vec, color_vec, label, is_group,
                              show_xaxis = FALSE) {
    df <- data.frame(ov = ow_vec, col = sq(color_vec), y = label,
                     stringsAsFactors = FALSE)
    df$y_fac <- factor(df$y, levels = label)
    ps  <- if (is_group) point_size_grp else point_size_feat
    ts  <- if (is_group) grp_axis_size  else feat_axis_size
    fw  <- if (is_group) "bold"         else "plain"
    bg  <- if (is_group) "grey92"       else "white"

    ggplot2::ggplot(df, ggplot2::aes(x = ov, y = y_fac, color = col)) +
      ggbeeswarm::geom_quasirandom(groupOnX = FALSE, width = 0.35,
                                    size = ps, alpha = point_alpha) +
      ggplot2::geom_vline(xintercept = 0, color = "grey50",
                          linewidth = 0.3, linetype = "dashed") +
      ggplot2::scale_color_gradient(
        low    = col_lo_a,
        high   = col_hi_a,
        limits = c(clim_lo, clim_hi),
        oob    = .oob,
        guide  = "none") +
      ggplot2::scale_y_discrete() +
      ggplot2::scale_x_continuous(
        limits = xlim,
        expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(x = if (show_xaxis) xlab else NULL, y = NULL) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(
        axis.text.y        = ggplot2::element_text(size = ts, face = fw, hjust = 1),
        axis.text.x        = if (show_xaxis)
                               ggplot2::element_text(size = 7)
                             else ggplot2::element_blank(),
        axis.ticks.x       = if (show_xaxis) ggplot2::element_line()
                             else ggplot2::element_blank(),
        axis.title.x       = ggplot2::element_text(size = 8),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.background   = ggplot2::element_rect(fill = bg),
        legend.position    = "none",
        plot.margin        = ggplot2::margin(
                               t = if (is_group) 4 else 0, r = 10, b = 0, l = 4))
  }

  # ── 8. standalone colorbar strip ──────────────────────────────────────────
  .make_colorbar_strip <- function() {
    col_leg_name <- "Value"
    xs    <- seq(clim_lo, clim_hi, length.out = 100)
    df_cb <- data.frame(x = xs, y = 0.5)
    ggplot2::ggplot(df_cb, ggplot2::aes(x = x, y = y, fill = x)) +
      ggplot2::geom_tile(height = 0.2) +
      ggplot2::scale_fill_gradient(
        low    = col_lo_a,
        high   = col_hi_a,
        limits = c(clim_lo, clim_hi),
        oob    = .oob,
        guide  = "none") +
      ggplot2::scale_x_continuous(
        limits = c(clim_lo, clim_hi),
        expand = ggplot2::expansion(0),
        breaks = c(clim_lo, clim_hi),
        labels = c("Low", "High")) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(0)) +
      ggplot2::labs(x = col_leg_name, y = NULL) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(
        axis.text.y    = ggplot2::element_blank(),
        axis.ticks.y   = ggplot2::element_blank(),
        axis.text.x    = ggplot2::element_text(size = 7),
        axis.title.x   = ggplot2::element_text(size = 8, hjust = 0.5),
        panel.grid     = ggplot2::element_blank(),
        panel.border   = ggplot2::element_blank(),
        plot.margin    = ggplot2::margin(t = 2, r = 10, b = 2, l = 4))
  }

  # ── 9. flatten all rows ────────────────────────────────────────────────────
  # Labels are the raw group/feature names from the result object.
  # No string transformation is applied — users control naming upstream.
  all_plots <- list(); all_heights <- numeric(0)

  for (ii in seq_len(K)) {
    gn      <- grp_order[ii]
    feats_g <- feat_by_grp[[gn]]
    col_grp <- match(gn, grp_names_all)
    ow_grp  <- if (!is.na(col_grp)) as.numeric(grp_attr[, col_grp]) else rep(0, n)
    cs_grp  <- .grp_color_vals(gn)
    all_plots   <- c(all_plots, list(
      list(ow_vec = ow_grp, color_vec = cs_grp, label = gn, is_group = TRUE)))
    all_heights <- c(all_heights, row_h_grp)

    for (jj in seq_along(feats_g)) {
      fn      <- feats_g[jj]
      col_fn  <- match(fn, colnames(ow_mat))
      ow_fn   <- if (!is.na(col_fn)) as.numeric(ow_mat[, col_fn]) else rep(0, n)
      xcol_fn <- if (fn %in% colnames(x_mat)) as.numeric(x_mat[, fn]) else rep(0, n)
      all_plots   <- c(all_plots, list(
        list(ow_vec = ow_fn, color_vec = xcol_fn, label = fn, is_group = FALSE)))
      all_heights <- c(all_heights, row_h_feat)
    }
  }

  total_rows <- length(all_plots)
  if (isTRUE(verbose))
    message(sprintf("[hierarchical] %d groups, %d total rows", K, total_rows))

  # ── 10. page splitting ────────────────────────────────────────────────────
  pages <- list(); cur_rows <- integer(0); cur_h <- 0
  for (ri in seq_len(total_rows)) {
    rh <- all_heights[ri]
    if (cur_h + rh > max_h_per_page && length(cur_rows) > 0) {
      pages <- c(pages, list(cur_rows)); cur_rows <- integer(0); cur_h <- 0
    }
    cur_rows <- c(cur_rows, ri); cur_h <- cur_h + rh
  }
  if (length(cur_rows) > 0) pages <- c(pages, list(cur_rows))

  n_pages <- length(pages)
  if (isTRUE(verbose))
    message(sprintf("[hierarchical] %d total rows → %d page(s)", total_rows, n_pages))

  # ── 11. assemble pages ────────────────────────────────────────────────────
  n_col      <- max(1L, as.integer(n_col))
  page_plots <- vector("list", n_pages)

  for (pi in seq_len(n_pages)) {
    page_row_ids <- pages[[pi]]; n_rows_page <- length(page_row_ids)
    rows_per_col <- ceiling(n_rows_page / n_col)

    col_row_ids <- lapply(seq_len(n_col), function(cc) {
      start <- (cc - 1L) * rows_per_col + 1L
      end   <- min(cc * rows_per_col, n_rows_page)
      if (start > n_rows_page) return(integer(0))
      page_row_ids[start:end]
    })

    col_strips <- lapply(seq_len(n_col), function(cc) {
      rids <- col_row_ids[[cc]]
      if (length(rids) == 0L)
        return(list(pw = patchwork::plot_spacer(), h = 1))
      col_plots <- list(); col_heights <- numeric(0)
      for (jj in seq_along(rids)) {
        ri       <- rids[jj]; row_info <- all_plots[[ri]]
        show_x   <- (jj == length(rids))
        p <- .make_row_plot(
          ow_vec    = row_info$ow_vec,
          color_vec = row_info$color_vec,
          label     = row_info$label,
          is_group  = row_info$is_group,
          show_xaxis = show_x)
        col_plots   <- c(col_plots,   list(p))
        col_heights <- c(col_heights, all_heights[ri])
      }
      strip_pw <- patchwork::wrap_plots(col_plots, ncol = 1L) +
        patchwork::plot_layout(heights = col_heights)
      list(pw = strip_pw, h = sum(col_heights))
    })

    ph_body <- max(sapply(col_strips, `[[`, "h"))
    pw_body <- patchwork::wrap_plots(
      lapply(col_strips, `[[`, "pw"),
      ncol   = n_col,
      widths = rep(1, n_col))

    # Standalone colorbar strip — centered with empty spacers on both sides,
    # so the strip occupies only the middle ~70% of the page width
    # (i.e. width is 30% narrower than the beeswarm body above it).
    if (isTRUE(show_colorbar)) {
      cb         <- .make_colorbar_strip()
      cb_spacer  <- patchwork::plot_spacer()
      cb_row     <- patchwork::wrap_plots(
        list(cb_spacer, cb, cb_spacer),
        nrow   = 1L,
        widths = c(0.15, 0.70, 0.15))
      ph <- ph_body + colorbar_h
      pw <- patchwork::wrap_plots(
        list(pw_body, cb_row),
        ncol    = 1L,
        heights = c(ph_body, colorbar_h))
    } else {
      ph <- ph_body
      pw <- pw_body
    }

    page_plots[[pi]] <- pw

    if (!is.null(save_path)) {
      fname_pdf <- file.path(save_path,
        if (n_pages == 1L)
          sprintf("fig_%s_beeswarm_hierarchical.pdf", lname)
        else
          sprintf("fig_%s_beeswarm_hierarchical_%d.pdf", lname, pi))
      fname_png <- sub("\\.pdf$", ".png", fname_pdf)
      ggplot2::ggsave(fname_pdf, pw, width = width_in, height = ph,
                      limitsize = FALSE)
      ggplot2::ggsave(fname_png, pw, width = width_in, height = ph,
                      limitsize = FALSE, dpi = dpi)
      if (isTRUE(verbose))
        message(sprintf("  %s saved (%.1f x %.1f in)", fname_pdf, width_in, ph))
    }
  }

  if (!is.null(save_path)) invisible(page_plots) else page_plots
}
