\name{treeowen_beeswarm}
\alias{treeowen_beeswarm}
\title{Beeswarm Plots of Owen Values}
\description{
Produces SHAP-style beeswarm plots of Owen value distributions at the
feature level, the group level, or both panels side by side. Points are coloured by the feature's raw value (feature panel) or by
a group-level summary statistic (group panel).
}
\usage{
treeowen_beeswarm(
 result,
 level = c("feature", "group", "both"),
 group_agg = c("sum_abs", "sum_abs_feat"),
 normalize_imp = FALSE,
 top_n_feature = 10L,
 top_n_group = 10L,
 group_color_stat = c("sum", "mean", "custom"),
 group_color_fn = NULL,
 color_low = "#0052A5",
 color_high = "#DC2626",
 color_alpha = 0.6,
 qlims = c(0.05, 0.95),
 point_size = 1.6,
 point_alpha = 0.55,
 quasirandom_width = 0.4,
 title_feature = NULL,
 title_group = NULL,
 xlab_feature = "Owen Value",
 xlab_group = "Owen Value",
 legend_label_feat = "Value",
 legend_label_grp = NULL,
 plot_title_size = 16,
 axis_title_x_size = 11,
 axis_text_x_size = 10,
 axis_text_y_size = 12,
 legend_text_size = 9,
 legend_title_size = 9,
 margin_t = 6,
 margin_r = 6,
 margin_b = 6,
 margin_l = 8,
 legend_position = "bottom",
 legend_direction = "horizontal",
 legend_barwidth = grid::unit(170, "pt"),
 legend_barheight = grid::unit(10, "pt"),
 verbose = FALSE
)
}
\arguments{
\item{result}{An object of class \code{"treeowen_result"} from
 \code{\link{treeowen}}.}

\item{level}{Character. Panel(s) to produce:
 \describe{
 \item{\code{"feature"}}{(default) One beeswarm row per feature.}
 \item{\code{"group"}}{One beeswarm row per group.}
 \item{\code{"both"}}{\pkg{patchwork}: group | feature side-by-side.
 Requires \pkg{patchwork}.}
 }}

\item{group_agg}{Character. Aggregation for sort order:
 \code{"sum_abs"} (default, mean absolute group attribution) or
 \code{"sum_abs_feat"} (mean of within-group feature importances).
 \code{"both"} is not permitted (ambiguous sort order).}

\item{normalize_imp}{Logical. Normalise importance before sorting.
 Default \code{FALSE}.}

\item{top_n_feature}{Integer or \code{NULL}. Top features to show.
 Default \code{10L}.}

\item{top_n_group}{Integer or \code{NULL}. Top groups to show.
 Default \code{10L}.}

\item{group_color_stat}{Character. Colour statistic for group panel:
 \code{"sum"} (default), \code{"mean"}, or \code{"custom"}.}

\item{group_color_fn}{Function or \code{NULL}. Required when
 \code{group_color_stat = "custom"}: \code{function(x_mat, feat_idx) ->
 numeric(n)}.}

\item{color_low}{Hex colour for low end of gradient.
 Default \code{"#0052A5"} (blue).}

\item{color_high}{Hex colour for high end. Default \code{"#DC2626"} (red).}

\item{color_alpha}{Numeric in [0,1]. Gradient transparency.
 Default \code{0.6}.}

\item{qlims}{Numeric vector of length 2. Quantile limits for colour
 clipping. Default \code{c(0.05, 0.95)}.}

\item{point_size}{Numeric. Point size. Default \code{1.6}.}
\item{point_alpha}{Numeric in [0,1]. Point transparency.
 Default \code{0.55}.}
\item{quasirandom_width}{Numeric. Quasirandom jitter spread.
 Default \code{0.4}.}

\item{title_feature}{Character or \code{NULL}. Feature panel title.}
\item{title_group}{Character or \code{NULL}. Group panel title.}
\item{xlab_feature}{Character. x-axis label, feature panel.
 Default \code{"Owen Value"}.}
\item{xlab_group}{Character. x-axis label, group panel.
 Default \code{"Owen Value"}.}
\item{legend_label_feat}{Character. Colour-bar label, feature panel.
 Default \code{"Value"}.}
\item{legend_label_grp}{Character or \code{NULL}. Colour-bar label,
 group panel. \code{NULL} auto-generates from \code{group_color_stat}.}

\item{plot_title_size}{Numeric (pt). Default \code{16}.}
\item{axis_title_x_size}{Numeric (pt). Default \code{11}.}
\item{axis_text_x_size}{Numeric (pt). Default \code{10}.}
\item{axis_text_y_size}{Numeric (pt). Default \code{12}.}
\item{legend_text_size}{Numeric (pt). Default \code{9}.}
\item{legend_title_size}{Numeric (pt). Default \code{9}.}

\item{margin_t, margin_r, margin_b, margin_l}{Numeric (pt). Plot margins.
 Defaults: \code{6}, \code{6}, \code{6}, \code{8}.}

\item{legend_position}{Character. Default \code{"bottom"}.}
\item{legend_direction}{Character. Default \code{"horizontal"}.}
\item{legend_barwidth}{\code{grid::unit}.
 Default \code{grid::unit(170, "pt")}.}
\item{legend_barheight}{\code{grid::unit}.
 Default \code{grid::unit(10, "pt")}.}

\item{verbose}{Logical. Default \code{FALSE}.}
}
\value{
\describe{
 \item{\code{level = "feature"} or \code{"group"}}{A \code{ggplot2}
 object.}
 \item{\code{level = "both"}}{An object of class \code{"treeowen_plots"}
 (named list):
 \describe{
 \item{\code{combined}}{\pkg{patchwork}: group | feature side-by-side.}
 \item{\code{feature}}{Standalone feature \code{ggplot}.}
 \item{\code{group}}{Standalone group \code{ggplot}.}
 \item{\code{feat_ordered}}{Features in display order.}
 \item{\code{grp_ordered}}{Groups in display order.}
 }}
}
}
\details{
Requires \pkg{ggplot2} and \pkg{grid}. \pkg{ggbeeswarm} is strongly
recommended (falls back to \code{geom_jitter} if unavailable).
\pkg{patchwork} is required for \code{level = "both"}. The two panels
(group and feature) are placed side by side.
\pkg{scales} is used for colour clipping when available.
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{treeowen_importance}},
\code{\link{treeowen_hierarchical_beeswarm}}, \code{\link{print.treeowen_plots}}
}
\examples{
\dontrun{
## ── Synthetic data + XGBoost (see ?treeowen for data setup) ─────────────────
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=20), "G", rep(1:20, each=5))
X <- as.data.frame(matrix(rnorm(100*100), 100, 100, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2
Y <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:20, function(k) paste0("F",1:5,"G",k)), paste0("G",1:20))

library(treeowen)
library(xgboost)
library(ggplot2)
library(ggbeeswarm)  # for geom_quasirandom
library(patchwork)   # required for level = "both"
xgb_mod <- xgboost(data = xgb.DMatrix(as.matrix(X), label=Y),
 nrounds = 50L, max_depth = 3L, eta = 0.1,
 objective = "binary:logistic", verbose = 0L)
result <- treeowen(xgboost_unify_compat(xgb_mod, X), X, groups,
 method="auto", dp_progress=FALSE, verbose=FALSE)

## ── Feature-level beeswarm ──────────────────────────────────────────────────
p_feat <- treeowen_beeswarm(result, level = "feature", top_n_feature = 15L,
 xlab_feature = "Owen Value (XGBoost)")
print(p_feat)

## ── Group-level beeswarm ────────────────────────────────────────────────────
p_grp <- treeowen_beeswarm(result, level = "group", top_n_group = 5L,
 group_color_stat = "sum")
print(p_grp)

## ── Both panels (group | feature side-by-side) ─────────────────────────────
## Requires patchwork
library(patchwork)
out <- treeowen_beeswarm(result, level = "both",
 top_n_feature = 10L, top_n_group = 5L)
print(out$combined) # patchwork figure
print(out$feature) # feature panel alone
print(out$group) # group panel alone
}
}
