\name{treeowen_hierarchical_beeswarm}
\alias{treeowen_hierarchical_beeswarm}
\title{Hierarchical Beeswarm Plot of Owen Values (Group x Feature)}
\description{
Produces a hierarchical beeswarm plot in which each group and its nested
features are displayed as individual rows, assembled into a multi-column
layout using \pkg{patchwork}. Groups are ordered by decreasing group
importance; features within each group by decreasing feature importance.
An optional standalone colour-bar strip is appended at the bottom of each
page.
}
\usage{
treeowen_hierarchical_beeswarm(
 ow_result,
 imp,
 top_n_group = NULL,
 top_n_feature = NULL,
 n_col = 2L,
 group_color_stat = "sum",
 color_low = "#0052A5",
 color_high = "#DC2626",
 color_alpha = 0.65,
 qlims = c(0.05, 0.95),
 point_size_grp = 2.2,
 point_size_feat = 1.2,
 point_alpha = 0.50,
 row_h_grp = 0.55,
 row_h_feat = 0.22,
 width_in = 14,
 max_h_per_page = 25,
 feat_axis_size = 7,
 grp_axis_size = 8.5,
 xlab = "Owen Value",
 show_colorbar = TRUE,
 colorbar_h = 0.35,
 lname = "treeowen",
 save_path = NULL,
 dpi = 300L,
 verbose = FALSE
)
}
\arguments{
\item{ow_result}{An object of class \code{"treeowen_result"} from
 \code{\link{treeowen}}.}

\item{imp}{An object of class \code{"treeowen_importance"} from
 \code{\link{treeowen_importance}} computed with \code{type = "both"}
 (must contain both \code{$feature} and \code{$group}).}

\item{top_n_group}{Integer or \code{NULL}. Top groups by group importance.
 \code{NULL} (default) shows all.}

\item{top_n_feature}{Integer or \code{NULL}. Restrict features to those
 among the top \code{top_n_feature} by feature importance that also belong
 to the selected groups. \code{NULL} (default) shows all features in the
 selected groups.}

\item{n_col}{Integer. Layout columns. Default \code{2L}.}

\item{group_color_stat}{Character. Colour statistic for group rows:
 \code{"sum"} (default) or \code{"mean"} of raw feature values in the
 group.}

\item{color_low}{Hex colour for low gradient end.
 Default \code{"#0052A5"} (blue).}

\item{color_high}{Hex colour for high gradient end.
 Default \code{"#DC2626"} (red).}

\item{color_alpha}{Numeric in [0,1]. Gradient transparency.
 Default \code{0.65}.}

\item{qlims}{Numeric vector of length 2. Quantile limits for colour
 clipping. Default \code{c(0.05, 0.95)}.}

\item{point_size_grp}{Numeric. Group row point size. Default \code{2.2}.}
\item{point_size_feat}{Numeric. Feature row point size. Default \code{1.2}.}
\item{point_alpha}{Numeric in [0,1]. Point transparency.
 Default \code{0.50}.}
\item{row_h_grp}{Numeric. Group row height (inches). Default \code{0.55}.}
\item{row_h_feat}{Numeric. Feature row height (inches). Default \code{0.22}.}
\item{width_in}{Numeric. Plot width (inches). Default \code{14}.}
\item{max_h_per_page}{Numeric. Maximum height per page (inches).
 Default \code{25}.}
\item{feat_axis_size}{Numeric (pt). Feature row y-axis text.
 Default \code{7}.}
\item{grp_axis_size}{Numeric (pt). Group row y-axis text.
 Default \code{8.5}.}
\item{xlab}{Character. x-axis label. Default \code{"Owen Value"}.}

\item{show_colorbar}{Logical. Append a standalone horizontal colour-bar
 strip at the bottom of each page. Default \code{TRUE}.}

\item{colorbar_h}{Numeric. Colour-bar height (inches). Default \code{0.35}.}

\item{lname}{Character. File-name prefix (used only when \code{save_path}
 is not \code{NULL}). Default \code{"treeowen"}.}

\item{save_path}{Character or \code{NULL}. Directory for PDF/PNG output.
 \code{NULL} (default) returns the list visibly without saving.}

\item{dpi}{Integer. PNG resolution. Default \code{300L}.}

\item{verbose}{Logical. Print diagnostic messages. Default \code{FALSE}.}
}
\value{
A list of \pkg{patchwork} objects, one per page. Each element is the
assembled plot (beeswarm columns + optional colour-bar strip).

Returned visibly when \code{save_path = NULL}; invisibly when
\code{save_path} is set (files are written to disk as
\file{fig_<lname>_beeswarm_hierarchical[_<i>].pdf} and the matching
\file{.png}).
}
\details{
\subsection{Layout}{
Each page is divided into \code{n_col} columns. Rows (one per group or
feature) fill columns top-to-bottom, then left-to-right. Group rows use
a grey background and bold text; feature rows use white background and
plain text. The x-axis label appears only on the last row of each column.
}
\subsection{Colour scale}{
The global colour range is computed from quantiles \code{qlims} of all
raw feature values. Group rows are coloured by \code{group_color_stat}
applied to the raw feature values of the group's members. Feature rows
are coloured by the individual feature's raw value. All rows share the
same gradient.
}
\subsection{Page splitting}{
When cumulative row height exceeds \code{max_h_per_page} the output is
automatically split across pages (one list element per page).
}
\subsection{Dependencies}{
Requires \pkg{ggplot2}, \pkg{ggbeeswarm}, and \pkg{patchwork}.
\pkg{scales} is used for colour clipping when available.
}
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{treeowen_importance}},
\code{\link{treeowen_beeswarm}}
}
\examples{
\dontrun{
library(treeowen)
library(xgboost)
library(ggplot2)
library(ggbeeswarm)   # required
library(patchwork)    # required

## Synthetic data: 10 groups of 5 features each, 50 observations
set.seed(42)
feat_names <- paste0("F", rep(1:5, times = 10), "G", rep(1:10, each = 5))
X <- as.data.frame(matrix(rnorm(50 * 50), 50, 50,
                          dimnames = list(NULL, feat_names)))
Y <- as.integer(0.8 * X$F1G1 - 0.6 * X$F2G1 + 0.5 * X$F1G3 + 0.3 * X$F1G2 > 0)
groups <- setNames(lapply(1:10, function(k) paste0("F", 1:5, "G", k)),
                   paste0("G", 1:10))

xgb_mod <- xgboost(data = xgb.DMatrix(as.matrix(X), label = Y),
                   nrounds = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
result <- treeowen(xgboost_unify_compat(xgb_mod, X), X, groups,
                   method = "auto", dp_progress = FALSE, verbose = FALSE)

## Importance (type = "both" required)
imp <- treeowen_importance(result, type = "both")

## Hierarchical beeswarm: returns one patchwork per page
pages <- treeowen_hierarchical_beeswarm(
  ow_result     = result,
  imp           = imp,
  top_n_group   = 5L,       # top 5 groups
  top_n_feature = NULL,     # all features in those groups
  n_col         = 2L,
  group_color_stat = "sum",
  show_colorbar = TRUE,
  width_in      = 10,
  max_h_per_page = 15,
  verbose       = TRUE
)
cat("Pages:", length(pages), "\n")
print(pages[[1]])

## Save to disk (PDF + PNG)
tmp_dir <- tempdir()
treeowen_hierarchical_beeswarm(
  ow_result  = result,
  imp        = imp,
  top_n_group = 5L,
  n_col      = 2L,
  width_in   = 10,
  lname      = "xgb_synthetic",
  save_path  = tmp_dir,
  dpi        = 150L,
  verbose    = TRUE
)
cat("Saved to:", tmp_dir, "\n")
}
}
