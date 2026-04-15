\name{treeowen-package}
\alias{treeowen-package}
\alias{_PACKAGE}
\docType{package}
\title{treeowen: Game-Theoretic Explanations for Tree-Based Ensembles via Owen Values}
\description{
Provides exact and Monte Carlo algorithms for computing Owen values in
tree-based ensemble models (XGBoost, LightGBM, Ranger). The Owen value
extends the Shapley value to settings where features are organised into
\emph{a priori} groups, allocating attribution in two stages: first across
groups, then within each group. Hierarchy-guided group aggregation with
tree-aware dynamic programming avoids explicit enumeration of outer group
subsets. Global importance measures and visualisation tools (flat and
hierarchical beeswarm plots) are also provided.
}
\section{Main workflow}{
\enumerate{
  \item Fit a tree model and unify it with \pkg{treeshap}:
        \code{treeshap::xgboost.unify()},
        \code{treeshap::lightgbm.unify()}, or
        \code{treeshap::ranger.unify()}.
  \item Define a feature partition \code{groups} — a named list where each
        element is a character vector of feature names belonging to that
        group.
  \item Compute Owen values with \code{\link{treeowen}()}.
  \item Summarise importance with \code{\link{treeowen_importance}()}.
  \item Visualise with \code{\link{treeowen_beeswarm}()} (flat) or
        \code{\link{treeowen_hierarchical_beeswarm}()} (hierarchical).
}
}
\section{Key functions}{
\describe{
  \item{\code{\link{treeowen}}}{Compute Owen values (exact or Monte Carlo).}
  \item{\code{\link{treeowen_importance}}}{Feature and group importance scores.}
  \item{\code{\link{treeowen_beeswarm}}}{Flat beeswarm plots.}
  \item{\code{\link{treeowen_hierarchical_beeswarm}}}{Hierarchical beeswarm.}
  \item{\code{\link{treeowen_exact_enum}}}{Reference model-agnostic enumerator.}
  \item{\code{\link{treeowen_exact_tuvalues}}}{Reference tree-aware enumerator.}
  \item{\code{\link{build_hierarchy_tree_binary}}}{Custom binary auxiliary tree.}
  \item{\code{\link{build_hierarchy_tree_from_layers}}}{Auxiliary tree from layer spec.}
  \item{\code{\link{clear_inner_enum_cache}}}{Free the memoisation cache.}
  \item{\code{\link{TREEOWEN_AUTO_EXACT_MAX_M}}}{Exact/MC switching threshold.}
  \item{\code{\link{TREEOWEN_INNER_BITMASK_DEFAULT}}}{Bitmask enumerator limit.}
}
}
\section{Datasets}{
\describe{
  \item{\code{\link{immuno}}}{Gut microbiome relative abundances from 219
    cancer patients treated with immune checkpoint inhibitor therapy, pooled
    across five independent cohorts.  Includes a pre-built group partition
    mapping species to genera.}
}
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\author{Hyunwook Koh \email{hyunwook.koh@stonybrook.edu}}
\seealso{
\url{https://github.com/hk1785/treeowen}
}
