\name{treeowen-package}
\alias{treeowen-package}
\alias{_PACKAGE}
\docType{package}
\title{treeowen: Game-Theoretic Explanations for Tree-Based Ensembles via Owen Values}
\description{
treeowen computes Owen values for tree-based ensemble models including
XGBoost, LightGBM, and Ranger.

Owen values extend standard SHAP (Shapley) values to settings where
features are organized into groups. The attribution is allocated in
two stages: first, each group receives a score reflecting how much
the group as a whole contributes to the prediction; second, that score
is divided among the individual features within the group. The result
satisfies the efficiency property: the feature-level Owen values sum
to the model prediction minus the average prediction across all
observations.

The package includes exact and Monte Carlo algorithms, global
importance measures, and flat and hierarchical beeswarm visualizations.
Robust compatibility wrappers (.xgboost_unify_compat,
.lightgbm_unify_compat, .ranger_unify_compat) handle all known
version-compatibility issues across different releases of xgboost,
lightgbm, ranger, and treeshap.
}
\section{Main workflow}{
\enumerate{
 \item Fit a tree model and unify it with \pkg{treeshap}:
 \code{.xgboost_unify_compat()},
 \code{.lightgbm_unify_compat()}, or
 \code{.ranger_unify_compat()}.
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
 across five independent cohorts. Includes a pre-built group partition
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
