\name{treeowen_importance}
\alias{treeowen_importance}
\alias{print.treeowen_importance}
\title{Compute Feature and Group Importance from Owen Values}
\description{
Derives feature-level and group-level importance scores from a
\code{\link{treeowen}} result. Feature importance is the mean absolute
Owen value across observations. Group importance is the mean absolute
sum of the within-group Owen values.
}
\usage{
treeowen_importance(
 result,
 type      = c("both", "feature", "group"),
 group_agg = c("sum_abs", "sum_abs_feat", "both"),
 sort      = TRUE,
 normalize = FALSE
)
}
\arguments{
\item{result}{An object of class \code{"treeowen_result"} returned by
 \code{\link{treeowen}}.}

\item{type}{Character. Which scores to compute.
 \describe{
 \item{\code{"both"}}{(default) Feature and group importance.}
 \item{\code{"feature"}}{Feature importance only.}
 \item{\code{"group"}}{Group importance only.}
 }}

\item{group_agg}{Character. How to aggregate feature-level Owen values
 to the group level.
 \describe{
 \item{\code{"sum_abs"}}{(default) For each observation, sum the
  absolute Owen values of all features in the group, then take the
  mean across observations.}
 \item{\code{"sum_abs_feat"}}{Average of the within-group feature
  importances.}
 \item{\code{"both"}}{Compute both; the \code{importance} column
  uses \code{"sum_abs"} for sorting; an extra column
  \code{importance_sum_abs_feat} is appended.}
 }}

\item{sort}{Logical. Sort results by decreasing importance.
 Default \code{TRUE}.}

\item{normalize}{Logical. Divide all scores by their sum so they add
 up to 1. Default \code{FALSE}.}
}

\details{
Feature importance is the mean absolute Owen value for each feature
across all observations.

Group importance aggregates the feature-level Owen values within each
group. The default aggregation (\code{"sum_abs"}) first sums the
absolute Owen values of all features in the group for each observation,
then takes the mean across observations.

When \code{normalize = TRUE}, all importance values are divided by
their sum so they add up to 1.
}

\value{
An object of class \code{"treeowen_importance"} (a named list):
\describe{
 \item{\code{feature}}{Data frame with columns \code{feature} and
  \code{importance}, sorted by decreasing importance. Present when
  \code{type} is \code{"both"} or \code{"feature"}.}
 \item{\code{group}}{Data frame with columns \code{group} and
  \code{importance}, sorted by decreasing importance. Present when
  \code{type} is \code{"both"} or \code{"group"}.}
 \item{\code{group_attr}}{Numeric matrix (observations x groups)
  of group-level Owen attributions. Always present; used internally
  by the visualization functions.}
 \item{\code{group_names}}{Character vector of group names.}
 \item{\code{feat_names}}{Character vector of feature names.}
}
}

\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{treeowen_beeswarm}},
\code{\link{treeowen_hierarchical_beeswarm}},
\code{\link{print.treeowen_importance}}
}
\examples{
\dontrun{
library(treeowen)
library(xgboost)

set.seed(42)
feat_names <- paste0("F", rep(1:5, times = 20), "G", rep(1:20, each = 5))
X <- as.data.frame(matrix(rnorm(100 * 100), 100, 100,
                           dimnames = list(NULL, feat_names)))
Y <- as.integer(0.8 * X$F1G1 - 0.6 * X$F2G1 + 0.5 * X$F1G3 > 0)
groups <- setNames(lapply(1:20, function(k) paste0("F", 1:5, "G", k)),
                   paste0("G", 1:20))

library(xgboost)
xgb_mod <- xgboost(data = xgb.DMatrix(as.matrix(X), label = Y),
                   nrounds = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
unified <- xgboost_unify_compat(xgb_mod, X)
result  <- treeowen(unified, X, groups, method = "auto",
                    dp_progress = FALSE, verbose = FALSE)

## Feature and group importance
imp <- treeowen_importance(result, type = "both", group_agg = "sum_abs")
print(imp)
head(imp$feature, 5L)
imp$group

## Normalised group importance (sums to 1)
imp_norm <- treeowen_importance(result, type = "group", normalize = TRUE)
stopifnot(abs(sum(imp_norm$group$importance) - 1) < 1e-10)
}
}
