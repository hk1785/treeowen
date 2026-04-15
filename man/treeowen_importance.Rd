\name{treeowen_importance}
\alias{treeowen_importance}
\title{Compute Feature and Group Importance from Owen Values}
\description{
Derives feature-level and group-level importance scores from a
\code{\link{treeowen}} result. Feature importance is the mean absolute Owen
value across observations; group importance is the mean absolute sum of
within-group Owen values (group-level attribution).
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
\item{result}{An object of class \code{"treeowen_result"} from
  \code{\link{treeowen}}.}

\item{type}{Character.  Which scores to compute:
  \describe{
    \item{\code{"both"}}{(default) Feature and group importance.}
    \item{\code{"feature"}}{Feature importance only.}
    \item{\code{"group"}}{Group importance only.}
  }}

\item{group_agg}{Character.  How to compute group importance:
  \describe{
    \item{\code{"sum_abs"}}{(default)
      \eqn{I_k = (1/n)\sum_j |\Psi_k^{(j)}|},
      where \eqn{\Psi_k^{(j)} = \sum_{i\in G_k}\psi_i^{(j)}}.}
    \item{\code{"sum_abs_feat"}}{\eqn{I_k = (1/|G_k|)\sum_{i\in G_k}I_i},
      average of within-group feature importances.}
    \item{\code{"both"}}{Compute both; the \code{importance} column uses
      \code{"sum_abs"} for sorting; an extra column
      \code{importance_sum_abs_feat} is appended.}
  }}

\item{sort}{Logical.  Sort by decreasing importance.  Default \code{TRUE}.}

\item{normalize}{Logical.  Normalize scores to sum to 1.
  Default \code{FALSE}.}
}
\value{
An object of class \code{"treeowen_importance"} (a named list):
\describe{
  \item{\code{feature}}{Data frame: columns \code{feature} and
    \code{importance}, sorted by decreasing importance.  Present when
    \code{type \%in\% c("both","feature")}.}
  \item{\code{group}}{Data frame: columns \code{group} and
    \code{importance} (plus \code{importance_sum_abs_feat} when
    \code{group_agg = "both"}), sorted by decreasing importance.
    Present when \code{type \%in\% c("both","group")}.}
  \item{\code{group_attr}}{Numeric matrix (\eqn{n \times K}):
    group-level Owen attributions \eqn{\Psi_k^{(j)}}.  Always present;
    used internally by plotting functions.}
  \item{\code{group_pos}}{Named list: group name \eqn{\to} column indices
    in \code{result$owens}.}
  \item{\code{group_names}}{Character vector of group names.}
  \item{\code{feat_names}}{Character vector of feature names.}
  \item{\code{n}, \code{K}, \code{p}}{Integer scalars: observations,
    groups, features.}
  \item{\code{type}, \code{group_agg}}{Arguments as supplied.}
}
}
\details{
\subsection{Definitions}{
Feature importance:
\deqn{I_i = \frac{1}{n}\sum_{j=1}^n |\psi_i^{(j)}|}

Group attribution (per observation \eqn{j}):
\deqn{\Psi_k^{(j)} = \sum_{i \in G_k} \psi_i^{(j)}}

Group importance (\code{"sum_abs"}):
\deqn{I_k = \frac{1}{n}\sum_{j=1}^n |\Psi_k^{(j)}|}

Note: \eqn{I_k \ne \sum_{i\in G_k}I_i} in general because signs may
cancel within a group.
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
## ── Synthetic data + XGBoost (see ?treeowen for data setup) ─────────────────
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=5), "G", rep(1:5, each=5))
X <- as.data.frame(matrix(rnorm(50*25), 50, 25, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2
Y      <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:5, function(k) paste0("F",1:5,"G",k)), paste0("G",1:5))

library(xgboost); library(treeshap)
xgb_mod <- xgboost(data    = xgb.DMatrix(as.matrix(X), label=Y),
                   nrounds = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
unified <- xgboost.unify(xgb_mod, X)
result  <- treeowen(unified, X, groups, method="auto",
                    dp_progress=FALSE, verbose=FALSE)

## ── Feature and group importance ────────────────────────────────────────────
imp <- treeowen_importance(result, type = "both", group_agg = "sum_abs")
print(imp)

# Top 5 features
head(imp$feature, 5L)

# Group importance
imp$group

# Normalised (sums to 1)
imp_norm <- treeowen_importance(result, type="group", normalize=TRUE)
stopifnot(abs(sum(imp_norm$group$importance) - 1) < 1e-10)

# group_attr matrix: n x K group-level attributions
dim(imp$group_attr)       # 50 x 5
head(imp$group_attr)
}
}
