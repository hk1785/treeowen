\name{treeowen}
\alias{treeowen}
\title{Compute Owen Values for Tree-Based Ensembles}
\description{
Computes exact or Monte Carlo Owen values for every observation in \code{x}
using a trained tree-based ensemble and a feature partition \code{groups}.
Owen values allocate model attribution in two stages: first across feature
groups, then within each group.
}
\usage{
treeowen(
 unified_model,
 x,
 groups,
 method = c("auto", "exact", "approx"),
 hierarchy = NULL,
 n_inner_mc = 64L,
 inner_antithetic = TRUE,
 target_se_inner = 1e-3,
 min_inner_mc = 32L,
 max_inner_mc = 512L,
 chunk_size_inner = 131072L,
 max_bytes = 512 * 1024^2,
 inner_bitmask_max = TREEOWEN_INNER_BITMASK_DEFAULT,
 auto_exact_max_m = TREEOWEN_AUTO_EXACT_MAX_M,
 check_efficiency = FALSE,
 efficiency_tol = 1e-6,
 dp_progress = TRUE,
 dp_print_every = 1L,
 dp_newline_every = 10L,
 eta_alpha = 0.2,
 n_cores = 1L,
 use_cpp = TRUE,
 verbose = TRUE
)
}
\arguments{
\item{unified_model}{A unified model object produced by one of the
 compatibility wrappers: \code{xgboost_unify_compat()} for XGBoost,
 \code{lightgbm_unify_compat()} for LightGBM, or
 \code{ranger_unify_compat()} for Ranger.
 These wrappers handle all known version-compatibility issues across
 different releases of xgboost, lightgbm, ranger, and treeshap,
 and automatically fall back to direct tree parsing when treeshap
 calls fail.}

\item{x}{A data frame or matrix of observations (rows = observations,
 columns = features). Must contain the same features as the model.}

\item{groups}{A named list of character vectors. Each element names the
 features belonging to that group. Together they must form a complete,
 non-overlapping partition of all features in the model.}

\item{method}{Character string selecting the computation algorithm:
 \describe{
 \item{\code{"auto"}}{(default) Uses the exact method for groups
  with fewer features than \code{auto_exact_max_m}, and Monte Carlo
  for larger groups.}
 \item{\code{"exact"}}{Forces the exact method for all groups.
  Requires group size at most 60.}
 \item{\code{"approx"}}{Forces Monte Carlo for all groups.}
 }}

\item{hierarchy}{Controls the auxiliary binary tree over groups used
 internally for computation. \code{NULL} (default) builds a balanced
 binary tree over the K groups. A nested list, numeric vector, or
 integer matrix can specify a custom structure. The choice does not
 affect the correctness of results, only computation speed.
 See \code{\link{build_hierarchy_tree_binary}} and
 \code{\link{build_hierarchy_tree_from_layers}}.}

\item{n_inner_mc}{Integer. Initial number of Monte Carlo permutation
 pairs per group context. Default \code{64L}.}

\item{inner_antithetic}{Logical. Use antithetic sampling to reduce
 Monte Carlo variance. Default \code{TRUE}.}

\item{target_se_inner}{Numeric. Target standard error for adaptive
 Monte Carlo stopping. Default \code{1e-3}.}

\item{min_inner_mc}{Integer. Minimum Monte Carlo pairs per context.
 Default \code{32L}.}

\item{max_inner_mc}{Integer. Maximum Monte Carlo pairs per context.
 Default \code{512L}.}

\item{chunk_size_inner}{Integer. Batch size for inner forest
 evaluation. Default \code{131072L}.}

\item{max_bytes}{Numeric. Memory limit in bytes for batched
 evaluation. Default \code{512 * 1024^2} (512 MB).}

\item{inner_bitmask_max}{Integer. Maximum group size for the
 bitmask-based exact enumerator. Default
 \code{TREEOWEN_INNER_BITMASK_DEFAULT} (60).}

\item{auto_exact_max_m}{Integer. Groups with this many features or
 more use Monte Carlo when \code{method = "auto"}. Default
 \code{TREEOWEN_AUTO_EXACT_MAX_M} (30).}

\item{check_efficiency}{Logical. After computation, verify that Owen
 values sum to the prediction minus the average prediction for each
 observation. Warns if the deviation exceeds \code{efficiency_tol}.
 Default \code{FALSE}.}

\item{efficiency_tol}{Numeric. Tolerance for the efficiency check.
 Default \code{1e-6}.}

\item{dp_progress}{Logical. Show a progress bar. Default \code{TRUE}.}

\item{dp_print_every}{Integer. Progress bar update frequency in rows.
 Default \code{1L}.}

\item{dp_newline_every}{Integer. Deprecated; retained for backwards
 compatibility. Default \code{10L}.}

\item{eta_alpha}{Numeric. EMA smoothing factor for ETA estimation.
 Default \code{0.2}.}

\item{n_cores}{Integer. Number of parallel cores. Default \code{1L}.
 Values greater than 1 use \code{parallel::mclapply()} on Linux and
 macOS. On Windows, keep this at 1.}

\item{use_cpp}{Logical. Use the compiled C++ backend. Must be
 \code{TRUE}; a pure-R fallback is not available. Default \code{TRUE}.}

\item{verbose}{Logical. Print diagnostic messages. Default \code{TRUE}.}
}

\details{
\subsection{What are Owen values?}{
Owen values extend standard Shapley values to settings where features
are organized into groups. Attribution is allocated in two stages.
In the first stage, each group receives a score reflecting how much
the group as a whole shifts the prediction away from the average
prediction. In the second stage, that group-level score is divided
among the individual features in the group using the standard Shapley
formula. The final Owen values satisfy the efficiency property: they
sum to the prediction minus the average prediction across all
observations.
}
\subsection{Algorithms}{
\describe{
 \item{Exact (\code{"exact"})}{Enumerates all subsets of each group
  explicitly. Feasible for groups with up to 30 features using R
  integer arithmetic, or up to 60 features using the C++ path.}
 \item{Monte Carlo (\code{"approx"})}{Uses antithetic permutation
  pairs to estimate inner Shapley values. Sampling is controlled
  adaptively by \code{target_se_inner}, \code{min_inner_mc}, and
  \code{max_inner_mc}.}
 \item{Auto (\code{"auto"})}{Applies exact to groups with fewer
  features than \code{auto_exact_max_m} and Monte Carlo to the rest.}
}
}
\subsection{Efficiency axiom}{
For the exact algorithm, Owen values sum exactly to the prediction
minus the average, up to floating-point precision. Set
\code{check_efficiency = TRUE} to monitor this property at runtime.
For Monte Carlo results the deviation is controlled by
\code{target_se_inner}.
}
\subsection{Parallelism}{
Uses \code{parallel::mclapply()} (fork-based) on Linux and macOS.
Each worker independently rebuilds the C++ internal representation.
On Windows, use \code{n_cores = 1L}.
}
}

\value{
An object of class \code{"treeowen_result"} (a named list):
\describe{
 \item{\code{owens}}{Data frame of Owen values (rows = observations,
  columns = features). Column names match the features in \code{x}.}
 \item{\code{observations}}{Data frame of the input observations.}
 \item{\code{groups}}{The \code{groups} argument as supplied.}
 \item{\code{baseline}}{Numeric vector: the average model prediction
  used as the reference baseline.}
 \item{\code{v_all}}{Numeric vector: the model prediction for each
  observation.}
 \item{\code{stats}}{Named list of diagnostics including the number
  of observations, features, groups, trees, effective method, group
  methods, C++ usage, and parallelism settings.}
 \item{\code{note}}{Character string summarising the run.}
}
}

\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen_importance}}, \code{\link{treeowen_beeswarm}},
\code{\link{treeowen_hierarchical_beeswarm}},
\code{\link{build_hierarchy_tree_from_layers}},
\code{\link{print.treeowen_result}}
}
\examples{
\dontrun{
library(treeowen)
library(xgboost)

## Synthetic data: 5 groups of 5 features each, 50 observations
set.seed(42)
feat_names <- paste0("F", rep(seq_len(5), times = 20),
                     "G", rep(seq_len(20), each  = 5))
X <- as.data.frame(matrix(rnorm(100 * 100), 100, 100,
                           dimnames = list(NULL, feat_names)))
Y <- as.integer(0.8 * X$F1G1 - 0.6 * X$F2G1 + 0.5 * X$F1G3 > 0)
groups <- setNames(
  lapply(seq_len(20), function(k) paste0("F", seq_len(5), "G", k)),
  paste0("G", seq_len(20)))

## Fit XGBoost and compute Owen values
library(xgboost)
xgb_mat <- xgb.DMatrix(as.matrix(X), label = Y)
xgb_mod <- xgboost(data = xgb_mat, nrounds = 50L, max_depth = 3L,
                   eta = 0.1, objective = "binary:logistic", verbose = 0L)
unified <- xgboost_unify_compat(xgb_mod, X)
result  <- treeowen(unified, X, groups,
                    method = "auto", dp_progress = FALSE, verbose = TRUE)
print(result)

## Check efficiency axiom
eff_err <- abs(rowSums(result$owens) - (result$v_all - result$baseline))
cat(sprintf("Max efficiency error: %.2e\n", max(eff_err)))
}
}
