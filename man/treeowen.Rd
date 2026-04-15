\name{treeowen}
\alias{treeowen}
\title{Compute Owen Values for Tree-Based Ensembles}
\description{
Computes exact or Monte Carlo Owen values for every observation in \code{x}
using a trained tree-based ensemble and a feature partition \code{groups}.
The Owen value allocates attribution in two stages: an outer Shapley game
across groups and an inner Shapley game within each group.
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
 compatibility wrappers: \code{.xgboost_unify_compat()} for XGBoost,
 \code{.lightgbm_unify_compat()} for LightGBM, or
 \code{.ranger_unify_compat()} for Ranger.
 These wrappers handle all known version-compatibility issues across
 different releases of each library and automatically fall back to
 direct tree parsing when treeshap calls fail. Use them instead of
 calling \code{treeshap::xgboost.unify()},
 \code{treeshap::lightgbm.unify()}, or
 \code{treeshap::ranger.unify()} directly.}

\item{x}{A data frame or matrix of observations (n x p).
 Must contain the same features as \code{unified_model}.}

\item{groups}{A named list of character vectors. Each element names the
 features belonging to that group. Together they must form a complete,
 non-overlapping partition of all features in \code{unified_model}.}

\item{method}{Character. Algorithm selection:
 \describe{
 \item{\code{"auto"}}{(default) Exact for groups with
 the group size < \code{auto_exact_max_m}; Monte Carlo otherwise.}
 \item{\code{"exact"}}{Force exact for all groups
 (requires the group size <= 60).}
 \item{\code{"approx"}}{Force Monte Carlo for all groups.}
 }}

\item{hierarchy}{Controls the auxiliary binary tree mathcal{B} over
 groups. \code{NULL} (default) builds a balanced binary tree safe for any
 K. A nested list (see \code{\link{build_hierarchy_tree_binary}}),
 a numeric vector, or a K-column integer matrix (see
 \code{\link{build_hierarchy_tree_from_layers}}) are also accepted.
 The choice does not affect correctness, only runtime.}

\item{n_inner_mc}{Integer. Initial Monte Carlo permutation pairs per
 context when \code{method = "approx"} or \code{"auto"}.
 Default \code{64L}.}

\item{inner_antithetic}{Logical. Antithetic sampling to reduce Monte Carlo
 variance. Default \code{TRUE}.}

\item{target_se_inner}{Numeric. Target standard error for adaptive Monte
 Carlo stopping. Default \code{1e-3}.}

\item{min_inner_mc}{Integer. Minimum Monte Carlo pairs per context.
 Default \code{32L}.}

\item{max_inner_mc}{Integer. Maximum Monte Carlo pairs per context.
 Default \code{512L}.}

\item{chunk_size_inner}{Integer. Batch size for inner forest evaluation.
 Default \code{131072L}.}

\item{max_bytes}{Numeric. Memory limit (bytes) for batched evaluation.
 Default \code{512 * 1024^2} (512 MB).}

\item{inner_bitmask_max}{Integer. Maximum the group size for bitmask-based
 exact enumeration. Default \code{TREEOWEN_INNER_BITMASK_DEFAULT} (60).}

\item{auto_exact_max_m}{Integer. Groups with the group size >=
 \code{auto_exact_max_m} use Monte Carlo when \code{method = "auto"}.
 Default \code{TREEOWEN_AUTO_EXACT_MAX_M} (30).}

\item{check_efficiency}{Logical. Verify the efficiency axiom
 sum_i psi_i = f(x) - E[f(X)] after computation. Warns if
 violated beyond \code{efficiency_tol}. Default \code{FALSE}.}

\item{efficiency_tol}{Numeric. Tolerance for the efficiency check.
 Default \code{1e-6}.}

\item{dp_progress}{Logical. Show a progress bar. Default \code{TRUE}.}

\item{dp_print_every}{Integer. Progress bar update frequency (rows).
 Default \code{1L}.}

\item{dp_newline_every}{Integer. Deprecated; kept for backwards
 compatibility. Default \code{10L}.}

\item{eta_alpha}{Numeric. EMA smoothing factor for ETA estimation.
 Default \code{0.2}.}

\item{n_cores}{Integer. Parallel cores. Default \code{1L}. Values
 \code{> 1} use \code{parallel::mclapply()} on Linux/macOS.
 Windows requires \code{n_cores = 1L}.}

\item{use_cpp}{Logical. Use compiled C++ routines (strongly recommended).
 Default \code{TRUE}.}

\item{verbose}{Logical. Print diagnostic messages. Default \code{TRUE}.}
}
\value{
An object of class \code{"treeowen_result"} (a named list):
\describe{
 \item{\code{owens}}{Data frame (n x p) of Owen values.
 Column names match the features in \code{x}.}
 \item{\code{observations}}{Data frame (n x p) of input
 observations.}
 \item{\code{groups}}{The \code{groups} argument as supplied.}
 \item{\code{baseline}}{Numeric vector of length n: the empty-set
 forest value v_x(varnothing) = E[f(X)] per observation.}
 \item{\code{v_all}}{Numeric vector of length n: the full-set value
 v_x(N) = f(x) per observation.}
 \item{\code{stats}}{Named list of diagnostics: \code{n}, \code{p},
 \code{K}, tree counts, effective method, group methods, C++ usage,
 and parallelism settings.}
 \item{\code{note}}{Character string summarising the run configuration.}
}
}
\details{
\subsection{Owen value definition}{
Let N = {1,ldots,p} be partitioned into K non-overlapping
groups G_1,ldots,G_K. The Owen value of feature i in G_k
is
} w(|T|,K)
 \sum_{S \subseteq G_k\setminus\{i\}} w(|S|,|G_k|)
 \bigl[v_x(G_T \cup S \cup \{i\}) - v_x(G_T \cup S)\bigr]
}
where G_T = bigcup_{k'in TG_{k'}}, v_x(cdot) is the
tree ensemble's conditional expectation, and
w(s,m) = s!(m-s-1)!/m! is the Shapley weight.
}
\subsection{Algorithms}{
\describe{
 \item{Exact (\code{"exact"})}{Bitmask enumeration of all
 2^{the group size} inner subsets. Feasible for the group size <= 30
 (R integer path) or up to 60 (C++ path).}
 \item{Monte Carlo (\code{"approx"})}{Antithetic permutation pairs.
 Adaptive stopping via \code{target_se_inner}, \code{min_inner_mc},
 \code{max_inner_mc}.}
 \item{Auto (\code{"auto"})}{Per-group selection: exact if
 the group size < \code{auto_exact_max_m}, Monte Carlo otherwise.}
}
}
\subsection{Efficiency axiom}{
For the exact algorithm, sum_ipsi_i(x)=f(x)-E[f(X)] holds up to
floating-point precision. Set \code{check_efficiency = TRUE} to monitor.
}
\subsection{Parallelism}{
Uses \code{parallel::mclapply()} (fork-based) on Linux/macOS. Each
worker independently rebuilds the C++ XPtr. Windows: \code{n_cores = 1L}
only.
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
## ── Synthetic data (5 groups × 5 features, n = 50) ─────────────────────────
set.seed(42)
n_obs <- 50L; n_groups <- 5L; n_feat_g <- 5L
feat_names <- paste0("F", rep(seq_len(n_feat_g), times = n_groups),
 "G", rep(seq_len(n_groups), each = n_feat_g))
X <- as.data.frame(matrix(rnorm(n_obs * 25L), n_obs, 25L,
 dimnames = list(NULL, feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3 + 0.3*X$F1G2
Y <- as.integer(log_odds > 0)
groups <- setNames(lapply(seq_len(n_groups),
 function(k) paste0("F", seq_len(n_feat_g), "G", k)),
 paste0("G", seq_len(n_groups)))

## ── XGBoost model ───────────────────────────────────────────────────────────
library(xgboost); library(treeshap)
xgb_mat <- xgb.DMatrix(as.matrix(X), label = Y)
xgb_mod <- xgboost(data = xgb_mat, nrounds = 50L, max_depth = 3L,
 eta = 0.1, objective = "binary:logistic", verbose = 0L)
unified <- .xgboost_unify_compat(xgb_mod, X)

## ── Owen values ─────────────────────────────────────────────────────────────
result <- treeowen(unified, X, groups,
 method = "auto", dp_progress = FALSE, verbose = TRUE)
print(result)
dim(result$owens) # 50 x 25
head(result$owens[, 1:5]) # first 5 Owen values per observation

## ── Efficiency axiom (should be < 1e-6 for exact) ──────────────────────────
eff_err <- abs(rowSums(result$owens) - (result$v_all - result$baseline))
cat(sprintf("Max efficiency error: %.2e\n", max(eff_err)))
}
}
