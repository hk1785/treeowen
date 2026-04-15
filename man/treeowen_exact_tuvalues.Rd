\name{treeowen_exact_tuvalues}
\alias{treeowen_exact_tuvalues}
\title{Tree-Aware Owen Value Enumerator (Reference)}
\description{
Computes Owen values by exhaustive enumeration using tree-aware dynamic
programming, implemented directly from the Owen value formula without the
hierarchy-guided context precomputation used by \code{\link{treeowen}}.
Intended as a correctness reference for small problems; complexity is
sum_k the group size cdot 2^{K+the group size-1}.
}
\usage{
treeowen_exact_tuvalues(
 unified_model,
 x,
 groups,
 chunk_size_inner = 131072L,
 max_bytes = 512 * 1024^2,
 inner_bitmask_max = TREEOWEN_INNER_BITMASK_DEFAULT,
 dp_progress = FALSE,
 use_cpp = TRUE,
 verbose = FALSE
)
}
\arguments{
\item{unified_model}{A unified model object from \pkg{treeshap}.}
\item{x}{A data frame or matrix of observations (n x p).}
\item{groups}{A named list of character vectors partitioning all features.}
\item{chunk_size_inner}{Integer. Batch size for inner forest evaluation.
 Default \code{131072L}.}
\item{max_bytes}{Numeric. Memory limit (bytes).
 Default \code{512 * 1024^2}.}
\item{inner_bitmask_max}{Integer. Maximum the group size for bitmask exact
 enumeration. Default \code{TREEOWEN_INNER_BITMASK_DEFAULT} (60).}
\item{dp_progress}{Logical. Show a progress bar. Default \code{FALSE}.}
\item{use_cpp}{Logical. Use C++ routines. Default \code{TRUE}.}
\item{verbose}{Logical. Print diagnostics. Default \code{FALSE}.}
}
\value{
An object of class \code{"treeowen_result"} (same structure as
\code{\link{treeowen}}): \code{owens}, \code{baseline}, \code{v_all},
\code{observations}, \code{groups}, \code{stats}, \code{note}.
}
\details{
\subsection{Algorithm}{
For each observation, group k, outer coalition
T subseteq [K]setminus{k}, feature i in G_k, and inner
coalition S subseteq G_ksetminus{i}, the function evaluates
v_x(G_T cup S) and v_x(G_T cup S cup {i)} via
\code{eval_v_knownset_cpp}.
}
\subsection{Comparison with treeowen()}{
Results match \code{\link{treeowen}} with \code{method = "exact"} to
floating-point precision. Unlike \code{\link{treeowen_exact_enum}}, there
is no K <= 30 limit, but the function scales poorly with large
K and group sizes.
}
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{treeowen_exact_enum}}
}
\examples{
\dontrun{
## ── Small problem: K=2, n=5 ─────────────────────────────────────────────────
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=5), "G", rep(1:5, each=5))
X <- as.data.frame(matrix(rnorm(50*25), 50, 25, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3
Y <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:5, function(k) paste0("F",1:5,"G",k)), paste0("G",1:5))

library(xgboost); library(treeshap)
xgb_mod <- xgboost(data = xgb.DMatrix(as.matrix(X), label=Y),
 nrounds = 50L, max_depth = 3L, eta = 0.1,
 objective = "binary:logistic", verbose = 0L)
unified <- .xgboost_unify_compat(xgb_mod, X)

small_x <- X[1:5, ]
small_g <- groups[c("G1","G2")]

ref <- treeowen_exact_tuvalues(unified, small_x, small_g, verbose=FALSE)
main <- treeowen(unified, small_x, small_g, method="exact",
 dp_progress=FALSE, verbose=FALSE)

max_diff <- max(abs(as.matrix(ref$owens) - as.matrix(main$owens)))
cat(sprintf("Max |tuvalues - treeowen| diff: %.2e\n", max_diff))
stopifnot(max_diff < 1e-6)
}
}
