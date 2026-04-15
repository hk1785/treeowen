\name{treeowen_exact_enum}
\alias{treeowen_exact_enum}
\title{Model-Agnostic Owen Value Enumerator (Reference)}
\description{
Computes Owen values by exhaustive model-agnostic combinatorial enumeration:
all outer coalitions \eqn{T \subseteq [K]\setminus\{k\}} and all inner
coalitions \eqn{S \subseteq G_k} are evaluated.  Intended as a correctness
reference for small problems; complexity is exponential in \eqn{K} and
\eqn{\max_k|G_k|}.  Use only when \eqn{K \le 5} and \eqn{|G_k| \le 5}.
}
\usage{
treeowen_exact_enum(
  unified_model,
  x,
  groups,
  chunk_size_inner = 131072L,
  max_bytes        = 512 * 1024^2,
  dp_progress      = FALSE,
  use_cpp          = TRUE,
  verbose          = FALSE
)
}
\arguments{
\item{unified_model}{A unified model object from \pkg{treeshap}.}
\item{x}{A data frame or matrix of observations (\eqn{n \times p}).}
\item{groups}{A named list of character vectors partitioning all features.}
\item{chunk_size_inner}{Integer.  Batch size for inner forest evaluation.
  Default \code{131072L}.}
\item{max_bytes}{Numeric.  Memory limit (bytes).
  Default \code{512 * 1024^2}.}
\item{dp_progress}{Logical.  Show a progress bar.  Default \code{FALSE}.}
\item{use_cpp}{Logical.  Use C++ routines.  Default \code{TRUE}.}
\item{verbose}{Logical.  Print diagnostics.  Default \code{FALSE}.}
}
\value{
An object of class \code{"treeowen_result"} (same structure as
\code{\link{treeowen}}): \code{owens}, \code{baseline}, \code{v_all},
\code{observations}, \code{groups}, \code{stats}, \code{note}.
}
\details{
\subsection{Limit}{
\eqn{K} must not exceed \code{TREEOWEN_EXACT_ENUM_K_MAX} (30), which is the
C++ compile-time limit (\code{d_minus < 30}).
}
\subsection{Comparison with treeowen()}{
Uses the same forest evaluation kernel as \code{\link{treeowen}} but without
hierarchy-guided context precomputation.  For small problems, results should
match \code{\link{treeowen}} with \code{method = "exact"} up to
floating-point precision.
}
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{treeowen_exact_tuvalues}}
}
\examples{
\dontrun{
## ── Small problem: K=2, n=5 ─────────────────────────────────────────────────
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=5), "G", rep(1:5, each=5))
X <- as.data.frame(matrix(rnorm(50*25), 50, 25, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3
Y      <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:5, function(k) paste0("F",1:5,"G",k)), paste0("G",1:5))

library(xgboost); library(treeshap)
xgb_mod <- xgboost(data    = xgb.DMatrix(as.matrix(X), label=Y),
                   nrounds = 50L, max_depth = 3L, eta = 0.1,
                   objective = "binary:logistic", verbose = 0L)
unified <- xgboost.unify(xgb_mod, X)

## treeowen_exact_enum requires a complete partition,
## so we pass the full feature matrix but only 2 groups
small_x <- X[1:5, ]
small_g <- groups[c("G1","G2")]   # 2 groups, but x still has all 25 cols

ref  <- treeowen_exact_enum(unified, small_x, small_g, verbose=FALSE)
main <- treeowen(unified, small_x, small_g, method="exact",
                 dp_progress=FALSE, verbose=FALSE)

## Results should match to floating-point precision
max_diff <- max(abs(as.matrix(ref$owens) - as.matrix(main$owens)))
cat(sprintf("Max |enum - treeowen| diff: %.2e\n", max_diff))
stopifnot(max_diff < 1e-6)
}
}
