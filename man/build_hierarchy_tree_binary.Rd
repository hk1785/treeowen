\name{build_hierarchy_tree_binary}
\alias{build_hierarchy_tree_binary}
\title{Build a Custom Binary Auxiliary Tree}
\description{
Constructs the auxiliary binary tree mathcal{B} used by
\code{\link{treeowen}} to organise the outer Shapley-weighted summation,
from a user-supplied nested-list specification.
}
\usage{
build_hierarchy_tree_binary(hier)
}
\arguments{
\item{hier}{A nested list representing the binary tree. Each leaf must be a
 length-1 integer (group index). Each internal node must have exactly
 two children as named elements \code{left}/\code{right} or positional
 elements \code{[[1]]}/\code{[[2]]}.

 Example for 4 groups \code{((1,2),(3,4))}:
 \preformatted{list(left = list(left=1L, right=2L),
 right = list(left=3L, right=4L))}
}
}
\value{
A normalised nested list representing the binary tree. Each node
contains:
\describe{
 \item{\code{type}}{Character: \code{"leaf"} or \code{"node"}.}
 \item{\code{uid}}{Integer unique identifier.}
 \item{\code{id}}{Integer group index (leaves only; \code{NA} for
 internal nodes).}
 \item{\code{children}}{List of child nodes (internal nodes only).}
 \item{\code{leaves}}{Integer vector of all group indices in this
 subtree.}
}
Pass this object as the \code{hierarchy} argument of \code{\link{treeowen}}.
}
\details{
The choice of hierarchy does not affect correctness; it only influences
runtime. A balanced binary tree (\code{hierarchy = NULL} in
\code{\link{treeowen}}) is the default and is suitable for most use cases.
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{build_hierarchy_tree_from_layers}}
}
\examples{
\dontrun{
library(treeowen)
library(xgboost)

## Custom hierarchy for 5 groups: ((G1,G2), (G3,(G4,G5)))
hier <- list(
 left = list(left = 1L, right = 2L),
 right = list(left = 3L,
 right = list(left = 4L, right = 5L))
)
tree <- build_hierarchy_tree_binary(hier)
str(tree, max.level = 2) # nested list structure

## Use with treeowen()
set.seed(42)
feat_names <- paste0("F", rep(1:5, times=20), "G", rep(1:20, each=5))
X <- as.data.frame(matrix(rnorm(100*100), 100, 100, dimnames=list(NULL,feat_names)))
log_odds <- 0.8*X$F1G1 - 0.6*X$F2G1 + 0.5*X$F1G3
Y <- as.integer(log_odds > 0)
groups <- setNames(lapply(1:20, function(k) paste0("F",1:5,"G",k)), paste0("G",1:20))

library(xgboost); library(treeshap)
xgb_mod <- xgboost(data = xgb.DMatrix(as.matrix(X), label=Y),
 nrounds = 50L, max_depth = 3L, eta = 0.1,
 objective = "binary:logistic", verbose = 0L)
unified <- xgboost_unify_compat(xgb_mod, X)

result_custom <- treeowen(unified, X, groups, method="exact",
 hierarchy = tree,
 dp_progress=FALSE, verbose=FALSE)

## Hierarchy-invariance: balanced vs custom hierarchy give same Owen values
result_default <- treeowen(unified, X, groups, method="exact",
 dp_progress=FALSE, verbose=FALSE)
max_diff <- max(abs(as.matrix(result_custom$owens) -
 as.matrix(result_default$owens)))
cat(sprintf("Max |custom - default| diff: %.2e (should be ~0)\n", max_diff))
stopifnot(max_diff < 1e-8)
}
}
