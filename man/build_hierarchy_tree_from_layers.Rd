\name{build_hierarchy_tree_from_layers}
\alias{build_hierarchy_tree_from_layers}
\title{Build Auxiliary Tree from a Layer Specification}
\description{
Constructs the auxiliary binary tree mathcal{B} from a layer-based
merge specification, a balanced binary default, or a fully specified nested
list. This is the entry point called internally by \code{\link{treeowen}}
to parse the \code{hierarchy} argument.
}
\usage{
build_hierarchy_tree_from_layers(K, hierarchy)
}
\arguments{
\item{K}{Integer. Number of groups.}

\item{hierarchy}{Controls the tree structure:
 \describe{
 \item{\code{NULL}}{Balanced binary tree over all K groups
 (default in \code{\link{treeowen}}).}
 \item{Named nested list}{Compatible with the output of
 \code{\link{build_hierarchy_tree_binary}}.}
 \item{Numeric vector of length K}{Treated as a single-row
 matrix; groups sharing the same label are merged.}
 \item{L x K integer matrix}{Layer-based spec with L
 merge layers. Each row merges groups that share the same integer
 label in that row.}
 }
}
}
\value{
A normalised nested list (same format as
\code{\link{build_hierarchy_tree_binary}}) with fields \code{type},
\code{uid}, \code{id}, \code{children}, and \code{leaves}. Pass this
directly as the \code{hierarchy} argument of \code{\link{treeowen}}.
}
\details{
\subsection{Layer-matrix format}{
A layer matrix has K columns and L rows. Example with
K=6:
\preformatted{
hierarchy <- matrix(
 c(1,1,2,2,3,3, # row 1: merge {G1,G2}, {G3,G4}, {G5,G6}
 1,1,1,1,2,2), # row 2: merge {G1..G4} vs {G5,G6}
 nrow=2, byrow=TRUE)
}
}
\subsection{Context overflow guard}{
With \code{hierarchy = NULL} (balanced binary tree) each internal node has
at most 2 children, so context depth is bounded and overflow is impossible.
Custom hierarchies with > 29 sibling leaves at any node will raise an
error; use a deeper hierarchy in that case.
}
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{build_hierarchy_tree_binary}}
}
\examples{
\dontrun{
## Balanced binary tree (default)
tree1 <- build_hierarchy_tree_from_layers(K = 5L, hierarchy = NULL)
str(tree1, max.level = 2)

## Layer-matrix: 2-level hierarchy for 6 groups
hier_mat <- matrix(c(1,1,2,2,3,3,
 1,1,1,1,2,2), nrow=2, byrow=TRUE)
tree2 <- build_hierarchy_tree_from_layers(K = 6L, hierarchy = hier_mat)

## Verify hierarchy-invariance of Owen values
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

## Layer-matrix hierarchy for 5 groups
hier5 <- matrix(c(1,1,2,2,3), nrow=1) # {G1,G2}, {G3,G4}, {G5}
tree5 <- build_hierarchy_tree_from_layers(K=5L, hierarchy=hier5)
r_layer <- treeowen(unified, X, groups, method="exact",
 hierarchy=tree5, dp_progress=FALSE, verbose=FALSE)
r_default <- treeowen(unified, X, groups, method="exact",
 dp_progress=FALSE, verbose=FALSE)
max_diff <- max(abs(as.matrix(r_layer$owens) - as.matrix(r_default$owens)))
cat(sprintf("Max |layer - default| diff: %.2e\n", max_diff))
stopifnot(max_diff < 1e-8)
}
}
