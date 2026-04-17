\name{ranger_unify_compat}
\alias{ranger_unify_compat}
\title{Unify a Ranger Model for Use with treeowen}
\description{
A robust replacement for \code{treeshap::ranger.unify()} that handles
all known version-compatibility issues across ranger 0.11 to 0.16+ and
treeshap 0.3 to 0.4.

The function tries three strategies in order: (1) \code{treeshap::ranger.unify()},
(2) patched \code{treeInfo()} output with \code{treeshap:::ranger_unify.common()},
and (3) direct parsing of \code{model$forest} internals.

Handled issues: \code{treeInfo()} returning \code{pred.1} (older ranger)
versus \code{prediction} (ranger 0.14+); column name variations
(\code{nodeID}/\code{node}, \code{splitvarID}/\code{splitVarID},
\code{splitval}/\code{splitVal}); \code{child.nodeIDs} stored as a
list of length-2 vectors versus separate left/right child columns;
leaf predictions as scalar, matrix, or list of vectors.

The model must be trained with \code{probability = TRUE}.
}
\usage{
ranger_unify_compat(model, data, recalculate = FALSE)
}
\arguments{
\item{model}{A trained \code{ranger} object. Must have been trained
 with \code{probability = TRUE} and \code{write.forest = TRUE}.}
\item{data}{A data frame or matrix of the training features (no outcome
 column).}
\item{recalculate}{Logical. Default \code{FALSE}.}
}
\value{
A \code{model_unified} object suitable for \code{\link{treeowen}}.
}
\seealso{
\code{\link{treeowen}}, \code{\link{xgboost_unify_compat}},
\code{\link{lightgbm_unify_compat}}
}
\examples{
\dontrun{
library(treeowen)
library(ranger)

## Train Ranger (probability = TRUE required) and unify for treeowen()
df      <- cbind(X, .y = factor(Y))
model   <- ranger(.y ~ ., data = df, num.trees = 100,
                  probability = TRUE, keep.inbag = TRUE)
unified <- ranger_unify_compat(model, X)
result  <- treeowen(unified, X, groups)
}
}
