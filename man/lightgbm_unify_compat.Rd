\name{lightgbm_unify_compat}
\alias{lightgbm_unify_compat}
\title{Unify a LightGBM Model for Use with treeowen}
\description{
A robust replacement for \code{treeshap::lightgbm.unify()} that handles
all known version-compatibility issues across lightgbm 3.x and 4.x and
treeshap 0.3 to 0.4.

The function tries three strategies in order: (1) \code{treeshap::lightgbm.unify()},
(2) direct parsing of \code{lgb.model.dt.tree()} output, and (3) parsing
the raw JSON tree dump from \code{model$dump_model()}.

Handled issues: Booster class restructuring in lightgbm 4.0; column name
changes (\code{split_gain}/\code{value}, \code{threshold}/\code{split_point},
\code{count}/\code{internal_count}); three different child-node ID encodings
(negative leaf index, positive node index, string format); feature name
extraction from multiple fallback sources.
}
\usage{
lightgbm_unify_compat(model, data, recalculate = FALSE)
}
\arguments{
\item{model}{A trained \code{lgb.Booster} object.}
\item{data}{A data frame or matrix of the training features (no outcome
 column).}
\item{recalculate}{Logical. Default \code{FALSE}.}
}
\value{
A \code{model_unified} object suitable for \code{\link{treeowen}}.
}
\seealso{
\code{\link{treeowen}}, \code{\link{xgboost_unify_compat}},
\code{\link{ranger_unify_compat}}
}
\examples{
\dontrun{
library(treeowen)
library(lightgbm)

library(lightgbm)
ds      <- lgb.Dataset(as.matrix(X), label = Y)
model   <- lgb.train(list(objective = "binary"), ds, nrounds = 100,
                     verbose = -1L)
unified <- lightgbm_unify_compat(model, X)
result  <- treeowen(unified, X, groups)
}
}
