\name{.xgboost_unify_compat}
\alias{.xgboost_unify_compat}
\title{Unify an XGBoost Model for Use with treeowen}
\description{
A robust replacement for \code{treeshap::xgboost.unify()} that handles
all known version-compatibility issues across xgboost 1.6 to 2.x and
treeshap 0.3 to 0.4. Builds the unified model object directly from
\code{xgb.model.dt.tree()}, bypassing treeshap internals entirely.

Handled issues: the \code{ID} column absent in older xgboost; the
split-gain column named \code{Quality} instead of \code{Gain} in some
versions; \code{xgb.model.dt.tree()} argument renamed from \code{model}
to \code{xgb_model} in xgboost 2.x; \code{Yes}/\code{No}/\code{Missing}
stored as node-ID strings instead of row indices; \code{model$feature_names}
being NULL for models trained from unnamed matrices.
}
\usage{
.xgboost_unify_compat(model, data, recalculate = FALSE)
}
\arguments{
\item{model}{A trained \code{xgb.Booster} object.}
\item{data}{A data frame or matrix of the training features (no outcome
 column). Used to resolve feature names when \code{model$feature_names}
 is NULL.}
\item{recalculate}{Logical. If \code{TRUE} and treeshap is available,
 calls \code{treeshap::set_reference_dataset()} after building the
 unified object. Default \code{FALSE}.}
}
\value{
A \code{model_unified} object suitable for \code{\link{treeowen}}.
}
\seealso{
\code{\link{treeowen}}, \code{\link{.lightgbm_unify_compat}},
\code{\link{.ranger_unify_compat}}
}
\examples{
\dontrun{
library(xgboost)
dm      <- xgb.DMatrix(as.matrix(X), label = Y)
model   <- xgboost(dm, nrounds = 100, objective = "binary:logistic",
                   verbose = 0)
unified <- .xgboost_unify_compat(model, X)
result  <- treeowen(unified, X, groups)
}
}
