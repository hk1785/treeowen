\name{print.treeowen_result}
\alias{print.treeowen_result}
\title{Print Method for treeowen_result Objects}
\description{
Prints a concise summary of a \code{"treeowen_result"} object: dimensions,
computation method, baseline and prediction value ranges, and tree structure.
}
\usage{
\method{print}{treeowen_result}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"treeowen_result"}.}
\item{...}{Additional arguments (currently ignored).}
}
\value{Returns \code{x} invisibly.}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{\code{\link{treeowen}}}
\examples{
\dontrun{
library(treeowen)
library(xgboost)

## Minimal setup (see ?treeowen for a full example)
unified <- xgboost_unify_compat(xgb_model, X)
result  <- treeowen(unified, X, groups,
                    dp_progress = FALSE, verbose = FALSE)
print(result)
}
}
