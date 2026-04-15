\name{print.treeowen_plots}
\alias{print.treeowen_plots}
\title{Print Method for treeowen_plots Objects}
\description{
Renders a \code{"treeowen_plots"} object produced by
\code{\link{treeowen_beeswarm}} with \code{level = "both"}. Prints the
combined patchwork plot when available; otherwise prints individual panels.
}
\usage{
\method{print}{treeowen_plots}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"treeowen_plots"}.}
\item{...}{Additional arguments (currently ignored).}
}
\value{Returns \code{x} invisibly.}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{\code{\link{treeowen_beeswarm}}}
\examples{
\dontrun{
library(treeowen)
library(xgboost)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)

out <- treeowen_beeswarm(result, level = "both")
print(out) # combined patchwork
print(out$feature) # feature panel only
}
}
