\name{print.treeowen_importance}
\alias{print.treeowen_importance}
\title{Print Method for treeowen_importance Objects}
\description{
Prints a formatted summary of a \code{"treeowen_importance"} object:
top-10 feature importances and all group importances.
}
\usage{
\method{print}{treeowen_importance}(x, digits = 4, ...)
}
\arguments{
\item{x}{An object of class \code{"treeowen_importance"}.}
\item{digits}{Integer.  Significant digits.  Default \code{4}.}
\item{...}{Additional arguments passed to \code{print.data.frame}.}
}
\value{Returns \code{x} invisibly.}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{\code{\link{treeowen_importance}}}
\examples{
\dontrun{
imp <- treeowen_importance(result, type = "both")
print(imp)
print(imp, digits = 6)
}
}
