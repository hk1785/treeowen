\name{TREEOWEN_AUTO_EXACT_MAX_M}
\alias{TREEOWEN_AUTO_EXACT_MAX_M}
\title{Exact/Monte Carlo Switching Threshold}
\description{
An exported integer constant equal to \code{15L}. When
\code{\link{treeowen}} is called with \code{method = "auto"}, groups with
the group size <= \code{TREEOWEN_AUTO_EXACT_MAX_M} use exact enumeration;
larger groups use Monte Carlo.
}
\format{An integer scalar (\code{15L}).}
\details{
Under the default threshold of 15, groups with the group size <= 15 are
evaluated using the exact algorithm, whereas groups with the group size > 15
are evaluated using the Monte Carlo algorithm.

The exact path itself is not limited to 15: the R-integer bitmask path
operates on 30-bit integers and the C++ bitmask path supports up to
\code{\link{TREEOWEN_INNER_BITMASK_DEFAULT}} (60), so \code{method = "exact"}
can still be forced for 15 < the group size <= 60.

Override per call:
\preformatted{result <- treeowen(unified, X, groups,
 method = "auto",
 auto_exact_max_m = 15L)}
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{TREEOWEN_INNER_BITMASK_DEFAULT}}
}
\examples{
TREEOWEN_AUTO_EXACT_MAX_M # 15L
}
