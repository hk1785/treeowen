\name{TREEOWEN_AUTO_EXACT_MAX_M}
\alias{TREEOWEN_AUTO_EXACT_MAX_M}
\title{Exact/Monte Carlo Switching Threshold}
\description{
An exported integer constant equal to \code{30L}.  When
\code{\link{treeowen}} is called with \code{method = "auto"}, groups with
\eqn{|G_k| \ge} \code{TREEOWEN_AUTO_EXACT_MAX_M} use Monte Carlo; smaller
groups use exact enumeration.
}
\format{An integer scalar (\code{30L}).}
\details{
The value 30 reflects the R-integer bitmask exact path, which operates on
30-bit integers.  The C++ bitmask path supports up to
\code{\link{TREEOWEN_INNER_BITMASK_DEFAULT}} (60), so \code{method = "exact"}
can still be forced for \eqn{30 \le |G_k| \le 60}.

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
TREEOWEN_AUTO_EXACT_MAX_M   # 30L
}
