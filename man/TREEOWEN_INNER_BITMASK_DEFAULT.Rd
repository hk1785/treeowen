\name{TREEOWEN_INNER_BITMASK_DEFAULT}
\alias{TREEOWEN_INNER_BITMASK_DEFAULT}
\title{Bitmask Enumerator Size Limit}
\description{
An exported integer constant equal to \code{60L}.  It represents the
maximum group size \eqn{|G_k|} supported by the C++ bitmask-based exact
inner enumerator.  Groups larger than this value require Monte Carlo.
}
\format{An integer scalar (\code{60L}).}
\details{
The C++ enumerator uses 64-bit integer bitmasks to enumerate all
\eqn{2^{|G_k|}} subsets.  The hard limit is 63 bits; \code{60} provides a
conservative safety margin.

For \eqn{|G_k| > 60}, \code{method = "exact"} raises an error.  Use
\code{method = "approx"} or \code{method = "auto"} for large groups.
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{TREEOWEN_AUTO_EXACT_MAX_M}}
}
\examples{
TREEOWEN_INNER_BITMASK_DEFAULT   # 60L
}
