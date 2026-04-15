\name{clear_inner_enum_cache}
\alias{clear_inner_enum_cache}
\title{Clear the Inner Shapley Weight Cache}
\description{
Removes all memoised inner Shapley weight tables from the package-internal
cache environment \code{.INNER_ENUM_CACHE}.  Call this between large
analyses with different group sizes to free memory.
}
\usage{
clear_inner_enum_cache()
}
\value{
Returns \code{NULL} invisibly.
}
\details{
The exact inner Shapley enumerator pre-computes and caches per-group-size
weight tables (popcount lookup, mask-without-bit tables) the first time it
encounters each group size \eqn{|G_k|}.  These tables are reused across
observations and across calls within the same R session.

Memory usage scales with \eqn{2^m} for each cached size \eqn{m \le 30}.
For \eqn{m = 30} the tables occupy approximately 32 MB.
}
\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
\code{\link{treeowen}}, \code{\link{TREEOWEN_AUTO_EXACT_MAX_M}}
}
\examples{
# Clear cache at any time to free memory
clear_inner_enum_cache()

\dontrun{
## Verify cache is empty after clearing
pkg_ns    <- getNamespace("treeowen")
cache_env <- pkg_ns$.INNER_ENUM_CACHE
n_before  <- length(ls(cache_env))
clear_inner_enum_cache()
n_after   <- length(ls(cache_env))
stopifnot(n_after == 0L)
cat(sprintf("Cache cleared: %d -> %d entries\n", n_before, n_after))
}
}
