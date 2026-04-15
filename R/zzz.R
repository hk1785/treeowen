# treeowen/R/zzz.R
# Package-level constants, mutable state, and .onLoad hook.

# ── Exported constants ────────────────────────────────────────────────────────

#' Auto-Switching Threshold for Exact vs Monte Carlo
#'
#' @name TREEOWEN_AUTO_EXACT_MAX_M
#' @title Auto-Switching Threshold for Exact vs Monte Carlo
#'
#' Default threshold used by \code{\link{treeowen}} when
#' \code{method = "auto"}: groups with \eqn{|G_k| \geq}
#' \code{TREEOWEN_AUTO_EXACT_MAX_M} use Monte Carlo; smaller groups use exact.
#'
#' @format An integer scalar (\code{30L}).
#' @export
TREEOWEN_AUTO_EXACT_MAX_M <- 30L

#' Bitmask Limit for the Exact Inner Enumerator
#'
#' @name TREEOWEN_INNER_BITMASK_DEFAULT
#' @title Bitmask Limit for the Exact Inner Enumerator
#'
#' Maximum group size (\eqn{|G_k|}) supported by the bitmask-based exact
#' inner enumerator. Groups larger than this value must use Monte Carlo.
#'
#' @format An integer scalar (\code{60L}).
#' @export
TREEOWEN_INNER_BITMASK_DEFAULT <- 60L

# ── Internal limits (not exported) ───────────────────────────────────────────
.TREEOWEN_K_HARD_MAX <- 2000L
.TREEOWEN_M_HARD_MAX <- 500L

TREEOWEN_K_DEFAULT        <- 1000L
TREEOWEN_M_DEFAULT        <- 500L
TREEOWEN_EXACT_ENUM_K_MAX <- 30L   # C++ d_minus < 30 제한

# ── Package-environment mutable state ────────────────────────────────────────
.CPP_LOADED    <- FALSE
.CPP_AVAILABLE <- FALSE

# ── Inner Shapley weight cache ────────────────────────────────────────────────
.INNER_ENUM_CACHE <- new.env(hash = TRUE, parent = emptyenv())

# Bitset block width — must equal the value compiled into treeowen_cpp.cpp.
.BITS_PER_BLOCK <- 30L

# ── .onLoad ───────────────────────────────────────────────────────────────────
# useDynLib + .registration = TRUE registers all [[Rcpp::export]] routines.
# We probe a lightweight sentinel to confirm the shared library loaded.
# exists() is NOT reliable inside .onLoad (namespace not yet sealed);
# .Call() on a registered name is the correct test.
.onLoad <- function(libname, pkgname) {
  ok <- tryCatch({
    # R_useDynamicSymbols(FALSE) is set, so we must use the exact registered
    # name from CallEntries: "_treeowen_treeowen_ping".
    # The plain "treeowen_ping" is NOT registered and will always fail.
    .Call("_treeowen_treeowen_ping", PACKAGE = pkgname)
    TRUE
  }, error = function(e) FALSE)
  .CPP_LOADED    <<- TRUE
  .CPP_AVAILABLE <<- ok
  invisible(NULL)
}
