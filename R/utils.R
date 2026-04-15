# treeowen/R/utils.R
#
# All internal helper functions (all "." prefixed, not exported unless noted).
# Source: treeowen_main.R v5.14
#
# Sections
# ─────────
#  §1   Package / Rcpp helpers           L142–227
#  §2   Progress & timing               L228–300
#  §3   Shapley weight helpers           L302–309
#  §4   Popcount helpers                 L315–337
#  §5   Tree preparation                 L344–591
#  §6   Input validation                 L598–620
#  §7   Bitset helpers                   L627–665
#  §8   Binary-tree builders (exported)  L672–754
#  §9   Tree annotation & common prep    L744–792
#  §10  Inner Shapley weight cache       L799–816
#  §11  Inner Shapley solvers (exact)    L818–854, L1234–1259
#  §12  Inner Shapley solvers (MC)       L1182–1296
#  §13  Outer exact recursion            L860–915
#  §14  Leaf-context precomputation      L924–994
#  §15  Row-loop runners                 L1001–1149
#  §16  Result constructor               L1151–1167
#  §17  Miscellaneous                    L1174–1179
#  §18  C++ / XPtr helpers               L206–221

# ── §1  Package / Rcpp helpers ───────────────────────────────────────────────

.pkg_ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# In an installed package the shared library is always pre-built.
# This function is a runtime fallback for the case where .onLoad could not
# confirm C availability (e.g. during R CMD check without compiled code).
.try_load_cpp <- function() {
  if (.CPP_LOADED) return(.CPP_AVAILABLE)
  # Probe the sentinel registered by useDynLib
  ok <- tryCatch({
    .Call("treeowen_ping", PACKAGE = "treeowen")
    TRUE
  }, error = function(e) FALSE)
  .CPP_LOADED    <<- TRUE
  .CPP_AVAILABLE <<- ok
  ok
}

# ── §18  C++ / XPtr helpers ──────────────────────────────────────────────────
# (placed early so treeowen() can reference them)

.trees_for_cpp <- function(dp_model) {
  lapply(dp_model$trees, function(tr) {
    tr$fcol_local[is.na(tr$fcol_local)] <- 0L
    c(tr, list(has_missing = dp_model$has_missing))
  })
}

.prepare_trees_xptr_once <- function(dp_model) {
  tlist <- .trees_for_cpp(dp_model)
  prepare_trees_xptr(tlist)
}

# backwards-compat alias
.prepare_trees_xptr_cached <- .prepare_trees_xptr_once

# ── §2  Progress & timing ─────────────────────────────────────────────────────

.fmt_hms <- function(sec) {
  if (!is.finite(sec) || is.na(sec)) return("NA")
  sec <- max(0, sec); h <- floor(sec / 3600)
  m <- floor((sec - 3600 * h) / 60); s <- floor(sec - 3600 * h - 60 * m)
  if (h > 0) sprintf("%dh%02dm%02ds", h, m, s) else sprintf("%dm%02ds", m, s)
}

.make_progress_state <- function()
  list(t_last = proc.time()[[3]], done_last = 0, rate_ema = NA_real_)

.update_eta_ema <- function(st, done, total, alpha = 0.2) {
  total <- as.numeric(total); done <- as.numeric(done)
  if (is.finite(total) && is.finite(done) && done >= total) {
    st$t_last <- proc.time()[[3]]; st$done_last <- done
    if (!is.finite(st$rate_ema)) st$rate_ema <- NA_real_
    return(list(state = st, eta = 0, rate = st$rate_ema))
  }
  t_now <- proc.time()[[3]]; dt <- t_now - st$t_last; dd <- done - st$done_last
  inst_rate <- if (is.finite(dt) && dt > 1e-9) dd / dt else NA_real_
  if (is.finite(inst_rate) && inst_rate > 0) {
    if (!is.finite(st$rate_ema)) st$rate_ema <- inst_rate
    else                         st$rate_ema <- alpha * inst_rate + (1 - alpha) * st$rate_ema
  }
  st$t_last <- t_now; st$done_last <- done
  eta <- if (is.finite(st$rate_ema) && st$rate_ema > 0 &&
             is.finite(total) && is.finite(done)) (total - done) / st$rate_ema else NA_real_
  list(state = st, eta = eta, rate = st$rate_ema)
}

.print_phase <- function(phase, pct, msg, elapsed, eta, newline = FALSE) {
  line <- sprintf("[%s] %5.1f%% | %s | elapsed=%s | eta=%s",
                  phase, pct, msg, .fmt_hms(elapsed), .fmt_hms(eta))
  ide_like <- nzchar(Sys.getenv("RSTUDIO")) || nzchar(Sys.getenv("VSCODE_PID")) ||
              nzchar(Sys.getenv("JPY_PARENT_PID")) || !interactive()
  if (ide_like) cat(line, "\n", sep = "")
  else if (newline) cat("\r", line, "\n", sep = "")
  else cat("\r", line, sep = "")
  flush.console()
}

.as_time_list <- function(tm) {
  tm_num <- as.numeric(tm); tm_names <- names(tm)
  pick <- function(cands, fallback_idx = NA_integer_) {
    if (!is.null(tm_names)) for (nm in cands) if (nm %in% tm_names) return(unname(tm[[nm]]))
    if (is.finite(fallback_idx) && length(tm_num) >= fallback_idx) return(unname(tm_num[[fallback_idx]]))
    NA_real_
  }
  list(elapsed = as.numeric(pick(c("elapsed","elapsed.sec"),             3L)),
       user    = as.numeric(pick(c("user.self","user","user.self.sec"),   1L)),
       system  = as.numeric(pick(c("sys.self","system","sys.self.sec"),   2L)))
}

# ── §3  Shapley weight helpers ────────────────────────────────────────────────

.safe_shapley_w <- function(s, m)
  exp(lfactorial(s) + lfactorial(m - s - 1L) - lfactorial(m))

.safe_outer_w <- function(t, K) .safe_shapley_w(t, K)
.safe_inner_w <- function(s, m) .safe_shapley_w(s, m)
.outer_w      <- .safe_outer_w
.inner_w      <- .safe_inner_w

# ── §4  Popcount helpers ──────────────────────────────────────────────────────

.popcount_int <- function(max_mask) {
  max_mask <- as.integer(max_mask); pc <- integer(max_mask + 1L)
  if (max_mask > 0L)
    for (mm in seq_len(max_mask))
      pc[mm + 1L] <- pc[bitwShiftR(mm, 1L) + 1L] + bitwAnd(mm, 1L)
  pc
}

.popcount_num <- function(x) {
  count <- 0L; xi <- x
  while (xi > 0) { xi <- xi - 2^(floor(log2(xi))); count <- count + 1L }
  count
}

.mask_has_bit_num <- function(masks_num, bp) {
  pow <- 2^(as.integer(bp) - 1L)
  (floor(as.numeric(masks_num) / pow) %% 2) == 1
}

# ── §5  Tree preparation ──────────────────────────────────────────────────────

.prepare_dp_model <- function(unified_model, feat_names) {
  stopifnot(inherits(unified_model, "model_unified"))
  model   <- unified_model$model
  Node    <- model$Node; Yes0 <- model$Yes; No0 <- model$No
  Feature <- model$Feature; Split <- model$Split
  Pred    <- model$Prediction
  if (is.factor(Pred)) Pred <- as.numeric(as.character(Pred))
  Cover <- model$Cover
  has_missing <- "Missing" %in% names(model)
  Missing0    <- if (has_missing) model$Missing else rep(NA_integer_, nrow(model))
  has_tree    <- "Tree" %in% names(model)
  TreeID      <- if (has_tree) model$Tree else rep(NA_integer_, nrow(model))
  fcol  <- match(Feature, feat_names)
  roots <- which(Node == 0)
  N     <- nrow(model)

  .is_row_pointer <- function(v) {
    v <- v[is.finite(v) & !is.na(v)]; if (!length(v)) return(TRUE)
    all(v >= 1L & v <= N)
  }
  child_is_row <- isTRUE(.is_row_pointer(Yes0) && .is_row_pointer(No0) &&
                          (if (has_missing) .is_row_pointer(Missing0) else TRUE))
  unused_env <- new.env(parent = emptyenv()); unused_env$unused <- rep(TRUE, N)

  get_kids_factory <- function(root_row) {
    if (child_is_row) {
      return(function(r) {
        kids <- c(Yes0[r], No0[r])
        if (has_missing) { mr <- Missing0[r]; if (!is.na(mr) && mr != Yes0[r] && mr != No0[r]) kids <- c(kids, mr) }
        as.integer(kids[!is.na(kids)])
      })
    }
    if (has_tree) {
      tr_id <- TreeID[root_row]; rows_in_tree <- which(TreeID == tr_id)
      node_to_row <- setNames(rows_in_tree, as.character(Node[rows_in_tree]))
      return(function(r) {
        kids_id <- c(Yes0[r], No0[r])
        if (has_missing) { mr <- Missing0[r]; if (!is.na(mr) && mr != Yes0[r] && mr != No0[r]) kids_id <- c(kids_id, mr) }
        kids_id <- kids_id[!is.na(kids_id)]; if (!length(kids_id)) return(integer(0))
        out <- unname(node_to_row[as.character(kids_id)]); as.integer(out[!is.na(out)])
      })
    }
    node_to_row <- new.env(parent = emptyenv())
    pick_row_for_node <- function(node_id) {
      node_id <- as.integer(node_id); if (is.na(node_id)) return(NA_integer_)
      key <- as.character(node_id)
      hit <- get0(key, envir = node_to_row, inherits = FALSE); if (!is.null(hit)) return(as.integer(hit))
      cand <- which(unused_env$unused & Node == node_id); if (!length(cand)) return(NA_integer_)
      rr <- cand[1L]; assign(key, rr, envir = node_to_row); unused_env$unused[rr] <- FALSE; as.integer(rr)
    }
    unused_env$unused[root_row] <- FALSE
    assign(as.character(Node[root_row]), root_row, envir = node_to_row)
    return(function(r) {
      kids_id <- c(Yes0[r], No0[r])
      if (has_missing) { mr <- Missing0[r]; if (!is.na(mr) && mr != Yes0[r] && mr != No0[r]) kids_id <- c(kids_id, mr) }
      kids_id <- kids_id[!is.na(kids_id)]; if (!length(kids_id)) return(integer(0))
      out <- integer(0)
      for (kid_id in kids_id) { rr <- pick_row_for_node(kid_id); if (!is.na(rr)) out <- c(out, rr) }
      as.integer(out)
    })
  }

  max_nn <- 0L; trees <- vector("list", length(roots))
  for (tt in seq_along(roots)) {
    rt <- roots[tt]; get_kids <- get_kids_factory(rt)
    stk_cap <- 2L * N + 8L
    stk   <- integer(stk_cap); sp  <- 1L; stk[1L] <- rt
    seen  <- rep(FALSE, N)
    nodes_buf <- integer(N); nb <- 0L
    while (sp > 0L) {
      r  <- stk[sp]; sp <- sp - 1L
      if (r < 1L || r > N || seen[r]) next
      seen[r] <- TRUE
      nb <- nb + 1L; nodes_buf[nb] <- r
      j  <- fcol[r]
      if (!is.na(j)) {
        kids <- get_kids(r)
        for (k in kids) {
          if (sp >= stk_cap) { stk_cap <- stk_cap * 2L; stk <- c(stk, integer(stk_cap - length(stk))) }
          sp <- sp + 1L; stk[sp] <- k
        }
      }
    }
    nodes <- nodes_buf[seq_len(nb)]
    nn    <- nb; if (nn > max_nn) max_nn <- nn
    local_of_global <- integer(N); local_of_global[nodes] <- seq_len(nn)
    cap2   <- nn * 3L + 8L
    stk_n  <- integer(cap2); stk_e <- logical(cap2); sp2 <- 1L
    stk_n[1L] <- rt; stk_e[1L] <- FALSE
    seen2  <- rep(FALSE, N)
    po_buf <- integer(nn); pb <- 0L
    while (sp2 > 0L) {
      r        <- stk_n[sp2]; expanded <- stk_e[sp2]; sp2 <- sp2 - 1L
      if (!expanded) {
        if (seen2[r]) next; seen2[r] <- TRUE
        sp2 <- sp2 + 1L; stk_n[sp2] <- r; stk_e[sp2] <- TRUE
        j   <- fcol[r]
        if (!is.na(j)) {
          kids <- get_kids(r)
          for (k in rev(kids)) { sp2 <- sp2 + 1L; stk_n[sp2] <- k; stk_e[sp2] <- FALSE }
        }
      } else {
        pb <- pb + 1L; po_buf[pb] <- local_of_global[r]
      }
    }
    post_local <- po_buf[seq_len(pb)]
    nkids <- integer(nn); kid1 <- integer(nn); kid2 <- integer(nn)
    p1 <- numeric(nn); p2 <- numeric(nn)
    yesL <- integer(nn); noL <- integer(nn); missL <- integer(nn)
    has_miss_dir <- logical(nn)
    for (i in seq_len(nn)) {
      g <- nodes[i]; j <- fcol[g]
      if (is.na(j)) { nkids[i] <- 0L; next }
      yl_row <- NA_integer_; nr_row <- NA_integer_; ms_row <- NA_integer_
      if (child_is_row) {
        yl_row <- Yes0[g]; nr_row <- No0[g]; if (has_missing) ms_row <- Missing0[g]
      } else {
        map_one <- function(id) {
          if (is.na(id)) return(NA_integer_); kids <- get_kids(g); if (!length(kids)) return(NA_integer_)
          rr <- kids[which(Node[kids] == id)]; if (length(rr)) rr[1L] else NA_integer_
        }
        yl_row <- map_one(Yes0[g]); nr_row <- map_one(No0[g])
        if (has_missing) ms_row <- map_one(Missing0[g])
      }
      ylL <- if (!is.na(yl_row)) local_of_global[yl_row] else 0L
      nrL <- if (!is.na(nr_row)) local_of_global[nr_row] else 0L
      yesL[i] <- ylL; noL[i] <- nrL
      mrL <- 0L
      if (has_missing && !is.na(ms_row)) {
        mrL <- local_of_global[ms_row]; has_miss_dir[i] <- TRUE
      }
      missL[i] <- mrL
      kids <- as.integer(c(yl_row, nr_row)[!is.na(c(yl_row, nr_row))])
      if (!length(kids)) { nkids[i] <- 0L
      } else if (length(kids) == 1L) {
        nkids[i] <- 1L; kid1[i] <- local_of_global[kids[1]]; p1[i] <- 1
      } else {
        nkids[i] <- 2L; yl2 <- kids[1]; nr2 <- kids[2]
        kid1[i] <- local_of_global[yl2]; kid2[i] <- local_of_global[nr2]
        denom <- Cover[yl2] + Cover[nr2]
        if (!is.finite(denom) || denom <= 0) { p1[i] <- 0.5; p2[i] <- 0.5 }
        else { p1[i] <- Cover[yl2] / denom; p2[i] <- Cover[nr2] / denom }
      }
    }
    trees[[tt]] <- list(
      root_global = rt, root_local = local_of_global[rt],
      nodes_global = nodes, post_local = post_local,
      fcol_local = fcol[nodes], split_local = Split[nodes], pred_local = Pred[nodes],
      yesL = yesL, noL = noL, missL = missL, has_miss_dir = has_miss_dir,
      nkids = nkids, kid1 = kid1, kid2 = kid2, p1 = p1, p2 = p2
    )
  }
  list(feat_names = feat_names, N = N, has_missing = has_missing, trees = trees,
       used_feat_cols = sort(unique(na.omit(fcol))), max_nn = as.integer(max_nn))
}

.eval_v_masks_forest_chunked_prepared <- function(
    dp_model, x_row, chunk_len,
    outer_known_feat, group_bitpos, col_masks_num
) {
  chunk_len <- as.integer(chunk_len)
  if (chunk_len <= 0L) return(numeric(0))
  total <- rep(0, chunk_len); max_nn <- dp_model$max_nn
  if (!is.finite(max_nn) || max_nn <= 0L) return(total)
  Vbuf <- matrix(0, nrow = max_nn, ncol = chunk_len)
  maxj <- length(outer_known_feat)
  if (length(group_bitpos)  != maxj)      stop("group_bitpos length mismatch")
  if (length(col_masks_num) != chunk_len) stop("col_masks length mismatch")
  for (tr in dp_model$trees) {
    nn <- length(tr$nodes_global); if (nn <= 0L) next
    Vbuf[seq_len(nn), ] <- 0
    for (loc in tr$post_local) {
      j  <- tr$fcol_local[loc]
      nk <- tr$nkids[loc]
      if (is.na(j)) { Vbuf[loc, ] <- tr$pred_local[loc]; next }
      if      (nk == 0L) Vbuf[loc, ] <- 0
      else if (nk == 1L) Vbuf[loc, ] <- Vbuf[tr$kid1[loc], ]
      else               Vbuf[loc, ] <- tr$p1[loc] * Vbuf[tr$kid1[loc], ] +
                                        tr$p2[loc] * Vbuf[tr$kid2[loc], ]
      if (j < 1L || j > maxj) next
      is_outer <- isTRUE(outer_known_feat[j])
      bp       <- group_bitpos[j]
      if (!is_outer && bp == 0L) next
      xv <- x_row[[j]]
      if (is.na(xv)) {
        if (dp_model$has_missing && isTRUE(tr$has_miss_dir[loc])) {
          mrL <- tr$missL[loc]
          if (is_outer) { Vbuf[loc, ] <- Vbuf[mrL, ]
          } else { sel <- .mask_has_bit_num(col_masks_num, bp); if (any(sel)) Vbuf[loc, sel] <- Vbuf[mrL, sel] }
        }
      } else {
        det_child_local <- if (xv <= tr$split_local[loc]) tr$yesL[loc] else tr$noL[loc]
        if (det_child_local == 0L) next
        if (is_outer) { Vbuf[loc, ] <- Vbuf[det_child_local, ]
        } else { sel <- .mask_has_bit_num(col_masks_num, bp); if (any(sel)) Vbuf[loc, sel] <- Vbuf[det_child_local, sel] }
      }
    }
    total <- total + Vbuf[tr$root_local, ]
  }
  total
}

.model_tree_stats <- function(unified_model) {
  model <- unified_model$model; ntrees <- sum(model$Node == 0, na.rm = TRUE); total <- nrow(model)
  list(ntrees = ntrees, total_nodes = total, avg_nodes = total / max(1, ntrees))
}

# ── §6  Input validation ──────────────────────────────────────────────────────

.validate_inputs <- function(unified_model, x, groups,
                             K_max = TREEOWEN_K_DEFAULT, m_max = TREEOWEN_M_DEFAULT) {
  K_max <- as.integer(K_max); m_max <- as.integer(m_max)
  if (!inherits(unified_model, "model_unified")) stop("unified_model must be 'model_unified' class")
  if (!is.data.frame(x) && !is.matrix(x)) stop("x must be data.frame or matrix")
  if (nrow(x) <= 0) stop("x has zero rows"); if (ncol(x) <= 0) stop("x has zero columns")
  if (!is.list(groups) || !length(groups)) stop("groups must be non-empty list")
  groups <- Filter(length, groups); if (!length(groups)) stop("groups must contain at least one non-empty group")
  K <- length(groups)
  if (K > .TREEOWEN_K_HARD_MAX) stop(sprintf("K=%d exceeds hard maximum K=%d.", K, .TREEOWEN_K_HARD_MAX))
  if (K > K_max) stop(sprintf("K=%d exceeds K_max=%d (raise K_max up to %d if intended).", K, K_max, .TREEOWEN_K_HARD_MAX))
  all_feats <- unlist(groups, use.names = FALSE); feat_names <- colnames(x)
  if (length(all_feats) != length(unique(all_feats))) stop("Features appear in multiple groups.")
  unknown <- setdiff(all_feats, feat_names)
  if (length(unknown)) stop(sprintf("Groups contain unknown feature(s): %s", paste(head(unknown,10), collapse=", ")))
  uncovered <- setdiff(feat_names, all_feats)
  if (length(uncovered)) stop(sprintf("TreeOwen requires FULL partition: %d x column(s) not assigned: %s",
                                       length(uncovered), paste(head(uncovered,10), collapse=", ")))
  max_m <- max(vapply(groups, length, 1L))
  if (max_m > .TREEOWEN_M_HARD_MAX) stop(sprintf("Group size m=%d exceeds hard maximum m=%d.", max_m, .TREEOWEN_M_HARD_MAX))
  if (max_m > m_max) stop(sprintf("Group size m=%d exceeds m_max=%d (raise m_max up to %d if intended).", max_m, m_max, .TREEOWEN_M_HARD_MAX))
  invisible(TRUE)
}

# ── §7  Bitset helpers ────────────────────────────────────────────────────────

.bitset_nblocks <- function(K) as.integer((as.integer(K) + .BITS_PER_BLOCK - 1L) %/% .BITS_PER_BLOCK)
.bitset_empty   <- function(K) rep.int(0L, .bitset_nblocks(K))

.bitset_from_ids <- function(ids, K) {
  ids <- as.integer(ids); if (!length(ids)) return(.bitset_empty(K))
  K <- as.integer(K); m <- .bitset_empty(K)
  for (g in ids) {
    blk <- as.integer((g - 1L) %/% .BITS_PER_BLOCK) + 1L
    off <- as.integer((g - 1L) %%  .BITS_PER_BLOCK)
    m[[blk]] <- bitwOr(m[[blk]], bitwShiftL(1L, off))
  }
  m
}

.bitset_or <- function(a, b) { if (length(a) != length(b)) stop("bitset_or: length mismatch"); bitwOr(a, b) }

.bitset_count <- function(mask, K)
  sum(vapply(seq_len(K), function(g) .bitset_has(mask, g), logical(1L)))

.bitset_has <- function(mask, g) {
  g <- as.integer(g)
  blk <- as.integer((g - 1L) %/% .BITS_PER_BLOCK) + 1L
  off  <- as.integer((g - 1L) %%  .BITS_PER_BLOCK)
  bitwAnd(mask[[blk]], bitwShiftL(1L, off)) != 0L
}

.groups_to_mask       <- function(group_ids, K) .bitset_from_ids(group_ids, K)

.mask_to_feat_logical <- function(group_mask, group_pos, maxj, K) {
  out <- logical(maxj)
  for (g in seq_len(K)) {
    if (.bitset_has(group_mask, g)) {
      feats <- group_pos[[g]]; feats <- feats[feats >= 1L & feats <= maxj]
      if (length(feats)) out[feats] <- TRUE
    }
  }
  out
}

# ── §8  Binary-tree builders (exported) ─────────────────────────────────────

.build_balanced_binary_tree <- function(leaf_ids) {
  leaf_ids <- as.integer(leaf_ids); uid_env <- new.env(parent = emptyenv()); uid_env$n <- 0L
  build <- function(ids) {
    uid_env$n <- uid_env$n + 1L; uid <- uid_env$n
    if (length(ids) == 1L) return(list(type="leaf", uid=uid, id=ids, children=NULL, leaves=ids))
    mid <- length(ids) %/% 2L; left <- build(ids[seq_len(mid)]); right <- build(ids[(mid+1L):length(ids)])
    list(type="node", uid=uid, id=NA_integer_, children=list(left, right), leaves=c(left$leaves, right$leaves))
  }
  build(leaf_ids)
}

#' Build a Balanced Binary Tree over Group Indices
#'
#' Constructs the auxiliary balanced binary tree \eqn{\mathcal{B}} used by
#' \code{\link{treeowen}} to organise the outer Shapley-weighted summation.
#' Exposed for users who need a custom \code{hierarchy} argument.
#'
#' The function accepts a NESTED LIST where each leaf is a single integer
#' group index and each internal node is \code{list(left = ..., right = ...)}.
#'
#' @param hier A nested list (or single integer) representing the binary tree.
#'   Each leaf must be a length-1 integer; each internal node must have exactly
#'   two children named \code{left}/\code{right} or as positional elements.
#' @return A normalised nested list representing the annotated binary tree.
#' @seealso \code{\link{build_hierarchy_tree_from_layers}}, \code{\link{treeowen}}
#' @export
build_hierarchy_tree_binary <- function(hier) {
  uid_env <- new.env(parent = emptyenv()); uid_env$n <- 0L
  norm_node <- function(x) {
    if (is.numeric(x) && length(x) == 1L) {
      uid_env$n <- uid_env$n + 1L; g <- as.integer(x)
      return(list(type="leaf", uid=uid_env$n, id=g, children=NULL, leaves=g))
    }
    if (!is.list(x)) stop("binary hierarchy: leaf int or list(left/right)")
    kids <- if (!is.null(x$left) && !is.null(x$right)) list(x$left, x$right)
            else if (length(x) == 2L) list(x[[1]], x[[2]])
            else stop("binary hierarchy: internal node must have exactly 2 children")
    uid_env$n <- uid_env$n + 1L; uid <- uid_env$n
    children <- lapply(kids, norm_node)
    leaves   <- unlist(lapply(children, `[[`, "leaves"), use.names = FALSE)
    list(type="node", uid=uid, id=NA_integer_, children=children, leaves=as.integer(leaves))
  }
  norm_node(hier)
}

#' Build the Auxiliary Binary Tree from a Hierarchy Specification
#'
#' Constructs the auxiliary binary tree \eqn{\mathcal{B}} from either a
#' balanced default (\code{hierarchy = NULL}) or a user-supplied specification.
#' This is the entry point called by \code{\link{treeowen}}.
#'
#' @param K Integer. Number of groups.
#' @param hierarchy Controls the tree structure. Accepted values:
#'   \describe{
#'     \item{\code{NULL}}{(default) A balanced binary tree over the \eqn{K} groups.}
#'     \item{Named list}{A fully specified nested-list tree, compatible with
#'       the output of \code{build_hierarchy_tree_binary()}.}
#'     \item{Numeric vector or \eqn{K}-column integer matrix}{A layer-based
#'       merge specification: each row merges groups sharing the same integer label.}
#'   }
#' @return A nested list representing the annotated binary tree.
#' @seealso \code{\link{build_hierarchy_tree_binary}}, \code{\link{treeowen}}
#' @export
build_hierarchy_tree_from_layers <- function(K, hierarchy) {
  K <- as.integer(K); uid_env <- new.env(parent = emptyenv()); uid_env$n <- 0L
  mk_leaf <- function(g) {
    uid_env$n <- uid_env$n + 1L
    list(type="leaf", uid=uid_env$n, id=as.integer(g), children=NULL, leaves=as.integer(g))
  }
  if (is.null(hierarchy)) return(.build_balanced_binary_tree(seq_len(K)))
  if (is.list(hierarchy) && !is.matrix(hierarchy)) {
    norm_any <- function(x) {
      if (is.numeric(x) && length(x) == 1L) return(mk_leaf(as.integer(x)))
      if (!is.list(x)) stop("hierarchy tree: node must be leaf int or list")
      kids <- if (!is.null(x$children)) x$children
              else if (!is.null(x$left) && !is.null(x$right)) list(x$left, x$right)
              else unname(x)
      uid_env$n <- uid_env$n + 1L; uid <- uid_env$n
      children <- lapply(kids, norm_any)
      leaves   <- unlist(lapply(children, `[[`, "leaves"), use.names = FALSE)
      list(type="node", uid=uid, id=NA_integer_, children=children, leaves=as.integer(leaves))
    }
    return(norm_any(hierarchy))
  }
  h <- hierarchy
  if (is.vector(h) && is.numeric(h)) h <- matrix(as.integer(h), nrow = 1L)
  if (!is.matrix(h) || !is.numeric(h)) stop("hierarchy must be NULL/list/vector/matrix")
  h <- apply(h, 2, as.integer); if (ncol(h) != K) stop("hierarchy matrix must have K columns")
  current_nodes <- lapply(seq_len(K), mk_leaf)
  for (r in seq_len(nrow(h))) {
    labels <- as.integer(h[r, ]); ulabels <- sort(unique(labels)); parents <- vector("list", length(ulabels))
    for (ii in seq_along(ulabels)) {
      lab  <- ulabels[ii]; kids <- Filter(function(nd) any(labels[nd$leaves] == lab), current_nodes)
      uid_env$n <- uid_env$n + 1L; uid <- uid_env$n
      leaves <- unlist(lapply(kids, `[[`, "leaves"), use.names = FALSE)
      parents[[ii]] <- list(type="node", uid=uid, id=NA_integer_, children=kids, leaves=as.integer(leaves))
    }
    current_nodes <- Filter(Negate(is.null), parents)
  }
  if (length(current_nodes) == 1L) return(current_nodes[[1L]])
  uid_env$n <- uid_env$n + 1L; uid <- uid_env$n
  leaves <- unlist(lapply(current_nodes, `[[`, "leaves"), use.names = FALSE)
  list(type="node", uid=uid, id=NA_integer_, children=current_nodes, leaves=as.integer(leaves))
}

# ── §9  Tree annotation & common prep ────────────────────────────────────────

.annotate_tree <- function(tree, K) {
  nodes <- list()
  rec <- function(node, parent_uid = 0L, depth = 0L) {
    node$parent_uid <- as.integer(parent_uid); node$depth <- as.integer(depth)
    nodes[[as.character(node$uid)]] <<- node
    if (identical(node$type, "node")) for (ch in node$children) rec(ch, parent_uid = node$uid, depth = depth + 1L)
  }
  rec(tree)
  ann <- list(root = nodes[[as.character(tree$uid)]], by_uid = nodes)
  attr(ann, "K") <- as.integer(K); ann
}

.prepare_common <- function(unified_model, x, groups, max_bytes, chunk_size_inner,
                            K_max = TREEOWEN_K_DEFAULT, m_max = TREEOWEN_M_DEFAULT) {
  stopifnot(inherits(unified_model, "model_unified"))
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if (!is.null(unified_model$feature_names))
    x <- x[, intersect(colnames(x), unified_model$feature_names), drop = FALSE]
  feat_names <- colnames(x); p <- ncol(x); n <- nrow(x)
  if (!p) stop("No usable columns in x after aligning to model features.")
  groups <- Filter(length, groups); K <- length(groups)
  if (K > .TREEOWEN_K_HARD_MAX) stop(sprintf("K=%d exceeds hard max %d", K, .TREEOWEN_K_HARD_MAX))
  if (K > as.integer(K_max))    stop(sprintf("K=%d > K_max=%d", K, as.integer(K_max)))
  fpos      <- setNames(seq_len(p), feat_names)
  group_pos <- lapply(groups, function(g) as.integer(fpos[as.character(g)]))
  if (any(vapply(group_pos, anyNA, logical(1)))) stop("Some group feature names not found in x.")
  all_group_cols <- sort(unique(unlist(group_pos, use.names = FALSE)))
  if (!identical(all_group_cols, seq_len(p))) stop("Requires FULL partition.")
  max_m <- max(vapply(group_pos, length, 1L))
  if (max_m > .TREEOWEN_M_HARD_MAX) stop(sprintf("Group size m=%d > hard max %d", max_m, .TREEOWEN_M_HARD_MAX))
  if (max_m > as.integer(m_max))    stop(sprintf("Group size m=%d > m_max=%d", max_m, as.integer(m_max)))
  dp_model <- .prepare_dp_model(unified_model, feat_names)
  maxj     <- if (length(dp_model$used_feat_cols)) max(dp_model$used_feat_cols) else 0L
  ts       <- .model_tree_stats(unified_model)
  max_bytes <- as.numeric(max_bytes); if (!is.finite(max_bytes) || max_bytes <= 0) max_bytes <- 512 * 1024^2
  max_nn <- as.integer(dp_model$max_nn); if (!is.finite(max_nn) || max_nn <= 0L) max_nn <- 1L
  chunk_cap  <- max(1L, floor(max_bytes / (max_nn * 8)))
  chunk_req  <- as.integer(chunk_size_inner); if (!is.finite(chunk_req) || chunk_req <= 0L) chunk_req <- 131072L
  chunk_size_inner <- min(chunk_req, as.integer(chunk_cap))
  list(x = x, feat_names = feat_names, p = p, n = n, K = K, groups = groups,
       group_pos = group_pos, max_m = max_m, dp_model = dp_model, maxj = maxj, ts = ts,
       max_bytes = max_bytes, max_nn = max_nn, chunk_cap = chunk_cap,
       chunk_size_inner = chunk_size_inner, chunk_size_inner_req = chunk_req)
}

# ── §10  Inner Shapley weight cache ──────────────────────────────────────────

#' Clear the memoised inner Shapley weight cache.
#'
#' Call between large analyses (different group sizes) to free memory.
#' @export
clear_inner_enum_cache <- function() {
  rm(list = ls(.INNER_ENUM_CACHE), envir = .INNER_ENUM_CACHE)
  invisible(NULL)
}

.get_inner_cache <- function(m) {
  m <- as.integer(m); key <- as.character(m)
  hit <- get0(key, envir = .INNER_ENUM_CACHE, inherits = FALSE); if (!is.null(hit)) return(hit)
  if (m == 0L) { out <- list(w_by_s = numeric(0), use_table = FALSE); assign(key, out, envir = .INNER_ENUM_CACHE); return(out) }
  use_table <- (m <= 30L)
  w_by_s    <- vapply(0L:(m - 1L), function(s) .safe_shapley_w(s, m), numeric(1))
  out <- list(w_by_s = w_by_s, use_table = use_table)
  if (use_table) {
    M <- bitwShiftL(1L, m); pc <- .popcount_int(M - 1L)
    masks_wo  <- vector("list", m); all_masks <- 0L:(M - 1L)
    for (i in seq_len(m)) { bit_i <- bitwShiftL(1L, i - 1L); masks_wo[[i]] <- all_masks[bitwAnd(all_masks, bit_i) == 0L] }
    out$M <- M; out$pc <- pc; out$masks_wo <- masks_wo
  }
  assign(key, out, envir = .INNER_ENUM_CACHE); out
}

# ── §11  Inner Shapley solvers (exact) ───────────────────────────────────────

.inner_shapley_in_group_given_context <- function(
    dp_model, x_row, g, group_pos, outer_group_mask,
    maxj, p,
    chunk_size        = 131072L,
    max_bytes         = 512 * 1024^2,
    inner_bitmask_max = TREEOWEN_INNER_BITMASK_DEFAULT,
    use_cpp           = FALSE,
    trees_xptr        = NULL
) {
  feats_g <- group_pos[[g]]; feats_g <- feats_g[feats_g >= 1L & feats_g <= maxj]
  if (length(feats_g) == 0L) return(numeric(p))
  if (length(feats_g) > 60L) stop(sprintf("Inner exact hard limit: m=%d > 60.", length(feats_g)))
  K <- length(group_pos)
  outer_known_cols <- integer(0)
  for (gg in seq_len(K)) {
    if (.bitset_has(outer_group_mask, gg)) {
      fc <- group_pos[[gg]]; fc <- fc[fc >= 1L & fc <= maxj]
      if (length(fc)) outer_known_cols <- c(outer_known_cols, fc)
    }
  }
  outer_known_cols <- unique(outer_known_cols)
  if (!use_cpp || !.CPP_AVAILABLE || is.null(trees_xptr) ||
      !exists("inner_shapley_exact_cpp", mode = "function"))
    stop("[TreeOwen] inner_shapley_exact_cpp unavailable. Ensure C++ is loaded.")
  as.numeric(inner_shapley_exact_cpp(
    trees_xptr       = trees_xptr,
    x_row            = as.numeric(x_row),
    feats_g          = as.integer(feats_g),
    outer_known_cols = as.integer(outer_known_cols),
    maxj             = as.integer(maxj),
    p                = as.integer(p),
    chunk_size       = as.integer(chunk_size)
  ))
}

.inner_shapley_in_group_given_context_direct <- function(
    dp_model, x_row, g, group_pos, outer_known_cols,
    maxj, p,
    chunk_size        = 131072L,
    max_bytes         = 512 * 1024^2,
    inner_bitmask_max = TREEOWEN_INNER_BITMASK_DEFAULT,
    use_cpp           = FALSE,
    trees_xptr        = NULL
) {
  feats_g <- group_pos[[g]]; feats_g <- feats_g[feats_g >= 1L & feats_g <= maxj]
  if (length(feats_g) == 0L) return(numeric(p))
  if (length(feats_g) > 60L) stop(sprintf("Inner exact hard limit: m=%d > 60.", length(feats_g)))
  if (!use_cpp || !.CPP_AVAILABLE || is.null(trees_xptr) ||
      !exists("inner_shapley_exact_cpp", mode = "function"))
    stop("[TreeOwen] inner_shapley_exact_cpp unavailable. Ensure C++ is loaded.")
  as.numeric(inner_shapley_exact_cpp(
    trees_xptr       = trees_xptr,
    x_row            = as.numeric(x_row),
    feats_g          = as.integer(feats_g),
    outer_known_cols = as.integer(outer_known_cols),
    maxj             = as.integer(maxj),
    p                = as.integer(p),
    chunk_size       = as.integer(chunk_size)
  ))
}

# ── §12  Inner Shapley solvers (Monte Carlo) ──────────────────────────────────

.inner_perm_mc_one_group <- function(
    dp_model, x_row, g, group_pos, outer_group_mask, maxj,
    n_inner_mc       = 64L,
    inner_antithetic = TRUE,
    target_se        = NA_real_,
    min_mc           = 32L,
    max_mc           = 512L,
    use_cpp          = FALSE,
    trees_xptr       = NULL
) {
  feats <- group_pos[[g]]; feats <- feats[feats >= 1L & feats <= maxj]
  m     <- length(feats); if (m == 0L) return(numeric(0))
  K     <- length(group_pos)
  outer_known_cols <- integer(0)
  for (gg in seq_len(K)) {
    if (.bitset_has(outer_group_mask, gg)) {
      fc <- group_pos[[gg]]; fc <- fc[fc >= 1L & fc <= maxj]
      if (length(fc)) outer_known_cols <- c(outer_known_cols, fc)
    }
  }
  outer_known_cols <- unique(outer_known_cols)
  n_inner_mc <- as.integer(n_inner_mc); if (!is.finite(n_inner_mc) || n_inner_mc <= 0L) n_inner_mc <- 64L
  min_mc     <- as.integer(min_mc);     if (!is.finite(min_mc)     || min_mc     <= 0L) min_mc     <- 32L
  max_mc     <- as.integer(max_mc);     if (!is.finite(max_mc)     || max_mc < min_mc)  max_mc     <- max(min_mc, 512L)
  if (!use_cpp || !.CPP_AVAILABLE || is.null(trees_xptr) ||
      !exists("inner_mc_perm_cpp", mode = "function"))
    stop("[TreeOwen] inner_mc_perm_cpp unavailable. Ensure C++ is loaded.")
  as.numeric(inner_mc_perm_cpp(
    trees_xptr       = trees_xptr,
    x_row            = as.numeric(x_row),
    feats_g          = as.integer(feats),
    outer_known_cols = as.integer(outer_known_cols),
    maxj             = as.integer(maxj),
    p                = as.integer(length(x_row)),
    n_perm           = as.integer(n_inner_mc),
    antithetic       = isTRUE(inner_antithetic),
    target_se        = if (is.finite(target_se) && target_se > 0) target_se else -1.0,
    min_perm         = as.integer(min_mc)
  ))
}

.inner_perm_mc_one_group_direct <- function(
    dp_model, x_row, g, group_pos, outer_known_cols,
    maxj, p,
    n_inner_mc       = 64L,
    inner_antithetic = TRUE,
    target_se        = NA_real_,
    min_mc           = 32L,
    max_mc           = 512L,
    use_cpp          = FALSE,
    trees_xptr       = NULL
) {
  feats <- group_pos[[g]]; feats <- feats[feats >= 1L & feats <= maxj]
  m     <- length(feats); if (m == 0L) return(numeric(p))
  n_inner_mc <- as.integer(n_inner_mc); if (!is.finite(n_inner_mc) || n_inner_mc <= 0L) n_inner_mc <- 64L
  min_mc     <- as.integer(min_mc);     if (!is.finite(min_mc)     || min_mc     <= 0L) min_mc     <- 32L
  max_mc     <- as.integer(max_mc);     if (!is.finite(max_mc)     || max_mc < min_mc)  max_mc     <- max(min_mc, 512L)
  if (!use_cpp || !.CPP_AVAILABLE || is.null(trees_xptr) ||
      !exists("inner_mc_perm_cpp", mode = "function"))
    stop("[TreeOwen] inner_mc_perm_cpp unavailable. Ensure C++ is loaded.")
  as.numeric(inner_mc_perm_cpp(
    trees_xptr       = trees_xptr,
    x_row            = as.numeric(x_row),
    feats_g          = as.integer(feats),
    outer_known_cols = as.integer(outer_known_cols),
    maxj             = as.integer(maxj),
    p                = as.integer(p),
    n_perm           = as.integer(n_inner_mc),
    antithetic       = isTRUE(inner_antithetic),
    target_se        = if (is.finite(target_se) && target_se > 0) target_se else -1.0,
    min_perm         = as.integer(min_mc)
  ))
}

# ── §13  Outer exact recursion ────────────────────────────────────────────────

.hierowen_outer_exact_node <- function(
    node_uid,
    ann,
    outer_group_mask,
    leaf_solver,
    max_children_exact = 20L,
    ...
) {
  K    <- as.integer(attr(ann, "K"))
  node <- ann$by_uid[[as.character(node_uid)]]
  if (identical(node$type, "leaf")) {
    g <- as.integer(node$id)
    return(leaf_solver(g = g, outer_group_mask = outer_group_mask, ...))
  }
  kids <- node$children; d <- length(kids)
  if (d < 1L) stop("Internal node has no children.")
  max_children_exact <- as.integer(max_children_exact)
  if (!is.finite(max_children_exact) || max_children_exact < 2L) max_children_exact <- 20L
  if (d > max_children_exact)
    stop(sprintf("HIEROwen exact: node degree d=%d > max_children_exact=%d.", d, max_children_exact))
  child_masks <- lapply(kids, function(ch) .groups_to_mask(ch$leaves, K = K))
  out <- NULL
  for (ci in seq_len(d)) {
    others <- setdiff(seq_len(d), ci); d_minus <- length(others)
    if (d_minus == 0L) {
      phi_ci <- .hierowen_outer_exact_node(
        node_uid = kids[[ci]]$uid, ann = ann, outer_group_mask = outer_group_mask,
        leaf_solver = leaf_solver, max_children_exact = max_children_exact, ...)
      if (is.null(out)) out <- rep(0, length(phi_ci)); out <- out + phi_ci; next
    }
    max_tmask <- bitwShiftL(1L, d_minus) - 1L
    for (tmask in 0L:max_tmask) {
      bits       <- as.logical(intToBits(tmask))[seq_len(d_minus)]
      T_children <- if (!any(bits)) integer(0) else others[which(bits)]
      mask_T     <- if (length(T_children)) Reduce(.bitset_or, child_masks[T_children]) else .bitset_empty(K)
      wT         <- .safe_shapley_w(length(T_children), d)
      phi_ci <- .hierowen_outer_exact_node(
        node_uid = kids[[ci]]$uid, ann = ann,
        outer_group_mask = .bitset_or(outer_group_mask, mask_T),
        leaf_solver = leaf_solver, max_children_exact = max_children_exact, ...)
      if (is.null(out)) out <- rep(0, length(phi_ci)); out <- out + wT * phi_ci
    }
  }
  out
}

# ── §14  Leaf-context precomputation ─────────────────────────────────────────

.precompute_leaf_contexts <- function(ann, K, group_pos, maxj) {
  contexts <- list()
  recurse <- function(node_uid, outer_group_mask) {
    node <- ann$by_uid[[as.character(node_uid)]]
    if (identical(node$type, "leaf")) {
      g <- as.integer(node$id)
      n_outer <- .bitset_count(outer_group_mask, K)
      weight  <- .safe_shapley_w(n_outer, K)
      okc_parts <- vector("list", K); n_parts <- 0L
      for (gg in seq_len(K)) {
        if (.bitset_has(outer_group_mask, gg)) {
          fc <- group_pos[[gg]]; fc <- fc[fc >= 1L & fc <= maxj]
          if (length(fc)) { n_parts <- n_parts + 1L; okc_parts[[n_parts]] <- fc }
        }
      }
      okc <- if (n_parts == 0L) integer(0) else unique(unlist(okc_parts[seq_len(n_parts)], use.names = FALSE))
      contexts[[length(contexts) + 1L]] <<- list(g = g, outer_known_cols = okc, weight = weight)
      return(invisible(NULL))
    }
    kids <- node$children; d <- length(kids)
    if (d < 1L) return(invisible(NULL))
    for (ci in seq_len(d)) {
      sibling_leaves <- unlist(
        lapply(kids[setdiff(seq_len(d), ci)], function(ch) ch$leaves), use.names = FALSE)
      n_sib <- length(sibling_leaves)
      if (n_sib == 0L) { recurse(kids[[ci]]$uid, outer_group_mask); next }
      if (n_sib >= 30L)
        stop(sprintf(paste0(
          "treeowen: outer context enumeration overflow \u2014 n_sib=%d >= 30.\n",
          "  This means one node has %d sibling leaves, requiring 2^%d=%s contexts.\n",
          "  Use a deeper hierarchy (hierarchy=NULL = balanced binary tree)."),
          n_sib, n_sib, n_sib, format(2^n_sib, big.mark = ",", scientific = FALSE)))
      for (tmask in 0L:(bitwShiftL(1L, n_sib) - 1L)) {
        bits   <- as.logical(intToBits(tmask))[seq_len(n_sib)]
        T_grp  <- sibling_leaves[which(bits)]
        mask_T <- if (length(T_grp)) .groups_to_mask(as.integer(T_grp), K) else .bitset_empty(K)
        recurse(kids[[ci]]$uid, .bitset_or(outer_group_mask, mask_T))
      }
    }
  }
  recurse(ann$root$uid, .bitset_empty(K))
  contexts
}

# ── §15  Row-loop runners ─────────────────────────────────────────────────────

.run_row_loop_serial <- function(n, compute_one_row,
                                 verbose, dp_progress, dp_print_every, dp_newline_every,
                                 eta_alpha, K, chunk_size_inner, chunk_cap, phase_tag = "DP") {
  ow <- NULL; basev <- numeric(n); vallv <- numeric(n)
  t0 <- proc.time()[[3]]; st <- .make_progress_state()
  done <- 0L
  dp_print_every <- max(1L, as.integer(dp_print_every))
  pb <- NULL
  if (isTRUE(verbose) && isTRUE(dp_progress) && n > 0L) {
    pb <- utils::txtProgressBar(min = 0L, max = n, style = 3, width = 60L, char = "=")
  }
  for (irow in seq_len(n)) {
    rr <- compute_one_row(irow)
    if (is.null(ow)) ow <- matrix(0, n, length(rr$phi))
    ow[irow, ] <- rr$phi; basev[irow] <- rr$base; vallv[irow] <- rr$vall
    done <- done + 1L
    if (!is.null(pb) && (done %% dp_print_every == 0L || done == n))
      utils::setTxtProgressBar(pb, done)
  }
  if (!is.null(pb)) {
    close(pb)
    elapsed_sec <- proc.time()[[3]] - t0
    cat(sprintf("[%s] %d/%d observations | K=%d | elapsed=%s\n",
                phase_tag, n, n, K, .fmt_hms(elapsed_sec)))
  }
  list(ow = ow, basev = basev, vallv = vallv)
}

.run_row_loop_parallel <- function(n, compute_chunk_fork_1, compute_chunk_fork_rest,
                                   compute_chunk_win,
                                   n_cores, verbose,
                                   phase_tag = "DP",
                                   cpp_available = FALSE) {
  chunk_ids <- vector("list", n_cores)
  base_size <- n %/% n_cores; remainder <- n %% n_cores; cur <- 1L
  for (ci in seq_len(n_cores)) {
    sz <- base_size + if (ci <= remainder) 1L else 0L
    if (sz > 0L) chunk_ids[[ci]] <- seq(cur, cur + sz - 1L)
    cur <- cur + sz
  }
  chunk_ids <- Filter(Negate(is.null), chunk_ids)
  n_chunks  <- length(chunk_ids)
  t0 <- proc.time()[[3]]
  if (isTRUE(verbose))
    cat(sprintf("[%s] Parallel: %d obs | %d cores | ~%d obs/core\n",
                phase_tag, n, n_chunks, ceiling(n / n_chunks)))
  use_mclapply <- (.Platform$OS.type != "windows") && .pkg_ok("parallel")
  if (use_mclapply) {
    workers <- c(list(compute_chunk_fork_1),
                 rep(list(compute_chunk_fork_rest), max(0L, n_chunks - 1L)))
    chunk_results <- parallel::mclapply(
      seq_len(n_chunks),
      function(ci) workers[[ci]](chunk_ids[[ci]]),
      mc.cores    = min(n_cores, n_chunks),
      mc.set.seed = TRUE
    )
    failed <- vapply(chunk_results, inherits, logical(1L), "try-error")
    if (any(failed)) {
      warning(sprintf("[TreeOwen] %d/%d parallel chunks failed; retrying serially.",
                      sum(failed), n_chunks))
      for (ci in which(failed))
        chunk_results[[ci]] <- compute_chunk_fork_rest(chunk_ids[[ci]])
    }
  } else if (.pkg_ok("foreach") && .pkg_ok("doParallel")) {
    cl <- parallel::makeCluster(n_cores); doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() }, add = TRUE)
    `%dopar%` <- foreach::`%dopar%`
    chunk_results <- foreach::foreach(irows = chunk_ids) %dopar% compute_chunk_win(irows)
  } else {
    if (isTRUE(verbose)) message("[TreeOwen] parallel packages unavailable; using serial.")
    chunk_results <- lapply(chunk_ids, compute_chunk_fork_rest)
  }
  if (isTRUE(verbose)) {
    elapsed_sec <- proc.time()[[3]] - t0
    cat(sprintf("[%s] %d/%d observations complete | %d cores | elapsed=%s\n",
                phase_tag, n, n, n_chunks, .fmt_hms(elapsed_sec)))
  }
  ow <- NULL; basev <- numeric(n); vallv <- numeric(n)
  for (ci in seq_along(chunk_ids)) {
    chunk_rr <- chunk_results[[ci]]
    for (li in seq_along(chunk_ids[[ci]])) {
      irow <- chunk_ids[[ci]][[li]]; rr <- chunk_rr[[li]]
      if (is.null(ow)) ow <- matrix(0, n, length(rr$phi))
      ow[irow, ] <- rr$phi; basev[irow] <- rr$base; vallv[irow] <- rr$vall
    }
  }
  list(ow = ow, basev = basev, vallv = vallv)
}

.run_row_loop <- function(n,
                          compute_one_row_serial,
                          compute_chunk_fork_1    = NULL,
                          compute_chunk_fork_rest = NULL,
                          compute_chunk_win       = NULL,
                          n_cores,
                          verbose, dp_progress, dp_print_every, dp_newline_every,
                          eta_alpha, K, chunk_size_inner, chunk_cap,
                          phase_tag = "DP", cpp_available = FALSE) {
  n_cores <- as.integer(n_cores); if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  if (n_cores == 1L)
    .run_row_loop_serial(n, compute_one_row_serial, verbose, dp_progress, dp_print_every,
                         dp_newline_every, eta_alpha, K, chunk_size_inner, chunk_cap, phase_tag)
  else
    .run_row_loop_parallel(n,
                           compute_chunk_fork_1    = compute_chunk_fork_1,
                           compute_chunk_fork_rest = compute_chunk_fork_rest,
                           compute_chunk_win       = compute_chunk_win,
                           n_cores = n_cores, verbose = verbose,
                           phase_tag = phase_tag, cpp_available = cpp_available)
}

# ── §16  Result constructor ───────────────────────────────────────────────────

.make_treeowen_result <- function(ow, basev, vallv, feat_names,
                                  unified_model, x, groups, stats, note) {
  rownames(ow) <- rownames(x); colnames(ow) <- feat_names
  res <- list(owens = as.data.frame(ow), baseline = basev, v_all = vallv,
              unified_model = unified_model, observations = x, groups = groups,
              note = note, stats = stats)
  class(res) <- c("treeowen_result", "list"); res
}

# ── §17  Miscellaneous ────────────────────────────────────────────────────────

.max_degree_hier_tree <- function(node) {
  if (identical(node$type, "leaf")) return(0L)
  d <- length(node$children); mx <- d
  for (ch in node$children) mx <- max(mx, .max_degree_hier_tree(ch))
  as.integer(mx)
}


# ══════════════════════════════════════════════════════════════════════════════
# §19  XGBoost / treeshap compatibility wrapper
# ══════════════════════════════════════════════════════════════════════════════
#
# Known failure modes handled (validated across xgboost 1.6 – 2.x, treeshap 0.3–0.4):
#
#  [XGB-1]  treeshap::xgboost.unify() calls model$getinfo() in ways that break
#           on newer xgboost class structures → bypass treeshap entirely.
#  [XGB-2]  xgb.model.dt.tree() omits the "ID" column in older xgboost (< 1.7).
#  [XGB-3]  Split-gain column named "Quality" (xgboost >= 1.7) instead of "Gain".
#  [XGB-4]  "Yes"/"No"/"Missing" columns are stored as node-ID strings
#           ("0-0", "0-1", …) that must be converted to row indices.
#  [XGB-5]  model$feature_names may be NULL when the model was built from a
#           plain matrix without colnames; fall back to colnames(data).
#  [XGB-6]  Leaf nodes in xgboost 2.x may have Gain = 0 / Cover = 0; keep as-is.
#  [XGB-7]  xgb.model.dt.tree() signature changed in xgboost 2.x:
#           argument "model" was renamed to "xgb_model" in some builds.
#
# Users should call this wrapper instead of treeshap::xgboost.unify():
#   unified <- xgboost_unify_compat(xgb_mod, X)
# ──────────────────────────────────────────────────────────────────────────────
#' Unify an XGBoost model for use with treeowen
#'
#' A robust replacement for \code{treeshap::xgboost.unify()} that handles
#' all known version-compatibility issues across xgboost 1.6 to 2.x and
#' treeshap 0.3 to 0.4. Builds the unified model object directly from
#' \code{xgb.model.dt.tree()}, bypassing treeshap internals entirely.
#'
#' @param model A trained \code{xgb.Booster} object.
#' @param data A data frame or matrix of the training features (no outcome
#'   column). Used to resolve feature names when \code{model$feature_names}
#'   is NULL.
#' @param recalculate Logical. If \code{TRUE} and treeshap is available,
#'   calls \code{treeshap::set_reference_dataset()} after unification.
#'   Default \code{FALSE}.
#' @return A \code{model_unified} object suitable for \code{\link{treeowen}}.
#' @export
xgboost_unify_compat <- function(model, data, recalculate = FALSE) {
  if (!requireNamespace("xgboost",    quietly = TRUE))
    stop('Package "xgboost" needed. Please install it.', call. = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE))
    stop('Package "data.table" needed. Please install it.', call. = FALSE)

  # [XGB-7] Try both argument names for xgb.model.dt.tree()
  xgbtree <- tryCatch(
    xgboost::xgb.model.dt.tree(model = model),
    error = function(e1) tryCatch(
      xgboost::xgb.model.dt.tree(xgb_model = model),
      error = function(e2) stop(sprintf(
        "xgb.model.dt.tree() failed with both 'model' and 'xgb_model' arguments.\n  Error 1: %s\n  Error 2: %s",
        conditionMessage(e1), conditionMessage(e2)
      ))
    )
  )
  xgbtree <- data.table::as.data.table(xgbtree)

  # [XGB-2] Reconstruct "ID" when missing
  if (!"ID" %in% colnames(xgbtree))
    xgbtree[, ID := paste(Tree, Node, sep = "-")]

  # [XGB-3] Rename "Quality" → "Gain" when needed
  if ("Quality" %in% colnames(xgbtree) && !"Gain" %in% colnames(xgbtree))
    data.table::setnames(xgbtree, "Quality", "Gain")

  # Verify required columns
  needed <- c("Tree", "Node", "ID", "Feature", "Split", "Yes", "No", "Missing", "Gain", "Cover")
  miss <- setdiff(needed, colnames(xgbtree))
  if (length(miss) > 0)
    stop("xgb.model.dt.tree() is missing required columns: ", paste(miss, collapse = ", "))

  # [XGB-4] Convert node-ID strings to row indices
  id_str  <- as.character(xgbtree$ID)
  xgbtree[, Yes     := match(as.character(Yes),     id_str)]
  xgbtree[, No      := match(as.character(No),      id_str)]
  xgbtree[, Missing := match(as.character(Missing), id_str)]

  # Mark leaf nodes
  xgbtree[is.na(Split), Feature := NA_character_]

  xgbtree[, Decision.type := factor(rep("<=", .N), levels = c("<=", "<"))]
  xgbtree[is.na(Feature), Decision.type := NA]

  # Select & rename to treeshap's expected layout
  xgbtree <- xgbtree[, .(Tree, Node, Feature, Decision.type, Split,
                          Yes, No, Missing, Gain, Cover)]
  setnames(xgbtree, "Gain", "Prediction")

  # [XGB-6] Internal nodes: Prediction field not used → NA
  xgbtree[!is.na(Feature), Prediction := NA_real_]

  # [XGB-5] Feature names
  feature_names <- model$feature_names
  if (is.null(feature_names) || !length(feature_names))
    feature_names <- colnames(data)
  if (is.null(feature_names) || !length(feature_names))
    stop("Unable to determine XGBoost feature names from model or data.")

  if (is.matrix(data)) data <- as.data.frame(data)
  data <- data[, intersect(colnames(data), feature_names), drop = FALSE]
  data <- data[, feature_names[feature_names %in% colnames(data)], drop = FALSE]

  ret <- list(
    model         = as.data.frame(xgbtree),
    data          = as.data.frame(data),
    feature_names = feature_names
  )
  class(ret) <- "model_unified"
  attr(ret, "missing_support") <- TRUE
  attr(ret, "model")           <- "xgboost"

  if (recalculate && requireNamespace("treeshap", quietly = TRUE))
    ret <- tryCatch(
      treeshap::set_reference_dataset(ret, as.data.frame(data)),
      error = function(e) ret
    )
  ret
}


# ══════════════════════════════════════════════════════════════════════════════
# §20  LightGBM / treeshap compatibility wrapper
# ══════════════════════════════════════════════════════════════════════════════
#
# Known failure modes handled (validated across lightgbm 3.x – 4.x, treeshap 0.3–0.4):
#
#  [LGB-1]  treeshap::lightgbm.unify() internals break on lightgbm >= 4.0 due to
#           Booster class restructuring and lgb.model.dt.tree() output changes.
#  [LGB-2]  lgb.model.dt.tree() may be unavailable in very old lightgbm (< 3.2);
#           fall back to model$dump_model() JSON parsing.
#  [LGB-3]  Column name changes across versions:
#             "split_gain"  ↔  "value"         (gain / leaf value)
#             "threshold"   ↔  "split_point"   (split threshold)
#             "count"       ↔  "internal_count" / "leaf_count"  (cover)
#             "num_cat"     may be absent
#             "decision_type" may be absent or empty
#  [LGB-4]  Child node IDs use three different encodings depending on version:
#             (a) negative-1-indexed leaves:   left_child < 0  → leaf #(-left_child-1)
#             (b) positive row indices (rare older builds)
#             (c) string "N<int>" / "L<int>" encoding
#  [LGB-5]  Feature names may be NULL for some training configurations;
#           fall back to model$feature_name() → colnames(data).
#  [LGB-6]  lgb.cv() / lgb.train() best_iter slot name changed:
#           cv$best_iter  (older)  ↔  cv$best_iteration  (newer).
#  [LGB-7]  lgb.cv() verbose argument moved out of params in newer lightgbm.
#
# Strategy: try treeshap::lightgbm.unify() first; if it fails, build
# model_unified directly from lgb.model.dt.tree() output; if that also fails,
# fall back to model$dump_model() JSON tree parsing.
# ──────────────────────────────────────────────────────────────────────────────
#' Unify a LightGBM model for use with treeowen
#'
#' A robust replacement for \code{treeshap::lightgbm.unify()} that handles
#' all known version-compatibility issues across lightgbm 3.x / 4.x and
#' treeshap 0.3 to 0.4. Falls back through three strategies: treeshap,
#' \code{lgb.model.dt.tree()}, and \code{model$dump_model()} JSON parsing.
#'
#' @param model A trained \code{lgb.Booster} object.
#' @param data A data frame or matrix of the training features.
#' @param recalculate Logical. Default \code{FALSE}.
#' @return A \code{model_unified} object suitable for \code{\link{treeowen}}.
#' @export
lightgbm_unify_compat <- function(model, data, recalculate = FALSE) {
  if (!requireNamespace("lightgbm", quietly = TRUE))
    stop('Package "lightgbm" needed. Please install it.', call. = FALSE)

  if (is.matrix(data)) data <- as.data.frame(data)

  # ── Path 1: treeshap::lightgbm.unify() ──────────────────────────────────────
  if (requireNamespace("treeshap", quietly = TRUE)) {
    ret <- tryCatch(treeshap::lightgbm.unify(model, data), error = function(e) NULL)
    if (!is.null(ret) && inherits(ret, "model_unified")) {
      if (recalculate)
        ret <- tryCatch(
          treeshap::set_reference_dataset(ret, as.data.frame(data)),
          error = function(e) ret)
      return(ret)
    }
  }

  # ── Helper: resolve feature names ───────────────────────────────────────────
  .lgb_feature_names <- function(model, data) {
    nm <- tryCatch(model$.__enclos_env__$private$valid_feature_names,
                   error = function(e) NULL)
    if (!is.null(nm) && length(nm)) return(nm)
    nm <- tryCatch(model$feature_name(), error = function(e) NULL)
    if (!is.null(nm) && length(nm)) return(nm)
    nm <- colnames(data)
    if (!is.null(nm) && length(nm)) return(nm)
    stop("lightgbm: unable to determine feature names from model or data.")
  }
  feature_names <- .lgb_feature_names(model, data)

  # ── Path 2: lgb.model.dt.tree() ─────────────────────────────────────────────
  lgbtree <- tryCatch({
    dt <- lightgbm::lgb.model.dt.tree(model)
    if (is.null(dt) || nrow(dt) == 0L) stop("empty")
    as.data.frame(dt, stringsAsFactors = FALSE)
  }, error = function(e) NULL)

  # ── Path 3: dump_model() JSON fallback ──────────────────────────────────────
  if (is.null(lgbtree)) {
    lgbtree <- tryCatch({
      dm  <- model$dump_model()
      ti  <- dm[["tree_info"]]
      if (is.null(ti)) stop("tree_info absent")
      lgb_parse_dump_model(ti)
    }, error = function(e) {
      stop(sprintf(
        "lightgbm: all tree-dump methods failed.\n  lgb.model.dt.tree(): unavailable\n  dump_model(): %s",
        conditionMessage(e)))
    })
  }

  lgbtree <- as.data.frame(lgbtree, stringsAsFactors = FALSE)
  nms     <- colnames(lgbtree)

  # ── [LGB-3] Normalise column names ──────────────────────────────────────────
  .alias <- function(df, pref, alts) {
    if (pref %in% colnames(df)) return(df)
    for (a in alts) if (a %in% colnames(df)) { df[[pref]] <- df[[a]]; return(df) }
    df[[pref]] <- NA_real_; df
  }
  lgbtree <- .alias(lgbtree, "threshold",         c("split_point"))
  lgbtree <- .alias(lgbtree, "split_gain",        c("value"))
  lgbtree <- .alias(lgbtree, "count",             c("internal_count", "leaf_count", "node_count"))
  lgbtree <- .alias(lgbtree, "leaf_value",        c("split_gain", "value"))
  lgbtree <- .alias(lgbtree, "decision_type",     character(0)); if (is.na(lgbtree$decision_type[1])) lgbtree$decision_type <- "<="
  lgbtree <- .alias(lgbtree, "num_cat",           character(0)); if (is.na(lgbtree$num_cat[1]))       lgbtree$num_cat       <- 0L
  lgbtree <- .alias(lgbtree, "missing_direction", c("missing_type"))

  # Ensure tree_index column
  if (!"tree_index" %in% colnames(lgbtree)) {
    if ("Tree" %in% colnames(lgbtree))  lgbtree$tree_index <- lgbtree$Tree
    else stop("lightgbm: cannot find tree_index column in tree dump")
  }

  n_rows  <- nrow(lgbtree)
  tree_id <- as.integer(lgbtree$tree_index)

  # ── [LGB-4] Build child row-pointer lookup ───────────────────────────────────
  # Identify leaves: no left/right child columns → use split_feature / feature to detect
  is_leaf <- rep(FALSE, n_rows)
  if ("left_child" %in% colnames(lgbtree)) {
    lc_raw <- lgbtree$left_child
    rc_raw <- lgbtree$right_child
    # Leaves: both NA, or both 0 in some encodings
    is_leaf <- is.na(lc_raw) | (suppressWarnings(as.integer(lc_raw)) == 0L &
                                  is.na(lgbtree$threshold))
  } else {
    # Infer: rows without a feature split are leaves
    feat_raw <- if ("feature" %in% colnames(lgbtree)) lgbtree$feature
                else if ("split_feature" %in% colnames(lgbtree)) lgbtree$split_feature
                else rep(NA_character_, n_rows)
    is_leaf <- is.na(feat_raw) | feat_raw == ""
  }

  # Encode each row with its tree-local node/leaf ID for lookup
  # (See [LGB-4] for the three encoding schemes)
  has_node_idx <- "node_index"  %in% colnames(lgbtree)
  has_leaf_idx <- "leaf_index"  %in% colnames(lgbtree)

  if (has_node_idx || has_leaf_idx) {
    node_idx <- if (has_node_idx) suppressWarnings(as.integer(lgbtree$node_index)) else rep(NA_integer_, n_rows)
    leaf_idx <- if (has_leaf_idx) suppressWarnings(as.integer(lgbtree$leaf_index)) else rep(NA_integer_, n_rows)
    encoded_id <- ifelse(is_leaf, -(leaf_idx + 1L), node_idx)
  } else {
    # Reconstruct sequential IDs within each tree
    node_ctr <- integer(max(tree_id, na.rm = TRUE) + 2L)
    leaf_ctr <- integer(max(tree_id, na.rm = TRUE) + 2L)
    encoded_id <- integer(n_rows)
    for (r in seq_len(n_rows)) {
      tid <- tree_id[r] + 1L
      if (is_leaf[r]) {
        encoded_id[r] <- -(leaf_ctr[tid])
        leaf_ctr[tid] <- leaf_ctr[tid] + 1L
      } else {
        encoded_id[r] <- node_ctr[tid]
        node_ctr[tid] <- node_ctr[tid] + 1L
      }
    }
  }

  lookup_key <- paste(tree_id, encoded_id, sep = ":")
  lookup_map <- setNames(seq_len(n_rows), lookup_key)

  # Parse a single child column → row indices
  .resolve_lgb_child <- function(raw_col) {
    if (is.null(raw_col)) return(rep(NA_integer_, n_rows))
    raw_char <- as.character(raw_col)
    encoded  <- suppressWarnings(as.integer(raw_char))
    # [LGB-4c] "N<int>" / "L<int>" string encoding
    n_mask <- grepl("^N", raw_char, perl = TRUE)
    l_mask <- grepl("^L", raw_char, perl = TRUE)
    encoded[n_mask] <-  as.integer(sub("^N", "", raw_char[n_mask]))
    encoded[l_mask] <- -(as.integer(sub("^L", "", raw_char[l_mask])) + 1L)
    key <- paste(tree_id, encoded, sep = ":")
    as.integer(unname(lookup_map[key]))
  }

  yes_rows  <- .resolve_lgb_child(lgbtree[["left_child"]])
  no_rows   <- .resolve_lgb_child(lgbtree[["right_child"]])

  # Missing direction
  miss_rows <- rep(NA_integer_, n_rows)
  if ("missing_direction" %in% colnames(lgbtree)) {
    md        <- tolower(as.character(lgbtree$missing_direction))
    miss_rows <- ifelse(md == "left",  yes_rows,
                 ifelse(md == "right", no_rows, NA_integer_))
  }

  # ── Feature column ───────────────────────────────────────────────────────────
  feat_raw <- if ("feature" %in% colnames(lgbtree))       lgbtree$feature
              else if ("split_feature" %in% colnames(lgbtree)) lgbtree$split_feature
              else rep(NA_character_, n_rows)
  feat_col <- as.character(feat_raw)
  feat_col[is_leaf] <- NA_character_

  # Convert 0-based numeric index → feature name
  idx_mask <- !is.na(feat_col) & !is.na(suppressWarnings(as.integer(feat_col)))
  if (any(idx_mask)) {
    idx0 <- suppressWarnings(as.integer(feat_col[idx_mask]))
    feat_col[idx_mask] <- ifelse(
      idx0 + 1L <= length(feature_names),
      feature_names[idx0 + 1L],
      NA_character_)
  }

  # ── Prediction (leaf value) ──────────────────────────────────────────────────
  pred_val <- rep(NA_real_, n_rows)
  for (lv_col in c("leaf_value", "split_gain", "value")) {
    if (lv_col %in% colnames(lgbtree)) {
      vals <- suppressWarnings(as.numeric(lgbtree[[lv_col]]))
      pred_val[is_leaf] <- vals[is_leaf]
      break
    }
  }

  # ── Other columns ────────────────────────────────────────────────────────────
  split_val <- suppressWarnings(as.numeric(lgbtree$threshold))
  split_val[is_leaf] <- NA_real_

  cover_val <- suppressWarnings(as.numeric(lgbtree$count))
  if (all(is.na(cover_val))) cover_val <- rep(1, n_rows)

  dec_type <- factor(rep("<=", n_rows), levels = c("<=", "<"))
  dec_type[is_leaf] <- NA

  model_df <- data.frame(
    Tree          = tree_id,
    Node          = as.integer(lgbtree$node_depth),
    Feature       = feat_col,
    Decision.type = dec_type,
    Split         = split_val,
    Yes           = yes_rows,
    No            = no_rows,
    Missing       = miss_rows,
    Prediction    = pred_val,
    Cover         = cover_val,
    stringsAsFactors = FALSE
  )

  data <- data[, intersect(colnames(data), feature_names), drop = FALSE]
  data <- data[, feature_names[feature_names %in% colnames(data)], drop = FALSE]

  ret <- list(
    model         = model_df,
    data          = as.data.frame(data),
    feature_names = feature_names
  )
  class(ret) <- "model_unified"
  attr(ret, "missing_support") <- "missing_direction" %in% colnames(lgbtree)
  attr(ret, "model")           <- "lightgbm"

  if (recalculate && requireNamespace("treeshap", quietly = TRUE))
    ret <- tryCatch(
      treeshap::set_reference_dataset(ret, as.data.frame(data)),
      error = function(e) ret)
  ret
}

# Internal: parse dump_model()$tree_info into a flat data.frame  [LGB-2 fallback]
lgb_parse_dump_model <- function(trees_info) {
  rows <- list()
  parse_node <- function(node, tree_id) {
    if (is.null(node)) return(invisible(NULL))
    is_leaf <- !is.null(node[["leaf_value"]]) || is.null(node[["split_feature"]])
    rows[[length(rows) + 1L]] <<- list(
      tree_index        = tree_id,
      node_depth        = if (!is.null(node$depth))          node$depth          else NA_integer_,
      node_index        = if (!is.null(node$split_index))    node$split_index    else NA_integer_,
      leaf_index        = if (!is.null(node$leaf_index))     node$leaf_index     else NA_integer_,
      feature           = if (!is_leaf) as.character(node$split_feature)         else NA_character_,
      threshold         = if (!is_leaf) suppressWarnings(as.numeric(node$threshold)) else NA_real_,
      split_gain        = suppressWarnings(as.numeric(node$split_gain)),
      leaf_value        = if (is_leaf)  suppressWarnings(as.numeric(node$leaf_value)) else NA_real_,
      count             = suppressWarnings(as.numeric(
                            node$internal_count %||% node$leaf_count %||% NA_real_)),
      missing_direction = as.character(node$missing_direction %||% NA_character_),
      decision_type     = as.character(node$decision_type     %||% "<="),
      left_child        = NA_integer_,
      right_child       = NA_integer_
    )
    if (!is_leaf) {
      parse_node(node$left_child,  tree_id)
      parse_node(node$right_child, tree_id)
    }
  }
  # NULL-coalescing helper (local scope)
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  for (ti in seq_along(trees_info)) {
    tree <- trees_info[[ti]]
    root <- tree$tree_structure %||% tree
    parse_node(root, tree_id = ti - 1L)
  }
  if (!length(rows)) stop("no rows parsed from dump_model() tree_info")
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}


# ══════════════════════════════════════════════════════════════════════════════
# §21  Ranger / treeshap compatibility wrapper
# ══════════════════════════════════════════════════════════════════════════════
#
# Known failure modes handled (validated across ranger 0.11 – 0.16+, treeshap 0.3–0.4):
#
#  [RNG-1]  treeshap::ranger.unify() assumes ranger::treeInfo() returns a "pred.1"
#           column, but ranger >= 0.14 / 0.16 returns "prediction" instead.
#  [RNG-2]  treeshap internals call ranger_unify.common() with a hard-coded
#           column list; the function signature and expected input changed.
#  [RNG-3]  ranger::treeInfo() column name variations across versions:
#             "nodeID"     / "node"      (node identifier)
#             "splitvarID" / "splitVarID"(split variable 0-based index)
#             "splitval"   / "splitVal"  (split threshold)
#             "terminal"   / "isTerminal"(leaf flag)
#             "prediction" / "pred.1"    (leaf probability, class 1)
#           All are normalised before use.
#  [RNG-4]  child.nodeIDs structure changed:
#             Older ranger: list of length-2 vectors [left_id, right_id] per node.
#             Newer ranger: separate columns "leftChild" / "rightChild" in treeInfo().
#  [RNG-5]  Probability forests (probability=TRUE) store leaf predictions as
#           a matrix (n_nodes × n_classes); the column for class "1" (second level)
#           is extracted.  Classification forests give integer class IDs.
#  [RNG-6]  model$forest may be NULL when write.forest=FALSE was used at training.
#  [RNG-7]  keep.inbag may be irrelevant for treeInfo()-based extraction; the
#           wrapper works regardless.
#  [RNG-8]  treeshap >= 0.4 added a check on probability=TRUE that raises an error
#           for non-probability forests; we warn instead of stopping.
#
# Strategy:
#   1. Try treeshap::ranger.unify() directly.
#   2. If it fails, patch treeInfo() column names and call ranger_unify.common()
#      (the internal treeshap function) directly.
#   3. If that also fails, build model_unified from forest internals manually.
# ──────────────────────────────────────────────────────────────────────────────
#' Unify a Ranger model for use with treeowen
#'
#' A robust replacement for \code{treeshap::ranger.unify()} that handles
#' all known version-compatibility issues across ranger 0.11 to 0.16+ and
#' treeshap 0.3 to 0.4. Falls back through three strategies: treeshap,
#' patched \code{treeInfo()} with \code{ranger_unify.common()}, and direct
#' \code{model$forest} parsing.
#'
#' The model must be trained with \code{probability = TRUE}.
#'
#' @param model A trained \code{ranger} object (\code{probability = TRUE}).
#' @param data A data frame or matrix of the training features (no outcome
#'   column).
#' @param recalculate Logical. Default \code{FALSE}.
#' @return A \code{model_unified} object suitable for \code{\link{treeowen}}.
#' @export
ranger_unify_compat <- function(model, data, recalculate = FALSE) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop('Package "ranger" needed. Please install it.', call. = FALSE)

  if (is.matrix(data)) data <- as.data.frame(data)

  # [RNG-8] Warn for non-probability forest
  is_prob <- isTRUE(model$treetype %in% c("Probability estimation", "probability"))
  if (!is_prob)
    warning(paste0(
      "[ranger_unify_compat] ranger model was not trained with probability=TRUE. ",
      "Leaf predictions will be class labels, not probabilities. ",
      "Re-train with probability=TRUE and keep.inbag=TRUE for reliable results."
    ), call. = FALSE)

  # Feature names from forest internals
  feature_names <- model$forest$independent.variable.names
  if (is.null(feature_names) || !length(feature_names))
    feature_names <- setdiff(colnames(data), ".y")
  if (is.null(feature_names) || !length(feature_names))
    feature_names <- colnames(data)

  # Clean data: remove outcome column if present
  X_clean <- data[, intersect(colnames(data), feature_names), drop = FALSE]
  X_clean <- X_clean[, feature_names[feature_names %in% colnames(X_clean)], drop = FALSE]

  # ── Path 1: standard treeshap::ranger.unify() ───────────────────────────────
  if (requireNamespace("treeshap", quietly = TRUE)) {
    ret <- tryCatch(treeshap::ranger.unify(model, X_clean), error = function(e) NULL)
    if (!is.null(ret) && inherits(ret, "model_unified")) {
      if (recalculate)
        ret <- tryCatch(
          treeshap::set_reference_dataset(ret, as.data.frame(X_clean)),
          error = function(e) ret)
      return(ret)
    }
  }

  # ── Path 2: patched treeInfo() + ranger_unify.common() ──────────────────────
  if (requireNamespace("treeshap",   quietly = TRUE) &&
      requireNamespace("data.table", quietly = TRUE) &&
      requireNamespace("ranger",     quietly = TRUE)) {

    ret2 <- tryCatch({
      n_trees <- model$num.trees
      tree_list <- lapply(seq_len(n_trees), function(ti) {
        td <- data.table::as.data.table(ranger::treeInfo(model, tree = ti))
        # [RNG-3] Normalise column names to what ranger_unify.common() expects
        # nodeID
        for (cand in c("node", "Node")) if (!"nodeID" %in% names(td) && cand %in% names(td))
          data.table::setnames(td, cand, "nodeID")
        # leftChild / rightChild
        for (cand in c("leftChildID", "left_child", "left")) if (!"leftChild" %in% names(td) && cand %in% names(td))
          data.table::setnames(td, cand, "leftChild")
        for (cand in c("rightChildID", "right_child", "right")) if (!"rightChild" %in% names(td) && cand %in% names(td))
          data.table::setnames(td, cand, "rightChild")
        # splitvarName
        for (cand in c("splitVarName", "split_var_name", "splitvar")) if (!"splitvarName" %in% names(td) && cand %in% names(td))
          data.table::setnames(td, cand, "splitvarName")
        # splitval
        for (cand in c("splitVal", "split_val", "split_value", "threshold")) if (!"splitval" %in% names(td) && cand %in% names(td))
          data.table::setnames(td, cand, "splitval")
        # [RNG-1] prediction column: "pred.1" → "prediction"
        if ("pred.1" %in% names(td) && !"prediction" %in% names(td))
          data.table::setnames(td, "pred.1", "prediction")
        # If still no "prediction", find the last numeric column as fallback
        if (!"prediction" %in% names(td)) {
          num_cols <- names(td)[vapply(td, is.numeric, logical(1L))]
          if (length(num_cols))
            data.table::setnames(td, tail(num_cols, 1L), "prediction")
        }
        need <- c("nodeID", "leftChild", "rightChild", "splitvarName", "splitval", "prediction")
        td[, intersect(need, names(td)), with = FALSE]
      })

      # [RNG-2] Locate ranger_unify.common regardless of export status
      common_fn <- tryCatch(
        getFromNamespace("ranger_unify.common", "treeshap"),
        error = function(e) NULL
      )
      if (is.null(common_fn)) stop("ranger_unify.common not found in treeshap")
      common_fn(x = tree_list, n = model$num.trees,
                data = X_clean, feature_names = feature_names)
    }, error = function(e) NULL)

    if (!is.null(ret2) && inherits(ret2, "model_unified")) {
      if (recalculate)
        ret2 <- tryCatch(
          treeshap::set_reference_dataset(ret2, as.data.frame(X_clean)),
          error = function(e) ret2)
      return(ret2)
    }
  }

  # ── Path 3: manual build from forest internals ───────────────────────────────
  # Try treeInfo() first for each tree; fall back to model$forest list structure
  n_trees   <- model$num.trees
  all_rows  <- vector("list", n_trees)

  for (ti in seq_len(n_trees)) {
    ti_df <- tryCatch(
      ranger_treeinfo_norm(ranger::treeInfo(model, tree = ti), feature_names),
      error = function(e) NULL
    )
    if (is.null(ti_df)) {
      ti_df <- tryCatch(
        ranger_tree_from_forest(model, ti, feature_names),
        error = function(e2) stop(sprintf(
          "[ranger_unify_compat] tree %d: treeInfo() and forest fallback both failed.\n  treeInfo error: %s",
          ti, conditionMessage(e2)))
      )
    }
    all_rows[[ti]] <- data.frame(
      Tree          = ti - 1L,
      Node          = ti_df$nodeid,
      Feature       = ti_df$feature,
      Decision.type = {
        dt <- factor(rep("<=", nrow(ti_df)), levels = c("<=", "<"))
        dt[ti_df$terminal] <- NA
        dt
      },
      Split         = ti_df$splitval,
      Yes           = ti_df$leftchild,
      No            = ti_df$rightchild,
      Missing       = NA_integer_,
      Prediction    = ti_df$prediction,
      Cover         = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  model_df <- do.call(rbind, all_rows)
  rownames(model_df) <- NULL

  # Convert tree-local row indices to global row indices
  row_offsets <- c(0L, cumsum(vapply(all_rows, nrow, integer(1L))))
  for (ti in seq_len(n_trees)) {
    rng <- seq(row_offsets[ti] + 1L, row_offsets[ti + 1L])
    off <- row_offsets[ti]
    model_df$Yes[rng] <- ifelse(is.na(model_df$Yes[rng]), NA_integer_,
                                 model_df$Yes[rng] + off)
    model_df$No[rng]  <- ifelse(is.na(model_df$No[rng]),  NA_integer_,
                                 model_df$No[rng]  + off)
  }

  ret3 <- list(
    model         = model_df,
    data          = as.data.frame(X_clean),
    feature_names = feature_names
  )
  class(ret3) <- "model_unified"
  attr(ret3, "missing_support") <- FALSE
  attr(ret3, "model")           <- "ranger"

  if (recalculate && requireNamespace("treeshap", quietly = TRUE))
    ret3 <- tryCatch(
      treeshap::set_reference_dataset(ret3, as.data.frame(X_clean)),
      error = function(e) ret3)
  ret3
}

# Internal: normalise a single treeInfo() data.frame to standard columns
# Returns: data.frame with cols: nodeid, feature, splitval, terminal, leftchild, rightchild, prediction
ranger_treeinfo_norm <- function(df, feature_names) {
  df   <- as.data.frame(df, stringsAsFactors = FALSE)
  nms  <- tolower(gsub("[._]", "", colnames(df)))   # flatten: nodeID → nodeid etc.
  colnames(df) <- nms

  # nodeID
  for (c in c("nodeid", "node"))             if (!"nodeid"    %in% nms && c %in% nms) { df$nodeid    <- df[[c]]; break }
  if (!"nodeid" %in% colnames(df))            df$nodeid    <- seq_len(nrow(df)) - 1L

  # splitvarName / splitvarID → feature name
  feat_col <- NULL
  for (c in c("splitvarname", "splitvariablename", "splitvar")) if (c %in% nms) { feat_col <- df[[c]]; break }
  if (is.null(feat_col)) {
    for (c in c("splitvarid", "splitvariableid", "varid")) {
      if (c %in% nms) {
        ids <- suppressWarnings(as.integer(df[[c]]))
        feat_col <- ifelse(!is.na(ids) & ids >= 0L & ids + 1L <= length(feature_names),
                          feature_names[ids + 1L], NA_character_)
        break
      }
    }
  }
  if (is.null(feat_col)) feat_col <- rep(NA_character_, nrow(df))

  # terminal flag
  term <- NULL
  for (c in c("terminal", "isterminal", "isleaf", "leaf")) if (c %in% nms) { term <- as.logical(df[[c]]); break }
  if (is.null(term)) term <- is.na(feat_col) | feat_col == ""
  feat_col[term] <- NA_character_

  # splitval
  sv <- NULL
  for (c in c("splitval", "splitvalue", "splitpoint", "threshold", "split")) if (c %in% nms) { sv <- suppressWarnings(as.numeric(df[[c]])); break }
  if (is.null(sv)) sv <- rep(NA_real_, nrow(df))
  sv[term] <- NA_real_

  # leftChild / rightChild (0-based node IDs → 1-based row indices)
  lc <- rc <- NULL
  for (c in c("leftchild", "leftchildid", "leftchild")) if (c %in% nms) { lc <- suppressWarnings(as.integer(df[[c]])); break }
  for (c in c("rightchild", "rightchildid", "rightchild")) if (c %in% nms) { rc <- suppressWarnings(as.integer(df[[c]])); break }
  if (is.null(lc) || is.null(rc)) {
    for (c in c("childnodeids", "child.nodeids")) {
      if (c %in% nms) {
        cn <- df[[c]]
        lc <- vapply(cn, function(x) as.integer(x[[1L]]), integer(1L))
        rc <- vapply(cn, function(x) as.integer(x[[2L]]), integer(1L))
        break
      }
    }
  }
  if (is.null(lc)) lc <- rc <- rep(NA_integer_, nrow(df))

  node_ids <- as.integer(df$nodeid)
  n2r <- setNames(seq_len(nrow(df)), as.character(node_ids))
  resolve_r <- function(ids) {
    ids <- as.integer(ids)
    out <- unname(n2r[as.character(ids)])
    out[is.na(ids) | ids < 0L] <- NA_integer_
    as.integer(out)
  }
  left_r  <- resolve_r(lc)
  right_r <- resolve_r(rc)

  # prediction (leaf probability for class 1 in probability forest)
  pred <- NULL
  for (c in c("prediction", "pred.1", "pred1", "nodeprediction")) if (c %in% nms) { pred <- df[[c]]; break }
  if (is.null(pred)) pred <- rep(NA_real_, nrow(df))
  if (is.list(pred)) {
    pred <- vapply(pred, function(v) {
      v <- as.numeric(v)
      if (length(v) >= 2L) v[2L] else if (length(v) == 1L) v[1L] else NA_real_
    }, numeric(1L))
  }
  if (is.matrix(pred))
    pred <- if (ncol(pred) >= 2L) as.numeric(pred[, 2L]) else as.numeric(pred[, 1L])
  pred <- suppressWarnings(as.numeric(pred))
  pred[!term] <- NA_real_

  data.frame(nodeid     = node_ids,
             feature    = feat_col,
             splitval   = sv,
             terminal   = term,
             leftchild  = left_r,
             rightchild = right_r,
             prediction = pred,
             stringsAsFactors = FALSE)
}

# Internal: build treeInfo-like data.frame from model$forest list
# (fallback when ranger::treeInfo() is unavailable)  [RNG-4 fallback]
ranger_tree_from_forest <- function(model, tree_idx, feature_names) {
  forest <- model$forest
  if (is.null(forest))
    stop("model$forest is NULL; re-train with write.forest=TRUE")

  cn <- forest$child.nodeIDs
  if (is.null(cn)) stop("model$forest$child.nodeIDs is NULL")

  tree_cn <- cn[[tree_idx]]
  if (is.matrix(tree_cn)) {
    left_ids  <- as.integer(tree_cn[1L, ])
    right_ids <- as.integer(tree_cn[2L, ])
  } else if (is.list(tree_cn) && length(tree_cn) == 2L) {
    left_ids  <- as.integer(tree_cn[[1L]])
    right_ids <- as.integer(tree_cn[[2L]])
  } else stop("unrecognised child.nodeIDs structure in tree ", tree_idx)

  n_nodes   <- length(left_ids)
  node_ids  <- seq_len(n_nodes) - 1L
  is_term   <- (left_ids == 0L & right_ids == 0L)

  # Split variable IDs (0-based)
  sv_list <- forest$splitvarIDs %||% forest$split.varIDs
  split_var_ids <- if (!is.null(sv_list)) as.integer(sv_list[[tree_idx]]) else rep(NA_integer_, n_nodes)

  # Split values
  sv_list2  <- forest$splitValues %||% forest$split.values
  split_vals <- if (!is.null(sv_list2)) suppressWarnings(as.numeric(sv_list2[[tree_idx]])) else rep(NA_real_, n_nodes)

  # Feature names (0-based splitvarID → name)
  feat_col <- rep(NA_character_, n_nodes)
  for (ni in which(!is_term)) {
    fid <- split_var_ids[ni]
    if (!is.na(fid) && fid >= 0L && fid + 1L <= length(feature_names))
      feat_col[ni] <- feature_names[fid + 1L]
  }

  # Leaf predictions (probability for class "1")
  pred_vals <- rep(NA_real_, n_nodes)
  tcc <- forest$terminal.class.counts
  if (!is.null(tcc) && length(tcc) >= tree_idx) {
    ct <- tcc[[tree_idx]]
    if (is.matrix(ct) && nrow(ct) == n_nodes) {
      tots <- rowSums(ct); tots[tots == 0] <- 1
      pred_vals[is_term] <- if (ncol(ct) >= 2L) (ct / tots)[is_term, 2L]
                            else                  (ct / tots)[is_term, 1L]
    } else if (is.list(ct)) {
      for (ni in which(is_term)) {
        cv <- as.numeric(ct[[ni]])
        if (!is.null(cv) && length(cv) >= 2L) {
          tot <- sum(cv); if (tot == 0) tot <- 1
          pred_vals[ni] <- cv[2L] / tot
        }
      }
    }
  }

  # Convert 0-based child IDs to 1-based row indices
  node_map <- setNames(seq_len(n_nodes), as.character(node_ids))
  resolve  <- function(ids) {
    out <- unname(node_map[as.character(as.integer(ids))])
    out[is_term | is.na(ids)] <- NA_integer_
    as.integer(out)
  }

  data.frame(nodeid     = node_ids,
             feature    = feat_col,
             splitval   = ifelse(is_term, NA_real_, split_vals),
             terminal   = is_term,
             leftchild  = resolve(left_ids),
             rightchild = resolve(right_ids),
             prediction = pred_vals,
             stringsAsFactors = FALSE)
}

# NULL-coalescing operator (file-local, not exported)
`%||%` <- function(a, b) if (!is.null(a)) a else b
