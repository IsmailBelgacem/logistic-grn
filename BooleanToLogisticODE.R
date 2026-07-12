# ============================================================================
#  BooleanToLogisticODE.R
#  -------------------------------------------------------------------------
#  General De Morgan translator: Boolean network  ->  product-of-logistics ODE.
#  Native-R equivalent of the Mathematica routine BooleanToOdeSystem
#  (package BooleanToLogisticODE.m). It applies the SAME recursive soft
#  De Morgan map, so the two pipelines produce identical regulatory functions.
#
#     literal   v          ->  fpos(v)                  (increasing logistic)
#     literal  !v          ->  fneg(v)   ( = 1 - fpos(v) for logistic & Hill )
#     !f      (compound)   ->  1 - soft(f)              (De Morgan negation)
#     f1 & f2 (AND)        ->  soft(f1) * soft(f2)      (product t-norm)
#     f1 | f2 (OR)         ->  1 - (1-soft(f1))(1-soft(f2))
#                                          (probabilistic sum = De Morgan dual)
#     TRUE / FALSE         ->  1 / 0
#
#  Each gene then obeys     x_i'(t) = kappa_i * Phi_i(x) - gamma_i * x_i(t),
#  with the steepness matched to a Hill coefficient n via  lambda_i = n/theta_i.
#
## ---- logistic kernels (lambda = n / thr) -----------------------------------
logisticp <- function(x, thr, n) 1 / (1 + exp(-(n/thr) * (x - thr)))   # increasing / activation
logisticm <- function(x, thr, n) 1 / (1 + exp( (n/thr) * (x - thr)))   # decreasing / repression

## ---- variables occurring in a quoted Boolean formula -----------------------
ruleVars <- function(expr) {
  v <- character(0)
  rec <- function(e) {
    if (is.name(e)) v[[length(v) + 1]] <<- as.character(e)
    else if (is.call(e)) for (i in 2:length(e)) rec(e[[i]])
  }
  rec(expr); unique(v)
}


minimizeDNF <- function(expr, vars = ruleVars(expr)) {
  ## Quine-McCluskey prime-implicant generation (exact) followed by essential-
  ## prime selection and a GREEDY cover of any minterms left uncovered (lines
  ## below). The greedy step is a standard set-cover heuristic and is not
  ## guaranteed to return the globally smallest possible DNF in every case;
  ## it is, however, always a CORRECT DNF (logically equivalent to `expr`,
  ## verified by construction from the prime-implicant chart) regardless of
  ## whether that particular cover is of minimum cardinality. For every rule
  ## in this paper's networks this has been checked, by direct execution, to
  ## reproduce the intended minimal form (e.g. AKT's rule below).
  k <- length(vars)
  if (k == 0) return(isTRUE(eval(expr, list())))
  grid <- expand.grid(rep(list(c(FALSE, TRUE)), k), KEEP.OUT.ATTRS = FALSE)
  names(grid) <- vars
  onset <- which(vapply(seq_len(nrow(grid)),
                        function(i) isTRUE(eval(expr, as.list(grid[i, , drop = FALSE]))),
                        logical(1)))
  if (length(onset) == 0)          return(FALSE)
  if (length(onset) == nrow(grid)) return(TRUE)
  mt <- lapply(onset, function(i) as.integer(unlist(grid[i, ])))  # 0/1 vector, vars[1] first
  ## iterative merging of adjacent implicants (NA = "dash")
  imps <- mt; primes <- list()
  repeat {
    n <- length(imps); used <- rep(FALSE, n); merged <- list()
    if (n > 1) for (a in 1:(n - 1)) for (b in (a + 1):n) {
      x <- imps[[a]]; y <- imps[[b]]
      if (all(is.na(x) == is.na(y))) {
        d <- which(x != y & !is.na(x))
        if (length(d) == 1) {
          z <- x; z[d] <- NA
          merged[[length(merged) + 1]] <- z; used[a] <- TRUE; used[b] <- TRUE
        }
      }
    }
    for (i in seq_len(n)) if (!used[i]) primes[[length(primes) + 1]] <- imps[[i]]
    if (length(merged) == 0) break
    key  <- vapply(merged, function(z) paste(ifelse(is.na(z), "-", z), collapse = ""), character(1))
    imps <- merged[!duplicated(key)]
  }
  pkey   <- vapply(primes, function(z) paste(ifelse(is.na(z), "-", z), collapse = ""), character(1))
  primes <- primes[!duplicated(pkey)]
  ## prime-implicant chart -> essential primes + greedy cover of the rest
  covers <- function(p, m) all(is.na(p) | p == m)
  chart  <- matrix(FALSE, length(primes), length(mt))
  for (i in seq_along(primes)) for (j in seq_along(mt)) chart[i, j] <- covers(primes[[i]], mt[[j]])
  chosen <- rep(FALSE, length(primes))
  for (j in seq_along(mt)) { cols <- which(chart[, j]); if (length(cols) == 1) chosen[cols] <- TRUE }
  covd <- if (any(chosen)) apply(chart[chosen, , drop = FALSE], 2, any) else rep(FALSE, length(mt))
  while (any(!covd)) {
    gains <- vapply(seq_along(primes),
                    function(i) if (chosen[i]) -1L else sum(chart[i, !covd]), integer(1))
    chosen[which.max(gains)] <- TRUE
    covd <- apply(chart[chosen, , drop = FALSE], 2, any)
  }
  ## rebuild a DNF expression from the chosen prime implicants
  clause <- function(p) {
    lits <- list()
    for (i in seq_len(k)) if (!is.na(p[i])) {
      v <- as.name(vars[i]); lits[[length(lits) + 1]] <- if (p[i] == 1) v else bquote(!.(v))
    }
    Reduce(function(a, b) bquote(.(a) & .(b)), lits)
  }
  cls <- lapply(primes[chosen], clause)
  Reduce(function(a, b) bquote(.(a) | .(b)), cls)
}

## ---- recursive soft De Morgan evaluation of a quoted Boolean formula -------
softValue <- function(expr, x, theta, n, fpos, fneg) {
  if (is.logical(expr)) return(as.numeric(expr))                       # TRUE / FALSE
  if (is.name(expr)) {                                                 # bare literal v
    gname <- as.character(expr)
    return(fpos(x[[gname]], theta[[gname]], n))
  }
  if (is.call(expr)) {
    op <- as.character(expr[[1]])
    if (op == "(") return(softValue(expr[[2]], x, theta, n, fpos, fneg))
    if (op == "!") {
      arg <- expr[[2]]
      if (is.name(arg)) {                                             # negated literal !v
        gname <- as.character(arg)
        return(fneg(x[[gname]], theta[[gname]], n))
      }
      return(1 - softValue(arg, x, theta, n, fpos, fneg))             # De Morgan negation
    }
    if (op == "&" || op == "&&")
      return(softValue(expr[[2]], x, theta, n, fpos, fneg) *
             softValue(expr[[3]], x, theta, n, fpos, fneg))
    if (op == "|" || op == "||") {
      a <- softValue(expr[[2]], x, theta, n, fpos, fneg)
      b <- softValue(expr[[3]], x, theta, n, fpos, fneg)
      return(1 - (1 - a) * (1 - b))
    }
  }
  stop("Unrecognised Boolean expression: ", deparse(expr))
}

## ---- optional minimal-DNF preprocessing of a rule set ----------------------
.prepRules <- function(rules, minimize) {
  if (!minimize) return(rules)
  lapply(rules, function(r) minimizeDNF(r, ruleVars(r)))
}

## ---- regulatory functions Phi_i(x) for the whole network (for inspection) --
booleanToLogisticPhi <- function(rules, theta, n,
                                 fpos = logisticp, fneg = logisticm, minimize = FALSE) {
  rules <- .prepRules(rules, minimize)
  genes <- names(rules)
  function(x) {
    xl <- as.list(x)
    setNames(vapply(genes,
                    function(gi) softValue(rules[[gi]], xl, theta, n, fpos, fneg),
                    numeric(1)), genes)
  }
}

## ---- deSolve vector field  x_i' = kappa_i Phi_i(x) - gamma_i x_i -----------
booleanToLogisticField <- function(rules, kappa, gamma, theta, n,
                                   fpos = logisticp, fneg = logisticm, minimize = FALSE) {
  rules <- .prepRules(rules, minimize)        # minimal-DNF preprocessing done ONCE here
  genes <- names(rules)
  function(t, s, p) {
    x  <- as.list(s)
    dx <- numeric(length(genes)); names(dx) <- genes
    for (gi in genes)
      dx[gi] <- kappa[[gi]] * softValue(rules[[gi]], x, theta, n, fpos, fneg) -
                gamma[[gi]] * s[[gi]]
    list(dx[genes])
  }
}
