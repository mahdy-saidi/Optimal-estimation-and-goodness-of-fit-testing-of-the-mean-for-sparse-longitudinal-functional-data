##########################################################################
#        APPLICATION TO REAL DATA -- DIFFUSION TENSOR IMAGING
#
#   NOTATION OF THE PAPER. The estimator is mu^sp of the spacing scheme,
#   expanded on the candidate bases B_{r,K}, r = 0, 1, 2 (half-cosine,
#   augmented by t, augmented by t and t^2). THROUGHOUT THIS FILE, AS IN
#   THE PAPER, K COUNTS THE COSINE TERMS c_1, ..., c_K, THE r MONOMIALS
#   OF THE AUGMENTATION EXCLUDED, so that B_{r,K} has dimension r + K.
#   basis_r(r, dim) is the only function that takes a DIMENSION, and
#   B_{r,K} = basis_r(r, r + K); the criterion searches K over 1..KMAX in
#   each system, that is dimensions r+1 .. r+KMAX, and K is what the
#   storage and the tables report.
#
#   WHAT IT COMPUTES, AND NOTHING ELSE:
#
#   (1) ESTIMATION, BY NESTED CROSS-VALIDATION. The patients are split
#   into V folds; on the V-1 training folds the criterion (eq:cv) selects
#   the pair (r, K) of the spacing estimator, the bandwidth h_t of the
#   local-linear smoother and the penalty lambda of the P-spline over V
#   INNER folds; each estimator is refitted on the training patients at
#   its own tuning and scored on the held-out ones by
#       E_e = (1/(N_u M_test)) sum ( Y_ij(u_l) - mu-hat_e(u_l, T_ij) )^2 ,
#       Delta_e = log( E_e / E_sp ),   e in {LL, spl}.
#   The outer split is repeated R_OUT times, and the section reports on
#   how many of those blocks the spacing estimator has the smaller
#   held-out error. THE LOCAL-LINEAR SMOOTHER IS FITTED IN THE VISIT-TIME
#   DIRECTION ONLY, as in the simulation: no smoothing in u is needed,
#   since each profile is observed over the whole of U, so ll_crit() has
#   ONE grid, LL_HGRID, and h_u has disappeared from the selection, from
#   the storage and from the tables.
#
#   EVERY BLOCK IS REPRODUCIBLE ON ITS OWN. Repetition b draws its outer
#   partition from the seed SEEDOUT(b) = BASESEED + 100 b, and block
#   (b, v) draws its inner partitions from SEEDBLK(b, v) = SEEDOUT(b) + v,
#   so that one_block(b, v) depends on nothing but the data and those two
#   integers: it can be replayed alone, in any order, without running the
#   other blocks. replay_block(b, v) does exactly that, and check_blocks()
#   replays blocks of a stored run and compares them to what was stored.
#
#   (2) TESTING H_0 : mu(u,t) = m(u), TWICE. The statistic is Q^sp_n of
#   (eq:def_Hn_sp) and its multiplier bootstrap Q^{sp,*}_n of
#   (eq:def_Hn_sp_boot): pairs of visits of DIFFERENT patients, window
#   integrals Phi-hat_kl of the spacing scheme, one two-point multiplier
#   zeta_i per patient, and the scale factor (M+1)^2/{n(n-1)}. Since
#   phi_1 is constant, the null says that the coefficients of the
#   non-constant modes vanish, and the section reads it in two ways, both
#   on the half-cosine basis B_0 (the augmentation plays no role here):
#       centred : profiles centred by the pooled mean of the sample,
#                 m-hat(u) = M^{-1} sum_ij Y_ij(u), modes 1..KBAR ;
#       raw     : profiles left as they are, modes 2..KBAR+1, which are
#                 orthogonal to the constants, so no centring is needed.
#   Both give (Q^sp_n, q^{sp,*}_{1-a}, p-value). What is STORED is the
#   observed statistic and the N_b bootstrap draws themselves, so the
#   level ALPHA -- hence the reported quantile -- is a parameter of the
#   TABLE layer, not of the run.
#
#       run_dti()  ->  res_dti.rds   then   make_dti()  ->  .tex / .csv
#
#   (3) FIGURES: OBSERVED PROFILES AND FITTED MEAN. Three panels, one per
#   estimator, each fitted ONCE on the whole sample at the tuning
#   cv_select() selected there (S$sel, the very tuning of tab_dti): every
#   observed profile Y_ij(.) in light grey against u, and the fitted
#   curve mu-hat_e(u, t_j) overlaid in black at three visit times t_j
#   (10th, 50th, 90th percentile of the pooled visit times), dotted for
#   the earliest, dashed for the middle, solid for the latest.
#
#       run_dti()  ->  res_dti.rds   then   make_dti()  ->  .tex / .csv
#                                    and     make_dti_profiles() -> panels
#   CHECKPOINTS AND RESUMPTION, as in simul_estimation.R. res_dti.rds is
#   rewritten every CKPT_EVERY blocks, atomically (<name>.rds.tmp, then
#   <name>.rds.bak, then the rename), each block carrying its own seeds,
#   so a machine that stops loses at most CKPT_EVERY blocks and a restart
#   with RESUME = TRUE produces exactly what an uninterrupted run would
#   have produced. R_OUT may be raised between two runs: the blocks
#   already stored are kept and the new ones appended.
#
#   OUTPUTS. tab_dti (the selected tunings and the held-out error),
#   dti_numbers.tex (the macros used by the text, so that no number is
#   typed by hand), dti_blocks.csv, which reports the seed, the tuning
#   and the two log-ratios of EVERY block, and fig_dti_sp/fig_dti_ll/
#   fig_dti_spl (.tex and .pdf), the three panels of (3). The seeds of
#   the blocks on which the spacing estimator wins are printed to the
#   console and nowhere else.
##########################################################################

suppressPackageStartupMessages(library(splines))   # spline benchmark

# =========================================================================
#  1. USER PARAMETERS  --  everything meant to be edited sits here
# =========================================================================

## ---- 1.1 the data --------------------------------------------------------
ARM       <- 1L        # 1 = multiple sclerosis arm of DTI, 0 = controls
TRACT     <- "cca"     # "cca" (corpus callosum) or "rcst"
DROP_BASE <- TRUE      # drop the reference scans, whose visit time is 0

## ---- 1.2 estimation ------------------------------------------------------
## KMAX bounds the K of the paper, the r monomials excluded: the criterion
## searches the dimensions r+1 .. r+KMAX in each system. LL_HGRID must
## reach far enough to the right: on these data the criterion asks for the
## strongest smoothing available in the visit-time direction, and a
## selection at the last value is reported by make_dti() as a grid to
## widen.
KMAX      <- 12L                   # largest K considered, monomials excluded
V_FOLDS   <- 10L                   # folds, of the outer split and the inner one
R_OUT     <- 5L                    # repetitions of the outer split
R_SEL     <- 5L                    # partitions of the whole-sample selection
LL_HGRID  <- exp(seq(log(0.05), log(0.60), length.out = 20))   # h_t
SPL_NIK   <- 20L                                               # spline knots
SPL_LGRID <- exp(seq(log(1e-5), log(1e1), length.out = 20))    # penalties

## ---- 1.3 the test --------------------------------------------------------
## KBAR is the number of modes of the statistic: the centred reading uses
## K_M = {1, ..., KBAR}, the raw one K_M = {2, ..., KBAR+1}. ALPHA is read
## by make_dti() alone, on the stored bootstrap draws.
KBAR   <- 25L        # number of modes of the statistic
B_TEST <- 1999L      # bootstrap draws
LAW    <- "mammen"   # multipliers: "mammen" (two-point) or "normal"
ALPHA  <- 0.05       # nominal level; only make_dti() reads it

## ---- 1.4 checkpoints -----------------------------------------------------
CKPT_EVERY <- 5L     # blocks between two writes
RESUME     <- TRUE   # FALSE ignores an existing file and overwrites it

## ---- 1.5 output and seeds ------------------------------------------------
OUTDIR   <- Sys.getenv("SIM_OUTDIR", unset = "./output")
LATEX    <- Sys.getenv("SIM_LATEX",  unset = "pdflatex")   # compiles the panels
BASESEED <- 2026L
RDS_DTI  <- "res_dti"
CEXF     <- 1.2      # typography multiplier of the three profile panels

## One seed per outer repetition, one per block. The offsets keep them
## apart from the two seeds of the whole-sample stage, and from each other
## as long as there are fewer than 100 folds.
stopifnot(V_FOLDS < 100L)
SEEDOUT  <- function(b)    BASESEED + 100L * as.integer(b)
SEEDBLK  <- function(b, v) SEEDOUT(b) + as.integer(v)
SEEDSEL  <- BASESEED + 1L    # whole-sample selection
SEEDBOOT <- BASESEED + 2L    # multiplier bootstrap

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# =========================================================================
#  2. THE DATA
# =========================================================================

## dti_prepare(): the retained visits, as the model of the paper reads
## them. Every patient with at least one follow-up visit is kept.
dti_prepare <- function(arm = ARM) {
  if (!requireNamespace("refund", quietly = TRUE))
    stop("the DTI study ships with the package refund, which is not installed")
  e <- new.env(parent = emptyenv())
  utils::data("DTI", package = "refund", envir = e)
  D <- e$DTI
  D <- D[!is.na(D$case) & D$case == arm, , drop = FALSE]
  if (!nrow(D)) stop(sprintf("arm %d is not in the DTI data", arm))
  D <- D[stats::complete.cases(D[[TRACT]]) & !is.na(D$visit.time), , drop = FALSE]
  if (DROP_BASE) D <- D[D$visit.time > 0, , drop = FALSE]
  if (!nrow(D))
    stop(sprintf(paste("arm %d has no visit left once the reference scans are",
                       "dropped: in this study the controls are scanned once,",
                       "at visit time 0"), arm))
  Y    <- D[[TRACT]]
  days <- as.numeric(D$visit.time)
  tt   <- (days - min(days)) / (max(days) - min(days))
  subj <- as.integer(factor(D$ID))
  mi   <- as.integer(table(subj))
  list(T = tt, Y = t(Y), subj = subj, n = max(subj), M = length(tt),
       ug = seq(0, 1, length.out = ncol(Y)), Nu = ncol(Y),
       days = days, mi = mi, n_one = sum(mi == 1L))
}


# =========================================================================
#  3. CANDIDATE BASES  --  B_{r,K}, r = 0 (half-cosine), 1 (+ t), 2 (+ t, t^2)
#
#  basis_r(r, dim) orthonormalises, by closed-form Gram-Schmidt, the first
#  'dim' elements of the ordered dictionary
#      r = 0 : (c_1, c_2, c_3, ...)
#      r = 1 : (c_1, t, c_2, c_3, ...)
#      r = 2 : (c_1, t, t^2, c_2, c_3, ...)
#  where c_1 = 1 and c_k(t) = sqrt2 cos((k-1) pi t). Its argument is a
#  DIMENSION, and the paper's B_{r,K} is basis_r(r, r + K). Each phi_k is
#  a polynomial of degree at most 2 plus finitely many cosines, so that
#  both its values and its window integrals are closed-form.
# =========================================================================

basis_r <- function(r, dim) {
  K      <- dim                    # number of dictionary elements retained
  poly   <- matrix(0, K, 3)        # coefficients of phi_k over 1, t, t^2
  freqs  <- vector("list", K)      # cosine frequencies of phi_k
  thetas <- vector("list", K)      # matching cosine amplitudes
  poly[1, ] <- c(1, 0, 0)          # phi_1 = 1 in every candidate
  
  P2 <- c(-sqrt(3),  2 * sqrt(3), 0)              # phi_2 = sqrt3 (2t - 1)
  P3 <- c( sqrt(5), -6 * sqrt(5), 6 * sqrt(5))    # phi_3 = sqrt5 (6t^2-6t+1)
  
  if (r == 0L) {                                  # plain half-cosine system
    if (K >= 2) for (k in 2:K) { freqs[[k]] <- (k - 1) * pi; thetas[[k]] <- 1 }
    return(list(poly = poly, freqs = freqs, thetas = thetas, K = K, r = 0L))
  }
  
  ## projections of the cosines on phi_2 (odd block) and phi_3 (even block),
  ## and the partial sums driving the rank-one Gram-Schmidt recursion
  j <- seq_len(max(K, 4L))
  al <- -4 * sqrt(6)  / ((2 * j - 1)^2 * pi^2)
  be <-  3 * sqrt(10) / (j^2 * pi^2)
  A  <- c(1, 1 - cumsum(al^2))
  B  <- c(1, 1 - cumsum(be^2))
  
  if (r == 1L) {                                  # dictionary (c_1, t, c_2, ...)
    if (K >= 2) poly[2, ] <- P2
    if (K >= 3) for (m in seq_len(K - 2)) {
      k <- 2 + m
      if (m %% 2 == 0) {                          # symmetric: already orthogonal
        freqs[[k]] <- m * pi; thetas[[k]] <- 1
      } else {                                    # antisymmetric: remove phi_2
        nu <- (m + 1) %/% 2; th <- numeric(nu)
        th[nu] <- sqrt(A[nu] / A[nu + 1])
        if (nu >= 2) th[seq_len(nu - 1)] <-
          al[nu] * al[seq_len(nu - 1)] / sqrt(A[nu] * A[nu + 1])
        poly[k, ]   <- (-al[nu] / sqrt(A[nu] * A[nu + 1])) * P2
        freqs[[k]]  <- (2 * seq_len(nu) - 1) * pi
        thetas[[k]] <- th
      }
    }
    return(list(poly = poly, freqs = freqs, thetas = thetas, K = K, r = 1L))
  }
  
  if (K >= 2) poly[2, ] <- P2                     # dictionary (c_1, t, t^2, ...)
  if (K >= 3) poly[3, ] <- P3
  if (K >= 4) for (m in seq_len(K - 3)) {
    k <- 3 + m
    if (m %% 2 == 1) {                            # odd m  -> phi_2, odd freqs
      nu <- (m + 1) %/% 2; th <- numeric(nu)
      th[nu] <- sqrt(A[nu] / A[nu + 1])
      if (nu >= 2) th[seq_len(nu - 1)] <-
        al[nu] * al[seq_len(nu - 1)] / sqrt(A[nu] * A[nu + 1])
      poly[k, ]  <- (-al[nu] / sqrt(A[nu] * A[nu + 1])) * P2
      freqs[[k]] <- (2 * seq_len(nu) - 1) * pi
    } else {                                      # even m -> phi_3, even freqs
      nu <- m %/% 2; th <- numeric(nu)
      th[nu] <- sqrt(B[nu] / B[nu + 1])
      if (nu >= 2) th[seq_len(nu - 1)] <-
        be[nu] * be[seq_len(nu - 1)] / sqrt(B[nu] * B[nu + 1])
      poly[k, ]  <- (-be[nu] / sqrt(B[nu] * B[nu + 1])) * P3
      freqs[[k]] <- (2 * seq_len(nu)) * pi
    }
    thetas[[k]] <- th
  }
  list(poly = poly, freqs = freqs, thetas = thetas, K = K, r = 2L)
}

## basis_phi(): phi_1..phi_dim evaluated at t. Returns length(t) x dim.
basis_phi <- function(co, t) {
  Phi <- matrix(0, length(t), co$K)
  for (k in seq_len(co$K)) {
    v <- co$poly[k, 1] + co$poly[k, 2] * t + co$poly[k, 3] * t^2
    if (!is.null(co$freqs[[k]]))
      v <- v + sqrt(2) * as.vector(cos(outer(t, co$freqs[[k]])) %*% co$thetas[[k]])
    Phi[, k] <- v
  }
  Phi
}

## basis_int(): the window integrals int_lo^hi phi_k, in closed form.
basis_int <- function(co, lo, hi) {
  D  <- matrix(0, length(lo), co$K)
  d1 <- hi - lo; d2 <- hi^2 - lo^2; d3 <- hi^3 - lo^3
  for (k in seq_len(co$K)) {
    v <- co$poly[k, 1] * d1 + co$poly[k, 2] * d2 / 2 + co$poly[k, 3] * d3 / 3
    if (!is.null(co$freqs[[k]])) {
      om <- co$freqs[[k]]
      v <- v + sqrt(2) * as.vector(
        (sin(outer(hi, om)) - sin(outer(lo, om))) %*% (co$thetas[[k]] / om))
    }
    D[, k] <- v
  }
  D
}

## check_bases(): orthonormality and closed-form window integrals, on a
## Simpson grid. Cheap, and run once by the driver.
check_bases <- function(dim = 16L, ng = 20001L) {
  tq <- seq(0, 1, length.out = ng)
  wq <- c(1, rep(c(4, 2), length.out = ng - 2L), 1) / (3 * (ng - 1))
  for (r in 0:2) {
    Phi <- basis_phi(basis_r(r, dim), tq)
    e <- max(abs(crossprod(Phi * wq, Phi) - diag(dim)))
    cat(sprintf("  [B_%d] max|Gram - I| = %.2e   sup_k ||phi_k||_inf = %.4f\n",
                r, e, max(abs(Phi))))
    stopifnot(e < 1e-6)
  }
  co <- basis_r(2L, dim); lo <- c(0.03, 0.44); hi <- lo + c(0.11, 0.28)
  Dq <- t(sapply(seq_along(lo), function(i) {
    tt <- seq(lo[i], hi[i], length.out = 4001)
    ww <- c(1, rep(c(4, 2), length.out = 3999), 1) * (hi[i] - lo[i]) / (3 * 4000)
    as.vector(ww %*% basis_phi(co, tt))
  }))
  e <- max(abs(basis_int(co, lo, hi) - Dq))
  cat(sprintf("  [B_2] max|window integral - quadrature| = %.2e\n", e))
  stopifnot(e < 1e-9)
  invisible(TRUE)
}


# =========================================================================
#  4. THE SPACING ESTIMATOR
#
#  beta^sp_k(u) = sum_l Y_(l)(u) Phi-hat_kl, with
#  Phi-hat_kl = (1/2h) int_{T_(l-h)}^{T_(l+h)} phi_k and the window
#  endpoints obtained by symmetric reflection at the two ends of the
#  pooled design. It uses NO design density: the geometry of the order
#  statistics reproduces the correction 1/g from the data alone. The same
#  window integrals carry the estimator and the test statistic, hence
#  w_spacing(), which returns them in the order the observations are
#  given.
# =========================================================================

## h_of_M(): the window half-width, h = ceil((1 + ln ln(M+20))/2); it grows
## so slowly that h = 2 for every realistic sample size.
h_of_M <- function(M) ceiling(0.5 + 0.5 * log(log(M + 20)))

w_spacing <- function(Tv, co, h = h_of_M(length(Tv))) {
  M <- length(Tv); ord <- order(Tv); Ts <- Tv[ord]
  Text <- function(rk) {                         # reflected order statistics
    out <- numeric(length(rk)); lo <- rk <= 0; hi <- rk > M; in_ <- !lo & !hi
    out[in_] <- Ts[rk[in_]]
    out[lo]  <- 2 * Ts[1] - Ts[2 - rk[lo]]
    out[hi]  <- 2 * Ts[M] - Ts[2 * M - rk[hi]]
    out
  }
  l <- seq_len(M)
  W <- matrix(0, M, co$K)
  W[ord, ] <- basis_int(co, Text(l - h), Text(l + h)) / (2 * h)
  W
}

beta_spacing <- function(Y, Tv, co, h = h_of_M(length(Tv))) Y %*% w_spacing(Tv, co, h)

## sse_curve(): held-out sum of squares of the cross-validation criterion,
## as a function of the DIMENSION of the fit, entry r + K being the
## candidate B_{r,K}. The factor 1/(N_u M) is common to all candidates and
## is applied by the caller.
sse_curve <- function(beta, Phi_te, Y_te, Kup) {
  out <- numeric(Kup); Rk <- Y_te
  for (K in seq_len(Kup)) {
    Rk <- Rk - outer(beta[, K], Phi_te[, K])
    out[K] <- sum(Rk * Rk)
  }
  out
}


# =========================================================================
#  5. THE TWO BENCHMARKS
#
#  Each is evaluated by ONE function, used twice with the same code path:
#  with stat = "sum" and the held-out times and profiles it returns the
#  held-out sum of squares of every candidate tuning, which is the fold's
#  contribution to the criterion; with a single tuning it returns the sum
#  of squares of the refitted estimator on the held-out patients. The
#  criterion is therefore literally the same object as the one used for
#  the spacing estimator, up to the constant 1/(N_u M).
# =========================================================================

## ll_crit(): the pooled local-linear smoother of CM2012 in the visit-time
## direction, fitted on (Tp, Y) and evaluated at every (u_l, t_eval). NO
## SMOOTHING IN u: each profile is observed over the whole of U, so the
## fit is a section of the pooled scatterplot at each u, and the normal
## equations are those of one univariate local-linear fit shared by all
## the slices.
ll_crit <- function(Y, Tp, ug, t_eval, target, hgrid = LL_HGRID,
                    stat = c("mean", "sum")) {
  stat <- match.arg(stat)
  Nu <- length(ug)
  Dt <- outer(Tp, t_eval, "-")
  fb <- rowMeans(Y)                       # fallback: pooled profile mean
  vapply(hgrid, function(h) {
    Kt <- pmax(1 - (Dt / h)^2, 0)                      # Epanechnikov in t
    T0 <- colSums(Kt); T1 <- colSums(Kt * Dt); T2 <- colSums(Kt * Dt * Dt)
    A0 <- Y %*% Kt; A1 <- Y %*% (Kt * Dt)
    den <- T0 * T2 - T1 * T1
    Fh <- (A0 * rep(T2, each = Nu) - A1 * rep(T1, each = Nu)) /
      rep(pmax(den, 1e-300), each = Nu)
    bad <- which(den <= 1e-10)                         # design too degenerate
    if (length(bad)) Fh[, bad] <- (A0 / rep(pmax(T0, 1e-300), each = Nu))[, bad]
    empty <- which(T0 <= 1e-12)                        # no design point in view
    if (length(empty)) Fh[, empty] <- fb
    if (stat == "mean") mean((Fh - target)^2) else sum((Fh - target)^2)
  }, numeric(1))
}

## spl_crit(): penalised cubic B-spline smoother in t, fitted slice-wise in
## u on (Tp, Y) and evaluated at t_eval. The normal system is factorised
## once per lambda and applied to all slices at once.
spl_crit <- function(Y, Tp, t_eval, target, lgrid = SPL_LGRID, nik = SPL_NIK,
                     stat = c("mean", "sum")) {
  stat <- match.arg(stat)
  kn <- seq(0, 1, length.out = nik + 2L)[-c(1L, nik + 2L)]
  Bm <- splines::bs(Tp, knots = kn, degree = 3L, intercept = TRUE,
                    Boundary.knots = c(0, 1))
  Bt <- predict(Bm, t_eval); q <- ncol(Bm); M <- length(Tp)
  Pen <- crossprod(diff(diag(q), differences = 2L))   # 2nd-order difference
  BtB <- crossprod(Bm); YB <- Y %*% Bm
  vapply(lgrid, function(lam) {
    A <- tryCatch(chol2inv(chol(BtB + lam * M * Pen + 1e-10 * diag(q))),
                  error = function(e) NULL)
    if (is.null(A)) return(NA_real_)
    Rm <- YB %*% A %*% t(Bt) - target
    if (stat == "mean") mean(Rm^2) else sum(Rm^2)
  }, numeric(1))
}


# =========================================================================
#  6. THE STATISTIC AND ITS MULTIPLIER BOOTSTRAP
#
#  Q^sp_n of (eq:def_Hn_sp) is, up to the factor (M+1)^2/{n(n-1)}, the
#  sum over the pairs of visits of DIFFERENT patients of
#      sum_{k in K_M} Phi-hat_kl Phi-hat_kl' int_U Y_l Y_l' dnu ,
#  that is sum(G) with G_ii' = sum_k <v_ik, v_i'k>_{L^2(U)} and the
#  diagonal set to zero, where v_ik(u) = sum_j Phi-hat_{k,l(i,j)} Y_ij(u)
#  is the contribution of patient i to the mode k. The bootstrap version
#  Q^{sp,*}_n gives every visit of patient i the same multiplier zeta_i,
#  so it is the quadratic form zeta' G zeta.
# =========================================================================

## trapz_w(): quadrature weights on the u-grid.
trapz_w <- function(u) {
  n <- length(u); w <- numeric(n)
  w[1] <- (u[2] - u[1]) / 2; w[n] <- (u[n] - u[n - 1]) / 2
  if (n > 2) w[2:(n - 1)] <- (u[3:n] - u[1:(n - 2)]) / 2
  w
}

## level_hat(): the pooled mean profile m-hat(u) = M^{-1} sum_ij Y_ij(u),
## which the centred reading of the test removes once, on the whole sample.
level_hat <- function(Y) rowMeans(Y)

## subject_blocks(): the (dim * nu) x n matrix whose column i stacks, mode
## by mode, the contribution v_ik of patient i on the u-grid, already
## multiplied by sqrt(w_u) so that the Euclidean inner product of two
## columns is the L^2(U) inner product summed over the modes.
subject_blocks <- function(Y, W, subj, n, wu) {
  nu <- length(wu); K <- ncol(W)
  Ys <- sqrt(wu) * Y
  Vs <- matrix(0, K * nu, n)
  for (k in seq_len(K))
    Vs[((k - 1) * nu + 1):(k * nu), ] <- t(rowsum(t(Ys * rep(W[, k], each = nu)), subj))
  Vs
}

## gram_modes(): the n x n matrix of inner products restricted to a SET of
## modes, with the diagonal set to zero, so that sum(G) is the
## between-patient statistic and t(zeta) G zeta its bootstrap version.
gram_modes <- function(V, modes, nu) {
  rows <- as.vector(outer(seq_len(nu), (modes - 1L) * nu, "+"))
  G <- crossprod(V[rows, , drop = FALSE])
  diag(G) <- 0
  G
}

MAMMEN <- list(lo = -(sqrt(5) - 1) / 2, hi = (sqrt(5) + 1) / 2,
               p  = (sqrt(5) + 1) / (2 * sqrt(5)))
draw_xi <- function(law, n, B) switch(law,
                                      normal = matrix(rnorm(n * B), n, B),
                                      mammen = matrix(ifelse(runif(n * B) < MAMMEN$p, MAMMEN$lo, MAMMEN$hi), n, B),
                                      stop("unknown multiplier law"))

boot_quad <- function(G, Xi) colSums(Xi * (G %*% Xi))
pval      <- function(Qobs, Qstar) (1 + sum(Qstar >= Qobs)) / (length(Qstar) + 1)

## dti_test(): one reading of H_0 : mu(u,t) = m(u). 'modes' is the set
## K_M of the paper, 'centre' says whether the pooled mean profile is
## removed first, and the multipliers are passed in so that the two
## readings share them and are paired. Returns the observed statistic and
## the N_b bootstrap draws, both on the scale of (eq:def_Hn_sp).
##
## 'lead' is what the level m alone puts into the statistic: writing
## s_ik = sum_j Phi-hat_{k,l(i,j)} and Y_ij = m + (centred part), the
## profiles contribute ||m||^2 sum_{k in K_M} { (sum_i s_ik)^2 -
## sum_i s_ik^2 } to Q^sp_n. It vanishes with the centring and is only
## as small as the window quadrature is accurate without it, so it is
## printed next to the raw statistic and read together with it.
dti_test <- function(Y, Tv, subj, n, ug, co, modes, centre, Xi) {
  M <- length(Tv); nu <- length(ug); wu <- trapz_w(ug)
  W  <- w_spacing(Tv, co)
  mh <- level_hat(Y)
  Yc <- if (centre) Y - mh else Y
  G  <- gram_modes(subject_blocks(Yc, W, subj, n, wu), modes, nu)
  cst <- (M + 1)^2 / (n * (n - 1))
  s  <- t(rowsum(W, subj))[modes, , drop = FALSE]
  lead <- if (centre) 0 else
    cst * sum(wu * mh^2) * sum(rowSums(s)^2 - rowSums(s^2))
  list(Q = cst * sum(G), Qstar = cst * boot_quad(G, Xi), lead = lead,
       modes = as.integer(range(modes)), centre = centre, G = G)
}

## check_test(): the bootstrap law is centred, with variance
## 2 sum_{i != i'} G_ii'^2, whatever the law of the multipliers.
check_test <- function(G, B = 20000L) {
  TS <- boot_quad(G, draw_xi(LAW, ncol(G), B))
  cat(sprintf("  E*[Q*] = %+.3e (0)   sd*/target = %.4f\n",
              mean(TS), sd(TS) / sqrt(2 * sum(G^2))))
  stopifnot(abs(mean(TS)) < 6 * sd(TS) / sqrt(B),
            abs(sd(TS) / sqrt(2 * sum(G^2)) - 1) < 0.05)
  invisible(TRUE)
}


# =========================================================================
#  7. STORAGE  --  one .rds, written ATOMICALLY, every table being built
#                  afterwards from that file alone
# =========================================================================

rds_file <- function(name) file.path(OUTDIR, paste0(name, ".rds"))

## save_rds(): a write that a power cut cannot corrupt. The object goes to
## <name>.rds.tmp; only once that file is complete is the current
## <name>.rds renamed <name>.rds.bak and the temporary file moved in its
## place. At every instant, at least one of the two files is a complete
## .rds, and load_rds() falls back on the backup if the main file is
## unreadable.
save_rds <- function(name, obj, quiet = FALSE) {
  f <- rds_file(name); tmp <- paste0(f, ".tmp"); bak <- paste0(f, ".bak")
  saveRDS(obj, tmp)
  if (file.exists(f)) {
    unlink(bak)
    if (!file.rename(f, bak)) { file.copy(f, bak, overwrite = TRUE); unlink(f) }
  }
  if (!file.rename(tmp, f)) { file.copy(tmp, f, overwrite = TRUE); unlink(tmp) }
  if (!quiet)
    cat(sprintf("  wrote %s.rds (%.0f kB)\n", name, file.size(f) / 1024))
  invisible(obj)
}

load_rds <- function(name) {
  f <- rds_file(name); bak <- paste0(f, ".bak")
  obj <- if (file.exists(f)) tryCatch(readRDS(f), error = function(e) NULL) else NULL
  if (is.null(obj) && file.exists(bak)) {
    message(sprintf("  [!] %s.rds unreadable: falling back on %s.rds.bak", name, name))
    obj <- tryCatch(readRDS(bak), error = function(e) NULL)
  }
  if (is.null(obj))
    stop(sprintf("%s not found or unreadable: run the experiment first.", f))
  obj
}

has_rds <- function(name)
  file.exists(rds_file(name)) || file.exists(paste0(rds_file(name), ".bak"))

## clean_ckpt(): remove the two auxiliary files of the atomic write once
## they are of no further use, the backup only after <name>.rds has been
## READ BACK in full.
clean_ckpt <- function(name, quiet = FALSE) {
  f <- rds_file(name); tmp <- paste0(f, ".tmp"); bak <- paste0(f, ".bak")
  if (file.exists(tmp)) unlink(tmp)
  ok <- file.exists(f) && !is.null(tryCatch(readRDS(f), error = function(e) NULL))
  if (!ok) {
    message(sprintf("  [!] %s.rds is missing or unreadable: %s.rds.bak is kept",
                    name, name))
    return(invisible(FALSE))
  }
  if (file.exists(bak)) unlink(bak)
  if (!quiet) cat(sprintf("  %s.rds verified; .tmp and .bak removed\n", name))
  invisible(TRUE)
}


# =========================================================================
#  8. OUTPUT HELPERS  --  booktabs tables and LaTeX numbers
# =========================================================================

write_table <- function(name, align, header, body, caption, label) {
  writeLines(c("\\begin{table}[H]", "\\centering",
               sprintf("\\begin{tabular}{%s}", align), "\\toprule",
               header, "\\midrule", body, "\\bottomrule", "\\end{tabular}",
               sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
               "\\end{table}"),
             file.path(OUTDIR, paste0(name, ".tex")))
  cat(sprintf("  wrote %s.tex\n", name))
}

## num_tex(): a number in LaTeX, in scientific form when it is small.
num_tex <- function(x, dig = 3) {
  if (!is.finite(x)) return("---")
  if (x == 0) return("$0$")
  if (abs(x) < 1e-2 || abs(x) >= 1e3) {
    e <- floor(log10(abs(x)))
    sprintf("$%.1f\\times 10^{%d}$", x / 10^e, e)
  } else sprintf("$%.*f$", dig, x)
}
num_math <- function(x, dig = 3) gsub("\\$", "", num_tex(x, dig))

## seed_list(): a set of seeds, as the console prints them.
seed_list <- function(s) if (!length(s)) "none" else paste(s, collapse = ", ")


# =========================================================================
#  9. ONE SELECTION BY THE CRITERION OF THE PAPER
#
#  The pair (r, K) of the spacing estimator over the 3 x KMAX candidates,
#  the bandwidth h_t of the local-linear smoother over LL_HGRID and the
#  penalty lambda of the P-spline over SPL_LGRID, all minimising the same
#  criterion (eq:cv) on the same partition of the patients, repeated over
#  'reps' partitions. K is the K of the paper: the search over the
#  dimensions r+1 .. r+KMAX is read back on the slice r + 1..KMAX.
# =========================================================================

cv_select <- function(Y, Tv, subj, ug, co, V = V_FOLDS, reps = 1L) {
  ids <- sort(unique(subj)); n <- length(ids); lab <- match(subj, ids)
  nb <- length(co); nht <- length(LL_HGRID); nlam <- length(SPL_LGRID)
  A_sp <- matrix(0, nb, KMAX); A_ll <- numeric(nht); A_spl <- numeric(nlam)
  for (b in seq_len(reps)) {
    fold <- sample(rep_len(seq_len(V), n))[lab]
    for (v in seq_len(V)) {
      tr <- fold != v
      if (!any(tr) || all(tr)) next
      Ttr <- Tv[tr];  Ytr <- Y[, tr,  drop = FALSE]
      Tte <- Tv[!tr]; Yte <- Y[, !tr, drop = FALSE]
      htr <- h_of_M(length(Ttr))
      for (bq in seq_len(nb)) {
        r0 <- bq - 1L                      # B_{r0,K} has dimension r0 + K
        A_sp[bq, ] <- A_sp[bq, ] +
          sse_curve(beta_spacing(Ytr, Ttr, co[[bq]], htr),
                    basis_phi(co[[bq]], Tte), Yte, r0 + KMAX)[r0 + seq_len(KMAX)]
      }
      A_ll  <- A_ll  + ll_crit(Ytr, Ttr, ug, Tte, Yte, stat = "sum")
      A_spl <- A_spl + spl_crit(Ytr, Ttr, Tte, Yte, stat = "sum")
    }
  }
  ij <- which(A_sp == min(A_sp, na.rm = TRUE), arr.ind = TRUE)[1, ]
  il <- which.min(A_ll); is <- which.min(A_spl)
  list(basis = as.integer(ij[1]) - 1L, K = as.integer(ij[2]),
       h_t = LL_HGRID[il], lambda = SPL_LGRID[is])
}


# =========================================================================
#  10. ONE HELD-OUT BLOCK, REPRODUCIBLE ON ITS OWN
#
#  one_block(b, v) is a pure function of the data and of the two integers
#  (b, v): SEEDOUT(b) fixes the outer partition of repetition b and
#  SEEDBLK(b, v) the inner partition of that block. Nothing is carried
#  over from the blocks computed before it, so the blocks can be replayed
#  one by one, in any order.
# =========================================================================

one_block <- function(dat, co, b, v, V = V_FOLDS) {
  Tv <- dat$T; Y <- dat$Y; subj <- dat$subj; ug <- dat$ug
  n <- dat$n; Nu <- dat$Nu
  
  set.seed(SEEDOUT(b))
  fold <- sample(rep_len(seq_len(V), n))[subj]
  tr   <- fold != v
  Ttr <- Tv[tr];  Ytr <- Y[, tr,  drop = FALSE]; str <- subj[tr]
  Tte <- Tv[!tr]; Yte <- Y[, !tr, drop = FALSE]
  
  set.seed(SEEDBLK(b, v))
  s <- cv_select(Ytr, Ttr, str, ug, co, V = V, reps = 1L)
  
  den <- Nu * length(Tte)
  cs  <- co[[s$basis + 1L]]
  ds  <- s$basis + s$K                    # the selected system is B_{r,K}
  E <- c(sp  = sse_curve(beta_spacing(Ytr, Ttr, cs, h_of_M(length(Ttr))),
                         basis_phi(cs, Tte), Yte, ds)[ds] / den,
         ll  = unname(ll_crit(Ytr, Ttr, ug, Tte, Yte,
                              hgrid = s$h_t, stat = "sum")) / den,
         spl = unname(spl_crit(Ytr, Ttr, Tte, Yte,
                               lgrid = s$lambda, stat = "sum")) / den)
  list(E = E, sel = unlist(s), seed = SEEDBLK(b, v), n_test = length(Tte))
}

## replay_block(): one block from scratch, given nothing but (b, v).
replay_block <- function(b, v, dat = NULL, V = V_FOLDS) {
  if (is.null(dat)) dat <- dti_prepare()
  one_block(dat, lapply(0:2, basis_r, dim = 2L + KMAX), b, v, V)
}

## check_blocks(): replay stored blocks and compare them, number by number.
check_blocks <- function(S = load_rds(RDS_DTI), which = NULL, dat = NULL) {
  if (is.null(dat)) dat <- dti_prepare()
  co  <- lapply(0:2, basis_r, dim = 2L + KMAX)
  idx <- if (is.null(which)) which(!is.na(S$E[, "sp"])) else which
  d <- 0
  for (q in idx) {
    r <- one_block(dat, co, S$BLK[q, "b"], S$BLK[q, "v"], S$V)
    stopifnot(r$seed == S$BLK[q, "seed"])
    d <- max(d, max(abs(r$E - S$E[q, ])), max(abs(r$sel - S$SEL[q, ])))
  }
  cat(sprintf("  replayed %d block(s): max|replay - stored| = %.2e\n",
              length(idx), d))
  stopifnot(d == 0)
  invisible(TRUE)
}


# =========================================================================
#  11. THE EXPERIMENT
# =========================================================================

## dti_sig(): the parameters two runs must share for their blocks to be
## comparable, hence for a stored file to be continued. R_OUT is NOT in
## it: outer repetitions can be added. ALPHA is not either: it is read by
## the table layer alone.
dti_sig <- function(V) list(
  arm = ARM, tract = TRACT, drop_base = DROP_BASE,
  Kmax = KMAX, V = V,
  k_convention = "K excludes the r monomials; dimension is r + K",
  h_grid_t = LL_HGRID, spl_nik = SPL_NIK, lam_grid = SPL_LGRID,
  kbar = KBAR, B = B_TEST, law = LAW, baseseed = BASESEED)

## dti_empty(): the storage of a run, all blocks pending.
dti_empty <- function(sig, V, r_out) {
  nblk <- r_out * V
  list(sig = sig, V = V, r_out = r_out, r_sel = R_SEL, kbar = KBAR, B = B_TEST,
       seed = BASESEED, seed_sel = SEEDSEL, seed_boot = SEEDBOOT, secs = 0,
       E   = matrix(NA_real_, nblk, 3L,
                    dimnames = list(NULL, c("sp", "ll", "spl"))),
       SEL = matrix(NA_real_, nblk, 4L,
                    dimnames = list(NULL, c("basis", "K", "h_t", "lambda"))),
       BLK = matrix(NA_integer_, nblk, 4L,
                    dimnames = list(NULL, c("b", "v", "seed_out", "seed"))))
}

## dti_resize(): the same storage with room for r_out repetitions, the
## blocks already computed being carried over. The key of block (b, v) is
## (b-1)V + v and V sits in the signature, so the keys are stable.
dti_resize <- function(S, V, r_out) {
  w <- dti_empty(S$sig, V, r_out)
  keep <- seq_len(min(nrow(S$E), nrow(w$E)))
  w$E[keep, ] <- S$E[keep, ]; w$SEL[keep, ] <- S$SEL[keep, ]
  w$BLK[keep, ] <- S$BLK[keep, ]
  w$secs <- S$secs; w$r_sel <- R_SEL
  w
}

## dti_status(): what the stored file holds, without running anything.
dti_status <- function() {
  if (!has_rds(RDS_DTI)) {
    cat(sprintf("%s.rds: absent\n", RDS_DTI)); return(invisible(NULL))
  }
  S <- load_rds(RDS_DTI)
  done <- sum(!is.na(S$E[, "sp"]))
  cat(sprintf("%s.rds: %d / %d blocks done (%.1f min)\n",
              RDS_DTI, done, nrow(S$E), S$secs / 60))
  invisible(S)
}


## edge_report(): how often each tuning was selected at an end of its
## grid, over all the held-out blocks, read off the stored file alone.
## A count on the FIRST value of a grid calls for lowering that value,
## one on the LAST for raising it; 'sel' is the number of selections
## made per parameter, that is the number of blocks scored.
edge_report <- function() {
  S <- load_rds(RDS_DTI)
  nht <- length(S$sig$h_grid_t); nlam <- length(S$sig$lam_grid)
  idx <- which(!is.na(S$E[, "sp"]))
  if (!length(idx)) stop(sprintf("%s.rds holds no completed block yet.", RDS_DTI))
  cnt <- c(ht_lo  = sum(S$SEL[idx, "h_t"]    == S$sig$h_grid_t[1]),
           ht_hi  = sum(S$SEL[idx, "h_t"]    == S$sig$h_grid_t[nht]),
           lam_lo = sum(S$SEL[idx, "lambda"] == S$sig$lam_grid[1]),
           lam_hi = sum(S$SEL[idx, "lambda"] == S$sig$lam_grid[nlam]),
           K_hi   = sum(S$SEL[idx, "K"]      == S$sig$Kmax),
           sel    = length(idx))
  tab <- cbind(count = cnt, pct = round(100 * cnt / length(idx), 1))
  cat(sprintf("grids: h_t [%.3f, %.3f]   lambda [%.1e, %.1e]   Kmax %d\n",
              S$sig$h_grid_t[1], S$sig$h_grid_t[nht],
              S$sig$lam_grid[1], S$sig$lam_grid[nlam], S$sig$Kmax))
  cat("counts of selections at an end of a grid; 'sel' is the number of\n",
      "selections made per parameter (held-out blocks)\n", sep = "")
  print(tab)
  invisible(tab)
}

run_dti <- function(V = V_FOLDS, r_out = R_OUT, r_sel = R_SEL,
                    resume = RESUME, every = CKPT_EVERY) {
  cat("\n[dti] reading the data ...\n")
  dat <- dti_prepare()
  Tv <- dat$T; Y <- dat$Y; subj <- dat$subj; ug <- dat$ug
  n <- dat$n; M <- dat$M; Nu <- dat$Nu
  cat(sprintf("    %d patients (%d seen once, all kept), %d visits, m_i in [%d, %d], N_u = %d\n",
              n, dat$n_one, M, min(dat$mi), max(dat$mi), Nu))
  
  sig <- dti_sig(V)
  co  <- lapply(0:2, basis_r, dim = 2L + KMAX)   # largest dimension needed
  
  ## ---- the stored file, if it is to be continued ------------------------
  ST <- NULL
  if (resume && has_rds(RDS_DTI)) {
    ST <- load_rds(RDS_DTI)
    if (!isTRUE(all.equal(ST$sig, sig)))
      stop(sprintf("%s.rds was produced under other parameters: move it aside, ",
                   RDS_DTI),
           "or call run_dti(resume = FALSE) to overwrite it.")
    ST <- dti_resize(ST, V, r_out)
    cat(sprintf("    resuming: %d / %d blocks already stored\n",
                sum(!is.na(ST$E[, "sp"])), nrow(ST$E)))
  }
  if (is.null(ST)) ST <- dti_empty(sig, V, r_out)
  
  ## ---- selection on the training patients, score on the held-out ones ---
  cat(sprintf("[dti] nested cross-validation: %d x %d held-out blocks ...\n",
              r_out, V))
  ndone <- 0L
  for (b in seq_len(r_out)) for (v in seq_len(V)) {
    q <- (b - 1L) * V + v
    ST$BLK[q, ] <- c(b, v, SEEDOUT(b), SEEDBLK(b, v))
    if (!is.na(ST$E[q, "sp"])) next
    t_blk <- proc.time()[3]
    r <- one_block(dat, co, b, v, V)
    ST$E[q, ] <- r$E; ST$SEL[q, ] <- r$sel
    ST$secs <- ST$secs + (proc.time()[3] - t_blk); ndone <- ndone + 1L
    if (ndone %% every == 0L) save_rds(RDS_DTI, ST, quiet = TRUE)
    cat(sprintf("    block (%d, %2d) seed %d : E_sp %.3e, Delta_ll %+.3f, Delta_spl %+.3f\n",
                b, v, r$seed, r$E["sp"], log(r$E["ll"] / r$E["sp"]),
                log(r$E["spl"] / r$E["sp"])))
  }
  ST$D <- cbind(ll  = log(ST$E[, "ll"]  / ST$E[, "sp"]),
                spl = log(ST$E[, "spl"] / ST$E[, "sp"]))
  save_rds(RDS_DTI, ST, quiet = TRUE)
  
  ## ---- one selection on the whole sample, for the text ------------------
  cat(sprintf("[dti] whole-sample selection over %d partitions ...\n", r_sel))
  set.seed(SEEDSEL)
  ST$sel <- cv_select(Y, Tv, subj, ug, co, V = V, reps = r_sel)
  
  ## ---- the two readings of H_0 : mu(u,t) = m(u) -------------------------
  ## Both are computed on the half-cosine system B_0: testing the
  ## dependence on t calls for no augmentation, and the modes 2..KBAR+1
  ## are then orthogonal to the constants, so the raw reading needs no
  ## centring. The two readings share the multipliers.
  cat("[dti] multiplier bootstrap, centred and raw profiles ...\n")
  co0 <- basis_r(0L, KBAR + 1L)
  set.seed(SEEDBOOT)
  Xi <- draw_xi(LAW, n, B_TEST)
  tc <- dti_test(Y, Tv, subj, n, ug, co0, modes = seq_len(KBAR),
                 centre = TRUE,  Xi = Xi)
  tr <- dti_test(Y, Tv, subj, n, ug, co0, modes = 1L + seq_len(KBAR),
                 centre = FALSE, Xi = Xi)
  check_test(tc$G); check_test(tr$G)
  ST$test <- list(centred = tc[c("Q", "Qstar", "lead", "modes", "centre")],
                  raw     = tr[c("Q", "Qstar", "lead", "modes", "centre")])
  
  ST$data <- dat[c("n", "M", "Nu", "mi", "days", "n_one")]
  save_rds(RDS_DTI, ST)
  if (all(!is.na(ST$E[, "sp"]))) clean_ckpt(RDS_DTI)
  else cat(sprintf("  run still partial: %s.rds.bak kept\n", RDS_DTI))
  
  cat(sprintf("    whole sample: B_%d, K = %d | h_t = %.3f | lambda = %.1e\n",
              ST$sel$basis, ST$sel$K, ST$sel$h_t, ST$sel$lambda))
  cat(sprintf("    median Delta: LL %+.3f (%.0f%% won), spl %+.3f (%.0f%% won)\n",
              median(ST$D[, 1], na.rm = TRUE), 100 * mean(ST$D[, 1] > 0, na.rm = TRUE),
              median(ST$D[, 2], na.rm = TRUE), 100 * mean(ST$D[, 2] > 0, na.rm = TRUE)))
  for (e in c("centred", "raw")) {
    z <- ST$test[[e]]
    cat(sprintf("    test (%-7s, modes %d..%d) : Q = %.3e, q_%.2f = %.3e, p = %.3f%s\n",
                e, z$modes[1], z$modes[2], z$Q, 1 - ALPHA,
                unname(quantile(z$Qstar, 1 - ALPHA)), pval(z$Q, z$Qstar),
                if (z$centre) "" else
                  sprintf("   [level in Q: %.3e, sd* %.3e]",
                          z$lead, sd(z$Qstar))))
  }
  invisible(ST)
}


# =========================================================================
#  12. THE TABLES AND THE MACROS OF THE SECTION
#
#  Everything here is read off res_dti.rds alone, so a change of layout,
#  of caption or of nominal level costs no computation. ALPHA enters HERE
#  and nowhere else: the file stores the statistic and the bootstrap
#  draws, not the rejection.
# =========================================================================

make_dti <- function(alpha = ALPHA) {
  cat("\n[dti] table, macros and block file from res_dti.rds ...\n")
  S <- load_rds(RDS_DTI); sel <- S$sel
  ok <- which(!is.na(S$E[, "sp"]))
  if (!length(ok)) stop("res_dti.rds holds no completed block yet.")
  if (length(ok) < nrow(S$E))
    cat(sprintf("  [i] partial file: %d of %d blocks\n", length(ok), nrow(S$E)))
  E <- S$E[ok, , drop = FALSE]; SEL <- S$SEL[ok, , drop = FALSE]
  BLK <- S$BLK[ok, , drop = FALSE]
  D <- cbind(ll = log(E[, "ll"] / E[, "sp"]), spl = log(E[, "spl"] / E[, "sp"]))
  nblk <- nrow(D)
  
  ## ---- which blocks the spacing estimator wins, and with what seeds -----
  i_ll   <- which(D[, "ll"]  > 0)
  i_spl  <- which(D[, "spl"] > 0)
  i_both <- which(D[, "ll"]  > 0 & D[, "spl"] > 0)
  pc     <- function(i) sprintf("%.0f", 100 * length(i) / nblk)
  cat(sprintf("    smaller error on %d/%d blocks vs LL, %d/%d vs spl, %d/%d vs both\n",
              length(i_ll), nblk, length(i_spl), nblk, length(i_both), nblk))
  
  ## ---- tunings selected at an end of their grid -------------------------
  nht <- length(S$sig$h_grid_t); nlam <- length(S$sig$lam_grid)
  edge <- c(ht_hi  = sum(SEL[, "h_t"]    == S$sig$h_grid_t[nht]),
            ht_lo  = sum(SEL[, "h_t"]    == S$sig$h_grid_t[1]),
            lam_hi = sum(SEL[, "lambda"] == S$sig$lam_grid[nlam]),
            lam_lo = sum(SEL[, "lambda"] == S$sig$lam_grid[1]),
            K_hi   = sum(SEL[, "K"]      == S$sig$Kmax))
  if (any(edge > 0L)) {
    cat(sprintf("  [i] tunings selected at an end of their grid, out of %d blocks:\n",
                nblk))
    txt <- c(ht_hi  = sprintf("h_t at %.3f, LAST of LL_HGRID   -> raise its last value",
                              S$sig$h_grid_t[nht]),
             ht_lo  = sprintf("h_t at %.3f, FIRST of LL_HGRID  -> lower its first value",
                              S$sig$h_grid_t[1]),
             lam_hi = sprintf("lambda at %.1e, LAST of SPL_LGRID  -> raise its last value",
                              S$sig$lam_grid[nlam]),
             lam_lo = sprintf("lambda at %.1e, FIRST of SPL_LGRID -> lower its first value",
                              S$sig$lam_grid[1]),
             K_hi   = sprintf("K at %d, LAST of 1..KMAX          -> raise KMAX",
                              S$sig$Kmax))
    for (k in names(edge)[edge > 0L])
      cat(sprintf("      %4d  (%5.1f%%)   %s\n", edge[[k]],
                  100 * edge[[k]] / nblk, txt[[k]]))
  }
  
  ## ---- dti_blocks.csv : every block, whatever its outcome ---------------
  write.csv(data.frame(block = ok, b = BLK[, "b"], v = BLK[, "v"],
                       seed_out = BLK[, "seed_out"], seed = BLK[, "seed"],
                       basis = SEL[, "basis"], K = SEL[, "K"],
                       h_t = SEL[, "h_t"], lambda = SEL[, "lambda"],
                       E_sp = E[, "sp"], E_ll = E[, "ll"], E_spl = E[, "spl"],
                       Delta_ll = D[, "ll"], Delta_spl = D[, "spl"]),
            file.path(OUTDIR, "dti_blocks.csv"), row.names = FALSE)
  cat("  wrote dti_blocks.csv\n")
  
  ## ---- tab_dti ----------------------------------------------------------
  ex <- floor(log10(max(apply(E, 2, median))))
  f  <- function(x) sprintf("$%.2f$", x / 10^ex)
  bas <- SEL[, "basis"]; Kc <- SEL[, "K"]
  bmod <- as.integer(names(which.max(table(bas))))
  Kmod <- as.integer(names(which.max(table(Kc[bas == bmod]))))
  body <- c(
    paste("$\\widehat\\mu^{\\rm sp}$",
          sprintf("$\\mathcal{B}_%d$ (%.0f\\%%), $K=%d$", bmod,
                  100 * mean(bas == bmod), Kmod),
          f(median(E[, "sp"])), sep = " & "),
    paste("$\\widehat\\mu^{\\rm LL}$",
          sprintf("$\\widehat h_t = %s$", num_math(median(SEL[, "h_t"]))),
          f(median(E[, "ll"])), sep = " & "),
    paste("$\\widehat\\mu^{\\rm spl}$",
          sprintf("$\\widehat\\lambda = %s$", num_math(median(SEL[, "lambda"]))),
          f(median(E[, "spl"])), sep = " & "))
  body <- paste0(body, " \\\\")
  write_table("tab_dti", "llc",
              sprintf(paste0("estimator & tuning selected & $\\mathcal{E}_e$ ",
                             "($\\times 10^{%d}$) \\\\"), ex),
              body,
              sprintf(paste0("Medians over the $%d\\times%d$ held-out blocks of what the ",
                             "criterion~\\eqref{eq:cv} selected on the training patients ",
                             "of each block and of the resulting held-out error ",
                             "$\\mathcal{E}_e$ of~\\eqref{eq:dti_delta}."),
                      S$r_out, S$V),
              "tab:dti")
  
  ## ---- the two readings of the test -------------------------------------
  tst <- lapply(S$test, function(z) list(
    Q = z$Q, q = unname(quantile(z$Qstar, 1 - alpha)), p = pval(z$Q, z$Qstar),
    lead = if (is.null(z$lead)) NA_real_ else z$lead, sd = sd(z$Qstar),
    lo = z$modes[1], hi = z$modes[2]))
  for (e in names(tst))
    cat(sprintf("    test (%-7s, modes %d..%d) : Q = %.3e, q_%.2f = %.3e, p = %.3f   [level in Q: %.3e, sd* %.3e]\n",
                e, tst[[e]]$lo, tst[[e]]$hi, tst[[e]]$Q, 1 - alpha,
                tst[[e]]$q, tst[[e]]$p, tst[[e]]$lead, tst[[e]]$sd))
  
  ## ---- dti_numbers.tex : the macros the section calls, and only those ---
  mac <- function(nm, val) sprintf("\\newcommand{\\%s}{%s}", nm, val)
  writeLines(c(
    mac("dtiN",         sprintf("%d", S$data$n)),
    mac("dtinone",      sprintf("%d", S$data$n_one)),
    mac("dtiM",         sprintf("%d", S$data$M)),
    mac("dtiNu",        sprintf("%d", S$data$Nu)),
    mac("dtimmin",      sprintf("%d", min(S$data$mi))),
    mac("dtimmax",      sprintf("%d", max(S$data$mi))),
    mac("dtidmin",      sprintf("%d", as.integer(min(S$data$days)))),
    mac("dtidmax",      sprintf("%d", as.integer(max(S$data$days)))),
    mac("dtiV",         sprintf("%d", S$V)),
    mac("dtiRout",      sprintf("%d", S$r_out)),
    mac("dtiblk",       sprintf("%d", nblk)),
    mac("dtiwinll",     sprintf("%d", length(i_ll))),
    mac("dtiwinspl",    sprintf("%d", length(i_spl))),
    mac("dtiwinboth",   sprintf("%d", length(i_both))),
    mac("dtipctll",     pc(i_ll)),
    mac("dtipctspl",    pc(i_spl)),
    mac("dtipctboth",   pc(i_both)),
    mac("dtimedll",     num_math(median(D[, "ll"]))),
    mac("dtimedspl",    num_math(median(D[, "spl"]))),
    mac("dtiseed",      sprintf("%d", S$seed)),
    mac("dtibasis",     sprintf("\\mathcal{B}_%d", sel$basis)),
    mac("dtiK",         sprintf("%d", sel$K)),
    mac("dtiKbar",      sprintf("%d", S$kbar)),
    mac("dtiKmlo",      sprintf("%d", tst$centred$lo)),
    mac("dtiKmhi",      sprintf("%d", tst$centred$hi)),
    mac("dtiKmlor",     sprintf("%d", tst$raw$lo)),
    mac("dtiKmhir",     sprintf("%d", tst$raw$hi)),
    mac("dtiNb",        sprintf("%d", S$B)),
    mac("dtilevel",     sprintf("%.2f", 1 - alpha)),
    mac("dtipm",        sprintf("%.3f", tst$centred$p)),
    mac("dtip",         sprintf("%.3f", tst$raw$p)),
    mac("dtiQstatm",    num_math(tst$centred$Q)),
    mac("dtiQstat",     num_math(tst$raw$Q)),
    mac("dtibquantilem", num_math(tst$centred$q)),
    mac("dtibquantile",  num_math(tst$raw$q))),
    file.path(OUTDIR, "dti_numbers.tex"))
  cat("  wrote dti_numbers.tex\n")
  
  ## ---- the seeds of the blocks won, printed here and nowhere else ------
  cat("\n[dti] seeds of the blocks on which the spacing estimator has the",
      "smaller held-out error\n")
  cat("    vs LL   (", length(i_ll),   "/", nblk, ") : ",
      seed_list(BLK[i_ll,   "seed"]), "\n", sep = "")
  cat("    vs spl  (", length(i_spl),  "/", nblk, ") : ",
      seed_list(BLK[i_spl,  "seed"]), "\n", sep = "")
  cat("    vs both (", length(i_both), "/", nblk, ") : ",
      seed_list(BLK[i_both, "seed"]), "\n", sep = "")
  cat("    every block, won or lost, is in dti_blocks.csv\n")
  
  invisible(list(S = S, D = D, test = tst, edge = edge))
}

# =========================================================================
#  13. FIGURES  --  observed profiles and fitted mean, one panel per
#                    estimator, at the tuning of tab_dti
#
#  Reads only res_dti.rds (for S$sel, the tuning cv_select() picked on
#  the whole sample) and the DTI data (dti_prepare(), not kept in the
#  .rds to keep it light). No cross-validation is redone here: each
#  estimator is fitted ONCE, at ITS OWN tuning of tab_dti.
# =========================================================================

## ---- 13.1 graphics machinery, as in simul_estimation.R -------------------

TEX_HEAD <- c("\\documentclass[border=2pt]{standalone}",
              "\\usepackage[T1]{fontenc}",
              "\\usepackage{amsmath,amssymb}",
              "\\usepackage{tikz}",
              "\\begin{document}")
TEX_FOOT <- "\\end{document}"

## check_graphics(): the figure toolchain is needed by make_dti_profiles()
## only, so that run_dti()/make_dti() run on a machine without LaTeX.
GRAPHICS_OK <- FALSE
check_graphics <- function() {
  if (GRAPHICS_OK) return(invisible(TRUE))
  if (!requireNamespace("tikzDevice", quietly = TRUE))
    stop("package tikzDevice is required to write the figures")
  if (Sys.which(LATEX) == "")
    stop(sprintf("'%s' is not on the PATH: it compiles each figure to PDF", LATEX))
  options(tikzDefaultEngine = "pdftex",
          tikzMetricsDictionary = file.path(OUTDIR, "tikzMetrics"),
          tikzLatexPackages = c(getOption("tikzLatexPackages"),
                                "\\usepackage{amsmath}\n", "\\usepackage{amssymb}\n"))
  ok <- tryCatch({                               # one tiny metric test
    f <- tempfile(fileext = ".tex")
    tikzDevice::tikz(f, width = 2, height = 2); plot.new(); text(.5, .5, "$x_1$"); TRUE
  }, error = function(e) FALSE)
  while (length(dev.list())) try(dev.off(), silent = TRUE)
  if (!isTRUE(ok))
    stop("tikzDevice cannot compute LaTeX metrics here: check the LaTeX ",
         "installation on the PATH")
  GRAPHICS_OK <<- TRUE
  invisible(TRUE)
}

## compile_tex(): run LATEX on OUTDIR/<name>.tex and keep only the PDF.
compile_tex <- function(name) {
  owd <- setwd(OUTDIR); on.exit(setwd(owd), add = TRUE)
  out <- suppressWarnings(
    system2(LATEX, c("-interaction=nonstopmode", "-halt-on-error",
                     shQuote(paste0(name, ".tex"))),
            stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0L) {
    err <- grep("^(!|l\\.[0-9])", out, value = TRUE)
    message(sprintf("  [!] %s failed on %s.tex; see %s.log", LATEX, name, name))
    message(paste(if (length(err)) utils::head(err, 6L) else utils::tail(out, 12L),
                  collapse = "\n"))
    unlink(paste0(name, ".aux"))
    return(invisible(FALSE))
  }
  unlink(paste0(name, c(".aux", ".log")))
  invisible(TRUE)
}

## emit_fig(): draw 'plotfun' into OUTDIR/<name>.tex, as ONE panel wrapped
## in the standalone preamble above, then compile it into <name>.pdf.
emit_fig <- function(name, width, height, plotfun) {
  check_graphics()
  tmp <- tempfile(fileext = ".tex")
  tikzDevice::tikz(tmp, width = width, height = height, standAlone = FALSE)
  ok <- tryCatch({ plotfun(); TRUE },
                 error = function(e) { message("  [!] ", conditionMessage(e)); FALSE })
  while (length(dev.list())) try(dev.off(), silent = TRUE)
  if (!ok) return(invisible(FALSE))
  out <- file.path(OUTDIR, paste0(name, ".tex"))
  writeLines(c(TEX_HEAD, readLines(tmp, warn = FALSE), TEX_FOOT), out)
  if (!compile_tex(name)) {
    cat(sprintf("  wrote %s.tex (%.0f kB), NOT compiled\n", name, file.size(out) / 1024))
    return(invisible(FALSE))
  }
  cat(sprintf("  wrote %s.tex (%.0f kB) and %s.pdf (%.0f kB)\n",
              name, file.size(out) / 1024,
              name, file.size(file.path(OUTDIR, paste0(name, ".pdf"))) / 1024))
  invisible(TRUE)
}

## ---- 13.2 the three fitted curves, at a single tuning ---------------------
##
## ll_crit()/spl_crit() of Section 5 only ever return a SUM OF SQUARES
## over a grid of tunings, for the criterion; beta_spacing() already
## gives the spacing fit. What is missing is the FIT of the two
## benchmarks at a SINGLE tuning: ll_fit() and spl_fit() are the same
## normal equations as ll_crit()/spl_crit(), specialised to one
## bandwidth or one penalty, returning the fitted matrix itself.

## sp_fit(): mu-hat^sp(u, t_eval), at the basis co and window half-width
## h, directly from beta_spacing() and basis_phi().
sp_fit <- function(Y, Tv, co, h, t_eval) beta_spacing(Y, Tv, co, h) %*% t(basis_phi(co, t_eval))

## ll_fit(): the pooled local-linear smoother of CM2012 in the
## visit-time direction, at a SINGLE bandwidth h, evaluated at every
## (u_l, t_eval). Same normal equations as ll_crit(), stat = "mean",
## with hgrid reduced to the one value h.
ll_fit <- function(Y, Tp, ug, t_eval, h) {
  Nu <- length(ug)
  Dt <- outer(Tp, t_eval, "-")
  Kt <- pmax(1 - (Dt / h)^2, 0)                        # Epanechnikov in t
  T0 <- colSums(Kt); T1 <- colSums(Kt * Dt); T2 <- colSums(Kt * Dt * Dt)
  A0 <- Y %*% Kt; A1 <- Y %*% (Kt * Dt)
  den <- T0 * T2 - T1 * T1
  Fh <- (A0 * rep(T2, each = Nu) - A1 * rep(T1, each = Nu)) /
    rep(pmax(den, 1e-300), each = Nu)
  bad <- which(den <= 1e-10)
  if (length(bad)) Fh[, bad] <- (A0 / rep(pmax(T0, 1e-300), each = Nu))[, bad]
  empty <- which(T0 <= 1e-12)
  if (length(empty)) Fh[, empty] <- rowMeans(Y)
  Fh
}

## spl_fit(): penalised cubic B-spline smoother in t, fitted slice-wise
## in u, at a SINGLE penalty lambda, evaluated at t_eval. Same normal
## system as spl_crit(), specialised to one lambda.
spl_fit <- function(Y, Tp, t_eval, lambda, nik = SPL_NIK) {
  kn <- seq(0, 1, length.out = nik + 2L)[-c(1L, nik + 2L)]
  Bm <- splines::bs(Tp, knots = kn, degree = 3L, intercept = TRUE,
                    Boundary.knots = c(0, 1))
  Bt <- predict(Bm, t_eval); q <- ncol(Bm); M <- length(Tp)
  Pen <- crossprod(diff(diag(q), differences = 2L))
  A   <- chol2inv(chol(crossprod(Bm) + lambda * M * Pen + 1e-10 * diag(q)))
  Y %*% Bm %*% A %*% t(Bt)
}

## ---- 13.3 the three panels -------------------------------------------------

## pick_t3(): three representative visit times, the 10th, 50th and 90th
## percentile of the pooled, rescaled visit times -- inside the support
## by construction, so no panel extrapolates beyond the observed design.
pick_t3 <- function(Tv) unname(quantile(Tv, c(0.10, 0.50, 0.90), type = 7))

LTY3 <- c(3, 2, 1)   # dotted, dashed, solid: earliest, middle, latest t

## profile_panel(): one panel -- every observed profile Y_ij(.) in light
## grey against u, and the fitted curve mu-hat_e(u, t_j) overlaid in
## black for each of the three t_j of 'legend_t', dotted/dashed/solid
## from the earliest to the latest. 'Fh' is the N_u x 3 matrix of fitted
## values, one column per t_j, in that same order.
profile_panel <- function(name, ug, Y, Fh, legend_t, title,
                          width = 3.1, height = 3.0) {
  emit_fig(name, width, height, function() {
    par(mar = c(3.2, 3.6, 1.3, 0.5) * c(CEXF, CEXF, 1, 1),
        mgp = c(2.0, 0.6, 0) * CEXF, tcl = -0.3,
        cex.axis = 0.82 * CEXF, cex.lab = CEXF, cex.main = 0.95 * CEXF)
    matplot(ug, Y, type = "l", lty = 1, col = "grey82", lwd = 0.35,
            xlab = "$u$", ylab = "fractional anisotropy", main = title,
            ylim = range(Y, Fh), bty = "l")
    for (j in 1:3)
      lines(ug, Fh[, j], lty = LTY3[j], lwd = 1.5 * CEXF, col = "black")
    legend("topright", legend = sprintf("$t=%.2f$", legend_t), lty = LTY3,
           lwd = 1.3, seg.len = 2.2, bty = "n", cex = 0.72 * CEXF)
  })
}

## make_dti_profiles(): the three panels, each estimator fitted once on
## the whole sample at ITS OWN tuning of tab_dti -- S$sel, which
## cv_select() already computed in run_dti() and which res_dti.rds
## already stores, so nothing is re-selected here.
make_dti_profiles <- function() {
  cat("\n[dti] fitted-profile panels from res_dti.rds and the DTI data ...\n")
  S <- load_rds(RDS_DTI); sel <- S$sel
  if (is.null(sel))
    stop("res_dti.rds has no whole-sample selection yet: run run_dti() first.")
  
  dat <- dti_prepare()
  Tv <- dat$T; Y <- dat$Y; ug <- dat$ug; M <- dat$M
  t3 <- pick_t3(Tv)
  h  <- h_of_M(M)
  
  ## mu^sp, at the basis and truncation of the paper's notation,
  ## B_{r,K} = basis_r(r, r + K), the r monomials excluded from K
  co_sp <- basis_r(sel$basis, sel$basis + sel$K)
  Fh_sp <- sp_fit(Y, Tv, co_sp, h, t3)
  profile_panel("fig_dti_sp", ug, Y, Fh_sp, t3,
                sprintf("$\\widehat\\mu^{\\rm sp}$ ($\\mathcal{B}_%d$, $K=%d$)",
                        sel$basis, sel$K))
  
  Fh_ll <- ll_fit(Y, Tv, ug, t3, sel$h_t)
  profile_panel("fig_dti_ll", ug, Y, Fh_ll, t3,
                sprintf("$\\widehat\\mu^{\\rm LL}$ ($\\widehat h_t=%s$)",
                        num_math(sel$h_t)))
  
  Fh_spl <- spl_fit(Y, Tv, t3, sel$lambda)
  profile_panel("fig_dti_spl", ug, Y, Fh_spl, t3,
                sprintf("$\\widehat\\mu^{\\rm spl}$ ($\\widehat\\lambda=%s$)",
                        num_math(sel$lambda)))
  
  cat(sprintf("    t_1, t_2, t_3 (10th/50th/90th pct. of the pooled visit times) = %.3f, %.3f, %.3f\n",
              t3[1], t3[2], t3[3]))
  invisible(list(t3 = t3, Fh_sp = Fh_sp, Fh_ll = Fh_ll, Fh_spl = Fh_spl))
}

# =========================================================================
#  14. DRIVER
#
#  Sourcing the script with SIM_NO_RUN defined loads the functions and
#  runs nothing, so that make_dti()/make_dti_profiles() can be called
#  alone on a stored file, or one block replayed with replay_block(b, v).
# =========================================================================

if (!exists("SIM_NO_RUN")) {
  t0 <- proc.time()[3]
  cat("==== basis self-tests ====\n"); check_bases()
  if (has_rds(RDS_DTI)) { cat("==== checkpoint found ====\n"); dti_status() }
  run_dti(); make_dti(); make_dti_profiles()
  cat("selections at an end of a tuning grid:\n")
  edge_report()
  cat("==== replay self-test ====\n")
  check_blocks(which = c(1L, R_OUT * V_FOLDS))
  cat(sprintf("\nTotal wall time: %.1f s\nOutputs written to: %s\n",
              proc.time()[3] - t0, normalizePath(OUTDIR)))
}