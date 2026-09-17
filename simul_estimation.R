##########################################################################
#            MONTE CARLO STUDY -- Estimating the mean function
#
#   TWO EXPERIMENTS, each with its OWN .rds and its OWN figures. They share
#   the parameters, the mechanism and the criterion of this file. Both can
#   be run, resumed, extended or thrown away without touching the other's file.
#
#   [1] THE COMPETITORS. Every estimator is tuned BY THE PRACTITIONER'S OWN
#       MEANS: the pair (basis, truncation) of the spacing estimator, the
#       bandwidth of the local-linear smoother and the penalty of the
#       P-spline are selected on each sample by the same criterion, the
#       subject-level V-fold criterion of the paper, evaluated on ONE
#       partition of the subjects shared by the three families. No
#       estimator receives a tuning that minimises a Monte Carlo estimate
#       of its risk. The three means are treated on the same draws.
#           run_cvbox()  -> res_est.rds
#           make_cvbox() -> fig_cvbox_mu1, fig_cvbox_mu2, fig_cvbox_mu3,
#                           tab_cvbox
#
#   [2] BASES AND WEIGHTING SCHEMES, on mu_1 only.
#       (i)  what the augmentation buys: the spacing estimator on B_{0,K}
#            and on B_{1,K}, the truncation being selected inside each
#            system, against the same estimator on B_{2,K};
#       (ii) what the design density costs: the fully data-driven mu^sp,
#            which never uses g, against three estimators built on the
#            ratio phi_k(T_ij)/g(T_ij), namely the plug-in mu^OBS_ghat,
#            whose Gaussian-kernel bandwidth is FIXED at b = M^{-1/3},
#            and the INFEASIBLE mu^OBS and mu^MC, which are given the true
#            g.
#           run_schemes()  -> res_schemes.rds
#           make_schemes() -> fig_basis_mu1, tab_basis,
#                             fig_dens_mu1,  tab_dens
#
#       [2] opens each replication with the seed of [1] and rebuilds its
#       fold partition, so the design, the subject fields, the errors and
#       the blocks are the same and the two sets of boxes are paired. The
#       reference mu^sp of (ii) is therefore not refitted: make_schemes()
#       READS it from res_est.rds, which the run itself never opens.
#
#   Every experiment stores EVERYTHING it computed, so that a change of
#   style, of scale, of caption or of layout costs no simulation at all:
#   the make_* routines read the file and nothing else.
#
#   CHECKPOINTS AND RESUMPTION. Each .rds is rewritten every CKPT_EVERY
#   replications, not once at the end, so that a machine that stops loses
#   at most CKPT_EVERY replications. Three points make this safe and exact:
#     (i)   the write is atomic: the object goes to <name>.rds.tmp, the
#           current file is renamed <name>.rds.bak, and the temporary file
#           then takes its place, so a cut always leaves one COMPLETE file
#           behind, never a half-written one;
#     (ii)  each replication opens with its OWN seed, rep_seed(n, rep), so
#           the replication that follows a restart is the very one an
#           uninterrupted run would have produced: the results depend
#           neither on where the run was cut nor on the order of the
#           sample sizes;
#     (iii) each stored slot carries R_done, the number of replications
#           completed at that n, so the make_* routines run on a partial
#           file and the figures can be looked at while the study is still
#           running.
#   Restarting with RESUME = TRUE picks each file up where it stopped.
#   R_REP may be raised between two runs (the replications already stored
#   are kept and the new ones appended) and N_GRID may be extended (a
#   sample size already computed is not recomputed, and one dropped from
#   the grid is kept in the file). Once every replication asked for is in,
#   and only then, the file is read back in full and the two auxiliary
#   files .rds.tmp and .rds.bak are removed.
#
#   FIGURES. One panel per file, each file a self-contained standalone
#   LaTeX document (\documentclass[border=2pt]{standalone}) COMPILED on the
#   spot, so that every panel comes as a .tex kept for editing and as a
#   .pdf included by \includegraphics.
#
#   NOTATION OF THE PAPER. The estimator is mu^sp of the spacing scheme,
#   expanded on the candidate bases B_{r,K}, r = 0, 1, 2 (half-cosine,
#   augmented by t, augmented by t and t^2). THROUGHOUT THIS FILE, AS IN
#   THE PAPER, K COUNTS THE COSINE TERMS c_1, ..., c_K, THE r MONOMIALS OF
#   THE AUGMENTATION EXCLUDED, so that B_{r,K} has dimension r + K. The
#   criterion searches K over 1..Kmax in each system, that is dimensions
#   r+1 .. r+Kmax, and K is what the storage and the tables report. The
#   selected tunings are (r-hat, K_CV)
#   for the spacing estimator, h_t-hat for the local-linear smoother and
#   lambda-hat for the spline; they are collected in tab_cvbox, next to the
#   truncation K_n of the rate rule. The local-linear smoother is fitted in
#   the visit-time direction only, as in the paper: no smoothing in u is
#   needed, since each profile is observed over the whole of U.
#
#   THIRD COMPETITOR OF [1]. The same spacing estimator, but with the basis
#   FIXED to the augmented B_2 and the truncation fixed by the rate rule
#       r + K_n = ceil( C (lambda_m n)^{1/5} ),  i.e. C (7n)^{1/5} by
#   default with r = R_BASIS, so the rule fixes the dimension of the fit,
#   with C = C_RATE chosen BY HAND. It is the only member of the comparison
#   that is not data-driven, and it isolates the cost of cross-validating
#   (basis, K) against a rate-tuned choice.
##########################################################################

## ---- shape parameters of the three mean functions ----------------------
a  <- 1                     # coefficient of t
b  <- -7                    # coefficient of t^2
cc <- 4                     # frequency multiplier (the "c" of the paper)
P  <- 8L                    # number of Fourier terms

## ---- design of the experiments -----------------------------------------
N_GRID   <- c(100L, 200L, 500L)    # numbers of subjects
R_REP    <- 2000L                  # Monte Carlo replications per design point
lambda_m <- 7L                     # Poisson mean of the visit counts
m_cap    <- 15L                    # visit counts conditioned on [2, m_cap]
tau      <- 0.5                    # noise standard deviation (homoscedastic)

## ---- tuning grids -------------------------------------------------------
## The grids are searched by the SAME criterion on the SAME folds.
## LL_HGRID is the expensive one: the cost of [1] is proportional to its
## length.
Kmax      <- 25L                   # largest K considered, the r monomials
# excluded: dimensions r+1 .. r+Kmax
V_FOLDS   <- 10L                   # folds of the subject-level criterion
LL_HGRID  <- exp(seq(log(0.05), log(0.60), length.out = 20))   # h_t
SPL_NIK   <- 20L                                               # spline knots
SPL_LGRID <- exp(seq(log(1e-5), log(1e1), length.out = 20))    # penalties

## ---- rate-tuned spacing competitor of [1] -------------------------------
## Basis fixed to B_2, truncation fixed by the rate rule; C_RATE is the
## hand-chosen constant C of the paper. lambda_m n is the expected number
## of pooled observations, equal to 7n with the default lambda_m. The rule
## fixes the DIMENSION, R_BASIS + K_n = ceil(C (lambda_m n)^{1/5}), so the
## fitted system is B_{R_BASIS, K_n} in the sense of the paper.
C_RATE  <- 2
R_BASIS <- 2L                      # the basis that competitor is expanded on
K_rate  <- function(n)
  max(1L, as.integer(ceiling(C_RATE * (lambda_m * n)^0.2)) - R_BASIS)

## ---- the plug-in and MC schemes of [2] ----------------------------------
## b_of_M is the bandwidth of the kernel estimator of g. It is FIXED, not
## cross-validated: only the pair (r, K) of the coefficient estimators is
## selected by the criterion. The rule undersmooths, the MISE-optimal order
## for a Gaussian kernel being M^{-1/5}, so that the smoothing bias of ghat
## stays negligible against the accuracy of the coefficients, at the price
## of a larger variance. Undersmoothing matters here: g runs over a full
## period of a cosine on [0,1], so a bandwidth of order M^{-1/4} already
## flattens it and the bias goes straight through the ratio phi_k/ghat.
## GH_FLOOR floors ghat away from zero (the true g is bounded below by 0.4).
b_of_M   <- function(M) M^(-1/3)
GH_FLOOR <- 0.05

## ---- evaluation grid ----------------------------------------------------
NU_EVAL <- 50L; NT_EVAL <- 50L     # grid on which every ISE is computed, and
# on which the profiles are recorded

## ---- checkpoints --------------------------------------------------------
## CKPT_EVERY is the number of replications between two writes; at 1 the
## file is up to date after every replication and nothing is ever lost, at
## the price of one small write per replication. RESUME = FALSE ignores an
## existing file and overwrites it.
CKPT_EVERY <- 1L
RESUME     <- TRUE

## ---- what the driver runs -----------------------------------------------
RUN_COMP <- TRUE                   # experiment [1]
RUN_SCH  <- TRUE                   # experiment [2]

## ---- typography of the panels -------------------------------------------
## Multiplier applied to every cex of the figures: tick labels, axis titles
## and legends, together with the margins that must follow them.
CEX <- 1.35

## ---- output -------------------------------------------------------------
OUTDIR    <- Sys.getenv("SIM_OUTDIR", unset = "./output")
LATEX     <- Sys.getenv("SIM_LATEX",  unset = "pdflatex")   # compiles the panels
BASESEED  <- 2026L

## One file per experiment, independent of one another.
RDS_COMP <- "res_est"              # [1] -> fig_cvbox_mu1/mu2/mu3, tab_cvbox
RDS_SCH  <- "res_schemes"          # [2] -> fig_basis_mu1, fig_dens_mu1

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# =========================================================================
#  MEAN FUNCTIONS
#
#  Each function maps two grids (u, t) to the matrix mu(u[i], t[j]).
#      mu_1(u,t) = (a t + b t^2) sum_p ((-1)^p / p^2) cos(pi p t + c pi u)
#      mu_2(u,t) = (a t + b t^2) sum_p ((-1)^p / p^2) cos(pi p u + c pi t)
#      mu_3(u,t) = a t + b t^2 + cos(c pi t) sum_p ((-1)^p / p^2) cos(pi p u)
# =========================================================================

QUAD_T <- function(t) a * t + b * t^2                  # the factor a t + b t^2

MU_1 <- function(u, t) {
  S <- matrix(0, length(u), length(t))
  for (p in 1:P) S <- S + ((-1)^p / p^2) * cos(outer(cc * pi * u, pi * p * t, "+"))
  S * matrix(QUAD_T(t), length(u), length(t), byrow = TRUE)
}
MU_2 <- function(u, t) {
  S <- matrix(0, length(u), length(t))
  for (p in 1:P) S <- S + ((-1)^p / p^2) * cos(outer(pi * p * u, cc * pi * t, "+"))
  S * matrix(QUAD_T(t), length(u), length(t), byrow = TRUE)
}
MU_3 <- function(u, t) {
  Su <- rowSums(sapply(1:P, function(p) ((-1)^p / p^2) * cos(pi * p * u)))
  matrix(QUAD_T(t), length(u), length(t), byrow = TRUE) +
    outer(Su, cos(cc * pi * t))
}

MEANS <- list(
  mu1 = list(f = MU_1, tex = "$\\mu_1$"),
  mu2 = list(f = MU_2, tex = "$\\mu_2$"),
  mu3 = list(f = MU_3, tex = "$\\mu_3$"))


# =========================================================================
#  DATA-GENERATING PROCESS
# =========================================================================

## Visit-time density on [0,1]: Lipschitz, bounded away from zero.
G_DENSITY <- function(t) 1 + 0.6 * cos(2 * pi * t)   # g(t) = 1 + 3cos(2 pi t)/5
G_CDF     <- function(t) t + 0.3 * sin(2 * pi * t) / pi   # its cdf on [0,1]
G_MAX     <- 1.6                                   # any upper bound of g
L_MODES   <- 100L                                  # modes of the subject field
L_NOISE   <- 50L                                   # terms of the error field

## sample_T(): m visit times with density g, by rejection sampling.
sample_T <- function(m) {
  out <- numeric(0)
  while (length(out) < m) {
    x <- runif(2 * m); acc <- runif(2 * m)
    out <- c(out, x[acc <= G_DENSITY(x) / G_MAX])
  }
  sort(out[1:m])
}

## make_Xi(): ONE random draw of the subject field, returned as a function
## (u, t) -> X(u[i], t[j]). Calling it again gives an independent subject.
make_Xi <- function(L_modes = L_MODES) {
  C <- matrix(rnorm(L_modes^2), L_modes, L_modes) /
    outer(seq_len(L_modes), seq_len(L_modes), function(l, lp) l^2 * lp^2)
  function(u, t)
    outer(u, seq_len(L_modes), function(u, l) cos(l * pi * u)) %*% C %*%
    t(outer(t, seq_len(L_modes), function(t, lp) cos(lp * pi * t)))
}

## gen_noise(): 'nc' error curves on the grid u, with unit variance at every
## u; the scale tau is applied in gen_data().
gen_noise <- function(u, nc, Lp = L_NOISE) {
  B  <- matrix(rnorm(Lp * nc), Lp, nc) / (seq_len(Lp)^4)
  Cu <- outer(u, seq_len(Lp), function(u, l) cos(l * pi * u))
  v  <- as.vector((Cu^2) %*% (seq_len(Lp)^(-8)))
  (Cu %*% B) / sqrt(v)
}

## gen_data(): one data set of n subjects, the visit counts being Poisson
## conditioned on [2, m_cap], WITHOUT the mean: both experiments draw the
## design, the subject fields and the noise ONCE per replication and add
## the mean afterwards, so that the comparisons are paired and differ only
## through mu. Returns the pooled visit times, the (length(ug)) x M matrix
## of centred profiles and the subject label of each observation.
gen_data <- function(n, ug) {
  Tp <- vector("list", n); Yc <- vector("list", n); subj <- integer(0)
  for (i in seq_len(n)) {
    mi <- 0L
    while (mi < 2L || mi > m_cap) mi <- rpois(1L, lambda_m)   # 2 <= m_i <= m_cap
    Ti <- sample_T(mi)
    Yc[[i]] <- make_Xi()(ug, Ti) + tau * gen_noise(ug, mi)
    Tp[[i]] <- Ti; subj <- c(subj, rep(i, mi))
  }
  list(T = unlist(Tp), Y = do.call(cbind, Yc), subj = subj)
}

## rep_seed(): the seed of replication 'rep' at sample size 'n'. It depends
## on nothing else, so a replication computed after a restart, in another
## order, in a later run that only extends R_REP, or in the OTHER
## experiment, is the very one an uninterrupted run would have produced.
rep_seed <- function(n, rep) as.integer(BASESEED + 14L + 100000L * n + rep)


# =========================================================================
#  CANDIDATE BASES  --  B_{r,K}, r = 0 (half-cosine), 1 (+ t), 2 (+ t, t^2)
#
#  basis_r(r, dim) orthonormalises, by closed-form Gram-Schmidt, the first
#  'dim' elements of the ordered dictionary
#      r = 0 : (c_1, c_2, c_3, ...)
#      r = 1 : (c_1, t, c_2, c_3, ...)
#      r = 2 : (c_1, t, t^2, c_2, c_3, ...)
#  where c_1 = 1 and c_k(t) = sqrt2 cos((k-1) pi t). Its argument is a
#  DIMENSION, and the paper's B_{r,K} is basis_r(r, r + K): the first
#  r + K elements are c_1, the r monomials, and c_2, ..., c_K. Everywhere
#  else in this file the tuning is the K of the paper, and the conversion
#  is the single formula
#        dimension = r + K.
#  Each phi_k is a polynomial of degree at most 2 plus finitely many
#  cosines, so that both its values and its window integrals are
#  closed-form.
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

## basis_phi(): phi_1..phi_K evaluated at t. Returns length(t) x K.
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


# =========================================================================
#  THE SPACING ESTIMATOR
#
#  beta^sp_k(u) = sum_l Y_(l)(u) Phi_kl, with
#  Phi_kl = (1/2h) int_{T_(l-h)}^{T_(l+h)} phi_k and the window endpoints
#  obtained by symmetric reflection at the two ends of the pooled design.
#  It uses NO design density: the geometry of the order statistics
#  reproduces the correction 1/g from the data alone.
# =========================================================================

## h_of_M(): the window half-width, h = ceil((1 + ln ln(M+20))/2); it grows
## so slowly that h = 2 for every realistic sample size.
h_of_M <- function(M) ceiling(0.5 + 0.5 * log(log(M + 20)))

beta_spacing <- function(Y, Tp, co, h) {
  M <- length(Tp); ord <- order(Tp); Ts <- Tp[ord]; Yo <- Y[, ord, drop = FALSE]
  Text <- function(rk) {                         # reflected order statistics
    out <- numeric(length(rk)); lo <- rk <= 0; hi <- rk > M; in_ <- !lo & !hi
    out[in_] <- Ts[rk[in_]]
    out[lo]  <- 2 * Ts[1] - Ts[2 - rk[lo]]
    out[hi]  <- 2 * Ts[M] - Ts[2 * M - rk[hi]]
    out
  }
  l <- seq_len(M)
  Yo %*% (basis_int(co, Text(l - h), Text(l + h)) / (2 * h))
}

## ise_curve(): integrated squared error of a series fit as a function of
## the DIMENSION 1..Kup of the fit, given Phi_t = basis_phi(co, tg). One
## pass, by rank-one updates; entry r + K is the fit on B_{r,K}.
ise_curve <- function(beta, Phi_t, MUtrue, Kup) {
  out <- numeric(Kup); Fh <- matrix(0, nrow(MUtrue), ncol(MUtrue))
  for (K in seq_len(Kup)) {
    Fh <- Fh + outer(beta[, K], Phi_t[, K])
    out[K] <- mean((Fh - MUtrue)^2)
  }
  out
}

## sse_curve(): held-out sum of squares of the cross-validation criterion,
## as a function of the DIMENSION of the fit, entry r + K being the
## candidate B_{r,K}. The factor 1/(N_u M) is common to all candidates and
## is dropped.
sse_curve <- function(beta, Phi_te, Y_te, Kup) {
  out <- numeric(Kup); Rk <- Y_te
  for (K in seq_len(Kup)) {
    Rk <- Rk - outer(beta[, K], Phi_te[, K])
    out[K] <- sum(Rk * Rk)
  }
  out
}


# =========================================================================
#  THE WEIGHTED ESTIMATORS OF [2]
#
#  All three estimate beta_k(u) = E[Y(u) phi_k(T)/g(T)] by
#      beta_k(u) = sum_l w_l Y_l(u) phi_k(T_l) / g_l ,
#  and differ only through the pair (w, g): OBS weights with the true g,
#  OBS weights with a kernel plug-in ghat, and the control-neighbours MC
#  weights with the true g. Only the plug-in is feasible.
# =========================================================================

beta_weighted <- function(Y, Tp, co, w, gv) Y %*% (basis_phi(co, Tp) * (w / gv))

## ghat_at(): kernel estimator of the design density, Gaussian kernel
## reflected at the two ends of T = [0,1] so that the boundary bias does not
## contaminate the ratio, and floored at GH_FLOOR.
ghat_at <- function(Tp, teval, b) {
  v <- rowSums(dnorm(outer(teval, Tp, "-") / b) +          # sample
                 dnorm(outer(teval, -Tp, "-") / b) +         # reflection at 0
                 dnorm(outer(teval, 2 - Tp, "-") / b))       # reflection at 1
  pmax(v / (length(Tp) * b), GH_FLOOR)
}

## mc_weights_sorted(): the control-neighbours weights, in closed form on
## the line, for a SORTED sample. w_l = (1 + chat_l - dhat_l)/M, with dhat_l
## the number of points whose leave-one-out nearest neighbour is T_(l), and
## chat_l the cumulative Voronoi volume, a volume being the integral of g
## over the cell. In dimension one the cells are the intervals cut at the
## midpoints, and deleting T_(j) moves ONE boundary of each of its two
## neighbours, so chat is (M-1) times the full-sample volume plus two
## corrections. By construction sum(chat) = sum(dhat) = M, hence sum(w) = 1.
mc_weights_sorted <- function(Ts) {
  M <- length(Ts)
  if (M < 4L) return(rep(1 / M, M))
  bnd <- c(0, (Ts[-M] + Ts[-1]) / 2, 1)                 # cell boundaries
  vol <- G_CDF(bnd[-1]) - G_CDF(bnd[-(M + 1L)])
  dR <- numeric(M); dL <- numeric(M)
  ir <- seq_len(M - 1L)                                 # delete j = l+1
  newR <- ifelse(ir == M - 1L, 1, (Ts[ir] + Ts[pmin(ir + 2L, M)]) / 2)
  dR[ir] <- G_CDF(newR) - G_CDF(bnd[ir + 1L])
  il <- 2:M                                             # delete j = l-1
  newL <- ifelse(il == 2L, 0, (Ts[pmax(il - 2L, 1L)] + Ts[il]) / 2)
  dL[il] <- G_CDF(bnd[il]) - G_CDF(newL)
  chat <- (M - 1) * vol + dR + dL
  nb <- integer(M); nb[1] <- 2L; nb[M] <- M - 1L        # leave-one-out 1-NN
  j <- 2:(M - 1L)
  nb[j] <- ifelse(Ts[j] - Ts[j - 1L] <= Ts[j + 1L] - Ts[j], j - 1L, j + 1L)
  (1 + chat - tabulate(nb, M)) / M
}

## weight_set(): the three (w, g) pairs of one sample, training or full.
weight_set <- function(Tp) {
  M <- length(Tp); bw <- b_of_M(M)
  gh <- ghat_at(Tp, Tp, bw)
  ord <- order(Tp); wmc <- numeric(M); wmc[ord] <- mc_weights_sorted(Tp[ord])
  list(gh  = list(w = rep(1 / M, M), g = gh),
       obs = list(w = rep(1 / M, M), g = G_DENSITY(Tp)),
       mc  = list(w = wmc,           g = G_DENSITY(Tp)),
       b = bw, floor = mean(gh <= GH_FLOOR + 1e-12))
}

## check_schemes(): the identity the MC weights rest on, their quadrature
## accuracy against the OBS weights, and the accuracy of the plug-in.
check_schemes <- function(M = 400L, seed = BASESEED + 21L) {
  set.seed(seed)
  Ts <- sort(sample_T(M)); w <- mc_weights_sorted(Ts)
  f  <- function(t) cos(3 * t) + t^2
  tt <- seq(0, 1, length.out = 200001); dt <- diff(tt[1:2])
  fg <- f(tt) * G_DENSITY(tt)
  tr <- sum(fg[-1] + fg[-length(fg)]) * dt / 2
  cat(sprintf("  [MC] sum of weights - 1 = %.2e   quadrature error: MC %.2e, OBS %.2e\n",
              sum(w) - 1, sum(w * f(Ts)) - tr, mean(f(Ts)) - tr))
  stopifnot(abs(sum(w) - 1) < 1e-10)
  gh <- ghat_at(Ts, tt, b_of_M(M))
  cat(sprintf("  [gh] b = %.4f (MISE-optimal order M^(-1/5) = %.4f)   int ghat = %.4f   sup|ghat-g| = %.3f\n",
              b_of_M(M), M^(-0.2), sum(gh[-1] + gh[-length(gh)]) * dt / 2,
              max(abs(gh - G_DENSITY(tt)))))
  invisible(TRUE)
}


# =========================================================================
#  THE TWO EXTERNAL BENCHMARKS OF [1]
#
#  Each is evaluated by ONE function, used twice with the same code path:
#    - with t_eval = tg and target = mu on the evaluation grid, it returns
#      the ISE of the tunings passed to it ("mean"), used once the
#      selection has been made;
#    - with t_eval = the held-out visit times and target = the held-out
#      profiles, it returns the held-out sum of squares of every candidate
#      tuning ("sum"), which is the fold's contribution to the criterion.
#  The criterion is therefore literally the same object as the one used for
#  the spacing estimator, up to the constant 1/(N_u M).
# =========================================================================

## ll_crit(): the pooled local-linear smoother of CM2012 in the visit-time
## direction, fitted on (Tp, Y) and evaluated at every (u_l, t_eval). No
## smoothing in u: each profile is observed over the whole of U, so the fit
## is a section of the pooled scatterplot at each u, and the normal
## equations are those of one univariate local-linear fit shared by all the
## slices.
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
#  STORAGE  --  one .rds per experiment, written ATOMICALLY, every figure
#               and every table being built afterwards from that file alone
# =========================================================================

rds_file <- function(name) file.path(OUTDIR, paste0(name, ".rds"))

## save_rds(): a write that a power cut cannot corrupt. The object goes to
## <name>.rds.tmp; only once that file is complete is the current
## <name>.rds renamed <name>.rds.bak and the temporary file moved in its
## place. At every instant, at least one of the two files is a complete
## .rds, and load_rds() falls back on the backup if the main file is
## unreadable. Both moves are renames, so the cost does not grow with the
## number of checkpoints.
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
    stop(sprintf("%s not found or unreadable: run the corresponding experiment first.", f))
  obj
}

has_rds <- function(name)
  file.exists(rds_file(name)) || file.exists(paste0(rds_file(name), ".bak"))

## clean_ckpt(): remove the two auxiliary files of the atomic write once
## they are of no further use. <name>.rds.tmp is only ever present as the
## leftover of a write that was cut, and goes unconditionally. The backup
## goes only after <name>.rds has been READ BACK in full: nothing is
## dropped on the strength of a write that may not have landed, and a run
## that ends short of its replications keeps its backup.
clean_ckpt <- function(name) {
  f <- rds_file(name); tmp <- paste0(f, ".tmp"); bak <- paste0(f, ".bak")
  if (file.exists(tmp)) unlink(tmp)
  ok <- file.exists(f) && !is.null(tryCatch(readRDS(f), error = function(e) NULL))
  if (!ok) {
    message(sprintf("  [!] %s.rds is missing or unreadable: %s.rds.bak is kept",
                    name, name))
    return(invisible(FALSE))
  }
  if (file.exists(bak)) unlink(bak)
  cat(sprintf("  %s.rds verified; .tmp and .bak removed\n", name))
  invisible(TRUE)
}


# =========================================================================
#  OUTPUT HELPERS  --  one standalone TikZ document per panel, and
#                      booktabs tables
# =========================================================================

## Every figure file is a compilable document on its own, opened and closed
## by these two blocks, and is compiled by emit_fig() right after being
## written. The paper includes the resulting PDF, two panels on one line
## being obtained with, e.g.,
##   \includegraphics[width=0.48\linewidth]{fig_basis_mu1} ...
TEX_HEAD <- c("\\documentclass[border=2pt]{standalone}",
              "\\usepackage[T1]{fontenc}",
              "\\usepackage{tikz}",
              "\\begin{document}")
TEX_FOOT <- "\\end{document}"

## check_graphics(): the figure toolchain is needed by the make_* routines
## only, so that the Monte Carlo can be run on a machine without LaTeX.
GRAPHICS_OK <- FALSE
check_graphics <- function() {
  if (GRAPHICS_OK) return(invisible(TRUE))
  if (!requireNamespace("tikzDevice", quietly = TRUE))
    stop("package tikzDevice is required to write the figures")
  if (Sys.which(LATEX) == "")
    stop(sprintf("'%s' is not on the PATH: it compiles each figure to PDF", LATEX))
  options(tikzDefaultEngine = "pdftex",
          tikzMetricsDictionary = file.path(OUTDIR, "tikzMetrics"))
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
## The compilation is done from inside OUTDIR, so that a path with spaces
## or accents never reaches the command line. On failure the tail of the
## log is shown and the log is kept.
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
## in the standalone preamble above, then compile it into <name>.pdf. Both
## sizes are printed, since it is what makes these figures usable or not.
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

## write_table(): one booktabs table float per file, to be \input as before.
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

## box_panel(): the one figure both experiments draw, one group of boxes
## per competitor and one box per sample size, the reference being the zero
## line. D is a list, one entry per competitor, each a list of one vector
## per sample size. NOTE: as.vector() on a matrix of mode list does NOT
## drop the dim attribute and boxplot() then fails, hence the unlist().
box_panel <- function(name, D, labels, ylab, n_grid, width = 6.6) {
  ne <- length(D); nn <- length(n_grid)
  shade <- grey.colors(nn, start = 0.88, end = 0.55)
  emit_fig(name, width, 4.2, function() {
    yl <- range(sapply(seq_len(ne), function(e)
      sapply(seq_len(nn), function(jn) quantile(D[[e]][[jn]], c(.02, .98)))))
    yl <- yl + c(-0.05, 0.50) * max(diff(yl), 0.1)     # room for the legend
    at  <- as.vector(sapply(seq_len(ne), function(e) (e - 1) * (nn + 1) + seq_len(nn)))
    ctr <- sapply(seq_len(ne), function(e) mean((e - 1) * (nn + 1) + seq_len(nn)))
    par(mar = c(3.4, 4.4, 1.2, 0.6) * c(CEX, CEX, 1, 1),
        mgp = c(3.0, 0.7, 0) * CEX, tcl = -0.3,
        cex.axis = 0.85 * CEX, cex.lab = CEX)
    boxplot(unlist(D, recursive = FALSE), at = at, col = rep(shade, ne),
            boxwex = 0.7, outline = FALSE, xaxt = "n", las = 1, bty = "l",
            ylim = yl, ylab = ylab)
    abline(h = 0, lty = 2, col = "grey35")
    mtext(labels, side = 1, at = ctr, line = 1.0 * CEX, cex = 0.8 * CEX)
    legend("top", sprintf("$n=%d$", n_grid), fill = shade, horiz = TRUE,
           bty = "n", cex = 0.82 * CEX)
  })
}

EST_TEX <- c(ll   = "$\\widehat\\mu^{\\mathrm{LL}}$",
             spl  = "$\\widehat\\mu^{\\mathrm{spl}}$",
             rate = "$\\widehat\\mu^{\\mathrm{rate}}$")

EST_SCH <- c(sp0  = "$\\widehat\\mu^{\\rm sp}_{\\mathcal B_{0,K}}$",
             sp1  = "$\\widehat\\mu^{\\rm sp}_{\\mathcal B_{1,K}}$",
             gh   = "$\\widehat\\mu^{\\rm OBS}_{\\widehat g}$",
             obs  = "$\\widehat\\mu^{\\rm OBS}$",
             mc   = "$\\widehat\\mu^{\\rm MC}$")


# =========================================================================
#  [1] THE COMPETITORS  --  three fully data-driven estimators, one
#      criterion                             -> res_est.rds
#
#  On each replication, one partition of the subjects into V folds serves
#  the three families at once: the pair (basis, truncation) of the spacing
#  estimator over the 3 x Kmax candidates, the bandwidth h_t of the
#  local-linear smoother over LL_HGRID, and the penalty lambda of the
#  P-spline over SPL_LGRID. Each estimator is then refitted on the whole
#  sample AT ITS OWN SELECTED TUNING, and only that fit is scored, through
#  the paired log-ratios
#        Delta_e = log( ISE(mu^e) / ISE(mu^sp) ),  e in {LL, spl, rate}.
#  A positive Delta_e means that the spacing estimator was the more
#  accurate on that replication.
# =========================================================================

## est_sig(): the parameters two runs must share for their replications
## to be comparable, hence for a stored file to be continued. N_GRID and
## R_REP are NOT in it: sample sizes and replications can be added.
est_sig <- function(V) list(
  means = names(MEANS), V = V, Kmax = Kmax,
  k_convention = "K excludes the r monomials; dimension is r + K",
  a = a, b = b, cc = cc, P = P, lambda_m = lambda_m, m_cap = m_cap, tau = tau,
  L_modes = L_MODES, L_noise = L_NOISE,
  h_grid_t = LL_HGRID, spl_nik = SPL_NIK, lam_grid = SPL_LGRID,
  C_rate = C_RATE, r_basis = R_BASIS,
  nu_eval = NU_EVAL, nt_eval = NT_EVAL, baseseed = BASESEED)

## empty_slot(): the storage of one sample size, all replications pending.
empty_slot <- function(n, nm, R) list(
  n = as.integer(n), K_n = K_rate(n), R_done = 0L, secs = 0,
  M_rep = rep(NA_real_, R),
  ise = list(sp   = matrix(NA_real_, nm, R), ll   = matrix(NA_real_, nm, R),
             spl  = matrix(NA_real_, nm, R), rate = matrix(NA_real_, nm, R)),
  K_cv  = matrix(NA_integer_, nm, R), basis = matrix(NA_integer_, nm, R),
  h_ll  = matrix(NA_integer_, nm, R), lambda_spl = matrix(NA_integer_, nm, R))

## resize_slot(): the same slot with room for R replications, the completed
## ones being carried over. R may grow (the new columns are pending) or
## shrink (the last replications are dropped).
resize_slot <- function(z, nm, R) {
  keep <- seq_len(min(z$R_done, R))
  w <- empty_slot(z$n, nm, R)
  if (length(keep)) {
    w$M_rep[keep] <- z$M_rep[keep]
    for (e in names(w$ise)) w$ise[[e]][, keep] <- z$ise[[e]][, keep]
    w$K_cv[, keep]  <- z$K_cv[, keep]
    w$basis[, keep] <- z$basis[, keep]
    w$h_ll[, keep]  <- z$h_ll[, keep]
    w$lambda_spl[, keep] <- z$lambda_spl[, keep]
  }
  w$R_done <- length(keep); w$secs <- z$secs
  w
}

## cvbox_status(): what the stored file holds, without running anything.
cvbox_status <- function() {
  if (!has_rds(RDS_COMP)) {
    cat(sprintf("%s.rds: absent\n", RDS_COMP)); return(invisible(NULL))
  }
  S <- load_rds(RDS_COMP)
  cat(sprintf("%s.rds: R = %d asked for\n", RDS_COMP, S$R))
  for (z in S$by_n)
    cat(sprintf("   n = %5d : %4d / %4d replications done (%.1f min)\n",
                z$n, z$R_done, S$R, z$secs / 60))
  invisible(S)
}

## edge_report(): how often each tuning was selected at an end of its grid,
## per sample size, read off the stored file alone. make_cvbox() prints the
## same counts as percentages; this one breaks them down by n, which is
## what tells whether a grid is too narrow at every size or only at the
## largest. A count on the FIRST value of a grid calls for lowering that
## value, one on the LAST for raising it.
edge_report <- function() {
  S <- load_rds(RDS_COMP)
  nht <- length(S$h_grid_t); nlam <- length(S$lam_grid)
  tab <- sapply(S$by_n, function(z) {
    idx <- seq_len(z$R_done)
    if (!length(idx))
      return(c(ht_lo = 0, ht_hi = 0, lam_lo = 0, lam_hi = 0, K_hi = 0, sel = 0))
    c(ht_lo  = sum(z$h_ll[, idx] == 1L),
      ht_hi  = sum(z$h_ll[, idx] == nht),
      lam_lo = sum(z$lambda_spl[, idx] == 1L),
      lam_hi = sum(z$lambda_spl[, idx] == nlam),
      K_hi   = sum(z$K_cv[, idx] == S$Kmax),
      sel    = length(S$means) * length(idx))
  })
  colnames(tab) <- sprintf("n=%d", vapply(S$by_n, function(z) z$n, integer(1)))
  tab <- cbind(tab, total = rowSums(tab))
  cat(sprintf("grids: h_t [%.3f, %.3f]   lambda [%.1e, %.1e]   Kmax %d\n",
              S$h_grid_t[1], S$h_grid_t[nht], S$lam_grid[1], S$lam_grid[nlam],
              S$Kmax))
  cat("counts of selections at an end of a grid; 'sel' is the number of\n",
      "selections made per parameter (means x replications)\n", sep = "")
  print(tab)
  invisible(tab)
}

run_cvbox <- function(n_grid = N_GRID, R = R_REP, V = V_FOLDS,
                      resume = RESUME, every = CKPT_EVERY) {
  cat("\n[1] data-driven spacing, LL and spline, all tuned by the same CV ...\n")
  ug <- seq(0, 1, length.out = NU_EVAL); tg <- seq(0, 1, length.out = NT_EVAL)
  nm <- length(MEANS); nb <- 3L
  nht <- length(LL_HGRID); nlam <- length(SPL_LGRID)
  sig <- est_sig(V)
  
  ## ---- the stored file, if it is to be continued -------------------------
  ST <- NULL
  if (resume && has_rds(RDS_COMP)) {
    ST <- load_rds(RDS_COMP)
    if (!isTRUE(all.equal(ST$sig, sig)))
      stop(sprintf("%s.rds was produced under other parameters: move it aside, ",
                   RDS_COMP),
           "or call run_cvbox(resume = FALSE) to overwrite it.")
    old <- ST$by_n
    names(old) <- vapply(old, function(z) as.character(z$n), "")
    ## sample sizes already stored are kept, even if they left the grid
    all_n <- sort(unique(c(as.integer(n_grid), as.integer(names(old)))))
    ST$by_n <- lapply(all_n, function(nn) {
      k <- as.character(nn)
      if (!is.null(old[[k]])) resize_slot(old[[k]], nm, R) else empty_slot(nn, nm, R)
    })
    ST$n_grid <- all_n; ST$R <- R
    cat(sprintf("    resuming: %s\n", paste(vapply(ST$by_n, function(z)
      sprintf("n=%d %d/%d", z$n, z$R_done, R), ""), collapse = ", ")))
  }
  if (is.null(ST)) {
    all_n <- sort(unique(as.integer(n_grid)))
    ST <- list(means = names(MEANS), mean_tex = vapply(MEANS, `[[`, "", "tex"),
               n_grid = all_n, R = R, V = V, keys = c("ll", "spl", "rate"),
               Kmax = Kmax, C_rate = C_RATE, lambda_m = lambda_m,
               r_basis = R_BASIS, h_grid_t = LL_HGRID, lam_grid = SPL_LGRID,
               sig = sig, by_n = lapply(all_n, empty_slot, nm = nm, R = R))
  }
  
  ## the bases are built long enough for BOTH the CV grid, of largest
  ## dimension 2 + Kmax, and the rate rule, of dimension R_BASIS + K_n
  Dbig  <- 2L + max(Kmax, max(vapply(ST$n_grid, K_rate, integer(1))))
  co    <- lapply(0:2, basis_r, dim = Dbig)        # the three candidate bases
  Phi_t <- lapply(co, basis_phi, t = tg)
  MUev  <- lapply(MEANS, function(z) z$f(ug, tg))
  
  todo <- which(ST$n_grid %in% as.integer(n_grid))
  for (jn in todo) {
    n <- ST$n_grid[jn]; Kn <- K_rate(n); z <- ST$by_n[[jn]]
    if (z$R_done >= R) { cat(sprintf("    n=%5d  already complete\n", n)); next }
    
    for (rep in seq.int(z$R_done + 1L, R)) {
      t_rep <- proc.time()[3]
      ## the seed of THIS replication: the run is reproducible whatever the
      ## point at which it was interrupted
      set.seed(rep_seed(n, rep))
      
      ## one pool of design, subject fields and noise, shared by the 3 means
      pool <- gen_data(n, ug)
      Tv <- pool$T; z$M_rep[rep] <- length(Tv); h <- h_of_M(length(Tv))
      fold <- sample(rep_len(seq_len(V), n))       # folds are SUBJECT-level
      trn  <- lapply(seq_len(V), function(v) fold[pool$subj] != v)
      
      for (q in seq_len(nm)) {
        Y <- pool$Y + MEANS[[q]]$f(ug, Tv)         # same noise, another mean
        
        ## ---- the common subject-level V-fold criterion -------------------
        cv_sp  <- matrix(0, nb, Kmax)
        cv_ll  <- numeric(nht)
        cv_spl <- numeric(nlam)
        for (v in seq_len(V)) {
          tr  <- trn[[v]]
          Ttr <- Tv[tr];  Ytr <- Y[, tr,  drop = FALSE]
          Tte <- Tv[!tr]; Yte <- Y[, !tr, drop = FALSE]
          htr <- h_of_M(length(Ttr))
          for (bq in seq_len(nb)) {
            r0 <- bq - 1L                          # B_{r0,K} has dimension r0+K
            cv_sp[bq, ] <- cv_sp[bq, ] +
              sse_curve(beta_spacing(Ytr, Ttr, co[[bq]], htr),
                        basis_phi(co[[bq]], Tte), Yte, r0 + Kmax)[r0 + seq_len(Kmax)]
          }
          cv_ll  <- cv_ll  + ll_crit(Ytr, Ttr, ug, Tte, Yte, stat = "sum")
          cv_spl <- cv_spl + spl_crit(Ytr, Ttr, Tte, Yte, stat = "sum")
        }
        ij <- which(cv_sp == min(cv_sp, na.rm = TRUE), arr.ind = TRUE)[1, ]
        bs <- ij[1]; Ks <- ij[2]
        il <- which.min(cv_ll)
        is <- which.min(cv_spl)
        z$K_cv[q, rep] <- Ks; z$basis[q, rep] <- bs
        z$h_ll[q, rep] <- il; z$lambda_spl[q, rep] <- is
        
        ## ---- refit on the whole sample, at the selected tunings ----------
        ## the selected system is B_{bs-1, Ks}, of dimension bs-1+Ks
        ds <- (bs - 1L) + Ks
        Bfull <- beta_spacing(Y, Tv, co[[bs]], h)
        z$ise$sp[q, rep]  <- ise_curve(Bfull, Phi_t[[bs]], MUev[[q]], ds)[ds]
        z$ise$ll[q, rep]  <- ll_crit(Y, Tv, ug, tg, MUev[[q]],
                                     hgrid = LL_HGRID[il], stat = "mean")
        z$ise$spl[q, rep] <- spl_crit(Y, Tv, tg, MUev[[q]],
                                      lgrid = SPL_LGRID[is], stat = "mean")
        ## the rate-tuned competitor: same fit, system B_{R_BASIS, K_n}
        Brt <- if (bs == R_BASIS + 1L) Bfull else
          beta_spacing(Y, Tv, co[[R_BASIS + 1L]], h)
        dn <- R_BASIS + Kn
        z$ise$rate[q, rep] <- ise_curve(Brt, Phi_t[[R_BASIS + 1L]], MUev[[q]], dn)[dn]
      }
      
      ## ---- the replication is complete: it can now be checkpointed -------
      z$R_done <- rep; z$secs <- z$secs + (proc.time()[3] - t_rep)
      ST$by_n[[jn]] <- z
      ckpt <- (rep %% every == 0L) || (rep == R)
      if (ckpt) save_rds(RDS_COMP, ST, quiet = TRUE)
      if (rep %% 5 == 0 || rep == 1 || rep == R)
        cat(sprintf("    n=%5d  rep %4d / %4d   (%.1f s/rep%s)\n", n, rep, R,
                    z$secs / z$R_done, if (ckpt) ", saved" else ""))
    }
    
    ST$by_n[[jn]] <- z
    save_rds(RDS_COMP, ST, quiet = TRUE)
    idx <- seq_len(z$R_done)
    for (q in seq_len(nm))
      cat(sprintf(paste0("     %-3s  med log(ISE_LL/ISE_sp) %+.2f   ",
                         "spl %+.2f   rate %+.2f\n"), names(MEANS)[q],
                  median(log(z$ise$ll[q, idx]   / z$ise$sp[q, idx])),
                  median(log(z$ise$spl[q, idx]  / z$ise$sp[q, idx])),
                  median(log(z$ise$rate[q, idx] / z$ise$sp[q, idx]))))
  }
  
  save_rds(RDS_COMP, ST)
  ## An interrupted run never reaches this line, so its backup survives on
  ## disk, which is the whole point of keeping one. The test is a guard for
  ## the remaining case, a sample size left short of its R replications.
  if (all(vapply(ST$by_n[todo], function(z) z$R_done >= R, logical(1))))
    clean_ckpt(RDS_COMP)
  else
    cat(sprintf("  run still partial: %s.rds.bak kept\n", RDS_COMP))
  invisible(ST)
}


# =========================================================================
#  FIGURES AND TABLE OF [1], from res_est.rds
#      fig_cvbox_mu1, fig_cvbox_mu2, fig_cvbox_mu3, tab_cvbox
#
#  tab_cvbox collects, in one place, every tuning the criterion selected,
#  for the three data-driven estimators, next to the truncation K_n of the
#  rate-tuned competitor. Only the replications actually completed are
#  read, so the figures and the table can be built while the study is still
#  running, or after an interruption.
# =========================================================================

make_cvbox <- function() {
  cat(sprintf("\n[1] figures and table from %s.rds ...\n", RDS_COMP))
  S <- load_rds(RDS_COMP)
  slots <- Filter(function(z) z$R_done > 0L, S$by_n)
  if (!length(slots)) stop(sprintf("%s.rds holds no completed replication yet.",
                                   RDS_COMP))
  n_grid <- vapply(slots, function(z) as.integer(z$n), integer(1))
  Rd     <- vapply(slots, function(z) as.integer(z$R_done), integer(1))
  if (any(Rd < S$R))
    cat(sprintf("  [i] partial file: %s replications per sample size (asked: %d)\n",
                paste(sprintf("n=%d:%d", n_grid, Rd), collapse = ", "), S$R))
  nm <- length(S$means); nb <- 3L; nn <- length(slots); keys <- S$keys
  nht <- length(S$h_grid_t); nlam <- length(S$lam_grid)
  
  D    <- array(list(), c(nm, length(keys), nn))   # paired log-ratios
  freq <- array(0, c(nm, nn, nb))                  # basis selection frequencies
  Bmod <- matrix(NA_integer_, nm, nn)              # modal basis
  Bpct <- matrix(NA_real_, nm, nn)                 # its selection frequency
  Kmed <- matrix(NA_real_, nm, nn)                 # median K_CV IN that basis
  Hmed <- matrix(NA_real_, nm, nn)                 # median h_t selected
  Lmed <- matrix(NA_real_, nm, nn)                 # median lambda selected
  Krat <- vapply(slots, function(z) as.integer(z$K_n), integer(1))
  Mbar <- vapply(slots, function(z) mean(z$M_rep[seq_len(z$R_done)]), numeric(1))
  ## Selections at an end of a tuning grid, counted APART for each parameter
  ## and each end: a selection at the top of LL_HGRID and one at the bottom
  ## of SPL_LGRID call for opposite corrections, and "widen" is meaningless
  ## until the two are separated.
  edge <- c(ht_lo = 0L, ht_hi = 0L, lam_lo = 0L, lam_hi = 0L, K_hi = 0L)
  nsel <- 0L                                       # selections per parameter
  
  for (jn in seq_len(nn)) {
    z <- slots[[jn]]; idx <- seq_len(z$R_done)
    for (q in seq_len(nm)) {
      for (e in seq_along(keys))
        D[[q, e, jn]] <- log(z$ise[[keys[e]]][q, idx] / z$ise$sp[q, idx])
      freq[q, jn, ] <- 100 * tabulate(z$basis[q, idx], nb) / length(idx)
      ## the selected PAIR: modal basis, and the median truncation over the
      ## replications that selected THAT basis -- a median of K_CV pooled
      ## across bases would compare truncations of different systems.
      Bmod[q, jn] <- which.max(freq[q, jn, ]) - 1L
      Bpct[q, jn] <- max(freq[q, jn, ])
      Kmed[q, jn] <- median(z$K_cv[q, idx][z$basis[q, idx] == Bmod[q, jn] + 1L])
      Hmed[q, jn] <- median(S$h_grid_t[z$h_ll[q, idx]])
      Lmed[q, jn] <- median(S$lam_grid[z$lambda_spl[q, idx]])
      ht <- z$h_ll[q, idx]; la <- z$lambda_spl[q, idx]; Kc <- z$K_cv[q, idx]
      edge <- edge + c(sum(ht == 1L), sum(ht == nht),
                       sum(la == 1L), sum(la == nlam), sum(Kc == S$Kmax))
      nsel <- nsel + length(idx)
    }
  }
  if (any(edge > 0L)) {
    txt <- c(
      ht_lo  = sprintf("h_t at %.3f, FIRST of LL_HGRID     -> lower its first value",
                       S$h_grid_t[1]),
      ht_hi  = sprintf("h_t at %.3f, LAST of LL_HGRID      -> raise its last value",
                       S$h_grid_t[nht]),
      lam_lo = sprintf("lambda at %.1e, FIRST of SPL_LGRID -> lower its first value",
                       S$lam_grid[1]),
      lam_hi = sprintf("lambda at %.1e, LAST of SPL_LGRID  -> raise its last value",
                       S$lam_grid[nlam]),
      K_hi   = sprintf("K_CV at %d, LAST of 1..Kmax         -> raise Kmax",
                       S$Kmax))
    cat(sprintf("  [i] tunings selected at an end of their grid, out of %d selections per parameter:\n",
                nsel))
    for (k in names(edge)[edge > 0L])
      cat(sprintf("      %5d  (%5.1f%%)   %s\n", edge[[k]],
                  100 * edge[[k]] / nsel, txt[[k]]))
    warning(sprintf("selections at an end of a tuning grid: %s (out of %d per parameter)",
                    paste(sprintf("%s %d", names(edge)[edge > 0L], edge[edge > 0L]),
                          collapse = ", "), nsel))
  }
  
  ## ---- long-format summary ----------------------------------------------
  res <- NULL
  for (q in seq_len(nm)) for (e in seq_along(keys)) for (jn in seq_len(nn)) {
    v <- D[[q, e, jn]]
    res <- rbind(res, data.frame(
      mean = S$means[q], est = keys[e], n = n_grid[jn], R_done = Rd[jn],
      med = median(v), win = 100 * mean(v > 0),
      basis = Bmod[q, jn], basis_pct = Bpct[q, jn],
      K_cv = Kmed[q, jn],
      h_t = Hmed[q, jn], lambda = Lmed[q, jn], K_rate = Krat[jn]))
  }
  rownames(res) <- NULL
  
  ## ---- one figure per mean ----------------------------------------------
  ## Three competitor groups, one box per sample size. The tunings selected
  ## by the criterion are NOT printed on the figure: they are in tab_cvbox.
  for (q in seq_len(nm))
    box_panel(sprintf("fig_cvbox_%s", S$means[q]),
              lapply(seq_along(keys), function(e)
                lapply(seq_len(nn), function(jn) D[[q, e, jn]])),
              EST_TEX[keys], "$\\Delta_e$", n_grid, width = 8.6)
  
  ## ---- tab_cvbox: every tuning the criterion selected --------------------
  Rtex <- if (length(unique(Rd)) == 1L) sprintf("$R=%d$ replications", Rd[1]) else
    sprintf("between $%d$ and $%d$ replications", min(Rd), max(Rd))
  body <- character(0)
  for (q in seq_len(nm)) {
    if (q > 1L) body <- c(body, "\\midrule")
    for (jn in seq_len(nn)) {
      row <- paste(c(
        if (jn == 1L) sprintf("\\multirow{%d}{*}{%s}", nn, S$mean_tex[q]) else "",
        sprintf("$%d$", n_grid[jn]), sprintf("$%d$", as.integer(round(Mbar[jn]))),
        sprintf("$%d$ (%.0f\\%%)", Bmod[q, jn], Bpct[q, jn]),
        sprintf("$%g$", Kmed[q, jn]),
        num_tex(Hmed[q, jn]), num_tex(Lmed[q, jn]),
        sprintf("$%d$", Krat[jn])), collapse = " & ")
      body <- c(body, paste0(row, " \\\\"))
    }
  }
  write_table("tab_cvbox", "llrccccc",
              paste0("& & & \\multicolumn{2}{c}{$\\widehat\\mu^{\\rm sp}$} & ",
                     "$\\widehat\\mu^{\\rm LL}$ & $\\widehat\\mu^{\\rm spl}$ & ",
                     "$\\widehat\\mu^{\\rm rate}$ \\\\",
                     "\\cmidrule(lr){4-5}\\cmidrule(lr){6-6}\\cmidrule(lr){7-7}",
                     "\\cmidrule(lr){8-8}\n",
                     "mean & $n$ & $\\overline M$ & $\\widehat r$ & ",
                     "$K_{\\mathrm{CV}}$ & $\\widehat h_t$ & ",
                     "$\\widehat\\lambda$ & $K_n$ \\\\"),
              body,
              sprintf(paste0("Tunings selected by the criterion~\\eqref{eq:cv} ",
                             "over %s: for ",
                             "$\\widehat\\mu^{\\rm sp}$, the modal augmentation ",
                             "$\\widehat r$ with its ",
                             "selection frequency and the median truncation over ",
                             "the replications that selected that augmentation, ",
                             "the retained system being $\\mathcal{B}_{\\widehat r, ",
                             "K_{\\mathrm{CV}}}$, of dimension $\\widehat r + ",
                             "K_{\\mathrm{CV}}$; for the ",
                             "two smoothers, the median selected bandwidth and ",
                             "penalty. The last column is the truncation ",
                             "$K_n$ of the rate-tuned competitor, fixed by ",
                             "$r + K_n = \\lceil %g\\,(%d n)^{1/5} \\rceil$ with ",
                             "$r=%d$, which receives no selection."),
                      Rtex, S$C_rate, as.integer(S$lambda_m), S$r_basis),
              "tab:cvbox")
  
  print(res[, c("mean", "est", "n", "R_done", "med", "win", "basis", "basis_pct",
                "K_cv", "K_rate")], digits = 3)
  invisible(res)
}


# =========================================================================
#  [2] BASES AND WEIGHTING SCHEMES, on mu_1  -> res_schemes.rds
#
#  Each replication opens with the seed of [1] and rebuilds its fold
#  partition, so the design, the subject fields, the errors and the blocks
#  are the same: the boxes of [2] are paired with those of fig_cvbox_mu1.
#  Nothing is written but res_schemes.rds. The fully data-driven mu^sp,
#  the reference of the second panel, is not refitted here: [1] stored it
#  on these very draws, and make_schemes() reads its ISE and its selected
#  (r, K) from res_est.rds.
# =========================================================================

KEYS_W   <- c("gh", "obs", "mc")           # (w, g) schemes, in panel order
SCH_KEYS <- c("sp0", "sp1", "sp2", KEYS_W)

sch_sig <- function(V) list(
  V = V, Kmax = Kmax,
  k_convention = "K excludes the r monomials; dimension is r + K", a = a, b = b, cc = cc, P = P, lambda_m = lambda_m,
  m_cap = m_cap, tau = tau, L_modes = L_MODES, L_noise = L_NOISE,
  nu_eval = NU_EVAL, nt_eval = NT_EVAL, baseseed = BASESEED,
  gh_floor = GH_FLOOR, gh_kernel = "gaussian", b_rule = deparse(body(b_of_M)))

sch_slot <- function(n, R) list(
  n = as.integer(n), R_done = 0L, secs = 0,
  M_rep = rep(NA_real_, R), b_M = rep(NA_real_, R), gh_floor = rep(NA_real_, R),
  ise = setNames(replicate(length(SCH_KEYS), rep(NA_real_, R), simplify = FALSE),
                 SCH_KEYS),
  K_r = matrix(NA_integer_, 3L, R),               # K selected inside each B_r
  sel = setNames(replicate(length(KEYS_W), matrix(NA_integer_, 2L, R),
                           simplify = FALSE), KEYS_W))

sch_resize <- function(z, R) {
  keep <- seq_len(min(z$R_done, R)); w <- sch_slot(z$n, R)
  if (length(keep)) {
    w$M_rep[keep] <- z$M_rep[keep]; w$b_M[keep] <- z$b_M[keep]
    w$gh_floor[keep] <- z$gh_floor[keep]
    for (e in SCH_KEYS) w$ise[[e]][keep] <- z$ise[[e]][keep]
    w$K_r[, keep] <- z$K_r[, keep]
    for (e in names(w$sel)) w$sel[[e]][, keep] <- z$sel[[e]][, keep]
  }
  w$R_done <- length(keep); w$secs <- z$secs
  w
}

sch_status <- function() {
  if (!has_rds(RDS_SCH)) {
    cat(sprintf("%s.rds: absent\n", RDS_SCH)); return(invisible(NULL))
  }
  S <- load_rds(RDS_SCH)
  cat(sprintf("%s.rds: R = %d asked for\n", RDS_SCH, S$R))
  for (z in S$by_n)
    cat(sprintf("   n = %5d : %4d / %4d replications done (%.1f min)\n",
                z$n, z$R_done, S$R, z$secs / 60))
  invisible(S)
}

run_schemes <- function(n_grid = N_GRID, R = R_REP, V = V_FOLDS,
                        resume = RESUME, every = CKPT_EVERY) {
  cat("\n[2] augmentation and design density, on mu_1 ...\n")
  ug <- seq(0, 1, length.out = NU_EVAL); tg <- seq(0, 1, length.out = NT_EVAL)
  nb <- 3L; sig <- sch_sig(V)
  co    <- lapply(0:2, basis_r, dim = 2L + Kmax)   # largest dimension needed
  Phi_t <- lapply(co, basis_phi, t = tg)
  MUev  <- MU_1(ug, tg)
  
  ST <- NULL
  if (resume && has_rds(RDS_SCH)) {
    ST <- load_rds(RDS_SCH)
    if (!isTRUE(all.equal(ST$sig, sig)))
      stop(sprintf("%s.rds was produced under other parameters: move it aside, ",
                   RDS_SCH), "or call run_schemes(resume = FALSE) to overwrite it.")
    old <- ST$by_n; names(old) <- vapply(old, function(z) as.character(z$n), "")
    all_n <- sort(unique(c(as.integer(n_grid), as.integer(names(old)))))
    ST$by_n <- lapply(all_n, function(nn) {
      k <- as.character(nn)
      if (!is.null(old[[k]])) sch_resize(old[[k]], R) else sch_slot(nn, R) })
    ST$n_grid <- all_n; ST$R <- R
    cat(sprintf("    resuming: %s\n", paste(vapply(ST$by_n, function(z)
      sprintf("n=%d %d/%d", z$n, z$R_done, R), ""), collapse = ", ")))
  }
  if (is.null(ST)) {
    all_n <- sort(unique(as.integer(n_grid)))
    ST <- list(n_grid = all_n, R = R, V = V, Kmax = Kmax, keys = SCH_KEYS,
               keys_w = KEYS_W, sig = sig, by_n = lapply(all_n, sch_slot, R = R))
  }
  
  todo <- which(ST$n_grid %in% as.integer(n_grid))
  for (jn in todo) {
    n <- ST$n_grid[jn]; z <- ST$by_n[[jn]]
    if (z$R_done >= R) { cat(sprintf("    n=%5d  already complete\n", n)); next }
    
    for (rep in seq.int(z$R_done + 1L, R)) {
      t_rep <- proc.time()[3]
      ## the seed, the pool and the folds of [1], drawn in the same order
      set.seed(rep_seed(n, rep))
      pool <- gen_data(n, ug)
      Tv <- pool$T; M <- length(Tv); h <- h_of_M(M)
      fold <- sample(rep_len(seq_len(V), n))
      trn  <- lapply(seq_len(V), function(v) fold[pool$subj] != v)
      Y <- pool$Y + MU_1(ug, Tv)                   # the mean mu_1 only
      z$M_rep[rep] <- M
      
      ## ---- the common subject-level V-fold criterion --------------------
      cv <- array(0, c(1L + length(KEYS_W), nb, Kmax),
                  dimnames = list(c("sp", KEYS_W), NULL, NULL))
      for (v in seq_len(V)) {
        tr  <- trn[[v]]
        Ttr <- Tv[tr];  Ytr <- Y[, tr,  drop = FALSE]
        Tte <- Tv[!tr]; Yte <- Y[, !tr, drop = FALSE]
        htr <- h_of_M(length(Ttr)); wts <- weight_set(Ttr)
        for (bq in seq_len(nb)) {
          r0 <- bq - 1L                            # B_{r0,K} has dimension r0+K
          Pte <- basis_phi(co[[bq]], Tte)
          keep <- r0 + seq_len(Kmax)
          cv["sp", bq, ] <- cv["sp", bq, ] +
            sse_curve(beta_spacing(Ytr, Ttr, co[[bq]], htr),
                      Pte, Yte, r0 + Kmax)[keep]
          for (e in KEYS_W)
            cv[e, bq, ] <- cv[e, bq, ] +
            sse_curve(beta_weighted(Ytr, Ttr, co[[bq]], wts[[e]]$w, wts[[e]]$g),
                      Pte, Yte, r0 + Kmax)[keep]
        }
      }
      pick <- function(A) as.integer(which(A == min(A, na.rm = TRUE),
                                           arr.ind = TRUE)[1, ])
      
      ## ---- refit on the whole sample, at the selected tunings -----------
      ## The fully data-driven mu^sp, which selects (r, K) jointly, is NOT
      ## recomputed here: [1] already stored it, on these very draws, and
      ## make_schemes() reads it from res_est.rds. Only the per-system
      ## truncations are read off this CV surface.
      wfull <- weight_set(Tv)
      z$b_M[rep] <- wfull$b; z$gh_floor[rep] <- wfull$floor
      
      Kr <- apply(cv["sp", , ], 1, which.min)      # K inside each system
      z$K_r[, rep] <- as.integer(Kr)
      for (r0 in 0:2) {
        Br <- beta_spacing(Y, Tv, co[[r0 + 1L]], h)
        dr <- r0 + Kr[r0 + 1L]
        z$ise[[paste0("sp", r0)]][rep] <-
          ise_curve(Br, Phi_t[[r0 + 1L]], MUev, dr)[dr]
      }
      for (e in KEYS_W) {
        s <- pick(cv[e, , ]); z$sel[[e]][, rep] <- s
        Be <- beta_weighted(Y, Tv, co[[s[1]]], wfull[[e]]$w, wfull[[e]]$g)
        de <- (s[1] - 1L) + s[2]
        z$ise[[e]][rep] <- ise_curve(Be, Phi_t[[s[1]]], MUev, de)[de]
      }
      
      z$R_done <- rep; z$secs <- z$secs + (proc.time()[3] - t_rep)
      ST$by_n[[jn]] <- z
      ckpt <- (rep %% every == 0L) || (rep == R)
      if (ckpt) save_rds(RDS_SCH, ST, quiet = TRUE)
      if (rep %% 25 == 0 || rep == 1 || rep == R)
        cat(sprintf("    n=%5d  rep %4d / %4d   (%.2f s/rep%s)\n", n, rep, R,
                    z$secs / z$R_done, if (ckpt) ", saved" else ""))
    }
    
    ST$by_n[[jn]] <- z; save_rds(RDS_SCH, ST, quiet = TRUE)
    idx <- seq_len(z$R_done)
    cat(sprintf("     med log(ISE_e/ISE_B2): r=0 %+.2f  r=1 %+.2f  ghat %+.2f  OBS %+.2f  MC %+.2f\n",
                median(log(z$ise$sp0[idx] / z$ise$sp2[idx])),
                median(log(z$ise$sp1[idx] / z$ise$sp2[idx])),
                median(log(z$ise$gh[idx]  / z$ise$sp2[idx])),
                median(log(z$ise$obs[idx] / z$ise$sp2[idx])),
                median(log(z$ise$mc[idx]  / z$ise$sp2[idx]))))
  }
  
  save_rds(RDS_SCH, ST)
  if (all(vapply(ST$by_n[todo], function(z) z$R_done >= R, logical(1))))
    clean_ckpt(RDS_SCH)
  else
    cat(sprintf("  run still partial: %s.rds.bak kept\n", RDS_SCH))
  invisible(ST)
}


# =========================================================================
#  FIGURES AND TABLES OF [2]
#      fig_basis_mu1 / tab_basis : the augmentation, against B_{2,K}
#      fig_dens_mu1  / tab_dens  : the ratio-based schemes, against mu^sp
#
#  The reference of the second pair is the fully data-driven mu^sp, which
#  [1] already fitted on these very draws: its ISE and its selected (r, K)
#  are READ from res_est.rds rather than recomputed, and only the
#  replications completed in BOTH files are used.
# =========================================================================

make_schemes <- function() {
  cat(sprintf("\n[2] figures and tables from %s.rds and %s.rds ...\n",
              RDS_SCH, RDS_COMP))
  S <- load_rds(RDS_SCH)
  if (!has_rds(RDS_COMP))
    stop(sprintf("%s.rds is needed for the reference mu^sp: run run_cvbox() first.",
                 RDS_COMP))
  C <- load_rds(RDS_COMP)
  cn <- vapply(C$by_n, function(w) as.integer(w$n), integer(1))
  
  ## the sample sizes present in both files, and the replications common
  ## to the two: the pairing is exact, both having run on rep_seed(n, rep)
  slots <- Filter(function(z) z$R_done > 0L && z$n %in% cn, S$by_n)
  if (!length(slots))
    stop("no sample size with completed replications in both files yet.")
  cslots <- lapply(slots, function(z) C$by_n[[which(cn == z$n)]])
  n_grid <- vapply(slots, function(z) as.integer(z$n), integer(1))
  Rd <- as.integer(mapply(function(z, w) min(z$R_done, w$R_done), slots, cslots))
  if (any(Rd < S$R))
    cat(sprintf("  [i] %s replications common to both files (asked: %d)\n",
                paste(sprintf("n=%d:%d", n_grid, Rd), collapse = ", "), S$R))
  if (any(Rd == 0L)) stop("some sample size has no replication in common.")
  
  ise_sp <- lapply(seq_along(slots), function(j) cslots[[j]]$ise$sp[1, seq_len(Rd[j])])
  ise    <- function(e, j) slots[[j]]$ise[[e]][seq_len(Rd[j])]
  
  ## ---- the two panels ----------------------------------------------------
  DA <- lapply(c("sp0", "sp1"), function(e)
    lapply(seq_along(slots), function(j) log(ise(e, j) / ise("sp2", j))))
  DB <- lapply(S$keys_w, function(e)
    lapply(seq_along(slots), function(j) log(ise(e, j) / ise_sp[[j]])))
  
  box_panel("fig_basis_mu1", DA, EST_SCH[c("sp0", "sp1")],
            paste0("$\\log\\{\\mathrm{ISE}(\\cdot)/",
                   "\\mathrm{ISE}(\\widehat\\mu^{\\rm sp}_{\\mathcal B_{2,K}})\\}$"),
            n_grid, width = 5.4)
  box_panel("fig_dens_mu1", DB, EST_SCH[S$keys_w],
            paste0("$\\log\\{\\mathrm{ISE}(\\cdot)/",
                   "\\mathrm{ISE}(\\widehat\\mu^{\\rm sp})\\}$"),
            n_grid)
  
  ## ---- tab_basis: the truncation selected inside each system -------------
  ## Median and interquartile range of K_CV in B_{0,K}, B_{1,K}, B_{2,K}.
  qtex <- function(v) sprintf("$%g$ $[%g,%g]$", median(v),
                              unname(quantile(v, .25)), unname(quantile(v, .75)))
  body <- vapply(seq_along(slots), function(j) {
    Kr <- slots[[j]]$K_r[, seq_len(Rd[j]), drop = FALSE]
    paste0(paste(c(sprintf("$%d$", n_grid[j]),
                   vapply(1:3, function(i) qtex(Kr[i, ]), "")),
                 collapse = " & "),
           " \\\\")
  }, "")
  write_table("tab_basis", "lccc",
              paste0("$n$ & $\\mathcal{B}_{0,K}$ & $\\mathcal{B}_{1,K}$ & ",
                     "$\\mathcal{B}_{2,K}$ \\\\"),
              body,
              paste0("Truncation selected by the criterion~\\eqref{eq:cv} inside ",
                     "each candidate system, on $\\mu_1$: median and, in brackets, ",
                     "interquartile range of $K$ over the replications, the $r$ ",
                     "monomials of the augmentation excluded, so that the fitted ",
                     "space has dimension $r+K$."),
              "tab:basis")
  
  ## ---- tab_dens: what the criterion selected for the four estimators -----
  ## For each estimator, the modal augmentation with its frequency and the
  ## median truncation over the replications that selected it, as in
  ## Table~\ref{tab:cvbox}.
  rk_tex <- function(bas, Kcv) {
    fr <- tabulate(bas, 3L); m <- which.max(fr)
    c(sprintf("$%d$ (%.0f\\%%)", m - 1L, 100 * fr[m] / length(bas)),
      sprintf("$%g$", median(Kcv[bas == m])))
  }
  body <- vapply(seq_along(slots), function(j) {
    idx <- seq_len(Rd[j])
    cells <- c(rk_tex(cslots[[j]]$basis[1, idx], cslots[[j]]$K_cv[1, idx]),
               unlist(lapply(S$keys_w, function(e)
                 rk_tex(slots[[j]]$sel[[e]][1, idx], slots[[j]]$sel[[e]][2, idx]))))
    paste0(paste(c(sprintf("$%d$", n_grid[j]), cells), collapse = " & "), " \\\\")
  }, "")
  write_table("tab_dens", "lcccccccc",
              paste0("& \\multicolumn{2}{c}{$\\widehat\\mu^{\\rm sp}$} & ",
                     "\\multicolumn{2}{c}{$\\widehat\\mu^{\\rm OBS}_{\\widehat g}$} & ",
                     "\\multicolumn{2}{c}{$\\widehat\\mu^{\\rm OBS}$} & ",
                     "\\multicolumn{2}{c}{$\\widehat\\mu^{\\rm MC}$} \\\\",
                     "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
                     "\\cmidrule(lr){8-9}\n",
                     "$n$ & $\\widehat r$ & $K_{\\mathrm{CV}}$ & $\\widehat r$ & ",
                     "$K_{\\mathrm{CV}}$ & $\\widehat r$ & $K_{\\mathrm{CV}}$ & ",
                     "$\\widehat r$ & $K_{\\mathrm{CV}}$ \\\\"),
              body,
              paste0("Pairs $(r,K)$ selected by the criterion~\\eqref{eq:cv} on ",
                     "$\\mu_1$: the modal augmentation with its selection ",
                     "frequency, and the median truncation over the replications ",
                     "that selected it, the retained system being ",
                     "$\\mathcal{B}_{\\widehat r,K_{\\mathrm{CV}}}$, of dimension ",
                     "$\\widehat r + K_{\\mathrm{CV}}$."),
              "tab:dens")
  
  ## ---- console summary ---------------------------------------------------
  res <- NULL
  for (e in c("sp0", "sp1")) for (j in seq_along(slots)) {
    v <- log(ise(e, j) / ise("sp2", j))
    res <- rbind(res, data.frame(panel = "basis", est = e, n = n_grid[j],
                                 R_done = Rd[j], med = median(v),
                                 win = 100 * mean(v > 0)))
  }
  for (e in S$keys_w) for (j in seq_along(slots)) {
    v <- log(ise(e, j) / ise_sp[[j]])
    res <- rbind(res, data.frame(panel = "density", est = e, n = n_grid[j],
                                 R_done = Rd[j], med = median(v),
                                 win = 100 * mean(v > 0)))
  }
  rownames(res) <- NULL; print(res, digits = 3)
  for (j in seq_along(slots)) {
    z <- slots[[j]]; idx <- seq_len(Rd[j])
    cat(sprintf("  n=%5d : Mbar %5.0f   b = %.4f   ghat floored on %.2f%% of the design\n",
                z$n, mean(z$M_rep[idx]), mean(z$b_M[idx]), 100 * mean(z$gh_floor[idx])))
  }
  invisible(res)
}


# =========================================================================
#  DRIVER
#
#  The two halves of each experiment are independent: once a .rds is there,
#  its figures and tables can be rebuilt alone, by sourcing the script with
#  SIM_NO_RUN defined and calling make_cvbox() or make_schemes(). The
#  reports printed here, cvbox_status()/sch_status() before a run and
#  edge_report() after it, compute nothing: they read the .rds and can be
#  called on their own at any moment, including while a study is running in
#  another session. The panels that need no Monte Carlo are in
#  simul_viz.R.
# =========================================================================

if (!exists("SIM_NO_RUN")) {
  t0 <- proc.time()[3]
  cat("==== self-tests ====\n"); check_schemes()
  if (has_rds(RDS_COMP)) { cat("==== [1] checkpoint found ====\n"); cvbox_status() }
  if (has_rds(RDS_SCH))  { cat("==== [2] checkpoint found ====\n"); sch_status() }
  
  if (RUN_COMP) {
    run_cvbox(); E <- make_cvbox()
    cat("\n============== SUMMARY OF [1] ==============\n")
    cat("median Delta_e (positive = the data-driven spacing estimator wins):\n")
    print(reshape(E[, c("mean", "est", "n", "med")],
                  idvar = c("mean", "est"), timevar = "n", direction = "wide"),
          digits = 3)
    cat("tunings selected by the criterion, and the rate truncation K_n:\n")
    print(reshape(E[E$est == "ll",
                    c("mean", "n", "basis", "K_cv", "h_t", "lambda", "K_rate")],
                  idvar = "mean", timevar = "n", direction = "wide"), digits = 3)
    ## the same saturation counts as make_cvbox(), but broken down by sample
    ## size: it tells whether a grid is too narrow everywhere or only at the
    ## largest n, hence which end of which grid to move.
    cat("selections at an end of a tuning grid, per sample size:\n")
    edge_report()
  }
  
  if (RUN_SCH) {
    run_schemes(); A <- make_schemes()
    cat("\n============== SUMMARY OF [2] ==============\n")
    cat("median log-ratios (positive = the reference estimator wins):\n")
    print(reshape(A[, c("panel", "est", "n", "med")],
                  idvar = c("panel", "est"), timevar = "n", direction = "wide"),
          digits = 3)
  }
  
  cat(sprintf("\nTotal wall time: %.1f s\nOutputs written to: %s\n",
              proc.time()[3] - t0, normalizePath(OUTDIR)))
}