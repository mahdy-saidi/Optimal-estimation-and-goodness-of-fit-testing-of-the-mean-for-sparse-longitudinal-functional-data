###########################################################################
#            MULTIPLIER BOOTSTRAP TEST OF  H0 : mu = 0
#   
#   SELF-CONTAINED: the pieces of the data-generating mechanism
#   it needs are reproduced here, so that this file runs on its own.
#
#   TWO WEIGHT SCHEMES. The statistic is built from the same subject
#   blocks in both cases, only the weights carried by an observation
#   differ:
#
#     spacing (BENCHMARK, uses no design density)
#         w_ij,k = (1/2h) int_{T_(l-h)}^{T_(l+h)} phi_k ,  l = rank of T_ij,
#     g known (COMPETITOR)
#         w_ij,k = phi_k(T_ij) / (M g(T_ij))              (OBS weights),
#
#     v_ik = sum_j w_ij,k Y_ij  in L^2(U),
#     T_n   = sum_{k <= kbar} sum_{i != i'} <v_ik, v_i'k>,
#     T_n^* = sum_{k <= kbar} sum_{i != i'} xi_i xi_i' <v_ik, v_i'k>,
#
#   with ONE multiplier per SUBJECT, so that the within-subject dependence
#   carried by X_i is left intact. The diagonal blocks i = i' are dropped
#   because they carry the error and process variance. Both schemes use the
#   half-cosine system, so that the comparison isolates the quadrature.
#
#   NOTATION, as in the paper. The basis is B_{r,K}, obtained from the
#   dictionary (c_1, p_1, ..., p_r, c_2, ..., c_K): K counts the COSINES
#   only, the r monomials are excluded from it. The test needs no
#   augmentation, so this file works in B_0 = B_{0,K} throughout.
#   The truncation of the test is kbar, the K-bar of the paper:
#   the index set is K_M = {1, ..., kbar} and the statistic sums over the
#   first n_basis(kbar) functions of the system.
#
#   OPTIONAL CENTRING. With CENTER = TRUE every profile is replaced, before
#   the subject blocks are formed, by its deviation from the pooled mean
#   profile,
#       Y_ij(u) - Ybar(u),   Ybar(u) = (1/M) sum_{i,j} Y_ij(u),
#   so that T_n and all its bootstrap copies are built from centred data.
#   The map is linear, hence it is applied to the null part and to each
#   alternative direction separately and the amplitude grid is still
#   reached exactly. CENTER is part of the stored parameters, so a file
#   is only continued under the same choice.
#
#   THE EXPERIMENT. A power study: rejection frequencies against the
#   amplitude, for the two schemes, two alternative directions and three
#   truncations, stored in res_test.rds.
#
#   ORGANISATION. The experiment is run once and stores EVERYTHING it
#   computed in res_test.rds; the figures and the table are then built
#   from that file alone:
#       run_all() ->  res_test.rds  then  write_outputs()  ->  .tex/.pdf
#   What is stored is not the rejection indicator but the number of
#   bootstrap draws below the observed statistic, so that the nominal level
#   ALPHA is a parameter of the reporting layer only, and so is the
#   abscissa cut XCUT1 of the power figures.
#
#   CHECKPOINTS AND RESUMPTION. The .rds is rewritten every CKPT_EVERY
#   replications, not once at the end, so that a machine that stops loses
#   at most CKPT_EVERY replications. Three points make this safe and
#   exact:
#     (i)   the write is atomic: the object goes to <name>.rds.tmp, the
#           current file is renamed <name>.rds.bak, and the temporary
#           file then takes its place, so a cut always leaves one COMPLETE
#           file behind, never a half-written one;
#     (ii)  each replication opens with its OWN seed, seed(n, r), so the
#           replication that follows a restart is the very one an
#           uninterrupted run would have produced: the results do not
#           depend on where the run was cut, nor on the order of the
#           sample sizes;
#     (iii) each stored slot carries R_done, the number of replications
#           completed at that n, so the reporting functions run on a
#           partial file and the curves can be looked at while the study
#           is running.
#   Restarting the script with RESUME = TRUE picks the file up where it
#   stopped. R_REP may be raised between two runs (the replications
#   already stored are kept and the new ones are appended) and N_GRID may
#   be extended (a sample size already computed is not recomputed, and one
#   dropped from the grid is kept in the file). Once every replication
#   asked for is in, and only then, the file is read back in full and the
#   two auxiliary files .rds.tmp and .rds.bak are removed, so a completed
#   study leaves nothing but its .rds behind. A stored file produced under
#   other parameters (another mechanism, another bootstrap, other
#   truncations or amplitudes) is refused rather than continued; ALPHA and
#   XCUT1 are deliberately NOT among them, since they enter the reporting
#   only.
#
#   WHAT IT PRODUCES (and nothing else). One panel per file, each a
#   standalone LaTeX document COMPILED on the spot, so that every panel
#   comes as a .tex kept for editing and as the .pdf the paper includes:
#     fig_test_power_mu1          rejection frequency against the
#                                 amplitude for the SPACING statistic,
#                                 alternative mu_1, kbar = KBAR_MAX, one
#                                 curve per number of subjects, abscissa
#                                 cut at XCUT1;
#     fig_test_power_mu3_kbar<k>  the same for the alternative mu_3 with
#                                 a = b = 0, one file per kbar, full
#                                 amplitude grid;
#     fig_test_cmp_mu1_n<n>       the two weight schemes against each
#                                 other, one panel per number of subjects,
#                                 alternative mu_1 at kbar = KBAR_MAX,
#                                 abscissa cut at XCUT1;
#     tab_test_level.tex          empirical level (amplitude 0) of the
#                                 spacing statistic, one row per kbar, one
#                                 column per n.
#
###########################################################################


# =========================================================================
#  1. USER PARAMETERS  (everything meant to be edited sits here)
# =========================================================================

## ---- 1.1 the two alternative directions ---------------------------------
## mu_1(u,t) = (a t + b t^2) sum_p ((-1)^p/p^2) cos(pi p t + c pi u)
## mu_3(u,t) = a3 t + b3 t^2 + cos(c pi t) sum_p ((-1)^p/p^2) cos(pi p u)
a  <- 1;  b  <- -7;  cc <- 4;  P <- 8L   # shape parameters, as in the
# simulation study ("c" is cc here)
A3 <- 0;  B3 <- 0                       # the a and b of mu_3: zero, so that
# mu_3 carries the single mode c + 1

## ---- 1.2 design of the experiment ---------------------------------------
N_GRID   <- c(100L, 200L, 500L)   # numbers of subjects
R_REP    <- 2000L                 # replications per configuration
lambda_m <- 7L                    # Poisson mean of the visit counts
m_cap    <- 15L                   # visit counts truncated to [2, m_cap]
tau      <- 0.5                   # noise standard deviation

## ---- 1.3 the test --------------------------------------------------------
METHODS  <- c(sp = "spacing", g = "$g$ known")   # sp is the BENCHMARK
R_BASIS  <- 0L                    # augmentation order of the basis B_{r,K}:
# the test uses the plain half-cosine system B_0, the augmentation by t and
# t^2 controlling a boundary bias that plays no role here
KBAR_MAX <- 25L                   # the K-bar of the paper: omnibus index set
# K_M = {1, ..., KBAR_MAX}, hence n_basis(KBAR_MAX) = R_BASIS + KBAR_MAX
# basis functions in the statistic
KBAR_MU3 <- c(cc, cc + 1L, KBAR_MAX)   # the three panels of the mu_3 figure
B_TEST   <- 999L                  # bootstrap draws; (B+1)(1-alpha) = 950
ALPHA    <- 0.05                  # nominal level (reporting layer only)
LAW      <- "mammen"              # multipliers: "mammen" (two-point) or "normal"
CENTER   <- FALSE                 # TRUE: centre the profiles by their pooled
# mean before the statistic and its bootstrap are built (see the header)
AMP1 <- c(0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50)  # amplitudes, mu_1
AMP3 <- c(0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50)  # amplitudes, mu_3
XCUT1    <- 0.15                  # the mu_1 panels, fig_test_power_mu1 and
# fig_test_cmp_mu1_n<n>, are drawn on the amplitudes <= XCUT1 only; the
# larger ones stay in the file. Reporting layer only: NA draws them all.

## ---- 1.4 checkpoints -----------------------------------------------------
## CKPT_EVERY is the number of replications between two writes of the .rds
## file; at 1 the file is up to date after every replication and nothing is
## ever lost, at the price of one small write per replication.
## RESUME = TRUE makes run_all() continue an existing file instead of
## starting again; set it to FALSE to ignore it and overwrite.
CKPT_EVERY <- 1L
RESUME     <- TRUE

## ---- 1.5 grids, output, seed ---------------------------------------------
NU_EVAL  <- 50L                   # u-grid on which the profiles are recorded
SETTINGS <- Sys.getenv("SIM_SETTINGS", unset = "FULL")   # "QUICK" | "FULL"
OUTDIR   <- Sys.getenv("SIM_OUTDIR", unset = "./output")
LATEX    <- Sys.getenv("SIM_LATEX",  unset = "pdflatex")  # compiles the panels
BASESEED <- 2026L

if (SETTINGS != "FULL") {
  N_GRID <- c(50L, 100L); R_REP <- 20L; B_TEST <- 199L
}
## phi_mat() builds the half-cosine system only, so R_BASIS is fixed at 0
## here; it is kept as a name, and not written as a literal zero, so that
## the K / dimension distinction stays visible in the code.
stopifnot(R_BASIS == 0L)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# =========================================================================
#  2. MEAN FUNCTIONS
#     Each maps two grids (u, t) to the matrix mu(u[i], t[j]).
# =========================================================================

MU_1 <- function(u, t) {
  S <- matrix(0, length(u), length(t))
  for (p in 1:P) S <- S + ((-1)^p / p^2) * cos(outer(cc * pi * u, pi * p * t, "+"))
  S * matrix(a * t + b * t^2, length(u), length(t), byrow = TRUE)
}

MU_3 <- function(u, t) {
  Su <- rowSums(sapply(1:P, function(p) ((-1)^p / p^2) * cos(pi * p * u)))
  matrix(A3 * t + B3 * t^2, length(u), length(t), byrow = TRUE) +
    outer(Su, cos(cc * pi * t))
}


# =========================================================================
#  3. DATA-GENERATING PROCESS  (no mean: the alternatives are added later)
# =========================================================================

G_DENSITY <- function(t) 1 + 0.6 * cos(2 * pi * t)   # visit-time density g
G_MAX     <- 1.6                                     # any upper bound of g
L_MODES   <- 100L                                    # modes of the subject field
L_NOISE   <- 50L                                     # terms of the error field

## sample_T(): m visit times with density g, by rejection sampling.
sample_T <- function(m) {
  out <- numeric(0)
  while (length(out) < m) {
    x <- runif(2 * m); acc <- runif(2 * m)
    out <- c(out, x[acc <= G_DENSITY(x) / G_MAX])
  }
  sort(out[1:m])
}

## The cosine matrices in u and the decay matrix do not depend on the
## subject, so they are built once per u-grid and cached.
XCACHE <- new.env(parent = emptyenv())
x_cache <- function(ug) {
  key <- paste0("u", length(ug))
  if (is.null(XCACHE[[key]]))
    XCACHE[[key]] <- list(
      CU = outer(ug, seq_len(L_MODES), function(u, l) cos(l * pi * u)),
      D  = outer(seq_len(L_MODES), seq_len(L_MODES), function(l, lp) 1 / (l^2 * lp^2)),
      CE = outer(ug, seq_len(L_NOISE), function(u, l) cos(l * pi * u)))
  XCACHE[[key]]
}

## gen_pool(): one data set with NO mean, Y = X_i + tau eta_ij, together
## with the pooled visit times and the subject labels.
gen_pool <- function(n, ug) {
  ca <- x_cache(ug)
  ve <- as.vector((ca$CE^2) %*% (seq_len(L_NOISE)^(-8)))   # Var eta(u) = 1
  Tp <- vector("list", n); Yc <- vector("list", n); subj <- integer(0)
  for (i in seq_len(n)) {
    mi <- 0L
    while (mi < 2L || mi > m_cap) mi <- rpois(1L, lambda_m)  # 2 <= m_i <= m_cap
    Ti <- sample_T(mi)
    Ci <- matrix(rnorm(L_MODES^2), L_MODES, L_MODES) * ca$D
    Xi <- (ca$CU %*% Ci) %*%
      t(outer(Ti, seq_len(L_MODES), function(t, lp) cos(lp * pi * t)))
    Ei <- (ca$CE %*% (matrix(rnorm(L_NOISE * mi), L_NOISE, mi) /
                        seq_len(L_NOISE)^4)) / sqrt(ve)
    Yc[[i]] <- Xi + tau * Ei
    Tp[[i]] <- Ti; subj <- c(subj, rep(i, mi))
  }
  list(T = unlist(Tp), Y = do.call(cbind, Yc), subj = subj, n = n)
}


# =========================================================================
#  4. THE TWO WEIGHT SCHEMES, THE STATISTIC AND ITS BOOTSTRAP
# =========================================================================

## trapz_w(): quadrature weights on the u-grid.
trapz_w <- function(u) {
  n <- length(u); w <- numeric(n)
  w[1] <- (u[2] - u[1]) / 2; w[n] <- (u[n] - u[n - 1]) / 2
  if (n > 2) w[2:(n - 1)] <- (u[3:n] - u[1:(n - 2)]) / 2
  w
}

## center_profiles(): the profiles minus their pooled mean,
## Y_ij(u) - Ybar(u) with Ybar(u) = (1/M) sum_{i,j} Y_ij(u), Y being the
## nu x M matrix of profiles (rows: u-grid, columns: observations; rowMeans
## is recycled down the columns). With center = FALSE, Y is returned
## untouched, with no arithmetic done on it.
center_profiles <- function(Y, center = CENTER)
  if (isTRUE(center)) Y - rowMeans(Y) else Y

## n_basis(): the number of functions of B_{r,K}, that is r + K, the ONLY
## place in the file where the r monomials are added to the K cosines.
## Every argument below that is a number of basis functions comes from it,
## so that a truncation K, or a K-bar, is never read as a dimension.
n_basis <- function(K, r = R_BASIS) as.integer(r + K)

## phi_mat(): the first 'nb' functions of the half-cosine system, phi_1 = 1
## and phi_k(t) = sqrt2 cos((k-1) pi t). Its argument is a NUMBER OF BASIS
## FUNCTIONS, n_basis(K) at the truncation K of the paper, and not K.
## Returns length(t) x nb.
phi_mat <- function(t, nb) {
  Pm <- matrix(1, length(t), nb)
  if (nb >= 2) for (k in 2:nb) Pm[, k] <- sqrt(2) * cos((k - 1) * pi * t)
  Pm
}

## h_of_M(): the window half-width of the spacing scheme,
## h = ceil((1 + ln ln(M+20))/2); it grows so slowly that h = 2 for every
## realistic sample size.
h_of_M <- function(M) ceiling(0.5 + 0.5 * log(log(M + 20)))

## w_spacing(): the BENCHMARK weights, an M x nb matrix, nb a number of
## basis functions, whose entry (ij, k) is the window integral of phi_k
## attached to the observation of rank l, divided by 2h. The windows are
## defined by symmetric reflection at the two ends of the pooled design,
## and the integrals are closed-form: int phi_1 = hi - lo and
## int phi_k = sqrt2 (sin(lam hi) - sin(lam lo))/lam with lam = (k-1) pi.
## No design density is used.
w_spacing <- function(Tv, nb) {
  M <- length(Tv); h <- h_of_M(M); ord <- order(Tv); Ts <- Tv[ord]
  Text <- function(rk) {                          # reflected order statistics
    out <- numeric(length(rk)); lo <- rk <= 0; hi <- rk > M; in_ <- !lo & !hi
    out[in_] <- Ts[rk[in_]]
    out[lo]  <- 2 * Ts[1] - Ts[2 - rk[lo]]
    out[hi]  <- 2 * Ts[M] - Ts[2 * M - rk[hi]]
    out
  }
  l  <- seq_len(M); lo <- Text(l - h); hi <- Text(l + h)
  Wo <- matrix(0, M, nb); Wo[, 1] <- (hi - lo) / (2 * h)
  if (nb >= 2) for (k in 2:nb) {
    lam <- (k - 1) * pi
    Wo[, k] <- sqrt(2) * (sin(lam * hi) - sin(lam * lo)) / (lam * 2 * h)
  }
  W <- matrix(0, M, nb); W[ord, ] <- Wo              # back to the input order
  W
}

## w_obs(): the COMPETITOR weights, phi_k(T_ij) / (M g(T_ij)), i.e. the OBS
## weights of the paper with the design density taken as known. As above,
## nb is a number of basis functions.
w_obs <- function(Tv, nb)
  (1 / length(Tv)) * (phi_mat(Tv, nb) / G_DENSITY(Tv))

## weights_of(): the weight matrix of the scheme 'e', always carrying the
## n_basis(KBAR_MAX) functions of B_{R_BASIS, KBAR_MAX}, a truncation being
## imposed later by gram0().
weights_of <- function(e, Tv)
  if (e == "sp") w_spacing(Tv, n_basis(KBAR_MAX)) else
    w_obs(Tv, n_basis(KBAR_MAX))

## subject_blocks(): the (n_basis(KBAR_MAX) * nu) x n matrix whose column i
## stacks, function by function, the subject contribution v_ik on the
## u-grid, already multiplied by sqrt(w_u) so that the Euclidean inner
## product of two columns is the L^2(U) inner product summed over the basis
## functions. Restricting to the first n_basis(kbar) * nu rows restricts
## the statistic to the index set K_M = {1, ..., kbar}.
subject_blocks <- function(Y, W, subj, n, ug, wu) {
  nu <- length(ug)
  Ys <- sqrt(wu) * Y
  Vs <- matrix(0, n_basis(KBAR_MAX) * nu, n)
  for (k in seq_len(n_basis(KBAR_MAX)))
    Vs[((k - 1) * nu + 1):(k * nu), ] <-
    t(rowsum(t(Ys * rep(W[, k], each = nu)), subj))
  Vs
}

## gram0(): the n x n matrix of inner products <v_i, v_i'> restricted to
## the index set K_M = {1, ..., kbar}, i.e. to the first n_basis(kbar)
## functions of the system, or its polarised form when two block matrices
## are given, WITH THE DIAGONAL SET TO ZERO, so that sum(G) is the
## between-subject statistic and t(xi) G xi its bootstrap version.
gram0 <- function(A, Bm = NULL, kbar = KBAR_MAX, nu = NU_EVAL) {
  rows <- seq_len(n_basis(kbar) * nu)
  A <- A[rows, , drop = FALSE]
  G <- if (is.null(Bm)) crossprod(A) else {
    Bm <- Bm[rows, , drop = FALSE]
    crossprod(A, Bm) + crossprod(Bm, A)
  }
  diag(G) <- 0
  G
}

## draw_xi(): centred unit-variance multipliers, one per subject. The
## two-point law of Mammen also matches the third moment.
MAMMEN <- list(lo = -(sqrt(5) - 1) / 2, hi = (sqrt(5) + 1) / 2,
               p  = (sqrt(5) + 1) / (2 * sqrt(5)))        # P(xi = lo)
draw_xi <- function(law, n, B) switch(law,
                                      normal = matrix(rnorm(n * B), n, B),
                                      mammen = matrix(ifelse(runif(n * B) < MAMMEN$p, MAMMEN$lo, MAMMEN$hi), n, B),
                                      stop("unknown multiplier law"))

## boot_quad(): the N_b bootstrap quadratic forms t(xi) G xi.
boot_quad <- function(G, Xi) colSums(Xi * (G %*% Xi))

## boot_sd(): the EXACT conditional standard deviation of t(xi) G xi given
## the data, for any centred unit-variance multipliers and a symmetric G
## with a zero diagonal: Var* = 2 sum_{i != i'} G_ii'^2. Used only to check
## the Monte Carlo bootstrap in check_test (b).
boot_sd <- function(G) sqrt(2 * sum(G^2))

## below(): the number of bootstrap draws strictly below the observed
## statistic. This is what the experiment stores; a replication rejects
## when below() reaches floor((N_b+1)(1-alpha)), which is exactly
## T_n > the ((N_b+1)(1-alpha))-th order statistic of the bootstrap sample.
below <- function(Tobs, Tstar) sum(Tstar < Tobs)
crit_rank <- function(B, alpha) floor((B + 1) * (1 - alpha))


# =========================================================================
#  5. SELF-TESTS
#     (a) v_i is linear in the profiles, so the whole amplitude grid is
#         reached exactly from two block matrices, V(theta) = V0 + theta Vm;
#     (b) the bootstrap law is centred with variance 2 sum_{i!=i'} G_ii'^2,
#         which is what boot_sd() returns in closed form;
#     (c) the two weight schemes estimate the same coefficients: on a
#         noiseless pool both reproduce beta_k = int mu_1(u,.) phi_k;
#     (d) the Gram matrix built directly at one amplitude and the
#         three-term recomposition the experiment uses agree exactly;
#     (e) with CENTER = TRUE, the centred profiles have a zero pooled mean.
#     Checks (a), (b) and (d) run on the profiles as the experiment uses
#     them, centred or not; with centring, (a) checks that the centring
#     commutes with the decomposition Y + theta mu_1.
# =========================================================================

check_test <- function(n = 40L, B = 20000L, tol = 1e-9) {
  set.seed(BASESEED + 99L)
  ug <- seq(0, 1, length.out = NU_EVAL); wu <- trapz_w(ug)
  po <- gen_pool(n, ug); Tv <- po$T; th <- 0.37
  if (isTRUE(CENTER)) {
    er <- max(abs(rowMeans(center_profiles(po$Y + th * MU_1(ug, Tv)))))
    cat(sprintf("  (e) max_u |mean of the centred profiles| = %.2e\n", er))
    stopifnot(er < tol)
  }
  for (e in names(METHODS)) {
    W  <- weights_of(e, Tv)
    V0 <- subject_blocks(center_profiles(po$Y), W, po$subj, n, ug, wu)
    Vm <- subject_blocks(center_profiles(MU_1(ug, Tv)), W, po$subj, n, ug, wu)
    Vd <- subject_blocks(center_profiles(po$Y + th * MU_1(ug, Tv)),
                         W, po$subj, n, ug, wu)
    er <- max(abs(Vd - (V0 + th * Vm)))
    cat(sprintf("  (a) [%s] max|V(theta) - (V0 + theta Vm)| = %.2e\n", e, er))
    stopifnot(er < tol)
    ## (d) direct Gram at theta against the polarised recomposition
    Gd <- gram0(V0 + th * Vm)
    rc <- sum(gram0(V0)) + th * sum(gram0(V0, Vm)) + th^2 * sum(gram0(Vm))
    rl <- abs(sum(Gd) - rc) / abs(sum(Gd))
    cat(sprintf("  (d) [%s] |T(theta) - recomposed| / |T(theta)| = %.2e\n", e, rl))
    stopifnot(rl < 1e-10)
    if (e == "sp") {
      G  <- gram0(V0 + th * Vm)
      TS <- boot_quad(G, draw_xi(LAW, n, B))
      cat(sprintf("  (b) E*[T*] = %+.3e (0)   sd*/target = %.4f\n",
                  mean(TS), sd(TS) / boot_sd(G)))
      stopifnot(abs(mean(TS)) < 6 * sd(TS) / sqrt(B),
                abs(sd(TS) / boot_sd(G) - 1) < 0.05)
    }
  }
  ## (c) both schemes on a noiseless pool of 500 subjects, first modes
  set.seed(BASESEED + 98L)
  ug <- seq(0, 1, length.out = 201L)
  po <- gen_pool(500L, ug); Tv <- po$T
  Kb <- 6L; nb <- n_basis(Kb)          # truncation of the paper, and its
  tq <- seq(0, 1, length.out = 4001L)  # number of basis functions
  wq <- c(1, rep(c(4, 2), length.out = length(tq) - 2L), 1) / (3 * (length(tq) - 1))
  bt <- (MU_1(ug, tq) * rep(wq, each = length(ug))) %*% phi_mat(tq, nb)  # true beta_k
  Ym <- MU_1(ug, Tv)
  nrm <- sqrt(colSums(bt^2 * trapz_w(ug)))
  for (e in names(METHODS)) {
    W  <- if (e == "sp") w_spacing(Tv, nb) else w_obs(Tv, nb)
    bh <- Ym %*% W
    er <- sqrt(colSums((bh - bt)^2 * trapz_w(ug))) / max(nrm)
    cat(sprintf("  (c) [%s] error on beta_1..beta_%d, in units of max_k ||beta_k||: %s\n",
                e, nb, paste(sprintf("%.1e", er), collapse = " ")))
    stopifnot(max(er) < 0.05)
  }
  invisible(TRUE)
}


# =========================================================================
#  6. STORAGE  (one .rds, written ATOMICALLY, every figure and the table
#               being built afterwards from that file alone)
# =========================================================================

rds_file <- function(name) file.path(OUTDIR, paste0(name, ".rds"))

## save_rds(): a write that a power cut cannot corrupt. The object goes to
## <name>.rds.tmp; only once that file is complete is the current
## <name>.rds renamed <name>.rds.bak and the temporary file moved in its
## place. At every instant at least one of the two files is a complete
## .rds, and load_rds() falls back on the backup if the main one is
## unreadable. Both moves are renames, so the cost of a checkpoint does not
## grow with the number of checkpoints already written.
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
## they are of no further use. <name>.rds.tmp is only ever present as the
## leftover of a write that was cut, and goes unconditionally. The backup
## goes only after <name>.rds has been READ BACK in full: nothing is
## dropped on the strength of a write that may not have landed, and a run
## that ends short of its replications keeps its backup.
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

## rep_seed(): the seed of replication 'r' at sample size 'n'. It depends
## on nothing else, so a replication computed after a restart, or in
## another order, or in a later run that only extends the number of
## replications, is the very one an uninterrupted run would have produced.
rep_seed <- function(n, r) as.integer(BASESEED + 7L + 100000L * n + r)


# =========================================================================
#  7. POWER EXPERIMENT
#
#  Both T_n and T_n^* are QUADRATIC in the amplitude theta, since
#  v_i(theta) = v_i(0) + theta v_i(mu). The whole amplitude grid is
#  therefore reached from three Gram matrices, and the bootstrap from three
#  quadratic forms, exactly and not approximately. The same replication
#  serves the two weight schemes, the two alternatives and the three values
#  of kbar, so that every curve below is paired. At theta = 0 the statistic
#  does not depend on the alternative, so the level depends on the scheme
#  and on kbar only.
#
#  run_all() stores, for every replication, the number of bootstrap draws
#  below the observed statistic, and rewrites res_test.rds every
#  CKPT_EVERY replications; write_outputs() turns the counts into rejection
#  frequencies at the nominal level ALPHA, on the file, complete or not.
# =========================================================================

## test_sig(): the parameters two runs must share for their replications to
## be comparable, hence for a stored file to be continued, the centring
## included. N_GRID and R_REP are NOT in it: sample sizes and replications
## can be added. ALPHA and XCUT1 are not either: they are applied when the
## counts are read.
test_sig <- function() list(
  methods = METHODS, a = a, b = b, cc = cc, P = P, A3 = A3, B3 = B3,
  lambda_m = lambda_m, m_cap = m_cap, tau = tau,
  L_modes = L_MODES, L_noise = L_NOISE, nu_eval = NU_EVAL,
  kbar_max = KBAR_MAX, kbar_mu3 = KBAR_MU3, amp1 = AMP1, amp3 = AMP3,
  B = B_TEST, law = LAW, baseseed = BASESEED, center = isTRUE(CENTER))

## empty_slot(): the storage of one sample size, all replications pending.
empty_slot <- function(n, ne, nk, na1, na3, R) list(
  n = as.integer(n), R_done = 0L, secs = 0,
  M_rep = rep(NA_real_, R),
  cnt1 = array(NA_integer_, c(R, ne, na1)),        # counts, mu_1
  cnt3 = array(NA_integer_, c(R, ne, nk, na3)))    # counts, mu_3

## resize_slot(): the same slot with room for R replications, the completed
## ones being carried over. R may grow (the new rows are pending) or shrink
## (the last replications are dropped).
resize_slot <- function(z, ne, nk, na1, na3, R) {
  keep <- seq_len(min(z$R_done, R))
  w <- empty_slot(z$n, ne, nk, na1, na3, R)
  if (length(keep)) {
    w$M_rep[keep]   <- z$M_rep[keep]
    w$cnt1[keep, , ]   <- z$cnt1[keep, , , drop = FALSE]
    w$cnt3[keep, , , ] <- z$cnt3[keep, , , , drop = FALSE]
  }
  w$R_done <- length(keep); w$secs <- z$secs
  w
}

## test_status(): what the stored file holds, without running anything.
test_status <- function() {
  if (!has_rds("res_test")) { cat("res_test.rds: absent\n"); return(invisible(NULL)) }
  S <- load_rds("res_test")
  cat(sprintf("res_test.rds: R = %d asked for, N_b = %d\n", S$R, S$B))
  for (z in S$by_n)
    cat(sprintf("   n = %5d : %4d / %4d replications done (%.1f min)\n",
                z$n, z$R_done, S$R, z$secs / 60))
  invisible(S)
}

run_all <- function(n_grid = N_GRID, R = R_REP, resume = RESUME,
                    every = CKPT_EVERY) {
  ug <- seq(0, 1, length.out = NU_EVAL); wu <- trapz_w(ug)
  ne <- length(METHODS); nk <- length(KBAR_MU3)
  na1 <- length(AMP1); na3 <- length(AMP3)
  ik1 <- which(KBAR_MU3 == KBAR_MAX)               # panel carrying kbar max
  jc  <- crit_rank(B_TEST, ALPHA)                  # for the progress lines
  sig <- test_sig()
  
  ## ---- the stored file, if it is to be continued -------------------------
  ST <- NULL
  if (resume && has_rds("res_test")) {
    ST <- load_rds("res_test")
    if (!isTRUE(all.equal(ST$sig, sig)))
      stop("res_test.rds was produced under other parameters: move it aside, ",
           "or call run_all(resume = FALSE) to overwrite it.")
    old <- ST$by_n
    names(old) <- vapply(old, function(z) as.character(z$n), "")
    ## sample sizes already stored are kept, even if they left the grid
    all_n <- sort(unique(c(as.integer(n_grid), as.integer(names(old)))))
    ST$by_n <- lapply(all_n, function(nn) {
      k <- as.character(nn)
      if (!is.null(old[[k]])) resize_slot(old[[k]], ne, nk, na1, na3, R)
      else empty_slot(nn, ne, nk, na1, na3, R)
    })
    ST$n_grid <- all_n; ST$R <- R
    cat(sprintf("  resuming: %s\n", paste(vapply(ST$by_n, function(z)
      sprintf("n=%d %d/%d", z$n, z$R_done, R), ""), collapse = ", ")))
  }
  if (is.null(ST)) {
    all_n <- sort(unique(as.integer(n_grid)))
    ST <- list(
      methods = METHODS, n_grid = all_n, R = R, B = B_TEST, law = LAW,
      kbar_max = KBAR_MAX, kbar_mu3 = KBAR_MU3, amp1 = AMP1, amp3 = AMP3,
      lambda_m = lambda_m, m_cap = m_cap, tau = tau, nu_eval = NU_EVAL,
      sig = sig,
      by_n = lapply(all_n, empty_slot, ne = ne, nk = nk, na1 = na1,
                    na3 = na3, R = R))
  }
  
  todo <- which(ST$n_grid %in% as.integer(n_grid))
  for (jn in todo) {
    n <- ST$n_grid[jn]; z <- ST$by_n[[jn]]
    if (z$R_done >= R) { cat(sprintf("    n = %4d   already complete\n", n)); next }
    
    for (r in seq.int(z$R_done + 1L, R)) {
      t_rep <- proc.time()[3]
      ## the seed of THIS replication: the run is reproducible whatever the
      ## point at which it was interrupted
      set.seed(rep_seed(n, r))
      
      po <- gen_pool(n, ug); Tv <- po$T; z$M_rep[r] <- length(Tv)
      Xi <- draw_xi(LAW, n, B_TEST)                # one draw, shared
      ## the profiles as the statistic sees them: centred by their pooled
      ## mean if CENTER, each part separately since the centring is linear
      Y0  <- center_profiles(po$Y)
      Mu1 <- center_profiles(MU_1(ug, Tv)); Mu3 <- center_profiles(MU_3(ug, Tv))
      
      for (e in seq_len(ne)) {
        W  <- weights_of(names(METHODS)[e], Tv)
        V0 <- subject_blocks(Y0,   W, po$subj, n, ug, wu)   # null part
        V1 <- subject_blocks(Mu1,  W, po$subj, n, ug, wu)   # direction 1
        V3 <- subject_blocks(Mu3,  W, po$subj, n, ug, wu)   # direction 3
        
        for (l in seq_len(nk)) {
          kb <- KBAR_MU3[l]
          G0 <- gram0(V0, kbar = kb)               # null Gram at this kbar
          u0 <- boot_quad(G0, Xi); T0 <- sum(G0)
          ## --- alternative mu_3, at every kbar --------------------------
          G1 <- gram0(V0, V3, kbar = kb); G2 <- gram0(V3, kbar = kb)
          u1 <- boot_quad(G1, Xi); u2 <- boot_quad(G2, Xi)
          s1 <- sum(G1); s2 <- sum(G2)
          for (q in seq_along(AMP3))
            z$cnt3[r, e, l, q] <- below(T0 + AMP3[q] * s1 + AMP3[q]^2 * s2,
                                        u0 + AMP3[q] * u1 + AMP3[q]^2 * u2)
          ## --- alternative mu_1, at the omnibus kbar only ---------------
          if (l == ik1) {
            H1 <- gram0(V0, V1, kbar = kb); H2 <- gram0(V1, kbar = kb)
            w1 <- boot_quad(H1, Xi); w2 <- boot_quad(H2, Xi)
            t1 <- sum(H1); t2 <- sum(H2)
            for (q in seq_along(AMP1))
              z$cnt1[r, e, q] <- below(T0 + AMP1[q] * t1 + AMP1[q]^2 * t2,
                                       u0 + AMP1[q] * w1 + AMP1[q]^2 * w2)
          }
        }
      }
      
      ## ---- the replication is complete: it can now be checkpointed -------
      z$R_done <- r; z$secs <- z$secs + (proc.time()[3] - t_rep)
      ST$by_n[[jn]] <- z
      ckpt <- (r %% every == 0L) || (r == R)
      if (ckpt) save_rds("res_test", ST, quiet = TRUE)
      if (r %% 50 == 0 || r == 1 || r == R)
        cat(sprintf("    n = %4d   rep %4d / %4d   (%.1f s/rep%s)\n", n, r, R,
                    z$secs / z$R_done, if (ckpt) ", saved" else ""))
    }
    
    ST$by_n[[jn]] <- z
    save_rds("res_test", ST, quiet = TRUE)
    idx <- seq_len(z$R_done)
    ## counts at one method and one amplitude, as a (replications x kbar)
    ## matrix whatever the number of replications already done
    kb_mat <- function(e, q) matrix(z$cnt3[idx, e, , q], length(idx), nk)
    for (e in seq_len(ne))
      cat(sprintf("  n = %4d | %-8s | level %s | power at theta = %.2f: mu_1 %.3f, mu_3 %s\n",
                  n, names(METHODS)[e],
                  paste(sprintf("%.3f", colMeans(kb_mat(e, 1L) >= jc)), collapse = "/"),
                  tail(AMP1, 1), mean(z$cnt1[idx, e, na1] >= jc),
                  paste(sprintf("%.3f", colMeans(kb_mat(e, na3) >= jc)), collapse = "/")))
  }
  
  save_rds("res_test", ST)
  ## An interrupted run never reaches this line, so its backup survives on
  ## disk, which is the whole point of keeping one. The test is a guard for
  ## the remaining case, a sample size left short of its R replications.
  if (all(vapply(ST$by_n[todo], function(z) z$R_done >= R, logical(1))))
    clean_ckpt("res_test")
  else
    cat("  run still partial: res_test.rds.bak kept\n")
  invisible(ST)
}


# =========================================================================
#  8. OUTPUT HELPERS, FIGURES AND TABLE
# =========================================================================

## Every figure file is a compilable document on its own, opened and closed
## by these two blocks, and is compiled by emit_fig() right after being
## written. The paper includes the resulting PDF, e.g.
##   \includegraphics[width=0.32\linewidth]{fig_test_power_mu3_kbar4} ...
TEX_HEAD <- c("\\documentclass[border=2pt]{standalone}",
              "\\usepackage[T1]{fontenc}",
              "\\usepackage{amsmath,amssymb}",
              "\\usepackage{tikz}",
              "\\begin{document}")
TEX_FOOT <- "\\end{document}"

## check_graphics(): the figure toolchain is needed by write_outputs() only,
## and is therefore checked at the first panel rather than when the file is
## sourced, so that the Monte Carlo can be run, and resumed, on a machine
## that has neither tikzDevice nor LaTeX.
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

## compile_tex(): run LATEX on OUTDIR/<name>.tex and keep only the PDF. The
## compilation is done from inside OUTDIR, so that a path with spaces or
## accents never reaches the command line. On failure the errors of the log
## are shown and the log is kept.
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

VAL_PCH <- c(19, 17, 15, 18); VAL_LTY <- c(1, 2, 5, 4)   # one per n
MET_PCH <- c(19, 1); MET_LTY <- c(1, 2)                  # one per weight scheme
MET_COL <- c(sp = "black", g = "grey45")
PT_CEX <- 1.05; LN_LWD <- 1.5

## empty_panel(): the common frame, amplitude on the x-axis, rejection
## frequency on the y-axis, with the nominal level as a dotted rule.
empty_panel <- function(amp, ylab = "", main = NULL) {
  plot(NA, xlim = range(amp), ylim = c(0, 1.02), las = 1, bty = "l",
       xlab = "amplitude $\\theta$", ylab = ylab, cex.axis = 0.9)
  grid(col = "grey85", lty = 3, lwd = 0.6)
  abline(h = ALPHA, lty = 3, col = "grey45")
  if (!is.null(main)) mtext(main, side = 3, line = 0.4, cex = 0.85)
}

## power_panel(): the SPACING statistic alone, one curve per number of
## subjects.
power_panel <- function(amp, Y, n_grid, ylab = "", main = NULL,
                        legend = TRUE) {
  empty_panel(amp, ylab, main)
  for (v in seq_along(n_grid)) {
    lines(amp, Y[v, ], lty = VAL_LTY[v], lwd = LN_LWD)
    points(amp, Y[v, ], pch = VAL_PCH[v], cex = PT_CEX)
  }
  if (legend)
    legend("right", sprintf("$n=%d$", n_grid),
           pch = VAL_PCH[seq_along(n_grid)], lty = VAL_LTY[seq_along(n_grid)],
           lwd = LN_LWD, bty = "n", cex = 0.82, seg.len = 2.6)
}

## cmp_panel(): one number of subjects, the two weight schemes against each
## other. Y is (number of schemes) x (number of amplitudes).
cmp_panel <- function(amp, Y, met, ylab = "", main = NULL, legend = TRUE) {
  empty_panel(amp, ylab, main)
  for (e in seq_along(met)) {
    lines(amp, Y[e, ], lty = MET_LTY[e], lwd = LN_LWD, col = MET_COL[e])
    points(amp, Y[e, ], pch = MET_PCH[e], cex = PT_CEX, col = MET_COL[e])
  }
  if (legend)
    legend("right", unname(met), pch = MET_PCH[seq_along(met)],
           lty = MET_LTY[seq_along(met)], col = MET_COL[seq_along(met)],
           lwd = LN_LWD, bty = "n", cex = 0.82, seg.len = 2.6)
}

## amp_keep(): the amplitudes a panel is drawn on. XCUT1 cuts the abscissa
## of the mu_1 panels; NA, or a cut that would empty the panel, keeps them
## all. The stored counts are untouched: this is the reporting layer.
amp_keep <- function(amp, cut) {
  if (is.null(cut) || is.na(cut)) return(seq_along(amp))
  k <- which(amp <= cut)
  if (length(k) < 2L) seq_along(amp) else k
}

## write_outputs(): the counts become rejection frequencies here, at the
## level ALPHA. Only the R_done first replications of each sample size are
## read, so the figures and the table can be built while the experiment is
## still running, or after an interruption; a sample size with no
## replication yet is left out.
write_outputs <- function() {
  cat("\n[2] figures and table from res_test.rds\n")
  S <- load_rds("res_test")
  slots <- Filter(function(z) z$R_done > 0L, S$by_n)
  if (!length(slots)) stop("res_test.rds holds no completed replication yet.")
  n_grid <- vapply(slots, function(z) as.integer(z$n), integer(1))
  Rd     <- vapply(slots, function(z) as.integer(z$R_done), integer(1))
  if (any(Rd < S$R))
    cat(sprintf("  [i] partial file: %s replications (asked: %d)\n",
                paste(sprintf("n=%d:%d", n_grid, Rd), collapse = ", "), S$R))
  nn <- length(slots); ne <- length(S$methods); nk <- length(S$kbar_mu3)
  jc <- crit_rank(S$B, ALPHA)                    # ALPHA lives HERE only
  
  rate1 <- array(NA_real_, c(nn, ne, length(S$amp1)))       # power, mu_1
  rate3 <- array(NA_real_, c(nn, ne, nk, length(S$amp3)))   # power, mu_3
  lev   <- array(NA_real_, c(nn, ne, nk))                   # level, theta = 0
  lev_se <- array(NA_real_, c(nn, ne, nk))
  for (jn in seq_len(nn)) {
    z <- slots[[jn]]; idx <- seq_len(z$R_done)
    rate1[jn, , ] <- apply(z$cnt1[idx, , , drop = FALSE] >= jc, c(2, 3), mean)
    rate3[jn, , , ] <- apply(z$cnt3[idx, , , , drop = FALSE] >= jc, c(2, 3, 4), mean)
    lev[jn, , ]    <- rate3[jn, , , 1]
    lev_se[jn, , ] <- sqrt(lev[jn, , ] * (1 - lev[jn, , ]) / z$R_done)
  }
  
  ## --- the power of the SPACING statistic --------------------------------
  isp <- which(names(S$methods) == "sp")
  k1  <- amp_keep(S$amp1, XCUT1)                 # abscissa cut, mu_1 panels
  if (length(k1) < length(S$amp1))
    cat(sprintf("  [i] mu_1 panels drawn up to theta = %g (%d of %d amplitudes)\n",
                XCUT1, length(k1), length(S$amp1)))
  
  ## fig_test_power_mu1: one panel, kbar = KBAR_MAX
  emit_fig("fig_test_power_mu1", 4.6, 4.2, function() {
    par(mar = c(4.0, 4.4, 1.0, 0.9), mgp = c(2.5, 0.8, 0), tcl = -0.3)
    power_panel(S$amp1[k1], matrix(rate1[, isp, k1], nn, length(k1)), n_grid,
                ylab = "rejection frequency")
  })
  
  ## fig_test_power_mu3_kbar<k>: one file per kbar, full amplitude grid
  for (l in seq_len(nk)) local({
    ll <- l; kb <- S$kbar_mu3[ll]
    emit_fig(sprintf("fig_test_power_mu3_kbar%d", kb), 4.2, 4.2, function() {
      par(mar = c(4.0, 4.4, 1.8, 0.9), mgp = c(2.5, 0.8, 0), tcl = -0.3)
      power_panel(S$amp3, matrix(rate3[, isp, ll, ], nn, length(S$amp3)), n_grid,
                  ylab = "rejection frequency",
                  main = sprintf("$\\bar K = %d$", kb),
                  legend = TRUE #(ll == 1L)
      )
    })
  })
  
  ## --- fig_test_cmp_mu1_n<n>: the two schemes, one panel per n -----------
  ## Direction mu_1 at the omnibus truncation, the configuration of
  ## fig_test_power_mu1; the panels are meant to be read two per line.
  for (jn in seq_len(nn)) local({
    jj <- jn
    emit_fig(sprintf("fig_test_cmp_mu1_n%d", n_grid[jj]), 4.2, 4.0, function() {
      par(mar = c(4.0, 4.4, 1.8, 0.9), mgp = c(2.5, 0.8, 0), tcl = -0.3)
      cmp_panel(S$amp1[k1], matrix(rate1[jj, , k1], ne, length(k1)), S$methods,
                ylab = "rejection frequency",
                main = sprintf("$n = %d$", n_grid[jj]))
    })
  })
  
  ## --- tab_test_level: empirical level of the SPACING statistic, one row
  ## per truncation. The table is there to read the effect of Kbar, so the
  ## competitor is left to the comparison figure.
  Rtex <- if (length(unique(Rd)) == 1L) sprintf("$R=%d$", Rd[1]) else
    sprintf("$R$ between $%d$ and $%d$", min(Rd), max(Rd))
  hdr <- paste0(" & ", paste(sprintf("$%d$", n_grid), collapse = " & "), " \\\\")
  body <- sapply(seq_len(nk), function(l)
    paste0(sprintf("$%d$", S$kbar_mu3[l]), " & ",
           paste(sprintf("$%.3f\\,(%.3f)$", lev[, isp, l], lev_se[, isp, l]),
                 collapse = " & "), " \\\\"))
  writeLines(c("\\begin{table}[H]", "\\centering",
               sprintf("\\begin{tabular}{l%s}", strrep("c", nn)),
               "\\toprule",
               sprintf("\\multirow{2}{*}{$\\bar K$} & \\multicolumn{%d}{c}{$n$} \\\\", nn),
               sprintf("\\cmidrule(lr){2-%d}", nn + 1L),
               hdr,
               "\\midrule", body, "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Empirical level of the test at ",
                              "$\\mathfrak a=%.2f$, i.e.\\ the rejection frequency at ",
                              "$\\theta=0$ of the spacing statistic, one row per ",
                              "truncation $\\bar K$ and one column per number of ",
                              "subjects ($\\tau=1/2$, $N_{\\rm b}=%d$, %s; Monte ",
                              "Carlo s.e.\\ in parentheses).}"), ALPHA, S$B, Rtex),
               "\\label{tab:test_level}", "\\end{table}"),
             file.path(OUTDIR, "tab_test_level.tex"))
  cat("  wrote tab_test_level.tex\n")
  
  invisible(list(rate1 = rate1, rate3 = rate3, lev = lev, lev_se = lev_se,
                 n = n_grid, R_done = Rd, methods = S$methods,
                 kbar = S$kbar_mu3, amp1 = S$amp1, amp3 = S$amp3,
                 amp1_shown = S$amp1[k1]))
}


# =========================================================================
#  9. DRIVER
#
#  Once res_test.rds is there, the figures and the table can be rebuilt
#  alone, by sourcing the script with SIM_NO_RUN defined and calling
#  write_outputs(). Restarting the script after an interruption re-runs
#  the self-tests, which are cheap, and picks the experiment up where its
#  checkpoint left it; test_status() reports, without computing anything,
#  how many replications are stored at each sample size.
# =========================================================================

if (!exists("SIM_NO_RUN")) {
  t0 <- proc.time()[3]
  cat(sprintf("==== SETTINGS=%s | R=%d | N_b=%d | kbar in {%s} | law=%s ====\n",
              SETTINGS, R_REP, B_TEST, paste(KBAR_MU3, collapse = ", "), LAW))
  cat("\n[0] self-tests\n"); check_test()
  
  cat("\n[1] power experiment\n")
  run_all()
  E <- write_outputs()
  
  cat("\n==================== SUMMARY ====================\n")
  for (e in seq_along(E$methods)) {
    cat(sprintf("[%s] empirical level (rows: kbar, columns: n)\n",
                names(E$methods)[e]))
    print(round(`dimnames<-`(t(E$lev[, e, ]), list(E$kbar, E$n)), 3))
    cat(sprintf("[%s] power against mu_1 at kbar = %d (rows: n, columns: amplitude)\n",
                names(E$methods)[e], KBAR_MAX))
    print(round(`dimnames<-`(E$rate1[, e, ], list(E$n, sprintf("%.3g", E$amp1))), 3))
    for (l in seq_along(E$kbar)) {
      cat(sprintf("[%s] power against mu_3 at kbar = %d\n",
                  names(E$methods)[e], E$kbar[l]))
      print(round(`dimnames<-`(E$rate3[, e, l, ], list(E$n, sprintf("%.3g", E$amp3))), 3))
    }
  }
  cat(sprintf("mu_1 panels drawn on theta in {%s}\n",
              paste(sprintf("%g", E$amp1_shown), collapse = ", ")))
  cat(sprintf("\nTotal wall time: %.1f s\nOutputs written to: %s\n",
              proc.time()[3] - t0, normalizePath(OUTDIR)))
}