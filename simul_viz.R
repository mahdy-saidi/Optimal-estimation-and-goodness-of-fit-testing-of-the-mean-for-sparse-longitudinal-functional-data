##########################################################################
#                 FIGURES THAT NEED NO MONTE CARLO
#   Optimal estimation and goodness-of-fit testing of the mean for
#   sparse longitudinal functional data
#
#   This file is SELF-CONTAINED: it duplicates the few objects it needs
#   (parameters, mean functions, mechanism, candidate bases, figure
#   helpers) and depends on no .rds. It writes
#
#       fig_means_mu1     fig_means_mu2     fig_means_mu3
#       fig_data_field    fig_data_profiles
#       fig_augbasis_phi  fig_augbasis_sup
#
#   ONE PANEL PER FILE, each .tex being a standalone document compiled on
#   the spot by pdflatex.
#
#   NOTATION OF THE PAPER. The candidate bases are B_{r,K}, r = 0, 1, 2
#   (half-cosine; augmented by t; augmented by t and t^2), where K counts
#   the cosine terms c_1, ..., c_K, the r monomials of the augmentation
#   EXCLUDED, so that B_{r,K} has dimension r + K. basis_r() takes that
#   DIMENSION as its argument: B_{r,K} is basis_r(r, r + K). The panels
#   below show the phi_k of the augmented system of Section "Basis choice".
##########################################################################


# =========================================================================
#  USER PARAMETERS
# =========================================================================

## ---- shape parameters of the three mean functions ----------------------
a  <- 1                     # coefficient of t
b  <- -7                    # coefficient of t^2
cc <- 4                     # frequency multiplier (the "c" of the paper)
P  <- 8L                    # number of Fourier terms

## ---- mechanism ----------------------------------------------------------
tau      <- 0.5                    # noise standard deviation (homoscedastic)

## ---- grids of the panels ------------------------------------------------
NU_FINE <- 1001L                   # grid on which the profiles are drawn
NU_VIZ  <- 61L; NT_VIZ <- 61L      # grids of the surface panels: the size of
# a .tex panel is proportional to their
# product, about 140 bytes per facet
VIZ_THIN <- 4L                     # 1 point out of VIZ_THIN is DRAWN in the
# profile panel
VIZ_MI   <- 3L                     # observed profiles shown in fig_data
VIZ_MEAN <- 1L                     # mean used by fig_data (1, 2 or 3)
D_SUP    <- 50L                    # dimension shown in the sup-norm panel,
# i.e. the system B_{2,D_SUP-2}

## ---- typography ---------------------------------------------------------
CEX    <- 1.35                     # multiplier applied to every cex
LN_LWD <- 1.2                      # line width of the curve panels

## ---- output -------------------------------------------------------------
OUTDIR   <- Sys.getenv("SIM_OUTDIR", unset = "./output")
LATEX    <- Sys.getenv("SIM_LATEX",  unset = "pdflatex")
BASESEED <- 2026L

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# =========================================================================
#  MEAN FUNCTIONS
#      mu_1(u,t) = (a t + b t^2) sum_p ((-1)^p / p^2) cos(pi p t + c pi u)
#      mu_2(u,t) = (a t + b t^2) sum_p ((-1)^p / p^2) cos(pi p u + c pi t)
#      mu_3(u,t) = a t + b t^2 + cos(c pi t) sum_p ((-1)^p / p^2) cos(pi p u)
# =========================================================================

QUAD_T <- function(t) a * t + b * t^2

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
  matrix(QUAD_T(t), length(u), length(t), byrow = TRUE) + outer(Su, cos(cc * pi * t))
}

MEANS <- list(mu1 = MU_1, mu2 = MU_2, mu3 = MU_3)


# =========================================================================
#  DATA-GENERATING PROCESS  (only what fig_data needs)
# =========================================================================

G_DENSITY <- function(t) 1 + 0.6 * cos(2 * pi * t)   # g(t) = 1 + 3cos(2 pi t)/5
G_MAX     <- 1.6
L_MODES   <- 100L
L_NOISE   <- 50L

sample_T <- function(m) {                            # m visit times, density g
  out <- numeric(0)
  while (length(out) < m) {
    x <- runif(2 * m); acc <- runif(2 * m)
    out <- c(out, x[acc <= G_DENSITY(x) / G_MAX])
  }
  sort(out[1:m])
}

make_Xi <- function(L_modes = L_MODES) {             # one subject field
  C <- matrix(rnorm(L_modes^2), L_modes, L_modes) /
    outer(seq_len(L_modes), seq_len(L_modes), function(l, lp) l^2 * lp^2)
  function(u, t)
    outer(u, seq_len(L_modes), function(u, l) cos(l * pi * u)) %*% C %*%
    t(outer(t, seq_len(L_modes), function(t, lp) cos(lp * pi * t)))
}

gen_noise <- function(u, nc, Lp = L_NOISE) {         # unit variance at every u
  B  <- matrix(rnorm(Lp * nc), Lp, nc) / (seq_len(Lp)^4)
  Cu <- outer(u, seq_len(Lp), function(u, l) cos(l * pi * u))
  v  <- as.vector((Cu^2) %*% (seq_len(Lp)^(-8)))
  (Cu %*% B) / sqrt(v)
}


# =========================================================================
#  CANDIDATE BASES  --  B_{r,K}, r = 0 (half-cosine), 1 (+ t), 2 (+ t, t^2)
#  Identical to simul_estimation.R: basis_r(r, dim) orthonormalises the
#  first 'dim' elements of the dictionary, and B_{r,K} of the paper is
#  basis_r(r, r + K). Each phi_k is a polynomial of degree at most 2 plus
#  finitely many cosines, so its values are closed-form.
# =========================================================================

basis_r <- function(r, dim) {
  K      <- dim                    # number of dictionary elements retained
  poly   <- matrix(0, K, 3); freqs <- vector("list", K); thetas <- vector("list", K)
  poly[1, ] <- c(1, 0, 0)
  P2 <- c(-sqrt(3),  2 * sqrt(3), 0)
  P3 <- c( sqrt(5), -6 * sqrt(5), 6 * sqrt(5))
  
  if (r == 0L) {
    if (K >= 2) for (k in 2:K) { freqs[[k]] <- (k - 1) * pi; thetas[[k]] <- 1 }
    return(list(poly = poly, freqs = freqs, thetas = thetas, K = K, r = 0L))
  }
  
  j <- seq_len(max(K, 4L))
  al <- -4 * sqrt(6)  / ((2 * j - 1)^2 * pi^2)
  be <-  3 * sqrt(10) / (j^2 * pi^2)
  A  <- c(1, 1 - cumsum(al^2))
  B  <- c(1, 1 - cumsum(be^2))
  
  if (r == 1L) {
    if (K >= 2) poly[2, ] <- P2
    if (K >= 3) for (m in seq_len(K - 2)) {
      k <- 2 + m
      if (m %% 2 == 0) { freqs[[k]] <- m * pi; thetas[[k]] <- 1 } else {
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
  
  if (K >= 2) poly[2, ] <- P2
  if (K >= 3) poly[3, ] <- P3
  if (K >= 4) for (m in seq_len(K - 3)) {
    k <- 3 + m
    if (m %% 2 == 1) {
      nu <- (m + 1) %/% 2; th <- numeric(nu)
      th[nu] <- sqrt(A[nu] / A[nu + 1])
      if (nu >= 2) th[seq_len(nu - 1)] <-
        al[nu] * al[seq_len(nu - 1)] / sqrt(A[nu] * A[nu + 1])
      poly[k, ]  <- (-al[nu] / sqrt(A[nu] * A[nu + 1])) * P2
      freqs[[k]] <- (2 * seq_len(nu) - 1) * pi
    } else {
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


# =========================================================================
#  FIGURE HELPERS  --  one standalone TikZ document per panel
# =========================================================================

TEX_HEAD <- c("\\documentclass[border=2pt]{standalone}",
              "\\usepackage[T1]{fontenc}",
              "\\usepackage{tikz}",
              "\\begin{document}")
TEX_FOOT <- "\\end{document}"

GRAPHICS_OK <- FALSE
check_graphics <- function() {
  if (GRAPHICS_OK) return(invisible(TRUE))
  if (!requireNamespace("tikzDevice", quietly = TRUE))
    stop("package tikzDevice is required to write the figures")
  if (Sys.which(LATEX) == "")
    stop(sprintf("'%s' is not on the PATH: it compiles each figure to PDF", LATEX))
  options(tikzDefaultEngine = "pdftex",
          tikzMetricsDictionary = file.path(OUTDIR, "tikzMetrics"))
  ok <- tryCatch({
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

compile_tex <- function(name) {
  owd <- setwd(OUTDIR); on.exit(setwd(owd), add = TRUE)
  out <- suppressWarnings(
    system2(LATEX, c("-interaction=nonstopmode", "-halt-on-error",
                     shQuote(paste0(name, ".tex"))), stdout = TRUE, stderr = TRUE))
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

## wire(): the grey wireframe surface of the fig_means and fig_data panels.
wire <- function(ug, tg, Z, main)
  persp(ug, tg, Z, theta = 35, phi = 25, expand = 0.72,
        xlab = "$u$", ylab = "$t$", zlab = "", ticktype = "detailed",
        nticks = 4, cex.axis = 0.55 * CEX, cex.lab = 0.80 * CEX,
        border = NA, col = "grey85", shade = 0.55, main = main,
        cex.main = 0.95 * CEX)

## setup_panel(): an empty curve panel, axes and margins already scaled.
setup_panel <- function(xlab, ylab, xlim, ylim) {
  par(mar = c(3.6, 4.2, 0.9, 0.8) * c(CEX, CEX, 1, 1),
      mgp = c(2.5, 0.7, 0) * CEX, tcl = -0.3,
      cex.axis = 0.85 * CEX, cex.lab = CEX)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, las = 1, bty = "l")
}

legend_box <- function(pos, labels, pch, lty, horiz = FALSE)
  legend(pos, legend = labels, pch = pch, lty = lty, lwd = LN_LWD, col = "black",
         bty = "n", horiz = horiz, cex = 0.78 * CEX, seg.len = 1.5)


# =========================================================================
#  THE THREE MEAN SURFACES AND ONE DRAW OF THE MECHANISM
#      fig_means_mu1, fig_means_mu2, fig_means_mu3,
#      fig_data_field, fig_data_profiles
#
#  The seed and the ORDER of the draws are those of simul_estimation.R,
#  so the two data panels are unchanged.
# =========================================================================

fig_means_and_data <- function() {
  cat("\n[fig] mean surfaces and one draw of the mechanism ...\n")
  set.seed(BASESEED + 1L)
  u_viz  <- seq(0, 1, length.out = NU_VIZ)
  t_viz  <- seq(0, 1, length.out = NT_VIZ)
  u_fine <- seq(0, 1, length.out = NU_FINE)
  
  Z  <- lapply(MEANS, function(f) f(u_viz, t_viz))
  Xi <- make_Xi(); Ti <- sample_T(VIZ_MI)
  Zx <- Xi(u_viz, t_viz)
  Yi <- MEANS[[VIZ_MEAN]](u_fine, Ti) + Xi(u_fine, Ti) +
    tau * gen_noise(u_fine, VIZ_MI)
  
  for (q in seq_along(Z)) {
    nm <- names(Z)[q]
    emit_fig(sprintf("fig_means_%s", nm), 3.5, 3.5, function() {
      par(mar = c(0.4, 0.4, 1.5 * CEX, 0.4))
      wire(u_viz, t_viz, Z[[q]], main = sprintf("$\\mu_%s(u,t)$", substr(nm, 3, 3)))
    })
  }
  
  emit_fig("fig_data_field", 3.5, 3.5, function() {
    par(mar = c(0.4, 0.4, 1.5 * CEX, 0.4))
    wire(u_viz, t_viz, Zx, main = "$X_i(u,t)$")
  })
  
  emit_fig("fig_data_profiles", 3.9, 3.5, function() {
    keep <- seq(1L, length(u_fine), by = VIZ_THIN)     # drawn, not computed
    par(mar = c(4.0, 4.2, 1.5, 0.8) * c(CEX, CEX, CEX, 1),
        mgp = c(2.4, 0.8, 0) * CEX)
    ltys <- rep_len(c(1L, 2L, 3L), ncol(Yi))
    matplot(u_fine[keep], Yi[keep, , drop = FALSE], type = "l", lty = ltys,
            lwd = 1.1, col = "black", xlab = "$u$", ylab = "$Y_{ij}(u)$",
            las = 1, bty = "l", main = "observed profiles",
            cex.axis = CEX, cex.lab = CEX, cex.main = 0.95 * CEX)
    legend("bottomright", sprintf("$T_{ij}=%.2f$", Ti), col = "black",
           lty = ltys, lwd = 1.1, bty = "n", cex = 0.70 * CEX)
  })
  invisible(TRUE)
}


# =========================================================================
#  THE AUGMENTED BASIS B_2: functions and sup-norm
#      fig_augbasis_phi : phi_1, ..., phi_5 of B_{2,K}
#      fig_augbasis_sup : ||phi_k||_inf against k, with the 2 sqrt2 level
#
#  Lemma "uniform bound" proves ||phi_k||_inf < 4 sqrt2; the panel shows
#  that the numerical ceiling is 2 sqrt2, approached from below, the
#  supremum being attained at an endpoint of T for every k.
# =========================================================================

exp_augbasis <- function(dim = D_SUP) {
  cat("\n[fig] augmented basis: functions and sup-norm ...\n")
  co  <- basis_r(2L, dim); co1 <- basis_r(1L, dim)
  tt  <- seq(0, 1, length.out = 4001)          # fine grid, endpoints included
  Phi <- basis_phi(co, tt); Phi1 <- basis_phi(co1, tt)
  supn  <- apply(abs(Phi),  2, max)            # ||phi_k||_inf on the grid
  supn1 <- apply(abs(Phi1), 2, max)
  argt  <- tt[apply(abs(Phi), 2, which.max)]   # location of the supremum
  
  emit_fig("fig_augbasis_phi", 3.9, 3.5, function() {
    tp <- seq(0, 1, length.out = 400); Pp <- basis_phi(co, tp)
    setup_panel("$t$", "$\\phi_k(t)$", c(0, 1), range(Pp[, 1:5]) + c(-0.1, 0.7))
    for (k in 1:5) lines(tp, Pp[, k], lty = k, lwd = LN_LWD)
    legend_box("top", sprintf("$\\phi_%d$", 1:5), rep(NA, 5), 1:5, horiz = TRUE)
  })
  
  emit_fig("fig_augbasis_sup", 3.9, 3.5, function() {
    setup_panel("$k$", "$\\|\\phi_k\\|_\\infty$", c(1, dim),
                c(0.95, 2 * sqrt(2) + 0.2))
    abline(h = 2 * sqrt(2), col = "grey55", lty = 2, lwd = 1.0)
    lines(1:dim, supn, lwd = LN_LWD); points(1:dim, supn, pch = 20, cex = 0.55)
    text(dim, 2 * sqrt(2) + 0.07, "$2\\sqrt2$", pos = 2, cex = 0.78 * CEX,
         col = "grey35")
  })
  
  cat(sprintf("  max_k ||phi_k||_inf (k <= %d): B_2 %.6f, B_1 %.6f  (2 sqrt2 = %.6f)\n",
              dim, max(supn), max(supn1), 2 * sqrt(2)))
  cat(sprintf("  sup attained at an endpoint for every k: %s\n",
              all(argt < 1e-6 | argt > 1 - 1e-6)))
  cat(sprintf("  sup-norm nondecreasing in k: %s ; below 2 sqrt2: %s\n",
              all(diff(supn) > -1e-12), all(supn < 2 * sqrt(2))))
  invisible(TRUE)
}


# =========================================================================
#  DRIVER
# =========================================================================

if (!exists("FIG_NO_RUN")) {
  t0 <- proc.time()[3]
  fig_means_and_data()
  exp_augbasis()
  cat(sprintf("\nTotal wall time: %.1f s\nOutputs written to: %s\n",
              proc.time()[3] - t0, normalizePath(OUTDIR)))
}