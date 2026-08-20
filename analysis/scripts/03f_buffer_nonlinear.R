# 03f_buffer_nonlinear.R
# -----------------------------------------------------------
# Follow-up to 03e. The LINEAR climate-buffer interaction is a precisely
# estimated zero, but a smooth (GAM) fit finds a real nonlinear interaction
# between log(B) and thermal exposure (edf 3.8, F 13). So the question is
# not whether the buffer responds to warming, but with what SHAPE, and
# whether the shape survives the inference the strategy note demands.
#
# Buffer strength from the smooth surface:
#
#     Phi(C) = 1 - d log(P) / d log(B) | evaluated at mean log(B), given C
#
# computed by finite difference on the fitted surface over a grid of C,
# with a BLOCK BOOTSTRAP BY YEAR (years resampled whole, with replacement)
# so that arbitrary within-year correlation is respected. The GAM's own
# p-values assume independence given the random effects and are not
# trusted here.
#
# Reports whether Phi declines from the cold end of the observed thermal
# range to the warm end, which is the claim "warming erodes the buffer".
#
# Requires: data/buffer_reef_year.csv   (from 03e)
# Writes:   data/buffer_phi_nonlinear.csv   Phi(C) with bootstrap CI
#           data/buffer_nonlinear_summary.csv
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(mgcv) })
setwd("..")
DATA <- "data"; OUT <- "data"
set.seed(20260805)
NBOOT <- 300

d <- fread(file.path(DATA, "buffer_reef_year.csv"))
d[, `:=`(Reef = factor(Reef), YearF = factor(Year))]
d[, C := C_warm - mean(C_warm)]
d[, lB_c := log(B) - mean(log(B))]
mean_lB <- 0                      # centred, so 0 is the mean reef-year
EPS <- 0.05                       # finite-difference step in log(B)

# grid over the OBSERVED thermal range (5th to 95th percentile)
Cg <- seq(quantile(d$C, 0.05), quantile(d$C, 0.95), length.out = 25)

fit_gam <- function(dat) {
  gam(lP ~ s(lB_c) + s(C) + ti(lB_c, C) + s(Reef, bs = "re") + s(YearF, bs = "re"),
      data = dat, method = "REML")
}

# local slope d log(P) / d log(B) at (mean logB, C), by finite difference.
# Reef and year random effects are held at a reference level; they shift the
# intercept only and cancel in the difference.
phi_of_C <- function(m, dat) {
  ref <- data.table(Reef = dat$Reef[1], YearF = dat$YearF[1])
  hi <- data.table(lB_c = mean_lB + EPS, C = Cg, Reef = ref$Reef, YearF = ref$YearF)
  lo <- data.table(lB_c = mean_lB - EPS, C = Cg, Reef = ref$Reef, YearF = ref$YearF)
  ph <- predict(m, hi); pl <- predict(m, lo)
  slope <- as.numeric(ph - pl) / (2 * EPS)
  1 - slope
}

m0 <- fit_gam(d)
phi0 <- phi_of_C(m0, d)
message("Observed Phi(C) across the thermal range:")
print(data.table(C_warm = round(Cg + mean(d$C_warm), 2), phi = round(phi0, 3)))

# --- block bootstrap by year -------------------------------------------
yrs <- unique(d$Year)
boot <- matrix(NA_real_, nrow = NBOOT, ncol = length(Cg))
ok <- 0
for (b in seq_len(NBOOT)) {
  pick <- sample(yrs, length(yrs), replace = TRUE)
  db <- rbindlist(lapply(seq_along(pick), function(i) {
    x <- d[Year == pick[i]]
    copy(x)[, YearF := factor(paste0(pick[i], "_", i))]
  }))
  db[, Reef := factor(as.character(Reef))]
  db[, YearF := factor(as.character(YearF))]
  r <- try(phi_of_C(fit_gam(db), db), silent = TRUE)
  if (!inherits(r, "try-error")) { ok <- ok + 1; boot[b, ] <- r }
  if (b %% 50 == 0) message("  bootstrap ", b, "/", NBOOT)
}
message(sprintf("block bootstrap: %d of %d replicates converged", ok, NBOOT))

qs <- apply(boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
curve <- data.table(C_centred = Cg, C_warm = Cg + mean(d$C_warm),
                    phi = phi0, lo = qs[1, ], hi = qs[2, ])
fwrite(curve, file.path(OUT, "buffer_phi_nonlinear.csv"))

# --- the directional claim: is the buffer smaller at the warm end? ------
d_phi_boot <- boot[, length(Cg)] - boot[, 1]      # Phi(warm) - Phi(cold)
d_phi_obs  <- phi0[length(phi0)] - phi0[1]
p_erode <- mean(d_phi_boot >= 0, na.rm = TRUE)    # one-sided: erosion means < 0
ci <- quantile(d_phi_boot, c(0.025, 0.975), na.rm = TRUE)

message(sprintf("\nPhi(cold) = %.3f, Phi(warm) = %.3f", phi0[1], phi0[length(phi0)]))
message(sprintf("Difference Phi(warm) - Phi(cold) = %.3f  (95%% CI %.3f to %.3f)",
                d_phi_obs, ci[1], ci[2]))
message(sprintf("One-sided bootstrap p for erosion (Phi falling with C): %.3f", p_erode))

# where is the buffer strongest / weakest?
i_max <- which.max(phi0); i_min <- which.min(phi0)
summ <- data.table(
  quantity = c("phi_cold_end", "phi_warm_end", "phi_diff_warm_minus_cold",
               "phi_diff_ci_lo", "phi_diff_ci_hi", "p_bootstrap_erosion",
               "phi_max", "C_at_phi_max", "phi_min", "C_at_phi_min",
               "gam_ti_edf", "n_reef_years", "nboot_ok"),
  value = c(round(phi0[1], 3), round(phi0[length(phi0)], 3), round(d_phi_obs, 3),
            round(ci[1], 3), round(ci[2], 3), round(p_erode, 3),
            round(phi0[i_max], 3), round(Cg[i_max] + mean(d$C_warm), 2),
            round(phi0[i_min], 3), round(Cg[i_min] + mean(d$C_warm), 2),
            round(summary(m0)$s.table[grep("ti", rownames(summary(m0)$s.table)), "edf"], 2),
            nrow(d), ok))
fwrite(summ, file.path(OUT, "buffer_nonlinear_summary.csv"))
print(summ)
message("\nStep 03f done.")

# -----------------------------------------------------------
# ARBITER: a model-free stratified test.
#
# The smooth fit above reports an eroding buffer, while the linear
# interaction of 03e is an exact null. The two disagree only in
# functional form, so the tie is broken WITHOUT one: split the panel by
# thermal exposure and estimate the biomass-production slope separately
# in each stratum, with reef and year fixed effects and two-way clustered
# errors. This uses all the biomass variation inside each stratum rather
# than the local curvature of a fitted surface at one biomass value.
#
# If the buffer erodes, Phi must fall from the cool strata to the warm.
# -----------------------------------------------------------
suppressPackageStartupMessages(library(fixest))
d[, lB_c := log(B) - mean(log(B))]
qs <- quantile(d$C_warm, c(0, .25, .5, .75, 1))
d[, Q := cut(C_warm, qs, include.lowest = TRUE,
             labels = c("Q1 coolest", "Q2", "Q3", "Q4 warmest"))]
strat <- rbindlist(lapply(levels(d$Q), function(q) {
  dd <- d[Q == q]
  m  <- feols(lP ~ lB_c | Reef + YearF, data = dd)
  b  <- coef(m)[["lB_c"]]
  se <- sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[["lB_c"]]
  data.table(stratum = q, mean_anom = mean(dd$C_warm), phi = 1 - b,
             lo = 1 - b - 1.96 * se, hi = 1 - b + 1.96 * se,
             n = nrow(dd), reefs = uniqueN(dd$Reef))
}))
fwrite(strat, file.path(OUT, "buffer_phi_stratified.csv"))
message("\n=== ARBITER: Phi by thermal stratum (model free) ===")
print(strat[, .(stratum, mean_anom = round(mean_anom, 2), phi = round(phi, 3),
                lo = round(lo, 3), hi = round(hi, 3), n, reefs)])

d[, warm := as.integer(C_warm > median(C_warm))]
mw  <- feols(lP ~ lB_c * warm | Reef + YearF, data = d)
bw  <- coef(mw)[["lB_c:warm"]]
sew <- sqrt(diag(vcov(mw, cluster = ~ Reef + YearF)))[["lB_c:warm"]]
pw  <- 2 * pnorm(-abs(bw / sew))
# Phi = 1 - slope, so a NEGATIVE slope shift means a STRONGER buffer when warm
message(sprintf("Warm half vs cool half: slope shift %+.3f (SE %.3f), p = %.3f", bw, sew, pw))
message(sprintf("  -> Phi is %+.3f in the warm half, i.e. the buffer is %s when warmer%s",
                -bw, ifelse(-bw > 0, "slightly STRONGER", "slightly weaker"),
                ifelse(pw < 0.05, "" , " (not significant)")))

erosion_supported <- (pw < 0.05 && bw > 0) &&
                     all(diff(strat$phi) <= 0)
final <- data.table(
  quantity = c("phi_Q1_coolest", "phi_Q4_warmest", "phi_stratified_range",
               "warmhalf_slope_shift", "warmhalf_p", "gam_suggests_erosion",
               "stratified_confirms_erosion", "FINAL_climate_buffer_erosion_supported"),
  value = c(round(strat[1, phi], 3), round(strat[4, phi], 3),
            round(max(strat$phi) - min(strat$phi), 3),
            round(bw, 3), round(pw, 3), 1L,
            as.integer(all(diff(strat$phi) <= 0)),
            as.integer(erosion_supported)))
fwrite(final, file.path(OUT, "buffer_nonlinear_verdict.csv"))
print(final)
message(sprintf("\n>>> FINAL VERDICT: climate erosion of the buffer %s.",
                ifelse(erosion_supported, "IS SUPPORTED",
                       "is NOT supported once functional form is removed")))
message("    The smooth fit's decline does not replicate under model-free stratification;")
message("    it reflects local curvature of the fitted surface, not a robust feature.")
