# 03h_buffer_winter.R
# -----------------------------------------------------------
# Does the buffer answer to the winter rather than the summer?
#
# 03e and 03f test whether the buffer Phi = 1 - dlog(P)/dlog(B) changes
# with thermal exposure, using the warm season anomaly (May to October)
# of the reef's own OISST cell. That is the season of the heatwaves. But
# the reef's winter canopy (Sargassum and the other cool water macroalgae)
# grows from December to March and is suppressed by warm winters, and the
# canopy is the habitat in which the small, fast renewing fish that carry
# the buffer feed and shelter. If the buffer runs on that canopy, it should
# answer to the winter preceding the survey rather than to the summer.
#
# This step repeats the buffer tests with winter exposure defined as the
# mean daily anomaly over December (t-1) to March (t) in the same cell,
# against the same 1991 to 2020 climatology, for reef year t. Surveys fall
# mostly in July to October, so both the winter and the warm season
# precede them. The model, the clustering, the wild cluster bootstrap,
# the model-free stratified arbiter and the levels check are those of
# 03e, 03f and 03g, unchanged; only the exposure is swapped. Both
# exposures are fitted jointly as well, and the quartile estimates for
# both are written so that Figure 3d can show them side by side.
#
# Reads:  data/buffer_reef_year.csv           (03e)
#         ../data/env/sst_reef_daily.rds      (00b)
#         ../data/env/sst_reef_year.csv       (00b, reef -> cell)
# Writes: data/buffer_phi_by_season.csv       Phi by exposure quartile, both seasons
#         data/buffer_winter_models.csv       linear interaction, every specification
#         data/buffer_winter_summary.csv      in-text numbers and the verdict
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(fixest) })
setwd("..")
DATA <- "data"; OUT <- "data"; ENV <- "../data/env"
set.seed(20260823)

d <- fread(file.path(DATA, "buffer_reef_year.csv"))
d[, Reef := factor(Reef)][, YearF := factor(Year)]

# --- winter exposure per reef year ---------------------------------------
daily <- readRDS(file.path(ENV, "sst_reef_daily.rds"))
cells <- unique(fread(file.path(ENV, "sst_reef_year.csv"))[, .(Reef, cell)])
# winter of reef year t = Dec(t-1) .. Mar(t); Nov(t-1) .. Apr(t) as a check
daily[, wy := fifelse(mo >= 11L, yr + 1L, yr)]
win <- daily[mo %in% c(12L, 1L, 2L, 3L), .(C_win = mean(anom)), by = .(cell, Year = wy)]
win6 <- daily[mo %in% c(11L, 12L, 1L, 2L, 3L, 4L), .(C_win6 = mean(anom)), by = .(cell, Year = wy)]
win <- merge(win, win6, by = c("cell", "Year"))
win <- merge(cells, win, by = "cell", allow.cartesian = TRUE)[, .(Reef, Year, C_win, C_win6)]
d <- merge(d, win, by = c("Reef", "Year"))
stopifnot(!anyNA(d$C_win))
d[, `:=`(Cw = C_warm - mean(C_warm), Cn = C_win - mean(C_win), Cn6 = C_win6 - mean(C_win6))]

message(sprintf("Panel: %d reef years, %d reefs. Winter anomaly %.2f to %.2f C (warm season %.2f to %.2f)",
                nrow(d), uniqueN(d$Reef), min(d$C_win), max(d$C_win), min(d$C_warm), max(d$C_warm)))
r_seasons <- cor(d$C_warm, d$C_win)
sd_within <- function(x) d[, .(s = sd(get(x))), by = Year][, mean(s, na.rm = TRUE)]
message(sprintf("corr(warm season, winter) across reef years = %.2f; within-year between-reef sd: winter %.3f C, warm season %.3f C",
                r_seasons, sd_within("C_win"), sd_within("C_warm")))

# --- the linear interaction, every specification ---------------------------
grab <- function(fit, d, label, term, cl = "twoway") {
  V <- if (cl == "twoway") vcov(fit, cluster = ~ Reef + YearF)
       else if (cl == "reef") vcov(fit, cluster = ~ Reef)
       else vcov(fit, cluster = ~ YearF)
  cf <- coef(fit); se <- sqrt(diag(V))
  b2 <- cf[[term]]; s2 <- se[[term]]
  data.table(spec = label, exposure = term, cluster = cl,
             b1 = cf[["lB_c"]], b2 = b2, se_b2 = s2,
             t = b2 / s2, p = 2 * pnorm(-abs(b2 / s2)),
             lo = b2 - 1.96 * s2, hi = b2 + 1.96 * s2,
             phi_at_mean = 1 - cf[["lB_c"]], dphi_dC = -b2,
             n = nobs(fit), reefs = uniqueN(d$Reef))
}
m_w  <- feols(lP ~ lB_c * Cw  | Reef + YearF, data = d)
m_n  <- feols(lP ~ lB_c * Cn  | Reef + YearF, data = d)
m_n6 <- feols(lP ~ lB_c * Cn6 | Reef + YearF, data = d)
m_j  <- feols(lP ~ lB_c * Cw + lB_c * Cn | Reef + YearF, data = d)
models <- rbindlist(list(
  grab(m_w,  d, "Warm season (May to Oct), as 03e",           "lB_c:Cw"),
  grab(m_n,  d, "Winter (Dec to Mar)",                        "lB_c:Cn"),
  grab(m_n,  d, "Winter, clustered on year only",             "lB_c:Cn", "year"),
  grab(m_n,  d, "Winter, clustered on reef only",             "lB_c:Cn", "reef"),
  grab(m_n6, d, "Winter, Nov to Apr window",                  "lB_c:Cn6"),
  grab(m_j,  d, "Both seasons jointly: warm season term",     "lB_c:Cw"),
  grab(m_j,  d, "Both seasons jointly: winter term",          "lB_c:Cn")))
fwrite(models, file.path(OUT, "buffer_winter_models.csv"))
message("\nLinear interaction (positive b2 = buffer erodes as exposure rises):")
print(models[, .(spec, cluster, b2 = round(b2, 4), se = round(se_b2, 4), p = round(p, 3),
                 dPhi_per_C = round(dphi_dC, 3))])

# --- wild cluster bootstrap on year for the winter interaction -------------
wild_boot <- function(d, xvar, B = 999) {
  f0 <- feols(as.formula(paste0("lP ~ lB_c + ", xvar, " | Reef + YearF")), data = d)
  fit0 <- fitted(f0); r0 <- resid(f0)
  term <- paste0("lB_c:", xvar)
  tstat <- function(m) coef(m)[[term]] / sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[[term]]
  t_obs <- tstat(feols(as.formula(paste0("lP ~ lB_c * ", xvar, " | Reef + YearF")), data = d))
  yrs <- levels(d$YearF); ts <- numeric(B)
  for (b in seq_len(B)) {
    w <- setNames(sample(c(-1, 1), length(yrs), replace = TRUE), yrs)
    db <- copy(d)[, lP := fit0 + r0 * w[as.character(YearF)]]
    ts[b] <- tstat(feols(as.formula(paste0("lP ~ lB_c * ", xvar, " | Reef + YearF")), data = db))
  }
  list(t_obs = t_obs, p = mean(abs(ts) >= abs(t_obs)))
}
wb <- wild_boot(d, "Cn")
message(sprintf("Wild cluster bootstrap on year, winter interaction: t = %.3f, p = %.3f", wb$t_obs, wb$p))

# --- model-free arbiter: Phi within quartiles of each exposure ------------
strat_by <- function(xvar, label) {
  qs <- quantile(d[[xvar]], c(0, .25, .5, .75, 1))
  d[, Q := cut(get(xvar), qs, include.lowest = TRUE,
               labels = c("Q1 coolest", "Q2", "Q3", "Q4 warmest"))]
  out <- rbindlist(lapply(levels(d$Q), function(q) {
    dd <- d[Q == q]
    m  <- feols(lP ~ lB_c | Reef + YearF, data = dd)
    b  <- coef(m)[["lB_c"]]
    se <- sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[["lB_c"]]
    data.table(season = label, stratum = q, mean_anom = mean(dd[[xvar]]),
               phi = 1 - b, lo = 1 - b - 1.96 * se, hi = 1 - b + 1.96 * se,
               n = nrow(dd), reefs = uniqueN(dd$Reef))
  }))
  d[, Q := NULL]
  out
}
strat <- rbind(strat_by("C_warm", "Warm season (May to Oct)"),
               strat_by("C_win",  "Winter (Dec to Mar)"))
fwrite(strat, file.path(OUT, "buffer_phi_by_season.csv"))
message("\nPhi by exposure quartile (model free):")
print(strat[, .(season, stratum, mean_anom = round(mean_anom, 2), phi = round(phi, 3),
                lo = round(lo, 3), hi = round(hi, 3), n)])

# warm half against cool half of the winter record
d[, warmwin := as.integer(C_win > median(C_win))]
mw  <- feols(lP ~ lB_c * warmwin | Reef + YearF, data = d)
bw  <- coef(mw)[["lB_c:warmwin"]]
sew <- sqrt(diag(vcov(mw, cluster = ~ Reef + YearF)))[["lB_c:warmwin"]]
pw  <- 2 * pnorm(-abs(bw / sew))
message(sprintf("Warm winters vs cool winters: slope shift %+.3f (SE %.3f), p = %.3f -> Phi %+.3f when warm",
                bw, sew, pw, -bw))

# --- levels check, the 03g rule: does winter exposure move the levels -----
#     of P and B by plausible amounts? If a 1 C warmer winter "raises"
#     production or biomass by more than 15% within reef and year, the
#     exposure is picking up geography, and its interaction is not clean.
lev <- rbindlist(lapply(c("lP", "lB"), function(r) {
  m <- feols(as.formula(paste(r, "~ Cn | Reef + YearF")), data = d)
  data.table(response = r, est = coef(m)[["Cn"]],
             se = sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[["Cn"]])
}))
lev_max <- lev[, max(abs(est))]
message(sprintf("\nLevels on winter exposure: log P %+.3f (SE %.3f), log B %+.3f (SE %.3f) per C",
                lev[response == "lP", est], lev[response == "lP", se],
                lev[response == "lB", est], lev[response == "lB", se]))
identification <- ifelse(lev_max > 0.15,
  "levels implausible -> identification NOT clean, interaction reported with that caveat",
  "levels behave -> interaction may be interpreted")
message("  -> ", identification)

# --- verdict, by the same rule as 03f ------------------------------------
sn <- strat[season == "Winter (Dec to Mar)"]
b2n <- models[spec == "Winter (Dec to Mar)", b2]
p_n <- models[spec == "Winter (Dec to Mar)", p]
erosion_supported <- (p_n < 0.05 && b2n > 0) && (pw < 0.05 && bw > 0) && all(diff(sn$phi) <= 0)
message(sprintf("\nVERDICT winter erosion of the buffer supported: %s", erosion_supported))

summ <- data.table(quantity = c(
  "corr_warm_winter_reef_years", "sd_within_year_winter_C", "sd_within_year_warm_C",
  "winter_range_C_lo", "winter_range_C_hi",
  "b2_winter", "b2_winter_se", "b2_winter_p", "b2_winter_lo", "b2_winter_hi", "b2_winter_wild_p",
  "b2_winter_joint", "b2_winter_joint_p", "b2_warm_joint", "b2_warm_joint_p",
  "b2_winter_nov_apr", "b2_winter_nov_apr_p",
  "phi_winter_Q1", "phi_winter_Q2", "phi_winter_Q3", "phi_winter_Q4",
  "phi_warm_Q1", "phi_warm_Q2", "phi_warm_Q3", "phi_warm_Q4",
  "winter_warmhalf_slope_shift", "winter_warmhalf_p",
  "levels_lP_per_C_winter", "levels_lB_per_C_winter", "identification",
  "dPhi_winter_CI_across_range_lo", "dPhi_winter_CI_across_range_hi",
  "winter_erosion_supported"),
  value = c(
    round(r_seasons, 2), round(sd_within("C_win"), 3), round(sd_within("C_warm"), 3),
    round(min(d$C_win), 2), round(max(d$C_win), 2),
    round(b2n, 4), round(models[spec == "Winter (Dec to Mar)", se_b2], 4), round(p_n, 3),
    round(models[spec == "Winter (Dec to Mar)", lo], 4), round(models[spec == "Winter (Dec to Mar)", hi], 4),
    round(wb$p, 3),
    round(models[spec == "Both seasons jointly: winter term", b2], 4),
    round(models[spec == "Both seasons jointly: winter term", p], 3),
    round(models[spec == "Both seasons jointly: warm season term", b2], 4),
    round(models[spec == "Both seasons jointly: warm season term", p], 3),
    round(models[spec == "Winter, Nov to Apr window", b2], 4),
    round(models[spec == "Winter, Nov to Apr window", p], 3),
    round(sn$phi, 3),
    round(strat[season != "Winter (Dec to Mar)", phi], 3),
    round(bw, 3), round(pw, 3),
    round(lev[response == "lP", est], 3), round(lev[response == "lB", est], 3), identification,
    # change in Phi across the observed winter range implied by the CI on b2
    round(-models[spec == "Winter (Dec to Mar)", hi] * diff(range(d$C_win)), 3),
    round(-models[spec == "Winter (Dec to Mar)", lo] * diff(range(d$C_win)), 3),
    as.character(erosion_supported)))
fwrite(summ, file.path(OUT, "buffer_winter_summary.csv"))
message("\nStep 03h done. Winter buffer tables written to manuscript/data/.")
