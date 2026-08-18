# 03d_buffer.R
# -----------------------------------------------------------
# The closing analysis: the reef has been living off its
# buffer, and is now living off its savings.
#
# Two products, both on the balanced 26-reef panel and both
# from the transect table written by 03b:
#
#   (a) the buffer over time: standing biomass (the savings)
#       and biomass production (the flow that pays the catch),
#       each as a reef-adjusted index set to 100 at its 1998
#       to 2004 mean. The gap between the two lines is the
#       compensation being spent.
#
#   (b) the rate of degradation under each pressure. Reefs are
#       classed by whether extraction is actually excluded
#       (Cabo Pulmo, the region's one long-enforced no-take
#       reserve) and years by the heatwave era (2014 onward,
#       when heatwave days per warm season tripled), giving
#       three regimes: fishing only (fished reefs before
#       2014), climate only (enforced reefs from 2014), and
#       both (fished reefs from 2014). Within each regime the
#       annual rate of change of biomass and of production is
#       fitted with reef fixed effects and reef-clustered
#       errors, exactly as every other trend in this paper.
#       The era split is a year-level classification, so the
#       regime contrast is descriptive in the same sense as
#       the heatwave-era comparisons already in the text.
#
# Requires: data/productivity_transects.csv        (from 03b)
#           data/productivity_fig2cd_residuals.csv (from 03b)
#           ../data/ltem.parquet                   (protection fields)
# Writes:   data/buffer_timeseries.csv
#           data/buffer_regime_rates.csv
#           data/buffer_summary.csv
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(arrow); library(data.table) })
setwd("..")
DATA <- "data"; OUT <- "data"

tr <- fread(file.path(DATA, "productivity_transects.csv"))

# --- protection per reef (08's rule: enforced = Cabo Pulmo no-take) -----
ltem <- as.data.table(read_parquet("../data/ltem.parquet"))
modal <- function(x) names(sort(table(x), decreasing = TRUE))[1]
prot <- ltem[Label == "PEC" & Reef %in% unique(tr$Reef),
             .(mpa = modal(MPA), level = modal(Protection_level)), by = Reef]
prot[, enforced := mpa == "Cabo Pulmo" & level == "Prohibited"]
message(sprintf("Panel protection: %d enforced no-take reefs, %d fished reefs",
                prot[enforced == TRUE, .N], prot[enforced == FALSE, .N]))
tr <- merge(tr, prot[, .(Reef, enforced)], by = "Reef")
tr[, era := fifelse(Year >= 2014, "heatwave era", "before 2014")]

# --- (a) buffer over time: % of the 1998-2004 baseline ------------------
res <- fread(file.path(DATA, "productivity_fig2cd_residuals.csv"))[panel == "c"]
res[, base := mean(mean[Year %in% 1998:2004]), by = group]
res[, `:=`(idx    = 100 * exp(mean - base),
           idx_lo = 100 * exp(mean - ci95 - base),
           idx_hi = 100 * exp(mean + ci95 - base))]
fwrite(res[, .(Year, group, idx, idx_lo, idx_hi)],
       file.path(OUT, "buffer_timeseries.csv"))

# --- (b) regime rates ---------------------------------------------------
clust_beta <- function(fit, cluster, term) {
  cf <- coef(fit); k <- names(cf)[!is.na(cf)]
  X  <- model.matrix(fit)[, k, drop = FALSE]
  u  <- residuals(fit); g <- factor(cluster[as.integer(names(u))])
  br <- chol2inv(chol(crossprod(X))); mt <- matrix(0, ncol(X), ncol(X))
  for (l in levels(g)) { i <- which(g == l)
    xu <- crossprod(X[i, , drop = FALSE], u[i]); mt <- mt + tcrossprod(xu) }
  m <- nlevels(g)
  V <- br %*% mt %*% br * (m / (m - 1)) * ((nrow(X) - 1) / (nrow(X) - ncol(X)))
  se <- sqrt(diag(V)); names(se) <- k
  list(beta = cf[term], se = se[term],
       p = 2 * pt(-abs(cf[term] / se[term]), df = m - 1))
}

regimes <- list(
  "Fishing only"  = quote(enforced == FALSE & era == "before 2014"),
  "Climate only"  = quote(enforced == TRUE  & era == "heatwave era"),
  "Both"          = quote(enforced == FALSE & era == "heatwave era"),
  "Neither (reference)" = quote(enforced == TRUE & era == "before 2014"))

rate_one <- function(d, value) {
  d <- copy(d)[!is.na(get(value))]
  d[, l_v := log(get(value) + 0.01)]
  if (uniqueN(d$Reef) > 1) {
    fit <- lm(l_v ~ Year + Reef, data = d)
  } else {
    fit <- lm(l_v ~ Year, data = d)
  }
  cb <- clust_beta(fit, d$Reef, "Year")
  data.table(pct_per_decade = 100 * (exp(cb$beta * 10) - 1),
             lo = 100 * (exp((cb$beta - 1.96 * cb$se) * 10) - 1),
             hi = 100 * (exp((cb$beta + 1.96 * cb$se) * 10) - 1),
             p  = cb$p, n_reefs = uniqueN(d$Reef), n_transects = nrow(d))
}

rates <- rbindlist(lapply(names(regimes), function(nm) {
  d <- tr[eval(regimes[[nm]])]
  rbind(data.table(regime = nm, response = "Biomass (the savings)",    rate_one(d, "b100")),
        data.table(regime = nm, response = "Production (the buffer)",  rate_one(d, "p100")))
}))
fwrite(rates, file.path(OUT, "buffer_regime_rates.csv"))
message("\nRegime rates (% per decade, reef FE, reef-clustered):")
print(rates[, .(regime, response, pct_per_decade = round(pct_per_decade, 1),
                lo = round(lo, 1), hi = round(hi, 1), p = signif(p, 2), n_reefs)])

# --- (c) the projection: how long the buffer holds -----------------------
# Two log-linear extrapolations of measured rates, anchored at the fitted
# 2025 level of the production index. Status quo continues the full-record
# production trend (everything as it has been, fishing and climate as
# experienced). The second scenario adds the measured thermal sensitivity of
# production on fished reefs, applied to the basin's ongoing warming trend.
# These are extrapolations of measured rates for timing, not population
# forecasts, and are labelled as such in the caption and Methods.

# full-record production trend on the panel, with clustered SE
d_all <- copy(tr)[!is.na(p100)]
d_all[, l_v := log(p100 + 0.01)]
fit_all <- lm(l_v ~ Year + Reef, data = d_all)
cb_all  <- clust_beta(fit_all, d_all$Reef, "Year")
r1    <- cb_all$beta                       # log rate per year
r1_lo <- cb_all$beta - 1.96 * cb_all$se
r1_hi <- cb_all$beta + 1.96 * cb_all$se

# thermal sensitivity: mean per-reef production response on fished reefs
rts <- fread(file.path(DATA, "reef_thermal_slopes.csv"))
rts <- merge(rts, prot[, .(Reef, enforced)], by = "Reef")
sens_pct_perC <- rts[enforced == FALSE & !is.na(slope), mean(slope)]

# basin warming trend, C per decade, from the same monthly series as Fig 1
mon <- fread(file.path("../data/env", "sst_gulf_monthly_1981_2026.csv"))
mon[, `:=`(yr = as.integer(substr(month, 1, 4)), mo = as.integer(substr(month, 6, 7)))]
cl  <- mon[yr %between% c(1991, 2020), .(clim = mean(mean_sst)), by = mo]
mon <- merge(mon, cl, by = "mo")
mon[, anom := mean_sst - clim]
warm_dec <- coef(lm(anom ~ I(yr + (mo - 0.5) / 12), data = mon))[2] * 10

# climate penalty as a log rate per year
pen <- log(1 + (sens_pct_perC / 100) * warm_dec) / 10
r2 <- r1 + pen; r2_lo <- r1_lo + pen; r2_hi <- r1_hi + pen
message(sprintf("Projection rates: status quo %.1f %%/decade; + warming %.1f %%/decade (sens %.1f %%/C x %.2f C/decade)",
                100 * (exp(10 * r1) - 1), 100 * (exp(10 * r2) - 1), sens_pct_perC, warm_dec))

# anchor: fitted index level at 2025 from the annual index series
idx <- res[group == "Biomass production", .(Year, idx)]
af  <- lm(log(idx) ~ Year, data = idx)
idx0 <- exp(predict(af, data.frame(Year = 2025)))

proj_years <- 2025:2060
proj <- rbindlist(list(
  data.table(scenario = "Status quo", Year = proj_years,
             idx    = idx0 * exp(r1    * (proj_years - 2025)),
             idx_lo = idx0 * exp(r1_lo * (proj_years - 2025)),
             idx_hi = idx0 * exp(r1_hi * (proj_years - 2025))),
  data.table(scenario = "Fishing plus continued warming", Year = proj_years,
             idx    = idx0 * exp(r2    * (proj_years - 2025)),
             idx_lo = idx0 * exp(r2_lo * (proj_years - 2025)),
             idx_hi = idx0 * exp(r2_hi * (proj_years - 2025)))))
fwrite(proj, file.path(OUT, "buffer_projection.csv"))

cross_year <- function(r, level) 2025 + log(level / idx0) / r
crossings <- data.table(
  quantity = c("proj_rate_statusquo_pct_decade", "proj_rate_climate_pct_decade",
               "proj_sensitivity_pct_perC", "proj_warming_C_per_decade",
               "proj_anchor_idx_2025",
               "half_baseline_year_statusquo", "half_baseline_year_climate",
               "quarter_baseline_year_statusquo", "quarter_baseline_year_climate"),
  value = c(round(100 * (exp(10 * r1) - 1), 1), round(100 * (exp(10 * r2) - 1), 1),
            round(sens_pct_perC, 1), round(warm_dec, 2), round(idx0, 1),
            round(cross_year(r1, 50)), round(cross_year(r2, 50)),
            round(cross_year(r1, 25)), round(cross_year(r2, 25))))
fwrite(crossings, file.path(OUT, "buffer_projection_summary.csv"))
print(crossings)

gr <- function(rg, rs, col) rates[regime == rg & response == rs][[col]]
summ <- data.table(quantity = c(
  "fishingonly_biomass_pct_decade", "fishingonly_production_pct_decade",
  "climateonly_biomass_pct_decade", "climateonly_production_pct_decade",
  "both_biomass_pct_decade", "both_production_pct_decade",
  "both_production_p", "fishingonly_production_p", "climateonly_production_p"),
  value = c(round(gr("Fishing only", "Biomass (the savings)",   "pct_per_decade"), 1),
            round(gr("Fishing only", "Production (the buffer)", "pct_per_decade"), 1),
            round(gr("Climate only", "Biomass (the savings)",   "pct_per_decade"), 1),
            round(gr("Climate only", "Production (the buffer)", "pct_per_decade"), 1),
            round(gr("Both", "Biomass (the savings)",   "pct_per_decade"), 1),
            round(gr("Both", "Production (the buffer)", "pct_per_decade"), 1),
            signif(gr("Both", "Production (the buffer)", "p"), 2),
            signif(gr("Fishing only", "Production (the buffer)", "p"), 2),
            signif(gr("Climate only", "Production (the buffer)", "p"), 2)))
fwrite(summ, file.path(OUT, "buffer_summary.csv"))
print(summ)
message("\nStep 03d done.  Buffer tables written to manuscript/data/.")
