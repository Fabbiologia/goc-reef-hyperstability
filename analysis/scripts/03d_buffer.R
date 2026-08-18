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
