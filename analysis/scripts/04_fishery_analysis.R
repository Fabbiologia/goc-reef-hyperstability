# 04_fishery_analysis.R
# -----------------------------------------------------------
# Two-mode lagged-climate model for the ARTISANAL ROCKY REEF fishery.
#
# Scope. The industrial fleet was removed upstream in 01 (we keep only
# TIPO AVISO == "MENORES"), and this step keeps only species whose genus
# is recorded on the LTEM rocky reef transects. Small pelagics, large
# pelagics, sharks and squid are therefore no longer modelled: they are
# different fisheries on different grounds and they confounded the reef
# signal. The reef aggregate below is the fishery the manuscript is
# about.
#
#   * CPUE = annual reef landings / annual total artisanal trips
#   * Single-lag scan k = 0..8 yr (impulse response)
#   * Two-mode model: ws_short (lag 0-3) + ws_long (lag 4-8)
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })
setwd("..")
DATA <- "data"; OUT <- "data"

art  <- fread(file.path(DATA, "artisanal_bcs_annual.csv"))
yt   <- fread(file.path(DATA, "artisanal_bcs_yearly_totals.csv"))
warm <- fread(file.path(DATA, "warm_season_anomaly_annual.csv"))

# Effort: reef_folios, the distinct landing receipts (trips) that landed at
# least one reef species. Every receipt at the same offices (folios) is the
# wrong denominator because the squid fishery moves it, see 01. See the note
# below on why boat-days cannot be used and why 2008 is a break point.
effort <- yt[, .(year, total_folios = reef_folios, all_folios = folios)]

ws <- setNames(warm$ws_anom, warm$year)

add_lags <- function(d) {
  for (k in 0:8) d[, paste0("ws_lag", k) := ws[as.character(year - k)]]
  d[, ws_short := rowMeans(.SD), .SDcols = paste0("ws_lag", 0:3)]
  d[, ws_long  := rowMeans(.SD), .SDcols = paste0("ws_lag", 4:8)]
  d[]
}

two_mode <- function(d, by) {
  d[, {
    fit <- lm(l_cpue ~ year + ws_short + ws_long, data = .SD)
    co  <- summary(fit)$coefficients
    .(n = .N, r2 = summary(fit)$r.squared,
      beta_year  = co["year", "Estimate"],     p_year  = co["year", "Pr(>|t|)"],
      beta_short = co["ws_short", "Estimate"], p_short = co["ws_short", "Pr(>|t|)"],
      beta_long  = co["ws_long", "Estimate"],  p_long  = co["ws_long", "Pr(>|t|)"])
  }, by = by]
}

lag_scan <- function(d, by) {
  out <- list()
  for (g in unique(d[[by]])) {
    s <- d[get(by) == g]
    for (k in 0:8) {
      fit <- lm(as.formula(paste0("l_cpue ~ year + ws_lag", k)), data = s)
      co  <- summary(fit)$coefficients
      out[[length(out) + 1]] <- data.table(
        group = g, lag = k, n = nrow(s),
        beta = co[paste0("ws_lag", k), "Estimate"],
        p    = co[paste0("ws_lag", k), "Pr(>|t|)"],
        r2   = summary(fit)$r.squared)
    }
  }
  setnames(rbindlist(out), "group", by)
}

# -----------------------------------------------------------
# 1. Reef aggregate CPUE
# -----------------------------------------------------------
reef_yr <- art[is_reef == TRUE,
               .(landings_t = sum(landings_t, na.rm = TRUE)), by = year]
reef_yr <- merge(reef_yr, effort, by = "year")
reef_yr[, CPUE   := landings_t / pmax(total_folios, 1)]
reef_yr[, l_cpue := log(pmax(CPUE, 1e-4))]
reef_yr <- add_lags(merge(reef_yr, warm, by = "year"))
reef_yr[, group := "Rocky reef fish"]

reef_two <- two_mode(reef_yr, "group")
fwrite(reef_two, file.path(OUT, "reef_aggregate_two_mode.csv"))
fwrite(lag_scan(reef_yr, "group"), file.path(OUT, "lag_scan_reef_aggregate.csv"))
message("Reef aggregate two-mode fit:"); print(reef_two)

# NOTE ON EFFORT. We previously carried a boat-days denominator
# (NUMERO EMBARCACIONES x DIAS EFECTIVOS) as a sensitivity check. It is
# unusable: those two fields contain impossible values in the later years,
# giving annual boat-day totals up to 3e14, so the check has been removed
# rather than reported. Distinct landing receipts (trips) is the only
# workable effort measure in these files, and it carries its own problem:
# the receipt count rises by about a quarter in 2008 with no matching change
# in the fishery (the exact figure is written by 07 to decoupling_summary.csv),
# so effort is NOT comparable across that year. Any statement about effort or
# catch per trip in the manuscript is therefore restricted to 2008 onward.

# -----------------------------------------------------------
# 2. Per-species reef targets
# -----------------------------------------------------------
# The named reef targets, ranked by landed weight. Groups are built from
# the scientific name in 01 (see REEF_GROUPS), not the common-name field.
top <- art[is_reef == TRUE & !is.na(reef_group) & reef_group != "",
           .(landings_t = sum(landings_t, na.rm = TRUE)),
           by = reef_group][order(-landings_t)]
fwrite(top, file.path(OUT, "reef_species_landed_ranking.csv"))

# keep groups present in most years so the 8-year lag model is identifiable
yrs   <- art[, uniqueN(year)]
cover <- art[is_reef == TRUE & !is.na(reef_group) & reef_group != "", .(ny = uniqueN(year)), by = reef_group]
keep  <- cover[ny >= 0.8 * yrs, reef_group]
reef_sp <- intersect(top$reef_group, keep)
message(sprintf("Modelling %d reef species groups: %s",
                length(reef_sp), paste(reef_sp, collapse = ", ")))

sp_yr <- art[is_reef == TRUE & reef_group %in% reef_sp,
             .(landings_t = sum(landings_t, na.rm = TRUE)),
             by = .(year, species = reef_group)]
sp_yr <- merge(sp_yr, effort, by = "year")
sp_yr[, CPUE   := landings_t / pmax(total_folios, 1)]
sp_yr[, l_cpue := log(pmax(CPUE, 1e-4))]
sp_yr <- add_lags(merge(sp_yr, warm, by = "year"))

sp_two <- two_mode(sp_yr, "species")
fwrite(sp_two[order(beta_long)], file.path(OUT, "reef_species_two_mode.csv"))
fwrite(lag_scan(sp_yr, "species"), file.path(OUT, "lag_scan_per_species.csv"))
print(sp_two[order(beta_long)])

message("\nStep 04 done.  Reef two-mode + impulse-response tables in manuscript/data/.")
