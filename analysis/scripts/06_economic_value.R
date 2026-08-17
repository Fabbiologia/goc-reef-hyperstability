# 06_economic_value.R
# -----------------------------------------------------------
# Ex-vessel (first-sale "precio de playa", valor_pesos) economic
# exposure of the Gulf of California ARTISANAL reef fishery, Figure 3.
#
# What changed from the earlier version, and why it matters: the old
# five-state total summed every fleet, so the headline value included
# the industrial catch. It also classified "reef" partly by a
# common-name regex that swept in soft-bottom and estuarine families
# (Sciaenidae, Mugilidae, Gerreidae). Both are gone. Here:
#
#   (a) Figure 3a -- five-state ARTISANAL first-sale value across the
#       Gulf states (BC, BCS, Sonora, Sinaloa, Nayarit), Pacific
#       litoral, 2022-2025. Small-scale fleet only.
#   (b) Figure 3b -- La Paz + Loreto species-level first-sale value,
#       2022-2025, rocky reef species only, joined to each species'
#       warming-response class.
#
# Rocky reef membership is decided in 01 by genus against the LTEM reef
# transects, not by common name.
#
# valor_pesos is only populated consistently from the 2022 reporting
# reform onward, so all economic figures use 2022-2025.
#
# Requires: data/artisanal_5state_annual.csv   (from 01)
# Writes:   data/economic_5state_artisanal.csv
#           data/economic_lapaz_loreto_species.csv
#           data/economic_summary.csv          (in-text numbers)
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })
setwd("..")
DATA <- "data"; OUT <- "data"

goc <- fread(file.path(DATA, "artisanal_5state_annual.csv"))
ECON_YEARS <- 2022:2025

# -----------------------------------------------------------
# Currency. All reported values are US dollars, because a peso figure is
# not interpretable to most readers. Each year is converted at that year's
# own average rate rather than at a single blended rate, so the conversion
# does not smuggle in exchange-rate movement. Rates are the annual averages
# of daily noon buying rates published by the US Federal Reserve (H.10,
# series AEXMXUS).
FX <- data.table(year = 2022:2025,
                 mxn_per_usd = c(20.1208, 17.7334, 18.3062, 19.2052))
FX_REF <- FX[, mean(mxn_per_usd)]   # 18.84, used for the constant-price series
goc <- merge(goc, FX, by = "year", all.x = TRUE)
goc[is.na(mxn_per_usd), mxn_per_usd := FX_REF]
goc[, value_usd   := value_mxn   / mxn_per_usd]
# value_const is already expressed at fixed 2022-2025 prices, so it takes the
# 2022-2025 mean rate throughout.
goc[, value_const_usd := value_const / FX_REF]
message(sprintf("currency: per-year FX applied; 2022-2025 mean = %.2f MXN per USD", FX_REF))

# -----------------------------------------------------------
# (a) Five-state artisanal first-sale value  (Figure 3a)
# -----------------------------------------------------------
st <- goc[year %in% ECON_YEARS,
          .(value_usd  = sum(value_usd,  na.rm = TRUE),
            value_mxn  = sum(value_mxn,  na.rm = TRUE),
            landings_t = sum(landings_t, na.rm = TRUE)),
          by = .(year, state)]
fwrite(st[order(state, year)], file.path(OUT, "economic_5state_artisanal.csv"))

annual_tot <- st[, .(value_musd = sum(value_usd) / 1e6), by = year][order(year)]
state_mean <- st[, .(value_Myr = sum(value_usd) / length(ECON_YEARS) / 1e6),
                 by = state][order(-value_Myr)]
message("\nFive-state ARTISANAL value by year (million USD):"); print(annual_tot)
message("Mean per state (M USD/yr):"); print(state_mean)

# reef share of that artisanal total, for the manuscript
reef_share <- goc[year %in% ECON_YEARS,
                  .(value_usd = sum(value_usd, na.rm = TRUE)), by = is_reef]
reef_share[, pct := round(100 * value_usd / sum(value_usd), 1)]
message("Rocky reef share of the five-state artisanal value:"); print(reef_share)

# -----------------------------------------------------------
# (a2) The full record, 2000 to 2025, at constant prices  (Figure 3a)
# -----------------------------------------------------------
# The landing receipts carry a usable price in over 99 percent of rows in
# every year back to 2000, so the value series is not limited to the recent
# window. It cannot be read in nominal pesos across 25 years of inflation,
# so every year's catch is valued at a fixed reference price per species
# (the 2022 to 2025 average). What remains is the change in how much is
# landed and of what, with all price movement removed.
eff <- fread(file.path(DATA, "artisanal_5state_effort.csv"))

ts <- goc[, .(all_const_M  = sum(value_const_usd, na.rm = TRUE) / 1e6,
              reef_const_M = sum(value_const_usd[is_reef == TRUE], na.rm = TRUE) / 1e6,
              reef_t       = sum(landings_t[is_reef == TRUE], na.rm = TRUE)),
          by = year][order(year)]
ts <- merge(ts, eff[, .(year, trips)], by = "year")
ts[, reef_kg_per_trip := 1000 * reef_t / trips]
ts[, reef_pct_of_all  := 100 * reef_const_M / all_const_M]
fwrite(ts, file.path(OUT, "economic_timeseries_constant_price.csv"))
message("\nConstant-price series (2022-2025 pesos), five Gulf states:")
print(ts[, .(year, all_const_M = round(all_const_M), reef_const_M = round(reef_const_M),
             trips, reef_kg_per_trip = round(reef_kg_per_trip))])

pk <- ts[year <= 2019][which.max(reef_const_M)]
rec <- ts[year %in% 2021:2025, .(v = mean(reef_const_M), cpue = mean(reef_kg_per_trip))]
add_ts <- function(k, v) summ_ts[[length(summ_ts) + 1]] <<- data.table(quantity = k, value = v)
summ_ts <- list()
add_ts("reef_value_2000_MUSD_const",      round(ts[year == 2000, reef_const_M]))
add_ts("reef_value_peak_year",           pk$year)
add_ts("reef_value_peak_MUSD_const",      round(pk$reef_const_M))
add_ts("reef_value_2021_2025_MUSD_const", round(rec$v))
add_ts("reef_value_pct_below_peak",      round(100 * (rec$v / pk$reef_const_M - 1), 1))
add_ts("trips_2000",                     ts[year == 2000, trips])
add_ts("trips_2025",                     ts[year == 2025, trips])
add_ts("reef_cpue_peak_kg_per_trip",     round(ts[year == pk$year, reef_kg_per_trip]))
add_ts("reef_cpue_2021_2025_kg_per_trip", round(rec$cpue))
add_ts("reef_cpue_pct_below_peak",       round(100 * (rec$cpue / ts[year == pk$year, reef_kg_per_trip] - 1), 1))

# -----------------------------------------------------------
# (b) La Paz + Loreto reef species value  (Figure 3b)
# -----------------------------------------------------------
ll <- goc[grepl("LA PAZ|LORETO", oficina) & is_reef == TRUE & year %in% ECON_YEARS]
ll <- ll[!is.na(reef_group) & reef_group != ""]
sp_mean <- ll[, .(value_Myr  = sum(value_usd,  na.rm = TRUE) / length(ECON_YEARS) / 1e6,
                  landings_t = sum(landings_t, na.rm = TRUE) / length(ECON_YEARS)),
              by = reef_group][order(-value_Myr)]

# Each group is classed by its measured TIME TREND on the monitored reefs
# (reef_group_trend.csv from 03), not by its response to warm years. The trend
# is identified at the reef level and survives clustering; the warming response
# is a year-level regressor and does not (see the warning in 03). Class drives
# the Figure-3 colour:
#   decline_sig   significant decline on the monitored reefs (p < 0.05)
#   no_sig_trend  no significant trend either way
#   increase_sig  significant increase
gt <- fread(file.path(DATA, "reef_group_trend.csv"))
# The Pacific red snapper is a commercial split of genus Lutjanus, so it has
# no genus-level reef row of its own and inherits the snapper trend.
gt <- rbind(gt, copy(gt[reef_group == "Snappers"])[, reef_group := "Pacific red snapper"])
reef <- merge(sp_mean, gt[, .(reef_group, trend_class, pct_per_decade, p_year)],
              by = "reef_group", all.x = TRUE)
setnames(reef, "trend_class", "warming_class")
reef[is.na(warming_class), warming_class := "not_monitored"]
reef <- reef[order(-value_Myr)]
fwrite(reef, file.path(OUT, "economic_lapaz_loreto_species.csv"))

# The denominator is the whole reef catch, including the groups whose warming
# response we could not resolve. Dropping them would inflate the share.
reef_total <- reef[, sum(value_Myr)]
declining  <- reef[warming_class == "decline_sig", sum(value_Myr)]
message("\nLa Paz + Loreto reef-species first-sale value (M pesos/yr, 2022-2025):")
print(reef[, .(reef_group, value_Myr = round(value_Myr, 1),
               pct_per_decade = round(pct_per_decade, 1), warming_class)])
message(sprintf("Reef-species total: %.0f M/yr | SIGNIFICANTLY declining share: %.0f%% (%.0f M/yr)",
                reef_total, 100 * declining / reef_total, declining))

# -----------------------------------------------------------
# In-text summary numbers -> economic_summary.csv (traceable)
# -----------------------------------------------------------
val <- function(sp) { v <- reef[reef_group == sp, value_Myr]; if (length(v)) round(v, 1) else NA_real_ }
summ_ts[[length(summ_ts)+1]] <- data.table(quantity='mxn_per_usd_mean_2022_2025', value=round(FX_REF,2))
summ <- rbindlist(summ_ts)
summ <- rbind(summ, data.table(
  quantity = c("five_state_artisanal_mean_MUSD_per_yr",
               "five_state_artisanal_reef_pct",
               "lapaz_loreto_reef_total_MUSD_per_yr",
               "lapaz_loreto_reef_declining_pct",
               "lapaz_loreto_reef_declining_MUSD_per_yr",
               "pacific_red_snapper_MUSD_per_yr", "leopard_grouper_MUSD_per_yr"),
  value = c(round(mean(annual_tot$value_musd)),
            reef_share[is_reef == TRUE, pct],
            round(reef_total, 1),
            round(100 * declining / reef_total, 1),
            round(declining, 1),
            val("Pacific red snapper"), val("Leopard grouper"))))
fwrite(summ, file.path(OUT, "economic_summary.csv"))
print(summ)
message("\nStep 06 done.  Economic tables + summary written to manuscript/data/.")
