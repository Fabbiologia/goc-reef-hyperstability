# 02_climate_analysis.R
# -----------------------------------------------------------
# Climate diagnostics from the monthly SST anomaly series:
#   * Linear secular trend in monthly anomaly (whole Gulf)
#   * Per-1°-latitude decadal trend
#   * Warm-season anomaly per year (May–October) used downstream
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })
setwd("..")
DATA <- "data"
OUT  <- "data"

mon <- fread(file.path(DATA, "sst_gulf_monthly_1981_2026.csv"))
mon[, dt    := as.Date(paste0(month, "-01"))]
mon[, year  := year(dt)]
mon[, mo    := month(dt)]

# 1991–2020 monthly climatology
clim <- mon[year %between% c(1991, 2020),
            .(clim = mean(mean_sst)),
            by = mo]
mon <- merge(mon, clim, by = "mo")
mon[, anom := mean_sst - clim]

# Linear secular trend on monthly anomaly
mon[, yfrac := year + (mo - 0.5)/12]
fit_all <- lm(anom ~ yfrac, data = mon)
trend_dec <- coef(fit_all)["yfrac"] * 10
mon[, anom_detrended := anom - predict(fit_all)]
message(sprintf("Whole-Gulf trend: %+.3f °C / decade", trend_dec))

# Warm-season (May–October) anomaly per year — used by downstream scripts
warm <- mon[mo %between% c(5, 10),
            .(ws_anom = mean(anom)),
            by = year][order(year)]
fwrite(warm, file.path(OUT, "warm_season_anomaly_annual.csv"))

# Per-latitude trend (using by-lat file)
lat_path <- file.path("..","data","env","sst_gulf_by_lat_degree_daily.csv")
if (file.exists(lat_path)) {
  d <- fread(lat_path)
  d[, doy  := yday(date)]
  d[, year := year(date)]
  clim_lat <- d[year %between% c(1991, 2020),
                .(clim = mean(sst_mean)), by = .(lat_degree, doy)]
  d <- merge(d, clim_lat, by = c("lat_degree","doy"))
  d[, anom := sst_mean - clim]
  d[, yfrac := year + (yday(date) - 0.5)/365]
  trend_lat <- d[, .(trend_dec = coef(lm(anom ~ yfrac))[2] * 10),
                  by = lat_degree][order(lat_degree)]
  fwrite(trend_lat, file.path(OUT, "sst_trend_by_latitude.csv"))
  message("Per-latitude trends written.")
}

message("\nStep 02 done.")
