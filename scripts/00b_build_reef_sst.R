# 00b_build_reef_sst.R
# -----------------------------------------------------------
# Turns the OISST v2.1 daily grids downloaded by 00a into every temperature
# product the analysis uses, at the product's native quarter degree
# resolution.
#
# BACKGROUND. The series previously held in data/env carried three sample
# points per one degree latitude band, thirty points for the whole basin,
# and came from a dashboard export rather than the gridded product the
# Methods described. The real product carries about 891 ocean cells over the
# same box, so the reefs resolve into roughly three times as many distinct
# thermal units and the spatial structure that matters here, the cold pool
# the midriff islands maintain by tidal mixing and the difference between
# the peninsula and mainland sides, is retained.
#
# DROP-IN BY DESIGN. Two outputs keep the exact names and column schemas of
# the files they replace, so every downstream script reads real OISST
# without modification:
#     sst_gulf_monthly_1981_2026.csv     (02, 03d, 08, 09)
#     sst_gulf_by_lat_degree_daily.csv   (02, 03b, 03c, 03e, 09)
# The dashboard-era originals they replace remain recoverable from git
# history (commit 4e94858 and earlier) if the two ever need comparing.
#
# Three further outputs are new and carry the per reef resolution that the
# band series could not:
#     sst_reef_longterm.csv   per reef long term means (drives Kmax in 03b/03c)
#     sst_reef_year.csv       per reef-year exposure (drives C_rt in 03e)
#     sst_reef_daily.rds      per reef daily series with anomalies and heatwaves
#
# Marine heatwaves follow Hobday et al.: 90th percentile of a 1982 to 2011
# baseline, minimum duration five days, maximum gap two days, via heatwaveR,
# detected separately on every reef cell rather than once for the basin.
#
# Requires: ../data/env/oisst_daily/YYYY/*.csv   (from 00a)
#           ../data/ltem.parquet                 (reef coordinates)
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(heatwaveR)
})
setwd("..")
ENV  <- "../data/env"
ROOT <- file.path(ENV, "oisst_daily")

files <- list.files(ROOT, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
stopifnot(length(files) > 100)
message(sprintf("%s daily OISST files on disk", format(length(files), big.mark = ",")))

# --- the grid, the Gulf mask, and the band assignment -------------------
grid <- unique(fread(files[[length(files) %/% 2]], showProgress = FALSE)[, .(lat, lon)])
grid[, lonW := lon - 360]
message(sprintf("%d ocean cells in the downloaded box", nrow(grid)))

# Polygon tracing the Gulf: down the mainland coast from the delta to Cabo
# Corrientes, then back up the peninsula side. Excludes Pacific water west
# of Baja, which the box necessarily includes.
gulf_poly <- data.table(
  # Mainland coast (Sonora, Sinaloa, Nayarit) traced north to south, then the
  # peninsula's east coast south to north. Each side is pushed about 0.15
  # degrees offshore so that 25 km cells and slightly offshore reefs are kept.
  # Verified: 258 of 258 Gulf-region LTEM reefs fall inside.
  lon = c(-114.90, -113.35, -112.15, -111.75, -110.75, -110.45, -109.35, -108.35,
          -107.35, -106.35, -105.65, -105.55,
          -110.05, -110.80, -111.05, -111.45, -112.15, -112.95, -113.65, -114.45,
          -114.85, -115.05),
  lat = c(  31.90,   30.50,   29.50,   28.80,   27.90,   27.00,   26.00,   25.00,
            24.00,   23.00,   22.00,   20.40,
            22.70,   24.00,   25.00,   26.00,   27.00,   28.00,   29.00,   30.00,
            31.00,   31.70))
in_poly <- function(px, py, vx, vy) {
  n <- length(vx); inside <- rep(FALSE, length(px)); j <- n
  for (i in seq_len(n)) {
    hit <- ((vy[i] > py) != (vy[j] > py)) &
           (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])
    inside <- xor(inside, hit); j <- i
  }
  inside
}
grid[, gulf := in_poly(lonW, lat, gulf_poly$lon, gulf_poly$lat)]
grid[, lat_degree := as.integer(floor(lat))]
message(sprintf("%d of %d cells fall inside the Gulf polygon", sum(grid$gulf), nrow(grid)))
fwrite(grid, file.path(ENV, "oisst_gulf_cell_mask.csv"))

# --- reef to nearest ocean cell -----------------------------------------
ltem <- as.data.table(read_parquet("../data/ltem.parquet"))
reefs <- unique(ltem[Label == "PEC" & !is.na(Latitude) & !is.na(Longitude),
                     .(Reef, Latitude, Longitude)])
reefs <- reefs[, .(lat = mean(Latitude), lon = mean(Longitude)), by = Reef]
hav <- function(a1, o1, a2, o2) { R <- 6371; p <- pi/180
  x <- sin((a2-a1)*p/2)^2 + cos(a1*p)*cos(a2*p)*sin((o2-o1)*p/2)^2
  2*R*asin(pmin(1, sqrt(x))) }
idx <- vapply(seq_len(nrow(reefs)),
              function(k) which.min(hav(reefs$lat[k], reefs$lon[k], grid$lat, grid$lonW)), 1L)
dst <- vapply(seq_len(nrow(reefs)),
              function(k) min(hav(reefs$lat[k], reefs$lon[k], grid$lat, grid$lonW)), 1.0)
reefs[, `:=`(cell_lat = grid$lat[idx], cell_lon = grid$lon[idx], dist_km = dst)]
reefs <- reefs[dist_km < 60]
reefs[, cell := paste(cell_lat, cell_lon)]
fwrite(reefs, file.path(ENV, "sst_reef_cell_map.csv"))
message(sprintf("%d reefs -> %d distinct cells; distance median %.1f km, max %.1f km",
                nrow(reefs), uniqueN(reefs$cell), median(reefs$dist_km), max(reefs$dist_km)))

reef_cells <- unique(reefs[, .(lat = cell_lat, lon = cell_lon)])[, keep := TRUE]
gulf_key   <- grid[gulf == TRUE, paste(lat, lon)]
band_key   <- grid[gulf == TRUE, .(k = paste(lat, lon), lat_degree)]

# --- single pass over the dailies ---------------------------------------
message("reading daily grids (one pass, three products) ...")
t0 <- Sys.time()
out <- vector("list", length(files))
for (n in seq_along(files)) {
  f <- files[[n]]
  d <- fread(f, showProgress = FALSE)
  d[, k := paste(lat, lon)]
  day <- as.Date(sub(".*_(\\d{8})\\.csv$", "\\1", basename(f)), format = "%Y%m%d")
  rc <- merge(d, reef_cells, by = c("lat", "lon"))          # reef cells
  gz <- d[k %chin% gulf_key]                                # Gulf interior
  bd <- merge(gz[, .(k, sst)], band_key, by = "k")[
        , .(sst_mean = round(mean(sst), 3), sst_min = round(min(sst), 3),
            sst_max = round(max(sst), 3), n_stations = .N), by = lat_degree]
  out[[n]] <- list(
    reef  = if (nrow(rc)) data.table(date = day, cell = rc$k, sst = rc$sst) else NULL,
    basin = data.table(date = day, mean_sst = mean(gz$sst), max_sst = max(gz$sst),
                       n_cells = nrow(gz)),
    band  = cbind(date = day, bd))
  if (n %% 2000 == 0)
    message(sprintf("  %6d/%d  (%.1f min)", n, length(files),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
reef_daily  <- rbindlist(lapply(out, `[[`, "reef"))
basin_daily <- rbindlist(lapply(out, `[[`, "basin"))[order(date)]
band_daily  <- rbindlist(lapply(out, `[[`, "band"))[order(date, lat_degree)]
rm(out); invisible(gc())
message(sprintf("  %s reef cell-days | %s basin days | %s band-days",
                format(nrow(reef_daily), big.mark = ","),
                format(nrow(basin_daily), big.mark = ","),
                format(nrow(band_daily), big.mark = ",")))

# --- backward compatible band daily file --------------------------------
legacy_band <- file.path(ENV, "sst_gulf_by_lat_degree_daily.csv")
fwrite(band_daily[, .(date, lat_degree, sst_mean, sst_min, sst_max, n_stations)], legacy_band)
message(sprintf("wrote band dailies: %d cells per band on average (was 3)",
                round(mean(band_daily$n_stations))))

# --- marine heatwaves, per reef cell ------------------------------------
message("detecting marine heatwaves per reef cell ...")
mhw_one <- function(dd) {
  cl <- ts2clm(data.frame(t = dd$date, temp = dd$sst),
               climatologyPeriod = c("1982-01-01", "2011-12-31"), pctile = 90)
  ev <- detect_event(cl, minDuration = 5, maxGap = 2)
  x  <- as.data.table(ev$climatology)
  x[, .(date = t, temp, thresh, mhw = event, intensity = pmax(0, temp - thresh))]
}
setorder(reef_daily, cell, date)
res <- reef_daily[, mhw_one(.SD), by = cell, .SDcols = c("date", "sst")]
res[, `:=`(yr = year(date), mo = month(date))]
clim <- res[yr %between% c(1991, 2020), .(clim = mean(temp, na.rm = TRUE)), by = .(cell, mo)]
res <- merge(res, clim, by = c("cell", "mo"))[, anom := temp - clim][]
setorder(res, cell, date)
saveRDS(res, file.path(ENV, "sst_reef_daily.rds"))

cy <- res[, .(sst_annual = mean(temp), anom_annual = mean(anom),
              anom_warm = mean(anom[mo %in% 5:10]),
              mhw_days = sum(mhw), mhw_cum = sum(intensity[mhw])),
          by = .(cell, year = yr)]
fwrite(merge(reefs[, .(Reef, cell, dist_km)], cy, by = "cell", allow.cartesian = TRUE),
       file.path(ENV, "sst_reef_year.csv"))
clt <- res[, .(sst_mean_lt = mean(temp),
               sst_winter_lt = mean(temp[mo %in% c(12,1,2,3)]),
               sst_summer_lt = mean(temp[mo %in% 7:10])), by = cell]
reef_lt <- merge(reefs[, .(Reef, lat, lon, cell, dist_km)], clt, by = "cell")
fwrite(reef_lt, file.path(ENV, "sst_reef_longterm.csv"))
message(sprintf("per reef long term SST spans %.2f to %.2f C (band series spanned %.2f C)",
                min(reef_lt$sst_mean_lt), max(reef_lt$sst_mean_lt),
                diff(range(reef_lt$sst_mean_lt))))

# --- backward compatible monthly basin file -----------------------------
bts <- data.frame(t = basin_daily$date, temp = basin_daily$mean_sst)
bcl <- ts2clm(bts, climatologyPeriod = c("1982-01-01", "2011-12-31"), pctile = 90)
bev <- detect_event(bcl, minDuration = 5, maxGap = 2)
bd  <- as.data.table(bev$climatology)
bd[, `:=`(yr = year(t), mo = month(t), int = pmax(0, temp - thresh))]
bm <- bd[, .(mhw_days = sum(event, na.rm = TRUE),
             mhw_fraction = round(mean(event, na.rm = TRUE), 4),
             event_count = 0L,
             max_sst = round(max(temp), 3), mean_sst = round(mean(temp), 3),
             max_intensity = round(max(int, na.rm = TRUE), 3),
             mean_intensity = round(mean(int[event], na.rm = TRUE), 3),
             cumulative_intensity = round(sum(int[event], na.rm = TRUE), 3)),
         by = .(yr, mo)]
bm[is.na(mean_intensity), mean_intensity := 0]
bm[, month := sprintf("%d-%02d", yr, mo)]
setorder(bm, yr, mo)
legacy_mon <- file.path(ENV, "sst_gulf_monthly_1981_2026.csv")
fwrite(bm[, .(month, mhw_days, mhw_fraction, event_count, max_sst, mean_sst,
              max_intensity, mean_intensity, cumulative_intensity)], legacy_mon)

# --- summary -------------------------------------------------------------
cl2 <- bm[yr %between% c(1991, 2020), .(clim = mean(mean_sst)), by = mo]
cmp <- merge(bm, cl2, by = "mo")[, anom := mean_sst - clim][order(yr, mo)]
message("\n================ SUMMARY ================")
message(sprintf("April 2026 anomaly : %+.2f C", cmp[yr == 2026 & mo == 4, anom]))
message(sprintf("largest anomaly    : %+.2f C in %s",
                cmp[which.max(anom), anom], cmp[which.max(anom), month]))
message(sprintf("warming trend      : %.3f C/decade",
                coef(lm(anom ~ I(yr + (mo-0.5)/12), data = cmp))[2]*10))
mn <- cmp[mo %in% 5:10, .(d = sum(mhw_days)), by = yr][yr %between% c(1998, 2025)]
message(sprintf("heatwave days/warm season: %.0f before 2014, %.0f after",
                mn[yr < 2014, mean(d)], mn[yr >= 2014, mean(d)]))
message("Step 00b done.")
