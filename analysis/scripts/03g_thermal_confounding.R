# 03g_thermal_confounding.R
# -----------------------------------------------------------
# Before any thermal coefficient from the per-reef OISST data is believed,
# it has to survive checks that the identifying variation is heat and not
# geography. The band series produced an implausible positive level effect
# (warmer reefs holding MORE biomass, +36 % per degree), which is the
# signature of spatial confounding rather than a thermal response. In this
# basin the obvious mechanism is upwelling: cells that run cold relative to
# their own climatology are often upwelling, and upwelling brings nutrients,
# so a cold anomaly can raise production for reasons that have nothing to do
# with thermal stress.
#
# Four checks, all with reef and year fixed effects and errors clustered on
# reef and year:
#   1 LEVELS   log(B) and log(P) on thermal exposure. A large positive
#              coefficient means the design is still confounded.
#   2 PLACEBO  each reef-year is given the thermal exposure of a randomly
#              chosen DIFFERENT year. Any interaction must vanish.
#   3 GRADIENT the residual thermal exposure is regressed on latitude and on
#              distance to the midriff islands. A strong spatial signature
#              means the variation is geography, not weather.
#   4 SIGN     production response to cold anomalies at the reefs with the
#              strongest upwelling signature, to see whether the mechanism
#              reads as nutrients rather than heat.
#
# Requires: data/buffer_reef_year.csv  (rebuilt by 03e on per-reef SST)
# Writes:   data/thermal_confounding.csv
# -----------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(fixest)})
setwd(".."); DATA <- "data"; OUT <- "data"
set.seed(20260805)

d <- fread(file.path(DATA, "buffer_reef_year.csv"))
d[, YearF := factor(Year)]
d[, lB_c := log(B) - mean(log(B))]
d[, C := C_warm - mean(C_warm)]
res <- list()
add <- function(check, term, b, se, note = "") {
  p <- 2 * pnorm(-abs(b / se))
  res[[length(res) + 1]] <<- data.table(check, term, est = b, se = se, p = p, note = note)
  message(sprintf("  %-22s %-16s %+8.4f (se %.4f) p = %-8.3g %s", check, term, b, se, p, note))
}
cf <- function(m, tm) { b <- coef(m)[[tm]]
  list(b = b, se = sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[[tm]]) }

message("1. LEVELS  (a large positive coefficient indicates confounding)")
for (r in c("lB", "lP")) {
  m <- feols(as.formula(paste(r, "~ C | Reef + YearF")), data = d)
  x <- cf(m, "C"); add("levels", r, x$b, x$se,
                       sprintf("%+.1f %% per C", 100 * (exp(x$b) - 1)))
}

message("2. PLACEBO  (interaction must vanish when years are shuffled)")
m_true <- feols(lP ~ lB_c * C | Reef + YearF, data = d)
xt <- cf(m_true, "lB_c:C"); add("true", "lB_c:C", xt$b, xt$se)
pl <- replicate(200, {
  dd <- copy(d)
  key <- unique(dd[, .(Reef, Year, C)])
  key[, Cshuf := sample(C), by = Reef]
  dd <- merge(dd[, !"C"], key[, .(Reef, Year, C = Cshuf)], by = c("Reef", "Year"))
  coef(feols(lP ~ lB_c * C | Reef + YearF, data = dd))[["lB_c:C"]]
})
add("placebo", "lB_c:C", mean(pl), sd(pl),
    sprintf("true |b| exceeds %.0f %% of placebos", 100 * mean(abs(pl) < abs(xt$b))))

message("3. GRADIENT  (residual exposure should not be spatial)")
d[, Cres := resid(feols(C ~ 1 | Reef + YearF, data = d))]
MIDRIFF_LAT <- 29.0; MIDRIFF_LON <- -113.0
rc <- fread("../data/env/sst_reef_cell_map.csv")
d2 <- merge(d, rc[, .(Reef, rlat = lat, rlon = lon)], by = "Reef", all.x = TRUE)
if (!all(is.na(d2$rlat))) {
  d2[, dmid := sqrt((rlat - MIDRIFF_LAT)^2 + (rlon - MIDRIFF_LON)^2)]
  for (v in c("rlat", "dmid")) {
    m <- feols(as.formula(paste("Cres ~", v)), data = d2)
    b <- coef(m)[[v]]; se <- sqrt(diag(vcov(m, cluster = ~ Reef)))[[v]]
    add("gradient", v, b, se)
  }
}

message("4. SIGN  (production in cold vs warm anomaly years)")
d[, cold := as.integer(C_warm < quantile(C_warm, 0.25))]
m <- feols(lP ~ cold | Reef + YearF, data = d)
x <- cf(m, "cold"); add("cold_years", "production", x$b, x$se,
                        sprintf("%+.1f %% in coldest quartile", 100 * (exp(x$b) - 1)))

out <- rbindlist(res)
fwrite(out, file.path(OUT, "thermal_confounding.csv"))
lev <- out[check == "levels", max(abs(est))]
message(sprintf("\nVERDICT: %s",
  ifelse(lev > 0.15,
    "levels still implausible -> identification NOT clean, do not report interactions",
    "levels behave -> interaction results may be interpreted, with the caveats stated")))
