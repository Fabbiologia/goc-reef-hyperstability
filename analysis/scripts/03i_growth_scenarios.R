# 03i_growth_scenarios.R
# -----------------------------------------------------------
# How conservative is the growth model's temperature treatment?
#
# The production engine (03b, 03c, 03e) predicts each species' Kmax at the
# reef's 1991 to 2020 climatological mean temperature, so growth capacity
# never responds to warming and every production estimate in the warming
# half of the record is a floor: Kmax rises with temperature in the trait
# model, and holding temperature at the climatology withholds that rise.
# R. Morais (the model's author) suggested bracketing this choice with
# explicit scenarios and running temperature windows. This step does that.
#
# Three scenarios, identical fish, identical model, different temperature:
#   conservative  Kmax at the reef's climatological mean (the manuscript's
#                 primary estimate everywhere; growth capacity fixed, so
#                 production change is community change only)
#   tracking      Kmax at the mean SST of the reef's own OISST cell over
#                 the 365 days before that reef year's median survey date;
#                 the fitted spatial Kmax to temperature gradient is read
#                 as a fully expressed temporal response
#   extreme       Kmax at the survey year's warm season mean (May to Oct)
#                 of the same cell: the fish grows all year at the hottest
#                 season's temperature, the upper bracket
#
# Reality sits between conservative and tracking: fish acclimatise and
# adapt on the spatial gradient in ways an anomaly year does not allow,
# so the spatial coefficient overstates the instantaneous response. The
# envelope, not any single line, is the estimate ("hybridisation").
#
# The tracking scenario is also computed with 30, 90 and 180 day windows
# to show how the choice of temperature integration window matters.
#
# The scenarios feed the LEVELS results (production and turnover trends,
# the buffer strength Phi). They must NOT feed the climate interaction
# tests of 03e/03f/03h: putting observed reef-year temperature into the
# growth model bakes thermal exposure into the response, and the
# interaction would partly measure the model's own assumption. Those
# tests stay on the temperature-blind conservative engine by design.
#
# Reads:  ../data/ltem.parquet, ../data/ltem_name_lookup.csv,
#         ../data/ltem_fish_traits.csv, ../data/env/sst_reef_longterm.csv,
#         ../data/env/sst_reef_daily.rds, ../data/env/sst_reef_year.csv,
#         data/artisanal_bcs_annual.csv
# Writes: data/growth_scenario_trends.csv     trends per scenario, 26-reef panel
#         data/growth_scenario_phi.csv        Phi per scenario, >=6-year panel
#         data/growth_scenario_wedge.csv      annual production gap vs conservative
#         data/growth_scenario_windows.csv    window sensitivity of the tracking run
#         data/growth_scenario_summary.csv    in-text numbers
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(arrow); library(rfishprod); library(fixest)
})
setwd("..")
DATA <- "data"; OUT <- "data"; ENV <- "../data/env"
set.seed(20260805)                     # same seed as 03b
data(db, package = "rfishprod")

# -----------------------------------------------------------
# 1. Fish records, exactly as 03b builds them, for ALL reefs
# -----------------------------------------------------------
ltem <- as.data.table(read_parquet("../data/ltem.parquet"))
ltem[, Year := as.integer(Year)]
ltem <- ltem[Year != 2020]
trans_keys <- c("Year","Region","Reef","Habitat","Depth2","Transect",
                "Latitude","Degree","Area")
all_trans <- unique(ltem[Label == "PEC", ..trans_keys])
reef_yrs  <- all_trans[, .(yrs = uniqueN(Year)), by = Reef]
core      <- reef_yrs[yrs >= 15, Reef]     # the balanced 26-reef trend panel
wide      <- reef_yrs[yrs >= 6,  Reef]     # the 03e buffer panel
fish_core_t <- all_trans[Reef %in% core]

pec  <- ltem[Label == "PEC"]
fish <- pec[Taxa2 == "Actinopterygii" & !is.na(Size) & Size > 0 &
            !is.na(Quantity) & Quantity > 0 & !is.na(A_ord) & !is.na(B_pen)]
lookup <- fread(file.path("../data", "ltem_name_lookup.csv"))
traits <- fread(file.path("../data", "ltem_fish_traits.csv"))
fish[, species_std := Species]
fish[lookup, on = .(species_std = species_raw), species_std := i.species]
fish <- merge(fish, traits[, .(species_std = species, MaxSizeTL, Diet, Position, Method)],
              by = "species_std", all.x = TRUE)
modal <- function(x) names(sort(table(x), decreasing = TRUE))[1]
gen_tr <- traits[, .(MaxSizeTL = max(MaxSizeTL), Diet = modal(Diet),
                     Position = modal(Position), Method = modal(Method)), by = genus]
fam_tr <- traits[, .(MaxSizeTL = max(MaxSizeTL), Diet = modal(Diet),
                     Position = modal(Position), Method = modal(Method)), by = family]
fish[is.na(MaxSizeTL), c("MaxSizeTL","Diet","Position","Method") :=
       gen_tr[.SD, on = .(genus = Genus), .(MaxSizeTL, Diet, Position, Method)]]
fish[is.na(MaxSizeTL), c("MaxSizeTL","Diet","Position","Method") :=
       fam_tr[.SD, on = .(family = Family), .(MaxSizeTL, Diet, Position, Method)]]
fish <- fish[!is.na(MaxSizeTL)]
fish[, Size := pmin(Size, MaxSizeTL)]
fish <- fish[Reef %in% union(core, wide)]
message(sprintf("Fish table: %s rows, %d reefs (%d core, %d buffer-panel)",
                format(nrow(fish), big.mark = ","), uniqueN(fish$Reef),
                length(core), length(wide)))

# -----------------------------------------------------------
# 2. The temperature each reef year actually experienced
# -----------------------------------------------------------
reef_lt <- fread(file.path(ENV, "sst_reef_longterm.csv"))[, .(Reef, sst_mean_lt, cell)]
daily   <- readRDS(file.path(ENV, "sst_reef_daily.rds"))
setkey(daily, cell, date)

# median survey date per reef year
sd_med <- unique(ltem[Label == "PEC", .(Reef, Year, Month, Day)])[
  , .(sdate = median(as.IDate(sprintf("%d-%02d-%02d", Year, Month, pmax(Day, 1))))),
  by = .(Reef, Year)]
missing_cell <- setdiff(union(core, wide), reef_lt$Reef)
if (length(missing_cell)) message("  reefs without an OISST cell, excluded: ",
                                  paste(missing_cell, collapse = ", "))
core <- intersect(core, reef_lt$Reef); wide <- intersect(wide, reef_lt$Reef)
sd_med <- sd_med[Reef %in% union(core, wide)]
sd_med <- merge(sd_med, reef_lt[, .(Reef, cell)], by = "Reef")

# running-window means before the survey date, per reef year
win_mean <- function(W) {
  q <- copy(sd_med)[, `:=`(d0 = sdate - W, d1 = sdate - 1L)]
  j <- daily[q, on = .(cell, date >= d0, date <= d1),
             .(t = mean(temp, na.rm = TRUE), nd = .N), by = .EACHI]
  out <- cbind(q[, .(Reef, Year)], j[, .(t, nd)])
  stopifnot(out[, all(nd > 0.9 * W)])
  out[, .(Reef, Year, t)]
}
tw <- list(); for (W in c(30, 90, 180, 365)) tw[[as.character(W)]] <- win_mean(W)
temps <- Reduce(function(a, b) merge(a, b, by = c("Reef", "Year")),
                Map(function(d, w) setnames(copy(d), "t", paste0("t", w)),
                    tw, names(tw)))

# warm-season mean of the survey year, per reef
warm <- unique(daily[mo %in% 5:10, .(t_warm = mean(temp)), by = .(cell, Year = yr)])
temps <- merge(temps, merge(sd_med[, .(Reef, Year, cell)], warm,
                            by = c("cell", "Year"))[, .(Reef, Year, t_warm)],
               by = c("Reef", "Year"))
temps <- merge(temps, reef_lt[, .(Reef, t_clim = sst_mean_lt)], by = "Reef")
message(sprintf("Reef-year temperatures: %d reef years; tracking (365 d) runs %.2f to %.2f C, climatology %.2f to %.2f C",
                nrow(temps), min(temps$t365), max(temps$t365),
                min(temps$t_clim), max(temps$t_clim)))

# -----------------------------------------------------------
# 3. Kmax for every species x temperature the scenarios need
# -----------------------------------------------------------
fish <- merge(fish, temps, by = c("Reef", "Year"))
scenario_temp <- c(conservative = "t_clim", tracking = "t365", extreme = "t_warm",
                   track30 = "t30", track90 = "t90", track180 = "t180")
for (v in scenario_temp) fish[, (v) := round(get(v), 1)]

grid <- unique(rbindlist(lapply(scenario_temp, function(v)
  fish[, .(species_std, MaxSizeTL, Diet, Position, Method, sstmean = get(v))])))
setorder(grid, species_std, sstmean)
message(sprintf("Kmax grid: %s species x temperature combinations", format(nrow(grid), big.mark = ",")))
g <- tidytrait(as.data.frame(grid), db)
kmax <- as.data.table(predKmax(g, dataset = db,
                               fmod = formula(~ sstmean + MaxSizeTL + Diet + Position + Method),
                               niter = 100, return = "pred")$pred)
kmax <- kmax[, .(species_std, sstmean, Kmax)]

# implied temperature sensitivity of the trait model itself, for the text
ks <- kmax[, if (.N > 3 && diff(range(sstmean)) > 1)
             .(el = coef(lm(log(Kmax) ~ sstmean))[["sstmean"]]), by = species_std]
message(sprintf("Kmax temperature elasticity across species: median %+.1f%% per C",
                100 * (exp(median(ks$el)) - 1)))

# -----------------------------------------------------------
# 4. Production per transect under each scenario
# -----------------------------------------------------------
fish[, w := A_ord * Size^B_pen]
per_scenario <- function(v) {
  f <- merge(fish, kmax, by.x = c("species_std", v), by.y = c("species_std", "sstmean"))
  stopifnot(nrow(f) == nrow(fish))
  f[, growth_day := somaGain(a = A_ord, b = B_pen, Lmeas = Size, t = 1,
                             Lmax = MaxSizeTL, Kmax = Kmax, silent = TRUE)]
  f[, .(b = sum(w * Quantity, na.rm = TRUE),
        p = sum(growth_day * Quantity, na.rm = TRUE)), by = trans_keys]
}

zero_fill <- function(rec) {
  full <- merge(fish_core_t, rec, by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]; full[is.na(p), p := 0]
  full[, b100 := b / pmax(Area, 1) * 100]
  full[, p100 := p / pmax(Area, 1) * 100]
  full[, turnover := fifelse(b > 0, 100 * p / b, NA_real_)]
  full[]
}

# 03b's estimators, verbatim
clust_p <- function(fit, cluster, term) {
  cf <- coef(fit); k <- names(cf)[!is.na(cf)]
  X  <- model.matrix(fit)[, k, drop = FALSE]
  u  <- residuals(fit); g <- factor(cluster[as.integer(names(u))])
  br <- chol2inv(chol(crossprod(X))); mt <- matrix(0, ncol(X), ncol(X))
  for (l in levels(g)) { i <- which(g == l)
    xu <- crossprod(X[i, , drop = FALSE], u[i]); mt <- mt + tcrossprod(xu) }
  m <- nlevels(g)
  V <- br %*% mt %*% br * (m / (m - 1)) * ((nrow(X) - 1) / (nrow(X) - ncol(X)))
  se <- sqrt(diag(V)); names(se) <- k
  2 * pt(-abs(cf[term] / se[term]), df = m - 1)
}
trend_one <- function(d, value, offset = 0.01) {
  d <- copy(d)[!is.na(get(value))]
  d[, l_v := log(get(value) + offset)]
  fit <- lm(l_v ~ Year + Reef, data = d)
  pr  <- d[, .(e = mean(get(value)[Year <= 2013], na.rm = TRUE),
               l = mean(get(value)[Year >= 2014], na.rm = TRUE)), by = Reef
           ][is.finite(e) & is.finite(l) & e > 0 & l > 0]
  data.table(pct_per_decade = 100 * (exp(coef(fit)["Year"] * 10) - 1),
             p_reefclust    = clust_p(fit, d$Reef, "Year"),
             paired_pct     = 100 * (exp(mean(log(pr$l / pr$e))) - 1))
}

phi_one <- function(rec) {
  # identical construction to 03e: zero-fill the surveyed transects, then
  # reef-year MEANS of the per-100 m2 rates, then B > 0 & P > 0
  ry <- merge(all_trans[Reef %in% wide], rec, by = trans_keys, all.x = TRUE)
  ry[is.na(b), b := 0]; ry[is.na(p), p := 0]
  ry[, `:=`(b100 = b / pmax(Area, 1) * 100, p100 = p / pmax(Area, 1) * 100)]
  ry <- ry[, .(B = mean(b100, na.rm = TRUE), P = mean(p100, na.rm = TRUE)),
           by = .(Reef, Year)][B > 0 & P > 0]
  ry[, `:=`(lP = log(P), lB_c = log(B) - mean(log(B)),
            Reef = factor(Reef), YearF = factor(Year))]
  m  <- feols(lP ~ lB_c | Reef + YearF, data = ry)
  b  <- coef(m)[["lB_c"]]
  se <- sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[["lB_c"]]
  data.table(phi = 1 - b, lo = 1 - b - 1.96 * se, hi = 1 - b + 1.96 * se,
             n = nrow(ry), reefs = uniqueN(ry$Reef))
}

resid_traj <- function(d, value, offset = 0.01) {
  d <- copy(d)[!is.na(get(value))]
  d[, l_v := log(get(value) + offset)]
  d[, r := l_v - mean(l_v, na.rm = TRUE), by = Reef]
  d[, .(mean = mean(r, na.rm = TRUE)), by = Year][order(Year)]
}

scen_main <- c("conservative", "tracking", "extreme")
trends <- list(); phis <- list(); traj <- list(); winds <- list()
for (sc in names(scenario_temp)) {
  rec <- per_scenario(scenario_temp[[sc]])
  if (sc %in% scen_main) {
    d <- zero_fill(rec[Reef %in% core])
    trends[[sc]] <- rbind(
      data.table(scenario = sc, response = "production", trend_one(d, "p100")),
      data.table(scenario = sc, response = "turnover",   trend_one(d, "turnover", offset = 0.001)))
    phis[[sc]] <- data.table(scenario = sc, phi_one(rec))
    traj[[sc]] <- data.table(scenario = sc, resid_traj(d, "p100"))
    message(sprintf("  %s done", sc))
  }
  if (grepl("^track", sc) || sc == "tracking") {
    W <- if (sc == "tracking") 365 else as.integer(sub("track", "", sc))
    d <- zero_fill(rec[Reef %in% core])
    winds[[sc]] <- data.table(window_days = W,
                              data.table(response = "production", trend_one(d, "p100")))
  }
}
trends <- rbindlist(trends); phis <- rbindlist(phis)
traj <- rbindlist(traj);     winds <- rbindlist(winds)[order(window_days)]

fwrite(trends, file.path(OUT, "growth_scenario_trends.csv"))
fwrite(phis,   file.path(OUT, "growth_scenario_phi.csv"))
fwrite(winds,  file.path(OUT, "growth_scenario_windows.csv"))

message("\nTrends on the 26-reef panel by scenario:")
print(trends[, .(scenario, response, pct_per_decade = round(pct_per_decade, 1),
                 p = signif(p_reefclust, 2), paired_pct = round(paired_pct, 1))])
message("\nBuffer strength Phi by scenario (>=6-year panel):")
print(phis[, .(scenario, phi = round(phi, 3), lo = round(lo, 3), hi = round(hi, 3), n, reefs)])
message("\nTracking-scenario production trend by temperature window:")
print(winds[, .(window_days, pct_per_decade = round(pct_per_decade, 1), paired_pct = round(paired_pct, 1))])

# -----------------------------------------------------------
# 5. The wedge: how the scenario gap moves with temperature
# -----------------------------------------------------------
wide_traj <- dcast(traj, Year ~ scenario, value.var = "mean")
anom <- fread(file.path(ENV, "sst_reef_year.csv"))[
  Reef %in% core, .(anom = mean(anom_annual)), by = .(Year = year)]
wedge <- merge(wide_traj, anom, by = "Year")
wedge[, `:=`(gap_tracking = 100 * (exp(tracking - conservative) - 1),
             gap_extreme  = 100 * (exp(extreme  - conservative) - 1))]
fwrite(wedge, file.path(OUT, "growth_scenario_wedge.csv"))
r_wedge <- wedge[, cor(gap_tracking, anom)]
message(sprintf("\nWedge: tracking runs %+.1f%% above conservative in the warmest year, %+.1f%% in the coolest; corr with the annual anomaly r = %.2f",
                wedge[which.max(anom), gap_tracking], wedge[which.min(anom), gap_tracking], r_wedge))

summ <- data.table(quantity = c(
  "kmax_elasticity_pct_per_C",
  "prod_trend_conservative", "prod_trend_tracking", "prod_trend_extreme",
  "turn_trend_conservative", "turn_trend_tracking", "turn_trend_extreme",
  "phi_conservative", "phi_tracking", "phi_extreme",
  "phi_lo_tracking", "phi_hi_tracking",
  "wedge_warmest_year_pct", "wedge_coolest_year_pct", "wedge_anom_corr",
  "prod_trend_window30", "prod_trend_window90", "prod_trend_window180", "prod_trend_window365"),
  value = c(
    round(100 * (exp(median(ks$el)) - 1), 1),
    trends[response == "production", round(pct_per_decade, 1)],
    trends[response == "turnover",   round(pct_per_decade, 1)],
    phis[, round(phi, 3)],
    round(phis[scenario == "tracking", lo], 3), round(phis[scenario == "tracking", hi], 3),
    round(wedge[which.max(anom), gap_tracking], 1), round(wedge[which.min(anom), gap_tracking], 1),
    round(r_wedge, 2),
    winds[, round(pct_per_decade, 1)]))
fwrite(summ, file.path(OUT, "growth_scenario_summary.csv"))
message("\nStep 03i done. Growth scenario tables written to manuscript/data/.")
