# 03e_buffer_climate.R
# -----------------------------------------------------------
# THE DECISIVE TEST for "climate-compressed hyperstability".
#
# Buffer strength, following the strategy note:
#
#     Phi(B,C) = 1 - d log(P) / d log(B)
#
#   Phi > 0  production falls proportionally less than biomass: a buffer
#   Phi = 0  production and biomass move together
#   Phi < 0  production falls faster: the compensation has inverted
#
# Climate erodes the buffer if and only if  d Phi / d C < 0.
#
# Model, at the reef-year level:
#
#   log(P_rt) = a_r + d_t + b1*log(B_rt) + b2*log(B_rt):C_rt + b3*C_rt + e_rt
#
#   a_r  reef fixed effects   (permanent differences between reefs)
#   d_t  year fixed effects   (anything common to all reefs in a year)
#   C_rt reef-specific thermal exposure
#
# Then  d log(P)/d log(B) = b1 + b2*C,  so  Phi(C) = 1 - b1 - b2*C
# and   d Phi / d C = -b2.  The buffer erodes with warming iff b2 > 0.
#
# WHY THIS CAN BE IDENTIFIED WHERE THE EARLIER ATTRIBUTION COULD NOT.
# A Gulf-wide annual anomaly is collinear with year fixed effects and is
# therefore not identified (this is the attribution limit stated in the
# Methods). Here C varies BETWEEN reefs WITHIN a year, because thermal
# exposure is resolved to one-degree latitude bands. The interaction is
# identified from reefs that ran hotter than their contemporaries in the
# same year, not from the passage of time. That variation is real but
# modest: the panel spans a handful of latitude bands, so the test is
# reported with its power limits, and with the full inference battery the
# strategy note requires.
#
# Inference: two-way cluster-robust (reef and year) standard errors,
# a wild cluster bootstrap on year (27 clusters), leave-one-year-out and
# leave-one-reef-out, sensitivity to dropping the enforced reserve and the
# extreme thermal years, the whole community against commercial species,
# and a smooth (GAM) alternative compared with the linear interaction.
#
# Requires: ../data/ltem.parquet, ../data/ltem_fish_traits.csv,
#           ../data/ltem_name_lookup.csv,
#           ../data/env/sst_gulf_by_lat_degree_daily.csv
#           data/artisanal_bcs_annual.csv
# Writes:   data/buffer_phi_models.csv      coefficient table, all specs
#           data/buffer_phi_curve.csv       Phi(C) with CI, for plotting
#           data/buffer_phi_summary.csv     the headline numbers
#           data/buffer_reef_year.csv       the reef-year panel itself
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(rfishprod); library(fixest)
})
setwd("..")
DATA <- "data"; OUT <- "data"
set.seed(20260805)
data(db, package = "rfishprod")
MIN_YEARS <- 6      # widest panel: maximises spatial spread of thermal exposure

ltem <- as.data.table(read_parquet("../data/ltem.parquet"))
ltem[, Year := as.integer(Year)]
ltem <- ltem[Year != 2020]

trans_keys <- c("Year","Region","Reef","Habitat","Depth2","Transect",
                "Latitude","Degree","Area")
all_trans <- unique(ltem[Label == "PEC", ..trans_keys])
core <- all_trans[, .(yrs = uniqueN(Year)), by = Reef][yrs >= MIN_YEARS, Reef]
message(sprintf("Buffer panel: %d reefs surveyed in >= %d years", length(core), MIN_YEARS))

# --- production per transect (same construction as 03b) -----------------
fish <- ltem[Label == "PEC" & Reef %in% core & Taxa2 == "Actinopterygii" &
             !is.na(Size) & Size > 0 & !is.na(Quantity) & Quantity > 0 &
             !is.na(A_ord) & !is.na(B_pen)]
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

# long-term band temperature drives Kmax, NOT the temperature of the survey
# year, so production cannot track warming by construction (strategy note 9)
sstd <- fread(file.path("../data/env", "sst_gulf_by_lat_degree_daily.csv"))
sstd[, date := as.Date(date)]
sstd[, `:=`(yr = as.integer(format(date, "%Y")), mo = as.integer(format(date, "%m")))]
sst_lt <- sstd[, .(sst_mean_lt = mean(sst_mean, na.rm = TRUE)), by = lat_degree]
fish[, Degree := as.integer(Degree)]
fish[is.na(Degree), Degree := as.integer(floor(Latitude))]
fish <- merge(fish, sst_lt, by.x = "Degree", by.y = "lat_degree", all.x = TRUE)
fish <- fish[!is.na(sst_mean_lt)]

grid <- unique(fish[, .(species_std, MaxSizeTL, Diet, Position, Method,
                        sstmean = round(sst_mean_lt, 2))])
g <- tidytrait(as.data.frame(grid), db)
kmax <- as.data.table(predKmax(g, dataset = db,
                               fmod = formula(~ sstmean + MaxSizeTL + Diet + Position + Method),
                               niter = 100, return = "pred")$pred)
message(sprintf("Kmax predicted for %d species x temperature combinations", nrow(kmax)))
fish[, sstmean := round(sst_mean_lt, 2)]
fish <- merge(fish, kmax[, .(species_std, sstmean, Kmax)],
              by = c("species_std", "sstmean"), all.x = TRUE)
fish <- fish[!is.na(Kmax)]
fish[, w := A_ord * Size^B_pen]
fish[, growth_day := somaGain(a = A_ord, b = B_pen, Lmeas = Size, t = 1,
                              Lmax = MaxSizeTL, Kmax = Kmax, silent = TRUE)]

art_g  <- fread(file.path(OUT, "artisanal_bcs_annual.csv"))
landed <- unique(toupper(art_g[is_reef == TRUE & landings_t > 0, genus]))
landed <- landed[!is.na(landed) & landed != ""]
fish[, commercial := toupper(Genus) %in% landed]

agg_transect <- function(d) {
  rec <- d[, .(b = sum(w * Quantity, na.rm = TRUE),
               p = sum(growth_day * Quantity, na.rm = TRUE)), by = trans_keys]
  full <- merge(all_trans[Reef %in% core], rec, by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]; full[is.na(p), p := 0]
  full[, `:=`(b100 = b / pmax(Area, 1) * 100, p100 = p / pmax(Area, 1) * 100)]
  full[]
}
tr_all  <- agg_transect(fish)
tr_comm <- agg_transect(fish[commercial == TRUE])

# --- reef-year panel ----------------------------------------------------
to_reef_year <- function(d) {
  d[, .(B = mean(b100, na.rm = TRUE), P = mean(p100, na.rm = TRUE),
        n_transects = .N, Degree = first(Degree), Region = first(Region)),
    by = .(Reef, Year)][B > 0 & P > 0]
}
ry_all  <- to_reef_year(tr_all)
ry_comm <- to_reef_year(tr_comm)

# --- reef-specific thermal exposure C_rt --------------------------------
# Anomalies are computed within each latitude band against that band's own
# 1991-2020 monthly climatology, so C carries between-band variation within
# every year. MHW days use the band's 90th percentile over 1982-2011.
clim <- sstd[yr %between% c(1991, 2020),
             .(clim = mean(sst_mean, na.rm = TRUE)), by = .(lat_degree, mo)]
sstd <- merge(sstd, clim, by = c("lat_degree", "mo"))
sstd[, anom := sst_mean - clim]
thr <- sstd[yr %between% c(1982, 2011),
            .(p90 = quantile(sst_mean, 0.90, na.rm = TRUE)), by = .(lat_degree, mo)]
sstd <- merge(sstd, thr, by = c("lat_degree", "mo"))
env <- sstd[, .(C_ann  = mean(anom, na.rm = TRUE),
                C_warm = mean(anom[mo %in% 5:10], na.rm = TRUE),
                mhw_days = sum(sst_mean > p90, na.rm = TRUE)),
            by = .(Degree = lat_degree, Year = yr)]
message("Thermal exposure: ", uniqueN(env$Degree), " latitude bands x ",
        uniqueN(env$Year), " years")

attach_env <- function(ry) {
  d <- merge(ry, env, by = c("Degree", "Year"))
  d[, `:=`(lB = log(B), lP = log(P))]
  d[, `:=`(lB_c = lB - mean(lB), C = C_warm - mean(C_warm))]   # centred
  d[, Reef := factor(Reef)][, YearF := factor(Year)]
  d[]
}
d_all  <- attach_env(ry_all)
d_comm <- attach_env(ry_comm)
fwrite(d_all, file.path(OUT, "buffer_reef_year.csv"))
message(sprintf("Reef-year panel: %d rows, %d reefs, %d years; C range %.2f to %.2f C",
                nrow(d_all), uniqueN(d_all$Reef), uniqueN(d_all$Year),
                min(d_all$C_warm), max(d_all$C_warm)))
message(sprintf("Within-year spread of C (mean sd across years): %.3f C",
                d_all[, .(s = sd(C_warm)), by = Year][, mean(s, na.rm = TRUE)]))

# --- the model ----------------------------------------------------------
# b2 = coefficient on lB_c:C.  Phi(C) = 1 - (b1 + b2*C).  dPhi/dC = -b2.
fit_buffer <- function(d) feols(lP ~ lB_c * C | Reef + YearF, data = d)

grab <- function(fit, d, label, cl = "twoway") {
  V <- if (cl == "twoway") vcov(fit, cluster = ~ Reef + YearF)
       else if (cl == "reef") vcov(fit, cluster = ~ Reef)
       else vcov(fit, cluster = ~ YearF)
  cf <- coef(fit); se <- sqrt(diag(V))
  nm <- "lB_c:C"
  if (!nm %in% names(cf)) return(NULL)
  b2 <- cf[[nm]]; s2 <- se[[nm]]
  b1 <- cf[["lB_c"]]
  data.table(spec = label, cluster = cl,
             b1 = b1, b2 = b2, se_b2 = s2,
             t = b2 / s2, p = 2 * pnorm(-abs(b2 / s2)),
             lo = b2 - 1.96 * s2, hi = b2 + 1.96 * s2,
             phi_at_mean = 1 - b1,
             dphi_dC = -b2,
             n = nobs(fit), reefs = uniqueN(d$Reef), years = uniqueN(d$Year))
}

res <- list()
m_all <- fit_buffer(d_all)
res[[length(res)+1]] <- grab(m_all, d_all, "Whole community", "twoway")
res[[length(res)+1]] <- grab(m_all, d_all, "Whole community", "reef")
res[[length(res)+1]] <- grab(m_all, d_all, "Whole community", "year")
m_comm <- fit_buffer(d_comm)
res[[length(res)+1]] <- grab(m_comm, d_comm, "Commercial species", "twoway")

# sensitivity: drop the enforced reserve; drop extreme thermal years
d_nocp <- d_all[!grepl("CABO_PULMO", as.character(Reef))]
res[[length(res)+1]] <- grab(fit_buffer(d_nocp), d_nocp, "Excluding Cabo Pulmo", "twoway")
extreme <- d_all[, .(m = mean(C_warm)), by = Year][order(-m)][1:3, Year]
d_noext <- d_all[!Year %in% extreme]
res[[length(res)+1]] <- grab(fit_buffer(d_noext), d_noext, "Excluding 3 hottest years", "twoway")

# sensitivity: narrower, better-sampled panels
for (k in c(8, 10, 15)) {
  keep <- d_all[, .(y = uniqueN(Year)), by = Reef][y >= k, Reef]
  dk <- d_all[Reef %in% keep]
  if (uniqueN(dk$Reef) > 5)
    res[[length(res)+1]] <- grab(fit_buffer(dk), dk, sprintf("Reefs with >= %d years", k), "twoway")
}
# alternative thermal metrics
d_alt <- copy(d_all)[, C := C_ann - mean(C_ann)]
res[[length(res)+1]] <- grab(fit_buffer(d_alt), d_alt, "Annual anomaly instead of warm season", "twoway")
d_mhw <- copy(d_all)[, C := (mhw_days - mean(mhw_days)) / sd(mhw_days)]
res[[length(res)+1]] <- grab(fit_buffer(d_mhw), d_mhw, "Heatwave days instead of anomaly", "twoway")

models <- rbindlist(res, fill = TRUE)
fwrite(models, file.path(OUT, "buffer_phi_models.csv"))
message("\n=== Buffer-climate interaction (b2 > 0 means warming erodes the buffer) ===")
print(models[, .(spec, cluster, b2 = round(b2, 4), se = round(se_b2, 4),
                 p = signif(p, 3), lo = round(lo, 3), hi = round(hi, 3),
                 phi_mean = round(phi_at_mean, 3), n, reefs, years)])

# --- wild cluster bootstrap on year (Rademacher, 27 clusters) -----------
wild_boot <- function(d, B = 999) {
  f0 <- feols(lP ~ lB_c + C | Reef + YearF, data = d)      # restricted: b2 = 0
  r0 <- resid(f0); fit0 <- predict(f0)
  yrs <- unique(d$YearF)
  t_obs <- { m <- fit_buffer(d)
             coef(m)[["lB_c:C"]] / sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[["lB_c:C"]] }
  ts <- numeric(B)
  for (b in seq_len(B)) {
    w <- setNames(sample(c(-1, 1), length(yrs), TRUE), as.character(yrs))
    db <- copy(d)[, lP := fit0 + r0 * w[as.character(YearF)]]
    m <- feols(lP ~ lB_c * C | Reef + YearF, data = db)
    ts[b] <- coef(m)[["lB_c:C"]] / sqrt(diag(vcov(m, cluster = ~ Reef + YearF)))[["lB_c:C"]]
  }
  list(t_obs = t_obs, p = mean(abs(ts) >= abs(t_obs)))
}
wb <- wild_boot(d_all, B = 999)
message(sprintf("\nWild cluster bootstrap on year: t = %.3f, p = %.3f", wb$t_obs, wb$p))

# --- leave-one-out ------------------------------------------------------
loo <- function(d, by) {
  vals <- unique(d[[by]])
  sapply(vals, function(v) {
    dd <- d[get(by) != v]
    m <- fit_buffer(dd)
    coef(m)[["lB_c:C"]]
  })
}
loo_y <- loo(d_all, "Year"); loo_r <- loo(d_all, "Reef")
b2_main <- models[spec == "Whole community" & cluster == "twoway", b2]
message(sprintf("Leave-one-year-out b2: range %.4f to %.4f, sign stable: %s",
                min(loo_y), max(loo_y), all(sign(loo_y) == sign(b2_main))))
message(sprintf("Leave-one-reef-out b2: range %.4f to %.4f, sign stable: %s",
                min(loo_r), max(loo_r), all(sign(loo_r) == sign(b2_main))))

# --- smooth alternative -------------------------------------------------
gam_ok <- FALSE
try({
  suppressPackageStartupMessages(library(mgcv))
  gm <- bam(lP ~ s(lB_c) + s(C) + ti(lB_c, C) + s(Reef, bs = "re") + s(YearF, bs = "re"),
            data = d_all, discrete = TRUE)
  ti_p <- summary(gm)$s.table[grep("ti\\(", rownames(summary(gm)$s.table)), "p-value"]
  message(sprintf("GAM tensor interaction ti(logB, C): p = %.3g", ti_p))
  gam_ok <- TRUE
}, silent = TRUE)

# --- Phi(C) curve -------------------------------------------------------
Vm <- vcov(m_all, cluster = ~ Reef + YearF)
b1 <- coef(m_all)[["lB_c"]]; b2 <- coef(m_all)[["lB_c:C"]]
Cs <- seq(quantile(d_all$C, 0.05), quantile(d_all$C, 0.95), length.out = 60)
phi <- 1 - (b1 + b2 * Cs)
va  <- Vm["lB_c","lB_c"] + Cs^2 * Vm["lB_c:C","lB_c:C"] + 2 * Cs * Vm["lB_c","lB_c:C"]
curve <- data.table(C_centred = Cs, C_warm = Cs + mean(d_all$C_warm),
                    phi = phi, lo = phi - 1.96 * sqrt(va), hi = phi + 1.96 * sqrt(va))
fwrite(curve, file.path(OUT, "buffer_phi_curve.csv"))

survives <- models[spec == "Whole community" & cluster == "twoway", p] < 0.05 &&
            b2_main > 0 && wb$p < 0.05 && all(sign(loo_y) == sign(b2_main))
summ <- data.table(
  quantity = c("phi_at_mean_C", "dphi_dC", "dphi_dC_se", "dphi_dC_p_twoway",
               "dphi_dC_p_wildboot", "phi_cold", "phi_warm",
               "buffer_panel_reefs", "buffer_panel_years", "buffer_panel_n",
               "within_year_C_sd", "climate_buffer_erosion_supported"),
  value = c(round(1 - b1, 3), round(-b2, 4),
            round(models[spec == "Whole community" & cluster == "twoway", se_b2], 4),
            signif(models[spec == "Whole community" & cluster == "twoway", p], 3),
            round(wb$p, 3),
            round(curve[1, phi], 3), round(curve[.N, phi], 3),
            uniqueN(d_all$Reef), uniqueN(d_all$Year), nrow(d_all),
            round(d_all[, .(s = sd(C_warm)), by = Year][, mean(s, na.rm = TRUE)], 3),
            as.integer(survives)))
fwrite(summ, file.path(OUT, "buffer_phi_summary.csv"))
print(summ)
message(sprintf("\n>>> VERDICT: climate-buffer erosion %s under the full battery.",
                ifelse(survives, "SURVIVES", "is NOT supported")))
message("Step 03e done.")
