# 03b_productivity.R
# -----------------------------------------------------------
# Biomass production and turnover on the SAME balanced 26-reef
# panel as the biomass analysis in 03, so the paper can measure
# internally the engine of hyperstability: a reef that renews
# its biomass faster as its standing stock falls keeps supplying
# catch while it declines.
#
# Method follows Morais & Bellwood (2018, 2020) as implemented
# in rfishprod: Kmax (the size-standardised von Bertalanffy
# growth coefficient) is predicted for every taxon from maximum
# length, diet, position, ageing method and the LONG-TERM mean
# SST of the reef's one-degree latitude band; every individual
# is grown forward one day along its growth curve; the length
# increment is converted to mass with the record's own
# length-weight coefficients; daily rates are scaled by 365.
# Using the long-term mean temperature (not the survey year)
# keeps production from tracking temperature by construction.
# Gross production is used throughout (no mortality subtraction),
# matching the definition behind published turnover estimates.
#
# The trait table (data/ltem_fish_traits.csv) is a committed
# input of this repository: the monitoring programme's curated
# trait table first, then the rfishprod growth database and
# FishBase, with gaps filled from congeners then confamilials.
# Elasmobranchs are excluded (the growth model was not built
# for them); their biomass share on the panel is logged below.
#
# Requires: ../data/ltem.parquet
#           ../data/ltem_fish_traits.csv, ../data/ltem_name_lookup.csv
#           ../data/env/sst_gulf_by_lat_degree_daily.csv
#           data/artisanal_bcs_annual.csv        (from 01)
# Writes:   data/productivity_panel_trends.csv   (trend table)
#           data/productivity_fig2cd_residuals.csv (Figure 2c,d)
#           data/productivity_panel_summary.csv  (in-text numbers)
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(rfishprod)
})
setwd("..")
DATA <- "data"; OUT <- "data"
set.seed(20260805)          # fixed seed for the Kmax bootstrap
data(db, package = "rfishprod")

ltem <- as.data.table(read_parquet("../data/ltem.parquet"))
ltem[, Year := as.integer(Year)]
ltem <- ltem[Year != 2020]

trans_keys <- c("Year","Region","Reef","Habitat","Depth2","Transect",
                "Latitude","Degree","Area")

# Balanced panel: identical construction to 03
all_trans <- unique(ltem[Label == "PEC", ..trans_keys])
reef_yrs  <- all_trans[, .(yrs = uniqueN(Year)), by = Reef]
core      <- reef_yrs[yrs >= 15, Reef]
fish_core_t <- all_trans[Reef %in% core]
message(sprintf("Balanced reef panel: %d reefs, %d transects",
                length(core), nrow(fish_core_t)))

# -----------------------------------------------------------
# Fish records eligible for the growth model
# -----------------------------------------------------------
pec <- ltem[Label == "PEC" & Reef %in% core]
elasmo_share <- pec[, 100 * sum(Biomass[Taxa2 != "Actinopterygii"], na.rm = TRUE) /
                        sum(Biomass, na.rm = TRUE)]
message(sprintf("Elasmobranch share of panel biomass (excluded): %.2f%%", elasmo_share))

fish <- pec[Taxa2 == "Actinopterygii" & !is.na(Size) & Size > 0 &
            !is.na(Quantity) & Quantity > 0 & !is.na(A_ord) & !is.na(B_pen)]

# Harmonise names, then join traits; fall back genus -> family for the few
# taxa outside the companion's 186 (fallback shares are logged).
lookup <- fread(file.path("../data", "ltem_name_lookup.csv"))
traits <- fread(file.path("../data", "ltem_fish_traits.csv"))
fish[, species_std := Species]
fish[lookup, on = .(species_std = species_raw), species_std := i.species]

fish <- merge(fish, traits[, .(species_std = species, MaxSizeTL, Diet, Position, Method)],
              by = "species_std", all.x = TRUE)
n_ind <- fish[, sum(Quantity)]
cov_sp <- fish[!is.na(MaxSizeTL), sum(Quantity)] / n_ind * 100

modal <- function(x) names(sort(table(x), decreasing = TRUE))[1]
gen_tr <- traits[, .(MaxSizeTL = max(MaxSizeTL), Diet = modal(Diet),
                     Position = modal(Position), Method = modal(Method)), by = genus]
fam_tr <- traits[, .(MaxSizeTL = max(MaxSizeTL), Diet = modal(Diet),
                     Position = modal(Position), Method = modal(Method)), by = family]
fish[is.na(MaxSizeTL), c("MaxSizeTL","Diet","Position","Method") :=
       gen_tr[.SD, on = .(genus = Genus), .(MaxSizeTL, Diet, Position, Method)]]
cov_gen <- fish[!is.na(MaxSizeTL), sum(Quantity)] / n_ind * 100
fish[is.na(MaxSizeTL), c("MaxSizeTL","Diet","Position","Method") :=
       fam_tr[.SD, on = .(family = Family), .(MaxSizeTL, Diet, Position, Method)]]
cov_fam <- fish[!is.na(MaxSizeTL), sum(Quantity)] / n_ind * 100
dropped <- fish[is.na(MaxSizeTL)]
message(sprintf("Trait coverage of individuals: %.2f%% species, %.2f%% + genus, %.2f%% + family; dropped %.3f%%",
                cov_sp, cov_gen, cov_fam, 100 - cov_fam))
fish <- fish[!is.na(MaxSizeTL)]

# Size entry errors: cap at the published maximum, preserving abundance
trunc_pct <- fish[, 100 * sum(Quantity[Size > MaxSizeTL]) / sum(Quantity)]
fish[, Size := pmin(Size, MaxSizeTL)]
message(sprintf("Size records capped at published maximum: %.2f%% of individuals", trunc_pct))

# -----------------------------------------------------------
# Long-term mean SST per one-degree latitude band
# -----------------------------------------------------------
# Per-reef long-term mean from OISST v2.1 at quarter degree resolution
# (00b). Each reef takes its own nearest ocean cell rather than a one-degree
# band average, so reefs that sit in the midriff cold pool or on the warmer
# mainland side get the temperature they actually experience.
reef_lt <- fread(file.path("../data/env", "sst_reef_longterm.csv"))
fish <- merge(fish, reef_lt[, .(Reef, sst_mean_lt)], by = "Reef", all.x = TRUE)
# Any reef without a matched cell falls back to its latitude band.
if (anyNA(fish$sst_mean_lt)) {
  sstd <- fread(file.path("../data/env", "sst_gulf_by_lat_degree_daily.csv"))
  bl <- sstd[, .(band_lt = mean(sst_mean, na.rm = TRUE)), by = lat_degree]
  fish[, Degree := as.integer(Degree)]
  fish[is.na(Degree), Degree := as.integer(floor(Latitude))]
  fish <- merge(fish, bl, by.x = "Degree", by.y = "lat_degree", all.x = TRUE)
  n_fb <- fish[is.na(sst_mean_lt), uniqueN(Reef)]
  fish[is.na(sst_mean_lt), sst_mean_lt := band_lt]
  if (n_fb) message("  ", n_fb, " reefs fell back to their latitude band")
}
stopifnot(!any(is.na(fish$sst_mean_lt)))
message(sprintf("Per-reef long-term SST on the panel: %.2f to %.2f C across %d reefs",
                min(fish$sst_mean_lt), max(fish$sst_mean_lt), uniqueN(fish$Reef)))

# -----------------------------------------------------------
# Kmax per unique species x temperature combination
# -----------------------------------------------------------
grid <- unique(fish[, .(species_std, MaxSizeTL, Diet, Position, Method,
                        sstmean = round(sst_mean_lt, 2))])
g <- tidytrait(as.data.frame(grid), db)
kmax <- as.data.table(predKmax(g, dataset = db,
                               fmod = formula(~ sstmean + MaxSizeTL + Diet + Position + Method),
                               niter = 100, return = "pred")$pred)
message(sprintf("Kmax predicted for %d species x temperature combinations", nrow(kmax)))
fish[, sstmean := round(sst_mean_lt, 2)]
fish <- merge(fish,
              kmax[, .(species_std, sstmean, Kmax)],
              by = c("species_std", "sstmean"), all.x = TRUE)
stopifnot(!any(is.na(fish$Kmax)))

# -----------------------------------------------------------
# Growth over one day, mass increment, per-area rates
# -----------------------------------------------------------
fish[, w := A_ord * Size^B_pen]
fish[, growth_day := somaGain(a = A_ord, b = B_pen, Lmeas = Size, t = 1,
                              Lmax = MaxSizeTL, Kmax = Kmax, silent = TRUE)]
stopifnot(all(fish$growth_day >= 0))

# Transect-level totals per 100 m2 (matching the per100 convention in 03)
agg <- fish[, .(b = sum(w * Quantity, na.rm = TRUE),
                p = sum(growth_day * Quantity, na.rm = TRUE)),
            by = trans_keys]

# Commercial split: genera present in the artisanal reef landing receipts
art_g  <- fread(file.path(OUT, "artisanal_bcs_annual.csv"))
landed <- unique(toupper(art_g[is_reef == TRUE & landings_t > 0, genus]))
landed <- landed[!is.na(landed) & landed != ""]
agg_c <- fish[toupper(Genus) %in% landed,
              .(b = sum(w * Quantity, na.rm = TRUE),
                p = sum(growth_day * Quantity, na.rm = TRUE)), by = trans_keys]
agg_n <- fish[!toupper(Genus) %in% landed,
              .(b = sum(w * Quantity, na.rm = TRUE),
                p = sum(growth_day * Quantity, na.rm = TRUE)), by = trans_keys]

zero_fill <- function(rec) {
  full <- merge(fish_core_t, rec, by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]; full[is.na(p), p := 0]
  full[, b100 := b / pmax(Area, 1) * 100]
  full[, p100 := p / pmax(Area, 1) * 100]
  full[, turnover := fifelse(b > 0, 100 * p / b, NA_real_)]   # % of biomass per day
  full[]
}
panels <- list("Whole community"        = zero_fill(agg),
               "Commercial species"     = zero_fill(agg_c),
               "Non-commercial species" = zero_fill(agg_n))

# transect-level table for the buffer analysis (03d)
fwrite(panels[["Whole community"]][, .(Year, Region, Reef, Habitat, Depth2, Transect,
                                       Area, b100, p100, turnover)],
       file.path(OUT, "productivity_transects.csv"))

med <- panels[["Whole community"]][b > 0]
message(sprintf("Panel medians: biomass %.1f g/m2, production %.1f g/m2/yr, turnover %.3f%%/day",
                med[, median(b100 / 100)], med[, median(p100 / 100 * 365)],
                med[, median(turnover, na.rm = TRUE)]))

# -----------------------------------------------------------
# Trends: same estimators as 03 (reef FE + reef-clustered SE,
# and the reef-paired early/late comparison)
# -----------------------------------------------------------
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
  tt <- t.test(log(pr$l / pr$e))
  data.table(pct_per_decade = 100 * (exp(coef(fit)["Year"] * 10) - 1),
             p_reefclust    = clust_p(fit, d$Reef, "Year"),
             paired_pct     = 100 * (exp(mean(log(pr$l / pr$e))) - 1),
             paired_lo      = 100 * (exp(tt$conf.int[1]) - 1),
             paired_hi      = 100 * (exp(tt$conf.int[2]) - 1),
             paired_p       = tt$p.value,
             n_reefs        = nrow(pr))
}

trends <- rbindlist(lapply(names(panels), function(nm) {
  d <- panels[[nm]]
  rbind(data.table(group = nm, response = "biomass",    trend_one(d, "b100")),
        data.table(group = nm, response = "production", trend_one(d, "p100")),
        data.table(group = nm, response = "turnover",   trend_one(d, "turnover", offset = 0.001)))
}))
fwrite(trends, file.path(OUT, "productivity_panel_trends.csv"))
message("\nTrends on the balanced panel (reef FE, reef-clustered p):")
print(trends[, .(group, response, pct_per_decade = round(pct_per_decade, 1),
                 p = signif(p_reefclust, 2), paired_pct = round(paired_pct, 1),
                 n_reefs)])

# -----------------------------------------------------------
# Reef-adjusted annual residual trajectories for Figure 2c,d
# (same construction as the panels a,b residuals in 03)
# -----------------------------------------------------------
resid_traj <- function(d, value, offset = 0.01) {
  d <- copy(d)[!is.na(get(value))]
  d[, l_v := log(get(value) + offset)]
  d[, r := l_v - mean(l_v, na.rm = TRUE), by = Reef]
  out <- d[, .(mean = mean(r, na.rm = TRUE),
               sem = sd(r, na.rm = TRUE) / sqrt(sum(!is.na(r))),
               count = sum(!is.na(r))), by = Year][order(Year)]
  out[, ci95 := 1.96 * sem][]
}
wc <- panels[["Whole community"]]; cm <- panels[["Commercial species"]]
fig2cd <- rbindlist(list(
  data.table(panel = "c", group = "Standing biomass",   resid_traj(wc, "b100")),
  data.table(panel = "c", group = "Biomass production", resid_traj(wc, "p100")),
  data.table(panel = "d", group = "Whole community",    resid_traj(wc, "turnover", offset = 0.001)),
  data.table(panel = "d", group = "Commercial species", resid_traj(cm, "turnover", offset = 0.001))))
fwrite(fig2cd, file.path(OUT, "productivity_fig2cd_residuals.csv"))

# -----------------------------------------------------------
# In-text summary numbers
# -----------------------------------------------------------
gv <- function(g, r, col) trends[group == g & response == r][[col]]
summ <- data.table(quantity = c(
  "panel_biomass_median_g_m2", "panel_production_median_g_m2_yr",
  "panel_turnover_median_pct_day",
  "community_biomass_pct_decade", "community_biomass_p",
  "community_production_pct_decade", "community_production_p",
  "community_turnover_pct_decade", "community_turnover_p",
  "commercial_biomass_pct_decade", "commercial_production_pct_decade",
  "commercial_turnover_pct_decade", "commercial_turnover_p",
  "noncommercial_turnover_pct_decade",
  "community_production_paired_pct", "community_biomass_paired_pct",
  "community_turnover_paired_pct",
  "elasmobranch_biomass_share_pct", "trait_coverage_individuals_pct",
  "size_cap_individuals_pct"),
  value = c(
    round(med[, median(b100 / 100)], 1), round(med[, median(p100 / 100 * 365)], 1),
    round(med[, median(turnover, na.rm = TRUE)], 3),
    round(gv("Whole community", "biomass", "pct_per_decade"), 1),
    signif(gv("Whole community", "biomass", "p_reefclust"), 2),
    round(gv("Whole community", "production", "pct_per_decade"), 1),
    signif(gv("Whole community", "production", "p_reefclust"), 2),
    round(gv("Whole community", "turnover", "pct_per_decade"), 1),
    signif(gv("Whole community", "turnover", "p_reefclust"), 2),
    round(gv("Commercial species", "biomass", "pct_per_decade"), 1),
    round(gv("Commercial species", "production", "pct_per_decade"), 1),
    round(gv("Commercial species", "turnover", "pct_per_decade"), 1),
    signif(gv("Commercial species", "turnover", "p_reefclust"), 2),
    round(gv("Non-commercial species", "turnover", "pct_per_decade"), 1),
    round(gv("Whole community", "production", "paired_pct"), 1),
    round(gv("Whole community", "biomass", "paired_pct"), 1),
    round(gv("Whole community", "turnover", "paired_pct"), 1),
    round(elasmo_share, 2), round(cov_fam, 2), round(trunc_pct, 2)))
fwrite(summ, file.path(OUT, "productivity_panel_summary.csv"))
print(summ)
message("\nStep 03b done.  Production and turnover tables written to manuscript/data/.")
