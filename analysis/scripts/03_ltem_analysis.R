# 03_ltem_analysis.R
# -----------------------------------------------------------
# LTEM reef-fixed-effect biomass analyses (2020 excluded):
#   * Reef-balanced panel (≥15 of 28 years)
#   * Fish functional-group warming response
#   * Foundation invertebrate (Asteroidea, Echinoidea,
#     Scleractinia, Holaxonia) reef-corrected trajectories
#   * Top commercial reef-fish species sensitivities
#
# Requires:  ../data/ltem.parquet   (from convert_ltem.R)
#            data/warm_season_anomaly_annual.csv
# -----------------------------------------------------------

# =====================================================================
# WARNING, UNRESOLVED (added after adversarial review). The warm-season
# anomaly is a YEAR-level variable with only 27 distinct values, but every
# model below fits it across thousands of TRANSECTS and clusters standard
# errors on REEF. Reef clustering does not fix this: the correct cluster for
# a year-level regressor is the year. Re-tested with year-clustered errors:
#
#     quantity                    naive p   reef p    YEAR p
#     year trend (the decline)    3e-15     2e-08     0.023   survives
#     community warming response  0.006     0.031     0.444   DIES
#     threshold (quadratic)       0.017     0.040     0.452   DIES
#     mean body size vs anomaly   0.003     0.020     0.486   DIES
#     per-group warming classes   various   various   all >0.13 except Wrasses
#
# At the year level (27 annual reef-adjusted indices) the warming response is
# p = 0.32 alone and p = 0.93 once a year trend is included, and the anomaly
# is correlated with year at r = 0.46, so warming cannot be separated from the
# secular trend. Using latitude-resolved SST with year fixed effects does not
# rescue it either: the coefficient flips to +0.75, which is spatial
# confounding across the four latitude degrees in the balanced panel.
#
# CONSEQUENCE: the biomass DECLINE is real and survives; its ATTRIBUTION to
# warming is not identified by these data. Any beta_anom or threshold result
# from this script must not be reported as established until this is resolved.
# =====================================================================

suppressPackageStartupMessages({
  library(arrow); library(data.table)  # models are plain lm(); lmerTest not required
})
setwd("..")
DATA <- "data"
OUT  <- "data"

ltem_path <- "../data/ltem.parquet"
if (!file.exists(ltem_path)) stop("LTEM parquet missing")
ltem <- as.data.table(read_parquet(ltem_path))
warm <- fread(file.path(DATA, "warm_season_anomaly_annual.csv"))

ltem[, Year := as.integer(Year)]
ltem <- ltem[Year != 2020]

trans_keys <- c("Year","Region","Reef","Habitat","Depth2","Transect",
                "Latitude","Degree","Area")

# Reef-balanced panel (≥15 surveyed years)
all_trans <- unique(ltem[Label=="PEC", ..trans_keys])
reef_yrs  <- all_trans[, .(yrs = uniqueN(Year)), by = Reef]
core      <- reef_yrs[yrs >= 15, Reef]
message(sprintf("Balanced reef panel: %d reefs", length(core)))

# -----------------------------------------------------------
# Fish functional-group warming model
# log(biomass + 0.01) ~ year + ws_anom + Reef
# -----------------------------------------------------------
fish <- ltem[Label == "PEC" & Reef %in% core & !is.na(Functional_groups)]
groups <- unique(fish$Functional_groups)
fg_res <- data.table()
for (g in groups) {
  rec <- fish[Functional_groups == g,
              .(b = sum(Biomass, na.rm = TRUE)),
              by = trans_keys]
  # ensure zero-fill for transects without this group
  full <- merge(all_trans[Reef %in% core], rec,
                by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]
  full[, per100 := b / pmax(Area, 1) * 100]
  full[, l_b := log(per100 + 0.01)]
  full <- merge(full, warm, by.x = "Year", by.y = "year")
  fit <- lm(l_b ~ Year + ws_anom + Reef, data = full)
  co <- summary(fit)$coefficients
  fg_res <- rbind(fg_res, data.table(
    Functional_group = g,
    n_obs    = nrow(full),
    n_reefs  = uniqueN(full$Reef),
    beta_year = co["Year","Estimate"],
    p_year    = co["Year","Pr(>|t|)"],
    beta_anom = co["ws_anom","Estimate"],
    p_anom    = co["ws_anom","Pr(>|t|)"]
  ))
}
fwrite(fg_res, file.path(OUT, "fish_functional_group_warming.csv"))

# -----------------------------------------------------------
# Foundation invertebrates: per-class reef-corrected fit
# -----------------------------------------------------------
inv <- ltem[Label == "INV" & Reef %in% core]
inv[, EchClass := fifelse(Taxa2 == "Asteroidea", "sea_star",
                  fifelse(Taxa2 == "Echinoidea", "urchin", NA_character_))]

class_specs <- list(
  "Holaxonia (gorgonians)"     = quote(Taxa3 == "Holaxonia"),
  "Scleractinia (hard corals)" = quote(Taxa3 == "Scleractinia"),
  "Asteroidea (sea stars)"     = quote(EchClass == "sea_star"),
  "Echinoidea (sea urchins)"   = quote(EchClass == "urchin")
)
all_inv_trans <- unique(inv[, ..trans_keys])
inv_res <- data.table()
for (nm in names(class_specs)) {
  rec <- inv[eval(class_specs[[nm]]),
             .(q = sum(Quantity, na.rm = TRUE)),
             by = trans_keys]
  full <- merge(all_inv_trans, rec, by = trans_keys, all.x = TRUE)
  full[is.na(q), q := 0]
  full[, per100 := q / pmax(Area, 1) * 100]
  full[, l_b := log(per100 + 0.01)]
  full <- merge(full, warm, by.x = "Year", by.y = "year")
  fit <- lm(l_b ~ Year + ws_anom + Reef, data = full)
  co <- summary(fit)$coefficients
  inv_res <- rbind(inv_res, data.table(
    group = nm, n_obs = nrow(full), n_reefs = uniqueN(full$Reef),
    beta_year = co["Year","Estimate"], p_year = co["Year","Pr(>|t|)"],
    beta_anom = co["ws_anom","Estimate"], p_anom = co["ws_anom","Pr(>|t|)"]
  ))
}
fwrite(inv_res, file.path(OUT, "invert_warming_reef_corrected.csv"))

# -----------------------------------------------------------
# Top commercial reef-fish species (Loreto + La Paz + Corredor)
# -----------------------------------------------------------
focus_sp <- c(
  "Mycteroperca rosacea","Lutjanus argentiventris",
  "Hoplopagrus guentherii","Anisotremus interruptus",
  "Lutjanus viridis","Haemulon maculicauda",
  "Cephalopholis panamensis","Lutjanus novemfasciatus",
  "Epinephelus labriformis","Paranthias colonus",
  "Haemulon sexfasciatum","Caranx caballus",
  "Microlepidotus inornatus"
)

fish_bcs <- fish[Region %in% c("Loreto","La Paz","Corredor")]
core_sub <- fish_bcs[, .(yrs = uniqueN(Year)),
                     by = Reef][yrs >= 10, Reef]
trans_id_sub <- unique(fish_bcs[Reef %in% core_sub, ..trans_keys])

sp_res <- data.table()
for (sp in focus_sp) {
  rec <- fish_bcs[Reef %in% core_sub & Species == sp,
                  .(b = sum(Biomass, na.rm = TRUE)),
                  by = trans_keys]
  full <- merge(trans_id_sub, rec, by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]
  if (nrow(full[b > 0]) < 100) next
  full[, per100 := b / pmax(Area, 1) * 100]
  full[, l_b := log(per100 + 0.01)]
  full <- merge(full, warm, by.x = "Year", by.y = "year")
  fit <- lm(l_b ~ Year + ws_anom + Reef, data = full)
  co <- summary(fit)$coefficients
  sp_res <- rbind(sp_res, data.table(
    Species  = sp,
    n_obs    = nrow(full),
    n_reefs  = uniqueN(full$Reef),
    beta_year = co["Year","Estimate"],
    p_year    = co["Year","Pr(>|t|)"],
    beta_anom = co["ws_anom","Estimate"],
    p_anom    = co["ws_anom","Pr(>|t|)"]
  ))
}
fwrite(sp_res[order(beta_anom)], file.path(OUT, "top_commercial_species_warming.csv"))

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

# -----------------------------------------------------------
# Reef-paired change: each reef against itself, early vs late
# -----------------------------------------------------------
# This is the loss estimate the manuscript leads on, because it does not
# depend on the log offset (the fitted log trend swings between -17 and -46
# percent depending on that choice) and it does not depend on the CONAPESCA
# landing statistics at all. Every reef is compared with itself, so changes
# in which reefs were surveyed cannot drive it.
rp_rec  <- ltem[Label == "PEC" & Reef %in% core, .(b = sum(Biomass, na.rm = TRUE)),
                by = trans_keys]
rp_full <- merge(all_trans[Reef %in% core], rp_rec, by = trans_keys, all.x = TRUE)
rp_full[is.na(b), b := 0]
rp_full[, per100 := b / pmax(Area, 1) * 100]
paired <- rp_full[, .(early = mean(per100[Year <= 2013]),
                      late  = mean(per100[Year >= 2014])), by = Reef
                  ][!is.na(early) & !is.na(late) & early > 0 & late > 0]
paired[, log_ratio := log(late / early)]
paired[, pct_change := 100 * (late / early - 1)]
tt <- t.test(paired$log_ratio)
fwrite(paired[order(pct_change)], file.path(OUT, "reef_paired_change.csv"))
rp_summary <- data.table(
  n_reefs   = nrow(paired),
  n_down    = paired[log_ratio < 0, .N],
  pct_mean  = 100 * (exp(mean(paired$log_ratio)) - 1),
  ci_lo     = 100 * (exp(tt$conf.int[1]) - 1),
  ci_hi     = 100 * (exp(tt$conf.int[2]) - 1),
  p_value   = tt$p.value)
fwrite(rp_summary, file.path(OUT, "reef_paired_summary.csv"))
message(sprintf("Reef-paired change: %+.1f%% [%.1f to %.1f], p=%.4f, %d of %d reefs down",
        rp_summary$pct_mean, rp_summary$ci_lo, rp_summary$ci_hi,
        rp_summary$p_value, rp_summary$n_down, rp_summary$n_reefs))

# -----------------------------------------------------------
# Warming response per REEF GROUP (the groups used by the fishery and
# economic steps). Fitting the same model at the group level is what
# lets Figure 3b colour each commercial group by a measured LTEM
# response instead of an asserted one.
# -----------------------------------------------------------
rg <- fread(file.path(OUT, "reef_group_genera.csv"))
grp_res <- data.table()
for (g in unique(rg$reef_group)) {
  gen <- rg[reef_group == g, genus]
  rec <- fish_bcs[Reef %in% core_sub & toupper(Genus) %in% gen,
                  .(b = sum(Biomass, na.rm = TRUE)), by = trans_keys]
  full <- merge(trans_id_sub, rec, by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]
  if (nrow(full[b > 0]) < 100) next
  full[, per100 := b / pmax(Area, 1) * 100]
  full[, l_b := log(per100 + 0.01)]
  full <- merge(full, warm, by.x = "Year", by.y = "year")
  fit <- lm(l_b ~ Year + ws_anom + Reef, data = full)
  co <- summary(fit)$coefficients
  grp_res <- rbind(grp_res, data.table(
    reef_group = g, n_obs = nrow(full), n_reefs = uniqueN(full$Reef),
    beta_year = co["Year","Estimate"],   p_year = co["Year","Pr(>|t|)"],
    beta_anom = co["ws_anom","Estimate"], p_anom = co["ws_anom","Pr(>|t|)"]))
}
fwrite(grp_res[order(beta_anom)], file.path(OUT, "reef_group_warming.csv"))

# Per reef group, the TIME TREND with reef-clustered errors. This is what
# Figure 3 uses to colour each commercial group, because the trend is
# identified at the reef level and survives clustering, whereas the warming
# response is a year-level regressor that does not (see the warning above).
grp_trend <- data.table()
for (g in unique(rg$reef_group)) {
  gen <- rg[reef_group == g, genus]
  rec <- fish_bcs[Reef %in% core_sub & toupper(Genus) %in% gen,
                  .(b = sum(Biomass, na.rm = TRUE)), by = trans_keys]
  full <- merge(trans_id_sub, rec, by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]
  if (nrow(full[b > 0]) < 100) next
  full[, per100 := b / pmax(Area, 1) * 100]
  full[, l_b := log(per100 + 0.01)]
  fit <- lm(l_b ~ Year + Reef, data = full)
  pr <- full[, .(e = mean(per100[Year <= 2013]), l = mean(per100[Year >= 2014])), by = Reef
             ][!is.na(e) & !is.na(l) & e > 0 & l > 0]
  grp_trend <- rbind(grp_trend, data.table(
    reef_group = g,
    beta_year  = coef(fit)["Year"],
    p_year     = clust_p(fit, full$Reef, "Year"),
    pct_per_decade = 100 * (exp(coef(fit)["Year"] * 10) - 1),
    paired_pct = if (nrow(pr) > 3) 100 * (exp(mean(log(pr$l / pr$e))) - 1) else NA_real_,
    n_reefs    = nrow(pr)))
}
grp_trend[, trend_class := fcase(
  p_year < 0.05 & beta_year < 0, "decline_sig",
  p_year < 0.05 & beta_year > 0, "increase_sig",
  default = "no_sig_trend")]
fwrite(grp_trend[order(pct_per_decade)], file.path(OUT, "reef_group_trend.csv"))
message("Reef-group time trends (reef-clustered):")
print(grp_trend[order(pct_per_decade), .(reef_group, pct_per_decade = round(pct_per_decade, 1),
                                          p_year = signif(p_year, 3), trend_class)])
message("Reef-group warming response:"); print(grp_res[order(beta_anom)])

# -----------------------------------------------------------
# Per-year mean +/- 95% CI trajectories for Figure 2
#   community + Mycteroperca: balanced fish panel (`core`), per 100 m^2
#   sea stars + urchins:      all invertebrate reefs, per 100 m^2
# (ci95 = 1.96 * sem)
# -----------------------------------------------------------
traj <- function(per100_dt, with_median = FALSE) {
  out <- per100_dt[, .(mean = mean(per100, na.rm = TRUE), median = median(per100, na.rm = TRUE),
                       sem = sd(per100, na.rm = TRUE)/sqrt(sum(!is.na(per100))),
                       count = sum(!is.na(per100))), by = Year][order(Year)]
  out[, ci95 := 1.96 * sem]
  if (with_median) setcolorder(out, c("Year","mean","median","sem","count","ci95"))
  else { out[, median := NULL]; setcolorder(out, c("Year","mean","sem","count","ci95")) }
  out[]
}
per100_panel <- function(dt, all_t, vcol, filt = NULL) {
  sub <- if (is.null(filt)) dt else dt[eval(filt)]
  rec <- sub[, .(v = sum(get(vcol), na.rm = TRUE)), by = trans_keys]
  full <- merge(all_t, rec, by = trans_keys, all.x = TRUE)
  full[is.na(v), v := 0]
  full[, per100 := v / pmax(Area, 1) * 100][]
}

# Whole community, commercial and non-commercial fish on the balanced panel
# (2020 already dropped). "Commercial" means the genus actually appears in the
# artisanal reef landing receipts, so the split is made by what the fishery
# takes rather than by an assumed target list. Splitting this way lets the
# reef panel be read against the landings directly: the non-commercial fish
# are the internal control that fishing does not remove.
fish_core_t <- all_trans[Reef %in% core]
comm <- per100_panel(ltem[Label == "PEC" & Reef %in% core], fish_core_t, "Biomass")
fwrite(traj(comm, with_median = TRUE), file.path(OUT, "ltem_balanced_community_biomass.csv"))

art_g  <- fread(file.path(OUT, "artisanal_bcs_annual.csv"))
landed <- unique(toupper(art_g[is_reef == TRUE & landings_t > 0, genus]))
landed <- landed[!is.na(landed) & landed != ""]
message(sprintf("commercial reef genera (present in the receipts): %d", length(landed)))

com_f <- per100_panel(ltem[Label == "PEC" & Reef %in% core], fish_core_t, "Biomass",
                      quote(toupper(Genus) %in% landed))
non_f <- per100_panel(ltem[Label == "PEC" & Reef %in% core], fish_core_t, "Biomass",
                      quote(!toupper(Genus) %in% landed))
fwrite(traj(com_f), file.path(OUT, "ltem_balanced_commercial_biomass.csv"))
fwrite(traj(non_f), file.path(OUT, "ltem_balanced_noncommercial_biomass.csv"))

# sea stars + urchins on all invertebrate reefs (per-transect Quantity per 100 m^2)
inv_all   <- ltem[Label == "INV"]
inv_all_t <- unique(inv_all[, ..trans_keys])
ss <- per100_panel(inv_all, inv_all_t, "Quantity", quote(Taxa2 == "Asteroidea"))
fwrite(traj(ss), file.path(OUT, "ltem_seastars_annual.csv"))
ur <- per100_panel(inv_all, inv_all_t, "Quantity", quote(Taxa2 == "Echinoidea"))
fwrite(traj(ur), file.path(OUT, "ltem_urchins_annual.csv"))
co <- per100_panel(inv_all, inv_all_t, "Quantity", quote(Taxa3 == "Scleractinia"))
fwrite(traj(co), file.path(OUT, "ltem_corals_annual.csv"))
go <- per100_panel(inv_all, inv_all_t, "Quantity", quote(Taxa3 == "Holaxonia"))
fwrite(traj(go), file.path(OUT, "ltem_gorgonians_annual.csv"))

# Reef-adjusted residual trajectories for Figure 2 (deviation from reef mean
# log-density / log-biomass; the r_ijt of the Supplementary Methods), on the
# balanced core panel. ci95 = 1.96 * sem.
resid_traj <- function(p100) {
  d <- copy(p100); d[, l_b := log(per100 + 0.01)]
  d[, r := l_b - mean(l_b, na.rm = TRUE), by = Reef]
  out <- d[, .(mean = mean(r, na.rm = TRUE), sem = sd(r, na.rm = TRUE)/sqrt(sum(!is.na(r))),
               count = sum(!is.na(r))), by = Year][order(Year)]
  out[, ci95 := 1.96 * sem][]
}
inv_core   <- ltem[Label == "INV" & Reef %in% core]
inv_core_t <- unique(inv_core[, ..trans_keys])
ssc <- per100_panel(inv_core, inv_core_t, "Quantity", quote(Taxa2 == "Asteroidea"))
urc <- per100_panel(inv_core, inv_core_t, "Quantity", quote(Taxa2 == "Echinoidea"))
coc <- per100_panel(inv_core, inv_core_t, "Quantity", quote(Taxa3 == "Scleractinia"))
goc <- per100_panel(inv_core, inv_core_t, "Quantity", quote(Taxa3 == "Holaxonia"))
fig2 <- rbindlist(list(
  data.table(panel = "a", group = "Sea stars (Asteroidea)",   resid_traj(ssc)),
  data.table(panel = "a", group = "Sea urchins (Echinoidea)", resid_traj(urc)),
  data.table(panel = "a", group = "Hard corals (Scleractinia)", resid_traj(coc)),
  data.table(panel = "a", group = "Gorgonians (Holaxonia)",   resid_traj(goc)),
  data.table(panel = "b", group = "Whole reef-fish community",  resid_traj(comm)),
  data.table(panel = "b", group = "Commercial species",         resid_traj(com_f)),
  data.table(panel = "b", group = "Non-commercial species",     resid_traj(non_f))))
fwrite(fig2, file.path(OUT, "ltem_fig2_residuals.csv"))

# Reef-fixed-effect trends for the whole community and the snappers (Fig 2b).
# p_year is clustered on reef, which is the correct unit for a time trend
# fitted across repeatedly surveyed reefs. p_anom is reported with the same
# clustering but see the warning at the top of this file: the anomaly is a
# year-level regressor and does not survive year-level clustering, so it must
# not be reported as an established result.
fitFE <- function(p100) {
  d <- merge(p100, warm, by.x = "Year", by.y = "year")
  d[, l_b := log(per100 + 0.01)]
  fit <- lm(l_b ~ Year + ws_anom + Reef, data = d)
  cf  <- coef(fit)
  data.table(beta_year = cf["Year"],    p_year_reefclust = clust_p(fit, d$Reef, "Year"),
             pct_per_decade = 100 * (exp(cf["Year"] * 10) - 1),
             beta_anom = cf["ws_anom"],
             p_anom_reefclust = clust_p(fit, d$Reef, "ws_anom"),
             p_anom_yearclust = clust_p(fit, d$Year, "ws_anom"))
}
fish_fe <- rbind(data.table(group = "Whole reef-fish community", fitFE(comm)),
                 data.table(group = "Commercial species",        fitFE(com_f)),
                 data.table(group = "Non-commercial species",    fitFE(non_f)))
fwrite(fish_fe, file.path(OUT, "reef_fish_community_warming.csv"))
print(fish_fe)

message("\nStep 03 done.  Intermediates + Figure-2 trajectories written to manuscript/data/.")
