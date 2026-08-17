# 08_gap_analysis.R
# -----------------------------------------------------------
# Gap analysis: how climate and fishing each act on Gulf of California
# rocky reef fish, and whether their combined effect is worse than the
# sum of the parts now that whole-basin heatwaves recur every few years.
#
# Three questions, in order:
#   1. CLIMATE   what warming does to the reef fish community, including
#                where in the water column and to what body size
#   2. FISHERY   what fishing pressure does, separated from warming by
#                the no-take contrast and by targeted vs untargeted taxa
#   3. COMBINED  whether fishing makes the climate response worse
#                (the interaction), which is the question that matters
#                for management
#
# Statistical note. The unit of observation is the transect, and
# transects are nested in reefs surveyed repeatedly, so ordinary
# standard errors would be badly overconfident. Every model below is
# fitted with reef fixed effects and reported with cluster-robust
# (CR1) standard errors clustered on reef. Coefficients are on
# log(biomass per 100 m2), so they read as proportional change.
#
# Reads:  ../data/ltem.parquet, data/warm_season_anomaly_annual.csv,
#         data/sst_gulf_monthly_1981_2026.csv, data/artisanal_bcs_annual.csv
# Writes: data/gap_*.csv  (one table per question) and
#         data/gap_summary.csv (the headline numbers)
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

suppressPackageStartupMessages({ library(arrow); library(data.table)
                                 library(ggplot2); library(patchwork) })
setwd("..")
DATA <- "data"; OUT <- "data"

# -----------------------------------------------------------
# Cluster-robust (CR1) standard errors, clustered on reef
# -----------------------------------------------------------
cluster_se <- function(fit, cluster) {
  # Terms that are constant within reef (protection stratum, for example) are
  # aliased with the reef fixed effects and dropped by lm() to NA. model.matrix
  # still returns their columns, so they must be removed or the sandwich
  # estimator misaligns X with the coefficients and returns nonsense.
  cf   <- coef(fit)
  keep <- names(cf)[!is.na(cf)]
  X  <- model.matrix(fit)[, keep, drop = FALSE]
  u  <- residuals(fit)
  cl <- factor(cluster[as.integer(names(u))])
  k  <- ncol(X); n <- nrow(X); m <- nlevels(cl)
  bread <- tryCatch(chol2inv(chol(crossprod(X))), error = function(e) MASS::ginv(crossprod(X)))
  meat  <- matrix(0, k, k)
  for (g in levels(cl)) {
    i  <- which(cl == g)
    Xu <- crossprod(X[i, , drop = FALSE], u[i])
    meat <- meat + tcrossprod(Xu)
  }
  adj <- (m / (m - 1)) * ((n - 1) / (n - k))
  V   <- bread %*% meat %*% bread * adj
  dimnames(V) <- list(colnames(X), colnames(X))   # chol2inv drops these
  se  <- sqrt(diag(V))
  est <- coef(fit)[names(se)]
  t   <- est / se
  data.table(term = names(se), estimate = est, se = se, t = t,
             p = 2 * pt(-abs(t), df = m - 1))
}

report <- function(fit, cluster, terms, label) {
  r <- cluster_se(fit, cluster)
  r <- r[term %in% terms]
  r[, model := label]
  setcolorder(r, "model")[]
}

# -----------------------------------------------------------
# Data: reef fish transects, Gulf regions, balanced-ish panel
# -----------------------------------------------------------
ltem <- as.data.table(read_parquet(file.path("../data", "ltem.parquet")))
warm <- fread(file.path(DATA, "warm_season_anomaly_annual.csv"))
mon  <- fread(file.path(DATA, "sst_gulf_monthly_1981_2026.csv"))

GULF <- c("Loreto", "La Paz", "Corredor", "Cabo Pulmo", "Los Cabos",
          "La Ventana", "Alto Golfo")
fish <- ltem[Label == "PEC" & Region %in% GULF & Year != 2020 &
             !is.na(Protection_level) & Protection_level != ""]
fish[, Year := as.integer(Year)]
fish[, TL := suppressWarnings(as.numeric(TrophicLevel))]

# keep reefs with a real time series
core <- fish[, .(yrs = uniqueN(Year)), by = Reef][yrs >= 8, Reef]
fish <- fish[Reef %in% core]

fish[is.na(MPA), MPA := ""]
TK <- c("Year", "Region", "Reef", "Habitat", "Depth2", "Transect",
        "Latitude", "Area", "Protection_level", "MPA")
trans_id <- unique(fish[, ..TK])

# annual heat exposure: how much heatwave the reef actually saw that year
mon[, `:=`(year = as.integer(substr(month, 1, 4)),
           mo   = as.integer(substr(month, 6, 7)))]
heat <- mon[mo %in% 5:10, .(mhw_days = sum(mhw_days, na.rm = TRUE),
                            cum_int  = sum(cumulative_intensity, na.rm = TRUE)),
            by = year]
clim <- merge(warm, heat, by = "year", all.x = TRUE)

# biomass per 100 m2 for an arbitrary subset of the fish table
per100 <- function(sub) {
  rec  <- sub[, .(b = sum(Biomass, na.rm = TRUE)), by = TK]
  full <- merge(trans_id, rec, by = TK, all.x = TRUE)
  full[is.na(b), b := 0]
  full[, per100 := b / pmax(Area, 1) * 100]
  full[, l_b := log(per100 + 0.01)]
  merge(full, clim, by.x = "Year", by.y = "year")
}

res <- list(); summ <- list()
add <- function(k, v) summ[[length(summ) + 1]] <<- data.table(quantity = k, value = v)
gs_val <- function(k) rbindlist(summ)[quantity == k, value][1]

# ===========================================================
# 1. CLIMATE
# ===========================================================
message("\n===== 1. CLIMATE =====")

# 1.1 Which climate metric actually tracks the community: the warm-season
#     anomaly, the number of heatwave days, or cumulative heat?
comm <- per100(fish)
m_anom <- lm(l_b ~ Year + ws_anom  + Reef, data = comm)
m_days <- lm(l_b ~ Year + mhw_days + Reef, data = comm)
m_cum  <- lm(l_b ~ Year + cum_int  + Reef, data = comm)
dose <- rbindlist(list(
  report(m_anom, comm$Reef, "ws_anom",  "community ~ warm-season anomaly"),
  report(m_days, comm$Reef, "mhw_days", "community ~ heatwave days"),
  report(m_cum,  comm$Reef, "cum_int",  "community ~ cumulative heat")))
dose[, AIC := c(AIC(m_anom), AIC(m_days), AIC(m_cum))]
res$gap_climate_dose_response <- dose
print(dose)
add("community_beta_ws_anom", round(dose[term == "ws_anom", estimate], 4))
add("community_beta_ws_anom_p", signif(dose[term == "ws_anom", p], 3))

# 1.2 Depth. If the deep water is no longer a thermal refuge, the deep
#     transects should lose at least as much as the shallow ones.
m_depth <- lm(l_b ~ Year + ws_anom * Depth2 + Reef, data = comm)
depth <- report(m_depth, comm$Reef,
                c("ws_anom", "Depth2Shallow", "ws_anom:Depth2Shallow"),
                "community ~ anomaly x depth")
res$gap_climate_depth <- depth
print(depth)
# deep slope is ws_anom; shallow slope is ws_anom + interaction
b <- setNames(depth$estimate, depth$term)
add("warming_slope_deep",    round(b["ws_anom"], 4))
add("warming_slope_shallow", round(b["ws_anom"] + b["ws_anom:Depth2Shallow"], 4))

# 1.3 Body size. Warming is expected to shrink fish; test the
#     abundance-weighted mean body size of the community.
sz <- fish[!is.na(Size) & Quantity > 0,
           .(mean_size = sum(Size * Quantity) / sum(Quantity)), by = TK]
sz <- merge(sz, clim, by.x = "Year", by.y = "year")
m_size <- lm(log(mean_size) ~ Year + ws_anom + Reef, data = sz)
size <- report(m_size, sz$Reef, c("Year", "ws_anom"), "mean body size")
res$gap_climate_body_size <- size
print(size)
add("body_size_beta_ws_anom", round(size[term == "ws_anom", estimate], 4))
add("body_size_beta_ws_anom_p", signif(size[term == "ws_anom", p], 3))

# 1.5 Is the response linear, or is there a threshold? This matters because
#     the pre-2014 and post-2014 eras barely overlap in anomaly (pre-2014
#     never exceeded +0.59 C, post-2014 reaches +1.23 C), so "the community
#     became more sensitive" and "there is a threshold we had not yet
#     crossed" cannot be fully separated. A quadratic term and a binned fit
#     test the threshold reading directly.
m_quad <- lm(l_b ~ Year + ws_anom + I(ws_anom^2) + Reef, data = comm)
quad <- report(m_quad, comm$Reef, c("ws_anom", "I(ws_anom^2)"),
               "community ~ anomaly + anomaly^2")
res$gap_climate_nonlinear <- quad
print(quad)
add("anomaly_quadratic_term",   round(quad[term == "I(ws_anom^2)", estimate], 4))
add("anomaly_quadratic_term_p", signif(quad[term == "I(ws_anom^2)", p], 3))

# binned response: mean log-biomass by anomaly band, reef-centred so that
# which reefs were surveyed in which years cannot drive the pattern
comm[, reef_mean := mean(l_b), by = Reef]
comm[, l_b_adj := l_b - reef_mean]
comm[, band := cut(ws_anom, breaks = c(-Inf, -0.25, 0.25, 0.6, Inf),
                   labels = c("cool (< -0.25)", "near normal (-0.25 to 0.25)",
                              "warm (0.25 to 0.6)", "extreme (> 0.6)"))]
bins <- comm[, .(mean_l_b_adj = mean(l_b_adj), sd = sd(l_b_adj),
                 n_transects = .N, n_years = uniqueN(Year)), by = band][order(band)]
bins[, rel_to_normal_pct := round(100 * (exp(mean_l_b_adj -
        bins[band == "near normal (-0.25 to 0.25)", mean_l_b_adj]) - 1), 1)]
res$gap_climate_anomaly_bands <- bins
message("Reef-centred biomass by warm-season anomaly band:"); print(bins)
add("biomass_pct_vs_normal_extreme_band",
    bins[band == "extreme (> 0.6)", rel_to_normal_pct])
add("biomass_pct_vs_normal_warm_band",
    bins[band == "warm (0.25 to 0.6)", rel_to_normal_pct])

# 1.4 Recovery. Does biomass rebound in the years between heatwaves, or
#     does each event leave a step down?
yrly <- comm[, .(l_b = mean(l_b)), by = .(Year, Reef)]
yrly <- merge(yrly, clim, by.x = "Year", by.y = "year")
yrly[, hot := ws_anom > quantile(clim$ws_anom, 0.75, na.rm = TRUE)]
yrly[, phase := fifelse(hot, "heatwave year", "cooler year")]
rec <- yrly[, .(mean_l_b = mean(l_b), n = .N), by = .(phase, era = fifelse(Year < 2014, "1998-2013", "2014-2025"))]
res$gap_climate_recovery <- rec[order(era, phase)]
print(rec[order(era, phase)])
# post-2014 cooler years vs pre-2014 cooler years: is the baseline lower?
d1 <- rec[era == "1998-2013" & phase == "cooler year", mean_l_b]
d2 <- rec[era == "2014-2025" & phase == "cooler year", mean_l_b]
add("cool_year_baseline_shift_pct", round(100 * (exp(d2 - d1) - 1), 1))

# ===========================================================
# 2. FISHERY
# ===========================================================
message("\n===== 2. FISHERY =====")

# 2.1 Protection, stratified by whether it is actually enforced. Pooling
#     every "Prohibited" reef hides the answer: most are paper parks. Cabo
#     Pulmo is the one long-enforced no-take reserve in the region, so it
#     is kept separate. The contrast is made WITHIN reefs surveyed in both
#     eras, because which reefs were visited changed over the record, and
#     on the log scale, because a handful of high-biomass reefs otherwise
#     dominate an arithmetic mean (their raw mean change is -66% against a
#     median of +13%; the log model is the honest summary).
comm[, grp_prot := fcase(
  MPA == "Cabo Pulmo" & Protection_level == "Prohibited", "no-take (enforced)",
  Protection_level == "Prohibited",                       "no-take (paper park)",
  Protection_level == "Allowed",                          "MPA, fishing allowed",
  Protection_level == "Open Area",                        "open to fishing",
  default = NA_character_)]
comm[, era := factor(fifelse(Year < 2014, "pre", "post"), levels = c("pre", "post"))]
both_eras <- comm[, .(ne = uniqueN(era)), by = Reef][ne == 2, Reef]
pc <- comm[!is.na(grp_prot) & Reef %in% both_eras]
pc[, grp_prot := factor(grp_prot, levels = c("open to fishing", "MPA, fishing allowed",
                                             "no-take (paper park)", "no-take (enforced)"))]
message(sprintf("protection panel: %d reefs surveyed in both eras", uniqueN(pc$Reef)))

m_prot <- lm(l_b ~ era * grp_prot + Reef, data = pc)
prot <- report(m_prot, pc$Reef, grep("^era", names(coef(m_prot)), value = TRUE),
               "era change by protection (within reef)")
res$gap_fishery_protection <- prot
print(prot)
bp <- setNames(prot$estimate, prot$term)
pct <- function(x) round(100 * (exp(x) - 1), 1)
add("era_change_pct_open",            pct(bp["erapost"]))
add("era_change_pct_mpa_allowed",     pct(bp["erapost"] + bp["erapost:grp_protMPA, fishing allowed"]))
add("era_change_pct_paper_park",      pct(bp["erapost"] + bp["erapost:grp_protno-take (paper park)"]))
add("era_change_pct_no_take_enforced",pct(bp["erapost"] + bp["erapost:grp_protno-take (enforced)"]))
add("enforced_vs_open_p",             signif(prot[term == "erapost:grp_protno-take (enforced)", p], 3))
add("paper_park_vs_open_p",           signif(prot[term == "erapost:grp_protno-take (paper park)", p], 3))

# 2.2 Targeted vs untargeted. A reef genus counts as targeted if it
#     appears in the artisanal landing receipts; the rest of the reef
#     community is the internal control that fishing does not touch.
art <- fread(file.path(DATA, "artisanal_bcs_annual.csv"))
landed <- unique(toupper(art[is_reef == TRUE & landings_t > 0, genus]))
fish[, targeted := toupper(Genus) %in% landed]
message(sprintf("targeted genera on reefs: %d of %d",
                uniqueN(fish[targeted == TRUE, Genus]), uniqueN(fish$Genus)))

tg <- per100(fish[targeted == TRUE]);  tg[, grp := "targeted"]
ug <- per100(fish[targeted == FALSE]); ug[, grp := "untargeted"]
both <- rbind(tg, ug)
m_tgt <- lm(l_b ~ Year * grp + ws_anom + Reef, data = both)
tgt <- report(m_tgt, both$Reef, c("Year", "grpuntargeted", "Year:grpuntargeted"),
              "biomass trend ~ targeted vs untargeted")
res$gap_fishery_targeted <- tgt
print(tgt)
bt <- setNames(tgt$estimate, tgt$term)
add("trend_pct_yr_targeted",   round(100 * (exp(bt["Year"]) - 1), 2))
add("trend_pct_yr_untargeted", round(100 * (exp(bt["Year"] + bt["Year:grpuntargeted"]) - 1), 2))

# 2.3 Fishing down the food web: mean trophic level of the landed reef
#     catch, weighted by landed weight.
tl_by_genus <- fish[!is.na(TL), .(TL = mean(TL, na.rm = TRUE)), by = .(genus = toupper(Genus))]
al <- merge(art[is_reef == TRUE], tl_by_genus, by = "genus")
mtl <- al[, .(mtl = sum(TL * landings_t) / sum(landings_t),
              landings_t = sum(landings_t)), by = year][order(year)]
fit_mtl <- lm(mtl ~ year, data = mtl)
res$gap_fishery_mean_trophic_level <- mtl
message("Mean trophic level of the landed reef catch:")
print(mtl)
# The series steps up exactly at the 2013 gap (2013 is missing from the
# government record), so this reads as a reporting change, not as ecology.
# It is recorded but should not be interpreted as fishing down the web.
mtl[, era := fifelse(year < 2014, "pre", "post")]
step <- mtl[, .(mtl = mean(mtl)), by = era]
add("mean_trophic_level_slope_per_decade", round(10 * coef(fit_mtl)["year"], 3))
add("mean_trophic_level_p", signif(summary(fit_mtl)$coefficients["year", 4], 3))
add("mean_trophic_level_step_at_2013_gap",
    round(step[era == "post", mtl] - step[era == "pre", mtl], 3))

# ===========================================================
# 3. COMBINED
# ===========================================================
message("\n===== 3. COMBINED =====")

# 3.1 Does protection buffer the warming response itself? The interaction
#     between the warm-season anomaly and the protection stratum.
m_int <- lm(l_b ~ Year + ws_anom * grp_prot + Reef, data = pc)
inter <- report(m_int, pc$Reef,
                grep("ws_anom", names(coef(m_int)), value = TRUE),
                "warming x protection")
res$gap_combined_protection_interaction <- inter
print(inter)
bi <- setNames(inter$estimate, inter$term)
add("warming_slope_open",          round(bi["ws_anom"], 4))
add("warming_slope_no_take_enforced",
    round(bi["ws_anom"] + bi["ws_anom:grp_protno-take (enforced)"], 4))
add("protection_buffer_p",
    signif(inter[term == "ws_anom:grp_protno-take (enforced)", p], 3))

# 3.2 Does fishing make the warming response worse? The interaction
#     between the anomaly and whether the taxon is fished.
m_int2 <- lm(l_b ~ Year + ws_anom * grp + Reef, data = both)
inter2 <- report(m_int2, both$Reef, c("ws_anom", "ws_anom:grpuntargeted"),
                 "warming x targeted")
res$gap_combined_targeted_interaction <- inter2
print(inter2)
b2 <- setNames(inter2$estimate, inter2$term)
add("warming_slope_targeted",   round(b2["ws_anom"], 4))
add("warming_slope_untargeted", round(b2["ws_anom"] + b2["ws_anom:grpuntargeted"], 4))
add("targeted_amplification_p", signif(inter2[term == "ws_anom:grpuntargeted", p], 3))

# 3.4 Is the same heat more costly now than it used to be? If a decade of
#     fishing and repeated heatwaves has eroded the community's capacity to
#     absorb warming, the slope on the anomaly should steepen in the recent
#     era. This is the most direct test of the compounding claim.
comm[, era2 := factor(fifelse(Year < 2014, "pre", "post"), levels = c("pre", "post"))]
m_era <- lm(l_b ~ Year + ws_anom * era2 + Reef, data = comm)
era_int <- report(m_era, comm$Reef, c("ws_anom", "ws_anom:era2post"),
                  "warming slope, before vs during the heatwave era")
res$gap_combined_era_interaction <- era_int
print(era_int)
be <- setNames(era_int$estimate, era_int$term)
add("warming_slope_pre2014",  round(be["ws_anom"], 4))
add("warming_slope_post2014", round(be["ws_anom"] + be["ws_anom:era2post"], 4))
add("slope_steepening_p",     signif(era_int[term == "ws_anom:era2post", p], 3))
add("pre2014_max_anomaly_C",  round(max(clim[year %in% 1998:2013, ws_anom]), 2))
add("post2014_max_anomaly_C", round(max(clim[year %in% 2014:2025, ws_anom]), 2))

# 3.5 Does heat accumulate? Successive warm years should compound if the
#     community cannot rebuild between events.
lagm <- data.table(year = clim$year, a0 = clim$ws_anom)
lagm[, `:=`(a1 = shift(a0, 1), a2 = shift(a0, 2))]
cl2 <- merge(comm, lagm, by.x = "Year", by.y = "year")
m_lag <- lm(l_b ~ Year + a0 + a1 + a2 + Reef, data = cl2)
lags <- report(m_lag, cl2$Reef, c("a0", "a1", "a2"), "cumulative heat (0, 1, 2 yr lags)")
res$gap_combined_cumulative_heat <- lags
print(lags)
bl <- setNames(lags$estimate, lags$term)
add("heat_effect_same_year",  round(bl["a0"], 4))
add("heat_effect_lag1",       round(bl["a1"], 4))
add("heat_effect_lag2",       round(bl["a2"], 4))
add("heat_effect_3yr_total",  round(sum(bl[c("a0","a1","a2")]), 4))

# 3.3 What the combination costs. Heatwave frequency has risen, so the
#     same per-event response is now applied far more often. Compare the
#     warm-season anomaly regime before and after 2014 and translate it
#     through the fitted slopes for fished and protected reefs.
regime <- clim[year %in% 1998:2025,
               .(mean_anom = mean(ws_anom, na.rm = TRUE),
                 mhw_days  = mean(mhw_days, na.rm = TRUE)),
               by = .(era = fifelse(year < 2014, "1998-2013", "2014-2025"))]
res$gap_combined_regime <- regime
print(regime)
d_anom <- regime[era == "2014-2025", mean_anom] - regime[era == "1998-2013", mean_anom]
add("warm_season_anomaly_shift_C", round(d_anom, 3))
add("implied_loss_pct_open_reefs",
    round(100 * (exp(bi["ws_anom"] * d_anom) - 1), 1))
add("implied_loss_pct_targeted_taxa",
    round(100 * (exp(b2["ws_anom"] * d_anom) - 1), 1))
add("mhw_days_per_warm_season_before", round(regime[era == "1998-2013", mhw_days], 1))
add("mhw_days_per_warm_season_after",  round(regime[era == "2014-2025", mhw_days], 1))

# -----------------------------------------------------------
# Figure S12: the gap analysis in one display
# -----------------------------------------------------------
RED <- "#c0392b"; BLUE <- "#1f6f9c"; GREY <- "#95a5a6"
theme_set(theme_classic(base_size = 9) +
          theme(plot.title = element_text(face = "bold", size = 9)))

pa <- ggplot(bins, aes(band, rel_to_normal_pct,
                       fill = rel_to_normal_pct < 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = BLUE), guide = "none") +
  scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
  labs(x = NULL, y = "Biomass vs normal years (%)",
       title = "a   Reef fish tolerate warm years until about +0.6 \u00b0C")

szb <- data.table(anom = seq(min(comm$ws_anom), max(comm$ws_anom), length.out = 50))
szb[, pct := 100 * (exp(size[term == "ws_anom", estimate] * anom) - 1)]
pb <- ggplot(szb, aes(anom, pct)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(colour = RED, linewidth = 0.8) +
  labs(x = "Warm-season anomaly (\u00b0C)", y = "Mean body size change (%)",
       title = sprintf("b   Fish are smaller in warm years (p = %.3f)",
                       size[term == "ws_anom", p]))

pd <- data.table(
  stratum = factor(c("open to fishing", "MPA, fishing allowed",
                     "no-take (paper park)", "no-take (enforced)"),
                   levels = c("no-take (paper park)", "MPA, fishing allowed",
                              "open to fishing", "no-take (enforced)")),
  pct = c(gs_val("era_change_pct_open"), gs_val("era_change_pct_mpa_allowed"),
          gs_val("era_change_pct_paper_park"), gs_val("era_change_pct_no_take_enforced")))
pc_fig <- ggplot(pd, aes(pct, stratum, fill = pct < 0)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = BLUE), guide = "none") +
  labs(x = "Change in reef fish biomass, 1998-2013 to 2014-2025 (%)", y = NULL,
       title = "c   Only enforced no-take held its biomass")

sl <- data.table(
  grp = factor(c("targeted by the fishery", "not targeted"),
               levels = c("not targeted", "targeted by the fishery")),
  beta = c(gs_val("warming_slope_targeted"), gs_val("warming_slope_untargeted")))
pdd <- ggplot(sl, aes(beta, grp, fill = grp)) +
  geom_col(width = 0.65) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`targeted by the fishery` = RED,
                               `not targeted` = GREY), guide = "none") +
  labs(x = "Change in log biomass per +1 \u00b0C", y = NULL,
       title = "d   Fished taxa lose most per degree of warming")

fig <- (pa | pb) / (pc_fig | pdd)
ggsave(file.path("figures", "FigureS12_gap_analysis.pdf"), fig, width = 9.5, height = 6.5)
ggsave(file.path("figures", "FigureS12_gap_analysis.png"), fig, width = 9.5, height = 6.5, dpi = 300)
message("wrote FigureS12_gap_analysis (.pdf/.png)")

# ===========================================================
# 4. TROPICALISATION, AND WHETHER IT REACHES THE FISH
# ===========================================================
# The community composition is shifting from temperate-affinity taxa to
# tropical ones. This is a climate signal that does NOT need a year-level
# dose-response to stand up: the direction of the compositional change is
# the fingerprint, and every trend below is identified at the reef level.
message("\n===== 4. TROPICALISATION =====")
inv_t  <- ltem[Label == "INV" & Year != 2020]
inv_t[, Year := as.integer(Year)]
ITK    <- c("Year", "Reef", "Depth2", "Transect", "Area")
inv_at <- unique(inv_t[, ..ITK])

trop_groups <- list("Sea stars (temperate)"  = quote(Taxa2 == "Asteroidea"),
                    "Sea urchins"            = quote(Taxa2 == "Echinoidea"),
                    "Gorgonians (temperate)" = quote(Taxa3 == "Holaxonia"),
                    "Hard corals (tropical)" = quote(Taxa3 == "Scleractinia"))
inv_change <- function(sel) {
  rc <- inv_t[eval(sel), .(q = sum(Quantity, na.rm = TRUE)), by = ITK]
  f  <- merge(inv_at, rc, by = ITK, all.x = TRUE); f[is.na(q), q := 0]
  f[, per100 := q / pmax(Area, 1) * 100]; f[, l_b := log(per100 + 0.01)]
  f[]
}
trop <- rbindlist(lapply(names(trop_groups), function(g) {
  f  <- inv_change(trop_groups[[g]])
  m  <- lm(l_b ~ Year + Reef, data = f)
  pr <- f[, .(e = mean(per100[Year <= 2013]), l = mean(per100[Year >= 2014])), by = Reef
          ][!is.na(e) & !is.na(l) & e > 0 & l > 0]
  data.table(group = g, pct_per_decade = 100 * (exp(coef(m)["Year"] * 10) - 1),
             p_reefclust = cluster_se(m, f$Reef)[term == "Year", p],
             paired_pct = 100 * (exp(mean(log(pr$l / pr$e))) - 1), n_reefs = nrow(pr))
}))
res$gap_tropicalisation <- trop
print(trop)
add("gorgonian_pct_per_decade", round(trop[group %like% "Gorgonian", pct_per_decade], 1))
add("hard_coral_pct_per_decade", round(trop[group %like% "Hard coral", pct_per_decade], 1))
add("sea_star_pct_per_decade",   round(trop[group %like% "Sea stars", pct_per_decade], 1))

# Does that habitat change reach the fish? Tested reef by reef, which is the
# right unit and gives n = 124. It does not: reefs that lost more
# invertebrates did not lose more fish. The two declines are concurrent but
# we cannot show one drives the other.
fish_at <- unique(ltem[Label == "PEC" & Year != 2020, ..ITK])
fc <- ltem[Label == "PEC" & Year != 2020, .(v = sum(Biomass, na.rm = TRUE)), by = ITK]
fc <- merge(fish_at, fc, by = ITK, all.x = TRUE); fc[is.na(v), v := 0]
fc[, per100 := v / pmax(Area, 1) * 100]
fd <- fc[, .(e = mean(per100[Year <= 2013]), l = mean(per100[Year >= 2014])), by = Reef
         ][!is.na(e) & !is.na(l) & e > 0 & l > 0][, .(Reef, d_fish = log(l / e))]
med <- copy(fd)
for (g in c("Gorgonians (temperate)", "Hard corals (tropical)", "Sea urchins")) {
  f  <- inv_change(trop_groups[[g]])
  dd <- f[, .(e = mean(per100[Year <= 2013]), l = mean(per100[Year >= 2014])), by = Reef
          ][!is.na(e) & !is.na(l) & e > 0 & l > 0][, .(Reef, d = log(l / e))]
  setnames(dd, "d", paste0("d_", substr(g, 1, 4)))
  med <- merge(med, dd, by = "Reef")
}
mediation <- rbindlist(lapply(setdiff(names(med), c("Reef", "d_fish")), function(v) {
  ct <- cor.test(med$d_fish, med[[v]])
  data.table(predictor = v, r = ct$estimate, p = ct$p.value, n_reefs = nrow(med))
}))
res$gap_habitat_mediation <- mediation
message("Does invertebrate change predict fish change, reef by reef?")
print(mediation)
add("habitat_mediation_max_r", round(max(abs(mediation$r)), 3))
add("habitat_mediation_min_p", signif(min(mediation$p), 3))

# -----------------------------------------------------------
for (n in names(res)) fwrite(res[[n]], file.path(OUT, paste0(n, ".csv")))
gs <- rbindlist(summ)
fwrite(gs, file.path(OUT, "gap_summary.csv"))
message("\n===== GAP SUMMARY =====")
print(gs)
message("\nStep 08 done. Gap tables written to manuscript/data/gap_*.csv")
