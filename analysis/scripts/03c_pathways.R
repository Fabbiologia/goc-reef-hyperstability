# 03c_pathways.R
# -----------------------------------------------------------
# Energy pathways on the balanced 26-reef panel: who carries
# the production that props up the catch, and whether the
# pathway that buffers other reef systems is doing so here.
#
# Three quantities feed main-text Figure 3:
#   (a) the share of individuals, biomass and production each
#       consumer pathway carries (the planktivore mismatch:
#       many mouths, little production);
#   (b) the trend of each pathway's production over the record
#       (reef FE + reef-clustered errors, as everywhere else);
#   (c) whether reefs drawing a larger share of production
#       from the water column respond less to warm years.
#       With year fixed effects the anomaly main effect is
#       absorbed, so the interaction is identified from
#       cross-reef contrast within years and does not rest on
#       the year-level regressor the attribution limit warns
#       about.
#
# Pathways follow the consumer classification of Parravicini
# et al. (2020) collapsed to five, mapped from the rfishprod
# diet codes in the trait table:
#   Plktiv          -> pelagic (planktivory)
#   HerMac, HerDet  -> herbivory and detritivory
#   InvMob, InvSes  -> benthic invertivory
#   FisCep          -> piscivory
#   Omnivr          -> omnivory
# The pelagic pathway indexes the allochthonous subsidy; the
# pathways beginning on the reef index autochthonous support.
#
# Requires: ../data/ltem.parquet, ../data/ltem_fish_traits.csv,
#           ../data/ltem_name_lookup.csv,
#           ../data/env/sst_gulf_by_lat_degree_daily.csv,
#           data/warm_season_anomaly_annual.csv (from 02)
# Writes:   data/pathway_shares.csv
#           data/pathway_trends.csv
#           data/pathway_fig3_series.csv
#           data/subsidy_interaction.csv
#           data/pathway_summary.csv
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(rfishprod)
})
setwd("..")
DATA <- "data"; OUT <- "data"
set.seed(20260805)
data(db, package = "rfishprod")

ltem <- as.data.table(read_parquet("../data/ltem.parquet"))
ltem[, Year := as.integer(Year)]
ltem <- ltem[Year != 2020]

trans_keys <- c("Year","Region","Reef","Habitat","Depth2","Transect",
                "Latitude","Degree","Area")

all_trans <- unique(ltem[Label == "PEC", ..trans_keys])
core      <- all_trans[, .(yrs = uniqueN(Year)), by = Reef][yrs >= 15, Reef]
fish_core_t <- all_trans[Reef %in% core]

# --- same record preparation as 03b ------------------------------------
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

# Per-reef long-term mean from OISST v2.1 at quarter degree resolution
# (00b): each reef takes its own nearest ocean cell rather than a one-degree
# band average. Long-term, not the survey year, so production still cannot
# track warming by construction.
reef_lt <- fread(file.path("../data/env", "sst_reef_longterm.csv"))
fish <- merge(fish, reef_lt[, .(Reef, sst_mean_lt)], by = "Reef", all.x = TRUE)
fish <- fish[!is.na(sst_mean_lt)]
message(sprintf("Per-reef long-term SST: %.2f to %.2f C across %d reefs",
                min(fish$sst_mean_lt), max(fish$sst_mean_lt), uniqueN(fish$Reef)))

grid <- unique(fish[, .(species_std, MaxSizeTL, Diet, Position, Method,
                        sstmean = round(sst_mean_lt, 2))])
g <- tidytrait(as.data.frame(grid), db)
kmax <- as.data.table(predKmax(g, dataset = db,
                               fmod = formula(~ sstmean + MaxSizeTL + Diet + Position + Method),
                               niter = 100, return = "pred")$pred)
fish[, sstmean := round(sst_mean_lt, 2)]
fish <- merge(fish, kmax[, .(species_std, sstmean, Kmax)],
              by = c("species_std", "sstmean"), all.x = TRUE)
stopifnot(!any(is.na(fish$Kmax)))
fish[, w := A_ord * Size^B_pen]
fish[, growth_day := somaGain(a = A_ord, b = B_pen, Lmeas = Size, t = 1,
                              Lmax = MaxSizeTL, Kmax = Kmax, silent = TRUE)]

# --- pathway assignment -------------------------------------------------
path_map <- c(Plktiv = "Pelagic (planktivory)",
              HerMac = "Herbivory and detritivory",
              HerDet = "Herbivory and detritivory",
              InvMob = "Benthic invertivory",
              InvSes = "Benthic invertivory",
              FisCep = "Piscivory",
              Omnivr = "Omnivory")
fish[, pathway := path_map[as.character(Diet)]]
stopifnot(!any(is.na(fish$pathway)))

# --- (a) shares of individuals, biomass and production ------------------
sh <- fish[, .(individuals = sum(Quantity),
               biomass     = sum(w * Quantity),
               production  = sum(growth_day * Quantity)), by = pathway]
sh[, `:=`(pct_individuals = 100 * individuals / sum(individuals),
          pct_biomass     = 100 * biomass     / sum(biomass),
          pct_production  = 100 * production  / sum(production))]
setorder(sh, -pct_production)
fwrite(sh, file.path(OUT, "pathway_shares.csv"))
message("Pathway shares on the panel:")
print(sh[, .(pathway, pct_individuals = round(pct_individuals, 1),
             pct_biomass = round(pct_biomass, 1),
             pct_production = round(pct_production, 1))])
onreef_pct  <- sh[pathway %in% c("Herbivory and detritivory", "Benthic invertivory"),
                  sum(pct_production)]
pelagic_row <- sh[pathway == "Pelagic (planktivory)"]

# --- (b) pathway production trends --------------------------------------
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

path_panel <- function(pw) {
  rec <- fish[pathway %in% pw,
              .(p = sum(growth_day * Quantity, na.rm = TRUE)), by = trans_keys]
  full <- merge(fish_core_t, rec, by = trans_keys, all.x = TRUE)
  full[is.na(p), p := 0]
  full[, p100 := p / pmax(Area, 1) * 100][]
}
resid_traj <- function(d, value = "p100", offset = 0.01) {
  d <- copy(d)[!is.na(get(value))]
  d[, l_v := log(get(value) + offset)]
  d[, r := l_v - mean(l_v, na.rm = TRUE), by = Reef]
  out <- d[, .(mean = mean(r, na.rm = TRUE),
               sem = sd(r, na.rm = TRUE) / sqrt(sum(!is.na(r))),
               count = sum(!is.na(r))), by = Year][order(Year)]
  out[, ci95 := 1.96 * sem][]
}

pathway_sets <- list(
  "Pelagic (planktivory)"     = "Pelagic (planktivory)",
  "Herbivory and detritivory" = "Herbivory and detritivory",
  "Benthic invertivory"       = "Benthic invertivory",
  "Piscivory"                 = "Piscivory",
  "Omnivory"                  = "Omnivory",
  "Autochthonous (reef-based)" = c("Herbivory and detritivory", "Benthic invertivory"))

trends <- rbindlist(lapply(names(pathway_sets), function(nm) {
  d <- path_panel(pathway_sets[[nm]])
  d[, l_v := log(p100 + 0.01)]
  fit <- lm(l_v ~ Year + Reef, data = d)
  data.table(pathway = nm,
             pct_per_decade = 100 * (exp(coef(fit)["Year"] * 10) - 1),
             p_reefclust    = clust_p(fit, d$Reef, "Year"),
             n_reefs        = uniqueN(d$Reef))
}))
fwrite(trends, file.path(OUT, "pathway_trends.csv"))
message("\nPathway production trends (reef FE, reef-clustered p):")
print(trends[, .(pathway, pct_per_decade = round(pct_per_decade, 1),
                 p = signif(p_reefclust, 2))])

# annual residual series for Figure 3b: the subsidy against the reef-based base
fig3b <- rbindlist(list(
  data.table(group = "Pelagic (planktivory)",
             resid_traj(path_panel("Pelagic (planktivory)"))),
  data.table(group = "Autochthonous (reef-based)",
             resid_traj(path_panel(c("Herbivory and detritivory", "Benthic invertivory"))))))
fwrite(fig3b, file.path(OUT, "pathway_fig3_series.csv"))

# --- (c) does the subsidy damp the thermal response? --------------------
# Reef-level pelagic share of production across the record, then the
# interaction share x warm-season anomaly with YEAR fixed effects, so the
# question is asked across reefs within the same year.
warm <- fread(file.path(OUT, "warm_season_anomaly_annual.csv"))

reef_share <- fish[, .(pel = sum(growth_day * Quantity * (pathway == "Pelagic (planktivory)")),
                       tot = sum(growth_day * Quantity)), by = Reef]
reef_share[, share := 100 * pel / tot]
message(sprintf("\nReef pelagic share: median %.1f%%, p10 %.1f%%, p90 %.1f%%",
                median(reef_share$share), quantile(reef_share$share, 0.1),
                quantile(reef_share$share, 0.9)))

tot <- path_panel(unique(unlist(pathway_sets[1:5])))   # total production
tot <- merge(tot, warm, by.x = "Year", by.y = "year")
tot <- merge(tot, reef_share[, .(Reef, share)], by = "Reef")
tot[, l_v := log(p100 + 0.01)]
tot[, share_c := share - mean(share)]
tot[, YearF := factor(Year)]

fit_int <- lm(l_v ~ Reef + YearF + ws_anom:share_c, data = tot)
b_int <- coef(fit_int)["ws_anom:share_c"]
p_int <- clust_p(fit_int, tot$Reef, "ws_anom:share_c")

# descriptive companion without year FE (anomaly main effect estimable but
# subject to the attribution limit; reported only in the supplement)
fit_desc <- lm(l_v ~ Reef + ws_anom * share_c, data = tot)
b_anom  <- coef(fit_desc)["ws_anom"]
b_intd  <- coef(fit_desc)["ws_anom:share_c"]
q       <- quantile(reef_share$share, c(0.1, 0.9)) - mean(reef_share$share)
slope_lo <- 100 * (exp(b_anom + b_intd * q[1]) - 1)   # response at p10 share
slope_hi <- 100 * (exp(b_anom + b_intd * q[2]) - 1)   # response at p90 share

inter <- data.table(
  quantity = c("interaction_beta_yearFE", "interaction_p_yearFE_reefclust",
               "descriptive_resp_pct_per_C_p10share", "descriptive_resp_pct_per_C_p90share",
               "reef_share_median", "reef_share_p10", "reef_share_p90"),
  value = c(b_int, p_int, slope_lo, slope_hi,
            median(reef_share$share), quantile(reef_share$share, 0.1),
            quantile(reef_share$share, 0.9)))
fwrite(inter, file.path(OUT, "subsidy_interaction.csv"))
fwrite(reef_share[, .(Reef, share)], file.path(OUT, "reef_pelagic_share.csv"))
message(sprintf("Interaction (year FE, reef-clustered): beta = %.4f, p = %.3g", b_int, p_int))
message(sprintf("Descriptive response at p10/p90 share: %+.1f%% / %+.1f%% per degree",
                slope_lo, slope_hi))

# per-reef descriptive thermal slopes for Figure 3c
reef_slopes <- tot[, {
  if (uniqueN(Year) >= 10) {
    ann <- .SD[, .(r = mean(l_v)), by = .(Year, ws_anom)]
    f <- lm(r ~ ws_anom, data = ann)
    .(slope = 100 * (exp(coef(f)["ws_anom"]) - 1))
  } else .(slope = NA_real_)
}, by = Reef]
reef_slopes <- merge(reef_slopes, reef_share[, .(Reef, share)], by = "Reef")
fwrite(reef_slopes, file.path(OUT, "reef_thermal_slopes.csv"))

# --- summary numbers ----------------------------------------------------
summ <- data.table(quantity = c(
  "pelagic_pct_individuals", "pelagic_pct_production", "pelagic_pct_biomass",
  "onreef_pct_production",
  "pelagic_trend_pct_decade", "pelagic_trend_p",
  "onreef_trend_pct_decade", "onreef_trend_p",
  "interaction_p"),
  value = c(round(pelagic_row$pct_individuals, 1), round(pelagic_row$pct_production, 1),
            round(pelagic_row$pct_biomass, 1), round(onreef_pct, 1),
            round(trends[pathway == "Pelagic (planktivory)", pct_per_decade], 1),
            signif(trends[pathway == "Pelagic (planktivory)", p_reefclust], 2),
            round(trends[pathway == "Autochthonous (reef-based)", pct_per_decade], 1),
            signif(trends[pathway == "Autochthonous (reef-based)", p_reefclust], 2),
            signif(p_int, 2)))
fwrite(summ, file.path(OUT, "pathway_summary.csv"))
print(summ)
message("\nStep 03c done.  Pathway tables written to manuscript/data/.")
