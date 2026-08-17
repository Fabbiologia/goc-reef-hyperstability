# 09_decoupling.R
# -----------------------------------------------------------
# Why the landing statistics did not show what the reefs were doing.
#
# Landings-based management implicitly assumes catch per unit effort is
# proportional to stock: CPUE = q * B, so log(CPUE) = log(q) + b*log(B)
# with b = 1. A value b < 1 is hyperstability, where catch rates hold up
# while the stock falls. This step estimates b against an INDEPENDENT
# stock index, the LTEM survey biomass, which is what makes the test
# possible: the survey never sees the fishery and the fishery never sees
# the survey.
#
# Three questions:
#   1. Does catch per trip track survey biomass at all?
#   2. Can proportionality (b = 1) be rejected?
#   3. Does the relationship degrade as heat exposure rises?
#
# Reads:  ../data/ltem.parquet, data/artisanal_bcs_annual.csv,
#         data/artisanal_bcs_yearly_totals.csv, data/warm_season_anomaly_annual.csv
# Writes: data/decoupling_series.csv, data/decoupling_models.csv,
#         data/decoupling_summary.csv  (Figure 3c is drawn in 09_figures.R)
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(arrow); library(data.table) })
setwd("..")
DATA <- "data"; OUT <- "data"

BREAK_YEAR <- 2008   # receipt filing changes here; effort not comparable across it

ltem <- as.data.table(read_parquet(file.path("../data", "ltem.parquet")))
ltem[, Year := as.integer(Year)]
ltem <- ltem[Year != 2020]
art  <- fread(file.path(DATA, "artisanal_bcs_annual.csv"))
yt   <- fread(file.path(DATA, "artisanal_bcs_yearly_totals.csv"))
warm <- fread(file.path(DATA, "warm_season_anomaly_annual.csv"))

# -----------------------------------------------------------
# Independent stock index: reef-adjusted survey biomass
# -----------------------------------------------------------
tk   <- c("Year", "Reef", "Depth2", "Transect", "Area")
allt <- unique(ltem[Label == "PEC", ..tk])
core <- allt[, .(y = uniqueN(Year)), by = Reef][y >= 15, Reef]
rec  <- ltem[Label == "PEC" & Reef %in% core, .(b = sum(Biomass, na.rm = TRUE)), by = tk]
f    <- merge(allt[Reef %in% core], rec, by = tk, all.x = TRUE)
f[is.na(b), b := 0]
f[, per100 := b / pmax(Area, 1) * 100]
f[, l_b := log(per100 + 0.01)]
f[, reef_mean := mean(l_b), by = Reef]
f[, adj := l_b - reef_mean]
B <- f[, .(B = mean(adj)), by = Year][order(Year)]

# Fishery side: reef landings and trips
C <- art[is_reef == TRUE, .(landings_t = sum(landings_t, na.rm = TRUE)), by = .(Year = year)]
d <- merge(merge(B, C, by = "Year"), yt[, .(Year = year, trips = folios)], by = "Year")
d <- merge(d, warm[, .(Year = year, ws_anom)], by = "Year")[order(Year)]
d[, cpue := landings_t / trips]
d[, l_cpue := log(cpue)]
fwrite(d, file.path(OUT, "decoupling_series.csv"))

summ <- list()
add  <- function(k, v) summ[[length(summ) + 1]] <<- data.table(quantity = k, value = v)

# -----------------------------------------------------------
# 1. Does catch per trip track the survey at all?
# -----------------------------------------------------------
ct <- cor.test(d$B, d$l_cpue)
message(sprintf("corr(survey biomass, log CPUE) = %+.2f (p = %.3f, n = %d)",
                ct$estimate, ct$p.value, nrow(d)))
add("corr_biomass_cpue", round(ct$estimate, 3))
add("corr_biomass_cpue_p", signif(ct$p.value, 3))

# -----------------------------------------------------------
# 2. Test proportionality. H0: b = 1 is the assumption that lets a
#    manager read stock status off catch rates.
# -----------------------------------------------------------
test_b <- function(dat, label) {
  m  <- lm(l_cpue ~ B, data = dat)
  cf <- summary(m)$coefficients
  b  <- cf["B", 1]; se <- cf["B", 2]; df <- nrow(dat) - 2
  p1 <- 2 * pt(-abs((b - 1) / se), df)
  data.table(window = label, n = nrow(dat), beta = b,
             lo = b - qt(.975, df) * se, hi = b + qt(.975, df) * se,
             p_beta_eq_0 = cf["B", 4], p_beta_eq_1 = p1)
}
prop <- rbind(test_b(d, "full record"),
              test_b(d[Year >= BREAK_YEAR], "post-2008 (no reporting break)"))
message("\nProportionality test (b = 1 means CPUE tracks biomass):")
print(prop[, .(window, n, beta = round(beta, 3),
               CI = sprintf("%+.2f to %+.2f", lo, hi),
               p_vs_1 = signif(p_beta_eq_1, 3))])
add("hyperstability_beta_post2008", round(prop[2, beta], 3))
add("hyperstability_beta_lo", round(prop[2, lo], 3))
add("hyperstability_beta_hi", round(prop[2, hi], 3))
add("p_proportionality_rejected", signif(prop[2, p_beta_eq_1], 3))

# -----------------------------------------------------------
# 3. Does the relationship degrade as heat exposure rises?
#    Low power (11 vs 14 years), so reported as a direction.
# -----------------------------------------------------------
era <- rbind(test_b(d[Year <  2014], "before 2014"),
             test_b(d[Year >= 2014], "2014 onward"))
message("\nBy era:"); print(era[, .(window, n, beta = round(beta, 3),
                                    CI = sprintf("%+.2f to %+.2f", lo, hi))])
mi <- lm(l_cpue ~ B * I(Year >= 2014), data = d)
add("era_interaction_p", signif(summary(mi)$coefficients[4, 4], 3))
message(sprintf("interaction p = %.3f (low power: %d vs %d years)",
                summary(mi)$coefficients[4, 4], era[1, n], era[2, n]))

fwrite(rbind(prop, era), file.path(OUT, "decoupling_models.csv"))
fwrite(rbindlist(summ), file.path(OUT, "decoupling_summary.csv"))

message("Step 07 done.  Decoupling tables written to manuscript/data/.")
