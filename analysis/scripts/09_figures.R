# 07_figures.R
# -----------------------------------------------------------
# Single figure generator for the manuscript. Replaces the
# absent external (Python) generator: every figure is rebuilt
# from the pipeline CSVs in manuscript/data/.
#
# Outputs (figures/): Figure1_climate_signal, Figure2_reef_community,
# Figure3_hyperstability (4 panels incl. the buffer projection), FigureS1..S13
# (PDF + PNG; S14 is written by 08_gap_analysis.R), and a
# traceable in_text_statistics.csv with every number quoted in
# the manuscript.
#
# NOTE: Figure 1 panels (c) and (d) (diving-smartwatch temperature
# profiles) require the raw profile dataset, which is NOT in the
# repository. If data/smartwatch_profiles.csv is present they are
# drawn; otherwise Figure 1 is rendered with panels (a) and (b)
# only and a message is printed. The headline SST value is the
# raw anomaly above the 1991-2020 average (+2.91 C), used
# consistently here and in the text.
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
setwd("..")
DATA <- "data"; OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_minimal(base_size = 9) +
          theme(plot.title = element_text(face = "bold", size = 10),
                panel.grid.minor = element_blank()))
RED <- "#c0392b"; BLUE <- "#2c7fb8"; GREEN <- "#27ae60"; ORANGE <- "#e67e22"; GREY <- "#7f8c8d"
save_fig <- function(p, name, w, h) {
  for (ext in c(".pdf", ".png")) {
    f <- file.path(OUT, paste0(name, ext))
    if (file.exists(f)) unlink(f)   # overwrite of a synced file can time out; create fresh
    dev <- if (ext == ".pdf" && capabilities("cairo")) grDevices::cairo_pdf else NULL
    if (is.null(dev)) ggsave(f, p, width = w, height = h, units = "in", dpi = 300)
    else ggsave(f, p, width = w, height = h, units = "in", device = dev)
  }
  message("wrote ", name, " (.pdf/.png)")
}
mhw_windows <- data.table(
  x1 = as.Date(c("2014-01-01","2019-01-01","2023-01-01","2026-01-01")),
  x2 = as.Date(c("2016-12-31","2020-12-31","2024-12-31","2026-04-30")))
stat_rows <- list()
addstat <- function(q, v) stat_rows[[length(stat_rows)+1]] <<- data.table(quantity = q, value = v)

# ===========================================================
# FIGURE 1 -- climate signal
# ===========================================================
mon <- fread(file.path(DATA, "sst_gulf_monthly_1981_2026.csv"))
mon[, dt := as.Date(paste0(month, "-01"))]
mon[, mo := as.integer(format(dt, "%m"))]
mon[, yr := as.integer(format(dt, "%Y"))]
clim <- mon[yr %between% c(1991, 2020), .(clim = mean(mean_sst)), by = mo]
mon <- merge(mon, clim, by = "mo")[order(dt)]
mon[, anom := mean_sst - clim]                     # raw anomaly above 1991-2020 average
apr26 <- mon[yr == 2026 & mo == 4, anom]
trend_dec <- coef(lm(anom ~ I(yr + (mo-0.5)/12), data = mon))[2] * 10
addstat("april_2026_anomaly_above_1991_2020_C", round(apr26, 3))
addstat("whole_gulf_trend_C_per_decade", round(trend_dec, 3))
prev <- mon[yr != 2026][which.max(anom)]   # previous record, excluding the 2026 event
addstat("previous_record_anomaly_C", round(prev$anom, 3))
addstat("previous_record_month", prev$month)
latt <- fread(file.path(DATA, "sst_trend_by_latitude.csv"))
addstat("northern_gulf_max_trend_C_per_decade", round(max(latt$trend_dec), 3))
ssd <- fread(file.path(DATA, "ltem_seastars_annual.csv")); setnames(ssd, "mean", "dens")
addstat("seastar_post2014_vs_pre2014_pct",
        round(100 * (mean(ssd[Year >= 2014, dens], na.rm = TRUE) /
                     mean(ssd[Year < 2014, dens], na.rm = TRUE) - 1), 1))

p1a <- ggplot(mon, aes(dt, anom)) +
  geom_col(aes(fill = anom > 0), width = 31, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = BLUE)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  annotate("point", x = mon[yr==2026 & mo==4, dt], y = apr26, colour = "black", size = 1.2) +
  annotate("text", x = as.Date("2012-06-01"), y = apr26 * 0.78,
           label = sprintf("April 2026\n+%.2f°C (record)", apr26),
           colour = RED, fontface = "bold", size = 3, hjust = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = "SST anomaly (°C)",
       title = "a")

p1b <- ggplot(mon, aes(dt, cumulative_intensity)) +
  geom_col(fill = RED, width = 31, show.legend = FALSE) +
  labs(x = NULL, y = "MHW intensity (°C·days)",
       title = "b")

prof_csv <- file.path(DATA, "smartwatch_profiles.csv")
if (file.exists(prof_csv)) {
  pr <- fread(prof_csv)   # expects: phase, depth_m, temp_C  (+ optional anom_C)
  p1c <- ggplot(pr, aes(temp_C, depth_m, colour = phase)) +
    geom_path(linewidth = 0.8) + scale_y_reverse() +
    scale_colour_manual(values = c(Before = BLUE, During = RED, After = "grey50")) +
    labs(x = "Temperature (°C)", y = "Depth (m)", colour = NULL,
         title = "c")
  p1d <- ggplot(pr[phase == "During"], aes(depth_m)) +
    labs(title = "d")
  fig1 <- (p1a | p1b) / (p1c | p1d)
  save_fig(fig1, "Figure1_climate_signal", 9, 7)
} else {
  message("NOTE: data/smartwatch_profiles.csv absent -> Figure 1 c,d not regenerated; ",
          "writing Figure1ab_climate_regenerated and leaving the deposited 4-panel Figure1 intact.")
  fig1 <- p1a / p1b
  save_fig(fig1, "Figure1ab_climate_regenerated", 7, 5)
}

# ===========================================================
# FIGURE 2 -- reef community reorganisation (reef-adjusted residuals)
# ===========================================================
f2 <- fread(file.path(DATA, "ltem_fig2_residuals.csv"))
mhwyr <- data.table(x1 = c(2014, 2019, 2023, 2026), x2 = c(2016.99, 2020.99, 2024.99, 2026.4))

reef_panel <- function(d, lv, col, lty, shp, ttl, ylab, brk = 2) {
  d <- copy(d); d[, group := factor(group, levels = lv)]
  ggplot(d, aes(Year, mean, colour = group, fill = group)) +
    geom_rect(data = mhwyr, inherit.aes = FALSE,
              aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf), fill = "grey80", alpha = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_ribbon(aes(ymin = mean - ci95, ymax = mean + ci95), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = group), linewidth = 0.7) +
    geom_point(aes(shape = group), size = 1.5) +
    scale_colour_manual(values = col) + scale_fill_manual(values = col) +
    scale_linetype_manual(values = lty) + scale_shape_manual(values = shp) +
    scale_x_continuous(breaks = seq(1998, 2024, brk)) +
    labs(x = NULL, y = ylab, title = ttl, colour = NULL, fill = NULL, linetype = NULL, shape = NULL) +
    guides(colour = guide_legend(ncol = 2)) +
    theme(legend.position = c(0.015, 0.04), legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
          legend.key.size = unit(0.8, "lines"), legend.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
inv_lv  <- c("Sea stars (Asteroidea)", "Sea urchins (Echinoidea)",
             "Hard corals (Scleractinia)", "Gorgonians (Holaxonia)")
p2a <- reef_panel(f2[panel == "a"], inv_lv,
                  setNames(c("#1f77b4","#e8a33d","#2ca25f","#d6649c"), inv_lv),
                  setNames(c("solid","dashed","dotted","twodash"), inv_lv),
                  setNames(c(16,15,17,18), inv_lv),
                  "a",
                  "Reef-adjusted log-density\n(deviation from reef mean)", brk = 4)
fish_lv <- c("Whole reef-fish community", "Commercial species", "Non-commercial species")
p2b <- reef_panel(f2[panel == "b"], fish_lv,
                  setNames(c("#1f77b4","#7b241c","#2ca25f"), fish_lv),
                  setNames(c("solid","dashed","dotted"), fish_lv),
                  setNames(c(16,15,17), fish_lv),
                  "b",
                  "Reef-adjusted log-biomass\n(deviation from reef mean)", brk = 4)
# Panels c,d: production and turnover on the same panel (from 03b). Panel c
# shows the stock and the flow diverging; panel d shows that divergence as the
# rising pace of renewal, the engine that keeps catch coming off a shrinking
# stock (the precondition for hyperstability tested in Figure 3c).
fcd <- fread(file.path(DATA, "productivity_fig2cd_residuals.csv"))
sf_lv <- c("Standing biomass", "Biomass production")
p2c <- reef_panel(fcd[panel == "c"], sf_lv,
                  setNames(c("#1f77b4","#c0392b"), sf_lv),
                  setNames(c("solid","dashed"), sf_lv),
                  setNames(c(16,15), sf_lv),
                  "c",
                  "Reef-adjusted log-value\n(deviation from reef mean)", brk = 4)
to_lv <- c("Whole community", "Commercial species")
p2d <- reef_panel(fcd[panel == "d"], to_lv,
                  setNames(c("#1f77b4","#7b241c"), to_lv),
                  setNames(c("solid","dashed"), to_lv),
                  setNames(c(16,15), to_lv),
                  "d",
                  "Reef-adjusted log-turnover\n(deviation from reef mean)", brk = 4)
save_fig((p2a | p2b) / (p2c | p2d), "Figure2_reef_community", 8.2, 7.6)

# Headline production and turnover numbers, traceable
pt <- fread(file.path(DATA, "productivity_panel_trends.csv"))
gp <- function(g, r, col) pt[group == g & response == r][[col]]
addstat("panel_turnover_median_pct_day",
        fread(file.path(DATA, "productivity_panel_summary.csv"))[
          quantity == "panel_turnover_median_pct_day", value])
addstat("community_turnover_pct_per_decade", round(gp("Whole community", "turnover", "pct_per_decade"), 1))
addstat("community_turnover_p", signif(gp("Whole community", "turnover", "p_reefclust"), 2))
addstat("commercial_turnover_pct_per_decade", round(gp("Commercial species", "turnover", "pct_per_decade"), 1))
addstat("commercial_turnover_p", signif(gp("Commercial species", "turnover", "p_reefclust"), 2))
addstat("community_production_paired_pct", round(gp("Whole community", "production", "paired_pct"), 1))
addstat("community_production_paired_lo", round(gp("Whole community", "production", "paired_lo"), 1))
addstat("community_production_paired_hi", round(gp("Whole community", "production", "paired_hi"), 1))
addstat("community_production_pct_per_decade", round(gp("Whole community", "production", "pct_per_decade"), 1))
addstat("community_turnover_paired_pct", round(gp("Whole community", "turnover", "paired_pct"), 1))
addstat("commercial_production_paired_pct", round(gp("Commercial species", "production", "paired_pct"), 1))
addstat("commercial_biomass_pct_per_decade", round(gp("Commercial species", "biomass", "pct_per_decade"), 1))
addstat("subset_biomass_pct_per_decade", round(gp("Whole community", "biomass", "pct_per_decade"), 1))

# ===========================================================
# FIGURE 3 -- energy pathways: the compensation has no buffer
# ===========================================================
psh <- fread(file.path(DATA, "pathway_shares.csv"))
ptr <- fread(file.path(DATA, "pathway_trends.csv"))
pfs <- fread(file.path(DATA, "pathway_fig3_series.csv"))
rts <- fread(file.path(DATA, "reef_thermal_slopes.csv"))
pin <- fread(file.path(DATA, "subsidy_interaction.csv"))

psh[, pathway := factor(pathway, levels = rev(psh[order(-pct_production), pathway]))]
shm <- melt(psh[, .(pathway, `Share of individuals` = pct_individuals,
                    `Share of production` = pct_production)],
            id.vars = "pathway", variable.name = "what", value.name = "pct")
p3n_a <- ggplot(shm, aes(pct, pathway, fill = what)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f", pct)), hjust = -0.25, size = 2.6,
            position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c(`Share of individuals` = "grey70",
                               `Share of production` = "#1f77b4")) +
  scale_x_continuous(limits = c(0, 68), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Percent of panel total", y = NULL, fill = NULL, title = "a") +
  theme(legend.position = c(0.98, 0.92), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.size = unit(0.8, "lines"), legend.text = element_text(size = 7))

pw_lv <- c("Autochthonous (reef-based)", "Pelagic (planktivory)")
p3n_b <- reef_panel(pfs[, .(panel = "x", group, Year, mean, sem, count, ci95)],
                    pw_lv,
                    setNames(c("#2ca25f", "#1f77b4"), pw_lv),
                    setNames(c("solid", "dashed"), pw_lv),
                    setNames(c(16, 15), pw_lv),
                    "b",
                    "Reef-adjusted log-production\n(deviation from reef mean)", brk = 4)

p3n_c <- ggplot(rts[!is.na(slope)], aes(share, slope)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_smooth(method = "lm", colour = RED, fill = RED, alpha = 0.15, linewidth = 0.7) +
  geom_point(size = 1.8, colour = "black") +
  labs(x = "Reef pelagic share of production (%)",
       y = "Production response to warm years\n(% per °C, per reef)", title = "a")

# shares + pathway series move to the Supplementary (Fig. S5)
save_fig(p3n_a / p3n_b + patchwork::plot_layout(heights = c(0.8, 1)),
         "FigureS5_energy_pathways", 7.5, 7.2)

addstat("pelagic_pct_individuals", psh[pathway == "Pelagic (planktivory)", round(pct_individuals, 1)])
addstat("pelagic_pct_production",  psh[pathway == "Pelagic (planktivory)", round(pct_production, 1)])
addstat("onreef_pct_production",   psh[pathway %in% c("Herbivory and detritivory", "Benthic invertivory"),
                                       round(sum(pct_production), 1)])
addstat("pelagic_trend_pct_decade", ptr[pathway == "Pelagic (planktivory)", round(pct_per_decade, 1)])
addstat("pelagic_trend_p",          ptr[pathway == "Pelagic (planktivory)", signif(p_reefclust, 2)])
addstat("piscivory_trend_pct_decade", ptr[pathway == "Piscivory", round(pct_per_decade, 1)])
addstat("onreef_trend_pct_decade",  ptr[pathway == "Autochthonous (reef-based)", round(pct_per_decade, 1)])
addstat("onreef_trend_p",           ptr[pathway == "Autochthonous (reef-based)", signif(p_reefclust, 2)])
addstat("subsidy_interaction_p",    pin[quantity == "interaction_p_yearFE_reefclust", signif(value, 2)])
addstat("resp_pct_perC_p10_share",  pin[quantity == "descriptive_resp_pct_per_C_p10share", round(value, 1)])
addstat("resp_pct_perC_p90_share",  pin[quantity == "descriptive_resp_pct_per_C_p90share", round(value, 1)])

# ===========================================================
# FIGURE 4 -- economic exposure
# ===========================================================
ts  <- fread(file.path(DATA, "economic_timeseries_constant_price.csv"))
pr  <- fread(file.path(DATA, "reef_paired_change.csv"))
rps <- fread(file.path(DATA, "reef_paired_summary.csv"))

# (a) What the reefs have already lost. Each reef is compared with itself,
#     early years against late years, so this cannot be driven by which reefs
#     were surveyed, and it does not use the landing statistics at all.
pr <- pr[order(pct_change)][, Reef := factor(Reef, levels = Reef)]
p3a <- ggplot(pr, aes(y = Reef)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = pct_change, yend = Reef,
                   colour = pct_change < 0), linewidth = 1.1) +
  geom_point(aes(x = pct_change, colour = pct_change < 0), size = 1.6) +
  geom_vline(xintercept = rps$pct_mean, linetype = 2, colour = "grey25", linewidth = 0.4) +
  scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = BLUE), guide = "none") +
  annotate("label", x = min(pr$pct_change) * 0.96,
           y = nrow(pr) - 3.2, hjust = 0, size = 2.6, fontface = "bold",
           colour = "#7b241c", fill = "#fdecea", label.size = 0.25,
           label = sprintf("%d of %d reefs lost fish\nmean %.0f%% (95%% CI %.0f to %.0f)",
                           rps$n_down, rps$n_reefs, rps$pct_mean, rps$ci_lo, rps$ci_hi)) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Change in reef fish biomass, 1998 to 2013 against 2014 to 2025", y = NULL,
       title = "a") +
  theme(axis.text.y = element_text(size = 5.5))

# (b) The landings gave no warning. Reef catch rose while the reef itself
#     was losing fish. The 2008 break is a change in how receipts were filed,
#     not in the fishery, so effort is not comparable across it.
p3b_cpue <- ggplot(ts, aes(year)) +
  geom_rect(data = mhw_windows, inherit.aes = FALSE,
            aes(xmin = as.integer(format(x1, "%Y")), xmax = as.integer(format(x2, "%Y")),
                ymin = -Inf, ymax = Inf), fill = RED, alpha = 0.10) +
  geom_col(aes(y = reef_t), fill = GREY, alpha = 0.55, width = 0.75) +
  geom_vline(xintercept = 2007.5, linetype = 3, linewidth = 0.5, colour = "grey20") +
  annotate("text", x = 2007.3, y = max(ts$reef_t) * 0.97, hjust = 1, size = 2.4,
           colour = "grey20", label = "receipt filing\nchanged", lineheight = 0.9) +
  labs(x = NULL, y = "Reef fish landed (tonnes)",
       title = "b")

# (c) Where the money sits, and what can actually be said about it. Most of
#     the value is in species whose response to warming is not resolved.
sp <- fread(file.path(DATA, "economic_lapaz_loreto_species.csv"))
sp[, status := fcase(
  warming_class == "decline_sig",  "Significant decline on the reefs",
  warming_class == "increase_sig", "Significant increase",
  warming_class == "not_monitored","Not monitored on the reefs",
  default =                        "No significant trend")]
sp <- sp[value_Myr > 0][order(value_Myr)]
sp[, lab := factor(reef_group, levels = reef_group)]
declining_val <- sp[status == "Significant decline on the reefs", sum(value_Myr)]
p3c <- ggplot(sp, aes(value_Myr, lab, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value_Myr)), hjust = -0.15, size = 2.7) +
  scale_fill_manual(values = c(
    "Significant decline on the reefs" = "#7b241c",
    "No significant trend"             = "#95a5a6",
    "Significant increase"             = "#2ca25f",
    "Not monitored on the reefs"       = "#d5d8dc"), name = NULL) +
  scale_x_continuous(limits = c(0, max(sp$value_Myr) * 1.22), expand = expansion(mult = c(0, 0))) +
  labs(x = "Mean annual ex-vessel value 2022 to 2025 (million USD)", y = NULL) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.7, "lines"))

# (d) THE CLOSURE. What the reef holds against what the fishery takes, on one
#     axis, with heat exposure behind it. The reef index is the monitored
#     biomass (reef-adjusted, so survey composition cannot drive it); the
#     landings index is the reef catch. Both are set to 100 at their
#     1998-2004 mean. The gap between them is the compound picture: the
#     fishery took steadily more from a reef that steadily held less, and it
#     did so as heat exposure moved outside anything in the earlier record.
f2c  <- fread(file.path(DATA, "ltem_fig2_residuals.csv"))[group == "Whole reef-fish community"]
base <- f2c[Year %in% 1998:2004, mean(mean)]
f2c[, idx := 100 * exp(mean - base)]

lnd <- ts[, .(Year = year, reef_t)]
lb  <- lnd[Year %in% 2000:2004, mean(reef_t)]
lnd[, idx := 100 * reef_t / lb]

hd <- fread(file.path(DATA, "sst_gulf_monthly_1981_2026.csv"))
hd[, `:=`(Year = as.integer(substr(month, 1, 4)), mo = as.integer(substr(month, 6, 7)))]
hd <- hd[mo %in% 5:10 & Year %in% 1998:2025,
         .(mhw_days = sum(mhw_days, na.rm = TRUE)), by = Year]

ymax <- max(c(f2c$idx, lnd$idx), na.rm = TRUE) * 1.05
hd[, bar := mhw_days / max(mhw_days) * ymax * 0.30]

p3d <- ggplot() +
  geom_col(data = hd, aes(Year, bar, fill = "Marine heatwave days per warm season"),
           alpha = 0.30, width = 0.85) +
  scale_fill_manual(values = c("Marine heatwave days per warm season" = RED), name = NULL) +
  geom_hline(yintercept = 100, linetype = 3, linewidth = 0.4, colour = "grey35") +
  geom_line(data = lnd, aes(Year, idx, colour = "Reef landings"),
            linewidth = 0.95) +
  geom_line(data = f2c, aes(Year, idx, colour = "Reef fish biomass (surveys)"),
            linewidth = 0.95) +
  geom_point(data = f2c, aes(Year, idx, colour = "Reef fish biomass (surveys)"), size = 1) +
  scale_colour_manual(values = c(
    "Reef landings" = "#1f6f9c",
    "Reef fish biomass (surveys)"   = "#7b241c"), name = NULL) +
  scale_y_continuous(limits = c(0, ymax)) +
  scale_x_continuous(breaks = seq(1998, 2024, 4)) +
  labs(x = NULL, y = "% of the 1998 to 2004 baseline",
       title = "b") +
  guides(colour = guide_legend(nrow = 2, order = 1),
         fill   = guide_legend(nrow = 1, order = 2)) +
  theme(legend.position = c(0.02, 0.97), legend.justification = c(0, 1),
        legend.spacing.y = unit(0, "cm"), legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key.size = unit(0.8, "lines"), legend.text = element_text(size = 6.8))

# (c) The formal test behind panel a: catch per trip against the independent
#     survey index. The dashed line is proportionality, the assumption that
#     lets a manager read stock status off catch rates.
dc <- fread(file.path(DATA, "decoupling_series.csv"))
p3e <- ggplot(dc, aes(B, l_cpue)) +
  geom_abline(slope = 1, intercept = mean(dc$l_cpue) - mean(dc$B),
              linetype = 2, colour = "grey40", linewidth = 0.5) +
  geom_smooth(method = "lm", se = TRUE, colour = RED, fill = RED,
              alpha = 0.12, linewidth = 0.7) +
  geom_point(aes(fill = Year), shape = 21, colour = "grey20", size = 2.2, stroke = 0.3) +
  scale_fill_gradient(low = "#d6eaf8", high = "#154360", name = NULL) +
  labs(x = "Survey biomass index (reef adjusted)",
       y = "log catch per trip", title = "c") +
  theme(legend.position = "right", legend.key.width = unit(0.35, "cm"),
        legend.text = element_text(size = 6.5))

# (d) how long the buffer holds: Worm-style extrapolation of measured rates
proj <- fread(file.path(DATA, "buffer_projection.csv"))
pjs  <- fread(file.path(DATA, "buffer_projection_summary.csv"))
gv <- function(q) pjs[quantity == q, value]
sc_lv <- c("Status quo", "Fishing plus continued warming")
proj[, scenario := factor(scenario, levels = sc_lv)]
p3f <- ggplot() +
  geom_hline(yintercept = 100, linetype = 3, linewidth = 0.4, colour = "grey60") +
  geom_hline(yintercept = c(50, 25), linetype = 3, linewidth = 0.4, colour = "grey45") +
  geom_ribbon(data = proj, aes(Year, ymin = idx_lo, ymax = idx_hi, fill = scenario), alpha = 0.13) +
  geom_line(data = proj, aes(Year, idx, colour = scenario, linetype = scenario), linewidth = 0.95) +
  scale_colour_manual(values = setNames(c("#1f6f9c", RED), sc_lv), name = NULL) +
  scale_fill_manual(values = setNames(c("#1f6f9c", RED), sc_lv), name = NULL) +
  scale_linetype_manual(values = setNames(c("solid", "longdash"), sc_lv), name = NULL) +
  annotate("text", x = 2025.6, y = 103.5, label = "the 1998 to 2004 baseline", hjust = 0,
           size = 2.4, colour = "grey45") +
  annotate("text", x = 2025.6, y = 53.5, label = "half the baseline", hjust = 0,
           size = 2.4, colour = "grey30") +
  annotate("text", x = 2025.6, y = 28.5, label = "a quarter", hjust = 0,
           size = 2.4, colour = "grey30") +
  annotate("label", x = gv("half_baseline_year_climate"), y = 50,
           label = gv("half_baseline_year_climate"), size = 2.6, colour = RED,
           fontface = "bold", label.padding = unit(0.12, "lines"), label.size = 0) +
  annotate("label", x = gv("half_baseline_year_statusquo"), y = 50,
           label = gv("half_baseline_year_statusquo"), size = 2.6, colour = "#1f6f9c",
           fontface = "bold", label.padding = unit(0.12, "lines"), label.size = 0) +
  scale_x_continuous(breaks = seq(2025, 2060, 5), limits = c(2025, 2060.5),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 108)) +
  labs(x = NULL, y = "Production (% of 1998 to 2004 baseline)", title = "d") +
  theme(legend.position = c(0.97, 0.97), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.size = unit(0.8, "lines"), legend.text = element_text(size = 6.8))

save_fig((p3n_c | p3d) / (p3e | p3f), "Figure3_hyperstability", 8.6, 7.8)
save_fig(p3c, "FigureS9_reef_value", 7, 4.6)

for (q in c("proj_rate_statusquo_pct_decade", "proj_rate_climate_pct_decade",
            "proj_sensitivity_pct_perC", "proj_warming_C_per_decade",
            "half_baseline_year_statusquo", "half_baseline_year_climate",
            "quarter_baseline_year_statusquo", "quarter_baseline_year_climate"))
  addstat(q, gv(q))

# the per-reef paired changes move to the Supplementary
save_fig(p3a, "FigureS2_reef_paired_change", 7, 6.5)

# economic in-text numbers
es <- fread(file.path(DATA, "economic_summary.csv"))
for (i in seq_len(nrow(es))) addstat(es$quantity[i], es$value[i])

# ===========================================================
# SUPPLEMENTARY FIGURES S1-S9
# ===========================================================
coefbar <- function(dt, ycol, title, xlab) {
  dt <- copy(dt); setnames(dt, ycol, "lab")
  dt[, sig := beta_anom < 0]; dt[, signif := p_anom < 0.05]
  dt <- dt[order(beta_anom)]; dt[, lab := factor(lab, levels = lab)]
  ggplot(dt, aes(beta_anom, lab, fill = sig, colour = signif)) +
    geom_col(linewidth = 0.5) +
    scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = GREEN), guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(x = xlab, y = NULL, title = title)
}

# S1 -- detrended monthly SST anomaly by 1-degree latitude band
latf <- "../data/env/sst_gulf_by_lat_degree_daily.csv"
if (file.exists(latf)) {
  d <- fread(latf); d[, date := as.Date(date)]
  d[, `:=`(yr = as.integer(format(date,"%Y")), mo = as.integer(format(date,"%m")))]
  dm <- d[, .(sst = mean(sst_mean, na.rm = TRUE)), by = .(lat_degree, yr, mo)]
  cl <- dm[yr %between% c(1991,2020), .(clim = mean(sst)), by = .(lat_degree, mo)]
  dm <- merge(dm, cl, by = c("lat_degree","mo"))[, anom := sst - clim]
  dm[, yf := yr + (mo-0.5)/12]
  dm[, anomd := anom - predict(lm(anom ~ yf)), by = lat_degree]   # detrend per band
  # Decimal-year x gives perfectly even monthly spacing, so the raster is
  # seamless (Date x leaves hairline gaps because months differ in length).
  # Each band's removed warming trend is annotated in its axis label.
  tr <- fread(file.path(DATA, "sst_trend_by_latitude.csv"))
  lat_lab <- setNames(sprintf("%d  (+%.2f)", tr$lat_degree, tr$trend_dec),
                      as.character(tr$lat_degree))
  s1 <- ggplot(dm, aes(yf, factor(lat_degree), fill = anomd)) +
    geom_raster() +
    scale_fill_gradient2(low = BLUE, mid = "white", high = RED, midpoint = 0, name = "Anomaly (°C)") +
    scale_x_continuous(breaks = seq(1985, 2025, 5), expand = c(0, 0)) +
    scale_y_discrete(labels = lat_lab, expand = c(0, 0)) +
    labs(x = NULL, y = "Latitude (°N), removed warming trend in parentheses (°C per decade)") +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7.5))
  save_fig(s1, "FigureS1_sst_by_latitude", 8, 5)
}

# S2 -- fish functional-group warming sensitivity
fg <- fread(file.path(DATA, "fish_functional_group_warming.csv"))
save_fig(coefbar(fg, "Functional_group", "Functional-group warming sensitivity",
                 "β on warm-season anomaly (Δlog-biomass per +1 °C)") +
           theme(plot.title = element_text(size = 10, face = "bold")),
         "FigureS4_funcgroup_warming", 7, 4.5)

# S3 -- per reef-group decadal biomass trend. This shows the TREND, which is
# identified at the reef level and is what the main text cites, rather than the
# warm-season sensitivity, which is a year-level regressor and unresolved.
gt <- fread(file.path(DATA, "reef_group_trend.csv"))
gt <- gt[order(pct_per_decade)][, lab := factor(reef_group, levels = reef_group)]
gt[, sig := p_year < 0.05]
s3 <- ggplot(gt, aes(pct_per_decade, lab, fill = sig)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = GREY), name = NULL,
                    labels = c(`TRUE` = "p < 0.05", `FALSE` = "not significant")) +
  labs(x = "Biomass trend (% per decade, reef-clustered)", y = NULL) +
  theme(legend.position = "bottom")
save_fig(s3, "FigureS3_top_species", 7, 5)

# Reconstruct reef CPUE from 01's artisanal tables (matches 04).
# Industrial vessels were dropped in 01 and the non-reef blocs (small
# pelagics, large pelagics, sharks, squid) are no longer analysed: this
# manuscript is about the small-scale rocky reef fishery.
art  <- fread(file.path(DATA, "artisanal_bcs_annual.csv"))
yt   <- fread(file.path(DATA, "artisanal_bcs_yearly_totals.csv"))
warm <- fread(file.path(DATA, "warm_season_anomaly_annual.csv"))
trips <- yt[, .(year, folios)]

reef_yr <- merge(art[is_reef == TRUE, .(landings_t = sum(landings_t)), by = year],
                 trips, by = "year")
reef_yr[, CPUE := landings_t / pmax(folios, 1)]
reef_yr <- merge(reef_yr, warm, by = "year")

# S6 -- total artisanal activity (landings + folios)
s6a <- ggplot(yt, aes(year, landings_t)) + geom_col(fill = ORANGE) +
  labs(x = NULL, y = "Landings (t)", title = "a")
s6b <- ggplot(yt, aes(year, folios)) + geom_col(fill = BLUE) +
  geom_vline(xintercept = 2022, linetype = 2) +
  labs(x = NULL, y = "Trip notices (folios)", title = "b")
save_fig(s6a / s6b, "FigureS10_fishery_activity", 7, 5)

# S7 -- rocky reef aggregate CPUE with marine heatwave windows
s7 <- ggplot(reef_yr, aes(year)) +
  geom_rect(data = mhw_windows, inherit.aes = FALSE,
            aes(xmin = as.integer(format(x1,"%Y")), xmax = as.integer(format(x2,"%Y")),
                ymin = -Inf, ymax = Inf), fill = RED, alpha = 0.08) +
  geom_line(aes(y = CPUE), colour = "grey20") +
  geom_smooth(aes(y = CPUE), method = "loess", se = FALSE, colour = RED,
              linewidth = 0.6, span = 0.6) +
  labs(x = NULL, y = "CPUE (t per trip)")
save_fig(s7, "FigureS6_reef_cpue_trend", 7, 4)

# S8 -- per reef-species CPUE (the species modelled in 04)
reef_sp <- fread(file.path(DATA, "reef_species_two_mode.csv"))$species
sp_yr <- merge(art[is_reef == TRUE & reef_group %in% reef_sp,
                   .(landings_t = sum(landings_t)), by = .(year, species = reef_group)],
               trips, by = "year")
sp_yr[, CPUE := landings_t / pmax(folios, 1)]
s8 <- ggplot(sp_yr, aes(year, CPUE)) +
  geom_rect(data = mhw_windows, inherit.aes = FALSE,
            aes(xmin = as.integer(format(x1,"%Y")), xmax = as.integer(format(x2,"%Y")),
                ymin = -Inf, ymax = Inf), fill = RED, alpha = 0.08) +
  geom_line(colour = "grey20") +
  geom_smooth(method = "loess", se = FALSE, colour = RED, linewidth = 0.6, span = 0.7) +
  facet_wrap(~species, scales = "free_y") +
  labs(x = NULL, y = "CPUE (t per trip)")
save_fig(s8, "FigureS7_reef_species_cpue", 9, 6)

# S4 -- two-timescale climate response: per-species lag scan + per-bloc two-mode
lagsp <- fread(file.path(DATA, "lag_scan_per_species.csv"))
focus <- head(reef_sp, 4)
s4a <- ggplot(lagsp[species %in% focus], aes(lag, beta, colour = species)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_line() +
  geom_point(data = lagsp[species %in% focus & p < 0.05], size = 1.8) +
  scale_x_continuous(breaks = 0:8) + scale_colour_brewer(palette = "Set1") +
  labs(x = "Lag (years)", y = "β per +1 °C", colour = NULL,
       title = "a")
bm <- fread(file.path(DATA, "reef_species_two_mode.csv"))
setnames(bm, "species", "bloc")
bml <- melt(bm[, .(bloc, beta_short, beta_long)], id.vars = "bloc", value.name = "beta")
bml[, mode := ifelse(variable == "beta_short", "Short lag 0–3", "Long lag 4–8")]
s4b <- ggplot(bml, aes(beta, bloc, fill = mode)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c(`Short lag 0–3` = ORANGE, `Long lag 4–8` = BLUE), name = NULL) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(x = "β on warm-season anomaly", y = NULL, title = "b") +
  theme(legend.position = "bottom")
save_fig(s4a / s4b, "FigureS11_lag_effects", 7, 7)

# S5 -- per-bloc long-run CPUE trend (%/yr) + two-mode sensitivity
bm[, year_pct := (exp(beta_year) - 1) * 100]
s5a <- ggplot(bm[order(year_pct)], aes(year_pct, factor(bloc, levels = bloc))) +
  geom_col(fill = BLUE) + geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(x = "CPUE trend (%/yr)", y = NULL, title = "a")
s5 <- s5a / s4b
save_fig(s5, "FigureS12_reef_fishery", 7, 7)

# S9 -- Borcard variance partitioning
vp <- fread(file.path(DATA, "varpart_two_mode.csv"))
vpl <- melt(vp[, .(species, `Pure fishery` = pure_fishery, `Pure climate` = pure_climate,
                   Shared = shared, Unexplained = unexplained)],
            id.vars = "species", variable.name = "component", value.name = "v")
vpl[, pct := 100 * abs(v) / sum(abs(v)), by = species]
s9 <- ggplot(vpl, aes(pct, species, fill = component)) +
  geom_col() +
  scale_fill_manual(values = c(`Pure fishery` = ORANGE, `Pure climate` = BLUE,
                               Shared = GREEN, Unexplained = GREY), name = NULL) +
  labs(x = "Share of explained + unexplained variance (%)", y = NULL) +
  theme(legend.position = "bottom")
save_fig(s9, "FigureS13_variance_partitioning", 7, 4)

# S11 -- five-state artisanal value composition (moved out of Figure 3)
e5 <- fread(file.path(DATA, "economic_5state_artisanal.csv"))
e5[, value_M := value_usd / 1e6]
st_disp <- c("SINALOA"="Sinaloa","BAJA CALIFORNIA SUR"="B.C. Sur","SONORA"="Sonora",
             "BAJA CALIFORNIA"="B.C. Norte","NAYARIT"="Nayarit")
e5[, st := factor(st_disp[state], levels = c("Sinaloa","B.C. Sur","Sonora","B.C. Norte","Nayarit"))]
tot <- e5[, .(tot_M = sum(value_M)), by = year]
s11 <- ggplot(e5, aes(factor(year), value_M, fill = st)) +
  geom_col(width = 0.7) +
  geom_text(data = tot, inherit.aes = FALSE,
            aes(factor(year), tot_M, label = sprintf("%.0f M US$", tot_M)),
            vjust = -0.6, fontface = "bold", size = 2.9) +
  scale_fill_manual(values = c("Sinaloa"="#1f6f9c","B.C. Sur"="#7fb8e0","Sonora"="#2ca25f",
                               "B.C. Norte"="#f5d33f","Nayarit"="#e8762b"), name = NULL) +
  scale_y_continuous(limits = c(0, max(tot$tot_M) * 1.12), expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Year", y = "Small-scale market value (million US$)") +
  theme(legend.position = "bottom", legend.key.size = unit(0.8, "lines"))
save_fig(s11, "FigureS8_state_value", 7, 4.5)

# ===========================================================
fwrite(rbindlist(stat_rows), file.path(OUT, "in_text_statistics.csv"))
message("\nAll figures done. in_text_statistics.csv written.")
