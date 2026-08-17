# 05_variation_partitioning.R
# -----------------------------------------------------------
# Borcard (1992) variance partitioning of LTEM reef-fish
# biomass into pure fishery, pure climate, shared, and
# unexplained components.
#
# Climate side uses the two-mode (short + long lag) model.
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(arrow); library(data.table)
})
setwd("..")
DATA <- "data"; OUT <- "data"

ltem <- as.data.table(read_parquet("../data/ltem.parquet"))
warm <- fread(file.path(DATA, "warm_season_anomaly_annual.csv"))
art  <- fread(file.path(DATA, "artisanal_bcs_annual.csv"))

ltem[, Year := as.integer(Year)]
ltem <- ltem[Year != 2020 & Label == "PEC" &
             Region %in% c("Loreto","La Paz","Corredor","Santa Rosalia","San Basilio")]

trans_keys <- c("Year","Region","Reef","Habitat","Depth2","Transect","Latitude","Degree","Area")
all_trans <- unique(ltem[, ..trans_keys])
reef_yrs <- all_trans[, .(yrs = uniqueN(Year)), by = Reef]
core <- reef_yrs[yrs >= 10, Reef]
trans_id <- all_trans[Reef %in% core]

ltem_annual <- function(genera) {
  rec <- ltem[Reef %in% core & Genus %in% genera,
              .(b = sum(Biomass, na.rm = TRUE)),
              by = trans_keys]
  full <- merge(trans_id, rec, by = trans_keys, all.x = TRUE)
  full[is.na(b), b := 0]
  full[, per100 := b / pmax(Area, 1) * 100]
  full[, l_b    := log(per100 + 0.01)]
  full[, reef_mean := mean(l_b), by = Reef]
  full[, l_b_adj  := l_b - reef_mean]
  full[, .(ltem_l_b_adj = mean(l_b_adj)), by = Year]
}

# Reef groups (shared definition written by 01) -> LTEM genera. Landings
# are the artisanal reef catch for the same group, so the fishery and
# climate predictors describe the same taxa.
rg <- fread(file.path(DATA, "reef_group_genera.csv"))
guild_map <- split(rg$genus, rg$reef_group)
guild_map <- lapply(guild_map, function(g)
  paste0(substr(g, 1, 1), tolower(substr(g, 2, nchar(g)))))   # GENUS -> Genus

# Add two-mode warm-season lag features to a per-year frame
ws <- setNames(warm$ws_anom, warm$year)
add_modes <- function(d) {
  for (k in 0:8) d[, paste0("ws_lag", k) := ws[as.character(year - k)]]
  d[, ws_short := rowMeans(.SD), .SDcols = paste0("ws_lag", 0:3)]
  d[, ws_long  := rowMeans(.SD), .SDcols = paste0("ws_lag", 4:8)]
  d
}

results <- list()
for (cn in names(guild_map)) {
  lt <- ltem_annual(guild_map[[cn]])
  setnames(lt, "Year", "year")
  catch <- art[reef_group == cn & is_reef == TRUE,
               .(landings_t = sum(landings_t)), by = year]
  catch[, l_catch := log(pmax(landings_t, 0.01))]
  catch[, l_catch_lag1 := shift(l_catch, 1)]
  s <- merge(lt, catch, by = "year", all.x = TRUE)
  s <- merge(s, warm, by = "year", all.x = TRUE)
  s <- add_modes(s)
  fit <- s[complete.cases(s$ltem_l_b_adj, s$l_catch, s$l_catch_lag1,
                          s$ws_short, s$ws_long)]
  if (nrow(fit) < 10) next
  r2_full <- summary(lm(ltem_l_b_adj ~ l_catch + l_catch_lag1 + ws_short + ws_long, data = fit))$r.squared
  r2_cli  <- summary(lm(ltem_l_b_adj ~ ws_short + ws_long, data = fit))$r.squared
  r2_fis  <- summary(lm(ltem_l_b_adj ~ l_catch + l_catch_lag1, data = fit))$r.squared
  results[[cn]] <- data.table(
    species = cn, n = nrow(fit),
    R2_full = r2_full,
    pure_fishery = r2_full - r2_cli,
    pure_climate = r2_full - r2_fis,
    shared       = r2_cli + r2_fis - r2_full,
    unexplained  = 1 - r2_full
  )
}
vp <- rbindlist(results)
fwrite(vp, file.path(OUT, "varpart_two_mode.csv"))
print(vp)
message("\nStep 05 done.")
