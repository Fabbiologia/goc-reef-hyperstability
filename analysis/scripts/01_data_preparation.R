# 01_data_preparation.R
# -----------------------------------------------------------
# Load and clean the three primary datasets:
#   * NOAA OISST v2.1 daily SST for the Gulf of California
#   * LTEM visual census (parquet)
#   * CONAPESCA avisos de arribo, ARTISANAL fleet only
#
# Fishery source: the official CONAPESCA "avisos de arribo" landing
# receipts, 2000-2026, ingested from the raw government files by
# scripts/01a_ingest_conapesca_raw.py. Every filter below is applied
# here, in the open, so the chain from the government file to the
# manuscript number is fully auditable.
#
# Two decisions define the fishery used throughout the manuscript:
#
#   1. ARTISANAL ONLY. We keep TIPO AVISO == "MENORES", the small-scale
#      panga fleet, and drop "MAYORES" (industrial vessels) and
#      "COSECHA" (aquaculture harvest). The industrial fleet works
#      different grounds, gear and species and would confound the reef
#      signal. Its share is written to artisanal_fleet_shares.csv so
#      the manuscript can state exactly what was set aside.
#
#   2. ROCKY REEF SPECIES ONLY (flagged, not yet filtered). A landed
#      species counts as a rocky reef species when its genus is
#      recorded on the LTEM rocky reef transects. This ties the catch
#      directly to the reefs we monitor, instead of relying on the
#      common-name field or a hand-written family list.
#
# Produces analysis-ready intermediate CSVs in manuscript/data/.
# -----------------------------------------------------------

PROJECT_ROOT <- ".."
setwd(PROJECT_ROOT)

suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(arrow)
})

DATA <- "../data"
OUT  <- "data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------
# 1. Climate: copy SST monthly + MHW catalog from env folder
# -----------------------------------------------------------
file.copy(file.path("../data/env/sst_gulf_monthly_1981_2026.csv"),
          file.path(OUT, "sst_gulf_monthly_1981_2026.csv"),
          overwrite = TRUE)
file.copy(file.path("../data/env/sst_gulf_mhw_events.csv"),
          file.path(OUT, "sst_gulf_mhw_events.csv"),
          overwrite = TRUE)
message("SST + MHW catalogue copied to manuscript/data/")

# -----------------------------------------------------------
# 2. CONAPESCA avisos de arribo -> parquet
# -----------------------------------------------------------
gz  <- file.path(DATA, "conapesca_avisos_pacifico_2000_2026.csv.gz")
pq  <- file.path(DATA, "conapesca_avisos_pacifico_2000_2026.parquet")
if (!file.exists(gz) && !file.exists(pq))
  stop("Run scripts/01a_ingest_conapesca_raw.py first (builds ", gz, ")")

if (!file.exists(pq) || file.mtime(pq) < file.mtime(gz)) {
  message("Reading raw ingest and writing parquet ...")
  av <- fread(gz, showProgress = FALSE)
  write_parquet(av, pq, compression = "zstd")
} else {
  message("Reading ", pq)
  av <- as.data.table(read_parquet(pq))
}
message(sprintf("  %s rows, %d-%d", format(nrow(av), big.mark = ","),
                min(av$anio_corte), max(av$anio_corte)))

av[, `:=`(estado  = toupper(trimws(estado)),
          oficina = toupper(trimws(oficina)),
          fleet   = toupper(trimws(tipo_aviso)),
          year    = as.integer(anio_corte))]

# -----------------------------------------------------------
# 3. Rocky reef species, defined from the LTEM reef transects
# -----------------------------------------------------------
# Genus is the join key: it survives the synonymy and "SPP" entries in
# the government species field, and it is the level at which the LTEM
# transects and the landing receipts actually agree.
ltem <- as.data.table(read_parquet(file.path(DATA, "ltem.parquet")))

# Being recorded on a transect is not enough on its own: the transects also
# pick up jacks, mackerel, sharks, rays and mullet passing over the reef, and
# those are separate fisheries on separate grounds. A genus must also belong
# to a rocky-reef-resident family. This positive list is what defines "rocky
# reef species" throughout the manuscript.
REEF_FAMILIES <- c("Serranidae", "Lutjanidae", "Haemulidae", "Scaridae",
                   "Labridae", "Balistidae", "Sparidae", "Malacanthidae",
                   "Kyphosidae", "Pomacanthidae", "Acanthuridae",
                   "Holocentridae", "Muraenidae", "Scorpaenidae",
                   "Cirrhitidae", "Chaetodontidae", "Pomacentridae",
                   "Diodontidae", "Ostraciidae", "Tetraodontidae")

reef_genera <- sort(unique(toupper(trimws(
  ltem[Label == "PEC" & Family %in% REEF_FAMILIES &
       !is.na(Genus) & Genus != "", Genus]))))
message(sprintf("LTEM rocky reef fish: %d genera in %d reef families (of %d fish genera recorded)",
                length(reef_genera), length(REEF_FAMILIES),
                uniqueN(ltem[Label == "PEC", Genus])))

# genus from the scientific name; strips the non-breaking spaces and
# mojibake that appear in some year files (e.g. "LUTJANUS SPP")
av[, sci := toupper(trimws(gsub("[^A-Za-z ]", " ", nombre_cientifico)))]
av[, genus := tstrsplit(trimws(gsub("\\s+", " ", sci)), " ", keep = 1)[[1]]]
av[, is_reef := !is.na(genus) & genus %in% reef_genera]

message(sprintf("  scientific name present on %.1f%% of rows; %.1f%% match a reef genus",
                100 * mean(!is.na(av$nombre_cientifico)), 100 * mean(av$is_reef)))

# Reef groups are built from the scientific name, not from the government's
# common-name field. That field is too coarse for this fishery: all triggerfish
# catch (Balistes, ~19,000 t) is filed under "OTRAS" ("other"), so a
# common-name analysis simply cannot see it. Genus is recorded per landing and
# maps one-to-one onto the LTEM reef transects.
REEF_GROUPS <- list(
  "Leopard grouper"  = "MYCTEROPERCA",
  "Cabrillas"        = c("PARALABRAX", "CEPHALOPHOLIS"),
  "Groupers"         = c("EPINEPHELUS", "HYPORTHODUS", "DERMATOLEPIS"),
  "Snappers"         = c("LUTJANUS", "HOPLOPAGRUS"),
  "Triggerfishes"    = c("BALISTES", "SUFFLAMEN", "PSEUDOBALISTES"),
  "Parrotfishes"     = c("SCARUS", "NICHOLSINA"),
  "Grunts"           = c("HAEMULON", "ANISOTREMUS", "MICROLEPIDOTUS", "ORTHOPRISTIS"),
  "Ocean whitefish"  = "CAULOLATILUS",
  "Porgies"          = "CALAMUS",
  "Wrasses"          = c("BODIANUS", "HALICHOERES"),
  "Chubs"            = "KYPHOSUS")
av[, reef_group := NA_character_]
for (g in names(REEF_GROUPS)) av[genus %in% REEF_GROUPS[[g]], reef_group := g]

# The Pacific red snapper is the most valuable reef species and is reported
# separately. It cannot be split off by scientific name: the receipts record
# it as "LUTJANUS SPP" or, wrongly, as Lutjanus campechanus, the Atlantic red
# snapper. The genus is reliable, the epithet is not, so the split uses the
# agency's own commercial category (GUACHINANGO) within genus Lutjanus.
av[reef_group == "Snappers" & nombre_principal == "GUACHINANGO",
   reef_group := "Pacific red snapper"]

cov <- av[is_reef == TRUE, .(landings_t = sum(peso_desembarcado_kg, na.rm = TRUE) / 1000),
          by = reef_group][order(-landings_t)]
message("Reef groups (Pacific, all fleets, t):"); print(cov)

fwrite(unique(av[is_reef == TRUE, .(genus, nombre_cientifico, nombre_principal)])[order(genus)],
       file.path(OUT, "reef_species_classification.csv"))

# one shared definition of the reef groups, used by 03 (LTEM warming
# response) and 06 (economic value) so the two cannot drift apart
fwrite(rbindlist(lapply(names(REEF_GROUPS), function(g)
         data.table(reef_group = g, genus = REEF_GROUPS[[g]]))),
       file.path(OUT, "reef_group_genera.csv"))

# -----------------------------------------------------------
# 3b. Price quality control and a constant-price value
# -----------------------------------------------------------
# valor_pesos equals precio_pesos x peso_desembarcado_kg in essentially every
# row, so where the value is wrong the error is in the unit price. A small
# number of receipts carry impossible prices (up to 91,000 pesos/kg against a
# median of 20), and two of them alone inflate the 2015 and 2016 Pacific
# totals from about 4-5 to 10 and 14.5 billion pesos. We flag any price above
# 20 times that species' median across the whole record, which touches 0.35
# percent of rows, and rebuild those rows from the species-year median price.
av[, med_sp := median(precio_pesos[precio_pesos > 0], na.rm = TRUE), by = nombre_principal]
av[, price_bad := !is.na(precio_pesos) & precio_pesos > 20 * med_sp]
av[, med_spyr := median(precio_pesos[!price_bad & precio_pesos > 0], na.rm = TRUE),
   by = .(nombre_principal, anio_corte)]
av[, value_mxn := fifelse(price_bad & !is.na(med_spyr),
                          med_spyr * peso_desembarcado_kg, valor_pesos)]
message(sprintf("price QC: %d of %d rows rebuilt (%.2f%%)",
                sum(av$price_bad, na.rm = TRUE), nrow(av),
                100 * mean(av$price_bad, na.rm = TRUE)))

# Nominal pesos cannot be compared across 25 years of inflation. We therefore
# also value every year's catch at a fixed reference price per species (the
# 2022 to 2025 average ex-vessel price), which removes all price movement and
# leaves a series that reflects only how much is landed and of what. This is
# the series used for the time trend; nominal value is kept for the recent
# window where prices are directly comparable.
refp <- av[price_bad == FALSE & anio_corte %in% 2022:2025 & peso_desembarcado_kg > 0,
           .(pref = sum(value_mxn, na.rm = TRUE) / sum(peso_desembarcado_kg, na.rm = TRUE)),
           by = nombre_principal]
av <- merge(av, refp, by = "nombre_principal", all.x = TRUE)
av[is.na(pref), pref := med_sp]
av[, value_const := pref * peso_desembarcado_kg]

# -----------------------------------------------------------
# 4. Fleet shares: what the artisanal filter sets aside
# -----------------------------------------------------------
shares <- av[year %in% 2001:2025,
             .(landings_t = sum(peso_desembarcado_kg, na.rm = TRUE) / 1000,
               value_mxn  = sum(valor_pesos, na.rm = TRUE)),
             by = .(year, fleet)][order(year, fleet)]
fwrite(shares, file.path(OUT, "artisanal_fleet_shares.csv"))
fl <- shares[, .(landings_t = sum(landings_t)), by = fleet][order(-landings_t)]
fl[, pct := round(100 * landings_t / sum(landings_t), 1)]
message("Fleet split of Pacific landings, 2001-2025:"); print(fl)

# -----------------------------------------------------------
# 5. ARTISANAL fleet
# -----------------------------------------------------------
art <- av[fleet == "MENORES"]

# Effort. A trip is a distinct landing receipt (folio). The CPUE
# denominator throughout is reef_folios, the receipts that landed at least
# one reef species: it is the trip count that actually fished the reef.
# folios, every artisanal receipt at the same offices, is kept alongside as
# the comparison series. The two diverge because of fisheries that have
# nothing to do with the reef: the jumbo squid fishery at Santa Rosalia
# filed ~1,800 receipts a year in 2006 to 2008, fewer than 50 a year from
# 2017 to 2024, and 1,951 in 2025. Dividing reef catch by every receipt
# therefore inflated reef CPUE as squid left the denominator and crashed it
# when squid returned, neither of which is a change on the reef. boat-days
# uses the vessel and effective-day fields the raw receipts carry, and is
# written alongside as a sensitivity check (see 04_fishery_analysis.R).
art[, boat_days := pmax(n_embarcaciones, 1, na.rm = TRUE) *
                   pmax(dias_efectivos,  1, na.rm = TRUE)]

# 2013 is present and complete in the raw government files (the earlier
# analysis excluded it because the derived snapshot we then used was missing
# it). Only 2020 is genuinely anomalous, being COVID-suppressed, and it is
# dropped from the models while remaining in the descriptive value series.
BAD_YEARS <- c(2020)
art_all <- art[!is.na(year) & year %in% 2000:2025]
art <- art_all[!(year %in% BAD_YEARS)]

# --- (a) the four Gulf-coast offices of Baja California Sur, for the reef
#         fishery analysis. CONAPESCA files them under litoral PACIFICO,
#         which is the administrative label for the whole west coast of
#         Mexico including the Gulf of California; geographically every one
#         of them faces the Gulf.
bcs <- art[estado %in% c("BAJA CALIFORNIA SUR", "BAJA CALIFORNIA") &
           grepl("LA PAZ|LORETO|SANTA ROSAL|MULEG", oficina)]
bcs[, species := toupper(trimws(nombre_especie))]

annual <- bcs[, .(landings_kg = sum(peso_desembarcado_kg, na.rm = TRUE),
                  value_mxn   = sum(value_mxn,   na.rm = TRUE),
                  value_const = sum(value_const, na.rm = TRUE),
                  records     = .N,
                  folios      = uniqueN(folio_aviso),
                  boat_days   = sum(boat_days, na.rm = TRUE)),
              by = .(year, office = oficina, species, genus, reef_group, is_reef,
                     nombre_principal)]
annual[, landings_t := landings_kg / 1000]
fwrite(annual, file.path(OUT, "artisanal_bcs_annual.csv"))

yr_tot <- bcs[, .(landings_t   = sum(peso_desembarcado_kg, na.rm = TRUE) / 1000,
                  reef_t       = sum(peso_desembarcado_kg[is_reef == TRUE], na.rm = TRUE) / 1000,
                  records      = .N,
                  folios       = uniqueN(folio_aviso),
                  reef_folios  = uniqueN(folio_aviso[is_reef == TRUE]),
                  boat_days    = sum(boat_days, na.rm = TRUE),
                  n_species  = uniqueN(species),
                  n_offices  = uniqueN(oficina)),
              by = year][order(year)]
fwrite(yr_tot, file.path(OUT, "artisanal_bcs_yearly_totals.csv"))
message(sprintf("Artisanal BCS Gulf-coast offices: %s rows, %d-%d, %.0f%% of catch from reef genera",
                format(nrow(annual), big.mark = ","), min(annual$year), max(annual$year),
                100 * annual[is_reef == TRUE, sum(landings_t)] / annual[, sum(landings_t)]))

# --- (b) the five Gulf states, for the economic step
GOC_STATES <- c("BAJA CALIFORNIA", "BAJA CALIFORNIA SUR", "SONORA", "SINALOA", "NAYARIT")
goc <- art_all[estado %in% GOC_STATES,
           .(landings_t  = sum(peso_desembarcado_kg, na.rm = TRUE) / 1000,
             value_mxn   = sum(value_mxn,   na.rm = TRUE),
             value_const = sum(value_const, na.rm = TRUE),
             folios      = uniqueN(folio_aviso)),
           by = .(year, state = estado, oficina, species = nombre_especie,
                  genus, reef_group, is_reef, nombre_principal)]
goc_effort <- art_all[estado %in% GOC_STATES,
                      .(trips      = uniqueN(folio_aviso),
                        reef_trips = uniqueN(folio_aviso[is_reef == TRUE]),
                        landings_t = sum(peso_desembarcado_kg, na.rm = TRUE) / 1000),
                      by = year][order(year)]
fwrite(goc_effort, file.path(OUT, "artisanal_5state_effort.csv"))

fwrite(goc, file.path(OUT, "artisanal_5state_annual.csv"))
message(sprintf("Artisanal 5-state: %s rows", format(nrow(goc), big.mark = ",")))

# -----------------------------------------------------------
# 6. LTEM check
# -----------------------------------------------------------
message("LTEM parquet ready at ", file.path(DATA, "ltem.parquet"))
message("\nStep 01 done. Intermediates in manuscript/data/.")
