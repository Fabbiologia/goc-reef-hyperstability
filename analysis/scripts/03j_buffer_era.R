# 03j_buffer_era.R
# -----------------------------------------------------------
# Is the buffer itself weakening over time?
#
# The savings account has two ways to fail: the capital runs out, or the
# rate of compensation falls. 03d shows the capital running down; this
# step asks whether the compensation coupling Phi = 1 - dlog(P)/dlog(B)
# is itself lower in the heatwave era. The split at 2014 is the
# manuscript's pre-registered era boundary, not chosen from these data.
# Same panel, model and clustering as 03e; interaction and stratified
# estimates, like the winter test in 03h.
#
# Reads:  data/buffer_reef_year.csv (03e)
# Writes: data/buffer_phi_era.csv, data/buffer_era_summary.csv
# -----------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(fixest) })
setwd("..")
DATA <- "data"; OUT <- "data"

d <- fread(file.path(DATA, "buffer_reef_year.csv"))
d[, `:=`(Reef = factor(Reef), YearF = factor(Year), era = as.integer(Year >= 2014))]

m  <- feols(lP ~ lB_c * era | Reef + YearF, data = d)
V  <- vcov(m, cluster = ~ Reef + YearF)
cf <- coef(m); se <- sqrt(diag(V))
b1 <- cf[["lB_c"]]; bi <- cf[["lB_c:era"]]; si <- se[["lB_c:era"]]
p_int <- 2 * pnorm(-abs(bi / si))

strat <- rbindlist(lapply(0:1, function(e) {
  me  <- feols(lP ~ lB_c | Reef + YearF, data = d[era == e])
  be  <- coef(me)[["lB_c"]]
  see <- sqrt(diag(vcov(me, cluster = ~ Reef + YearF)))[["lB_c"]]
  data.table(era = ifelse(e == 0, "1998 to 2013", "2014 to 2025"),
             phi = 1 - be, lo = 1 - be - 1.96 * see, hi = 1 - be + 1.96 * see,
             n = nobs(me), reefs = uniqueN(d[era == e, Reef]))
}))
fwrite(strat, file.path(OUT, "buffer_phi_era.csv"))

message(sprintf("Phi by era (stratified): %.3f (%.2f to %.2f) before 2014, %.3f (%.2f to %.2f) after",
                strat$phi[1], strat$lo[1], strat$hi[1],
                strat$phi[2], strat$lo[2], strat$hi[2]))
message(sprintf("Interaction: shift in Phi %+.3f (SE %.3f), p = %.2f -> %s",
                -bi, si, p_int,
                ifelse(p_int < 0.05, "significant", "a direction, not a result")))

summ <- data.table(quantity = c(
  "phi_era_pre", "phi_era_pre_lo", "phi_era_pre_hi",
  "phi_era_post", "phi_era_post_lo", "phi_era_post_hi",
  "phi_era_shift", "phi_era_shift_se", "phi_era_shift_p"),
  value = c(round(strat$phi[1], 3), round(strat$lo[1], 3), round(strat$hi[1], 3),
            round(strat$phi[2], 3), round(strat$lo[2], 3), round(strat$hi[2], 3),
            round(-bi, 3), round(si, 3), round(p_int, 2)))
fwrite(summ, file.path(OUT, "buffer_era_summary.csv"))
message("Step 03j done.")
