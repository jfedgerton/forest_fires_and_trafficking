################################################################################
# 17_channel_triage.R
#
# Three questions, none of which requires new data.
#
# TEST 1 - CROP YEAR ALIGNMENT. Everything in 16_seasonal_windows.R assumes
# SIDRA keys crops to the HARVEST year, so that for outcome year y the exposure
# window is Q4(y-1) + Q1(y) and the fire season Q3(y) falls after the harvest.
# If SIDRA keys to the PLANTING year instead, every window shifts and the
# reading collapses. The data can settle it: crop stress must show up as a
# NEGATIVE coefficient in the true exposure window, so run both alignments and
# see which one puts the negative where agronomy says it belongs.
#
# TEST 2 - DOES THE VIOLATION MATTER FOR THE LATE? 2SLS identifies an effect
# over compliers - places where VPD actually moves fire. Municipalities that
# never burn have a zero first stage and contribute zero weight, so a direct
# effect measured THERE may badly overstate the direct effect where the estimate
# actually comes from. If the VPD-to-agriculture effect shrinks as fire exposure
# rises, the bias in the estimation sample is smaller than the zero-first-stage
# calibration implies. If it does not shrink, that defence is unavailable.
#
# TEST 3 - WHAT SIGN DOES THE AGRICULTURAL CHANNEL CARRY? The seasonal run hinted
# that crop-window drying LOWERS rescues while fire-window drying RAISES them.
# If that holds, contamination from the agricultural channel pushes the fire
# estimate toward zero and the headline is a conservative bound rather than an
# inflated one. That is the difference between a caveat and a retraction, so it
# is worth establishing carefully rather than asserting from two coefficients.
#
# RESULTS WHEN RUN, 29 September 2026
#   1. Harvest-year keying wins. Exposure-window coefficients are negative
#      (Q1 -0.315, Q2 -0.471) and the fit is 3.7x better than planting-year.
#   2. The contamination-to-signal ratio falls 67x from the low to the high
#      fire-exposure tercile - first stage rises 0.28 -> 8.44 while the direct
#      agricultural effect halves. CAVEAT: the alternative reading is simply
#      that frontier municipalities have less agriculture to damage, which
#      makes the estimand a frontier LATE rather than a national one.
#   3. Opposite signs. Q2 crops -0.471 and rescues -0.0049; Q3 fire +3.99 and
#      rescues +0.0095. UNRESOLVED: Q2 drying RAISES burned area (+2.20) yet
#      LOWERS rescues, which cannot both run through the fire channel. Resolve
#      before any of this goes in a paper.
#
# Output: output/results/channel_triage_{alignment,gradient,sign}.csv
# Seed 123. Run from code/.
################################################################################

source("00_config.R")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

Z <- VPD_MAIN
AG <- c("log_crop_value", "log_crop_area_harvested", "log_crop_area_planted")

d <- load_panel()[year <= TRUNC_YEAR]

# quarter-specific VPD, keyed to an arbitrary reference year via `shift`
grab <- function(qtr, shift, nm) {
  x <- d[quarter == qtr, .(code_muni_7, year = year + shift, v = get(Z))]
  setnames(x, "v", nm); x
}
ann <- d[, c(lapply(.SD, function(x) x[1]),
             .(burn_yr = sum(percent_area_burned, na.rm = TRUE),
               tf_yr = sum(tf_count, na.rm = TRUE),
               ue_yr = sum(ue_count, na.rm = TRUE))),
         by = .(code_muni_7, year), .SDcols = AG]
ann[, `:=`(tf_bin = as.numeric(tf_yr > 0), ue_bin = as.numeric(ue_yr > 0))]

################################################################################
# TEST 1 - which alignment puts crop stress where agronomy says it belongs?
################################################################################
say("")
say("=== TEST 1: crop-year alignment ===")
say("  Harvest-year keying -> exposure is Q4(y-1) + Q1(y)")
say("  Planting-year keying -> exposure is Q4(y)   + Q1(y+1)")
say("  Crop stress must be NEGATIVE in the true exposure window.")

align <- function(tag, s_q4, s_q1, s_q2, s_q3) {
  w <- Reduce(function(a, b) merge(a, b, by = c("code_muni_7", "year")),
              list(grab(4L, s_q4, "vQ4"), grab(1L, s_q1, "vQ1"),
                   grab(2L, s_q2, "vQ2"), grab(3L, s_q3, "vQ3")))
  y <- merge(ann, w, by = c("code_muni_7", "year"))
  y <- y[stats::complete.cases(y[, .(vQ4, vQ1, vQ2, vQ3)])]
  y[, muni := factor(code_muni_7)][, yr := factor(year)]
  rbindlist(lapply(AG, function(dv) {
    m <- feols(as.formula(sprintf("%s ~ vQ4 + vQ1 + vQ2 + vQ3 | muni + yr", dv)),
               y, cluster = ~muni, notes = FALSE)
    ct <- coeftable(m)
    data.table(alignment = tag, outcome = dv,
               b_Q4 = ct["vQ4", 1], b_Q1 = ct["vQ1", 1],
               b_Q2 = ct["vQ2", 1], b_Q3 = ct["vQ3", 1],
               p_Q1 = ct["vQ1", 4], p_Q2 = ct["vQ2", 4],
               n = nobs(m), r2 = fitstat(m, "wr2")$wr2)
  }))
}
A1 <- align("harvest_year", 1L, 0L, 0L, 0L)
A2 <- align("planting_year", 0L, -1L, -1L, -1L)
AL <- rbind(A1, A2)
print(AL[, .(alignment, outcome, b_Q4 = signif(b_Q4, 3), b_Q1 = signif(b_Q1, 3),
             b_Q2 = signif(b_Q2, 3), b_Q3 = signif(b_Q3, 3), n)])
fwrite(AL, file.path(DIR_RES, "channel_triage_alignment.csv"))

say("")
say("  Which alignment puts the most negative coefficient in its exposure window?")
for (tag in c("harvest_year", "planting_year")) {
  s <- AL[alignment == tag]
  say(sprintf("   %-14s mean b in exposure window (Q4,Q1) = %+.4f | mean fit R2 = %.4f",
              tag, mean(c(s$b_Q4, s$b_Q1)), mean(s$r2)))
}

################################################################################
# TEST 2 - does the direct effect shrink where fire exposure is higher?
################################################################################
say("")
say("=== TEST 2: does VPD -> agriculture weaken where fire is stronger? ===")
w <- Reduce(function(a, b) merge(a, b, by = c("code_muni_7", "year")),
            list(grab(4L, 1L, "vQ4"), grab(1L, 0L, "vQ1"),
                 grab(2L, 0L, "vQ2"), grab(3L, 0L, "vQ3")))
y <- merge(ann, w, by = c("code_muni_7", "year"))
y <- y[stats::complete.cases(y[, .(vQ4, vQ1, vQ2, vQ3)])]
y[, muni := factor(code_muni_7)][, yr := factor(year)]

# fire exposure is a PREDETERMINED municipality characteristic - its long-run
# mean burned area - so the split is not itself moved by the instrument
expo <- y[, .(fire_bar = mean(burn_yr, na.rm = TRUE)), by = code_muni_7]
expo[, tercile := cut(fire_bar, quantile(fire_bar, c(0, 1/3, 2/3, 1)),
                      include.lowest = TRUE, labels = c("Low", "Mid", "High"))]
y <- merge(y, expo[, .(code_muni_7, fire_bar, tercile)], by = "code_muni_7")

G <- rbindlist(lapply(levels(y$tercile), function(tt) {
  s <- y[tercile == tt]
  rbindlist(lapply(AG, function(dv) {
    m <- feols(as.formula(sprintf("%s ~ vQ4 + vQ1 + vQ2 + vQ3 | muni + yr", dv)),
               s, cluster = ~muni, notes = FALSE)
    mf <- feols(burn_yr ~ vQ4 + vQ1 + vQ2 + vQ3 | muni + yr, s,
                cluster = ~muni, notes = FALSE)
    ct <- coeftable(m)
    data.table(tercile = tt, outcome = dv,
               crop_stress_Q2 = ct["vQ2", 1], p_Q2 = ct["vQ2", 4],
               fire_window_Q3 = ct["vQ3", 1], p_Q3 = ct["vQ3", 4],
               first_stage_Q3 = coeftable(mf)["vQ3", 1],
               n_muni = uniqueN(s$code_muni_7))
  }))
}))
print(G[, .(tercile, outcome, crop_stress_Q2 = signif(crop_stress_Q2, 3),
            fire_window_Q3 = signif(fire_window_Q3, 3),
            first_stage_Q3 = signif(first_stage_Q3, 3), n_muni)])
fwrite(G, file.path(DIR_RES, "channel_triage_gradient.csv"))

################################################################################
# TEST 3 - what sign does each channel carry into rescues?
################################################################################
say("")
say("=== TEST 3: sign of each channel on rescues ===")
S <- rbindlist(lapply(c("tf_bin", "ue_bin", "burn_yr", "log_crop_value"),
  function(dv) {
    m  <- feols(as.formula(sprintf("%s ~ vQ4 + vQ1 + vQ2 + vQ3 | muni + yr", dv)),
                y, cluster = ~muni, notes = FALSE)
    m2 <- feols(as.formula(sprintf("%s ~ vQ4 + vQ1 + vQ2 + vQ3 | muni + yr", dv)),
                y, cluster = ~code_muni_7 + yr, notes = FALSE)
    ct <- coeftable(m); ct2 <- coeftable(m2)
    data.table(outcome = dv, window = c("vQ4", "vQ1", "vQ2", "vQ3"),
               b = ct[c("vQ4","vQ1","vQ2","vQ3"), 1],
               p = ct[c("vQ4","vQ1","vQ2","vQ3"), 4],
               p_2way = ct2[c("vQ4","vQ1","vQ2","vQ3"), 4])
  }))
print(dcast(S, outcome ~ window, value.var = "b")[, lapply(.SD, function(x)
  if (is.numeric(x)) signif(x, 3) else x)])
say("  p-values, municipality clustered:")
print(dcast(S, outcome ~ window, value.var = "p")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 4) else x)])
fwrite(S, file.path(DIR_RES, "channel_triage_sign.csv"))

say("")
say("=== READ-OUT ===")
cs <- S[outcome == "log_crop_value" & window == "vQ2", b]
rs <- S[outcome == "tf_bin" & window == "vQ2", b]
fs <- S[outcome == "burn_yr" & window == "vQ3", b]
rf <- S[outcome == "tf_bin" & window == "vQ3", b]
bq2 <- S[outcome == "burn_yr" & window == "vQ2", b]
say(sprintf("  crop window Q2:  crops %+.4f, fire %+.4f, rescues %+.5f", cs, bq2, rs))
say(sprintf("  fire window Q3:  fire  %+.4f, rescues %+.5f", fs, rf))
if (sign(rs) != sign(rf)) {
  say("  OPPOSITE SIGNS. The agricultural channel pulls rescues the other way")
  say("  from the fire channel, so contamination biases the fire estimate")
  say("  TOWARD ZERO. The headline would be a conservative bound.")
} else {
  say("  SAME SIGN. The agricultural channel reinforces the fire channel, so")
  say("  contamination inflates the fire estimate. That is the bad case and it")
  say("  must be reported as such.")
}
# Internal consistency check that the 29 September run FAILED. Flag it loudly
# rather than leaving it for a referee.
if (sign(bq2) == sign(fs) && sign(rs) != sign(rf))
  say("  WARNING, UNRESOLVED: Q2 drying raises burned area but lowers rescues. ",
      "If fire raises rescues those cannot both run through the fire channel. ",
      "Resolve before using the conservative-bound argument.")
say("")
say("done")
