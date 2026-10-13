################################################################################
# 15_plausibly_exogenous.R
#
# Conley, Hansen & Rossi (2012) sensitivity analysis, with the direct effect
# calibrated from a zero-first-stage sample as in van Kippersluis & Rietveld
# (2018). This is the standard remedy when an instrument has a plausible second
# channel, and it is what a referee will ask for by name.
#
# THE SETUP. The exclusion restriction says VPD enters the outcome equation only
# through burned area. Relax it to
#       y = b*fire + g*z + controls + FE + e
# where g is the direct effect. Conventional IV assumes g = 0. Here g is
# estimated rather than assumed.
#
# THE ZERO-FIRST-STAGE SAMPLE. In municipalities that never record a detected
# burn, percent_area_burned is identically zero, so the first stage is not weak
# there - it is exactly zero. Any reduced-form effect of VPD on rescues in that
# subsample cannot be running through fire, so it identifies g directly. This is
# the cleanest version of the van Kippersluis-Rietveld design; most applications
# have to settle for a subsample where the first stage is merely small.
#
# WHAT IS ASSUMED, and it is not innocuous. Transporting g from the never-burn
# subsample to the full population requires the direct effect to be the SAME in
# both. Never-burn municipalities are wetter, more urban and carry 1.3% of
# task-force cases, so that is a real assumption and the reason the breakdown
# analysis matters more than any single adjusted point estimate.
#
# HOW TO READ THE 2026 RESULTS, which were over-interpreted when first produced.
# Every headline estimate lost significance under the LTZ adjustment. That was
# NOT because the direct effect was shown to be large - g-hat is statistically
# indistinguishable from zero on every outcome (p .05 to .65). The adjustment
# inflates the standard error by g's own imprecision, and the never-burn sample
# is small. The honest conclusion is that this analysis is INCONCLUSIVE here,
# not that it refutes the design. Report it that way.
#
# OUTPUTS
#   beta_ltz     local-to-zero adjusted estimate and its inflated standard error
#   uci_lo/hi    union of confidence intervals over the plausible support of g
#   g_breakdown  the value of g at which the adjusted estimate reaches zero,
#                expressed as a multiple of g-hat
#
# Output: output/results/plausibly_exogenous.csv
# Seed 123. Run from code/.
################################################################################
source("00_config.R")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

d <- load_panel()[year <= TRUNC_YEAR]
CT <- build_control_sets()$K4_ag10_gdp
need <- c(unname(DVS), FIRE, VPD_MAIN, CT); d <- d[stats::complete.cases(d[, ..need])]

burn_tot <- d[, .(burn = sum(get(FIRE), na.rm = TRUE)), by = code_muni_7]
never <- burn_tot[burn == 0, code_muni_7]
zfs <- d[code_muni_7 %in% never]
say("zero-first-stage sample: ", length(never), " municipalities, ", nrow(zfs), " rows")
stopifnot("the zero-first-stage sample is not actually zero-first-stage" =
            max(abs(zfs[[FIRE]])) == 0)
say("  verified: burned area is identically zero there, so the first stage is 0")

out <- list(); k <- 0L
for (dvn in names(DVS)) for (tr in c("binary", "log1p", "count")) {
  full <- apply_transform(copy(d),   DVS[[dvn]], tr)
  zero <- apply_transform(copy(zfs), DVS[[dvn]], tr)

  m_fs <- feols(as.formula(sprintf("%s ~ %s + %s | %s", FIRE, VPD_MAIN, rhs(CT), FE_MAIN)),
                full, cluster = CLUSTER, notes = FALSE)
  m_rf <- feols(as.formula(sprintf(".y ~ %s + %s | %s", VPD_MAIN, rhs(CT), FE_MAIN)),
                full, cluster = CLUSTER, notes = FALSE)
  m_iv <- feols(as.formula(sprintf(".y ~ %s | %s | %s ~ %s", rhs(CT), FE_MAIN, FIRE, VPD_MAIN)),
                full, cluster = CLUSTER, notes = FALSE)
  m_g  <- feols(as.formula(sprintf(".y ~ %s + %s | %s", VPD_MAIN, rhs(CT), FE_MAIN)),
                zero, cluster = CLUSTER, notes = FALSE)

  fs <- coeftable(m_fs)[VPD_MAIN, 1]
  rf <- coeftable(m_rf)[VPD_MAIN, 1]
  iv <- coeftable(m_iv)[paste0("fit_", FIRE), 1:2]
  g  <- coeftable(m_g)[VPD_MAIN, 1:2]

  # LOCAL TO ZERO. b(g) = (rf - g)/fs, and the variance picks up g's own
  # uncertainty scaled by the first stage. With one instrument this is the exact
  # Conley LTZ expression, not an approximation.
  b_ltz  <- (rf - g[1]) / fs
  se_ltz <- sqrt(iv[2]^2 + (g[2] / fs)^2)

  # UNION OF CONFIDENCE INTERVALS over g in its own 95% interval. Conservative:
  # it takes g's sampling uncertainty as the support rather than assuming a
  # tighter prior.
  g_lo <- g[1] - 1.96 * g[2]; g_hi <- g[1] + 1.96 * g[2]
  gg <- seq(g_lo, g_hi, length.out = 201)
  bs <- (rf - gg) / fs
  uci_lo <- min(bs - 1.96 * iv[2]); uci_hi <- max(bs + 1.96 * iv[2])

  # BREAKDOWN. b(g) = 0 exactly when g = rf: the direct channel would have to
  # account for the ENTIRE reduced form.
  g_break <- rf
  g_sig   <- rf - 1.96 * se_ltz * fs

  k <- k + 1L
  out[[k]] <- data.table(
    dv = dvn, trans = tr, fs = fs, rf = rf,
    iv_b = iv[1], iv_se = iv[2],
    g_hat = g[1], g_se = g[2], g_p = coeftable(m_g)[VPD_MAIN, 4],
    g_share_of_rf = g[1] / rf,
    b_ltz = b_ltz, se_ltz = se_ltz, p_ltz = 2 * pnorm(-abs(b_ltz / se_ltz)),
    uci_lo = uci_lo, uci_hi = uci_hi,
    g_break = g_break, break_mult = g_break / g[1],
    g_sig = g_sig, sig_mult = g_sig / g[1])
}
r <- rbindlist(out)
fwrite(r, file.path(DIR_RES, "plausibly_exogenous.csv"))

say("")
say("DIRECT EFFECT g, estimated in the zero-first-stage sample")
print(r[, .(dv, trans, g_hat = signif(g_hat, 3), g_se = signif(g_se, 3),
            g_p = round(g_p, 3), pct_of_reduced_form = round(100 * g_share_of_rf))])
say("  NOTE: if g_p is large on every row, the adjustment below is propagating")
say("  noise rather than demonstrating bias. Say so when reporting it.")
say("")
say("ADJUSTED ESTIMATES")
print(r[, .(dv, trans, iv = signif(iv_b, 3), iv_p = round(2*pnorm(-abs(iv_b/iv_se)), 4),
            ltz = signif(b_ltz, 3), ltz_p = round(p_ltz, 4),
            uci = sprintf("[%.5f, %.5f]", uci_lo, uci_hi))])
say("")
say("BREAKDOWN: how large the direct effect must be to zero the estimate,")
say("  as a multiple of what the zero-first-stage sample actually measures")
print(r[, .(dv, trans, g_hat = signif(g_hat, 3), g_break = signif(g_break, 3),
            break_mult = round(break_mult, 2))])
say("")
say("done")
