################################################################################
# 25_first_stage_diagnostic.R
#
# RECONSTRUCTED 13 October 2026 from the session transcript after the working
# directory was lost. Never re-executed since. Verify before trusting.
#
# How damaging is the variation in first-stage strength across subsamples?
#
# Three separate questions get conflated under "the first stage is a problem",
# and they have different answers:
#
#   Q1. WEAK INSTRUMENT. Is any subsample's first stage so weak that the 2SLS
#       point estimate and its Wald confidence interval are unreliable? Only the
#       Q3 split is anywhere near the line (F = 10.6). The answer is an
#       Anderson-Rubin confidence set, which is exact in the just-identified
#       case regardless of instrument strength, so it settles the question
#       rather than arguing about F thresholds.
#
#   Q2. LATE HETEROGENEITY. The first-stage coefficient ranges over a factor of
#       32 across splits (0.24 in low-fuel municipalities to 7.72 in the Legal
#       Amazon). That is not a weak-instrument problem. It means the full-sample
#       2SLS estimate is a first-stage-weighted average of very different local
#       effects, and we should be able to say whose effect it mostly is.
#
#   Q3. EXCLUSION RESTRICTION. Under a valid exclusion restriction the reduced
#       form is the product of the first stage and a structural parameter, so
#       the reduced form must SCALE with the first stage. If the first stage
#       varies by 32x and the reduced form is flat, either the structural
#       parameter varies inversely or VPD reaches the outcome through a channel
#       other than burned area. Computing the reduced form by subsample is what
#       distinguishes these, and nothing else in the pipeline does it.
#
# WHAT THIS FOUND IN SEPTEMBER, and the caveat that matters. The reduced form
# does not scale cleanly with the first stage. But the implied beta also varies
# 14x across subsamples, which means treatment effects are genuinely
# heterogeneous - and heterogeneous beta produces the SAME signature as a direct
# channel. With nine overlapping, non-independent subsamples this cannot tell
# them apart. Do not read the intercept of rf-on-fs as a direct effect without
# saying so.
#
# Output: output/results/first_stage_diagnostic.csv   (per-split estimates)
#         output/results/ar_confidence_sets.csv       (weak-IV robust intervals)
# Seed 123. Run from code/.
################################################################################

source("00_config.R")

d <- load_panel()[year <= TRUNC_YEAR]
CT   <- build_control_sets()$K4_ag10_gdp
need <- c(unname(DVS), FIRE, VPD_MAIN, CT, "uf_code", "quarter")
d <- d[stats::complete.cases(d[, ..need])]
say("estimation sample: ", nrow(d), " municipality-quarters")

## ---- the same splits the subsample robustness uses, rebuilt identically -----
clim <- d[, .(vpd_bar = mean(get(VPD_MAIN), na.rm = TRUE),
              fire_bar = mean(get(FIRE), na.rm = TRUE)), by = code_muni_7]
clim[, vpd_tercile := cut(vpd_bar, quantile(vpd_bar, c(0, 1/3, 2/3, 1)),
                          include.lowest = TRUE, labels = c("Low", "Mid", "High"))]
clim[, fuel := fifelse(fire_bar > median(fire_bar), "Above median", "Below median")]
d <- merge(d, clim[, .(code_muni_7, vpd_tercile, fuel)], by = "code_muni_7")
LEGAL_AMAZON <- c(11, 12, 13, 14, 15, 16, 17, 21)
d[, amazon := fifelse(uf_code %in% LEGAL_AMAZON, "Legal Amazon", "Rest of Brazil")]
d[, season := fifelse(quarter == 3L, "Q3 fire season", "Q1, Q2, Q4")]

## ---- first stage, reduced form and 2SLS on identical rows -------------------
# The reduced form is the piece the pipeline never computed. It is the only one
# of the three invariant to how strong the first stage is, so it is the quantity
# to compare across splits.
fit_one <- function(dd, split, level, dvn) {
  yv <- DVS[[dvn]]
  dd <- apply_transform(copy(dd), yv, "binary")

  f_fs <- as.formula(sprintf("%s ~ %s + %s | %s", FIRE, VPD_MAIN, rhs(CT), FE_MAIN))
  f_rf <- as.formula(sprintf(".y ~ %s + %s | %s",  VPD_MAIN, rhs(CT), FE_MAIN))
  f_iv <- as.formula(sprintf(".y ~ %s | %s | %s ~ %s",
                             rhs(CT), FE_MAIN, FIRE, VPD_MAIN))
  m_fs <- feols(f_fs, data = dd, cluster = CLUSTER, notes = FALSE)
  m_rf <- feols(f_rf, data = dd, cluster = CLUSTER, notes = FALSE)
  m_iv <- feols(f_iv, data = dd, cluster = CLUSTER, notes = FALSE)

  fs <- coeftable(m_fs)[VPD_MAIN, ]
  rf <- coeftable(m_rf)[VPD_MAIN, ]
  iv <- coeftable(m_iv)[paste0("fit_", FIRE), ]

  data.table(split = split, level = level, dv = dvn,
             n = nobs(m_iv), n_muni = uniqueN(dd$code_muni_7),
             fs_b = fs[1], fs_se = fs[2], fs_F = first_stage_F(m_fs, VPD_MAIN),
             rf_b = rf[1], rf_se = rf[2], rf_t = rf[3], rf_p = rf[4],
             iv_b = iv[1], iv_se = iv[2], iv_p  = iv[4])
}

rows <- list(); k <- 0L
add <- function(dd, split, level)
  for (dvn in names(DVS)) { k <<- k + 1L; rows[[k]] <<- fit_one(dd, split, level, dvn) }

add(d, "Full sample", "All")
for (v in unique(d$season))      add(d[season == v],      "Season", v)
for (v in levels(d$vpd_tercile)) add(d[vpd_tercile == v], "Baseline VPD", v)
for (v in unique(d$fuel))        add(d[fuel == v],        "Fuel", v)
for (v in unique(d$amazon))      add(d[amazon == v],      "Legal Amazon", v)
r <- rbindlist(rows)

# Under a valid exclusion restriction rf_b == fs_b * iv_b exactly (an algebraic
# identity for just-identified 2SLS on the same rows, not a test). Check it so
# we know the three fits really are on identical samples.
r[, implied := fs_b * iv_b]
stopifnot(max(abs(r$implied - r$rf_b) / pmax(abs(r$rf_b), 1e-12)) < 1e-6)
say("identity rf = fs x iv holds to ",
    signif(max(abs(r$implied - r$rf_b) / abs(r$rf_b)), 3))
fwrite(r, file.path(DIR_RES, "first_stage_diagnostic.csv"))

## ---- Q2 / Q3: dispersion of each quantity across splits ---------------------
say("")
say("DISPERSION ACROSS SPLITS (full sample excluded)")
cv <- r[split != "Full sample",
        .(n_splits = .N,
          fs_min = min(fs_b), fs_max = max(fs_b), fs_cv = sd(fs_b) / mean(fs_b),
          rf_min = min(rf_b), rf_max = max(rf_b), rf_cv = sd(rf_b) / mean(rf_b),
          iv_min = min(iv_b), iv_max = max(iv_b), iv_cv = sd(iv_b) / mean(iv_b)),
        by = dv]
cv[, `:=`(fs_range = fs_max / fs_min, rf_range = rf_max / rf_min,
          iv_range = iv_max / iv_min)]
print(cv[, .(dv, n_splits,
             fs_range = round(fs_range, 1), fs_cv = round(fs_cv, 2),
             rf_range = round(rf_range, 1), rf_cv = round(rf_cv, 2),
             iv_range = round(iv_range, 1), iv_cv = round(iv_cv, 2))])

## ---- is the reduced form flat, or just noisily estimated? -------------------
# Splits overlap each other (a municipality is in both a fuel half and a VPD
# tercile), so the Q statistic is run WITHIN each split family, where levels are
# disjoint and the estimates are therefore independent.
say("")
say("HOMOGENEITY OF THE REDUCED FORM WITHIN EACH SPLIT FAMILY")
qtest <- r[split != "Full sample", {
  w <- 1 / rf_se^2
  mu <- sum(w * rf_b) / sum(w)
  Q  <- sum(w * (rf_b - mu)^2)
  .(k = .N, common_rf = mu, Q = Q, df = .N - 1L,
    p = pchisq(Q, .N - 1L, lower.tail = FALSE))
}, by = .(dv, split)]
print(qtest[, .(dv, split, k, common_rf = signif(common_rf, 3),
                Q = round(Q, 2), df, p = round(p, 4))])

qtest_iv <- r[split != "Full sample", {
  w <- 1 / iv_se^2
  mu <- sum(w * iv_b) / sum(w)
  Q  <- sum(w * (iv_b - mu)^2)
  .(k = .N, Q = Q, df = .N - 1L, p = pchisq(Q, .N - 1L, lower.tail = FALSE))
}, by = .(dv, split)]
say("")
say("SAME TEST ON THE 2SLS COEFFICIENTS")
print(qtest_iv[, .(dv, split, k, Q = round(Q, 2), df, p = round(p, 4))])

## ---- Q1: Anderson-Rubin confidence sets, weak-instrument robust -------------
# For a just-identified model, testing H0: beta = b0 by regressing (y - b0*x) on
# z with the same controls and fixed effects is exact whatever the first-stage
# strength. Invert it over a grid to get a confidence set. Where the first stage
# is weak the set is wide or unbounded; where it is strong it nearly matches the
# Wald interval. This replaces arguing about whether F = 10.6 clears a threshold.
#
# SLOW. Each target is ~800 feols fits on the full panel, roughly 5 minutes.
# This was the step that kept getting killed when run alongside other jobs.
say("")
say("ANDERSON-RUBIN 95% CONFIDENCE SETS (weak-instrument robust)")

ar_set <- function(dd, dvn, lo, hi, ngrid = 801L) {
  yv <- DVS[[dvn]]
  dd <- apply_transform(copy(dd), yv, "binary")
  grid <- seq(lo, hi, length.out = ngrid)
  f <- as.formula(sprintf(".ar ~ %s + %s | %s", VPD_MAIN, rhs(CT), FE_MAIN))
  keep <- vapply(grid, function(b0) {
    dd[, .ar := .y - b0 * get(FIRE)]
    m <- feols(f, data = dd, cluster = CLUSTER, notes = FALSE)
    coeftable(m)[VPD_MAIN, 4] > 0.05
  }, logical(1))
  if (!any(keep)) return(c(NA_real_, NA_real_, 0))
  c(min(grid[keep]), max(grid[keep]), mean(keep))
}

ar_rows <- list(); j <- 0L
ar_targets <- list(
  list(dd = d,                             split = "Full sample",  level = "All"),
  list(dd = d[season == "Q3 fire season"], split = "Season",       level = "Q3 fire season"),
  list(dd = d[season == "Q1, Q2, Q4"],     split = "Season",       level = "Q1, Q2, Q4"),
  list(dd = d[fuel == "Below median"],     split = "Fuel",         level = "Below median"),
  list(dd = d[vpd_tercile == "High"],      split = "Baseline VPD", level = "High")
)
for (tg in ar_targets) for (dvn in names(DVS)) {
  w <- r[split == tg$split & level == tg$level & dv == dvn]
  # search a window 12 Wald standard errors wide, far enough out that an
  # unbounded set shows up as one running to the edge of the grid
  lo <- w$iv_b - 12 * w$iv_se; hi <- w$iv_b + 12 * w$iv_se
  s <- ar_set(tg$dd, dvn, lo, hi)
  j <- j + 1L
  ar_rows[[j]] <- data.table(
    split = tg$split, level = tg$level, dv = dvn, fs_F = w$fs_F,
    iv_b = w$iv_b, wald_lo = w$iv_b - 1.96 * w$iv_se,
    wald_hi = w$iv_b + 1.96 * w$iv_se,
    ar_lo = s[1], ar_hi = s[2], ar_frac = s[3],
    grid_lo = lo, grid_hi = hi)
  say(sprintf("  %-14s %-16s %s  F=%6.1f  Wald [%8.5f,%8.5f]  AR [%8.5f,%8.5f]%s",
              tg$split, tg$level, dvn, w$fs_F,
              w$iv_b - 1.96 * w$iv_se, w$iv_b + 1.96 * w$iv_se, s[1], s[2],
              if (isTRUE(s[1] <= lo + 1e-12 || s[2] >= hi - 1e-12))
                "  <- runs to grid edge, treat as unbounded" else ""))
}
fwrite(rbindlist(ar_rows), file.path(DIR_RES, "ar_confidence_sets.csv"))

## ---- whose effect is the pooled estimate? -----------------------------------
# 2SLS on the pooled sample is a first-stage-weighted average of the subsample
# LATEs, so report the weights rather than leaving the reader to assume the
# estimate is national.
say("")
say("WEIGHT ON EACH GROUP IN THE POOLED 2SLS (binary splits only)")
wt <- r[split %in% c("Fuel", "Legal Amazon") & dv == "TF"]
wt[, raw_w := n * fs_b]
wt[, share := raw_w / sum(raw_w), by = split]
print(wt[, .(split, level, n, fs_b = round(fs_b, 3),
             iv_b = signif(iv_b, 3), weight = paste0(round(100 * share), "%"))])

say("")
say("done")
