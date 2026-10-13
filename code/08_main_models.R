################################################################################
# 08_main_models.R
#
# The main-text models. Deliberately short and separate from the full grid in
# 09, so Table 1 can be reproduced without running 2,520 specifications; 09
# re-estimates the same cells and this script asserts the two agree.
#
# INSTRUMENT: vpd_max_q, the quarterly mean of daily-max VPD, used alone.
# SPECIFICATION: contemporaneous fire, municipality and year-quarter fixed
# effects, standard errors two-way clustered on municipality and year-quarter.
#
#   First stage      percent_area_burned ~ vpd_max_q | muni + yq
#                    estimated separately for each DV, because the estimation
#                    samples differ (FE-Poisson-free here, but complete-case
#                    and control-set coverage still differ between TF and UE)
#
#   Second stages    binary      1(y > 0)        2SLS
#                    count       y               2SLS
#                    log1p       log(1 + y)      2SLS
#                    each with no controls and with the full control set
#
# All four are 2SLS. A count outcome poses no problem for linear IV - the
# coefficient is a level effect rather than a semi-elasticity, and the clustered
# analytic SE is valid. FE-Poisson would force a control function and a
# bootstrap; that route is kept in 13 as a multiplicative robustness check only.
#
# Output: output/results/main_models.rds  (fixest objects, consumed by 10 and 11)
#         output/results/main_models_tidy.csv
################################################################################

source("00_config.R")

d <- load_panel()
say("panel: ", nrow(d), " rows")

CONTROLS <- build_control_sets()
FULL_CTRL <- CONTROLS$K4_ag10_gdp          # "full controls" = ag(.10) + GDP
say("full control set: ", paste(FULL_CTRL, collapse = ", "))
stopifnot(length(FULL_CTRL) > 0)

Z <- VPD_MAIN
models <- list(); tidy <- list(); k <- 0L

for (dvn in names(DVS)) {
  yv <- DVS[[dvn]]

  # Three columns, not two. The full control set draws on SIDRA 5938 and is only
  # estimable through 2021, so "no controls" is run on BOTH windows: once on the
  # full sample (the headline), and once truncated to 2002-2021 so that the
  # with-versus-without-controls comparison differs by controls alone and not by
  # sample. Reporting only the full-window no-controls column against a
  # truncated with-controls column confounds the two.
  for (ctrl_name in c("none_fullwindow", "none", "full")) {
    ctrl <- if (ctrl_name == "full") FULL_CTRL else character(0)
    truncate <- ctrl_name != "none_fullwindow"
    dd0 <- if (truncate) d[year <= TRUNC_YEAR] else d
    win <- if (truncate) "2002-2021" else "2002-2024"

    for (trans in c("binary", "count", "log1p")) {
      need <- c(yv, FIRE, Z, ctrl)
      dd <- dd0[stats::complete.cases(dd0[, ..need])]
      dd <- apply_transform(copy(dd), yv, trans)

      f_fs <- as.formula(sprintf("%s ~ %s + %s | %s", FIRE, Z, rhs(ctrl), FE_MAIN))
      m_fs <- safely(feols(f_fs, data = dd, cluster = CLUSTER, notes = FALSE))
      f_iv <- as.formula(sprintf(".y ~ %s | %s | %s ~ %s", rhs(ctrl), FE_MAIN, FIRE, Z))
      m_iv <- safely(feols(f_iv, data = dd, cluster = CLUSTER, notes = FALSE))
      f_rf <- as.formula(sprintf(".y ~ %s + %s | %s", Z, rhs(ctrl), FE_MAIN))
      m_rf <- safely(feols(f_rf, data = dd, cluster = CLUSTER, notes = FALSE))
      f_ols <- as.formula(sprintf(".y ~ %s + %s | %s", FIRE, rhs(ctrl), FE_MAIN))
      m_ols <- safely(feols(f_ols, data = dd, cluster = CLUSTER, notes = FALSE))

      key <- sprintf("%s__%s__%s", dvn, trans, ctrl_name)
      models[[key]] <- list(fs = m_fs, iv = m_iv, rf = m_rf, ols = m_ols,
                            dv = dvn, trans = trans, controls = ctrl_name,
                            window = win, n = nobs(m_iv))

      iv <- pick(m_iv, paste0("fit_", FIRE))
      fs <- pick(m_fs, Z); rf <- pick(m_rf, Z); ols <- pick(m_ols, FIRE)
      k <- k + 1L
      tidy[[k]] <- data.table(
        dv = dvn, transformation = trans, controls = ctrl_name, window = win,
        n_obs = one(nobs(m_iv)),
        fs_b = one(fs[1]), fs_se = one(fs[2]), fs_p = one(fs[4]),
        fs_F = first_stage_F(m_fs, Z),
        ols_b = one(ols[1]), ols_se = one(ols[2]), ols_p = one(ols[4]),
        rf_b  = one(rf[1]),  rf_se  = one(rf[2]),  rf_p  = one(rf[4]),
        iv_b  = one(iv[1]),  iv_se  = one(iv[2]),  iv_p  = one(iv[4]))
      say(sprintf("%-4s %-9s %-4s | n=%7d | FS F=%6.1f | IV %+0.5f (p=%.4f)",
                  dvn, trans, ctrl_name, nobs(m_iv), first_stage_F(m_fs, Z),
                  iv[1], iv[4]))
    }
  }
}

tidy <- rbindlist(tidy)

# A fixest object keeps per-observation vectors (residuals, fitted values,
# scores). With 18 specifications x 4 models each on a 500k-row panel that is
# ~1.9 GB on disk, most of it useless downstream: 10_tables.R only needs the
# coefficient table, the fit statistics and the call. Strip the heavy slots
# before saving, and check that etable still renders afterwards.
HEAVY <- c("residuals", "fitted.values", "scores", "working_residuals",
           "linear.predictors", "sumFE", "obs_selection", "y", "X",
           "irls_weights", "weights", "iv_residuals", "fixef_id",
           "model_matrix_info")
slim <- function(m) {
  if (is.null(m) || !is.list(m)) return(m)
  for (nm in HEAVY) m[[nm]] <- NULL
  # iv_first_stage is a LIST of nested fixest objects and is the single largest
  # component (~25 MB per model). etable needs it for the ivf fit statistic, so
  # slim it recursively rather than dropping it - dropping it makes etable fail
  # with "$ operator is invalid for atomic vectors".
  if (!is.null(m$iv_first_stage))
    m$iv_first_stage <- lapply(m$iv_first_stage, slim)
  m
}
models <- lapply(models, function(x) {
  x[c("fs", "iv", "rf", "ols")] <- lapply(x[c("fs", "iv", "rf", "ols")], slim)
  x
})
stopifnot(!is.null(fixest::coeftable(models[[1]]$iv)))   # still usable
saveRDS(models, file.path(DIR_RES, "main_models.rds"))
say("main_models.rds: ",
    round(file.size(file.path(DIR_RES, "main_models.rds")) / 1e6, 1), " MB")
fwrite(tidy, file.path(DIR_RES, "main_models_tidy.csv"))

## ---- validation --------------------------------------------------------------
stopifnot(nrow(tidy) == length(DVS) * 3L * 3L)
stopifnot(all(is.finite(tidy$iv_b)))
stopifnot(all(tidy$fs_F > 10, na.rm = TRUE))   # no weak first stage in the main table
say("all ", nrow(tidy), " main models estimated; min first-stage F = ",
    round(min(tidy$fs_F, na.rm = TRUE), 1))

cat("\n===== MAIN MODELS =====\n")
print(tidy[, .(dv, transformation, controls, window, n_obs,
               fs_F = round(fs_F, 1),
               ols = round(ols_b, 5), ols_p = round(ols_p, 3),
               rf = round(rf_b, 5), rf_p = round(rf_p, 3),
               iv = round(iv_b, 5), iv_se = round(iv_se, 5), iv_p = round(iv_p, 4))])

## ---- cross-check against the full grid, if it has been run -------------------
grid_file <- file.path(DIR_RES, "spec_grid_2sls.csv")
if (file.exists(grid_file)) {
  g <- fread(grid_file)
  g <- g[zset == "Z2_max" & lagid == "L0" & kset == "K4_ag10_gdp"]
  chk <- merge(tidy[controls == "full", .(dv, transformation, iv_b)],
               g[, .(dv = source, transformation, grid_b = iv_b)],
               by = c("dv", "transformation"))
  if (nrow(chk)) {
    worst <- max(abs(chk$iv_b - chk$grid_b), na.rm = TRUE)
    say("cross-check vs 09 grid: max |difference| = ", signif(worst, 3))
    stopifnot(worst < 1e-8)
  }
} else {
  say("09 grid not yet run - skipping cross-check")
}
say("done")
