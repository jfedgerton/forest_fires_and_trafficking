################################################################################
# 09_specification_grid.R
#
# The robustness grid: every combination of analytic choice, estimated by 2SLS,
# with nothing filtered on significance. Feeds the specification curves in 11.
#
# GRID
#   DV source        TF (workers rescued), UE (seguro-desemprego victims)
#   transformation   binary, count, log(1+y), intensity (log(1+y) | y > 0)
#   instrument set   3 singletons, 3 pairs, 1 triple
#   lag structure    L0, L1, L2, LALL
#   control x window 9 estimable combinations (see below)
#   = 2 x 4 x 7 x 4 x 9 = 2,016 specifications
#
# WINDOW IS ITS OWN DIMENSION. Control sets drawn from SIDRA 5938 stop at 2021.
# If the window were merely a side-effect of the control set, "controls change
# the answer" would be inseparable from "the last three years change the answer".
# The no-GDP sets therefore run on both windows; GDP-containing sets run
# truncated only, since they are not estimable after 2021.
#
# EVERYTHING IS 2SLS, including the count outcome. Linear IV takes a count on
# the left-hand side without difficulty - the coefficient is a level effect and
# the clustered analytic standard error is valid. FE-Poisson would force a
# control function whose analytic SE ignores the first stage and would need a
# bootstrap; that is 13, a robustness check, not this.
#
# TWO CAVEATS THE CURVE SCRIPTS ACT ON, recorded here so they travel with the
# results rather than living in someone's head:
#   1. L1 and L2 enter a lagged fire term WITHOUT the contemporaneous one. Within
#      municipality and net of these fixed effects, corr(fire_t, fire_t-2) is
#      about -0.22, so an L2-only term proxies for LOW current fire and flips
#      sign mechanically. They are estimated for completeness and excluded from
#      the headline curves; `lag_only` marks them.
#   2. In earlier versions the binary outcome and the two hurdle participation
#      equations were byte-identical regressions counted three times, which
#      inflated every pooled significance rate. There is now one binary row.
#
# Checkpointed one .rds per cell. Re-run to resume.
# Output: output/results/spec_grid_2sls.csv / .rds
################################################################################

source("00_config.R")

CELLS <- file.path(DIR_RES, "grid_cells"); dir.create(CELLS, showWarnings = FALSE)
PILOT <- as.logical(Sys.getenv("PILOT", "FALSE"))

d <- load_panel()
say("panel: ", nrow(d), " rows")

CONTROLS <- build_control_sets()
WINDOWS  <- c(Wfull = YEAR_MAX, Wtrunc = TRUNC_YEAR)

grid <- CJ(dv = names(DVS), trans = TRANS, zset = names(INSTSETS),
           lagid = names(LAGSTRUCT), kset = names(CONTROLS),
           wid = names(WINDOWS), sorted = FALSE)
grid <- grid[!(sapply(kset, function(k) needs_trunc(CONTROLS[[k]])) & wid == "Wfull")]
grid[, cell := sprintf("%s__%s__%s__%s__%s__%s", dv, trans, zset, lagid, kset, wid)]
if (PILOT) grid <- grid[sample(.N, min(20L, .N))]
say("grid cells: ", nrow(grid), if (PILOT) "  [PILOT]" else "")

t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  g  <- grid[i]
  fp <- file.path(CELLS, paste0(g$cell, ".rds"))
  if (file.exists(fp) && !PILOT) next

  yv    <- DVS[[g$dv]]
  ctrl  <- CONTROLS[[g$kset]]
  Ls    <- LAGSTRUCT[[g$lagid]]
  fires <- paste0(FIRE, sapply(Ls, sfx))
  zs    <- as.vector(outer(INSTSETS[[g$zset]], sapply(Ls, sfx), paste0))

  dd <- d[year <= WINDOWS[[g$wid]]]
  need <- c(yv, fires, zs, ctrl)
  dd <- dd[stats::complete.cases(dd[, ..need])]
  dd <- apply_transform(copy(dd), yv, g$trans)
  if (nrow(dd) < 500L || length(zs) < length(fires)) {
    saveRDS(data.table(cell = g$cell, status = "not_estimable", n_obs = nrow(dd)), fp); next
  }

  m_ols <- safely(feols(as.formula(sprintf(".y ~ %s + %s | %s", rhs(fires), rhs(ctrl), FE_MAIN)),
                        data = dd, cluster = CLUSTER, notes = FALSE))
  m_fs  <- safely(feols(as.formula(sprintf("%s ~ %s + %s | %s", fires[1], rhs(zs), rhs(ctrl), FE_MAIN)),
                        data = dd, cluster = CLUSTER, notes = FALSE))
  m_rf  <- safely(feols(as.formula(sprintf(".y ~ %s + %s | %s", rhs(zs), rhs(ctrl), FE_MAIN)),
                        data = dd, cluster = CLUSTER, notes = FALSE))
  m_iv  <- safely(feols(as.formula(sprintf(".y ~ %s | %s | %s ~ %s",
                                           rhs(ctrl), FE_MAIN, rhs(fires), rhs(zs))),
                        data = dd, cluster = CLUSTER, notes = FALSE))

  ols <- pick(m_ols, fires[1]); fs <- pick(m_fs, zs[1])
  rf  <- pick(m_rf, zs[1]);     iv <- pick(m_iv, paste0("fit_", fires[1]))

  saveRDS(data.table(
    cell = g$cell, status = "ok",
    source = g$dv, transformation = g$trans, zset = g$zset, lagid = g$lagid,
    kset = g$kset, wid = g$wid,
    window = if (g$wid == "Wtrunc") "2002-2021" else "2002-2024",
    controls = paste(ctrl, collapse = "+"),
    lag_only = as.integer(g$lagid %in% c("L1", "L2")),
    n_endog = length(fires), n_inst = length(zs), n_obs = one(nobs(m_iv)),
    ols_b = one(ols[1]), ols_se = one(ols[2]), ols_p = one(ols[4]),
    fs_b  = one(fs[1]),  fs_se  = one(fs[2]),  fs_p  = one(fs[4]),
    fs_F  = first_stage_F(m_fs, zs),
    rf_b  = one(rf[1]),  rf_se  = one(rf[2]),  rf_p  = one(rf[4]),
    iv_b  = one(iv[1]),  iv_se  = one(iv[2]),  iv_p  = one(iv[4])), fp)

  if (i %% 100 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say(sprintf("cell %d/%d | %.1f min | eta %.0f min", i, nrow(grid), el,
                el / i * (nrow(grid) - i)))
  }
}

res <- rbindlist(lapply(list.files(CELLS, "\\.rds$", full.names = TRUE), readRDS), fill = TRUE)
fwrite(res, file.path(DIR_RES, "spec_grid_2sls.csv"))
saveRDS(res, file.path(DIR_RES, "spec_grid_2sls.rds"))

ok <- res[status == "ok"]
say("estimated ", nrow(ok), " of ", nrow(grid), " cells")
stopifnot(nrow(ok) > 0.9 * nrow(grid))

cat("\n===== IV ESTIMATE BY LAG STRUCTURE (lag-only setups flagged) =====\n")
print(ok[, .(n = .N, median = round(median(iv_b, na.rm = TRUE), 5),
             pct_pos = round(100 * mean(iv_b > 0, na.rm = TRUE), 1),
             pct_p05 = round(100 * mean(iv_p < 0.05, na.rm = TRUE), 1),
             median_F = round(median(fs_F, na.rm = TRUE), 1)),
         by = .(lagid, lag_only)][order(lagid)])

cat("\n===== CONTEMPORANEOUS ONLY, BY DV AND TRANSFORMATION =====\n")
print(ok[lagid == "L0", .(n = .N, median = round(median(iv_b, na.rm = TRUE), 5),
                          pct_pos = round(100 * mean(iv_b > 0, na.rm = TRUE), 1),
                          pct_p05 = round(100 * mean(iv_p < 0.05, na.rm = TRUE), 1)),
         by = .(source, transformation)][order(source, transformation)])
say("done")
