################################################################################
# 14_tables_meeting.R
#
# Four sets of tables for the manuscript / meeting handout, plus a compiled
# LaTeX document. Reads the estimated grid from 09; estimates nothing new, so
# every number here is traceable to a cell of spec_grid_2sls.csv.
#
#   Set 1  Main IV results. Both DVs in one table, two panels. Columns are the
#          outcome transformation crossed with the number of controls.
#          Instrument: vpd_max_q.
#   Set 2  One table per instrument operationalization (mean VPD, mean daily-max
#          VPD, degrees above 3 kPa), same layout, so the three can be read
#          against each other.
#   Set 3  First-stage statistics: instrument coefficient, standard error, joint
#          F, and sample size, by instrument and control count.
#   Set 4  Lag structures: contemporaneous, lag 1 alone, lag 2 alone, and all
#          three entered jointly.
#
# EVERY TABLE IS ON THE 2002-2021 WINDOW. The 6-control set draws on SIDRA 5938,
# which ends in 2021. Holding the window fixed across all columns means the
# with-versus-without-controls comparison differs by controls alone. A
# full-window column would confound controls with the three extra years, and in
# this project that confound was large enough to flip a sign.
#
# Output: output/tables/tab*_meeting.tex  and  output/meeting/fire_trafficking_results.pdf
################################################################################

source("00_config.R")

DIR_MTG <- file.path(DIR_OUT, "meeting"); dir.create(DIR_MTG, showWarnings = FALSE)

g <- fread(file.path(DIR_RES, "spec_grid_2sls.csv"))[status == "ok"]
g <- g[window == "2002-2021"]
say("grid cells on the common window: ", nrow(g))

## ---- presentation choices ----------------------------------------------------
# Three control counts, not six sets: 0, 4 (agriculture only), 6 (agriculture
# plus population and GDP per capita). K1/K3/K5 are dropped from the tables as
# intermediate cases; they are in the appendix grid.
KSETS <- c(K0_none = "0", K2_ag10 = "4", K4_ag10_gdp = "6")
KLAB  <- c(K0_none = "No controls", K2_ag10 = "Agricultural",
           K4_ag10_gdp = "Agricultural + economic")
TSETS <- c(binary = "Binary", count = "Count", log1p = "Log(1+y)")
ZMAIN <- "Z2_max"
ZSETS <- c(Z1_mean = "Mean VPD", Z2_max = "Mean daily-max VPD",
           Z3_exc3 = "Degrees above 3 kPa")
LSETS <- c(L0 = "Contemporaneous", L1 = "Lag 1 alone", L2 = "Lag 2 alone",
           LALL = "Contemporaneous $+$ lags 1, 2")

stars <- function(p) if (is.na(p)) "" else
  if (p < 0.01) "$^{***}$" else if (p < 0.05) "$^{**}$" else
  if (p < 0.10) "$^{*}$" else ""
fmt  <- function(x, d = 4) if (is.na(x)) "" else formatC(x, format = "f", digits = d)
fmtN <- function(x) if (is.na(x)) "" else formatC(round(x), format = "d", big.mark = ",")

# one coefficient / standard-error pair of rows across a set of cells
coef_rows <- function(cells, label) {
  c(sprintf("%s & %s \\\\", label,
            paste(sapply(seq_len(nrow(cells)), function(i)
              paste0(fmt(cells$iv_b[i]), stars(cells$iv_p[i]))), collapse = " & ")),
    sprintf(" & %s \\\\[2pt]",
            paste(sprintf("(%s)", sapply(cells$iv_se, fmt)), collapse = " & ")))
}

write_tex <- function(lines, file) {
  writeLines(lines, file.path(DIR_TAB, file))
  say("wrote ", file)
}

################################################################################
# SET 1 and SET 2 - main results, one table per instrument
################################################################################
# argument is named zid, not zset: inside data.table's `[` a bare `zset` would
# resolve to the COLUMN, so `zset == zset` is trivially true for every row
make_iv_table <- function(zid, num, main = FALSE) {
  sub <- g[zset == zid & lagid == "L0" &
             kset %in% names(KSETS) & transformation %in% names(TSETS)]
  ncol_body <- length(TSETS) * length(KSETS)
  hdr_t <- paste(sprintf("\\multicolumn{%d}{c}{%s}", length(KSETS), TSETS),
                 collapse = " & ")
  cmids <- paste(sapply(seq_along(TSETS), function(j) {
    a <- 2 + (j - 1) * length(KSETS); sprintf("\\cmidrule(lr){%d-%d}", a, a + length(KSETS) - 1)
  }), collapse = " ")
  hdr_k <- paste(rep(unname(KSETS), length(TSETS)), collapse = " & ")

  # 3 transformations x 3 control counts = 9 body columns, which overflows a
  # portrait text block at any readable font size. Rotate rather than shrink.
  L <- c(sprintf("\\begin{sidewaystable}[p]\\centering"),
         sprintf("\\caption{%s}",
                 if (main) sprintf("Effect of area burned on forced-labour outcomes (instrument: %s)", ZSETS[[zid]])
                 else sprintf("Robustness to instrument operationalization: %s", ZSETS[[zid]])),
         sprintf("\\label{tab:iv_%s}", tolower(zid)),
         "\\footnotesize",
         sprintf("\\begin{tabular}{l%s}", paste(rep("c", ncol_body), collapse = "")),
         "\\toprule",
         sprintf(" & %s \\\\", hdr_t), cmids,
         sprintf("Number of controls & %s \\\\", hdr_k),
         "\\midrule")

  for (dvn in names(DVS)) {
    L <- c(L, sprintf("\\multicolumn{%d}{l}{\\emph{Panel %s: %s}} \\\\[2pt]",
                      ncol_body + 1, if (dvn == "TF") "A" else "B", LAB_DV[[dvn]]))
    cells <- rbindlist(lapply(names(TSETS), function(tr)
      rbindlist(lapply(names(KSETS), function(k)
        sub[source == dvn & transformation == tr & kset == k]))))
    L <- c(L, coef_rows(cells, "Area burned (\\% of municipality)"))
    L <- c(L, sprintf("First-stage $F$ & %s \\\\",
                      paste(sapply(cells$fs_F, function(x) fmt(x, 1)), collapse = " & ")),
           sprintf("Observations & %s \\\\",
                   paste(sapply(cells$n_obs, fmtN), collapse = " & ")))
    if (dvn == "TF") L <- c(L, "\\midrule")
  }

  L <- c(L, "\\bottomrule", "\\end{tabular}",
         sprintf("\\begin{minipage}{%s\\textwidth}\\vspace{4pt}\\footnotesize",
                 "0.98"),
         "Two-stage least squares. Area burned is instrumented by the quarterly",
         sprintf("%s. All specifications include municipality and year-quarter fixed effects;",
                 ZSETS[[zid]]),
         "standard errors, in parentheses, are two-way clustered on municipality and year-quarter.",
         "The binary column is a linear probability model for whether the municipality-quarter",
         "records any case; the count column is a level effect; log(1+y) is approximately",
         "proportional. Agricultural controls are cattle density, crop area share, log cattle",
         "herd and log crop value; economic controls add log population and log GDP per capita.",
         "The sample is 2002--2021 in every column, because the economic controls are only",
         "available through 2021 and holding the window fixed keeps the comparison across",
         "columns about controls rather than about sample.",
         "$^{*}p<0.1$; $^{**}p<0.05$; $^{***}p<0.01$.",
         "\\end{minipage}\\end{sidewaystable}")
  L
}

write_tex(make_iv_table(ZMAIN, 1, main = TRUE), "tab1_main_iv_meeting.tex")
k <- 0L
for (z in c("Z1_mean", "Z3_exc3")) {
  k <- k + 1L
  write_tex(make_iv_table(z, k), sprintf("tab2%s_iv_%s_meeting.tex", letters[k], z))
}

################################################################################
# SET 3 - first-stage statistics
################################################################################
# The first stage is percent_area_burned ~ instrument | FE. It does NOT depend on
# which outcome is on the left-hand side, and both DVs use the same complete-case
# sample, so a TF panel and a UE panel would be identical row for row. Reporting
# it once, with that stated, is honest; reporting it twice implies a distinction
# that does not exist. Verified below rather than assumed.
fs <- g[lagid == "L0" & kset %in% names(KSETS) & zset %in% names(ZSETS)]
chk <- dcast(fs[transformation == "binary"], zset + kset ~ source, value.var = "fs_b")
stopifnot(isTRUE(all.equal(chk$TF, chk$UE)))
fs <- fs[source == "TF" & transformation == "binary"]

L <- c("\\begin{table}[htbp]\\centering",
       "\\caption{First stage: vapour pressure deficit and area burned}",
       "\\label{tab:first_stage_meeting}", "\\footnotesize",
       "\\begin{tabular}{llccccc}", "\\toprule",
       "Instrument & Controls & Coefficient & Std. error & $p$ & Joint $F$ & Obs. \\\\",
       "\\midrule")
for (z in names(ZSETS)) {
  for (i in seq_along(KSETS)) {
    kk <- names(KSETS)[i]
    r <- fs[zset == z & kset == kk]
    if (!nrow(r)) next
    L <- c(L, sprintf("%s & %s & %s%s & (%s) & %s & %s & %s \\\\",
                      if (i == 1L) ZSETS[[z]] else "", KLAB[[kk]],
                      fmt(r$fs_b, 4), stars(r$fs_p), fmt(r$fs_se, 4),
                      if (r$fs_p < 0.001) "$<$0.001" else fmt(r$fs_p, 3),
                      fmt(r$fs_F, 1), fmtN(r$n_obs)))
  }
  if (z != names(ZSETS)[length(ZSETS)]) L <- c(L, "\\addlinespace")
}
L <- c(L, "\\bottomrule", "\\end{tabular}",
       "\\begin{minipage}{0.92\\textwidth}\\vspace{4pt}\\footnotesize",
       "Each row regresses area burned on the named instrument with municipality and",
       "year-quarter fixed effects and two-way clustered standard errors. The first stage",
       "does not depend on which outcome is used, and both outcomes are estimated on the",
       "same complete-case sample, so these rows are common to the task force and",
       "seguro-desemprego models rather than being reported twice. The joint $F$ is computed",
       "on the excluded instrument only, not on all covariates. Conventional weak-instrument",
       "thresholds are cleared by a wide margin throughout, including for the",
       "degrees-above-3-kPa measure, which is strong here but weak in the second stage",
       "because it is zero for most municipality-quarters.",
       "$^{*}p<0.1$; $^{**}p<0.05$; $^{***}p<0.01$.",
       "\\end{minipage}\\end{table}")
write_tex(L, "tab3_first_stage_meeting.tex")

################################################################################
# SET 4 - lag structures
################################################################################
lg <- g[zset == ZMAIN & kset == "K4_ag10_gdp" & transformation %in% names(TSETS)]
L <- c("\\begin{table}[htbp]\\centering",
       "\\caption{Lag structure}", "\\label{tab:lags_meeting}", "\\footnotesize",
       "\\begin{tabular}{lccc}", "\\toprule",
       sprintf("Lag structure & %s \\\\", paste(TSETS, collapse = " & ")),
       "\\midrule")
for (dvn in names(DVS)) {
  L <- c(L, sprintf("\\multicolumn{4}{l}{\\emph{Panel %s: %s}} \\\\[2pt]",
                    if (dvn == "TF") "A" else "B", LAB_DV[[dvn]]))
  for (ll in names(LSETS)) {
    cells <- rbindlist(lapply(names(TSETS), function(tr)
      lg[source == dvn & transformation == tr & lagid == ll]))
    if (!nrow(cells)) next
    L <- c(L, coef_rows(cells, LSETS[[ll]]))
  }
  L <- c(L, sprintf("First-stage $F$ (contemporaneous) & \\multicolumn{3}{c}{%s} \\\\",
                    fmt(lg[source == dvn & lagid == "L0" & transformation == "binary", fs_F], 1)))
  if (dvn == "TF") L <- c(L, "\\midrule")
}
L <- c(L, "\\bottomrule", "\\end{tabular}",
       "\\begin{minipage}{0.75\\textwidth}\\vspace{4pt}\\footnotesize",
       "Instrument: quarterly mean of daily-max VPD; six controls; 2002--2021.",
       "In the final row of each panel all three fire terms are entered together and each is",
       "instrumented by its own correspondingly-lagged VPD; the reported coefficient is on",
       "contemporaneous fire. \\textbf{The two middle rows are misspecified by omission} and",
       "are shown only because they are the obvious thing to try: they enter a lagged fire",
       "term without the contemporaneous one, and within municipality after these fixed",
       "effects $\\mathrm{corr}(\\text{fire}_t, \\text{fire}_{t-2}) = -0.22$, so a lag-2 term",
       "proxies for \\emph{low} current fire and changes sign mechanically.",
       "$^{*}p<0.1$; $^{**}p<0.05$; $^{***}p<0.01$.",
       "\\end{minipage}\\end{table}")
write_tex(L, "tab4_lags_meeting.tex")

################################################################################
# numbers quoted in the executive summary, written out so the prose cannot drift
################################################################################
m <- g[zset == ZMAIN & lagid == "L0" & kset == "K4_ag10_gdp" & transformation == "binary"]
sd_fire <- sd(load_panel()[[FIRE]], na.rm = TRUE)
key <- data.table(
  dv = m$source, coef = m$iv_b, se = m$iv_se, p = m$iv_p, F = m$fs_F, n = m$n_obs,
  baseline = BASELINE_RATE[m$source],
  pct_of_baseline_per_sd = 100 * m$iv_b * sd_fire / BASELINE_RATE[m$source])
fwrite(key, file.path(DIR_RES, "meeting_key_numbers.csv"))
cat("\n===== NUMBERS FOR THE EXECUTIVE SUMMARY =====\n"); print(key)
say("sd of area burned: ", round(sd_fire, 3), " pp")

# how often the contemporaneous binary result is significant, across the grid
rob <- g[lagid == "L0" & transformation == "binary",
         .(n = .N, pct_pos = round(100 * mean(iv_b > 0), 1),
           pct_p05 = round(100 * mean(iv_p < 0.05), 1)), by = source]
print(rob); fwrite(rob, file.path(DIR_RES, "meeting_robustness.csv"))
say("done")
