################################################################################
# 16_robustness_sweep.R
#
# Three sweeps over functional form and fire measurement.
#
#   SWEEP 1  OLS only. Every fire measure x raw and log(1+x) x every outcome
#            transformation x with and without the full control set.
#   SWEEP 2  IV. The same, crossed with all seven instrument sets, and with the
#            SECOND-STAGE treatment entered either raw or as log(1+fire) - in
#            which case the logged treatment is what gets instrumented.
#   SWEEP 3  Is folded into 1 and 2 rather than run separately: `fire_measure`
#            is a dimension throughout, so every cell exists for all three
#            measures with and without logging.
#
# THE THREE FIRE MEASURES, and why they are not interchangeable:
#   percent_area_burned         no-double-count area / municipal area x 100.
#                               Scale-free. A 1-unit change means the same thing
#                               in Altamira and in a small Sao Paulo municipality.
#   area_burned_no_double_sqkm  km2 burned, a cell counted once per quarter no
#                               matter how many of the three months it burned.
#   total_burned_area_sqkm      km2 summed over the quarter's three monthly
#                               grids, so a cell that burns twice counts twice.
#                               This is the REPEAT-FIRE measure.
#
# The two km2 measures correlate at 0.997, so repeat burning is rare enough that
# they will rarely diverge; the useful contrast is either of them against the
# percentage. Both are levels and therefore confounded with municipality size
# (r = 0.20 with municipal area), which the municipality fixed effect absorbs in
# levels but not in logs - worth remembering when reading the logged rows.
#
# WHY LOGGING IS ON THE MENU AT ALL. percent_area_burned is 73.4% exact zeros
# with skewness 5.8; a linear-in-levels treatment fits one slope dominated by the
# top few percent. Logging also strengthens the first stage substantially
# (F rises from ~102 to ~272), which is a reason to prefer it that does not
# depend on what it does to any second-stage result.
#
# WINDOW IS FIXED AT 2002-2021 THROUGHOUT. The full control set draws on SIDRA
# 5938, which ends in 2021. Letting the window move with the control set would
# confound "controls change the answer" with "the last three years change the
# answer", which in this project was large enough to flip a sign.
#
# Output: output/results/robustness_sweep.csv
#         output/tables/tab_sweep{1_ols,2_iv}.tex
# Seed 123 (set in 00_config.R; nothing here draws random numbers).
################################################################################

source("00_config.R")

d <- load_panel()[year <= TRUNC_YEAR]
CTRL_FULL <- build_control_sets()$K4_ag10_gdp
need <- c(unname(DVS), FIRE, VPD_ALL, CTRL_FULL,
          "total_burned_area_sqkm", "area_burned_no_double_sqkm")
d <- d[stats::complete.cases(d[, ..need])]
say("sample: ", nrow(d), " rows, ", uniqueN(d$code_muni_7), " municipalities")

FIRE_MEASURES <- c(pct   = "percent_area_burned",
                   nodbl = "area_burned_no_double_sqkm",
                   total = "total_burned_area_sqkm")
FIRE_LAB <- c(pct = "\\% of area burned", nodbl = "km$^2$ burned (no double count)",
              total = "km$^2$ burned (repeat fire counted)")
FIRE_TRANS <- c(raw = "Raw", log = "Log(1+x)")
OUT_TRANS  <- c(binary = "Binary", log1p = "Log(1+y)", count = "Raw count")
CTRL_OPT   <- list(none = character(0), full = CTRL_FULL)
CTRL_LAB   <- c(none = "No controls", full = "Full controls")

# build every transformed treatment column once
for (fm in names(FIRE_MEASURES)) {
  v <- FIRE_MEASURES[[fm]]
  d[, (paste0("T_", fm, "_raw")) := get(v)]
  d[, (paste0("T_", fm, "_log")) := log1p(get(v))]
}

rows <- list(); k <- 0L
t0 <- Sys.time()

################################################################################
# SWEEP 1 - OLS
################################################################################
say("sweep 1: OLS")
for (dvn in names(DVS)) for (ot in names(OUT_TRANS)) for (fm in names(FIRE_MEASURES))
  for (ft in names(FIRE_TRANS)) for (cn in names(CTRL_OPT)) {
    yv <- DVS[[dvn]]; tv <- paste0("T_", fm, "_", ft); ctrl <- CTRL_OPT[[cn]]
    dd <- apply_transform(copy(d), yv, ot)
    m <- safely(feols(as.formula(sprintf(".y ~ %s + %s | %s", tv, rhs(ctrl), FE_MAIN)),
                      data = dd, cluster = CLUSTER, notes = FALSE))
    s <- pick(m, tv)
    k <- k + 1L
    rows[[k]] <- data.table(sweep = "1_OLS", estimator = "OLS", dv = dvn,
      out_trans = ot, fire_measure = fm, fire_trans = ft, controls = cn,
      zset = NA_character_, n_obs = one(tryCatch(nobs(m), error = function(e) NA)),
      b = one(s[1]), se = one(s[2]), p = one(s[4]), fs_F = NA_real_)
  }

################################################################################
# SWEEP 2 and 3 - IV, all instrument sets, both treatment transforms,
#                  all three fire measures
################################################################################
say("sweep 2/3: IV")
for (dvn in names(DVS)) for (ot in names(OUT_TRANS)) for (fm in names(FIRE_MEASURES))
  for (ft in names(FIRE_TRANS)) for (cn in names(CTRL_OPT)) for (zn in names(INSTSETS)) {
    yv <- DVS[[dvn]]; tv <- paste0("T_", fm, "_", ft)
    ctrl <- CTRL_OPT[[cn]]; zs <- INSTSETS[[zn]]
    dd <- apply_transform(copy(d), yv, ot)
    m_fs <- safely(feols(as.formula(sprintf("%s ~ %s + %s | %s", tv, rhs(zs), rhs(ctrl), FE_MAIN)),
                         data = dd, cluster = CLUSTER, notes = FALSE))
    m_iv <- safely(feols(as.formula(sprintf(".y ~ %s | %s | %s ~ %s",
                                            rhs(ctrl), FE_MAIN, tv, rhs(zs))),
                         data = dd, cluster = CLUSTER, notes = FALSE))
    s <- pick(m_iv, paste0("fit_", tv))
    k <- k + 1L
    rows[[k]] <- data.table(sweep = "2_IV", estimator = "2SLS", dv = dvn,
      out_trans = ot, fire_measure = fm, fire_trans = ft, controls = cn,
      zset = zn, n_obs = one(tryCatch(nobs(m_iv), error = function(e) NA)),
      b = one(s[1]), se = one(s[2]), p = one(s[4]),
      fs_F = first_stage_F(m_fs, zs))
    if (k %% 100 == 0) say(sprintf("  %d models, %.1f min", k,
                                   as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

r <- rbindlist(rows, fill = TRUE)
r[, sig05 := as.integer(p < 0.05)]
fwrite(r, file.path(DIR_RES, "robustness_sweep.csv"))
say("models estimated: ", nrow(r))
stopifnot(nrow(r) == 2*3*3*2*2 + 2*3*3*2*2*length(INSTSETS))
stopifnot(mean(is.finite(r$b)) > 0.95)

################################################################################
# READOUTS
################################################################################
cat("\n============ SWEEP 1: OLS, EVERY TRANSFORM ============\n")
o <- r[sweep == "1_OLS"]
print(dcast(o, dv + fire_measure + fire_trans + controls ~ out_trans,
            value.var = "p")[order(dv, fire_measure, fire_trans, controls)])
cat("\n-- OLS: how often significant, by treatment transform --\n")
print(o[, .(n = .N, pct_sig = round(100*mean(sig05), 1),
            median_b = signif(median(b, na.rm = TRUE), 3)),
        by = .(dv, fire_trans)][order(dv, fire_trans)])

cat("\n============ SWEEP 2: IV ============\n")
iv <- r[sweep == "2_IV"]
print(iv[, .(n = .N, pct_pos = round(100*mean(b > 0, na.rm = TRUE), 1),
             pct_sig = round(100*mean(sig05, na.rm = TRUE), 1),
             median_F = round(median(fs_F, na.rm = TRUE), 1)),
         by = .(dv, fire_trans, out_trans)][order(dv, fire_trans, out_trans)])

cat("\n-- IV by instrument set --\n")
print(iv[, .(n = .N, pct_sig = round(100*mean(sig05, na.rm = TRUE), 1),
             median_F = round(median(fs_F, na.rm = TRUE), 1)),
         by = .(zset, fire_trans)][order(zset, fire_trans)])

cat("\n============ SWEEP 3: FIRE MEASURE ============\n")
print(r[, .(n = .N, pct_sig = round(100*mean(sig05, na.rm = TRUE), 1),
            median_b = signif(median(b, na.rm = TRUE), 3)),
        by = .(estimator, fire_measure, fire_trans)][order(estimator, fire_measure, fire_trans)])

cat("\n-- does repeat-fire counting matter? total vs no-double, same cells --\n")
cmp <- merge(r[fire_measure == "total", .(sweep,dv,out_trans,fire_trans,controls,zset,b_total=b,p_total=p)],
             r[fire_measure == "nodbl", .(sweep,dv,out_trans,fire_trans,controls,zset,b_nodbl=b,p_nodbl=p)],
             by = c("sweep","dv","out_trans","fire_trans","controls","zset"))
cmp[, rel_diff := (b_total - b_nodbl)/abs(b_nodbl)]
print(cmp[, .(n = .N, median_rel_diff = round(median(rel_diff, na.rm = TRUE), 4),
              max_abs_rel_diff = round(max(abs(rel_diff), na.rm = TRUE), 3),
              sign_flips = sum(sign(b_total) != sign(b_nodbl), na.rm = TRUE),
              sig_disagree = sum((p_total<0.05) != (p_nodbl<0.05), na.rm = TRUE)),
          by = .(sweep, fire_trans)])

################################################################################
# LATEX
################################################################################
stars <- function(p) if (is.na(p)) "" else if (p<.01) "$^{***}$" else
  if (p<.05) "$^{**}$" else if (p<.10) "$^{*}$" else ""
fm4 <- function(x) if (is.na(x)) "" else formatC(x, format="g", digits=3)

mk <- function(dt, caption, label, file, showF = TRUE) {
  L <- c("\\begin{table}[htbp]\\centering", sprintf("\\caption{%s}", caption),
         sprintf("\\label{%s}", label), "\\footnotesize",
         sprintf("\\begin{tabular}{lll%s}", paste(rep("c", length(OUT_TRANS)), collapse="")),
         "\\toprule",
         sprintf("Fire measure & Transform & Controls & %s \\\\",
                 paste(OUT_TRANS, collapse = " & ")), "\\midrule")
  for (dvn in names(DVS)) {
    L <- c(L, sprintf("\\multicolumn{%d}{l}{\\emph{Panel %s: %s}} \\\\[2pt]",
                      3 + length(OUT_TRANS), if (dvn=="TF") "A" else "B", LAB_DV[[dvn]]))
    for (fm in names(FIRE_MEASURES)) for (ft in names(FIRE_TRANS)) for (cn in names(CTRL_OPT)) {
      cells <- rbindlist(lapply(names(OUT_TRANS), function(ot)
        dt[dv==dvn & fire_measure==fm & fire_trans==ft & controls==cn & out_trans==ot]))
      if (!nrow(cells)) next
      L <- c(L, sprintf("%s & %s & %s & %s \\\\",
                        FIRE_LAB[[fm]], FIRE_TRANS[[ft]], CTRL_LAB[[cn]],
                        paste(sapply(seq_len(nrow(cells)), function(i)
                          paste0(fm4(cells$b[i]), stars(cells$p[i]))), collapse=" & ")),
             sprintf(" & & & %s \\\\[1pt]",
                     paste(sprintf("(%s)", sapply(cells$se, fm4)), collapse=" & ")))
      if (showF) L <- c(L, sprintf(" & & \\footnotesize $F$ & %s \\\\[2pt]",
                        paste(sapply(cells$fs_F, function(x) fm4(round(x,1))), collapse=" & ")))
    }
    if (dvn=="TF") L <- c(L, "\\midrule")
  }
  c(L, "\\bottomrule", "\\end{tabular}",
    "\\begin{minipage}{0.95\\textwidth}\\vspace{4pt}\\footnotesize",
    "Municipality and year-quarter fixed effects; standard errors two-way clustered.",
    "Sample 2002--2021 throughout so that the with- and without-controls rows differ by",
    "controls alone. Full controls are cattle density, crop area share, log cattle herd,",
    "log crop value, log population and log GDP per capita.",
    "$^{*}p<0.1$; $^{**}p<0.05$; $^{***}p<0.01$.",
    "\\end{minipage}\\end{table}")
}

writeLines(mk(r[sweep=="1_OLS"], "OLS across fire measures, transformations and controls",
              "tab:sweep_ols", showF = FALSE),
           file.path(DIR_TAB, "tab_sweep1_ols.tex"))
writeLines(mk(r[sweep=="2_IV" & zset=="Z2_max"],
              "Two-stage least squares across fire measures, transformations and controls (instrument: mean daily-max VPD)",
              "tab:sweep_iv"),
           file.path(DIR_TAB, "tab_sweep2_iv.tex"))
say("wrote tab_sweep1_ols.tex, tab_sweep2_iv.tex")
say("done")
