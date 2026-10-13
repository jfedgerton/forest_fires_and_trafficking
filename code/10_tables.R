################################################################################
# 10_tables.R
#
# Main-text tables. Reads the fitted objects from 08 and writes LaTeX (for the
# manuscript) and CSV (for reading and for the appendix workbook).
#
#   Table 1  First stage. percent_area_burned ~ vpd_max_q | muni + yq,
#            estimated separately in the TF and UE samples and with and without
#            the full control set. Reports the Kleibergen-Paap-style joint F on
#            the excluded instrument.
#
#   Table 2  Second stage, TF. Binary / count / log(1+y), each without controls
#            and with the full control set.
#   Table 3  Second stage, UE. Same layout.
#
#   Table 4  Instrument balance. Each candidate control regressed on
#            vpd_max_q | muni + yq. This is the exclusion-restriction evidence:
#            if the instrument moves a control, conditioning on that control
#            blocks part of the pathway the IV runs through.
#
# A WINDOW NOTE that belongs in the manuscript, not only here. The full control
# set draws on SIDRA 5938, which ends in 2021. Every with-controls column is
# therefore 2002-2021. Each table carries BOTH a no-controls column on the full
# 2002-2024 window and a no-controls column truncated to 2002-2021, so the
# with-versus-without comparison is like-for-like and the reader can see the
# window's contribution separately from the controls'.
#
# Output: output/tables/*.tex, output/tables/*.csv
################################################################################

source("00_config.R")

models <- readRDS(file.path(DIR_RES, "main_models.rds"))
tidy   <- fread(file.path(DIR_RES, "main_models_tidy.csv"))
d      <- load_panel()
CONTROLS  <- build_control_sets()
FULL_CTRL <- CONTROLS$K4_ag10_gdp

DICT <- c(
  percent_area_burned      = "Area burned (\\% of municipality)",
  fit_percent_area_burned  = "Area burned (\\% of municipality)",
  vpd_max_q                = "Mean daily-max VPD (kPa)",
  log_pop                  = "Log population",
  log_gdp_pc               = "Log GDP per capita",
  log_cattle               = "Log cattle herd",
  cattle_density           = "Cattle per km$^2$",
  crop_area_share          = "Crop area share",
  log_crop_value           = "Log crop value",
  .y                       = "Outcome",
  muni                     = "Municipality",
  yq                       = "Year-quarter"
)
setFixest_dict(DICT)

CTRL_ORDER <- c("none_fullwindow", "none", "full")
CTRL_LAB   <- c(none_fullwindow = "No controls, 2002-2024",
                none            = "No controls, 2002-2021",
                full            = "Full controls, 2002-2021")

## =============================================================================
## Table 1 - first stage
## =============================================================================
fs_models <- list(); fs_names <- character(0)
for (dvn in names(DVS)) for (cn in CTRL_ORDER) {
  m <- models[[sprintf("%s__binary__%s", dvn, cn)]]$fs   # first stage is common
  fs_models <- c(fs_models, list(m))                      # across transformations
  fs_names  <- c(fs_names, sprintf("%s: %s", dvn, CTRL_LAB[[cn]]))
}
names(fs_models) <- fs_names

etable(fs_models, tex = TRUE, file = file.path(DIR_TAB, "table1_first_stage.tex"),
       replace = TRUE, title = "First stage: vapour pressure deficit and area burned",
       label = "tab:first_stage",
       fitstat = ~ n + r2 + wr2,
       notes = paste("Municipality-quarter panel, Brazil.",
                     "Standard errors two-way clustered on municipality and year-quarter.",
                     "Full controls are the agriculture set selected at p<.10 plus log population",
                     "and log GDP per capita; that set is drawn from SIDRA 5938 and is only",
                     "available through 2021."))

t1 <- tidy[transformation == "binary",
           .(dv, controls, window, n_obs,
             first_stage_coef = round(fs_b, 5), se = round(fs_se, 5),
             p = round(fs_p, 4), joint_F = round(fs_F, 1))]
t1[, controls := CTRL_LAB[controls]]
fwrite(t1, file.path(DIR_TAB, "table1_first_stage.csv"))
cat("\n===== TABLE 1: FIRST STAGE =====\n"); print(t1)

## =============================================================================
## Tables 2 and 3 - second stage, one per DV
## =============================================================================
for (dvn in names(DVS)) {
  iv_models <- list(); iv_names <- character(0)
  for (trans in c("binary", "count", "log1p")) for (cn in CTRL_ORDER) {
    iv_models <- c(iv_models, list(models[[sprintf("%s__%s__%s", dvn, trans, cn)]]$iv))
    iv_names  <- c(iv_names, sprintf("%s | %s", LAB_TRAN[[trans]], CTRL_LAB[[cn]]))
  }
  names(iv_models) <- iv_names
  num <- if (dvn == "TF") 2 else 3
  etable(iv_models, tex = TRUE,
         file = file.path(DIR_TAB, sprintf("table%d_second_stage_%s.tex", num, dvn)),
         replace = TRUE,
         title = sprintf("Second stage: effect of fire on %s", LAB_DV[[dvn]]),
         label = sprintf("tab:second_stage_%s", tolower(dvn)),
         fitstat = ~ n + ivf,
         notes = paste("Two-stage least squares. Area burned instrumented by the quarterly",
                       "mean of daily-max vapour pressure deficit. Municipality and",
                       "year-quarter fixed effects; standard errors two-way clustered.",
                       "The binary column is a linear probability model; the count column is",
                       "a level effect; log(1+y) is approximately proportional."))
}

t23 <- tidy[, .(dv, transformation, controls, window, n_obs,
                iv_coef = round(iv_b, 5), iv_se = round(iv_se, 5), iv_p = round(iv_p, 4),
                reduced_form = round(rf_b, 5), rf_p = round(rf_p, 4),
                ols = round(ols_b, 5), ols_p = round(ols_p, 4),
                first_stage_F = round(fs_F, 1))]
t23[, controls := CTRL_LAB[controls]]
t23[, transformation := LAB_TRAN[transformation]]
fwrite(t23, file.path(DIR_TAB, "table23_second_stage.csv"))
cat("\n===== TABLES 2-3: SECOND STAGE =====\n"); print(t23)

## ---- effect sizes, expressed as a percentage of baseline ---------------------
sd_fire <- sd(d[[FIRE]], na.rm = TRUE)
eff <- tidy[transformation == "binary",
            .(dv, controls,
              per_pp = iv_b,
              per_sd = iv_b * sd_fire,
              baseline = BASELINE_RATE[dv],
              pct_of_baseline_per_sd = round(100 * iv_b * sd_fire / BASELINE_RATE[dv], 1))]
eff[, controls := CTRL_LAB[controls]]
fwrite(eff, file.path(DIR_TAB, "table_effect_sizes.csv"))
cat("\n===== EFFECT SIZES (binary outcome; SD of area burned = ",
    round(sd_fire, 3), " pp) =====\n", sep = ""); print(eff)

## =============================================================================
## Table 4 - instrument balance on the controls
## =============================================================================
# Each candidate control on the instrument, same fixed effects and clustering as
# the main models. A significant coefficient is a warning about that control,
# not a result about fire.
bal <- rbindlist(lapply(names(CANDIDATE_CONTROLS), function(v) {
  dd <- d[stats::complete.cases(d[, c(v, VPD_MAIN), with = FALSE])]
  if (nrow(dd) < 1000L) return(NULL)
  m <- safely(feols(as.formula(sprintf("%s ~ %s | %s", v, VPD_MAIN, FE_MAIN)),
                    data = dd, cluster = CLUSTER, notes = FALSE))
  s <- pick(m, VPD_MAIN)
  data.table(variable = v, label = CANDIDATE_CONTROLS[[v]],
             block = if (v %in% AG_CANDIDATES) "Agricultural" else "Economic",
             n_obs = one(nobs(m)),
             coef = one(s[1]), se = one(s[2]), t = one(s[3]), p = one(s[4]),
             # standardised so magnitudes are comparable across differently scaled vars
             coef_sd = one(s[1]) * sd(dd[[VPD_MAIN]], na.rm = TRUE) / sd(dd[[v]], na.rm = TRUE))
}))
bal[, sig_05 := as.integer(p < 0.05)]
bal <- bal[order(block, p)]
fwrite(bal, file.path(DIR_TAB, "table4_instrument_balance.csv"))

# LaTeX by hand: etable cannot stack models with different dependent variables
tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Instrument balance: does the instrument move the controls?}",
  "\\label{tab:instrument_balance}",
  "\\begin{tabular}{llrrrr}",
  "\\hline\\hline",
  "Block & Candidate control & Coef. & Std. err. & $p$ & SD per SD \\\\",
  "\\hline")
for (b in unique(bal$block)) {
  sub <- bal[block == b]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i]
    stars <- if (is.na(r$p)) "" else if (r$p < 0.01) "$^{***}$" else
             if (r$p < 0.05) "$^{**}$" else if (r$p < 0.10) "$^{*}$" else ""
    tex <- c(tex, sprintf("%s & %s & %.4f%s & %.4f & %.3f & %.4f \\\\",
                          if (i == 1L) b else "", r$label, r$coef, stars,
                          r$se, r$p, r$coef_sd))
  }
  tex <- c(tex, "\\hline")
}
tex <- c(tex, "\\end{tabular}",
         "\\begin{minipage}{\\textwidth}\\footnotesize",
         "Each row is a separate regression of the named control on the quarterly mean of",
         "daily-max vapour pressure deficit, with municipality and year-quarter fixed effects",
         "and standard errors two-way clustered. The final column rescales the coefficient to",
         "standard deviations of the control per standard deviation of the instrument.",
         sprintf("%d of %d rows are significant at the 5\\%% level.",
                 sum(bal$sig_05, na.rm = TRUE), nrow(bal)),
         "$^{*}p<0.1$; $^{**}p<0.05$; $^{***}p<0.01$.",
         "\\end{minipage}\\end{table}")
writeLines(tex, file.path(DIR_TAB, "table4_instrument_balance.tex"))

cat("\n===== TABLE 4: INSTRUMENT BALANCE =====\n")
print(bal[, .(block, label, coef = round(coef, 5), se = round(se, 5),
              p = round(p, 4), coef_sd = round(coef_sd, 4), sig_05)])
say(sprintf("balance: %d of %d controls significant at 5%% (chance expectation %.1f)",
            sum(bal$sig_05, na.rm = TRUE), nrow(bal), 0.05 * nrow(bal)))

say("wrote ", length(list.files(DIR_TAB)), " files to ", DIR_TAB)
