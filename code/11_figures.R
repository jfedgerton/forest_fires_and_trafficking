################################################################################
# 11_figures.R
#
# Every figure except the maps (those are 12_maps.py, which needs geopandas).
#
#   Fig 1  First-stage binscatter. VPD against area burned, both residualised on
#          the fixed effects. The single most persuasive panel in an IV paper and
#          the one this project did not have.
#   Fig 2  Seasonality. Fire, rescues and victims by quarter of year. The first
#          referee objection to a national year-quarter fixed effect is "this is
#          seasonality"; this figure is where that argument is won or lost.
#   Fig 3  Effect of a one-standard-deviation increase in area burned, expressed
#          expressed as a percentage change. Binary and count are a percentage of
#          their own baseline; log(1+y) is an approximate proportional change.
#          These are NOT a common scale, so each panel carries its own axis.
#   Fig 4  Specification curves, faceted by DV and transformation.
#   Fig 5  Specification curves faceted by DV x instrument set, split by whether
#          controls are included.
#
# Fig 3 is labelled an effect size rather than a marginal effect on purpose. In a
# linear probability model the marginal effect IS the coefficient; calling it
# something else invites the reader to think a transformation has been applied.
#
# The curves use the contemporaneous and contemporaneous-plus-lags structures
# only. L1 and L2 enter a lagged fire term without the contemporaneous one and
# flip sign mechanically (corr(fire_t, fire_t-2) = -0.22 within municipality
# after the fixed effects); they are estimated in 09 and excluded here.
################################################################################

source("00_config.R")

d    <- load_panel()
tidy <- fread(file.path(DIR_RES, "main_models_tidy.csv"))
GRID <- file.path(DIR_RES, "spec_grid_2sls.csv")

## =============================================================================
## Fig 1 - first-stage binscatter
## =============================================================================
say("figure 1: first-stage binscatter")
fs_resid <- function(dd) {
  mx <- feols(as.formula(sprintf("%s ~ 1 | %s", VPD_MAIN, FE_MAIN)), data = dd, notes = FALSE)
  my <- feols(as.formula(sprintf("%s ~ 1 | %s", FIRE, FE_MAIN)), data = dd, notes = FALSE)
  data.table(z = resid(mx), y = resid(my))
}
r <- fs_resid(d[stats::complete.cases(d[, c(VPD_MAIN, FIRE), with = FALSE])])
r[, bin := cut(z, breaks = quantile(z, probs = seq(0, 1, length.out = 41), na.rm = TRUE),
               include.lowest = TRUE, labels = FALSE)]
bs <- r[, .(z = mean(z), y = mean(y), n = .N), by = bin][order(z)]
slope <- coef(lm(y ~ z, data = r))[2]

g1 <- ggplot(bs, aes(z, y)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_smooth(data = r, aes(z, y), method = "lm", formula = y ~ x,
              colour = "#B2182B", fill = "#B2182B", alpha = 0.18, linewidth = 0.7) +
  geom_point(size = 2, colour = "grey20") +
  annotate("text", x = min(bs$z), y = max(bs$y), hjust = 0, vjust = 1,
           label = sprintf("slope = %.3f", slope), size = 3.5) +
  labs(title = "First stage: vapour pressure deficit and area burned",
       subtitle = "40 equal-sized bins; both axes residualised on municipality and year-quarter fixed effects",
       x = "Mean daily-max VPD (kPa), residualised",
       y = "Area burned (% of municipality), residualised") +
  theme_ff()
save_fig(g1, "fig1_first_stage_binscatter", 8, 6)

## =============================================================================
## Fig 2 - seasonality
## =============================================================================
say("figure 2: seasonality")
seas <- d[, .(`Area burned (%)` = mean(percent_area_burned, na.rm = TRUE),
              `Mean daily-max VPD (kPa)` = mean(vpd_max_q, na.rm = TRUE),
              `TF: any rescue (%)` = 100 * mean(tf_count > 0),
              `UE: any victim (%)` = 100 * mean(ue_count > 0)), by = quarter]
seas <- melt(seas, id.vars = "quarter", variable.name = "series", value.name = "value")
seas[, quarter_lab := factor(paste0("Q", quarter), levels = paste0("Q", 1:4))]

g2 <- ggplot(seas, aes(quarter_lab, value, group = 1)) +
  geom_col(fill = "grey75", width = 0.65) +
  geom_point(colour = "#B2182B", size = 1.8) +
  geom_line(colour = "#B2182B", linewidth = 0.5) +
  facet_wrap(~ series, scales = "free_y", nrow = 1) +
  labs(title = "Seasonality of fire, the instrument, and both outcomes",
       subtitle = "Municipality-quarter means, 2002-2024. Fire peaks in Q3; the outcomes do not track it mechanically.",
       x = NULL, y = NULL) +
  theme_ff()
save_fig(g2, "fig2_seasonality", 11, 4)

## =============================================================================
## Fig 3 - effect of a one-SD increase in area burned
## =============================================================================
say("figure 3: one-SD effect sizes")
sd_fire <- sd(d[[FIRE]], na.rm = TRUE)
BASE_LEVEL <- c(TF = mean(d$tf_count), UE = mean(d$ue_count))

eff <- copy(tidy)[transformation %in% c("binary", "count", "log1p")]
eff[, `:=`(b_sd = iv_b * sd_fire, lo_sd = (iv_b - 1.96 * iv_se) * sd_fire,
           hi_sd = (iv_b + 1.96 * iv_se) * sd_fire)]
# put all three on one axis: percent change relative to that outcome's baseline
eff[, denom := fifelse(transformation == "binary", BASELINE_RATE[dv],
                fifelse(transformation == "count", BASE_LEVEL[dv], 1))]
eff[, `:=`(pct = 100 * b_sd / denom, pct_lo = 100 * lo_sd / denom,
           pct_hi = 100 * hi_sd / denom)]
# log(1+y) is already proportional; express it as a percentage directly
eff[transformation == "log1p", `:=`(pct = 100 * b_sd, pct_lo = 100 * lo_sd,
                                    pct_hi = 100 * hi_sd)]
eff[, trans_lab := factor(LAB_TRAN[transformation], levels = unname(LAB_TRAN))]
eff[, dv_lab := factor(LAB_DV[dv], levels = unname(LAB_DV))]
eff[, ctrl_lab := factor(c(none_fullwindow = "No controls\n2002-2024",
                           none = "No controls\n2002-2021",
                           full = "Full controls\n2002-2021")[controls],
                         levels = c("No controls\n2002-2024", "No controls\n2002-2021",
                                    "Full controls\n2002-2021"))]

# facet_wrap, not facet_grid: facet_grid ties the y scale within a row, which
# crushes the log(1+y) panel against zero. The three transformations are NOT on a
# common scale - binary and count are percentages of their own baseline, log(1+y)
# is an approximate proportional change - so each panel gets its own axis and its
# own unit label rather than a shared one that flatters the comparison.
eff[, panel := factor(paste0(dv_lab, "\n", trans_lab))]
g3 <- ggplot(eff, aes(ctrl_lab, pct, colour = dv_lab)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_pointrange(aes(ymin = pct_lo, ymax = pct_hi),
                  linewidth = 0.6, size = 0.5) +
  facet_wrap(~ panel, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = unname(COL_DV), guide = "none") +
  labs(title = "Effect of a one-standard-deviation increase in area burned",
       subtitle = sprintf(paste0("One SD = %.2f pp of municipal area. Binary and count are a percentage of ",
                                 "that outcome's own baseline;\nlog(1+y) is an approximate proportional ",
                                 "change. Scales differ by panel; 95%% confidence intervals."), sd_fire),
       x = NULL, y = "Percent change") +
  theme_ff()
save_fig(g3, "fig3_one_sd_effects", 11, 7)
fwrite(eff[, .(dv, transformation, controls, window, iv_b, iv_se,
               effect_per_sd = b_sd, pct_change = pct, ci_lo = pct_lo, ci_hi = pct_hi)],
       file.path(DIR_RES, "one_sd_effects.csv"))

## =============================================================================
## Figs 4 and 5 - specification curves
## =============================================================================
if (!file.exists(GRID)) {
  say("spec grid not found - run 09_specification_grid.R for figures 4 and 5")
} else {
  res <- fread(GRID)
  s <- res[status == "ok" & lag_only == 0]
  s[, trans_lab := factor(LAB_TRAN[transformation], levels = unname(LAB_TRAN))]
  s[, dv_lab := factor(LAB_DV[source], levels = unname(LAB_DV))]
  s[, inst_lab := factor(LAB_Z[zset], levels = unname(LAB_Z))]
  s[, lag_lab := factor(LAB_LAG[lagid], levels = unname(LAB_LAG))]
  s[, has_ctrl := factor(fifelse(kset == "K0_none", "No controls", "With controls"),
                         levels = c("No controls", "With controls"))]
  s[, sig := factor(fifelse(iv_p < 0.05, "p < 0.05",
                     fifelse(iv_p < 0.10, "0.05 <= p < 0.10", "p >= 0.10")),
                    levels = names(COL_SIG))]
  s[, `:=`(lo = iv_b - 1.96 * iv_se, hi = iv_b + 1.96 * iv_se)]
  say("spec curve specifications: ", nrow(s))

  curve <- function(dd, facet, stub, ttl, w = 11, h = 8) {
    dd <- copy(dd)[order(iv_b)]
    dd[, x := seq_len(.N), by = eval(all.vars(facet))]
    med <- dd[, .(m = median(iv_b, na.rm = TRUE)), by = eval(all.vars(facet))]
    g <- ggplot(dd, aes(x, iv_b)) +
      geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.4) +
      geom_hline(data = med, aes(yintercept = m), linetype = 2,
                 colour = "#2166AC", linewidth = 0.4) +
      geom_linerange(aes(ymin = lo, ymax = hi, colour = sig),
                     linewidth = 0.3, alpha = 0.6) +
      geom_point(aes(colour = sig), size = 0.7) +
      facet_grid(facet, scales = "free") +
      scale_colour_manual(values = COL_SIG, name = NULL, drop = FALSE) +
      labs(title = ttl,
           subtitle = "Each point is one specification; 95% CI; blue dashed line is the panel median",
           x = "Specification, ordered by estimate within panel", y = "IV estimate") +
      theme_ff() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    save_fig(g, stub, w, h)
  }

  curve(s[lagid == "L0"], trans_lab ~ dv_lab, "fig4_speccurve_by_outcome",
        "Specification curve, contemporaneous fire, by outcome and DV")
  curve(s, trans_lab ~ dv_lab + lag_lab, "fig4b_speccurve_by_outcome_and_lag",
        "Specification curve by outcome, DV and lag structure", 13, 8)
  curve(s[lagid == "L0" & transformation == "binary"], inst_lab ~ dv_lab + has_ctrl,
        "fig5_speccurve_by_instrument_controls",
        "Specification curve by instrument set, with and without controls (binary outcome)",
        11, 10)
  curve(s[lagid == "L0" & transformation == "count"], inst_lab ~ dv_lab + has_ctrl,
        "fig5b_speccurve_by_instrument_controls_count",
        "Specification curve by instrument set, with and without controls (count outcome)",
        11, 10)

  tab <- s[lagid == "L0", .(n = .N, median = round(median(iv_b, na.rm = TRUE), 5),
                            pct_pos = round(100 * mean(iv_b > 0, na.rm = TRUE), 1),
                            pct_p05 = round(100 * mean(iv_p < 0.05, na.rm = TRUE), 1),
                            median_F = round(median(fs_F, na.rm = TRUE), 1)),
           by = .(dv_lab, trans_lab, inst_lab, has_ctrl)]
  fwrite(tab, file.path(DIR_RES, "speccurve_summary.csv"))
  cat("\n===== SPEC CURVE SUMMARY (contemporaneous) =====\n")
  print(tab[order(dv_lab, trans_lab, inst_lab, has_ctrl)])
}

say("figures written to ", DIR_FIG)
