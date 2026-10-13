################################################################################
# 19_figures_deck.R
#
# Slide-sized variants of two figures from 11_figures.R.
#
# The paper figures are drawn on an 11 x 7 inch canvas and placed on a slide at
# about 5.5 inches wide. That is a 0.5 scale factor, so a 10 pt axis label lands
# on the projector at 5 pt and nobody in the room can read it. These variants are
# the same plots on a small canvas with a large base size, so that after the
# slide scales them down the text is still 9-10 pt on screen.
#
# Nothing here is a different result. If a number moves between these and the
# paper figures, something is wrong.
#
# Output: output/figures/deck_fig3_one_sd.{pdf,png}
#         output/figures/deck_fig4_speccurve.{pdf,png}
################################################################################

source("00_config.R")

d    <- load_panel()
tidy <- fread(file.path(DIR_RES, "main_models_tidy.csv"))

LAB_TRAN_D <- c(binary = "Binary", count = "Count", log1p = "Log(1+y)")
LAB_DV_D   <- c(TF = "TF - rescues", UE = "UE - victims")

## =============================================================================
## deck fig 3 - one-SD effects
## =============================================================================
say("deck figure 3")
sd_fire <- sd(d[[FIRE]], na.rm = TRUE)
BASE_LEVEL <- c(TF = mean(d$tf_count), UE = mean(d$ue_count))

eff <- copy(tidy)[transformation %in% c("binary", "count", "log1p")]
eff[, `:=`(b_sd = iv_b * sd_fire, lo_sd = (iv_b - 1.96 * iv_se) * sd_fire,
           hi_sd = (iv_b + 1.96 * iv_se) * sd_fire)]
eff[, denom := fifelse(transformation == "binary", BASELINE_RATE[dv],
                fifelse(transformation == "count", BASE_LEVEL[dv], 1))]
eff[, `:=`(pct = 100 * b_sd / denom, pct_lo = 100 * lo_sd / denom,
           pct_hi = 100 * hi_sd / denom)]
eff[transformation == "log1p", `:=`(pct = 100 * b_sd, pct_lo = 100 * lo_sd,
                                    pct_hi = 100 * hi_sd)]
eff[, trans_lab := factor(LAB_TRAN_D[transformation], levels = unname(LAB_TRAN_D))]
eff[, dv_lab    := factor(LAB_DV_D[dv], levels = unname(LAB_DV_D))]
# Short two-line labels, horizontal, aligned against the plot. The long form
# ("No controls / 2002-2024") collides with its neighbours once six panels are
# squeezed onto half a slide.
eff[, ctrl_lab := factor(c(none_fullwindow = "No ctrl\n02-24",
                           none            = "No ctrl\n02-21",
                           full            = "Full ctrl\n02-21")[controls],
                         levels = c("No ctrl\n02-24", "No ctrl\n02-21",
                                    "Full ctrl\n02-21"))]
eff[, panel := factor(paste0(dv_lab, "  |  ", trans_lab))]

g <- ggplot(eff, aes(ctrl_lab, pct, colour = dv_lab)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_pointrange(aes(ymin = pct_lo, ymax = pct_hi), linewidth = 0.7, size = 0.5) +
  facet_wrap(~ panel, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = unname(COL_DV), guide = "none") +
  labs(title = "One SD more burned area",
       subtitle = sprintf(paste0("1 SD = %.1f pp of municipal area. Binary and count as a percent\n",
                                 "of their own baseline; log(1+y) as a proportional change.\n",
                                 "Free y scale by panel; 95%% CI."), sd_fire),
       x = NULL, y = "Percent change") +
  theme_ff(13) + theme(axis.text.x = element_text(size = 9.5))
save_fig(g, "deck_fig3_one_sd", width = 7.4, height = 5.2)

## =============================================================================
## deck fig 4 - specification curve, contemporaneous
## =============================================================================
say("deck figure 4")
GRID <- file.path(DIR_RES, "spec_grid_2sls.csv")
stopifnot(file.exists(GRID))
s <- fread(GRID)[status == "ok" & lag_only == 0 & lagid == "L0"]
s[, trans_lab := factor(LAB_TRAN_D[transformation], levels = unname(LAB_TRAN_D))]
s <- s[!is.na(trans_lab)]        # drops the intensity form, which has its own panel
s[, dv_lab := factor(LAB_DV_D[source], levels = unname(LAB_DV_D))]
s[, sig := factor(fifelse(iv_p < 0.05, "p < 0.05",
                   fifelse(iv_p < 0.10, "0.05 <= p < 0.10", "p >= 0.10")),
                  levels = names(COL_SIG))]
s[, `:=`(lo = iv_b - 1.96 * iv_se, hi = iv_b + 1.96 * iv_se)]
setorder(s, iv_b)
s[, x := seq_len(.N), by = .(trans_lab, dv_lab)]
med <- s[, .(m = median(iv_b, na.rm = TRUE)), by = .(trans_lab, dv_lab)]
say("deck spec curve specifications: ", nrow(s))

g4 <- ggplot(s, aes(x, iv_b)) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.4) +
  geom_hline(data = med, aes(yintercept = m), linetype = 2,
             colour = "#2166AC", linewidth = 0.45) +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = sig), linewidth = 0.3,
                 alpha = 0.6) +
  geom_point(aes(colour = sig), size = 0.75) +
  facet_grid(trans_lab ~ dv_lab, scales = "free") +
  scale_colour_manual(values = COL_SIG, name = NULL, drop = FALSE) +
  labs(title = "Specification curve, contemporaneous fire",
       subtitle = "One point per specification; 95% CI; dashed line is the panel median",
       x = "Specification, ordered by estimate within panel", y = "IV estimate") +
  theme_ff(13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(g4, "deck_fig4_speccurve", width = 7.8, height = 5.6)

say("done")
