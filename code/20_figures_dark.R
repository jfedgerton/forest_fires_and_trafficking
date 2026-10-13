################################################################################
# 20_figures_dark.R
#
# Dark-background versions of the figures that go on slides, matching the
# black Edgerton deck template. A white-panelled ggplot dropped onto a black
# slide reads as a hole punched in the slide; these are drawn on the slide's own
# background so the figure sits in the page rather than on top of it.
#
# Same estimates as 11_figures.R and 19_figures_deck.R - only the ink changes.
# If a point moves between the light and dark versions, something is wrong.
#
# Canvas sizes are deliberately small. A figure drawn at 11 x 7 in and placed on
# a slide at 5.5 in scales to 0.5, which puts a 10 pt axis label on the
# projector at 5 pt. Each figure here is drawn at roughly the size it will
# occupy, so the type survives.
#
# Output: output/figures/dk_*.{pdf,png}
################################################################################

source("00_config.R")

BG   <- "#000000"
INK  <- "#FFFFFF"
GRID <- "#3A3A3A"
PINK <- "#FD738A"     # the template's emphasis colour
AMBR <- "#FFC000"
BLUE <- "#7FB3E8"
COL_DV_DK  <- c(TF = PINK, UE = BLUE)
COL_SIG_DK <- c("p < 0.05" = PINK, "0.05 <= p < 0.10" = AMBR, "p >= 0.10" = "#8A8A8A")

# Axis text horizontal throughout; strips on top and right, never on the left.
theme_dk <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = BG, colour = NA),
      panel.background  = element_rect(fill = BG, colour = NA),
      legend.background = element_rect(fill = BG, colour = NA),
      legend.key        = element_rect(fill = BG, colour = NA),
      panel.grid.major  = element_line(colour = GRID, linewidth = 0.25),
      panel.grid.minor  = element_blank(),
      text              = element_text(colour = INK),
      axis.text         = element_text(colour = "#C9C9C9"),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      axis.text.y       = element_text(angle = 0, hjust = 1),
      axis.title        = element_text(colour = INK),
      strip.background  = element_rect(fill = "#1E1E1E", colour = NA),
      strip.text        = element_text(colour = INK, angle = 0, lineheight = 1.05),
      strip.text.y      = element_text(colour = INK, angle = 0),
      plot.title        = element_text(colour = INK, face = "bold"),
      plot.subtitle     = element_text(colour = "#B8B8B8"),
      legend.position   = "bottom",
      legend.text       = element_text(colour = INK))
}

save_dk <- function(g, stub, width, height) {
  ggsave(file.path(DIR_FIG, paste0(stub, ".pdf")), g, width = width, height = height,
         bg = BG)
  ggsave(file.path(DIR_FIG, paste0(stub, ".png")), g, width = width, height = height,
         dpi = 200, bg = BG)
  say("  wrote ", stub)
}

d    <- load_panel()
tidy <- fread(file.path(DIR_RES, "main_models_tidy.csv"))
LAB_TRAN_D <- c(binary = "Binary", count = "Count", log1p = "Log(1+y)")
LAB_DV_D   <- c(TF = "TF - rescues", UE = "UE - victims")

## =============================================================================
## dk 1 - first-stage binscatter
## =============================================================================
say("dark figure 1: first stage")
Z <- VPD_MAIN
fs <- d[stats::complete.cases(d[, c(FIRE, Z), with = FALSE])]
m_y <- feols(as.formula(sprintf("%s ~ 1 | %s", FIRE, FE_MAIN)), data = fs, notes = FALSE)
m_x <- feols(as.formula(sprintf("%s ~ 1 | %s", Z, FE_MAIN)), data = fs, notes = FALSE)
bs <- data.table(ry = resid(m_y), rx = resid(m_x))
bs[, bin := cut(rx, breaks = quantile(rx, seq(0, 1, length.out = 41)),
                include.lowest = TRUE, labels = FALSE)]
pts <- bs[, .(x = mean(rx), y = mean(ry)), by = bin]
fit <- coef(lm(ry ~ rx, data = bs))

g1 <- ggplot(pts, aes(x, y)) +
  geom_abline(intercept = fit[1], slope = fit[2], colour = PINK, linewidth = 0.8) +
  geom_point(colour = INK, size = 2) +
  labs(title = "First stage: dryness predicts burning",
       subtitle = paste0("40 equal-sized bins. Both axes residualised on the\n",
                         "municipality and year-quarter fixed effects."),
       x = "Vapour pressure deficit, residualised (kPa)",
       y = "Share of municipality burned, residualised (pp)") +
  theme_dk(14)
save_dk(g1, "dk_fig1_first_stage", 7.6, 4.6)

## =============================================================================
## dk 2 - one-SD effects
## =============================================================================
say("dark figure 2: one-SD effects")
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
eff[, ctrl_lab := factor(c(none_fullwindow = "No ctrl\n02-24",
                           none            = "No ctrl\n02-21",
                           full            = "Full ctrl\n02-21")[controls],
                         levels = c("No ctrl\n02-24", "No ctrl\n02-21", "Full ctrl\n02-21"))]
eff[, panel := factor(paste0(dv_lab, "  |  ", trans_lab))]

g2 <- ggplot(eff, aes(ctrl_lab, pct, colour = dv_lab)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.4) +
  geom_pointrange(aes(ymin = pct_lo, ymax = pct_hi), linewidth = 0.8, size = 0.55) +
  facet_wrap(~ panel, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = unname(COL_DV_DK), guide = "none") +
  labs(title = "One SD more burned area",
       subtitle = sprintf(paste0("1 SD = %.1f pp of municipal area. Binary and count as a percent of their own\n",
                                 "baseline; log(1+y) as a proportional change. Free y scale by panel; 95%% CI."),
                          sd_fire),
       x = NULL, y = "Percent change") +
  scale_y_continuous(expand = expansion(mult = 0.10)) +
  theme_dk(13) + theme(axis.text.x = element_text(size = 9.5))
save_dk(g2, "dk_fig2_one_sd", 7.6, 5.0)

## =============================================================================
## dk 3 - specification curve, contemporaneous
## =============================================================================
say("dark figure 3: specification curve")
GR <- file.path(DIR_RES, "spec_grid_2sls.csv")
stopifnot(file.exists(GR))
s <- fread(GR)[status == "ok" & lag_only == 0 & lagid == "L0"]
s[, trans_lab := factor(LAB_TRAN_D[transformation], levels = unname(LAB_TRAN_D))]
s <- s[!is.na(trans_lab)]        # drops the intensity form, which has its own panel
s[, dv_lab := factor(LAB_DV_D[source], levels = unname(LAB_DV_D))]
s[, sig := factor(fifelse(iv_p < 0.05, "p < 0.05",
                   fifelse(iv_p < 0.10, "0.05 <= p < 0.10", "p >= 0.10")),
                  levels = names(COL_SIG_DK))]
s[, `:=`(lo = iv_b - 1.96 * iv_se, hi = iv_b + 1.96 * iv_se)]
setorder(s, iv_b)
s[, x := seq_len(.N), by = .(trans_lab, dv_lab)]
med <- s[, .(m = median(iv_b, na.rm = TRUE)), by = .(trans_lab, dv_lab)]
say("  specifications: ", nrow(s))

g3 <- ggplot(s, aes(x, iv_b)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.4) +
  geom_hline(data = med, aes(yintercept = m), linetype = 2, colour = BLUE,
             linewidth = 0.45) +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = sig), linewidth = 0.3, alpha = 0.7) +
  geom_point(aes(colour = sig), size = 0.75) +
  facet_grid(trans_lab ~ dv_lab, scales = "free") +
  scale_colour_manual(values = COL_SIG_DK, name = NULL, drop = FALSE) +
  labs(title = "Specification curve, contemporaneous fire",
       subtitle = "One point per specification; 95% CI; dashed line is the panel median",
       x = "Specification, ordered by estimate within panel", y = "IV estimate") +
  # four breaks, not the default: at slide size the bottom tick of one panel and
  # the top tick of the next run together
  scale_y_continuous(n.breaks = 4, expand = expansion(mult = 0.12)) +
  theme_dk(13) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 10))
save_dk(g3, "dk_fig3_speccurve", 6.8, 5.0)

## =============================================================================
## dk 4 - OLS by treatment transform (from the 16 sweep)
## =============================================================================
say("dark figure 4: OLS by transform")
r <- fread(file.path(DIR_RES, "robustness_sweep.csv"))
r <- r[is.finite(b) & is.finite(se) & se > 0 & se < 1e3]
LAB_FIRE_M <- c(pct = "% of area burned", nodbl = "km2, repeat once",
                total = "km2, repeat counted")
# Short labels: at the width this figure gets on a slide, "Raw count" and
# "Log(1 + y)" run into each other.
LAB_OUT    <- c(binary = "Binary", count = "Count", log1p = "Log(1+y)")
o <- r[sweep == "1_OLS"]
o[, lo := b - 1.96 * se][, hi := b + 1.96 * se]
o[, fm_lab  := factor(LAB_FIRE_M[fire_measure], levels = LAB_FIRE_M)]
o[, out_lab := factor(LAB_OUT[out_trans], levels = LAB_OUT)]
o[, dv_lab  := factor(LAB_DV_D[dv], levels = unname(LAB_DV_D))]
o[, ft_lab  := factor(c(raw = "Raw fire", log = "Log(1 + fire)")[fire_trans],
                      levels = c("Raw fire", "Log(1 + fire)"))]

g4 <- ggplot(o, aes(x = out_lab, y = b, colour = fm_lab, shape = controls)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.4) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.72), size = 0.4, linewidth = 0.55) +
  facet_grid(ft_lab ~ dv_lab, scales = "free_y") +
  scale_colour_manual(values = c(PINK, AMBR, BLUE), name = NULL) +
  scale_shape_manual(values = c(full = 16, none = 1), name = NULL,
                     labels = c(full = "Full controls", none = "No controls")) +
  labs(title = "OLS: null on the raw scale,\npositive on the log scale",
       subtitle = "Burned area is 73% exact zeros, skew 5.8.\nFree y scale by row.",
       x = NULL, y = "Coefficient on fire") +
  theme_dk(13) + theme(legend.box = "vertical", legend.spacing.y = unit(1, "pt"))
save_dk(g4, "dk_fig4_ols_transform", 6.2, 4.6)

## =============================================================================
## dk 5 - seasonality, for the "is this just seasonality" question
## =============================================================================
say("dark figure 5: seasonality")
sea <- d[, .(Fire = mean(get(FIRE), na.rm = TRUE),
             TF   = mean(tf_count > 0, na.rm = TRUE),
             UE   = mean(ue_count > 0, na.rm = TRUE)), by = quarter]
sea <- melt(sea, id.vars = "quarter")
sea[, variable := factor(variable, levels = c("Fire", "TF", "UE"),
      labels = c("Share of area burned", "P(any rescue)", "P(any victim)"))]
g5 <- ggplot(sea, aes(factor(quarter), value, group = 1)) +
  geom_line(colour = PINK, linewidth = 0.9) +
  geom_point(colour = INK, size = 2) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Seasonality is real, and the fixed effects absorb it",
       subtitle = "Quarter-of-year means. The year-quarter fixed effect removes any national seasonal pattern.",
       x = "Quarter of year", y = NULL) +
  theme_dk(13)
save_dk(g5, "dk_fig5_seasonality", 8.0, 3.2)

say("done")
