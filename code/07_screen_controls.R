################################################################################
# 07_screen_controls.R
#
# Screens candidate controls before they are used, and exports the coefficient
# plots. Defines the agricultural control sets that 08 and 09 consume, so that
# no downstream script hard-codes the screen's result.
#
# TWO SCREENS, because "which controls matter" has two readings that answer
# different questions.
#
#   A - BALANCE / EXCLUSION.   control ~ instrument | muni + yq
#       Does the instrument move the candidate control? If VPD predicts crop
#       area or cattle stocking, that variable is POST-INSTRUMENT, and
#       conditioning on it blocks part of the pathway the IV runs through - the
#       textbook bad-control problem. A significant coefficient here is a reason
#       for caution about that variable, not a reason to include it.
#
#   B - OUTCOME RELEVANCE.     outcome ~ instrument + control | muni + yq
#       Conditional on the instrument, does the control predict the outcome?
#       This is the screen that defines the ag@.05 and ag@.10 sets.
#
# A CAVEAT THAT BELONGS IN THE APPENDIX TEXT. Selecting controls on their own
# significance and then reporting the resulting models is pre-test selection:
# the nominal p-values in the selected-control specifications are conditional on
# the screen having been run and are not exact. That is defensible inside a
# specification curve, where the selected sets are two branches among many. It
# is not defensible as a route to one preferred model - which is why 08 uses a
# fixed full control set rather than a screened one.
#
# Output: output/results/ag_control_selection.rds  (consumed by 00_config.R)
#         output/results/ag_control_screen.csv
#         output/figures/fig0_screen_{balance,relevance}.{pdf,png}
################################################################################

source("00_config.R")

d <- load_panel()
say("panel: ", nrow(d), " rows")

SCREEN_OUTCOMES <- DVS
rows <- list(); k <- 0L

## ---- Screen A: does the instrument move the control? -------------------------
say("screen A: balance / exclusion")
for (v in names(CANDIDATE_CONTROLS)) {
  for (zn in VPD_ALL) {
    for (L in 0:2) {
      z  <- paste0(zn, sfx(L))
      dd <- d[stats::complete.cases(d[, c(v, z), with = FALSE])]
      if (nrow(dd) < 1000L) next
      m <- safely(feols(as.formula(sprintf("%s ~ %s | %s", v, z, FE_MAIN)),
                        data = dd, cluster = CLUSTER, notes = FALSE))
      s <- pick(m, z)
      if (all(is.na(s))) next
      k <- k + 1L
      rows[[k]] <- data.table(
        screen = "A_balance", variable = v, var_label = CANDIDATE_CONTROLS[[v]],
        is_ag = v %in% AG_CANDIDATES, instrument = zn, lag = L,
        outcome = NA_character_, estimate = s[1], se = s[2], p = s[4],
        n = one(nobs(m)),
        est_sd = s[1] * sd(dd[[z]], na.rm = TRUE) / sd(dd[[v]], na.rm = TRUE))
    }
  }
}

## ---- Screen B: does the control predict the outcome? -------------------------
say("screen B: outcome relevance")
for (v in names(CANDIDATE_CONTROLS)) {
  for (zn in VPD_ALL) {
    for (L in 0:2) {
      for (yn in names(SCREEN_OUTCOMES)) {
        z  <- paste0(zn, sfx(L)); yv <- SCREEN_OUTCOMES[[yn]]
        dd <- d[stats::complete.cases(d[, c(v, z, yv), with = FALSE])]
        if (nrow(dd) < 1000L) next
        m <- safely(fepois(as.formula(sprintf("%s ~ %s + %s | %s", yv, z, v, FE_MAIN)),
                           data = dd, cluster = CLUSTER, notes = FALSE))
        s <- pick(m, v)
        if (all(is.na(s))) next
        k <- k + 1L
        rows[[k]] <- data.table(
          screen = "B_outcome", variable = v, var_label = CANDIDATE_CONTROLS[[v]],
          is_ag = v %in% AG_CANDIDATES, instrument = zn, lag = L,
          outcome = yv, estimate = s[1], se = s[2], p = s[4], n = one(nobs(m)),
          est_sd = s[1] * sd(dd[[v]], na.rm = TRUE))
      }
    }
  }
}

screen <- rbindlist(rows, fill = TRUE)
screen[, `:=`(sig_05 = as.integer(p < 0.05), sig_10 = as.integer(p < 0.10),
              ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)]
fwrite(screen, file.path(DIR_RES, "ag_control_screen.csv"))
say("screen rows: ", nrow(screen))

## ---- selected sets -----------------------------------------------------------
# A variable qualifies if it is significant in a majority of its screen-B cells
# for EITHER outcome. A both-outcomes rule collapses the .05 and .10 sets onto
# the same two cattle variables and makes the threshold a dead dimension in the
# specification curve; the either-outcome rule keeps the branches distinct.
selB <- screen[screen == "B_outcome" & is_ag == TRUE,
               .(cells = .N, frac_sig_05 = mean(sig_05, na.rm = TRUE),
                 frac_sig_10 = mean(sig_10, na.rm = TRUE)),
               by = .(variable, outcome)]
sel_wide <- dcast(selB, variable ~ outcome,
                  value.var = c("frac_sig_05", "frac_sig_10"))
balA <- screen[screen == "A_balance" & is_ag == TRUE,
               .(frac_moved_05 = mean(sig_05, na.rm = TRUE)), by = variable]

saveRDS(list(selection_table = sel_wide, balance_table = balA,
             endogenous_flag = balA[frac_moved_05 > 0.5, variable]),
        SELECTION)

sets <- build_control_sets()
cat("\n===== SELECTED CONTROL SETS =====\n")
for (nm in names(sets))
  cat(sprintf("  %-12s %s\n", nm,
              if (length(sets[[nm]])) paste(sets[[nm]], collapse = ", ") else "(none)"))
cat("\nAgricultural variables the instrument moves in >50% of balance cells: ",
    if (nrow(balA[frac_moved_05 > 0.5])) paste(balA[frac_moved_05 > 0.5, variable],
                                               collapse = ", ") else "(none)", "\n")
cat("\n--- screen A balance ---\n"); print(balA[order(-frac_moved_05)])
cat("\n--- screen B selection ---\n"); print(sel_wide)

nA <- screen[screen == "A_balance", .N]
say(sprintf("screen A: %d of %d tests significant at 5%% (chance expectation %.1f)",
            screen[screen == "A_balance", sum(sig_05, na.rm = TRUE)], nA, 0.05 * nA))

## ---- figures -----------------------------------------------------------------
pA <- screen[screen == "A_balance"]
pA[, var_label := factor(var_label, levels = rev(unname(CANDIDATE_CONTROLS)))]
pA[, lag_lab := factor(paste("Lag", lag), levels = paste("Lag", 0:2))]
pA[, inst_lab := factor(LAB_Z[match(instrument, VPD_ALL)],
                        levels = LAB_Z[1:3])]

gA <- ggplot(pA, aes(est_sd, var_label, colour = is_ag)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 1.5) +
  facet_grid(lag_lab ~ inst_lab, scales = "free_x") +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                      labels = c(`TRUE` = "Agricultural", `FALSE` = "Economic"),
                      name = NULL) +
  labs(title = "Screen A: does the instrument move the candidate control?",
       subtitle = "Standardised coefficient; municipality and year-quarter FE; two-way clustered SEs",
       x = "SD of the control per SD of the instrument", y = NULL) +
  theme_ff()
save_fig(gA, "fig0_screen_balance", 11, 8)

pB <- screen[screen == "B_outcome"]
pB[, var_label := factor(var_label, levels = rev(unname(CANDIDATE_CONTROLS)))]
pB[, lag_lab := factor(paste("Lag", lag), levels = paste("Lag", 0:2))]
pB[, dv_lab := factor(fifelse(outcome == "tf_count", "TF\nworkers rescued",
                                                     "UE\nvictims"))]
gB <- ggplot(pB, aes(estimate, var_label, colour = instrument)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  position = position_dodge(width = 0.7), size = 0.25) +
  facet_grid(lag_lab ~ dv_lab, scales = "free_x") +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(title = "Screen B: does the control predict the outcome, given the instrument?",
       subtitle = "FE-Poisson coefficient with 95% CI",
       x = "Coefficient", y = NULL) +
  theme_ff()
save_fig(gB, "fig0_screen_relevance", 11, 9)

say("done")
