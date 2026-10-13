################################################################################
# 16_seasonal_windows.R
#
# Two jobs, two windows: does the agricultural channel live in a different
# quarter from the fire channel?
#
# THE ARGUMENT. The instrument should sit in the window that moves fire and does
# NOT move agriculture. The Brazilian calendar suggests that window exists:
#
#     Q3 (Jul-Sep)     fire season peak; main crop already harvested, next crop
#                      not yet planted
#     Q4 (Oct-Dec)     main crop planted
#     Q1 (Jan-Mar)     main crop grows and is harvested
#     Q2 (Apr-Jun)     safrinha (second maize) grows
#
# So for harvest year y the crop is exposed in Q4(y-1) and Q1(y), while the fire
# season Q3(y) falls AFTER that harvest and BEFORE the next planting.
#
# NO AVERAGING. An annual mean VPD has no agronomic meaning - it mixes the
# window that stresses the crop with windows where nothing is in the ground, and
# the separability run showed Q3 and Q4+Q1 push burned area in OPPOSITE
# directions (+4.85 vs -8.40), so averaging can cancel real effects. Even a
# growing-season average imposes that Q4 and Q1 matter equally. Every quarter
# therefore enters as its own regressor and the data says which window carries
# each channel. This is how the climate-agriculture literature specifies it.
#
# WHERE THE AGGREGATION PROBLEM ACTUALLY LIVES. Only here. The main model is
# quarterly on both sides - quarterly VPD, quarterly fire, quarterly rescues -
# and needs no aggregation at all. The agricultural controls are the only annual
# object in the panel, and testing them against a quarterly instrument is what
# produced the uninformative balance nulls in the first place.
#
# FIVE WINDOWS, covering one full agricultural cycle for harvest year y:
#     vpd_Q3_lag   Q3(y-1)   fire season before planting
#     vpd_Q4_lag   Q4(y-1)   planting
#     vpd_Q1       Q1(y)     growing and harvest
#     vpd_Q2       Q2(y)     safrinha
#     vpd_Q3       Q3(y)     fire season after harvest  <- THE INSTRUMENT WINDOW
#
# WHAT WOULD VINDICATE THE DESIGN
#   (a) agricultural outcomes load on Q4_lag and Q1 and NOT on Q3, and
#   (b) burned area loads on Q3.
# Then the instrument window is clean and the agricultural channel lives
# somewhere the instrument does not.
#
# WHAT WOULD NOT. If Q3 predicts the agricultural outcomes too, the windows do
# not separate and this route is closed. That is a real possible outcome and the
# script is written to report it plainly rather than to search for the other one.
#
# RESULT WHEN RUN, 29 September 2026. The mechanical verdict printed NOT CLEAN -
# Q3 loaded on 5 of 6 agricultural outcomes. But the verdict counts significance
# and ignores SIGN, and the sign structure says something different: crop stress
# is NEGATIVE in Q1 and Q2 exactly as agronomy requires, while every Q3
# coefficient is POSITIVE, which is the wrong sign for crop stress and sits
# after the harvest in time. Read the signs, not the verdict.
#
# Output: output/results/seasonal_windows.csv
#         output/figures/fig_seasonal_windows.{pdf,png}
# Seed 123. Run from code/.
################################################################################

source("00_config.R")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

Z <- VPD_MAIN          # vpd_max_q, the instrument as the paper uses it
AG_OUT <- c("log_crop_area_planted", "log_crop_area_harvested",
            "log_crop_value", "log_cattle", "crop_area_share", "cattle_density")

d <- load_panel()[year <= TRUNC_YEAR]

## ---- build the five quarter-specific VPD columns on a municipality-year key --
# Each is pulled from the quarterly panel and keyed to the HARVEST year it
# belongs to, so a lag is a relabelled year rather than a shift operation - the
# construction cannot silently misalign.
grab <- function(qtr, shift, nm) {
  x <- d[quarter == qtr, .(code_muni_7, year = year + shift, v = get(Z))]
  setnames(x, "v", nm)
  x
}
w <- Reduce(function(a, b) merge(a, b, by = c("code_muni_7", "year")),
            list(grab(3L, 1L, "vpd_Q3_lag"),
                 grab(4L, 1L, "vpd_Q4_lag"),
                 grab(1L, 0L, "vpd_Q1"),
                 grab(2L, 0L, "vpd_Q2"),
                 grab(3L, 0L, "vpd_Q3")))
WIN <- c("vpd_Q3_lag", "vpd_Q4_lag", "vpd_Q1", "vpd_Q2", "vpd_Q3")
GROW <- c("vpd_Q4_lag", "vpd_Q1")          # the crop's exposure window
FIRE <- "vpd_Q3"                            # the instrument's window

## ---- annual outcomes ---------------------------------------------------------
ag <- d[, c(lapply(.SD, function(x) x[1]),
            .(burn_yr = sum(percent_area_burned, na.rm = TRUE),
              tf_yr   = sum(tf_count, na.rm = TRUE),
              ue_yr   = sum(ue_count, na.rm = TRUE),
              uf_code = uf_code[1])),
        by = .(code_muni_7, year), .SDcols = AG_OUT]

y <- merge(ag, w, by = c("code_muni_7", "year"))
y <- y[stats::complete.cases(y[, ..WIN])]
y[, muni := factor(code_muni_7)][, yr := factor(year)][, uf := factor(uf_code)]
y[, `:=`(tf_bin = as.numeric(tf_yr > 0), ue_bin = as.numeric(ue_yr > 0))]

say("municipality-year frame: ", nrow(y), " rows, ", uniqueN(y$code_muni_7),
    " municipalities, ", min(y$year), "-", max(y$year))
stopifnot(nrow(y) > 90000)

## ---- how distinct are the windows? ------------------------------------------
say("")
say("Within municipality + year, correlation between the five windows:")
R <- copy(y)
for (v in WIN) R[, (v) := resid(feols(as.formula(paste(v, "~ 1 | muni + yr")),
                                      R, notes = FALSE))]
print(round(cor(R[, ..WIN]), 2))

## ---- the regressions ---------------------------------------------------------
# All five windows enter together, so each coefficient is conditional on the
# others. Clustered on municipality, with state clustering reported alongside
# because the two disagreed badly in the annual balance tests and the
# disagreement should stay visible.
fit <- function(dv, dat = y) {
  f <- as.formula(sprintf("%s ~ %s | muni + yr", dv, paste(WIN, collapse = " + ")))
  m  <- feols(f, dat, cluster = ~muni, notes = FALSE)
  m2 <- feols(f, dat, cluster = ~uf,   notes = FALSE)
  ct <- coeftable(m); ct2 <- coeftable(m2)
  data.table(outcome = dv, window = WIN,
             b = ct[WIN, 1], se = ct[WIN, 2], p = ct[WIN, 4],
             p_state = ct2[WIN, 4],
             p_joint_grow = one(safely(fixest::wald(
               m, keep = paste0("^(", paste(GROW, collapse = "|"), ")$"),
               print = FALSE)$p)),
             p_fire = ct[FIRE, 4],
             n = nobs(m))
}

say("")
say("=== BLOCK 1: agricultural outcomes ===")
A <- rbindlist(lapply(AG_OUT, fit))
print(dcast(A, outcome ~ window, value.var = "b")[, lapply(.SD, function(x)
  if (is.numeric(x)) signif(x, 3) else x)])
say("")
say("p-values (municipality clustered):")
print(dcast(A, outcome ~ window, value.var = "p")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 4) else x)])

say("")
say("=== BLOCK 2: burned area, the first stage ===")
B <- fit("burn_yr")
print(B[, .(window, b = signif(b, 3), se = signif(se, 3),
            p = round(p, 4), p_state = round(p_state, 4))])

say("")
say("=== BLOCK 3: rescues, annual (PREVIEW ONLY - not the main specification, ")
say("    which is quarterly on both sides) ===")
C <- rbindlist(lapply(c("tf_bin", "ue_bin"), fit))
print(dcast(C, outcome ~ window, value.var = "b")[, lapply(.SD, function(x)
  if (is.numeric(x)) signif(x, 3) else x)])
print(dcast(C, outcome ~ window, value.var = "p")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 4) else x)])

out <- rbindlist(list(cbind(block = "agriculture", A),
                      cbind(block = "fire", B),
                      cbind(block = "rescues_annual", C)), fill = TRUE)
fwrite(out, file.path(DIR_RES, "seasonal_windows.csv"))

## ---- the verdict -------------------------------------------------------------
# NOTE, added after the 29 September run: this verdict counts significance and
# ignores sign, which made it misleading. Keep it for the mechanical summary but
# read the sign structure in BLOCK 1 before drawing any conclusion.
say("")
say("=== DOES THE INSTRUMENT WINDOW SEPARATE? (counts only - read signs too) ===")
ag_fire_sig  <- A[window == FIRE & p < 0.05, .N]
ag_grow_sig  <- A[window %in% GROW & p < 0.05, .N]
burn_fire_p  <- B[window == FIRE, p]
say("  agricultural outcomes loading on the FIRE window (Q3): ",
    ag_fire_sig, " of ", length(AG_OUT))
say("  agricultural outcomes loading on a CROP window (Q4lag, Q1): ",
    ag_grow_sig, " of ", 2 * length(AG_OUT))
say("  burned area on the fire window: p = ", signif(burn_fire_p, 3))
say("  SIGN CHECK - mean Q3 coefficient on crop outcomes: ",
    signif(A[window == FIRE & grepl("crop", outcome), mean(b)], 3),
    "  (positive = wrong sign for crop stress, and Q3 is post-harvest)")

## ---- figure ------------------------------------------------------------------
pl <- rbindlist(list(cbind(A, blk = "Agricultural outcomes"),
                     cbind(B, blk = "Burned area")), fill = TRUE)
pl[, win := factor(window, levels = WIN,
                   labels = c("Q3 (y-1)\nfire season", "Q4 (y-1)\nplanting",
                              "Q1 (y)\ngrow + harvest", "Q2 (y)\nsafrinha",
                              "Q3 (y)\nfire season"))]
pl[, sig := p < 0.05]
g <- ggplot(pl, aes(x = win, y = b)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se),
                width = 0, linewidth = 0.35, colour = "grey40") +
  geom_point(aes(shape = sig), size = 1.8) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("p >= .05", "p < .05"), name = NULL) +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(x = NULL, y = "Coefficient on quarter-specific VPD, with 95% CI",
       title = "Which quarter's drying moves each outcome?",
       subtitle = paste0("All five windows entered together. Municipality and ",
                         "year fixed effects, clustered on municipality, ",
                         min(y$year), "-", max(y$year), ".")) +
  theme_ff(9)
save_fig(g, "fig_seasonal_windows", width = 11, height = 7)

say("")
say("done")
