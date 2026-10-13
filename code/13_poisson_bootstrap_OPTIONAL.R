################################################################################
# 41_count_models_bootstrap.R
#
# Count-outcome IV models with valid standard errors.
#
# WHY THIS EXISTS. The count outcomes are estimated by FE-Poisson, and their IV
# is a control function: the first-stage residual is entered as a regressor in
# the Poisson. The analytic standard error that fixest reports for that model
# treats the residual as data rather than as an estimate, so it is wrong - which
# is why every count row in scripts 38-40 carries iv_se_valid = FALSE and is
# dropped from the curves. This script replaces those standard errors with a
# cluster bootstrap over municipalities, which is the unit the models cluster on
# and the unit the panel's dependence runs through.
#
# ESTIMANDS (contemporaneous specification only)
#   count      FE-Poisson on the count, full sample                  [bootstrapped]
#   log1p      OLS on log(1 + y), full sample                        [analytic SE, valid]
# The log1p column is carried alongside so the two are readable together; it is
# the same estimate already in the main tables, not a re-estimation.
#
# BOOTSTRAP. Municipalities are resampled with replacement, B = 199. Each drawn
# municipality gets a fresh fixed-effect id, so a municipality drawn twice
# contributes two separate fixed effects rather than being silently pooled.
# Both stages - the first-stage OLS and the Poisson - are re-estimated inside
# every replication, which is the point: the whole two-step procedure has to be
# bootstrapped, not just the second stage.
# Reported: bootstrap SE (sd of the replication estimates), a percentile 95%
# interval, and a normal-approximation p-value. Where the two disagree, trust
# the percentile interval - the Poisson CF estimates are not symmetric.
#
# Bootstrap draws are shared across instrument sets within a (DV, control set)
# context, so instruments are compared on the same resamples rather than on
# independent noise.
#
# Checkpointed one .rds per (DV, control set). Re-run to resume.
# Seed 123.
################################################################################

set.seed(123)
suppressPackageStartupMessages({library(data.table); library(fixest)})
setFixest_nthreads(max(1, parallel::detectCores()))
options(width = 200)

B <- as.integer(Sys.getenv("NBOOT", "199"))
OUT <- "count_bootstrap"; dir.create(OUT, showWarnings = FALSE)
LOGF <- file("41_count_bootstrap.log", open = "at")
say <- function(...) { m <- sprintf("%s %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); writeLines(m, LOGF); flush(LOGF) }

d <- fread("final_panel_muni_quarter.csv.gz")
d[, yq := factor(year_quarter)]
say("panel: ", nrow(d), " rows | B = ", B)

sel <- readRDS("ag_control_selection.rds")
st  <- as.data.table(sel$selection_table)
ag05 <- st[frac_sig_05_tf_count > 0.5 | frac_sig_05_ue_count > 0.5, variable]
ag10 <- st[frac_sig_10_tf_count > 0.5 | frac_sig_10_ue_count > 0.5, variable]
GDP <- c("log_pop", "log_gdp_pc")
CONTROLS <- list(K0_none = character(0), K1_ag05 = ag05, K2_ag10 = ag10,
                 K3_ag05_gdp = c(ag05, GDP), K4_ag10_gdp = c(ag10, GDP),
                 K5_gdp_only = GDP)
NEEDS_TRUNC <- sapply(CONTROLS, function(v) any(v %in% c("log_gdp_pc", GDP[2])))

VPD <- c("vpd_mean_q", "vpd_max_q", "vpd_excess3_q")
INSTSETS <- list(Z1_mean = VPD[1], Z2_max = VPD[2], Z3_exc3 = VPD[3],
                 Z4_mean_max = VPD[1:2], Z5_mean_exc3 = VPD[c(1, 3)],
                 Z6_max_exc3 = VPD[2:3], Z7_all = VPD)
FIRE <- "percent_area_burned"
DVS <- c(TF = "tf_count", UE = "ue_count")
rhs <- function(v) if (length(v)) paste(v, collapse = " + ") else "1"

# one control-function fit; returns the coefficient on fire, or NA
cf_fit <- function(DD, yv, zs, ctrl, fecol) {
  f1 <- as.formula(sprintf("%s ~ %s + %s | %s + yq", FIRE, rhs(zs), rhs(ctrl), fecol))
  m1 <- tryCatch(feols(f1, data = DD, notes = FALSE), error = function(e) NULL)
  if (is.null(m1)) return(NA_real_)
  DD[, .cf := NA_real_]
  DD[fixest::obs(m1), .cf := residuals(m1)]
  f2 <- as.formula(sprintf("%s ~ %s + .cf + %s | %s + yq", yv, FIRE, rhs(ctrl), fecol))
  m2 <- tryCatch(fepois(f2, data = DD, notes = FALSE), error = function(e) NULL)
  if (is.null(m2)) return(NA_real_)
  ct <- tryCatch(coeftable(m2), error = function(e) NULL)
  if (is.null(ct) || !(FIRE %in% rownames(ct))) return(NA_real_)
  as.numeric(ct[FIRE, 1])
}

for (dvn in names(DVS)) {
  yv <- DVS[[dvn]]
  for (kn in names(CONTROLS)) {
    fp <- file.path(OUT, sprintf("%s__%s.rds", dvn, kn))
    if (file.exists(fp)) next
    ctrl <- CONTROLS[[kn]]
    base <- if (NEEDS_TRUNC[[kn]]) d[year <= 2021L] else d
    need <- c(yv, FIRE, VPD, ctrl)
    base <- base[stats::complete.cases(base[, ..need])]
    base[, muni := factor(code_muni_7)]
    munis <- unique(base$code_muni_7)
    setkey(base, code_muni_7)
    say(sprintf("%s / %s : n = %d rows, %d munis", dvn, kn, nrow(base), length(munis)))

    ## ---- point estimates, all instrument sets --------------------------------
    pt <- sapply(names(INSTSETS), function(zn)
      cf_fit(copy(base), yv, INSTSETS[[zn]], ctrl, "muni"))

    ## ---- log1p companion (analytic SE is valid here) -------------------------
    base[, .ylog := log1p(get(yv))]
    lg <- rbindlist(lapply(names(INSTSETS), function(zn) {
      zs <- INSTSETS[[zn]]
      m <- tryCatch(feols(as.formula(sprintf(".ylog ~ %s | muni + yq | %s ~ %s",
                                             rhs(ctrl), FIRE, rhs(zs))),
                          data = base, cluster = ~muni + yq, notes = FALSE),
                    error = function(e) NULL)
      tm <- paste0("fit_", FIRE)
      if (is.null(m) || !(tm %in% rownames(coeftable(m))))
        return(data.table(zset = zn, log1p_b = NA_real_, log1p_se = NA_real_, log1p_p = NA_real_))
      ct <- coeftable(m)
      data.table(zset = zn, log1p_b = ct[tm, 1], log1p_se = ct[tm, 2], log1p_p = ct[tm, 4])
    }))

    ## ---- cluster bootstrap ----------------------------------------------------
    t0 <- Sys.time()
    boot <- matrix(NA_real_, nrow = B, ncol = length(INSTSETS),
                   dimnames = list(NULL, names(INSTSETS)))
    for (b in seq_len(B)) {
      drawn <- sample(munis, length(munis), replace = TRUE)
      idx <- data.table(code_muni_7 = drawn, .bid = seq_along(drawn))
      BD <- base[idx, on = "code_muni_7", allow.cartesian = TRUE]
      BD[, bmuni := factor(.bid)]          # a muni drawn twice gets two FEs
      for (zn in names(INSTSETS))
        boot[b, zn] <- cf_fit(BD, yv, INSTSETS[[zn]], ctrl, "bmuni")
      rm(BD)
      if (b %% 25 == 0) {
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        say(sprintf("  %s/%s boot %d/%d (%.1f min, eta %.0f min)",
                    dvn, kn, b, B, el, el / b * (B - b)))
      }
    }

    res <- rbindlist(lapply(names(INSTSETS), function(zn) {
      bs <- boot[, zn]; bs <- bs[is.finite(bs)]
      se <- if (length(bs) > 10) sd(bs) else NA_real_
      z  <- if (is.finite(se) && se > 0) pt[[zn]] / se else NA_real_
      data.table(source = dvn, kset = kn, zset = zn,
                 window = if (NEEDS_TRUNC[[kn]]) "2002-2021" else "2002-2024",
                 n_obs = nrow(base), n_boot_ok = length(bs),
                 count_b = pt[[zn]], count_se_boot = se,
                 count_ci_lo = if (length(bs) > 10) quantile(bs, .025) else NA_real_,
                 count_ci_hi = if (length(bs) > 10) quantile(bs, .975) else NA_real_,
                 count_t = z, count_p = if (is.finite(z)) 2 * pnorm(-abs(z)) else NA_real_)
    }))
    res <- merge(res, lg, by = "zset", sort = FALSE)
    saveRDS(res, fp)
    say(sprintf("  saved %s (%.1f min)", fp,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

fs <- list.files(OUT, pattern = "\\.rds$", full.names = TRUE)
all <- rbindlist(lapply(fs, readRDS), fill = TRUE)
fwrite(all, "count_models_bootstrapped.csv")
say("contexts complete: ", length(fs), " / 12 | rows ", nrow(all))
if (length(fs) == 12L) {
  cat("\n===== COUNT (FE-Poisson, bootstrapped SE) vs LOG1P, contemporaneous =====\n")
  print(all[, .(n = .N,
                count_median = round(median(count_b, na.rm = TRUE), 5),
                count_pct_pos = round(100 * mean(count_b > 0, na.rm = TRUE), 1),
                count_pct_p05 = round(100 * mean(count_p < 0.05, na.rm = TRUE), 1),
                log1p_median = round(median(log1p_b, na.rm = TRUE), 5),
                log1p_pct_p05 = round(100 * mean(log1p_p < 0.05, na.rm = TRUE), 1)),
            by = source])
  cat("\n-- by instrument --\n")
  print(all[, .(count_b = round(median(count_b, na.rm = TRUE), 5),
                count_p05 = round(100 * mean(count_p < 0.05, na.rm = TRUE), 1),
                log1p_b = round(median(log1p_b, na.rm = TRUE), 5),
                log1p_p05 = round(100 * mean(log1p_p < 0.05, na.rm = TRUE), 1)),
            by = .(source, zset)][order(source, zset)])
}
close(LOGF)
