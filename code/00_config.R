################################################################################
# 00_config.R
#
# Shared configuration for every R script in this folder. Source it first:
#     source("00_config.R")
#
# Everything that more than one script needs lives here - paths, the seed, the
# instrument and control definitions, the plotting theme, and the small helpers
# that were previously copy-pasted between scripts. If a constant is used twice,
# it belongs in this file.
################################################################################

set.seed(123)

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
})
if (requireNamespace("patchwork", quietly = TRUE)) suppressPackageStartupMessages(library(patchwork))

setFixest_nthreads(max(1L, parallel::detectCores()))
options(width = 200, scipen = 999)

## ---- paths ------------------------------------------------------------------
DIR_DATA    <- "."            # where the built panels live
DIR_OUT     <- "output"
DIR_TAB     <- file.path(DIR_OUT, "tables")
DIR_FIG     <- file.path(DIR_OUT, "figures")
DIR_RES     <- file.path(DIR_OUT, "results")
for (p in c(DIR_OUT, DIR_TAB, DIR_FIG, DIR_RES)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

PANEL_FILE  <- file.path(DIR_DATA, "final_panel_muni_quarter.csv.gz")
MESH_FILE   <- file.path(DIR_DATA, "malha_inter.json")
SELECTION   <- file.path(DIR_RES, "ag_control_selection.rds")

## ---- panel constants --------------------------------------------------------
N_MUNIS     <- 5570L          # fixed 2022 IBGE mesh, qualidade intermediaria
YEAR_MIN    <- 2002L
YEAR_MAX    <- 2024L
N_QUARTERS  <- (YEAR_MAX - YEAR_MIN + 1L) * 4L
TRUNC_YEAR  <- 2021L          # SIDRA 5938 (GDP, VAB) stops here

## ---- treatment and instruments ----------------------------------------------
FIRE        <- "percent_area_burned"
FIRE_ALT    <- "total_burned_area_sqkm"

VPD_MAIN    <- "vpd_max_q"    # quarterly mean of daily-max VPD - preferred measure
VPD_ALL     <- c("vpd_mean_q", "vpd_max_q", "vpd_excess3_q")

INSTSETS <- list(
  Z1_mean      = VPD_ALL[1],
  Z2_max       = VPD_ALL[2],
  Z3_exc3      = VPD_ALL[3],
  Z4_mean_max  = VPD_ALL[c(1, 2)],
  Z5_mean_exc3 = VPD_ALL[c(1, 3)],
  Z6_max_exc3  = VPD_ALL[c(2, 3)],
  Z7_all       = VPD_ALL
)
LAB_Z <- c(Z1_mean = "Mean VPD", Z2_max = "Mean daily-max VPD",
           Z3_exc3 = "Degrees above 3 kPa", Z4_mean_max = "Mean + max",
           Z5_mean_exc3 = "Mean + excess3", Z6_max_exc3 = "Max + excess3",
           Z7_all = "All three")

LAGSTRUCT <- list(L0 = 0L, L1 = 1L, L2 = 2L, LALL = c(0L, 1L, 2L))
LAB_LAG   <- c(L0 = "Contemporaneous", L1 = "Lag 1 only", L2 = "Lag 2 only",
               LALL = "Contemporaneous + lags 1, 2")

## ---- outcomes ---------------------------------------------------------------
DVS     <- c(TF = "tf_count", UE = "ue_count")
LAB_DV  <- c(TF = "TF - workers rescued", UE = "UE - trafficking victims")
# baseline P(y > 0) in the full panel, used to express effects as % of baseline
BASELINE_RATE <- c(TF = 0.00543, UE = 0.02927)

TRANS    <- c("binary", "count", "log1p", "intensity")
LAB_TRAN <- c(binary = "Binary: Pr(y > 0)", count = "Count (level)",
              log1p = "Log(1 + y)", intensity = "Intensity: log(1+y) | y > 0")

## ---- controls ---------------------------------------------------------------
# The agricultural sets are chosen by 07_screen_controls.R and read from disk so
# that no script hard-codes the screen's outcome. GDP_CTRL is fixed.
GDP_CTRL  <- c("log_pop", "log_gdp_pc")
# any control drawn from SIDRA 5938 forces the truncated window
TRUNC_VARS <- c("log_gdp_pc", "log_pib", "share_agro", "share_industry",
                "share_services", "share_public")

CANDIDATE_CONTROLS <- c(
  log_crop_area_planted   = "Log crop area planted",
  log_crop_area_harvested = "Log crop area harvested",
  log_crop_value          = "Log crop production value",
  log_cattle              = "Log cattle herd",
  crop_area_share         = "Crop area / municipal area",
  cattle_density          = "Cattle per km2",
  share_agro              = "Agriculture share of VAB",
  share_industry          = "Industry share of VAB",
  share_services          = "Services share of VAB",
  share_public            = "Public admin share of VAB",
  log_pop                 = "Log population",
  log_gdp_pc              = "Log GDP per capita"
)
AG_CANDIDATES <- names(CANDIDATE_CONTROLS)[1:7]

build_control_sets <- function(selection_file = SELECTION) {
  if (!file.exists(selection_file))
    stop("run 07_screen_controls.R first - ", selection_file, " not found")
  sel <- readRDS(selection_file)
  st  <- as.data.table(sel$selection_table)
  ag05 <- st[frac_sig_05_tf_count > 0.5 | frac_sig_05_ue_count > 0.5, variable]
  ag10 <- st[frac_sig_10_tf_count > 0.5 | frac_sig_10_ue_count > 0.5, variable]
  list(
    K0_none     = character(0),
    K1_ag05     = ag05,
    K2_ag10     = ag10,
    K3_ag05_gdp = c(ag05, GDP_CTRL),
    K4_ag10_gdp = c(ag10, GDP_CTRL),
    K5_gdp_only = GDP_CTRL
  )
}
LAB_K <- c(K0_none = "No controls", K1_ag05 = "Ag (p<.05)", K2_ag10 = "Ag (p<.10)",
           K3_ag05_gdp = "Ag (.05) + GDP", K4_ag10_gdp = "Ag (.10) + GDP",
           K5_gdp_only = "GDP only")
needs_trunc <- function(v) any(v %in% TRUNC_VARS)

## ---- panel loader -----------------------------------------------------------
load_panel <- function(file = PANEL_FILE) {
  if (!file.exists(file)) stop("run 06_build_analysis_panel.py first - ", file, " not found")
  d <- fread(file)
  d[, muni := factor(code_muni_7)]
  d[, yq   := factor(year_quarter)]
  d[, uf   := factor(uf_code)]
  stopifnot(nrow(d) == N_MUNIS * N_QUARTERS)
  stopifnot(!anyDuplicated(d[, .(code_muni_7, year, quarter)]))
  d[]
}

## ---- estimation helpers -----------------------------------------------------
sfx <- function(L) if (L == 0L) "" else paste0("_lag", L)
rhs <- function(v) if (length(v)) paste(v, collapse = " + ") else "1"
safely <- function(expr) tryCatch(expr, error = function(e) NULL)

# force a length-1 numeric; data.table() given a zero-length column silently
# returns a zero-row table and loses the whole specification
one <- function(x) {
  if (is.null(x) || length(x) != 1L ||
      !is.finite(suppressWarnings(as.numeric(x)))) return(NA_real_)
  as.numeric(x)
}
pick <- function(m, term) {
  if (is.null(m)) return(rep(NA_real_, 4))
  ct <- safely(coeftable(m))
  if (is.null(ct) || !(term %in% rownames(ct))) return(rep(NA_real_, 4))
  as.numeric(ct[term, 1:4])
}

# joint F on the excluded instruments only. fitstat(., "wald") tests every
# non-FE coefficient, so with controls in the model it is not the weak-
# instrument statistic.
first_stage_F <- function(m, zs) {
  one(if (!is.null(m)) safely(
    fixest::wald(m, keep = paste0("^", zs, "$"), print = FALSE)$stat) else NA_real_)
}

# apply an outcome transformation in place; returns the (possibly subset) data
apply_transform <- function(dd, yv, trans) {
  if (trans == "count")          dd[, .y := as.numeric(get(yv))]
  else if (trans == "log1p")     dd[, .y := log1p(get(yv))]
  else if (trans == "binary")    dd[, .y := as.numeric(get(yv) > 0)]
  else if (trans == "intensity") { dd <- dd[get(yv) > 0]; dd[, .y := log1p(get(yv))] }
  else stop("unknown transformation: ", trans)
  dd
}

CLUSTER <- ~ muni + yq
FE_MAIN <- "muni + yq"

## ---- plotting ---------------------------------------------------------------
# Axis text stays horizontal throughout; facet strips sit on top and right,
# never on the left.
theme_ff <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor  = element_blank(),
          axis.text.x       = element_text(angle = 0, hjust = 0.5),
          axis.text.y       = element_text(angle = 0, hjust = 1),
          strip.background  = element_rect(fill = "grey92", colour = NA),
          strip.text        = element_text(angle = 0, lineheight = 1.05),
          strip.text.y      = element_text(angle = 0),
          legend.position   = "bottom",
          plot.title        = element_text(face = "bold"))
}
COL_SIG <- c("p < 0.05" = "#B2182B", "0.05 <= p < 0.10" = "#EF8A62",
             "p >= 0.10" = "grey60")
COL_DV  <- c(TF = "#B2182B", UE = "#2166AC")

save_fig <- function(g, stub, width = 10, height = 8) {
  ggsave(file.path(DIR_FIG, paste0(stub, ".pdf")), g, width = width, height = height)
  ggsave(file.path(DIR_FIG, paste0(stub, ".png")), g, width = width, height = height, dpi = 200)
  invisible(file.path(DIR_FIG, stub))
}

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), paste0(...), "\n")
