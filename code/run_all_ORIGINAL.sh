#!/usr/bin/env bash
# =============================================================================
# run_all.sh - forest fires and forced labour, Brazil
#
#   ./run_all.sh              every stage in order
#   ./run_all.sh 08 10 11     named stages only
#   FORCE=1 ./run_all.sh 11   re-run a stage whose output already exists
#
# Stages are checkpointed and idempotent; an interrupted run resumes by
# re-invoking this script. Acquisition stages (01-03) are rate-limited by
# Copernicus and IBGE, not by this machine, and are written to be re-run until
# they print ALL DONE.
#
# Rough wall-clock with downloads already cached:
#   04 fire + VPD      ~2 h      09 specification grid  ~35 min
#   05 outcomes        ~2 min    10 tables              ~1 min
#   06 panel           ~2 min    11 figures             ~2 min
#   07 screen          ~8 min    12 maps                ~2 min
#   08 main models     ~3 min    13 bootstrap (opt)     ~90 min
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
STAGES=("$@")
run()  { [[ ${#STAGES[@]} -eq 0 ]] || [[ " ${STAGES[*]} " == *" $1 "* ]]; }
have() { [[ -z "${FORCE:-}" && -e "$1" ]]; }
banner() { printf '\n\033[1m=== %s ===\033[0m\n' "$*"; }

run 01 && { banner "01 ERA5 (re-run until ALL DONE)";      python3 01_acquire_era5.py; }
run 02 && { banner "02 MCD64A1 monthly grids";              python3 02_acquire_fire_mcd64.py; }
run 03 && { banner "03 IBGE SIDRA";                         python3 03_acquire_sidra.py; }

run 04 && { banner "04 fire measures and VPD instruments"
            have vpd_instruments_muni_quarter.csv || python3 04_build_fire_and_vpd.py; }
run 05 && { banner "05 TF/UE outcomes"
            have ue_muni_quarter.csv               || python3 05_build_outcomes.py; }
run 06 && { banner "06 analysis panel"
            have final_panel_muni_quarter.csv.gz   || python3 06_build_analysis_panel.py; }

run 07 && { banner "07 screen controls"
            have output/results/ag_control_selection.rds || Rscript 07_screen_controls.R; }
run 08 && { banner "08 main models";                        Rscript 08_main_models.R; }
run 09 && { banner "09 specification grid";                 Rscript 09_specification_grid.R; }
run 10 && { banner "10 tables";                             Rscript 10_tables.R; }
run 11 && { banner "11 figures";                            Rscript 11_figures.R; }
run 12 && { banner "12 maps";                               python3 12_maps.py; }
run 13 && { banner "13 Poisson bootstrap (optional)";       Rscript 13_poisson_bootstrap_OPTIONAL.R; }

banner "pipeline complete"
