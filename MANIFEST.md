# MANIFEST — what is in this repository

**Written 13 October 2026, after the cloud working directory was reclaimed and
the pipeline was recovered from Dropbox.**

---

## What happened

The analysis container holding the working copy was wiped some time between 29
September and 6 October 2026. Everything in it was lost — code, built panel,
climate extracts, outputs, slide deck.

The pipeline was then recovered intact from
`Dropbox/forest_fires_trafficking/code/`, last synced 14–16 September 2026. File
sizes and content hashes match the pre-loss state exactly.

---

## Contents

### Recovered from Dropbox, byte-identical

The September 2026 pipeline, in its **pre-consolidation numbering**.

| | |
|---|---|
| `00_config.R` | Constants, control sets, `theme_ff`, `save_fig`, `load_panel`, estimation helpers. Everything sources this. |
| `01_acquire_era5.py` | ERA5 daily temperature and dewpoint |
| `02_acquire_fire_mcd64.py` | MODIS MCD64A1 burned area |
| `03_acquire_sidra.py` | IBGE SIDRA 5457 / 3939 / 5938 / 6579 |
| `04_build_fire_and_vpd.py` | Fire aggregation and VPD construction |
| `05_build_outcomes.py` | TF and UE outcome build |
| `06_build_analysis_panel.py` | Municipality-quarter panel assembly |
| `07_screen_controls.R` | Agricultural control screen |
| `08_main_models.R` | Main 2SLS models |
| `09_specification_grid.R` | Specification grid |
| `10_tables.R`, `11_figures.R`, `12_maps.py` | Output generation |
| `13_poisson_bootstrap_OPTIONAL.R` | FE-Poisson control function, cluster bootstrap |
| `14_tables_meeting.R`, `15_diagnostics_power.R`, `16_robustness_sweep.R` | Meeting outputs, power, transform sweep |
| `19_figures_deck.R`, `20_figures_dark.R`, `22_render_equations.py` | Deck assets |
| `run_all_ORIGINAL.sh` | The original stage runner |

### Reconstructed from the session transcript — NOT byte-verified

Written after the 14–16 September Dropbox sync, so they exist nowhere else.
Treat as drafts until executed.

| | |
|---|---|
| `15_plausibly_exogenous.R` | Conley-Hansen-Rossi bounds, γ calibrated from a zero-first-stage sample (van Kippersluis & Rietveld 2018) |
| `16_seasonal_windows.R` | Five quarter-specific VPD windows against agriculture, fire and rescues |
| `17_channel_triage.R` | Crop-year alignment, contamination gradient, channel sign |

**Numbering collision.** The reconstructed `15` and `16` collide with the
recovered `15_diagnostics_power.R` and `16_robustness_sweep.R`. Renumber before
running anything. This is deliberate — a silent overwrite would be worse.

One bug already found and fixed by comparing against the recovered config: the
reconstructions called `save_fig(g, stub, w =, h =)` but the real signature is
`save_fig(g, stub, width, height)`.

### Lost and not recovered

Written in the container after the last Dropbox sync.

| | Status |
|---|---|
| The consolidated 00–11 pipeline (25 scripts merged into 12, ~40 name collisions resolved) | Gone. The pre-consolidation versions above are the same science in a messier layout. |
| `12_climate.py`, `degree_days.py` | TerraClimate pull and Schlenker-Roberts degree days. Full text is in the session transcript and can be rewritten. |
| `13_weather_placebo.R`, `14_poisson_iv_gmm.R` | Recoverable from transcript but heavily patched; reconstruction would be approximate. |
| `patch_metrics_500m.py` | Native-resolution burn-scar metrics. Recoverable from transcript. |
| `first_stage_diagnostic.R` | Recoverable from transcript. |
| `10_deck.js` | The 36-slide deck builder. Not recoverable. |

### Data — deliberately absent

No data in this repository, by `.gitignore`. Source files live in Dropbox.

| | Where | Cost to rebuild |
|---|---|---|
| TF / UE outcome CSVs, 1995–2025 | `raw_data/` | **Irreplaceable.** Ministry of Labour records. |
| `rescues_by_source.dta`, `inspections_by_mun_month.dta` | project root | Irreplaceable |
| ERA5 daily | — | ~44 h of Copernicus downloads |
| MODIS MCD64A1 | — | ~2 h |
| SIDRA tables | — | ~10 min |
| `final_panel_muni_quarter.csv.gz` | — | ~7 h on top of the above |
| TerraClimate, degree days, 500 m scars | — | ~75 min total |

---

## Findings

`docs/` holds the September analysis notes carried from the claude.ai project.
The numbers in them were produced by code that is only partly runnable now, so
read them as a record of what was found rather than as reproducible output.

`exclusion-restriction-failure-2026-09-23.md` is the important one — the
frequency comparison (median |t| 8.09 at municipality-year against 0.68 at
municipality-quarter), the clustering scorecard, and the channel checks.

---

## Before running anything

1. Renumber the three reconstructed scripts to avoid the 15/16 collision.
2. Diff `run_all_ORIGINAL.sh` against the stage list and decide whether to
   re-consolidate or stay with the original numbering.
3. Point `DIR_DATA` at wherever the panel is restored to.
4. The reconstructed scripts have never executed. Expect them to break.
