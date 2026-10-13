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
| `23_climate.py` | TerraClimate pull, 13 variables, four seasonal windows, municipality-year |
| `degree_days.py` | Schlenker-Roberts GDD/EDD from the on-disk daily ERA5. Helper to 23. |
| `24_patch_metrics_500m.py` | Native-resolution burn-scar metrics, 463 m MODIS sinusoidal |
| `25_first_stage_diagnostic.R` | Reduced form vs first stage by subsample, Anderson-Rubin sets |
| `26_plausibly_exogenous.R` | Conley-Hansen-Rossi bounds, γ calibrated from a zero-first-stage sample (van Kippersluis & Rietveld 2018) |
| `27_seasonal_windows.R` | Five quarter-specific VPD windows against agriculture, fire and rescues |
| `28_channel_triage.R` | Crop-year alignment, contamination gradient, channel sign |

Numbered from 23 so nothing collides with the recovered pipeline, which runs to
22. The original container numbering was different; this is a renumbering, not a
reordering.

One bug was already caught by diffing against the recovered config: the
reconstructions called `save_fig(g, stub, w =, h =)` but the real signature is
`save_fig(g, stub, width, height)`. Fixed. Expect more of this kind — none of
these seven files has ever executed.

### Lost and not recovered

Written in the container after the last Dropbox sync.

| | Status |
|---|---|
| The consolidated 00–11 pipeline (25 scripts merged into 12, ~40 name collisions resolved) | Gone. The pre-consolidation versions above are the same science in a messier layout. |
| `13_weather_placebo.R` | The positive-control battery and rainfall specification ladder. Recoverable from the transcript but patched many times by string replacement, so a reconstruction would be approximate. **Not attempted.** |
| `14_poisson_iv_gmm.R` | Mullahy multiplicative-error GMM with Wooldridge quasi-differencing, plus two synthetic validation harnesses. Same caveat. **Not attempted.** |
| `10_deck.js` | The 36-slide deck builder. Not recoverable. |

The two un-attempted scripts are the ones where an approximate reconstruction
would be actively dangerous — both carry estimators whose correctness depends on
details that were patched in place. Rewriting them from scratch against the
recovered `00_config.R` is safer than half-remembering them.

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

1. Point `DIR_DATA` in `00_config.R` at wherever the panel is restored to.
2. Rebuild the panel, or restore it from a surviving copy. `00_config.R` asserts
   `nrow == 5570 * 92`, so a partial panel fails loudly rather than quietly.
3. Run the recovered pipeline first and confirm it reproduces the September
   numbers before trusting anything numbered 23 or above.
4. The seven reconstructed scripts have never executed. Expect them to break.
