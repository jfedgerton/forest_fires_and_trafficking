# MANIFEST — what is in this repository and what is not

**Written 13 October 2026, after the cloud working directory was reclaimed.**

The analysis container holding the working copy of this project was wiped some
time between 29 September and 6 October 2026. This repository is a
reconstruction. Read this file before trusting anything in `code/`.

---

## Status of every file

### Recovered verbatim

These were written during the September 2026 sessions and their full text was
still in the conversation transcript, so they are exact rather than
paraphrased. They have NOT been re-executed since recovery.

| File | What it does | Confidence |
|---|---|---|
| `code/12_climate.py` | TerraClimate pull, 13 variables, municipality-year, four seasonal windows | High — full text plus all four later edits |
| `code/degree_days.py` | Schlenker-Roberts GDD/EDD from on-disk daily ERA5 | High — full text plus two later edits |
| `code/13_weather_placebo.R` | Positive-control battery, rainfall specification ladder, three clusterings | Medium — several patches applied by string replacement; diff against any surviving copy |
| `code/14_poisson_iv_gmm.R` | Mullahy multiplicative-error GMM, Wooldridge quasi-differencing, two synthetic tests | Medium — same caveat |
| `code/15_plausibly_exogenous.R` | Conley-Hansen-Rossi bounds with van Kippersluis-Rietveld zero-first-stage calibration | High — single version |
| `code/16_seasonal_windows.R` | Five quarter-specific VPD windows against agriculture, fire and rescues | High — single version |
| `code/17_channel_triage.R` | Crop-year alignment, contamination gradient, channel sign | High — single version |
| `code/first_stage_diagnostic.R` | Reduced form vs first stage by subsample, Anderson-Rubin sets | High — single version |
| `code/patch_metrics_500m.py` | Native-resolution burn-scar metrics, 463m MODIS sinusoidal | High — read in full early in the session |
| `code/run_all.sh` | Stage runner, final numbering | High — final version |

### NOT recovered — must come from your own copy

**This is the bulk of the pipeline and it is the part that cannot be
reconstructed from the transcript.** Only fragments were ever read into
context, and a half-correct `00_config.R` is more dangerous than none.

| File | Why it matters |
|---|---|
| `00_config.R` | Every constant, control set, the `theme_ff` plotting theme, `first_stage_F`, `apply_transform`, `load_panel`. Everything else sources it. |
| `01_acquire.py` | ERA5, MODIS MCD64A1, SIDRA acquisition |
| `02_build_panel.py` | Panel construction, VPD computation, outcome merge, patch metrics |
| `03_screen_controls.R` | Agricultural control screen |
| `04_models.R` | Main models and the 2,016-cell grid |
| `05_robustness.R` | Transform sweep, power, subsample, jackknife, ignition-vs-escape |
| `06_tables.R`, `07_figures.R`, `08_maps.py`, `09_latex_assets.py`, `10_deck.js` | Output generation |
| `11_poisson_bootstrap_OPTIONAL.R` | FE-Poisson control function with cluster bootstrap |

### Data — all gone, none of it belongs here anyway

| | Cost to rebuild |
|---|---|
| ERA5 daily temperature and dewpoint, 2002-2024 | ~44 hours of Copernicus downloads |
| MODIS MCD64A1 monthly burned area | ~2 hours |
| SIDRA 5457 / 3939 / 5938 / 6579 | ~10 minutes |
| `final_panel_muni_quarter.csv.gz` | ~7 hours of processing on top of the above |
| `climate_muni_year.csv` | 31 minutes |
| `degree_days_muni_year.csv` | 4 minutes |
| `patch_metrics_500m_muni_quarter.csv` | ~40 minutes, 3.1M scars |
| All `output/results/*.csv`, figures, the slide deck | derived |

---

## The findings survived

Everything substantive from September is in `docs/`, carried over from the
claude.ai project. The numbers in those documents were produced by code that no
longer exists in runnable form, so treat them as a record of what was found
rather than as reproducible output until the pipeline is restored.

Most important of these is `exclusion-restriction-failure-2026-09-23.md`, which
holds the frequency comparison (median |t| 8.09 at municipality-year against
0.68 at municipality-quarter), the clustering scorecard, and the channel checks.

---

## What has to happen next, in order

1. **Locate the original pipeline.** ROAR, GitHub, Dropbox, a local machine.
   Scripts 00-11 have to come from there; they are not recoverable here.
2. **Reconcile.** Diff the recovered files above against whatever surviving copy
   exists. Where they disagree, the surviving copy wins — these are
   reconstructions.
3. **Re-authenticate GitHub.** The token in this environment is rejected, so
   nothing can be pushed yet.
4. **Then** decide whether to re-acquire the data or work from a surviving panel.

---

## Note on what this repository is for

Everything here is code and prose. No data files, no outputs, no intermediate
artifacts. That separation is the point — the thing that was lost was recoverable
in principle precisely to the extent it had been written down rather than only
computed.
