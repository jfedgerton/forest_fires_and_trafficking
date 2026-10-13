# Forest fires and forced labour in Brazil — replication code

Municipality-quarter panel, 2002–2024, 5,570 municipalities on a fixed 2022 IBGE
mesh. Estimates the effect of area burned on forced-labour rescues (TF) and on
seguro-desemprego trafficking victims (UE), instrumenting fire with vapour
pressure deficit.

This folder replaces the earlier `replication/` folder, which is archived.

```bash
./run_all.sh              # everything, in order
./run_all.sh 08 10 11     # named stages only
FORCE=1 ./run_all.sh 11   # re-run a stage whose output already exists
```

Every stage is checkpointed and idempotent: re-running skips completed work, so
an interrupted run resumes by re-invoking the script.

## Layout

| Stage | # | Script | Does |
|---|---|---|---|
| — | 00 | `00_config.R` | Shared paths, seed, instrument and control definitions, plot theme, estimation helpers. Sourced by every R script. |
| **Data setup** | 01 | `01_acquire_era5.py` | All ERA5 downloads via the CDS async queue. |
| | 02 | `02_acquire_fire_mcd64.py` | MCD64A1 burned area → monthly presence grids. |
| | 03 | `03_acquire_sidra.py` | All IBGE SIDRA pulls: GDP and VAB by sector, crops, cattle, population. |
| **Data merging** | 04 | `04_build_fire_and_vpd.py` | Zonal aggregation → fire measures and the three VPD instruments. |
| | 05 | `05_build_outcomes.py` | Full TF/UE harmonization chain, five passes. |
| | 06 | `06_build_analysis_panel.py` | Reclaim leftovers, merge everything, attach controls, build lags. |
| **Main analyses** | 07 | `07_screen_controls.R` | Balance and relevance screens; defines the ag control sets. |
| | 08 | `08_main_models.R` | The main-text models. |
| **Robustness** | 09 | `09_specification_grid.R` | 2,016-cell specification grid. |
| | 13 | `13_poisson_bootstrap_OPTIONAL.R` | FE-Poisson control function with a cluster bootstrap. |
| **Tables** | 10 | `10_tables.R` | Main-text tables and the instrument-balance table. |
| **Visualizations** | 11 | `11_figures.R` | Binscatter, seasonality, effect sizes, specification curves. |
| | 12 | `12_maps.py` | Choropleths. |

Outputs land in `output/{tables,figures,results}/`.

## Design decisions that matter

**Fixed boundary vintage.** Every panel uses the 2022 IBGE mesh for all years. A
year-switching vintage would make the unit of analysis a function of time and
break the fixed-effect panel.

**VPD is built from daily values, then averaged** — never `es()` of an averaged
temperature. `es()` is convex in T, so averaging temperature first understates
mean VPD by an amount scaling with Var(T), which varies by biome and season.
That would be a non-classical measurement error correlated with treatment.

| Instrument | Definition |
|---|---|
| `vpd_mean_q` | quarterly mean of daily `es(t2m_mean) − es(d2m_mean)` |
| `vpd_max_q` | quarterly mean of daily `es(t2m_max) − es(d2m_mean)` — **preferred** |
| `vpd_excess3_q` | quarterly mean of daily `max(0, VPDmax − 3.0 kPa)` |

**Everything is 2SLS, including the count outcome.** Linear IV takes a count on
the left-hand side without difficulty — the coefficient is a level effect rather
than a semi-elasticity, and the clustered analytic SE is valid. FE-Poisson would
force a control function whose analytic SE ignores the first stage; that route is
stage 13, a robustness check, not the main specification.

**Controls are annual and held constant within the year**, not interpolated.

## Four traps this code exists to avoid

**1. SIDRA 5938 variable codes are mis-documented in the old pipeline.**
`collect_controls_fresh.R` had 513 as "PIB per capita", 517 as "agriculture VA",
525 as "total VA". All three are wrong:

| Code | Actually is |
|---|---|
| 37 | PIB a preços correntes |
| 498 | VAB total |
| **513** | **VAB agropecuária** |
| 516 | agropecuária share of total VAB (%) |
| 517 | VAB indústria |
| 525 | VAB administração pública |

Anyone re-running the old script gets an **industry** control labelled
agriculture. Stage 03 cross-checks its derived share against IBGE's published
variable 516 and aborts on disagreement above 0.5pp (observed worst case
0.019pp).

**2. SIDRA 6579 returns zero rows, not an error, for census years.** 2007, 2010,
2022 and 2023 are missing from the population estimates table. A naive loop drops
them silently. Stage 03 fills 2022 and 2010 from the censuses and interpolates
the two remaining single-year gaps, recording provenance in `pop_source`. The
2022 census produced genuine level shifts against prior estimates in 91
municipalities — real, not errors.

**3. Lag-only models are misspecified by omission.** `L1` and `L2` enter a lagged
fire term *without* the contemporaneous one. Within municipality and net of the
fixed effects, `corr(fire_t, fire_t−2) = −0.22`, so an L2-only term proxies for
*low* current fire and flips sign mechanically. Stage 09 estimates them for
completeness and flags them `lag_only = 1`; stages 11 and the headline tables
exclude them.

**4. Three grid rows used to be the same regression.** `binary`,
`hurdle_count`'s participation equation and `hurdle_log1p`'s participation
equation were byte-identical linear probability models — verified, maximum
coefficient difference exactly 0.0 across 336 matched pairs. Pooling without
deduplicating triple-counted the participation margin, which is precisely the
margin that is significant, and inflated every pooled rate. There is now one
binary row.

## The window / controls confound, and why every table has three columns

The full control set draws on SIDRA 5938, which ends in 2021, so every
with-controls column is 2002–2021. If the no-controls column were run on
2002–2024, "controls change the answer" would be inseparable from "the last three
years change the answer". Every main table therefore carries **no controls
2002–2024**, **no controls 2002–2021**, and **full controls 2002–2021**.

Decomposed, the answer is unambiguous: holding the window fixed, controls move
the TF binary estimate from 0.00062 to 0.00062. The window does all the work
that had appeared to be controls.

## Headline results

Contemporaneous fire, `vpd_max_q`, municipality and year-quarter FE, two-way
clustered SEs. First-stage F ≈ 99–102 throughout.

| DV | Outcome | No controls 2002–24 | No controls 2002–21 | Full controls 2002–21 |
|---|---|---|---|---|
| TF | Binary | +0.00079 (p=.0002) | +0.00062 (p=.0016) | +0.00062 (p=.0009) |
| TF | Count | +0.01734 (p=.022) | +0.01644 (p=.040) | +0.01607 (p=.030) |
| TF | Log(1+y) | +0.00205 (p=.0004) | +0.00170 (p=.0025) | +0.00170 (p=.0012) |
| UE | Binary | +0.00255 (p=.015) | +0.00255 (p=.016) | +0.00246 (p=.015) |
| UE | Count | +0.01116 (p=.027) | +0.01119 (p=.031) | +0.01093 (p=.032) |
| UE | Log(1+y) | +0.00317 (p=.016) | +0.00319 (p=.016) | +0.00310 (p=.016) |

A one-standard-deviation increase in area burned (5.77pp) raises the probability
of any TF rescue by about 66% of its baseline and any UE victim by about 48%.

**Fire acts on the extensive margin only.** The intensity margin — log(1+y)
conditional on y>0 — is never significant, in either DV, under any instrument or
control set. Fire changes *whether* a municipality-quarter records a case, not
*how many* once it does.

**Two things that do not hold.** Lag specifications: TF keeps its sign but loses
significance in most specifications; UE is a null at every window and control set
(all p > 0.47). And the OLS/IV gap: the naive estimate is null in 126 of 126
binary specifications while the IV is 23× larger for TF and sign-flipped for UE.
Measurement error in MCD64A1 attenuates OLS, but attenuation of that magnitude
implies a signal-to-noise ratio around 4%. That gap needs an explanation in the
manuscript.

## Verification status

| Stage | Status |
|---|---|
| 03, 07, 08, 09, 10, 11, 12 | Run end-to-end in this environment; outputs checked |
| 06 stage B | Run and verified |
| 01, 02, 04, 05, 06 stage A | Structural consolidations of code that produced the current panel. Not re-run end-to-end here — 01 would re-download ~44 hours of ERA5 and 02 about 90 minutes. The merged files parse and preserve the original sequence, but **verify against a known-good panel before trusting a fresh build.** |

`ENVIRONMENT.md` has pinned versions. Seed is 123 throughout.
