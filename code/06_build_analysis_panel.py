"""
06_build_analysis_panel.py
==========================
Assembles the modelling frame: reclaims the TF rows the five harmonization
passes could not place, merges fire, instruments and both outcomes onto the
fixed 2022 IBGE mesh, attaches the controls, and builds the lag structure.

RECLAIM RULES, applied to the 28 leftover rows. Both are flagged in the panel so
every result can be re-run without them.
  Coordinates over UF - nine rows whose parseable GPS falls outside the stated
    UF. Assigned by point-in-polygon, ignoring the UF label. They land somewhere
    sensible (a PB row at -11.21,-41.09 is Morro do Chapeu, BA), which is itself
    evidence the GPS is the reliable field. 32 rescued workers. Transcribed by
    hand rather than regex-parsed: six mutually inconsistent encodings across
    nine rows, where a general parser is likelier to introduce a silent error
    than to save effort.
  Locality names - five rows naming a place that is not a municipality but
    resolves to exactly one: Rio bairros (Tijuca, Jacarepagua -> Rio de
    Janeiro), pre-reform spellings (PARATI -> Paraty, ARES -> Arez) and a
    renamed municipality (Fortaleza do Tabocao -> Tabocao, 2019). 15 workers.
  NOT reclaimed: eight coordinate strings whose encoding admits more than one
    reading, and three names with no unique resolution.

CONTROLS ARE ANNUAL AND HELD CONSTANT WITHIN THE YEAR, not interpolated across
quarters. Interpolation would manufacture within-year variation the source does
not contain and put a smooth trend on the right-hand side of a quarterly model.
All models carry municipality fixed effects, so only between-year movement
identifies anything regardless.

COVERAGE, and why there are two control families:
  2002-2024  crop area, crop value, cattle (SIDRA 5457/3939), population
  2002-2021  everything from SIDRA 5938 - PIB, GDP per capita, VAB shares
`sample_2002_2021` marks the common window so with/without-controls comparisons
can be run on identical rows.

Outputs: analysis_panel_muni_quarter.csv.gz, final_panel_muni_quarter.csv.gz,
         tf_reclaim_log.csv, panel_coverage_report.txt
Seed: 123
"""


# =============================================================================
# STAGE A: reclaim leftover locations and merge the analysis panel
# =============================================================================
import logging
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

np.random.seed(123)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler("30_reclaim_and_merge.log"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

YEAR_MIN, YEAR_MAX = 2002, 2024
N_MUNIS_EXPECTED = 5570
N_QUARTERS_EXPECTED = (YEAR_MAX - YEAR_MIN + 1) * 4  # 92

# ----------------------------------------------------------------------------
# 0. Reference geography: the fixed 2022 IBGE mesh
# ----------------------------------------------------------------------------
mesh = gpd.read_file("malha_inter.json").set_crs("EPSG:4326", allow_override=True)
mesh["code_muni_7"] = mesh["codarea"].astype(int)
mesh = mesh[["code_muni_7", "geometry"]]
assert len(mesh) == N_MUNIS_EXPECTED, f"mesh has {len(mesh)} munis, expected {N_MUNIS_EXPECTED}"

ibge = pd.read_pickle("ibge.pkl")[["code7", "name", "uf"]].rename(columns={"code7": "code_muni_7"})
MESH_CODES = set(mesh.code_muni_7)
log.info("mesh loaded: %d municipalities", len(mesh))

# ----------------------------------------------------------------------------
# PART A1. Coordinate rows whose GPS contradicts the stated UF
# ----------------------------------------------------------------------------
# Parsed by hand from tf_unmatched_report.csv rather than by regex: there are
# 17 coordinate strings in 6 mutually inconsistent encodings, and a general
# parser on 17 rows is more likely to introduce a silent error than to save
# effort. Each entry below is (stated_UF, raw_string, lat, lon, n_rescued).
# Rows whose encoding is genuinely ambiguous are listed but set to None and
# are NOT reclaimed - see AMBIGUOUS below.

COORD_RECLAIM = [
    # stated_uf, raw location string,                                 lat,        lon,        rescued
    ("MG", "S 20° 02' 59,6\" LONG W 051° 03' 38.9\"",                 -20.049889, -51.060806, 18.0),
    ("PB", "-11.213878, -41.093041",                                  -11.213878, -41.093041,  7.0),
    ("PB", "8°29'46.0\"S 36°57'00.6\"W",                               -8.496111,  -36.950167,  4.0),
    ("PB", "11°15'14.6\"S 41°05'15.8\"W",                             -11.254056, -41.087722,  3.0),
    ("MG", "22°10'42.58\"S 53°52'46.52\"W",                           -22.178494, -53.879589,  0.0),
    ("DF", "17° 20' 3'' S 47° 49' 50'' W",                            -17.334167, -47.830556,  0.0),
    ("PA", "25°22'32.1\"S 52°18'31.4\"W",                             -25.375583, -52.308722,  0.0),
    ("PR", "4°29'40\" S 48°7'18\" W",                                  -4.494444,  -48.121667,  0.0),
    ("PA", "04°25'53.77\"S 47°14'40.00\"W",                            -4.431603,  -47.244444,  0.0),
]

# Encodings that admit more than one reading -> left unmatched on purpose.
AMBIGUOUS = [
    ("MA", "`-6,1828S, -434612W",                    17.0, "lon '-434612' reads as -43.4612 or -43°46'12\"=-43.7700"),
    ("MA", "-6,5166S/-444586W",                       2.0, "lon '-444586' reads as -44.4586 or -44°45'86\" (86s invalid)"),
    ("CE", "3º.35’55.39” e 38º58’39,40”",             1.0, "no hemisphere markers; stray '.' after degrees"),
    ("CE", "46' 19.538\" S, 38 33'9.981\" W",         1.0, "latitude degrees field missing entirely"),
    ("BA", "12O 99’72’’ S 38O 51’87’’ W",             0.0, "99 minutes / 72 seconds - not a valid sexagesimal"),
    ("ES", "19°24960=S 39°5949=O",                    0.0, "digits unseparated; no minute/second delimiters"),
    ("CE", "latitude 3,40,27,17839 Longitude 38’49’40,10879", 0.0, "four comma-separated fields, unclear order"),
    ("RJ", "-65.95739999999999",                      0.0, "single number - no latitude/longitude pair"),
]

log.info("PART A1: point-in-polygon for %d parseable coordinate rows", len(COORD_RECLAIM))
pts = gpd.GeoDataFrame(
    pd.DataFrame(COORD_RECLAIM, columns=["stated_uf", "location", "lat", "lon", "rescued"]),
    geometry=[Point(lo, la) for _, _, la, lo, _ in COORD_RECLAIM],
    crs="EPSG:4326",
)
hit = gpd.sjoin(pts, mesh, how="left", predicate="within").drop(columns="index_right")
hit = hit.merge(ibge, on="code_muni_7", how="left")
hit["landed_uf"] = hit["uf"]
hit["uf_conflict"] = hit.landed_uf.notna() & (hit.landed_uf != hit.stated_uf)
hit["rule"] = np.where(hit.code_muni_7.isna(), "coord_outside_brazil",
                       np.where(hit.uf_conflict, "A1_coord_over_uf", "A1_coord_uf_agrees"))

for r in hit.itertuples():
    log.info("  %-2s %-42s -> %s %s (%s) [%s] rescued=%.0f",
             r.stated_uf, r.location[:42],
             "None" if pd.isna(r.code_muni_7) else int(r.code_muni_7),
             r.name if isinstance(r.name, str) else "-",
             r.landed_uf if isinstance(r.landed_uf, str) else "-",
             r.rule, r.rescued)

# ----------------------------------------------------------------------------
# PART A2. Locality names that are not municipality names
# ----------------------------------------------------------------------------
# Each of these is a named place that does resolve to exactly one municipality
# in the 2022 mesh, but could not match on name because the raw string is a
# neighbourhood (bairro), a pre-reform spelling, or a former municipality name.
# Hard-coded with the reason; nothing here was inferred by string distance.

NAME_RECLAIM = [
    # stated_uf, raw location,        code_muni_7, ibge_name,          rescued, reason
    ("RJ", "Tijuca,",                 3304557, "Rio de Janeiro",  14.0, "bairro of Rio de Janeiro municipality"),
    ("RJ", "Jacarepaguá",             3304557, "Rio de Janeiro",   1.0, "bairro of Rio de Janeiro municipality"),
    ("RJ", "PARATI",                  3303807, "Paraty",           0.0, "pre-1970s orthography; fuzzy ratio 83 < 85 cutoff"),
    ("RN", "ARÊS",                    2401206, "Arez",             0.0, "s/z variant; only Are* municipality in RN besides Areia Branca"),
    ("TO", "FORTALEZA DO TABOCÃO",    1708254, "Tabocão",          0.0, "municipality renamed Tabocão (2019); no other Taboc* in TO"),
]

# Deliberately NOT reclaimed - no unique municipality resolution exists:
UNRESOLVED_NAMES = [
    ("RN", "Salina do Finado Laire", 4.0, "salt works, not a populated place; no municipality named"),
    ("RJ", "MARIA IZABEL",           0.0, "no municipality or known locality of this name in RJ"),
    ("RJ", "Taquara",                0.0, "bairro of Rio de Janeiro AND a municipality in RS - ambiguous"),
]

log.info("PART A2: %d name-based reclaims", len(NAME_RECLAIM))
for uf, loc, code, nm, resc, why in NAME_RECLAIM:
    assert code in MESH_CODES, f"{code} not in 2022 mesh"
    log.info("  %-2s %-24s -> %d %s | %s", uf, loc, code, nm, why)

# ----------------------------------------------------------------------------
# PART A3. Rebuild the TF municipality-quarter outcome with reclaimed rows
# ----------------------------------------------------------------------------
tf_rows = pd.read_csv("tf_cleaned_rowlevel.csv")
log.info("TF row-level input: %d rows", len(tf_rows))

loc_col = "location" if "location" in tf_rows.columns else tf_rows.columns[
    [c.lower().startswith("loc") for c in tf_rows.columns]].tolist()[0]
uf_col = "UF" if "UF" in tf_rows.columns else "uf"

# Build a lookup: (stated_uf, raw location) -> (code, reclaim_flag)
reclaim_map = {}
for r in hit.itertuples():
    if r.rule.startswith("A1") and pd.notna(r.code_muni_7):
        reclaim_map[(r.stated_uf, r.location)] = (int(r.code_muni_7),
                                                  "coord_over_uf" if r.uf_conflict else "coord_uf_agrees")
for uf, loc, code, nm, resc, why in NAME_RECLAIM:
    reclaim_map[(uf, loc)] = (code, "name_locality")

def apply_reclaim(row):
    key = (row[uf_col], row[loc_col])
    return pd.Series(reclaim_map.get(key, (np.nan, "")), index=["reclaim_code", "reclaim_flag"])

missing = tf_rows["code_muni_7"].isna()
rec = tf_rows.loc[missing].apply(apply_reclaim, axis=1)
tf_rows.loc[missing, "code_muni_7"] = rec["reclaim_code"]
tf_rows["reclaim_flag"] = ""
tf_rows.loc[missing, "reclaim_flag"] = rec["reclaim_flag"].fillna("")
n_rec = (tf_rows.reclaim_flag != "").sum()
log.info("reclaimed %d TF rows", n_rec)

tf_ok = tf_rows[tf_rows.code_muni_7.notna()].copy()
tf_ok["code_muni_7"] = tf_ok.code_muni_7.astype(int)
assert set(tf_ok.code_muni_7) <= MESH_CODES, "TF has codes outside the 2022 mesh"

resc_col = [c for c in tf_ok.columns if "rescued" in c.lower()][0]
tf_ok["is_reclaimed"] = (tf_ok.reclaim_flag != "").astype(int)
tf_ok["rescued_reclaimed"] = tf_ok[resc_col] * tf_ok.is_reclaimed

tf_q = (tf_ok.groupby(["code_muni_7", "year", "quarter"], as_index=False)
              .agg(workers_rescued_tf=(resc_col, "sum"),
                   n_operations=(resc_col, "size"),
                   workers_rescued_tf_reclaimed=("rescued_reclaimed", "sum"),
                   n_reclaimed_rows=("is_reclaimed", "sum")))
tf_q["workers_rescued_tf_strict"] = tf_q.workers_rescued_tf - tf_q.workers_rescued_tf_reclaimed
tf_q.to_csv("tf_muni_quarter_v2.csv", index=False)
log.info("tf_muni_quarter_v2.csv: %d rows | rescued total %.0f (strict %.0f, reclaimed %.0f)",
         len(tf_q), tf_q.workers_rescued_tf.sum(),
         tf_q.workers_rescued_tf_strict.sum(), tf_q.workers_rescued_tf_reclaimed.sum())

# audit trail
rec_log = hit[["stated_uf", "location", "lat", "lon", "landed_uf", "code_muni_7", "name", "rule", "rescued"]].copy()
amb = pd.DataFrame(AMBIGUOUS, columns=["stated_uf", "location", "rescued", "rule"])
unres = pd.DataFrame(UNRESOLVED_NAMES, columns=["stated_uf", "location", "rescued", "rule"])
nmr = pd.DataFrame(NAME_RECLAIM, columns=["stated_uf", "location", "code_muni_7", "name", "rescued", "rule"])
pd.concat([rec_log, nmr, amb, unres], ignore_index=True).to_csv("tf_reclaim_log.csv", index=False)

# ----------------------------------------------------------------------------
# PART B. Merge into the analysis panel
# ----------------------------------------------------------------------------
log.info("PART B: building the analysis panel")

fire = pd.read_csv("fire_measures_mcd64_muni_quarter_2002_2024.csv")
vpd = pd.read_csv("vpd_instruments_muni_quarter_2002_2024.csv")
ue = pd.read_csv("ue_muni_quarter.csv")

for nm, d in [("fire", fire), ("vpd", vpd)]:
    assert d.code_muni_7.nunique() == N_MUNIS_EXPECTED, f"{nm}: wrong muni count"
    assert set(d.code_muni_7) == MESH_CODES, f"{nm}: mesh mismatch"

# skeleton: every municipality x every quarter (balanced)
skel = pd.MultiIndex.from_product(
    [sorted(MESH_CODES), range(YEAR_MIN, YEAR_MAX + 1), range(1, 5)],
    names=["code_muni_7", "year", "quarter"]).to_frame(index=False)
assert len(skel) == N_MUNIS_EXPECTED * N_QUARTERS_EXPECTED

panel = skel.merge(
    fire[["code_muni_7", "year", "quarter", "municipality_area_km2",
          "total_burned_area_sqkm", "area_burned_no_double_sqkm", "percent_area_burned"]],
    on=["code_muni_7", "year", "quarter"], how="left")
panel = panel.merge(
    vpd[["code_muni_7", "year", "quarter", "vpd_mean_q", "vpd_max_q", "vpd_excess3_q"]],
    on=["code_muni_7", "year", "quarter"], how="left")

# outcomes: absent muni-quarter = no rescue / no victim recorded = 0
tfq = tf_q[(tf_q.year >= YEAR_MIN) & (tf_q.year <= YEAR_MAX)]
panel = panel.merge(tfq, on=["code_muni_7", "year", "quarter"], how="left")
panel = panel.merge(ue[["code_muni_7", "year", "quarter", "victims_ue"]],
                    on=["code_muni_7", "year", "quarter"], how="left")
for c in ["workers_rescued_tf", "workers_rescued_tf_strict", "workers_rescued_tf_reclaimed",
          "n_operations", "n_reclaimed_rows", "victims_ue"]:
    panel[c] = panel[c].fillna(0)

assert panel[["vpd_mean_q", "vpd_max_q", "vpd_excess3_q"]].notna().all().all(), "VPD gaps after merge"
assert panel["percent_area_burned"].notna().all(), "fire gaps after merge"

# ---- outcome forms: count / log1p / binary, for both sources -----------------
panel["tf_count"] = panel.workers_rescued_tf
panel["tf_log1p"] = np.log1p(panel.workers_rescued_tf)
panel["tf_binary"] = (panel.workers_rescued_tf > 0).astype(int)
panel["tf_count_strict"] = panel.workers_rescued_tf_strict
panel["tf_log1p_strict"] = np.log1p(panel.workers_rescued_tf_strict)
panel["tf_binary_strict"] = (panel.workers_rescued_tf_strict > 0).astype(int)
panel["ue_count"] = panel.victims_ue
panel["ue_log1p"] = np.log1p(panel.victims_ue)
panel["ue_binary"] = (panel.victims_ue > 0).astype(int)

# inspection-conditioned subsample flag (main + robustness, per 2026-09-12)
panel["inspected"] = (panel.n_operations > 0).astype(int)

# ---- lags: fire and every instrument, 1 and 2 quarters -----------------------
panel = panel.sort_values(["code_muni_7", "year", "quarter"]).reset_index(drop=True)
LAGGED = ["percent_area_burned", "total_burned_area_sqkm", "area_burned_no_double_sqkm",
          "vpd_mean_q", "vpd_max_q", "vpd_excess3_q"]
g = panel.groupby("code_muni_7", sort=False)
for c in LAGGED:
    panel[f"{c}_lag1"] = g[c].shift(1)
    panel[f"{c}_lag2"] = g[c].shift(2)

# lag sanity: first two quarters of each muni must be the only NaNs
nan_lag1 = panel[f"percent_area_burned_lag1"].isna().sum()
nan_lag2 = panel[f"percent_area_burned_lag2"].isna().sum()
assert nan_lag1 == N_MUNIS_EXPECTED, f"lag1 NaNs {nan_lag1} != {N_MUNIS_EXPECTED}"
assert nan_lag2 == 2 * N_MUNIS_EXPECTED, f"lag2 NaNs {nan_lag2} != {2*N_MUNIS_EXPECTED}"

panel["time_index"] = (panel.year - YEAR_MIN) * 4 + (panel.quarter - 1)
panel["code_muni_6"] = panel.code_muni_7 // 10
panel["uf_code"] = panel.code_muni_7 // 100000

panel.to_csv("analysis_panel_muni_quarter_2002_2024.csv", index=False)
panel.to_csv("analysis_panel_muni_quarter_2002_2024.csv.gz", index=False, compression="gzip")
log.info("analysis panel: %d rows x %d cols", len(panel), panel.shape[1])

# ----------------------------------------------------------------------------
# Validation report
# ----------------------------------------------------------------------------
lines = []
lines.append("PANEL BUILD REPORT")
lines.append(f"rows {len(panel):,} = {panel.code_muni_7.nunique():,} munis x {panel.time_index.nunique()} quarters")
lines.append(f"duplicate keys: {panel.duplicated(['code_muni_7','year','quarter']).sum()}")
lines.append("")
lines.append("OUTCOMES (muni-quarter totals over 2002-2024)")
lines.append(f"  TF workers rescued        {panel.tf_count.sum():>10,.0f}")
lines.append(f"     of which reclaimed     {panel.workers_rescued_tf_reclaimed.sum():>10,.0f}")
lines.append(f"  UE victims                {panel.ue_count.sum():>10,.0f}")
lines.append(f"  muni-quarters w/ TF>0     {panel.tf_binary.sum():>10,} ({100*panel.tf_binary.mean():.3f}%)")
lines.append(f"  muni-quarters w/ UE>0     {panel.ue_binary.sum():>10,} ({100*panel.ue_binary.mean():.3f}%)")
lines.append(f"  muni-quarters inspected   {panel.inspected.sum():>10,} ({100*panel.inspected.mean():.3f}%)")
lines.append("")
lines.append("TREATMENT / INSTRUMENTS (mean, sd)")
for c in LAGGED:
    lines.append(f"  {c:<28} {panel[c].mean():>9.4f}  {panel[c].std():>9.4f}")
report = "\n".join(lines)
open("panel_build_report.txt", "w").write(report + "\n")
print()
print(report)

# =============================================================================
# STAGE B: attach controls and write the final modelling frame
# =============================================================================
import logging
import sys

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    handlers=[logging.FileHandler("34_build_final_panel.log"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

N_MUNIS, N_QUARTERS = 5570, 92

panel = pd.read_csv("analysis_panel_muni_quarter_2002_2024.csv.gz")
econ = pd.read_csv("economic_controls_muni_year.csv")
# NOTE 2026-10-13: 03_acquire_sidra.py writes these WITHOUT the year suffixes and
# run_all_ORIGINAL.sh does not rename, so the pipeline could not run end to end.
# Row counts confirm the suffixes were descriptive, not filters (ag controls are
# already 2002-2021 because SIDRA 5938 stops in 2021; PAM is already 2002-2024).
# Reading whichever name is present.
import os
_ag = ("agricultural_controls_muni_year_2002_2021.csv"
       if os.path.exists("agricultural_controls_muni_year_2002_2021.csv")
       else "agricultural_controls_muni_year.csv")
_pam = ("pam_ppm_agriculture_muni_year_2002_2024.csv"
        if os.path.exists("pam_ppm_agriculture_muni_year_2002_2024.csv")
        else "pam_ppm_agriculture_muni_year.csv")
vab_ag = pd.read_csv(_ag)
pam = pd.read_csv(_pam)
log.info("panel %s | econ %s | vab_ag %s | pam %s",
         panel.shape, econ.shape, vab_ag.shape, pam.shape)

ECON_COLS = ["population", "pib", "gdp_pc", "share_agro", "share_industry",
             "share_services", "share_public", "log_pop", "log_gdp_pc", "log_pib"]
PAM_COLS = ["crop_area_planted_ha", "crop_area_harvested_ha",
            "crop_value_mil_reais", "cattle_head"]

panel = panel.merge(econ[["code_muni_7", "year"] + ECON_COLS],
                    on=["code_muni_7", "year"], how="left")
panel = panel.merge(pam[["code_muni_7", "year"] + PAM_COLS],
                    on=["code_muni_7", "year"], how="left")
# vab agro share from script 31 duplicates econ's share_agro; keep it under a
# distinct name and assert the two agree rather than silently preferring one
panel = panel.merge(
    vab_ag[["code_muni_7", "year", "agri_gdp_share"]].rename(
        columns={"agri_gdp_share": "share_agro_check"}),
    on=["code_muni_7", "year"], how="left")
both = panel.share_agro.notna() & panel.share_agro_check.notna()
maxdiff = (panel.loc[both, "share_agro"] - panel.loc[both, "share_agro_check"]).abs().max()
log.info("share_agro vs independent pull: n=%d  max abs diff %.6g", both.sum(), maxdiff)
assert maxdiff < 1e-6, "the two agricultural-share pulls disagree"
panel = panel.drop(columns="share_agro_check")

assert len(panel) == N_MUNIS * N_QUARTERS, f"row count changed to {len(panel)}"
assert panel.duplicated(["code_muni_7", "year", "quarter"]).sum() == 0

# ---- derived control transforms ---------------------------------------------
# log1p, not log: a municipality with no cropland or no cattle is a real zero,
# not missing, and dropping those rows would select on agricultural activity.
panel["log_crop_area_planted"] = np.log1p(panel.crop_area_planted_ha)
panel["log_crop_area_harvested"] = np.log1p(panel.crop_area_harvested_ha)
panel["log_crop_value"] = np.log1p(panel.crop_value_mil_reais)
panel["log_cattle"] = np.log1p(panel.cattle_head)
panel["crop_area_share"] = panel.crop_area_planted_ha / (panel.municipality_area_km2 * 100.0)
panel["cattle_density"] = panel.cattle_head / panel.municipality_area_km2

panel["sample_2002_2021"] = (panel.year <= 2021).astype(int)
panel["year_quarter"] = panel.year.astype(str) + "Q" + panel.quarter.astype(str)

# ---- coverage report ---------------------------------------------------------
lines = ["FINAL PANEL COVERAGE", f"rows {len(panel):,}", ""]
lines.append(f"{'column':<28}{'non-missing':>14}{'pct':>9}  {'first-last year fully covered'}")
for c in (ECON_COLS + PAM_COLS + ["log_crop_area_planted", "log_cattle",
                                  "crop_area_share", "cattle_density"]):
    nn = panel[c].notna()
    yrs = panel.loc[nn, "year"]
    span = f"{yrs.min():.0f}-{yrs.max():.0f}" if nn.any() else "none"
    lines.append(f"{c:<28}{nn.sum():>14,}{100*nn.mean():>8.1f}%  {span}")

lines.append("")
lines.append("CONTROL SETS (as used in script 35)")
lines.append("  C0 none                     - municipality + year-quarter FE only")
lines.append("  C1 economic                 - log_pop, log_gdp_pc              [2002-2021]")
lines.append("  C2 agricultural, full window- log_crop_area_planted, log_cattle [2002-2024]")
lines.append("  C3 agricultural, VAB share  - share_agro                       [2002-2021]")
lines.append("  C4 C1 + C2                                                     [2002-2021]")
lines.append("  C5 C1 + C2 + C3 + share_industry                               [2002-2021]")
lines.append("")
lines.append("NOTE: share_agro takes negative values in some municipality-years")
lines.append(f"      ({(panel.share_agro < 0).sum():,} rows, min {panel.share_agro.min():.4f}).")
lines.append("      SIDRA reports negative agricultural VAB where subsidies and")
lines.append("      inventory adjustments exceed gross output. Real, not a parse error,")
lines.append("      but it makes share_agro awkward to read as a proportion.")

report = "\n".join(lines)
open("final_panel_coverage_report.txt", "w").write(report + "\n")

# NOTE 2026-10-13: the float64 -> float32 downcast that used to sit here was a
# file-size optimization, but it truncated the treatment, the instrument and the
# log1p outcome transforms to ~7 significant digits and the CSV write then
# serialized them at that precision. Measured effect on headline estimates was
# at most 4.2e-08 relative with no significance flips, so no published number
# changes, but there is no reason to carry it. Integer downcasting is harmless
# and is kept. Cost of removal is 32 MB -> 48 MB gzipped.
for c in panel.columns:
    if panel[c].dtype == "int64":
        panel[c] = pd.to_numeric(panel[c], downcast="integer")

panel.to_csv("final_panel_muni_quarter.csv.gz", index=False,
             compression={"method": "gzip", "compresslevel": 9})
log.info("saved final_panel_muni_quarter.csv.gz: %d rows x %d cols", *panel.shape)
print()
print(report)
