"""
03_acquire_sidra.py
===================
Every IBGE SIDRA pull, in one script. Replaces four separate fetchers that
repeated the same request/parse/checkpoint pattern.

TABLES
  5938  PIB and gross value added by sector, municipality-year, 2002-2021.
  5457  Producao Agricola Municipal, crops (temporary + permanent combined,
        classification 782 category 0 = Total), 2002-2024.
  3939  Pesquisa da Pecuaria Municipal, cattle herd (classification 79
        category 2670 = Bovino), 2002-2024.
  6579  Population estimates, 2002-2024 - with holes, see below.
  4709  Censo 2022 population.
  1378  Censo 2010 population.

TWO TRAPS THIS SCRIPT EXISTS TO AVOID

1. SIDRA 5938 VARIABLE CODES. The pre-existing collect_controls_fresh.R
   documented 513 as "PIB per capita", 517 as "agriculture VA" and 525 as
   "total VA". All three are wrong. Verified against the table's own metadata
   endpoint:
       37    PIB a precos correntes
       498   VAB total
       513   VAB AGROPECUARIA          <- the agricultural one
       516   agropecuaria share of total VAB (%)
       517   VAB industria
       525   VAB administracao publica
       6575  VAB servicos
   Anyone re-running the old script gets an INDUSTRY control labelled
   agriculture. The derived agricultural share is cross-checked here against
   IBGE's own published share (516) and the script aborts if they disagree by
   more than 0.5 percentage points; observed worst case is 0.019pp.

2. SIDRA 6579 CENSUS-YEAR HOLES. The population-estimates table returns a
   header and ZERO RECORDS - not an error - for 2007, 2010, 2022 and 2023. A
   naive loop drops those years silently. 2022 and 2010 are filled from the
   censuses; 2007 and 2023 are single-year gaps bracketed by observations and
   are linearly interpolated within municipality. Provenance is recorded per row
   in `pop_source`. Note the 2022 census produced genuine level shifts against
   the prior estimates in 91 municipalities - those are real, not errors.

Checkpointed one CSV per (table, year) under sidra_cache/. Re-run until it
prints ALL DONE; SIDRA is rate-limited in practice to roughly one
municipality-level request every 25-40 seconds.

Outputs
  agricultural_controls_muni_year.csv   code_muni_7, year, vab_agro, vab_total,
                                        pib, agri_gdp_share, agri_share_sidra
  pam_ppm_agriculture_muni_year.csv     crop area planted/harvested, crop value,
                                        cattle head
  population_muni_year.csv              population, pop_source
  economic_controls_muni_year.csv       the assembled economic control set
Seed: 123
"""

import json
import logging
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

np.random.seed(123)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                    handlers=[logging.FileHandler("03_acquire_sidra.log"),
                              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

CACHE = "sidra_cache"
os.makedirs(CACHE, exist_ok=True)

GDP_YEARS = list(range(2002, 2022))    # 5938 stops at 2021
ALL_YEARS = list(range(2002, 2025))

VAB_MAP = {"37": "pib", "543": "taxes", "498": "vab_total", "513": "vab_agro",
           "517": "vab_industry", "6575": "vab_services", "525": "vab_public",
           "516": "agri_share_sidra"}
CROP_MAP = {"8331": "crop_area_planted_ha", "216": "crop_area_harvested_ha",
            "215": "crop_value_mil_reais"}
HERD_MAP = {"105": "cattle_head"}

URL = {
    "5938": "https://apisidra.ibge.gov.br/values/t/5938/n6/all/v/37,543,498,513,517,6575,525,516/p/{y}",
    "5457": "https://apisidra.ibge.gov.br/values/t/5457/n6/all/v/8331,216,215/p/{y}/c782/0",
    "3939": "https://apisidra.ibge.gov.br/values/t/3939/n6/all/v/105/p/{y}/c79/2670",
    "6579": "https://apisidra.ibge.gov.br/values/t/6579/n6/all/v/9324/p/{y}",
    "4709": "https://apisidra.ibge.gov.br/values/t/4709/n6/all/v/93/p/2022",
    "1378": "https://apisidra.ibge.gov.br/values/t/1378/n6/all/v/93/p/2010",
}


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "research-script"})
    with urllib.request.urlopen(req, timeout=240) as r:
        return json.loads(r.read().decode("utf-8"))[1:]   # [0] is the header dict


def fetch_wide(url, varmap, year):
    """Multi-variable tables: long -> wide on the variable code."""
    recs = get(url)
    if not recs:
        return None
    df = pd.DataFrame(recs)[["D1C", "D2C", "V"]]
    df.columns = ["code_muni_7", "var", "value"]
    df["code_muni_7"] = df.code_muni_7.astype(int)
    df["value"] = pd.to_numeric(df.value, errors="coerce")  # '...', '-', 'X' -> NaN
    w = (df.pivot_table(index="code_muni_7", columns="var", values="value",
                        aggfunc="first").rename(columns=varmap).reset_index())
    w["year"] = year
    return w


def fetch_single(url, colname, year):
    recs = get(url)
    if not recs:
        return None
    df = pd.DataFrame(recs)[["D1C", "V"]]
    df.columns = ["code_muni_7", colname]
    df["code_muni_7"] = df.code_muni_7.astype(int)
    df[colname] = pd.to_numeric(df[colname], errors="coerce")
    df["year"] = year
    return df


def cached(tag, year, fn):
    """Fetch once, then read from disk. Returns None if SIDRA has no rows."""
    p = os.path.join(CACHE, f"sidra_{tag}_{year}.csv")
    if os.path.exists(p) and os.path.getsize(p) > 5_000:
        return pd.read_csv(p)
    for attempt in range(3):
        try:
            w = fn()
            if w is None or len(w) < 5000:
                log.info("%s %d: no rows returned", tag, year)
                return None
            w.to_csv(p, index=False)
            log.info("%s %d: %d municipalities", tag, year, len(w))
            return w
        except Exception as e:
            log.error("%s %d attempt %d: %s", tag, year, attempt + 1, str(e)[:140])
            time.sleep(20)
    return None


# ---------------------------------------------------------------- 5938
log.info("=== SIDRA 5938: PIB and VAB by sector, %d-%d ===", GDP_YEARS[0], GDP_YEARS[-1])
gdp = [cached("5938", y, lambda y=y: fetch_wide(URL["5938"].format(y=y), VAB_MAP, y))
       for y in GDP_YEARS]
gdp = pd.concat([g for g in gdp if g is not None], ignore_index=True)

gdp["agri_gdp_share"] = gdp.vab_agro / gdp.vab_total
ok = gdp.agri_gdp_share.notna() & gdp.agri_share_sidra.notna()
diff = (100 * gdp.loc[ok, "agri_gdp_share"] - gdp.loc[ok, "agri_share_sidra"]).abs()
log.info("agricultural share cross-check vs IBGE variable 516: n=%d, max diff %.4f pp",
         ok.sum(), diff.max())
assert diff.max() < 0.5, "derived agricultural share disagrees with IBGE's published share"

gdp["code_muni_6"] = gdp.code_muni_7 // 10
gdp[["code_muni_7", "code_muni_6", "year", "vab_agro", "vab_total", "pib",
     "agri_gdp_share", "agri_share_sidra"]].to_csv(
    "agricultural_controls_muni_year.csv", index=False)

# ---------------------------------------------------------------- 5457 + 3939
log.info("=== SIDRA 5457 crops and 3939 cattle, %d-%d ===", ALL_YEARS[0], ALL_YEARS[-1])
crops = pd.concat([x for x in (
    cached("5457", y, lambda y=y: fetch_wide(URL["5457"].format(y=y), CROP_MAP, y))
    for y in ALL_YEARS) if x is not None], ignore_index=True)
herd = pd.concat([x for x in (
    cached("3939", y, lambda y=y: fetch_wide(URL["3939"].format(y=y), HERD_MAP, y))
    for y in ALL_YEARS) if x is not None], ignore_index=True)
pam = crops.merge(herd, on=["code_muni_7", "year"], how="outer")
pam["code_muni_6"] = pam.code_muni_7 // 10
assert pam.duplicated(["code_muni_7", "year"]).sum() == 0
pam.to_csv("pam_ppm_agriculture_muni_year.csv", index=False)
for c in CROP_MAP.values():
    log.info("  %-24s NaN %6d (%.1f%%)", c, pam[c].isna().sum(), 100 * pam[c].isna().mean())

# ---------------------------------------------------------------- population
log.info("=== population: 6579 estimates plus census patches ===")
frames = []
for y in ALL_YEARS:
    w = cached("6579", y, lambda y=y: fetch_single(URL["6579"].format(y=y), "population", y))
    if w is not None:
        w = w.copy(); w["pop_source"] = "sidra_6579_estimate"; frames.append(w)
have = sorted({int(f.year.iloc[0]) for f in frames})
log.info("6579 returned rows for: %s", have)
log.info("missing (census/count years): %s", [y for y in ALL_YEARS if y not in have])

for year, tag in ((2022, "4709"), (2010, "1378")):
    if year in have:
        continue
    w = cached(tag, year, lambda tag=tag, year=year:
               fetch_single(URL[tag], "population", year))
    if w is not None:
        w = w.copy(); w["pop_source"] = f"censo_{year}"; frames.append(w)

pop = pd.concat(frames, ignore_index=True).dropna(subset=["population"])
full = pd.MultiIndex.from_product([sorted(pop.code_muni_7.unique()), ALL_YEARS],
                                  names=["code_muni_7", "year"]).to_frame(index=False)
pop = full.merge(pop, on=["code_muni_7", "year"], how="left").sort_values(
    ["code_muni_7", "year"])
before = pop.population.isna().sum()
pop["population"] = pop.groupby("code_muni_7")["population"].transform(
    lambda s: s.interpolate(limit_direction="both"))
pop["pop_source"] = pop.pop_source.fillna("interpolated")
log.info("interpolated %d rows; %d still missing", before - pop.population.isna().sum(),
         pop.population.isna().sum())
log.info("source mix:\n%s", pop.pop_source.value_counts().to_string())
assert pop.population.isna().sum() == 0
pop[["code_muni_7", "year", "population", "pop_source"]].to_csv(
    "population_muni_year.csv", index=False)

# ---------------------------------------------------------------- economic set
econ = pop[["code_muni_7", "year", "population", "pop_source"]].merge(
    gdp.drop(columns=["agri_gdp_share", "agri_share_sidra", "code_muni_6"]),
    on=["code_muni_7", "year"], how="outer")
econ["code_muni_6"] = econ.code_muni_7 // 10
econ["gdp_pc"] = econ.pib / econ.population
for sec in ("agro", "industry", "services", "public"):
    econ[f"share_{sec}"] = econ[f"vab_{sec}"] / econ.vab_total
econ["log_pop"] = np.log(econ.population.where(econ.population > 0))
econ["log_gdp_pc"] = np.log(econ.gdp_pc.where(econ.gdp_pc > 0))
econ["log_pib"] = np.log(econ.pib.where(econ.pib > 0))

# Shares must sum to 1 where the decomposition exists. Restrict to complete rows:
# .sum(axis=1) skips NaN and would score a 2022-2024 row (no 5938 data at all) as
# a sum of 0, reporting a deviation of 1.0 that is not a decomposition failure.
SH = ["share_agro", "share_industry", "share_services", "share_public"]
complete = econ[SH].notna().all(axis=1)
s = econ.loc[complete, SH].sum(axis=1)
log.info("sector shares (n=%d complete): mean %.6f, max |dev from 1| %.2g",
         complete.sum(), s.mean(), (s - 1).abs().max())
# tolerance 1e-3, not machine epsilon: SIDRA rounds VAB to whole thousands of reais
assert (s - 1).abs().max() < 1e-3, "VAB sector shares do not sum to 1"

assert econ.duplicated(["code_muni_7", "year"]).sum() == 0
econ.to_csv("economic_controls_muni_year.csv", index=False)
log.info("rows %d | municipalities %d | years %d-%d",
         len(econ), econ.code_muni_7.nunique(), econ.year.min(), econ.year.max())
print("ALL DONE")
