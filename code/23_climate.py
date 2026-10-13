#!/usr/bin/env python3
"""
23_climate.py - TerraClimate monthly climate variables, municipality-year.

RECONSTRUCTED 13 October 2026 from the session transcript after the working
directory was lost. Never re-executed since. Verify before trusting.

WHY THIS EXISTS. The exclusion-restriction evidence in this paper is a set of
nulls: VPD does not move any of the agricultural controls. A null is only
evidence if the test could have detected something. Nothing in the pipeline
established that it could, because the panel contains no weather variable other
than VPD and its own two inputs.

This script adds the variables needed for a POSITIVE CONTROL: climate measures
that agronomy says must predict agricultural outcomes. If they do, the VPD null
is informative. If they do not, the test has no power and the exclusion
restriction is unsupported by it - which is a finding, not a failure.

THE CALIBRATION PROBLEM. Rainfall is not monotonic in crop outcomes: too little
is drought, too much is waterlogging and harvest loss. A linear precipitation
test can return a null purely from misspecification, which would be exactly the
wrong lesson. Three of the variables pulled here sidestep that entirely because
the non-monotonicity is already resolved by the physics:

  pdsi  Palmer Drought Severity Index - water supply netted against demand
  def   climatic water deficit, PET - AET: water the plants wanted and did not get
  aet   actual evapotranspiration: water actually transpired

Those are monotonic in plant water stress by construction. They are the
insurance against getting the rainfall functional form wrong. Raw ppt is pulled
too, and the placebo script runs it through a specification ladder.

SOURCE. TerraClimate via Microsoft Planetary Computer, the same access path
02_acquire_fire_mcd64.py already uses for MCD64A1. Monthly, 1/24 degree (~4 km),
and the MPC copy runs 1958-01 through 2021-12. That end date is fine for the
estimation window (TRUNC_YEAR = 2021) but means these variables CANNOT support
the 2022-2024 extension if the GDP control is later dropped. Noted, not worked
around.

INDEPENDENT VPD. TerraClimate computes its own VPD from its own inputs. It is
pulled deliberately: correlating it against the ERA5-derived vpd_max_q at
municipality-year is an external check on the instrument's construction. When
this was run in September the two correlated at r = 0.834 in levels but only
r = 0.593 within municipality, which is lower than two measures of one physical
quantity should be and remains unexplained.

AGGREGATION WINDOWS. Brazilian soy runs roughly October-March, so the
agricultural year straddles calendar years. Four windows are built per
municipality-year y:

  cal    Jan(y) - Dec(y)          calendar year, comparable to the SIDRA year
  agyr   Oct(y-1) - Sep(y)        agricultural year
  grow   Oct(y-1) - Mar(y)        main growing season, first crop
  harv   Feb(y) - Apr(y)          harvest window, when excess rain does damage

Flux variables (ppt, aet, pet, def, q, srad) are SUMMED within a window; state
variables (pdsi, soil, tmax, tmin, vpd, vap, ws) are AVERAGED. Summing a state
variable or averaging a flux would both be wrong and neither would be obvious
downstream, so the split is declared explicitly in FLUX below and asserted.

STAGES. Mirrors 06_build_analysis_panel.py's pattern:

    python3 23_climate.py --stage terraclimate   TerraClimate only
    python3 23_climate.py --stage degreedays     degree days only (delegates to
                                                 degree_days.py)
    python3 23_climate.py --stage all            both (default)

Output: climate_muni_year.csv, degree_days_muni_year.csv
Checkpointed one parquet per year under climate_tc/ and degdays/, so an
interrupted run resumes. Run from code/.
"""

import glob
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import geopandas as gpd
import planetary_computer
import pystac_client
import xarray as xr
from rasterio.features import rasterize
from rasterio.transform import from_origin

OUTDIR = "climate_tc"
OUT = "climate_muni_year.csv"
MESH = "malha_inter.json"

YEAR_MIN, YEAR_MAX = 2002, 2021      # estimation window
PULL_MIN = YEAR_MIN - 1              # need Oct-Dec of the prior year for agyr/grow

# Every TerraClimate variable. Pulling all of them costs little extra - the
# expensive part is the zonal-mean machinery, not the bytes - and it means a
# later question about radiation or wind does not require another pull.
VARS = ["ppt", "pdsi", "soil", "def", "aet", "pet", "tmax", "tmin",
        "vpd", "srad", "ws", "q", "vap"]

# Flux variables accumulate over time and are summed within a window; everything
# else is a state and is averaged. Getting this backwards is silent and fatal.
FLUX = {"ppt", "aet", "pet", "def", "q", "srad"}
STATE = set(VARS) - FLUX

LON0, LON1, LAT0, LAT1 = -74.05, -33.95, -34.10, 5.90
PX = 0.01            # municipality raster resolution, ~1.1 km

WINDOWS = {
    "cal":  [(0, m) for m in range(1, 13)],
    "agyr": [(-1, m) for m in range(10, 13)] + [(0, m) for m in range(1, 10)],
    "grow": [(-1, m) for m in range(10, 13)] + [(0, m) for m in range(1, 4)],
    "harv": [(0, m) for m in (2, 3, 4)],
}

os.makedirs(OUTDIR, exist_ok=True)


def say(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def retry(fn, what, tries=6, base=5):
    """Exponential backoff around a network call.

    Same policy as 24_patch_metrics_500m.py, and for the same reason: a
    transient proxy outage once burned through 220 months of that script in
    under a second because each failure was cheap. Failure must cost time, not
    data.
    """
    for i in range(tries):
        try:
            return fn()
        except Exception as e:
            if i == tries - 1:
                raise
            wait = base * (2 ** i)
            say(f"    {what}: {type(e).__name__}, retry {i + 1}/{tries - 1} in {wait}s")
            time.sleep(wait)


def municipality_raster():
    """Rasterise the fixed 2022 IBGE mesh at ~1 km over the Brazil bbox.

    Returns the code array plus the lon/lat coordinate vectors of its cell
    centres. Rasterising fine and then averaging whichever TerraClimate cell each
    fine cell falls in gives an AREA-WEIGHTED mean per municipality. Sampling
    TerraClimate at municipality centroids instead would drop every municipality
    smaller than a 4 km cell, and Brazil has a great many of those.
    """
    w = int(round((LON1 - LON0) / PX))
    h = int(round((LAT1 - LAT0) / PX))
    tr = from_origin(LON0, LAT1, PX, PX)
    mesh = gpd.read_file(MESH)
    code_col = next(c for c in mesh.columns
                    if c.lower() in ("codarea", "cd_mun", "code_muni_7"))
    mesh[code_col] = mesh[code_col].astype(np.int64)
    if mesh.crs is None:
        mesh = mesh.set_crs("EPSG:4326")
    say(f"rasterising {len(mesh)} municipalities at {h}x{w} (~1 km)")
    arr = rasterize(((g, c) for g, c in zip(mesh.geometry, mesh[code_col])),
                    out_shape=(h, w), transform=tr, fill=0, dtype="int64")
    n_cov = int((arr > 0).sum())
    n_muni = len(np.unique(arr)) - 1
    say(f"  {n_cov:,} cells assigned, {n_muni:,} municipalities present")
    assert n_muni > 5500, f"only {n_muni} municipalities rasterised, expected ~5570"
    lon = LON0 + (np.arange(w) + 0.5) * PX
    lat = LAT1 - (np.arange(h) + 0.5) * PX
    return arr, lon, lat


def open_terraclimate():
    cat = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace)
    a = cat.get_collection("terraclimate").assets["zarr-abfs"]
    kw = dict(a.extra_fields["xarray:open_kwargs"])
    kw.pop("engine", None)          # xarray>=2024 rejects engine= in open_zarr
    ds = retry(lambda: xr.open_zarr(a.href, **kw), "open terraclimate")
    say(f"terraclimate: {len(ds.data_vars)} vars, "
        f"{str(ds.time.values[0])[:7]} to {str(ds.time.values[-1])[:7]}")
    for v in VARS:
        assert v in ds.data_vars, f"{v} missing from TerraClimate"
    return ds


def build_index(tc_lon, tc_lat, muni, lon, lat):
    """Map every ~1 km municipality cell onto a TerraClimate cell.

    Returns the flat TerraClimate cell index and the municipality code for each
    covered fine cell, so a month reduces to two bincounts.
    """
    ci = np.searchsorted(tc_lon, lon)
    ci = np.clip(ci, 1, len(tc_lon) - 1)
    ci = np.where(np.abs(tc_lon[ci - 1] - lon) <= np.abs(tc_lon[ci] - lon), ci - 1, ci)
    # TerraClimate latitude runs north to south, same as the raster
    ri = np.searchsorted(-tc_lat, -lat)
    ri = np.clip(ri, 1, len(tc_lat) - 1)
    ri = np.where(np.abs(tc_lat[ri - 1] - lat) <= np.abs(tc_lat[ri] - lat), ri - 1, ri)

    ok = muni > 0
    rr = np.repeat(ri[:, None], len(lon), axis=1)[ok]
    cc = np.repeat(ci[None, :], len(lat), axis=0)[ok]
    flat = rr.astype(np.int64) * len(tc_lon) + cc.astype(np.int64)
    codes = muni[ok]
    say(f"  index built: {len(flat):,} fine cells -> "
        f"{len(np.unique(flat)):,} distinct TerraClimate cells")
    return flat, codes


def month_means(block, flat, code_index, n_muni):
    """One month, all variables -> mean per municipality.

    block is a dict var -> 2D array already subset to the Brazil window. Each
    fine cell takes the value of its TerraClimate cell, then values are averaged
    within municipality by bincount. NaN cells (ocean, gaps) are excluded from
    both numerator and denominator per variable, so a coastal municipality is
    averaged over its land cells rather than being poisoned to NaN.

    n_muni is the number of DISTINCT municipalities, not the number of fine
    cells - code_index runs over cells and indexes into the municipality list.
    """
    out = {}
    n = n_muni
    for v, arr2d in block.items():
        vals = arr2d.ravel()[flat]
        good = np.isfinite(vals)
        num = np.bincount(code_index[good], weights=vals[good], minlength=n)
        den = np.bincount(code_index[good], minlength=n)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[v] = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    return out


def pull_year(ds, y, flat, code_index, uniq, lat_sl, lon_sl):
    """Every month of calendar year y -> tidy municipality-month frame."""
    rows = []
    sub = ds.sel(time=slice(f"{y}-01-01", f"{y}-12-31"))
    nt = sub.sizes["time"]
    assert nt == 12, f"{y}: expected 12 months, got {nt}"
    loaded = {}
    for v in VARS:
        loaded[v] = retry(
            lambda v=v: sub[v].isel(lat=lat_sl, lon=lon_sl).values, f"{y} {v}")
    for mi in range(12):
        block = {v: loaded[v][mi] for v in VARS}
        mm = month_means(block, flat, code_index, len(uniq))
        df = pd.DataFrame(mm)
        df.insert(0, "month", mi + 1)
        df.insert(0, "year", y)
        df.insert(0, "code_muni_7", uniq)
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def aggregate(mon):
    """Municipality-month -> municipality-year, four windows.

    Flux variables are summed within a window, state variables averaged. A
    window is only reported when every one of its months is present; a partial
    window is set to NaN rather than quietly summing fewer months, which would
    make 2002's growing season look like a drought.
    """
    out = []
    for y in range(YEAR_MIN, YEAR_MAX + 1):
        parts = {}
        for wname, spec in WINDOWS.items():
            keys = set((y + dy, m) for dy, m in spec)
            sel = mon[[(yy, mm) in keys for yy, mm in zip(mon.year, mon.month)]]
            g = sel.groupby("code_muni_7")
            cnt = g.size()
            full = cnt[cnt == len(keys)].index
            for v in VARS:
                agg = g[v].sum() if v in FLUX else g[v].mean()
                agg = agg.where(agg.index.isin(full))
                parts[f"{v}_{wname}"] = agg
        df = pd.DataFrame(parts)
        df.insert(0, "year", y)
        out.append(df.reset_index())
    return pd.concat(out, ignore_index=True)


def stage_terraclimate():
    muni, lon, lat = municipality_raster()
    ds = open_terraclimate()

    tc_lon = ds.lon.values
    tc_lat = ds.lat.values
    lon_sl = slice(int(np.searchsorted(tc_lon, LON0)) - 1,
                   int(np.searchsorted(tc_lon, LON1)) + 2)
    lat_sl = slice(int(np.searchsorted(-tc_lat, -LAT1)) - 1,
                   int(np.searchsorted(-tc_lat, -LAT0)) + 2)
    lon_sl = slice(max(lon_sl.start, 0), min(lon_sl.stop, len(tc_lon)))
    lat_sl = slice(max(lat_sl.start, 0), min(lat_sl.stop, len(tc_lat)))
    sub_lon = tc_lon[lon_sl]
    sub_lat = tc_lat[lat_sl]
    say(f"brazil window: {len(sub_lat)} lat x {len(sub_lon)} lon TerraClimate cells")

    flat, codes = build_index(sub_lon, sub_lat, muni, lon, lat)
    uniq, code_index = np.unique(codes, return_inverse=True)
    say(f"{len(uniq):,} municipalities with at least one cell")

    years = list(range(PULL_MIN, YEAR_MAX + 1))
    todo = [y for y in years
            if not os.path.exists(os.path.join(OUTDIR, f"{y}.parquet"))]
    say(f"{len(todo)} of {len(years)} years to pull")

    t0 = time.time()
    failed = []
    for i, y in enumerate(todo, 1):
        try:
            df = pull_year(ds, y, flat, code_index, uniq, lat_sl, lon_sl)
            df.to_parquet(os.path.join(OUTDIR, f"{y}.parquet"), index=False)
            el = (time.time() - t0) / 60
            say(f"  {y} done ({i}/{len(todo)}, {el:.1f} min, "
                f"eta {el / i * (len(todo) - i):.0f} min)")
        except Exception as e:
            failed.append(y)
            say(f"  {y} FAILED after retries: {str(e)[:200]}")
            if len(failed) >= 3:
                sys.exit(f"aborting: {len(failed)} years failed even with backoff. "
                         f"Completed years are checkpointed in {OUTDIR}/.")

    files = sorted(glob.glob(os.path.join(OUTDIR, "*.parquet")))
    have = {int(os.path.basename(f)[:-8]) for f in files}
    missing = sorted(set(years) - have)
    if missing:
        sys.exit(f"REFUSING TO ASSEMBLE: years {missing} are missing. Re-run to "
                 f"fill them in; a partial climate panel would be silently wrong.")

    say(f"assembling {len(files)} years")
    mon = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    say(f"municipality-months: {len(mon):,}")
    assert mon.groupby(["year", "month"]).size().nunique() == 1, \
        "municipality count varies by month"

    yr = aggregate(mon).sort_values(["code_muni_7", "year"]).reset_index(drop=True)

    # ---- validation ---------------------------------------------------------
    n_exp = len(uniq) * (YEAR_MAX - YEAR_MIN + 1)
    assert len(yr) == n_exp, f"expected {n_exp} municipality-years, got {len(yr)}"

    # Physical sanity. These are ranges, not point checks: the aim is to catch a
    # unit error or a scrambled index, not to police the climatology.
    checks = [
        ("ppt_cal", 200, 6000, "annual precipitation, mm"),
        ("tmax_cal", 15, 40, "mean of monthly max temperature, C"),
        ("pdsi_cal", -8, 8, "Palmer index"),
        ("vpd_cal", 0.1, 4.0, "TerraClimate VPD, kPa"),
        ("def_cal", 0, 3000, "climatic water deficit, mm"),
    ]
    for col, lo, hi, lab in checks:
        m = yr[col].median()
        assert lo <= m <= hi, f"{col} median {m:.2f} outside [{lo}, {hi}] - {lab}"
        say(f"  {col:12s} median {m:8.2f}  [{lo}, {hi}]  {lab}")

    # PET must exceed AET everywhere, by definition of the deficit.
    bad = (yr["aet_cal"] > yr["pet_cal"] + 1e-6).sum()
    assert bad == 0, f"{bad} municipality-years have AET > PET, which is impossible"

    # def == pet - aet, to rounding. If this fails the variables are misaligned.
    resid = (yr["def_cal"] - (yr["pet_cal"] - yr["aet_cal"])).abs()
    say(f"  def vs pet-aet: max |diff| {resid.max():.3f} mm")
    assert resid.max() < 5.0, "def does not reconcile with pet - aet"

    # The growing-season window straddles years, so 2002 must be complete only
    # because 2001 was pulled. Verify it is not silently NaN.
    na_grow = yr.loc[yr.year == YEAR_MIN, "ppt_grow"].isna().mean()
    assert na_grow < 0.02, \
        f"{na_grow:.1%} of {YEAR_MIN} growing-season values are NaN - the " \
        f"{PULL_MIN} back-year did not merge"

    yr.to_csv(OUT, index=False)
    say(f"wrote {OUT}: {len(yr):,} municipality-years x {yr.shape[1]} columns")
    say(f"  municipalities {yr.code_muni_7.nunique():,}, "
        f"years {yr.year.min()}-{yr.year.max()}")


def stage_degreedays():
    """Delegate to degree_days.py, which owns its own checkpoints."""
    import degree_days
    degree_days.main()


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", default="all",
                    choices=["terraclimate", "degreedays", "all"])
    args = ap.parse_args()
    if args.stage in ("terraclimate", "all"):
        say("=== stage: terraclimate ===")
        stage_terraclimate()
    if args.stage in ("degreedays", "all"):
        say("=== stage: degreedays ===")
        stage_degreedays()


if __name__ == "__main__":
    main()
