#!/usr/bin/env python3
"""
degree_days.py - Schlenker-Roberts growing and extreme degree days.

RECONSTRUCTED 13 October 2026 from the session transcript after the working
directory was lost. Never re-executed since. Verify before trusting.

Helper to 23_climate.py, in the same relation 24_patch_metrics_500m.py has to
06_build_analysis_panel.py: standalone because it is slow and checkpointed,
called by the parent for the degreedays stage.

WHY. The weather placebo needs at least one variable whose expected sign is not
in dispute. Schlenker & Roberts (2009, PNAS 106:15594) is that benchmark: US
crop yields rise with temperature up to roughly 29C for maize and 30C for soy
and fall sharply above it, and the decline is the single most robust finding in
the climate-agriculture literature. If extreme heat does not hurt Brazilian crop
outcomes in this panel, the design cannot detect agricultural effects and the
VPD null that the exclusion restriction rests on is uninformative. That is the
whole point of computing this.

NO NEW DOWNLOAD. Daily 2m maximum temperature is already on disk in era5_hist/
and era5_2021_2024/, and daily mean in cds_data_hist/. Both cover 2002-2024 on
the same 0.25 degree grid the VPD instrument is built from, so degree days and
VPD are measured on identical footing - which matters, because the comparison
between them is the test.

DAILY MINIMUM IS NOT ON DISK, and is reconstructed as Tmin = 2*Tmean - Tmax.
That is exact only if the diurnal cycle is symmetric about its mean. It is not,
quite: real daily temperature spends longer near the minimum than the maximum,
so this OVERSTATES Tmin slightly and therefore understates the diurnal range.
The consequence for degree days is second-order (both GDD and EDD shift the same
direction) but it is a real approximation and it is the reason the EDD figures
here should be read as a calibration check, not as a quantity to report as an
estimate in its own right. Pulling true daily minimum from CDS would remove it.

METHOD. Single-sine degree days (Snyder 1985), the standard in this literature.
For threshold t, with M = (Tmin+Tmax)/2 and a = (Tmax-Tmin)/2:

    Tmin >= t                DD_above = M - t
    Tmax <= t                DD_above = 0
    otherwise                DD_above = [ (M-t)(pi/2 - th) + a*cos(th) ] / pi
                             where th = arcsin((t - M)/a)

Degree days in a band [t1, t2] are DD_above(t1) - DD_above(t2). Reported:

    gdd_10_29    growing degree days, the productive band
    edd_29       extreme degree days above 29C, the damaging tail
    edd_32       above 32C, a stricter tail for sensitivity

The formula was unit-tested when first written: tmin>=t returns M-t; tmax<=t
returns 0; a flat day at the threshold returns 0; a day symmetric about t
returns a/pi; and the band 10-29 for a 15-25C day returns exactly 10.

Output: degree_days_muni_year.csv, checkpointed one parquet per year under
degdays/. Run from code/.
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
import xarray as xr
from rasterio.features import rasterize
from rasterio.transform import from_origin

ROOT = ".."                      # daily ERA5 lives at the repository root
OUTDIR = "degdays"
OUT = "degree_days_muni_year.csv"
MESH = "malha_inter.json"

YEAR_MIN, YEAR_MAX = 2002, 2021
# The daily ERA5 archive on disk starts in 2002, so unlike TerraClimate there is
# no 2001 back-year to draw on. The calendar window is therefore complete from
# 2002, but every window that reaches back into October of the prior year is
# unavailable for 2002 and is left NaN rather than summed over nine months and
# passed off as a year. STRADDLE_MIN records that, and the validation at the end
# asserts it holds in both directions: 2002 must be missing, 2003 must be there.
PULL_MIN = YEAR_MIN
STRADDLE_MIN = YEAR_MIN + 1

BASE, UPPER, EXTREME = 10.0, 29.0, 32.0
KELVIN = 273.15

LON0, LON1, LAT0, LAT1 = -74.05, -33.95, -34.10, 5.90
PX = 0.01

WINDOWS = {
    "cal":  [(0, m) for m in range(1, 13)],
    "agyr": [(-1, m) for m in range(10, 13)] + [(0, m) for m in range(1, 10)],
    "grow": [(-1, m) for m in range(10, 13)] + [(0, m) for m in range(1, 4)],
    "harv": [(0, m) for m in (2, 3, 4)],
}

os.makedirs(OUTDIR, exist_ok=True)


def say(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def dd_above(tmin, tmax, t):
    """Single-sine degree days above threshold t. Vectorised over any shape."""
    M = (tmax + tmin) / 2.0
    a = (tmax - tmin) / 2.0
    out = np.zeros_like(M)

    hot = tmin >= t
    out[hot] = M[hot] - t

    mid = (~hot) & (tmax > t)
    if mid.any():
        Mm, am = M[mid], a[mid]
        # a can be ~0 when max and reconstructed min coincide; clamp so arcsin
        # is defined. Those cells are caught by the `hot` or zero branches
        # anyway, so the clamp never drives a reported value.
        am = np.maximum(am, 1e-6)
        r = np.clip((t - Mm) / am, -1.0, 1.0)
        th = np.arcsin(r)
        out[mid] = ((Mm - t) * (np.pi / 2.0 - th) + am * np.cos(th)) / np.pi

    return np.maximum(out, 0.0)


def municipality_raster():
    """Same ~1 km municipality raster 23_climate.py builds, for the same reason:
    area-weighted means rather than centroid sampling, so small municipalities
    are not dropped by a 0.25 degree grid."""
    w = int(round((LON1 - LON0) / PX))
    h = int(round((LAT1 - LAT0) / PX))
    tr = from_origin(LON0, LAT1, PX, PX)
    mesh = gpd.read_file(MESH)
    code_col = next(c for c in mesh.columns
                    if c.lower() in ("codarea", "cd_mun", "code_muni_7"))
    mesh[code_col] = mesh[code_col].astype(np.int64)
    if mesh.crs is None:
        mesh = mesh.set_crs("EPSG:4326")
    say(f"rasterising {len(mesh)} municipalities at {h}x{w}")
    arr = rasterize(((g, c) for g, c in zip(mesh.geometry, mesh[code_col])),
                    out_shape=(h, w), transform=tr, fill=0, dtype="int64")
    n_muni = len(np.unique(arr)) - 1
    assert n_muni > 5500, f"only {n_muni} municipalities rasterised"
    lon = LON0 + (np.arange(w) + 0.5) * PX
    lat = LAT1 - (np.arange(h) + 0.5) * PX
    return arr, lon, lat


def nearest_index(grid, pts, descending=False):
    g = -grid if descending else grid
    p = -pts if descending else pts
    i = np.clip(np.searchsorted(g, p), 1, len(g) - 1)
    return np.where(np.abs(g[i - 1] - p) <= np.abs(g[i] - p), i - 1, i)


def find_tmax(year):
    """Daily maximum: a yearly FULL file, or twelve monthly files.

    The archive is not uniform - 2002-2020 landed as yearly FULL files and 2021
    as twelve monthlies, because the two pulls were run at different times. Both
    shapes are accepted; a partial set of monthlies is refused rather than
    silently producing a short year.
    """
    for d in ("era5_hist", "era5_2021_2024"):
        p = os.path.join(ROOT, d, f"t2m_max_brazil_era5_{year}_FULL.nc")
        if os.path.exists(p):
            return [p]
    for d in ("era5_hist", "era5_2021_2024"):
        monthly = sorted(glob.glob(os.path.join(
            ROOT, d, f"t2m_max_brazil_era5_{year}_[0-9][0-9].nc")))
        if len(monthly) == 12:
            return monthly
        if monthly:
            raise FileNotFoundError(
                f"{year}: only {len(monthly)} of 12 monthly t2m_max files in "
                f"{d}. Refusing to build a partial year.")
    raise FileNotFoundError(f"no t2m_max file for {year}")


def find_tmean(year):
    """Daily mean: monthly files in cds_data_hist, or a yearly FULL file."""
    monthly = sorted(glob.glob(os.path.join(
        ROOT, "cds_data_hist", f"2m_temperature_brazil_era5_{year}_[0-9][0-9].nc")))
    if len(monthly) == 12:
        return monthly
    p = os.path.join(ROOT, "era5_2021_2024", f"t2m_mean_brazil_era5_{year}_FULL.nc")
    if os.path.exists(p):
        return [p]
    if monthly:
        raise FileNotFoundError(
            f"{year}: only {len(monthly)} of 12 monthly mean-temperature files "
            f"and no yearly file. Refusing to build a partial year.")
    raise FileNotFoundError(f"no t2m_mean source for {year}")


def load_year(year):
    """Daily Tmax and Tmean for one calendar year, in Celsius, aligned on time."""
    px = [xr.open_dataset(p) for p in find_tmax(year)]
    dx = xr.concat(px, dim="valid_time") if len(px) > 1 else px[0]
    parts = [xr.open_dataset(p) for p in find_tmean(year)]
    dm = xr.concat(parts, dim="valid_time") if len(parts) > 1 else parts[0]

    tmax = dx["t2m"].sortby("valid_time")
    tmean = dm["t2m"].sortby("valid_time")
    tmean = tmean.reindex(valid_time=tmax.valid_time)

    n = tmax.sizes["valid_time"]
    exp = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    assert n == exp, f"{year}: {n} days of Tmax, expected {exp}"
    miss = int(tmean.isnull().all(dim=("latitude", "longitude")).sum())
    assert miss == 0, f"{year}: {miss} days have no mean temperature"

    tx = tmax.values.astype("float32")
    tm = tmean.values.astype("float32")
    if np.nanmedian(tx) > 100:          # Kelvin
        tx = tx - KELVIN
        tm = tm - KELVIN
    assert 5 < np.nanmedian(tx) < 45, f"{year}: Tmax median {np.nanmedian(tx):.1f}C"

    tn = 2.0 * tm - tx                  # symmetric-diurnal reconstruction
    tn = np.minimum(tn, tx)             # cannot exceed the maximum
    months = pd.to_datetime(tmax.valid_time.values).month.values
    return tx, tn, months, tmax.longitude.values, tmax.latitude.values


def year_frame(year, muni, lon, lat):
    tx, tn, months, e_lon, e_lat = load_year(year)

    ci = nearest_index(e_lon, lon)
    ri = nearest_index(e_lat, lat, descending=True)
    ok = muni > 0
    rr = np.repeat(ri[:, None], len(lon), axis=1)[ok]
    cc = np.repeat(ci[None, :], len(lat), axis=0)[ok]
    flat = rr.astype(np.int64) * len(e_lon) + cc.astype(np.int64)
    uniq, code_index = np.unique(muni[ok], return_inverse=True)
    n = len(uniq)

    rows = []
    for m in range(1, 13):
        sel = months == m
        d10 = dd_above(tn[sel], tx[sel], BASE).sum(axis=0)
        d29 = dd_above(tn[sel], tx[sel], UPPER).sum(axis=0)
        d32 = dd_above(tn[sel], tx[sel], EXTREME).sum(axis=0)
        vals = {"gdd_10_29": (d10 - d29).ravel()[flat],
                "edd_29": d29.ravel()[flat],
                "edd_32": d32.ravel()[flat]}
        agg = {}
        for k, v in vals.items():
            good = np.isfinite(v)
            num = np.bincount(code_index[good], weights=v[good], minlength=n)
            den = np.bincount(code_index[good], minlength=n)
            agg[k] = np.where(den > 0, num / np.maximum(den, 1), np.nan)
        df = pd.DataFrame(agg)
        df.insert(0, "month", m)
        df.insert(0, "year", year)
        df.insert(0, "code_muni_7", uniq)
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def aggregate(mon):
    """Degree days are a flux - they accumulate - so every window sums."""
    cols = ["gdd_10_29", "edd_29", "edd_32"]
    out = []
    for y in range(YEAR_MIN, YEAR_MAX + 1):
        parts = {}
        for wname, spec in WINDOWS.items():
            keys = set((y + dy, m) for dy, m in spec)
            sel = mon[[(yy, mm) in keys for yy, mm in zip(mon.year, mon.month)]]
            g = sel.groupby("code_muni_7")
            cnt = g.size()
            full = cnt[cnt == len(keys)].index
            for v in cols:
                s = g[v].sum()
                parts[f"{v}_{wname}"] = s.where(s.index.isin(full))
        df = pd.DataFrame(parts)
        df.insert(0, "year", y)
        out.append(df.reset_index())
    return pd.concat(out, ignore_index=True)


def main():
    muni, lon, lat = municipality_raster()
    years = list(range(PULL_MIN, YEAR_MAX + 1))
    todo = [y for y in years
            if not os.path.exists(os.path.join(OUTDIR, f"{y}.parquet"))]
    say(f"{len(todo)} of {len(years)} years to build")

    t0 = time.time()
    for i, y in enumerate(todo, 1):
        df = year_frame(y, muni, lon, lat)
        df.to_parquet(os.path.join(OUTDIR, f"{y}.parquet"), index=False)
        el = (time.time() - t0) / 60
        say(f"  {y} done ({i}/{len(todo)}, {el:.1f} min, "
            f"eta {el / i * (len(todo) - i):.0f} min)")

    files = sorted(glob.glob(os.path.join(OUTDIR, "*.parquet")))
    have = {int(os.path.basename(f)[:-8]) for f in files}
    missing = sorted(set(years) - have)
    if missing:
        sys.exit(f"REFUSING TO ASSEMBLE: years {missing} missing.")

    mon = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    yr = aggregate(mon).sort_values(["code_muni_7", "year"]).reset_index(drop=True)

    # ---- validation ---------------------------------------------------------
    assert (yr.filter(like="gdd").min().min() >= 0), "negative growing degree days"
    assert (yr.filter(like="edd").min().min() >= 0), "negative extreme degree days"
    # EDD above 32 is a subset of EDD above 29 and must never exceed it.
    for w in WINDOWS:
        bad = (yr[f"edd_32_{w}"] > yr[f"edd_29_{w}"] + 1e-6).sum()
        assert bad == 0, f"{bad} rows have edd_32 > edd_29 in window {w}"
    m = yr["gdd_10_29_cal"].median()
    assert 3000 < m < 8000, f"annual GDD median {m:.0f} outside plausible range"
    say(f"  gdd_10_29_cal median {m:.0f}, "
        f"edd_29_cal median {yr['edd_29_cal'].median():.0f}, "
        f"edd_32_cal median {yr['edd_32_cal'].median():.0f}")

    # The straddling windows must be absent in the first year and present in the
    # second. Checking only one direction would let a silent off-by-one through.
    na_first = yr.loc[yr.year == YEAR_MIN, "gdd_10_29_grow"].isna().mean()
    assert na_first > 0.98, (
        f"{YEAR_MIN} growing-season GDD is {1 - na_first:.1%} populated, but the "
        f"daily archive has no Oct-Dec {YEAR_MIN - 1} - it should be entirely NaN")
    na_next = yr.loc[yr.year == STRADDLE_MIN, "gdd_10_29_grow"].isna().mean()
    assert na_next < 0.02, (
        f"{na_next:.1%} of {STRADDLE_MIN} growing-season GDD is NaN - the "
        f"{YEAR_MIN} back-months did not merge")
    say(f"  straddling windows: {YEAR_MIN} NaN as expected, "
        f"{STRADDLE_MIN} populated ({1 - na_next:.1%})")

    yr.to_csv(OUT, index=False)
    say(f"wrote {OUT}: {len(yr):,} municipality-years x {yr.shape[1]} columns")


if __name__ == "__main__":
    main()
