#!/usr/bin/env python3
"""
24_patch_metrics_500m.py - burn-scar size distribution at NATIVE MCD64A1 resolution.

RECONSTRUCTED 13 October 2026 from the session transcript after the working
directory was lost. Never re-executed since. Verify before trusting.

Replaces an earlier patch_metrics.py, which worked from the legacy 2 km
EPSG:4326 grids. That grid could not represent any scar smaller than 3.4 km2,
and a one-month test at native resolution found that 92% of scars fall below
that floor. The 2 km version was blind to nine fires out of ten. The full
native-resolution run confirmed it: 3,138,369 scars, median 0.429 km2, and
93.6% below the old floor.

WHAT CHANGED AND WHY IT MATTERS

  1. RESOLUTION. 463.3 m sinusoidal cells, 0.2147 km2 each. Sinusoidal is
     equal-area, so unlike the 2 km geographic grid there is no latitude
     correction to get wrong.
  2. BURN DATE. 02_acquire_fire_mcd64.py reads the Burn_Date band and
     immediately does `> 0`, discarding the day-of-year. Keeping it lets scars
     be separated in TIME as well as space. On the test month that took the
     scar count from 15,908 to 43,668: those are unrelated fires burning the
     same ground weeks apart that a monthly composite had glued into one blob.

  Neither is a new download. The pipeline already pulls this exact band and
  throws both properties away in the reproject.

METHOD
  For each month, for each MODIS tile: read Burn_Date, split the month into
  overlapping DAY_GAP-day windows, connected-component label within each
  window, and record every scar's area and centroid. Centroids are converted
  from sinusoidal to lon/lat and looked up in a municipality raster built once
  at ~500 m.

KNOWN LIMITATION, carried over and still unfixed: scars are labelled within
each MODIS tile, so a fire crossing a tile boundary is counted twice. Tiles are
~1200 km wide so this affects few scars, but it inflates counts slightly. A
mosaic pass would fix it; it is not free because the tiles are in different
sinusoidal blocks.

Checkpointed one file per month under patch500/, so an interrupted run resumes.

Output: patch_metrics_500m_muni_quarter.csv
Run from code/.
"""

import glob
import os
import sys
import time

import numpy as np
import pandas as pd
import geopandas as gpd
import planetary_computer
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.features import rasterize
from rasterio.transform import from_origin
from scipy import ndimage

OUTDIR = "patch500"
OUT = "patch_metrics_500m_muni_quarter.csv"
MESH = "malha_inter.json"
CELL_KM2 = (463.3127165279167 / 1000.0) ** 2
DAY_GAP = 4           # scars more than this many days apart are distinct fires
LARGE_SMALL, LARGE_BIG = 1.0, 10.0   # km2 thresholds, meaningful at 0.21 km2 cells

# municipality lookup raster, ~500 m in EPSG:4326 over the Brazil bbox
LON0, LON1, LAT0, LAT1 = -74.05, -33.95, -34.10, 5.90
PX = 0.004491          # degrees, about 500 m at the equator

os.makedirs(OUTDIR, exist_ok=True)


def say(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def retry(fn, what, tries=6, base=5):
    """Exponential backoff around a network call.

    The first version of this script had a bare try/except around each month
    and no retry. When the egress proxy dropped for a few minutes it burned
    through 220 consecutive months in under a second, logging each as FAILED,
    and then cheerfully tried to assemble a panel from the 56 that had already
    succeeded. Transient network failure must cost time, not data.
    """
    for i in range(tries):
        try:
            return fn()
        except Exception as e:
            if i == tries - 1:
                raise
            wait = base * (2 ** i)
            say(f"    {what}: {type(e).__name__}, retry {i + 1}/{tries - 1} "
                f"in {wait}s")
            time.sleep(wait)


def municipality_lookup():
    """Rasterise the fixed 2022 IBGE mesh once at ~500 m. Returns the array and
    the affine transform needed to turn lon/lat into row/col."""
    w = int(round((LON1 - LON0) / PX))
    h = int(round((LAT1 - LAT0) / PX))
    tr = from_origin(LON0, LAT1, PX, PX)
    mesh = gpd.read_file(MESH)
    code_col = next(c for c in mesh.columns
                    if c.lower() in ("codarea", "cd_mun", "code_muni_7"))
    mesh[code_col] = mesh[code_col].astype(np.int64)
    if mesh.crs is None:
        mesh = mesh.set_crs("EPSG:4326")
    say(f"rasterising {len(mesh)} municipalities at {h}x{w} (~500 m)")
    arr = rasterize(((g, c) for g, c in zip(mesh.geometry, mesh[code_col])),
                    out_shape=(h, w), transform=tr, fill=0, dtype="int32")
    covered = int((arr > 0).sum())
    say(f"  {covered:,} cells assigned ({100 * covered / arr.size:.1f}% of bbox)")
    assert covered > 10_000_000, "municipality raster looks empty"
    return arr, tr


def tile_scars(bd, src_crs, transform, to_wgs):
    """One tile's Burn_Date band -> (areas_km2, lon, lat) per scar."""
    pos = bd > 0
    if not pos.any():
        return None
    days = bd[pos]
    lo, hi = int(days.min()), int(days.max())
    areas, rows, cols = [], [], []
    for a in range(lo, hi + 1, DAY_GAP + 1):
        m = pos & (bd >= a) & (bd <= a + DAY_GAP)
        if not m.any():
            continue
        lab, n = ndimage.label(m, structure=np.ones((3, 3), int))
        if not n:
            continue
        idx = np.arange(1, n + 1)
        cnt = np.bincount(lab.ravel(), minlength=n + 1)[1:]
        cen = np.asarray(ndimage.center_of_mass(m, lab, index=idx))
        areas.append(cnt * CELL_KM2)
        rows.append(cen[:, 0])
        cols.append(cen[:, 1])
    if not areas:
        return None
    areas = np.concatenate(areas)
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    # pixel centre -> sinusoidal metres -> lon/lat
    x = transform.c + (cols + 0.5) * transform.a
    y = transform.f + (rows + 0.5) * transform.e
    lon, lat = to_wgs.transform(x, y)
    return areas, lon, lat


def do_month(cat, year, mon, muni, tr):
    def _search():
        q = cat.search(collections=["modis-64A1-061"], bbox=[-74, -34, -32, 5.5],
                       datetime=f"{year}-{mon:02d}-01/{year}-{mon:02d}-01")
        return list(q.items())
    items = retry(_search, f"search {year}-{mon:02d}")
    if len(items) < 15:
        say(f"  WARNING {year}-{mon:02d}: only {len(items)} tiles")
    A, LO, LA = [], [], []
    to_wgs = None
    for it in items:
        def _read(it=it):
            with rasterio.open(it.assets["Burn_Date"].href) as src:
                return src.read(1), src.crs, src.transform
        bd, crs, tf = retry(_read, f"tile {year}-{mon:02d}")
        if to_wgs is None:
            to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        r = tile_scars(bd, crs, tf, to_wgs)
        if r is None:
            continue
        A.append(r[0]); LO.append(r[1]); LA.append(r[2])
    if not A:
        return pd.DataFrame(columns=["code_muni_7", "area"])
    area = np.concatenate(A); lon = np.concatenate(LO); lat = np.concatenate(LA)
    col = ((lon - LON0) / PX).astype(int)
    row = ((LAT1 - lat) / PX).astype(int)
    ok = (row >= 0) & (row < muni.shape[0]) & (col >= 0) & (col < muni.shape[1])
    code = np.zeros(area.shape, dtype=np.int64)
    code[ok] = muni[row[ok], col[ok]]
    keep = code > 0
    return pd.DataFrame({"code_muni_7": code[keep], "area": area[keep]})


def aggregate(df):
    g = df.groupby(["code_muni_7", "year", "quarter"], sort=False)
    out = g.agg(n_patches=("area", "size"),
                patch_area_km2=("area", "sum"),
                max_patch_km2=("area", "max"),
                p90_patch_km2=("area", lambda s: float(np.percentile(s, 90))),
                median_patch_km2=("area", "median")).reset_index()
    ssq = g.apply(lambda d: float((d["area"] ** 2).sum()), include_groups=False)
    out["mean_patch_km2"] = ssq.values / out["patch_area_km2"].values
    for thr, name in ((LARGE_SMALL, "share_large_1"), (LARGE_BIG, "share_large_10")):
        big = df[df["area"] > thr].groupby(
            ["code_muni_7", "year", "quarter"], sort=False)["area"].sum()
        out = out.merge(big.rename(name).reset_index(),
                        on=["code_muni_7", "year", "quarter"], how="left")
        out[name] = (out[name].fillna(0.0) / out["patch_area_km2"]).astype(float)
    return out


def main():
    muni, tr = municipality_lookup()
    cat = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace)

    months = [(y, m) for y in range(2002, 2025) for m in range(1, 13)]
    todo = [(y, m) for (y, m) in months
            if not os.path.exists(os.path.join(OUTDIR, f"{y}_{m:02d}.parquet"))]
    say(f"{len(todo)} of {len(months)} months to do")

    t0 = time.time()
    failed = []
    for i, (y, m) in enumerate(todo, 1):
        try:
            d = do_month(cat, y, m, muni, tr)
            d.to_parquet(os.path.join(OUTDIR, f"{y}_{m:02d}.parquet"), index=False)
            if i % 12 == 0:
                el = (time.time() - t0) / 60
                say(f"  {i}/{len(todo)} months, {el:.1f} min, "
                    f"eta {el / i * (len(todo) - i):.0f} min")
        except Exception as e:
            failed.append((y, m))
            say(f"  {y}-{m:02d} FAILED after retries: {str(e)[:160]}")
            if len(failed) >= 5:
                sys.exit(f"aborting: {len(failed)} months failed even with "
                         f"backoff - the network or the API is down. Fix it and "
                         f"re-run; completed months are checkpointed in {OUTDIR}/.")

    files = sorted(glob.glob(os.path.join(OUTDIR, "*.parquet")))
    if len(files) != len(months):
        missing = sorted(set(months) - {
            tuple(int(x) for x in os.path.basename(f)[:-8].split("_")) for f in files})
        sys.exit(f"REFUSING TO ASSEMBLE: {len(missing)} of {len(months)} months "
                 f"are missing, first few {missing[:5]}. Re-run to fill them in; "
                 f"a panel built from a partial set would be silently wrong.")
    say(f"assembling {len(files)} months")
    parts = []
    for f in files:
        y, m = (int(x) for x in os.path.basename(f)[:-8].split("_"))
        d = pd.read_parquet(f)
        if not len(d):
            continue
        d["year"] = y
        d["quarter"] = (m - 1) // 3 + 1
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    say(f"total burn scars: {len(df):,}")
    say(f"  median {df.area.median():.3f} km2  p90 {df.area.quantile(.9):.2f}  "
        f"max {df.area.max():,.0f}")
    say(f"  under 3.4 km2 (the old 2 km floor): {100 * (df.area < 3.4).mean():.1f}%")

    out = aggregate(df).sort_values(
        ["code_muni_7", "year", "quarter"]).reset_index(drop=True)
    assert out["n_patches"].min() >= 1
    assert (out["patch_area_km2"] > 0).all()
    assert (out["mean_patch_km2"] <= out["max_patch_km2"] + 1e-6).all()
    for c in ("share_large_1", "share_large_10"):
        assert out[c].between(-1e-9, 1 + 1e-9).all(), c
        assert out[c].std() > 0.01, f"{c} is degenerate"
    out.to_csv(OUT, index=False)
    say(f"wrote {OUT}: {len(out):,} municipality-quarters, "
        f"{out.code_muni_7.nunique():,} municipalities")
    say(out[["n_patches", "median_patch_km2", "mean_patch_km2",
             "share_large_1", "share_large_10"]].describe().round(3).to_string())


if __name__ == "__main__":
    main()
