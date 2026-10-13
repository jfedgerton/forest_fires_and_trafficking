"""
04_build_fire_and_vpd.py
=======================
Both zonal-aggregation steps in one script: the MCD64A1 fire measures and the
three VPD instruments. They are merged because they share the whole machinery -
the fixed 2022 IBGE mesh, the raster profile, and the exact_extract
area-weighted aggregation - which was previously duplicated line-for-line.

STAGE A - FIRE MEASURES, from the monthly presence grids written by 02.
    presence_q  = max over the quarter's three months -> area_burned_no_double_sqkm
    frequency_q = sum over the quarter's three months -> total_burned_area_sqkm
    percent_area_burned = no_double / municipality area * 100

STAGE B - VPD INSTRUMENTS, from the ERA5 daily grids written by 01.
    vpd_mean_q     quarterly mean of daily es(t2m_mean) - es(d2m_mean)
    vpd_max_q      quarterly mean of daily es(t2m_max)  - es(d2m_mean)
    vpd_excess3_q  quarterly mean of daily max(0, VPDmax - 3.0 kPa)

  Tetens/FAO-56, es(T) = 0.6108 * exp(17.27*T/(T+237.3)) in kPa, VPD clamped at 0.
  Daily VPD is computed first and THEN averaged - never es() of an averaged
  temperature. es() is convex in T, so averaging temperature first understates
  mean VPD by an amount that scales with Var(T), which varies systematically by
  biome and season. That would be a non-classical measurement error correlated
  with the treatment.

FIXED BOUNDARY VINTAGE. Both stages use the 2022 IBGE mesh for every year. A
year-switching vintage would make the unit of analysis a function of time and
break the fixed-effect panel.

Checkpointed per quarter. Re-run until it prints ALL DONE.
Outputs: fire_measures_muni_quarter.csv.gz, vpd_instruments_muni_quarter.csv
Seed: 123
"""

# =============================================================================
# STAGE A - fire measures
# =============================================================================
import logging, os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.io import MemoryFile
from exactextract import exact_extract

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                    handlers=[logging.FileHandler("build_fire_panel_v2.log"),
                              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

MON = "mcd64_monthly"
START_YEAR, END_YEAR = 2002, 2024
CELL_KM2 = 4.0

gj = gpd.read_file("malha_inter.json").set_crs("EPSG:4326", allow_override=True)
gj["code_muni_7"] = gj["codarea"].astype(int)
gj["municipality_area_km2"] = gj.to_crs(5880).geometry.area / 1e6
munis = gj[["code_muni_7", "municipality_area_km2", "geometry"]]
log.info("munis: %d", len(munis))

ref = rasterio.open("tifs/burned_presence_2019_08.tif")
profile = ref.profile.copy()
profile.update(count=1, dtype="float32", nodata=np.nan)
ref.close()

def zonal_sum(arr):
    with MemoryFile() as mem:
        with mem.open(**profile) as ds:
            ds.write(arr.astype("float32"), 1)
        with mem.open() as ds:
            return exact_extract(ds, munis, "sum", output="pandas")["sum"].to_numpy()

rows = []
for y in range(START_YEAR, END_YEAR + 1):
    for q in range(1, 5):
        arrs = []
        for m in range((q - 1) * 3 + 1, q * 3 + 1):
            arrs.append(np.load(os.path.join(MON, f"presence_{y}_{m:02d}.npz"))["presence"].astype("float32"))
        pres_q = np.max(np.stack(arrs), 0)
        freq_q = np.sum(np.stack(arrs), 0)
        no_double = np.nan_to_num(zonal_sum(pres_q)) * CELL_KM2
        total = np.nan_to_num(zonal_sum(freq_q)) * CELL_KM2
        rows.append(pd.DataFrame({
            "code_muni_7": munis.code_muni_7.values,
            "municipality_area_km2": munis.municipality_area_km2.values,
            "year": y, "quarter": q,
            "total_burned_area_sqkm": total,
            "area_burned_no_double_sqkm": no_double,
        }))
        log.info("%d Q%d | burned munis %d | total km2 %.0f", y, q,
                 int((no_double > 0).sum()), total.sum())

panel = pd.concat(rows, ignore_index=True)
panel["percent_area_burned"] = (panel.area_burned_no_double_sqkm /
                                panel.municipality_area_km2 * 100)
panel["code_muni_6"] = panel.code_muni_7 // 10
panel.to_csv("fire_measures_mcd64_muni_quarter_2002_2024.csv", index=False)
log.info("saved fire_measures_mcd64_muni_quarter_2002_2024.csv: %d rows", len(panel))
print("PANEL DONE")


# =============================================================================
# STAGE B - VPD instruments
# =============================================================================
import logging, os, sys
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from exactextract import exact_extract

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                    handlers=[logging.FileHandler("build_instruments.log"),
                              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

OUT_DIR = "vpd_quarters"
os.makedirs(OUT_DIR, exist_ok=True)
FINAL = "vpd_instruments_muni_quarter_2002_2024.csv"

H, W = 159, 169
TRANSFORM = from_origin(-74.125, 5.625, 0.25, 0.25)
profile = dict(driver="GTiff", height=H, width=W, count=1, dtype="float32",
               crs="EPSG:4326", transform=TRANSFORM, nodata=np.nan)

gj = gpd.read_file("malha_inter.json").set_crs("EPSG:4326", allow_override=True)
gj["code_muni_7"] = gj["codarea"].astype(int)
munis = gj[["code_muni_7", "geometry"]]
log.info("munis: %d", len(munis))

def es_kpa(t_c):
    return 0.6108 * np.exp(17.27 * t_c / (t_c + 237.3))

def load_month(var, y, m):
    """Return (days, H, W) float32 array in Kelvin for var in
    {'t2m_mean','d2m_mean','t2m_max'}, or None if missing."""
    if y <= 2020:
        if var == "t2m_mean":
            p = f"cds_data_hist/2m_temperature_brazil_era5_{y}_{m:02d}.nc"
        elif var == "d2m_mean":
            p = f"cds_data_hist/2m_dewpoint_temperature_brazil_era5_{y}_{m:02d}.nc"
        else:
            p = f"era5_hist/t2m_max_brazil_era5_{y}_FULL.nc"
    else:
        p = f"era5_2021_2024/{var}_brazil_era5_{y}_FULL.nc"
        if not os.path.exists(p):
            p = f"era5_2021_2024/{var}_brazil_era5_{y}_{m:02d}.nc"
    if not os.path.exists(p):
        return None
    ds = xr.open_dataset(p)
    tdim = [d for d in ds.dims if "time" in d][0]
    da = ds[list(ds.data_vars)[0]]
    if "FULL" in p:
        da = da.sel({tdim: ds[tdim].dt.month == m})
    assert da.shape[1] == H and da.shape[2] == W, f"grid mismatch {p} {da.shape}"
    arr = da.values.astype("float32")
    ds.close()
    return arr

def zonal_mean(arr2d):
    with MemoryFile() as mem:
        with mem.open(**profile) as ds:
            ds.write(arr2d.astype("float32"), 1)
        with mem.open() as ds:
            return exact_extract(ds, munis, "mean", output="pandas")["mean"].to_numpy()

quarters = [(y, q) for y in range(2002, 2025) for q in range(1, 5)]
for y, q in quarters:
    out = os.path.join(OUT_DIR, f"vpd_{y}_q{q}.csv")
    if os.path.exists(out):
        continue
    months = range((q - 1) * 3 + 1, q * 3 + 1)
    tmean, dmean, tmax = [], [], []
    ok = True
    for m in months:
        a, b, c = (load_month("t2m_mean", y, m), load_month("d2m_mean", y, m),
                   load_month("t2m_max", y, m))
        if a is None or b is None or c is None:
            miss = [n for n, v in zip(["t2m_mean", "d2m_mean", "t2m_max"],
                                      [a, b, c]) if v is None]
            log.info("%d Q%d: missing %s for month %02d - skipping quarter",
                     y, q, ",".join(miss), m)
            ok = False
            break
        if not (a.shape[0] == b.shape[0] == c.shape[0]):
            log.error("%d-%02d day-count mismatch %s %s %s", y, m,
                      a.shape, b.shape, c.shape)
            ok = False
            break
        tmean.append(a); dmean.append(b); tmax.append(c)
    if not ok:
        continue
    t_mean = np.concatenate(tmean) - 273.15
    d_mean = np.concatenate(dmean) - 273.15
    t_max = np.concatenate(tmax) - 273.15
    ea = es_kpa(d_mean)
    vpd_mean = np.clip(es_kpa(t_mean) - ea, 0, None)
    vpd_max = np.clip(es_kpa(t_max) - ea, 0, None)
    excess3 = np.clip(vpd_max - 3.0, 0, None)
    g1 = np.nanmean(vpd_mean, axis=0)
    g2 = np.nanmean(vpd_max, axis=0)
    g3 = np.nanmean(excess3, axis=0)
    df = pd.DataFrame({
        "code_muni_7": munis.code_muni_7.values,
        "year": y, "quarter": q,
        "vpd_mean_q": zonal_mean(g1),
        "vpd_max_q": zonal_mean(g2),
        "vpd_excess3_q": zonal_mean(g3),
        "n_days": vpd_mean.shape[0],
    })
    df.to_csv(out, index=False)
    log.info("%d Q%d done | days %d | mean vpd %.3f | mean vpdmax %.3f | excess3 %.4f",
             y, q, vpd_mean.shape[0], np.nanmean(df.vpd_mean_q),
             np.nanmean(df.vpd_max_q), np.nanmean(df.vpd_excess3_q))

done = [f for f in os.listdir(OUT_DIR) if f.endswith(".csv")]
log.info("quarters built: %d / 92", len(done))
if len(done) == 92:
    panel = pd.concat([pd.read_csv(os.path.join(OUT_DIR, f)) for f in sorted(done)],
                      ignore_index=True)
    panel["code_muni_6"] = panel.code_muni_7 // 10
    panel = panel.sort_values(["code_muni_7", "year", "quarter"])
    panel.to_csv(FINAL, index=False)
    log.info("saved %s: %d rows", FINAL, len(panel))
    print("ALL QUARTERS DONE")
