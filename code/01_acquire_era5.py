"""
01_acquire_era5.py
==================
Every ERA5 download, in one script. Replaces three separate drips that differed
only in their request payload.

WHAT IS NEEDED, AND WHY EACH PIECE EXISTS
  2021-2024   daily mean t2m, daily max t2m, daily mean d2m. Extends the
              original series through 2024.
  2002-2020   daily MAX t2m only. The original project archived daily-mean
              grids but never the daily maxima, and `vpd_max_q` cannot be built
              without them, so they are re-pulled.
  2019-09     daily mean d2m for one month. Absent from the archive entirely -
              the folder holds eleven dewpoint files for 2019, not twelve.
  2024        daily mean t2m as twelve MONTHLY requests. The year-level request
              for 2024 failed three times with a bare {"status":"failed"} and an
              empty traceback; a monthly request for the same data succeeded.
              A transient CDS backend fault, not a malformed request, but the
              fallback is kept because it will be needed again if it recurs.

HOW THE QUEUE IS HANDLED. Copernicus allows two queued requests per user, and a
year-level request can sit queued for hours - during this project the queue
slowed from about one file an hour to one every seven. The script therefore
submits asynchronously, records request ids in a state file, polls them on the
next invocation, downloads whatever has completed, and tops the queue back up.
It is meant to be run repeatedly (from cron, or by hand) until it prints
ALL DONE. Nothing is lost if it is killed between runs.

Year-level requests are used wherever possible: they are the largest granularity
CDS allows and roughly twelve times more transfer-efficient than monthly ones.
Since the two-in-flight cap is the real constraint, monthly requests are not
faster - they are only a fallback for a year that will not build.

CREDENTIALS. Set CDS_API_KEY in the environment, or put it in ~/.cdsapirc.
Do not commit a key to this file.

Outputs: era5/*.nc  (one file per target)
Seed: 123
"""

import json
import logging
import os
import sys
import time

import numpy as np

np.random.seed(123)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                    handlers=[logging.FileHandler("01_acquire_era5.log"),
                              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

from ecmwf.datastores import Client   # noqa: E402

DATASET = "derived-era5-single-levels-daily-statistics"
AREA = [5.5, -74, -34, -32]          # Brazil bounding box, 0.25 deg -> 159 x 169
OUT = "era5"
STATE = "era5_requests.json"
MAX_INFLIGHT = 2                     # CDS hard cap per user
ALL_MONTHS = [f"{m:02d}" for m in range(1, 13)]
ALL_DAYS = [f"{d:02d}" for d in range(1, 32)]
os.makedirs(OUT, exist_ok=True)

API_KEY = os.environ.get("CDS_API_KEY")
if not API_KEY:
    sys.exit("set CDS_API_KEY in the environment (or configure ~/.cdsapirc)")


def req(variable, statistic, year, months):
    return {"product_type": "reanalysis", "variable": [variable],
            "year": str(year), "month": months, "day": ALL_DAYS,
            "daily_statistic": statistic, "time_zone": "utc+00:00",
            "frequency": "1_hourly", "area": AREA}


def good(path, min_bytes=5_000_000):
    if not os.path.exists(path) or os.path.getsize(path) < min_bytes:
        return False
    try:
        import xarray as xr
        with xr.open_dataset(path) as ds:
            return len(ds.data_vars) > 0
    except Exception:
        return False


# ---- the full target list ---------------------------------------------------
targets = []   # (key, output path, request dict, minimum plausible size)

for year in range(2021, 2025):
    for var, stat, tag in (("2m_temperature", "daily_mean", "t2m_mean"),
                           ("2m_temperature", "daily_maximum", "t2m_max"),
                           ("2m_dewpoint_temperature", "daily_mean", "d2m_mean")):
        targets.append((f"{tag}_{year}",
                        os.path.join(OUT, f"{tag}_{year}_FULL.nc"),
                        req(var, stat, year, ALL_MONTHS), 5_000_000))

for year in range(2002, 2021):
    targets.append((f"t2m_max_{year}",
                    os.path.join(OUT, f"t2m_max_{year}_FULL.nc"),
                    req("2m_temperature", "daily_maximum", year, ALL_MONTHS),
                    5_000_000))

targets.append(("d2m_mean_2019_09",
                os.path.join(OUT, "d2m_mean_2019_09.nc"),
                req("2m_dewpoint_temperature", "daily_mean", 2019, ["09"]),
                100_000))

# monthly fallback for 2024 t2m mean; only pursued if the year request has failed
if os.environ.get("ERA5_2024_MONTHLY_FALLBACK") == "1":
    for m in range(1, 13):
        targets.append((f"t2m_mean_2024_{m:02d}",
                        os.path.join(OUT, f"t2m_mean_2024_{m:02d}.nc"),
                        req("2m_temperature", "daily_mean", 2024, [f"{m:02d}"]),
                        100_000))

client = Client(url="https://cds.climate.copernicus.eu/api", key=API_KEY)
state = json.load(open(STATE)) if os.path.exists(STATE) else {}

# ---- poll what is already in flight -----------------------------------------
inflight = 0
for key, path, _, minb in targets:
    rid = state.get(key)
    if not rid:
        continue
    try:
        remote = client.get_remote(rid)
        status = remote.status
    except Exception as e:
        log.error("%s: poll failed (%s) - clearing so it resubmits", key, str(e)[:80])
        state.pop(key, None)
        continue
    if status == "successful":
        remote.download(path)
        log.info("%s DOWNLOADED (%.1f MB)", key, os.path.getsize(path) / 1e6)
        state.pop(key, None)
    elif status in ("failed", "rejected"):
        log.info("%s %s - will resubmit on a later run", key, status)
        state.pop(key, None)
    else:
        log.info("%s in flight: %s", key, status)
        inflight += 1

# ---- top the queue back up ---------------------------------------------------
for key, path, request, minb in targets:
    if inflight >= MAX_INFLIGHT:
        break
    if key in state or good(path, minb):
        continue
    try:
        remote = client.submit(DATASET, request)
        state[key] = remote.request_id
        log.info("submitted %s -> %s", key, remote.request_id)
        inflight += 1
        time.sleep(5)
    except Exception as e:
        log.error("submit %s failed: %s", key, str(e)[:120])

json.dump(state, open(STATE, "w"))
remaining = [k for k, p, _, mb in targets if not good(p, mb)]
if remaining:
    log.info("remaining: %d (%s%s)", len(remaining), ", ".join(remaining[:8]),
             " ..." if len(remaining) > 8 else "")
    log.info("re-run this script until it prints ALL DONE")
else:
    log.info("ALL DONE")
    print("ALL DONE")
