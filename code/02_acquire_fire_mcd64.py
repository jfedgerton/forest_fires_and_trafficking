"""Rebuild the full burned-area series 2002-2024 from MCD64A1 (Planetary Computer).

Same construction as the original pipeline, one consistent source:
  - monthly binary presence on the original 2km EPSG:4326 grid
    (2km cell = burned if ANY 500m MCD64A1 cell inside it burned that month)
  - checkpointed: one .npz per month in mcd64_monthly/, skip if present.
Run repeatedly until it prints ALL MONTHS DONE.
"""
import logging, os, sys
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
import pystac_client, planetary_computer

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                    handlers=[logging.FileHandler("mcd64_build_all.log"),
                              logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

OUT = "mcd64_monthly"
os.makedirs(OUT, exist_ok=True)

ref = rasterio.open("tifs/burned_presence_2019_08.tif")
dst_transform, dst_crs = ref.transform, ref.crs
H, W = ref.height, ref.width
ref.close()

cat = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace)

months = [(y, m) for y in range(2002, 2025) for m in range(1, 13)]
todo = [(y, m) for (y, m) in months
        if not os.path.exists(os.path.join(OUT, f"presence_{y}_{m:02d}.npz"))]
log.info("months remaining: %d", len(todo))

for (y, m) in todo:
    try:
        s = cat.search(collections=["modis-64A1-061"], bbox=[-74, -34, -32, 5.5],
                       datetime=f"{y}-{m:02d}-01/{y}-{m:02d}-01")
        items = list(s.items())
        if len(items) < 15:
            log.warning("%d-%02d: only %d tiles found", y, m, len(items))
        dest = np.zeros((H, W), dtype="float32")
        for it in items:
            with rasterio.open(it.assets["Burn_Date"].href) as src:
                burned = (src.read(1) > 0).astype("float32")
                tile_dest = np.zeros((H, W), dtype="float32")
                reproject(burned, tile_dest,
                          src_transform=src.transform, src_crs=src.crs,
                          dst_transform=dst_transform, dst_crs=dst_crs,
                          resampling=Resampling.max)
                dest = np.fmax(dest, tile_dest)
        np.savez_compressed(os.path.join(OUT, f"presence_{y}_{m:02d}.npz"),
                            presence=(dest > 0).astype("uint8"), n_tiles=len(items))
        log.info("%d-%02d done | tiles %d | burned cells %d", y, m, len(items),
                 int((dest > 0).sum()))
    except Exception as e:
        log.error("%d-%02d FAILED: %s", y, m, str(e)[:150])

remaining = [(y, m) for (y, m) in months
             if not os.path.exists(os.path.join(OUT, f"presence_{y}_{m:02d}.npz"))]
if not remaining:
    log.info("ALL MONTHS DONE")
else:
    log.info("still missing: %d", len(remaining))
