# Environment

Pinned from the container the 2026 rebuild was run in. Later versions will very
likely work; these are the ones the results were produced under.

## Python 3.11.15

| package | version |
|---|---|
| pandas | 3.0.2 |
| numpy | 2.4.4 |
| geopandas | 1.1.2 |
| shapely | 2.1.0 |
| xarray | 2026.7.0 |
| rasterio | 1.4.4 |
| exactextract | 0.3.0 |
| pyreadr | 0.5.6 |
| ecmwf-datastores-client | (CDS API client, stages 23-25) |

```bash
pip install pandas numpy geopandas shapely xarray rasterio exactextract \
            pyreadr ecmwf-datastores-client --break-system-packages
```

## R 4.3.3

| package | version |
|---|---|
| data.table | 1.18.6.1 |
| fixest | 0.14.2 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |
| R.utils | 2.13.0 |
| sandwich | 3.1.3 |
| lmtest | 0.9.40 |

```r
install.packages(c("data.table","fixest","ggplot2","patchwork",
                   "R.utils","sandwich","lmtest"))
```

`R.utils` is not optional — `data.table::fread` needs it to read the gzipped
panel directly.

## Credentials and network

- **Copernicus CDS** (stages 23-25): an API key is required. The key is
  currently hard-coded in the download scripts; move it to `~/.cdsapirc` or an
  environment variable before this leaves the working directory.
- **IBGE SIDRA** (stages 31-33): public, no key. Rate-limited in practice to
  roughly one municipality-level request every 25-40 seconds.
- **Microsoft Planetary Computer** (stage 26a): public STAC endpoint for
  MCD64A1.

## Seed

`123` throughout, set at the top of every script that has any stochastic step.
No stage in the current pipeline is actually stochastic; the seed is set so that
adding a bootstrap or a sampling step later does not silently change results.
