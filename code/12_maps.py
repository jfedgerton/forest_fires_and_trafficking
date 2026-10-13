"""
12_maps.py
==========
Municipality choropleths. Python rather than R because this is the only step
that needs a geospatial stack, and geopandas is already a dependency of the
panel-building scripts; adding `sf` would mean a GDAL/GEOS/PROJ build for one
figure.

FIGURES
  map1_overview        four panels on the fixed 2022 IBGE mesh:
                         (a) mean quarterly area burned
                         (b) mean daily-max VPD - the instrument
                         (c) TF workers rescued, total
                         (d) UE victims, total
  map2_enforcement     task-force operations, plotted separately and labelled
                       explicitly as ENFORCEMENT ACTIVITY rather than
                       prevalence - a map of where inspectors went is not a map
                       of where forced labour is, and an unlabelled version
                       invites exactly that misreading
  map3_periods         area burned and TF rescues by period (2002-09, 2010-16,
                       2017-24), which shows the frontier moving
  map4_first_stage     first-stage F statistic by state, i.e. where the
                       identifying variation actually comes from

WHY QUANTILE BINS. Brazilian municipal fire and rescue distributions are
extremely skewed - the top percentile of municipality-quarters carries orders of
magnitude more burned area than the median. A continuous colour ramp renders as
one dark blob over the arc of deforestation and white everywhere else, which
hides all the variation the models actually use. Every panel is therefore binned
by quantile, with the bin edges printed in the legend so nothing is concealed.

Inputs : final_panel_muni_quarter.csv.gz, malha_inter.json
Outputs: output/figures/map*.{pdf,png}
Seed   : 123
"""

import os
import sys

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap

np.random.seed(123)
mpl.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "figure.dpi": 110, "savefig.bbox": "tight", "pdf.fonttype": 42,
})

def say(*a):
    print(*a, flush=True)


PANEL = "final_panel_muni_quarter.csv.gz"
MESH = "malha_inter.json"
FIGDIR = os.path.join("output", "figures")
os.makedirs(FIGDIR, exist_ok=True)
N_BINS = 6


for f in (PANEL, MESH):
    if not os.path.exists(f):
        sys.exit(f"missing {f} - run the build scripts first")

say("loading panel and mesh")
d = pd.read_csv(PANEL, usecols=[
    "code_muni_7", "year", "quarter", "percent_area_burned", "vpd_max_q",
    "tf_count", "ue_count", "n_operations", "uf_code", "municipality_area_km2"])
mesh = gpd.read_file(MESH).set_crs("EPSG:4326", allow_override=True)
mesh["code_muni_7"] = mesh["codarea"].astype(int)
mesh = mesh[["code_muni_7", "geometry"]].copy()
assert len(mesh) == 5570, f"mesh has {len(mesh)} municipalities, expected 5570"

# The published IBGE mesh contains self-intersections that make GEOS refuse any
# dissolve ("side location conflict"). Repair once, here, rather than in each
# panel; make_valid preserves boundaries, unlike buffer(0).
_bad = (~mesh.geometry.is_valid).sum()
if _bad:
    mesh["geometry"] = mesh.geometry.make_valid()
    say(f"repaired {_bad} invalid municipality geometries")

# national outline, computed once and reused by every panel
OUTLINE = mesh.dissolve().boundary
STATES = mesh.assign(uf_code=mesh.code_muni_7 // 100000).dissolve(by="uf_code").reset_index()


def quantile_layer(ax, gdf, col, title, cmap, unit, zero_label=None):
    """Draw one choropleth with quantile bins. Municipalities at exactly zero get
    their own flat category rather than being folded into the lowest bin, because
    for the outcomes 'never any rescue' is a different statement from 'few'."""
    vals = gdf[col]
    nonzero = vals[vals > 0]
    if len(nonzero) < N_BINS:
        edges = np.unique(np.quantile(vals, np.linspace(0, 1, N_BINS + 1)))
    else:
        edges = np.unique(np.quantile(nonzero, np.linspace(0, 1, N_BINS + 1)))
    colors = plt.get_cmap(cmap)(np.linspace(0.12, 0.95, len(edges) - 1))
    cmap_d = ListedColormap(colors)
    norm = BoundaryNorm(edges, cmap_d.N)

    zero = gdf[vals <= 0]
    pos = gdf[vals > 0]
    if len(zero):
        zero.plot(ax=ax, color="#f0f0f0", edgecolor="none")
    pos.plot(ax=ax, column=col, cmap=cmap_d, norm=norm, edgecolor="none")
    OUTLINE.plot(ax=ax, color="grey", linewidth=0.25)

    ax.set_title(title)
    ax.set_axis_off()
    sm = plt.cm.ScalarMappable(cmap=cmap_d, norm=norm)
    cb = plt.colorbar(sm, ax=ax, fraction=0.030, pad=0.01,
                      boundaries=edges, ticks=edges)
    cb.ax.set_yticklabels([f"{e:,.3g}" for e in edges], fontsize=6.5)
    cb.set_label(unit, fontsize=7)
    if len(zero) and zero_label:
        ax.plot([], [], "s", color="#f0f0f0", markeredgecolor="grey",
                markersize=6, label=f"{zero_label} (n={len(zero):,})")
        ax.legend(loc="lower left", fontsize=6.5, frameon=False)


def save(fig, stub):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{stub}.{ext}"))
    plt.close(fig)
    say("  wrote", stub)


# ---------------------------------------------------------------- map 1
say("map 1: overview")
agg = d.groupby("code_muni_7", as_index=False).agg(
    burned=("percent_area_burned", "mean"),
    vpd=("vpd_max_q", "mean"),
    tf=("tf_count", "sum"),
    ue=("ue_count", "sum"),
    ops=("n_operations", "sum"))
g = mesh.merge(agg, on="code_muni_7", how="left")

fig, axes = plt.subplots(2, 2, figsize=(11, 10))
quantile_layer(axes[0, 0], g, "burned", "(a) Area burned",
               "YlOrRd", "% of municipality per quarter", "Never burned")
quantile_layer(axes[0, 1], g, "vpd", "(b) Instrument: mean daily-max VPD",
               "PuBu", "kPa")
quantile_layer(axes[1, 0], g, "tf", "(c) TF workers rescued, 2002-2024",
               "YlGnBu", "workers", "No rescues recorded")
quantile_layer(axes[1, 1], g, "ue", "(d) UE trafficking victims, 2002-2024",
               "BuPu", "victims", "No victims recorded")
fig.suptitle("Fire, the instrument, and both outcomes across Brazilian municipalities",
             fontsize=12, fontweight="bold", y=0.98)
fig.text(0.5, 0.015, "Quantile bins computed over municipalities with non-zero values; "
                     "bin edges shown on each scale. Fixed 2022 IBGE mesh, 5,570 municipalities.",
         ha="center", fontsize=7.5, style="italic")
save(fig, "map1_overview")

# ---------------------------------------------------------------- map 2
say("map 2: enforcement")
fig, ax = plt.subplots(figsize=(7, 7))
quantile_layer(ax, g, "ops", "Task-force operations, 2002-2024",
               "Greys", "operations", "Never inspected")
ax.set_title("ENFORCEMENT ACTIVITY, not prevalence", fontsize=10, color="#B2182B")
fig.text(0.5, 0.04,
         "This maps where labour inspectors went, not where forced labour is.\n"
         "Task forces are dispatched on tips and on prior cases, so coverage is\n"
         "selective; the blank areas are places nobody inspected, not places with\n"
         "no forced labour.",
         ha="center", fontsize=8, style="italic", color="#B2182B")
save(fig, "map2_enforcement")

# ---------------------------------------------------------------- map 3
say("map 3: by period")
periods = [(2002, 2009), (2010, 2016), (2017, 2024)]
fig, axes = plt.subplots(2, 3, figsize=(14, 9))
for j, (a, b) in enumerate(periods):
    sub = d[(d.year >= a) & (d.year <= b)].groupby("code_muni_7", as_index=False).agg(
        burned=("percent_area_burned", "mean"), tf=("tf_count", "sum"))
    gg = mesh.merge(sub, on="code_muni_7", how="left")
    quantile_layer(axes[0, j], gg, "burned", f"Area burned, {a}-{b}",
                   "YlOrRd", "% per quarter", "Never burned")
    quantile_layer(axes[1, j], gg, "tf", f"TF rescues, {a}-{b}",
                   "YlGnBu", "workers", "None")
fig.suptitle("Fire and rescues by period: the frontier moves",
             fontsize=12, fontweight="bold", y=0.97)
fig.text(0.5, 0.015, "Bins are computed within each panel, so colours are comparable "
                     "in rank but not in level across periods.",
         ha="center", fontsize=7.5, style="italic")
save(fig, "map3_periods")

# ---------------------------------------------------------------- map 4
say("map 4: first-stage strength by state")
rows = []
for uf, sub in d.groupby("uf_code"):
    sub = sub.dropna(subset=["percent_area_burned", "vpd_max_q"])
    if sub.code_muni_7.nunique() < 20:
        continue
    # within-state first stage, demeaned by municipality and by year-quarter
    x = sub.vpd_max_q.to_numpy()
    y = sub.percent_area_burned.to_numpy()
    for key in ("code_muni_7", "year"):
        gm = sub.groupby(key)
        x = x - gm.vpd_max_q.transform("mean").to_numpy()
        y = y - gm.percent_area_burned.transform("mean").to_numpy()
    if np.std(x) == 0:
        continue
    beta = np.sum(x * y) / np.sum(x * x)
    resid = y - beta * x
    n, k = len(x), sub.code_muni_7.nunique()
    se = np.sqrt(np.sum(resid ** 2) / max(n - k - 1, 1) / np.sum(x * x))
    rows.append({"uf_code": uf, "beta": beta, "F": (beta / se) ** 2 if se > 0 else np.nan,
                 "n_munis": sub.code_muni_7.nunique()})
fs = pd.DataFrame(rows)
fs.to_csv(os.path.join("output", "results", "first_stage_by_state.csv"), index=False)

states = STATES.merge(fs, on="uf_code", how="left")
fig, ax = plt.subplots(figsize=(7.5, 7.5))
quantile_layer(ax, states, "F", "First-stage F statistic by state",
               "Greens", "F on the excluded instrument")
fig.text(0.5, 0.04,
         "Within-state first stage, municipality and year demeaned. Shows where the\n"
         "identifying variation in the instrument actually comes from.",
         ha="center", fontsize=8, style="italic")
save(fig, "map4_first_stage_by_state")

say("\nfirst-stage F by state, weakest ten:")
say(fs.nsmallest(10, "F")[["uf_code", "beta", "F", "n_munis"]].round(3).to_string(index=False))
say("\ndone")
