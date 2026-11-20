"""Per-cell land cover and elevation from the DEM, the land-cover layers and the water layer (water wins).
Classes: water; upland = UPLAND layers; marsh = every other layer; unclassified = none. Writes MODEL/cell_stats.csv."""
import glob, os
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from config import CELLS, DEM, WATER, LANDCOVER, MODEL, UPLAND

os.makedirs(MODEL, exist_ok=True)
src = rasterio.open(DEM); dem = src.read(1).astype("float64"); valid = dem > -50
cells = gpd.read_file(CELLS); c4 = cells.to_crs(src.crs); NC = len(cells)
ras = lambda shapes, dtype: rasterize(shapes, out_shape=src.shape, transform=src.transform, fill=0, dtype=dtype)
cellr = ras(((g, i + 1) for i, g in enumerate(c4.geometry)), "uint16")
layers = [f for f in sorted(glob.glob(f"{LANDCOVER}/*.shp")) if os.path.basename(f) != "Water.shp"]
names = ["none"] + [os.path.basename(f)[:-4] for f in layers] + ["Water"]
cover = np.zeros(src.shape, dtype="uint8")
for k, f in enumerate(layers, 1):
    g = gpd.read_file(f).to_crs(src.crs); cover[ras(((gg, k) for gg in g.geometry if gg is not None), "uint8") > 0] = k
wu = gpd.read_file(WATER).to_crs(src.crs)
cover[ras([(wu.geometry.union_all(), 1)], "uint8") > 0] = len(names) - 1
m = (cellr > 0) & valid
ci = cellr[m].astype(np.int64) - 1; cl = cover[m].astype(np.int64); z = dem[m]; NK = len(names); key = ci * NK + cl
cnt = np.bincount(key, minlength=NC * NK).reshape(NC, NK); sz = np.bincount(key, weights=z, minlength=NC * NK).reshape(NC, NK)
szz = np.bincount(key, weights=z * z, minlength=NC * NK).reshape(NC, NK)
upl = [names.index(n) for n in UPLAND if n in names]; wat = [len(names) - 1]; non = [0]
mar = [k for k in range(NK) if k not in upl + wat + non]


def grp(idx):
    n = cnt[:, idx].sum(1); s = sz[:, idx].sum(1); ss = szz[:, idx].sum(1)
    mu = np.where(n > 0, s / np.maximum(n, 1), np.nan); sd = np.sqrt(np.maximum(ss / np.maximum(n, 1) - mu ** 2, 0)); sd[n == 0] = np.nan
    return n, mu, sd


tot = cnt.sum(1); nw, zw, _ = grp(wat); nm, zm, sm = grp(mar); nu, zu, _ = grp(upl); n0, _, _ = grp(non)
out = pd.DataFrame({"cell": range(NC), "area_ha": (cells.area / 1e4).round(1), "f_water": (nw / tot).round(3), "f_marsh": (nm / tot).round(3),
                    "f_upland": (nu / tot).round(3), "f_unclass": (n0 / tot).round(3), "bed_z_water": zw.round(3), "marsh_z_mean": zm.round(3),
                    "marsh_z_sd": sm.round(3), "upland_z": zu.round(3),
                    "dominant_cover": [names[mar[np.argmax(cnt[i, mar])]] if nm[i] > 0 else "-" for i in range(NC)]})
out.to_csv(f"{MODEL}/cell_stats.csv", index=False)
print(f"domain: water {nw.sum() / tot.sum():.3f}, marsh {nm.sum() / tot.sum():.3f}, upland {nu.sum() / tot.sum():.3f}")
print("cells without open water:", list(out.cell[out.f_water <= 0.005]), "->", f"{MODEL}/cell_stats.csv")
