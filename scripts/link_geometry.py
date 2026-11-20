"""Candidate links between neighbouring cells from the cell polygons, the water layer and the DEM.
For every pair: shared edge, open-water crossing, edge elevations, lowest ground 20-50 m inside each cell, distance.
Writes MODEL/link_geometry_candidates.csv."""
import os, warnings
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.mask import mask
from pyproj import Transformer
from config import CELLS, DEM, WATER, MODEL, BAY_CELL

warnings.filterwarnings("ignore")
os.makedirs(MODEL, exist_ok=True)
TOL = 3.0   # m: cell edges need not be snapped
cells = gpd.read_file(CELLS); cells["cell"] = range(len(cells)); N = len(cells)
water = gpd.read_file(WATER).to_crs(cells.crs).geometry.buffer(0).union_all()
src = rasterio.open(DEM); tr = Transformer.from_crs(cells.crs, src.crs, always_xy=True)


def sample(geom, step=1.0):
    pts = []
    for g in getattr(geom, "geoms", [geom]):
        if g.length == 0: continue
        n = max(2, int(g.length / step)); pts += [g.interpolate(t, normalized=True) for t in np.linspace(0, 1, n)]
    if not pts: return np.array([])
    xs, ys = tr.transform([p.x for p in pts], [p.y for p in pts])
    v = np.array([x[0] for x in src.sample(zip(xs, ys))]); return v[v > -100]


def band_low(poly, L, w0=20.0, w1=50.0, k=10):
    """lowest ground 20-50 m inside the cell (10th lowest pixel): a real channel carries on into the cell"""
    band = L.buffer(w1).difference(L.buffer(w0)).intersection(poly)
    if band.is_empty: return np.nan
    g = gpd.GeoSeries([band], crs=cells.crs).to_crs(src.crs).iloc[0]
    try: arr, _ = mask(src, [g.__geo_interface__], crop=True, filled=False)
    except ValueError: return np.nan
    v = np.sort(np.asarray(arr.compressed(), float)); v = v[(v > -50) & (v < 50)]
    return round(float(v[min(k - 1, v.size - 1)]), 3) if v.size else np.nan


pct = lambda v, q: round(float(np.percentile(v, q)), 3) if v.size else np.nan
rp = cells.geometry.representative_point(); rows = []
for i in range(N):
    for j in range(i + 1, N):
        if cells.geometry[i].distance(cells.geometry[j]) > TOL: continue
        L = cells.geometry[i].boundary.intersection(cells.geometry[j].boundary.buffer(TOL))
        if L.length < 5: continue
        Lw = L.intersection(water); Lm = L.difference(water)
        zw = sample(Lw); zm = sample(Lm); za = sample(L)
        zp2 = float(np.percentile(za, 2)) if za.size else np.nan
        low = za[za < zp2 + 0.30] if za.size else za
        rows.append(dict(a=i, b=j, shared_m=round(L.length, 1), water_cross_m=round(Lw.length, 1),
                         edge_z_p2=round(zp2, 3), edge_z_p10=pct(za, 10), edge_z_med=pct(za, 50), low_width_m=round(low.size * L.length / max(za.size, 1), 1),
                         invert_p10=pct(zw, 10), water_z_med=pct(zw, 50), marsh_cross_z_med=pct(zm, 50), marsh_cross_z_p25=pct(zm, 25),
                         band_low_a=band_low(cells.geometry[i], L), band_low_b=band_low(cells.geometry[j], L),
                         dist_m=round(rp[i].distance(rp[j]), 1), ax=round(rp[i].x, 1), ay=round(rp[i].y, 1), bx=round(rp[j].x, 1), by=round(rp[j].y, 1)))
df = pd.DataFrame(rows); df.to_csv(f"{MODEL}/link_geometry_candidates.csv", index=False)
print(f"pairs: {len(df)}, with an open-water crossing of 3 m or more: {(df.water_cross_m >= 3).sum()}")


def reach(sel):
    adj = {}
    for r in df[sel].itertuples(): adj.setdefault(r.a, set()).add(r.b); adj.setdefault(r.b, set()).add(r.a)
    seen = {BAY_CELL}; st = [BAY_CELL]
    while st:
        k = st.pop()
        for q in adj.get(k, set()) - seen: seen.add(q); st.append(q)
    return seen


for thr in (0.2, 0.4, 0.6):
    r = reach(df.edge_z_p2 < thr); print(f"edges lower than {thr} m: cells joined to the bay {len(r) - 1}/{N - 1}")
print("->", f"{MODEL}/link_geometry_candidates.csv")
