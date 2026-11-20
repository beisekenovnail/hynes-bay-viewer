"""Cumulative connectivity of the cells from the hourly water levels.
Each hour, two neighbouring cells are connected when the water in the higher cell stands above the lowest sill between them
and that cell is not dry (depth over 0.0501 m, as in ICM). The depth over the sill, times 1/24 day, is added to the cell
the water goes into (levels within 1 mm: half to each). Sums are in metre-days (m·d).
Needs RES/tables from collect_outputs.py. Env: YEARS as there. Writes RES/connectivity."""
import os, warnings
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from config import years, CELLS, MODEL, RUNS, RES, BAY_CELL as BAY

warnings.filterwarnings("ignore")
YEARS = years("") or sorted(int(d) for d in os.listdir(RUNS) if d.isdigit())
SPAN = f"{YEARS[0]}_{YEARS[-1]}"; OUT = f"{RES}/connectivity"; os.makedirs(OUT, exist_ok=True)
cells = gpd.read_file(CELLS); N = len(cells)
bed = pd.read_csv(f"{RUNS}/{YEARS[0]}/Cells.csv").bed_elev.values
lk = pd.read_csv(f"{MODEL}/links_table.csv"); lk["a"] = lk[["us", "ds"]].min(axis=1) - 1; lk["b"] = lk[["us", "ds"]].max(axis=1) - 1
pairs = lk.groupby(["a", "b"]).agg(inv=("inv", "min"), types=("type", lambda t: "+".join(sorted({str(x) for x in t})))).reset_index()
H = pd.read_csv(f"{RES}/tables/hourly_stage_m_navd88_by_cell_{SPAN}.csv", index_col=0, parse_dates=True)
S = H[[f"c{c:02d}" for c in range(N)]].values.astype("float64"); T = H.index
hday = (T - pd.Timedelta(hours=1)).normalize(); yrs = hday.year       # hour-ending: 01:00 .. 24:00 belong to that day
A, B = pairs.a.values, pairs.b.values; INV = pairs.inv.values; K = len(pairs)
sa, sb = S[:, A], S[:, B]; upa = sa >= sb; hi = np.where(upa, sa, sb); upbed = np.where(upa, bed[A], bed[B])
u4 = lambda x: np.rint(x * 1e4).astype(np.int64)                         # levels are written to 0.1 mm
conn = (u4(hi) > u4(INV)) & (u4(hi) - u4(upbed) > 501)
dep = np.where(conn, hi - INV, 0.0).astype("float32")
wa = np.where(np.abs(u4(sa) - u4(sb)) < 10, 0.5, upa.astype("float32")).astype("float32"); del sa, sb, hi, upbed
a2b = dep * wa / 24; b2a = dep * (1 - wa) / 24
days = pd.DatetimeIndex(sorted(set(hday))); di = days.get_indexer(hday); ND = len(days)


def daily(x):
    out = np.zeros((ND, x.shape[1]), "float64"); np.add.at(out, di, x); return out


DA, DB = daily(a2b), daily(b2a); HRS = daily(conn.astype("float32")).astype("uint8")
recv = np.zeros((ND, N)); sent = np.zeros((ND, N))
for k in range(K):
    recv[:, B[k]] += DA[:, k]; recv[:, A[k]] += DB[:, k]; sent[:, A[k]] += DA[:, k]; sent[:, B[k]] += DB[:, k]
np.savez_compressed(f"{OUT}/daily_pair_connection_depth_{SPAN}.npz", dates=np.array([d.strftime("%Y-%m-%d") for d in days]), a=A, b=B, inv=INV,
                    types=pairs.types.values, a_to_b_m_day=DA.astype("float32"), b_to_a_m_day=DB.astype("float32"), hours_connected=HRS)

# ---- tables
pairs["cum_a_to_b_m_day"] = DA.sum(0).round(3); pairs["cum_b_to_a_m_day"] = DB.sum(0).round(3)
pairs["cum_total_m_day"] = (DA.sum(0) + DB.sum(0)).round(3); pairs["days_connected"] = (HRS.sum(0) / 24).round(1)
pairs["pct_hours_connected"] = (100 * conn.mean(0)).round(1)
FL = pd.read_csv(f"{RES}/tables/daily_flow_cms_by_link_{SPAN}.csv", index_col=0, parse_dates=True); nyr = len(FL) / 365.25
lk["sgn"] = np.where(lk.us - 1 == lk.a, 1.0, -1.0)
q_ab = {(a, b): sum(FL[f"L{r.link}"].values * r.sgn for r in g.itertuples()) for (a, b), g in lk.groupby(["a", "b"])}
pairs["gross_Mm3_per_yr"] = [np.abs(q_ab[(a, b)]).sum() * 86400 / nyr / 1e6 for a, b in zip(A, B)]
pairs["net_a_to_b_Mm3_per_yr"] = [q_ab[(a, b)].sum() * 86400 / nyr / 1e6 for a, b in zip(A, B)]
pairs.round(4).to_csv(f"{OUT}/pair_cumulative_connectivity_{SPAN}.csv", index=False)
ycol = sorted(set(days.year)); dy = days.year.values
CR = pd.DataFrame({y: recv[dy == y].sum(0) for y in ycol}, index=pd.Index(range(N), name="cell"))
CR.insert(0, f"total_{SPAN}", recv.sum(0)); CR.insert(1, f"sent_total_{SPAN}", sent.sum(0)); CR.round(3).to_csv(f"{OUT}/cell_cumulative_connectivity_by_year_{SPAN}.csv")
cum = pd.DataFrame(np.cumsum(recv, 0), index=days, columns=[f"c{c:02d}" for c in range(N)]).drop(columns=f"c{BAY:02d}")
cum.index.name = "date"; cum.round(3).to_csv(f"{OUT}/daily_cumulative_connectivity_by_cell_{SPAN}.csv")
# hours each cell is joined to the bay through a chain of connected pairs
reach = np.zeros((len(T), N), bool); reach[:, BAY] = True
for _ in range(N):
    before = reach.sum()
    for k in range(K):
        c = conn[:, k]; a, b = A[k], B[k]; reach[:, a] |= reach[:, b] & c; reach[:, b] |= reach[:, a] & c
    if reach.sum() == before: break
area = cells.area.values; marsh = np.arange(N) != BAY
CY = pd.DataFrame({y: 100 * reach[yrs == y].mean(0) for y in sorted(set(yrs))}, index=pd.Index(range(N), name="cell")).round(1)
CY.insert(0, f"all_{SPAN}", (100 * reach.mean(0)).round(1)); CY.to_csv(f"{OUT}/cell_pct_hours_connected_to_bay_by_year_{SPAN}.csv")
Dd = pd.DataFrame({"marsh_area_connected_to_bay_pct": 100 * (reach[:, marsh] * area[marsh]).sum(1) / area[marsh].sum(), "bay_stage_m": S[:, BAY]}, index=hday).groupby(level=0).mean()
Dd["marsh_received_m_day"] = recv[:, marsh].sum(1); Dd.index.name = "date"; Dd.round(4).to_csv(f"{OUT}/daily_marsh_connectivity_{SPAN}.csv")
Yr = pd.DataFrame({"marsh_received_m_day": Dd.marsh_received_m_day.groupby(Dd.index.year).sum(),
                   "pct_marsh_area_connected_to_bay": Dd.marsh_area_connected_to_bay_pct.groupby(Dd.index.year).mean(),
                   "pct_pairs_hours_connected": pd.Series({y: 100 * conn[yrs == y].mean() for y in sorted(set(yrs))}),
                   "bay_stage_mean_m": Dd.bay_stage_m.groupby(Dd.index.year).mean()}).round(2)
Yr.index.name = "year"; Yr.to_csv(f"{OUT}/yearly_connectivity_{SPAN}.csv"); print(Yr.to_string())

# ---- map: connectivity received over the whole period, lines between connected neighbours
rp = cells.geometry.representative_point(); tot = CR[f"total_{SPAN}"].values; ptot = pairs.cum_total_m_day.values
NORM = LogNorm(vmin=0.1, vmax=10 ** np.ceil(np.log10(max(tot[marsh].max(), 1))))
fig, ax = plt.subplots(figsize=(8.8, 10))
cells.assign(v=np.maximum(tot, NORM.vmin)).plot(column="v", cmap="YlGnBu", norm=NORM, ax=ax, edgecolor="k", lw=0.25)
cells.iloc[[BAY]].plot(ax=ax, color="#c6dbef", edgecolor="k", lw=0.25)
for k in np.argsort(ptot):
    if ptot[k] <= 0 or BAY in (A[k], B[k]): continue
    ax.plot([rp[A[k]].x, rp[B[k]].x], [rp[A[k]].y, rp[B[k]].y], color="#c0392b", lw=0.3 + 3.2 * np.log10(1 + ptot[k]) / np.log10(1 + ptot.max()), alpha=0.8)
for c, p in enumerate(rp):
    if c != BAY: ax.text(p.x, p.y, str(c), fontsize=6, ha="center", va="center")
fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap="YlGnBu"), ax=ax, shrink=0.6, label="connectivity received (m·day)")
ax.set_title(f"Cumulative connectivity, {days[0]:%d %b %Y} – {days[-1]:%d %b %Y}", fontsize=10); ax.set_axis_off()
plt.tight_layout(); plt.savefig(f"{OUT}/connectivity_map_{SPAN}.png", dpi=110); plt.close()
print("->", OUT)
