"""Collect the yearly outputs into tables: columns c00, c01, ... = cells, L1, L2, ... = links.
Env: YEARS (default: every year folder; leave out a spin-up year). Writes RES/tables/*_<first>_<last>.csv."""
import os
import numpy as np, pandas as pd
from config import years, RUNS, MODEL, RES

YEARS = years("") or sorted(int(d) for d in os.listdir(RUNS) if d.isdigit())
SPAN = f"{YEARS[0]}_{YEARS[-1]}"; TB = f"{RES}/tables"; os.makedirs(TB, exist_ok=True)
LK = pd.read_csv(f"{MODEL}/links_table.csv")


def load(name, freq="D"):
    parts = []
    for y in YEARS:
        t0 = f"{y}-01-01 01:00" if freq == "h" else f"{y}-01-01"   # hourly values are written at the end of each hour
        a = np.atleast_2d(np.loadtxt(f"{RUNS}/{y}/{name}", delimiter=",")); parts.append(pd.DataFrame(a, index=pd.date_range(t0, periods=a.shape[0], freq=freq)))
    return pd.concat(parts)


sal, stg, tmp, flo, sth = load("SAL.out"), load("STG.out"), load("TMP.out"), load("FLO.out"), load("STGhr.out", "h")
cols = [f"c{c:02d}" for c in range(sal.shape[1])]
for d in (sal, stg, tmp, sth): d.columns = cols
flo.columns = [f"L{i}" for i in LK.link]
for d, n in ((sal, "salinity_ppt"), (stg, "stage_m_navd88"), (tmp, "temperature_C")):
    d.index.name = "date"; d.round(4).to_csv(f"{TB}/daily_{n}_by_cell_{SPAN}.csv"); d.resample("MS").mean().round(3).to_csv(f"{TB}/monthly_{n}_by_cell_{SPAN}.csv")
flo.index.name = "date"; flo.round(4).to_csv(f"{TB}/daily_flow_cms_by_link_{SPAN}.csv")
sth.index.name = "datetime"; sth.round(4).to_csv(f"{TB}/hourly_stage_m_navd88_by_cell_{SPAN}.csv")
print(f"{len(sal)} days, {len(cols)} cells, {len(flo.columns)} links ->", TB)
