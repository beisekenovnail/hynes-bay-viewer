"""Write the ICM-Hydro input files, one run folder per year (RUNS/<year>). ICM numbers cells from 1: ICM id = cell + 1.
Needs MODEL/cell_stats.csv, MODEL/link_geometry_candidates.csv and the forcing files. Env: YEARS (e.g. 2009-2026, the first
year is a cold start) and the parameters below."""
import os, shutil, warnings
import numpy as np, pandas as pd, geopandas as gpd
from config import E, years, CELLS, MODEL, FORCING, TAG, RUNS, TEMPLATE, INIT_SAL, BAY_CELL as BAY, INFLOW_CELL

warnings.filterwarnings("ignore")
YEARS = years("2009-2026")
N_CH = float(E("N_CH", 0.03)); N_MA = float(E("N_MA", 0.3))            # Manning n: channel / over-marsh links
EXY_CH = float(E("EXY_CH", 200)); EXY_MA = float(E("EXY_MA", 20))      # dispersion, m2/s
IADV = int(E("IADV", 1))                                                # transport: 1 uses the dispersion, 3 does not
DT = int(E("DT", 10)); MINWATER = float(E("MINWATER", 5000))            # time step (s); minimum open water per cell (m2)
CH_MAX = float(E("CH_MAX", 0.6))                                        # edge elevation (m) below which a channel link is built
BAY_BED = float(E("BAY_BED", -0.40))                                    # bed of the bay cell (m)
BAY_SILL = int(E("BAY_SILL", 1))          # 1: a bay link without a creek mouth gets the ground 20-50 m inside the marsh as its sill
NOPOND_BED = int(E("NOPOND_BED", 1))      # 1: a cell without ponds sits 5 cm below its lowest channel link
BAY_LINK_LEN = E("BAY_LINK_LEN", "shore")  # shore: bay link as long as the way to the cell's own shoreline; centroid: to the bay centre

cells = gpd.read_file(CELLS); cells["cell"] = range(len(cells)); N = len(cells)
cs = pd.read_csv(f"{MODEL}/cell_stats.csv").set_index("cell")
lg = pd.read_csv(f"{MODEL}/link_geometry_candidates.csv")
daily = pd.read_csv(f"{FORCING}/daily_forcing_{TAG}.csv", parse_dates=["date"], index_col="date")
hourly = pd.read_csv(f"{FORCING}/hourly_waterlevel_wind_{TAG}.csv", parse_dates=["datetime"], index_col="datetime")
if "trib_Q_cms" not in daily: daily["trib_Q_cms"] = 0.0
s0 = pd.read_csv(INIT_SAL).set_index("cell").sal_ppt if INIT_SAL else pd.Series(dtype=float)

# ---- cells
rows = []
for c in range(N):
    r = cs.loc[c]; A = cells.geometry[c].area; P = cells.geometry[c].length
    if c == BAY:
        fw, fu, fm = 1.0, 0.0, 0.0; bed = mz = BAY_BED; msd = 0.0
    else:
        fw = float(r.f_water); fu = float(r.f_upland); fm = max(0.0, 1.0 - fw - fu)
        mz = float(r.marsh_z_mean) if np.isfinite(r.marsh_z_mean) else float(r.upland_z)
        msd = float(r.marsh_z_sd) if np.isfinite(r.marsh_z_sd) else 0.1
        bed = float(r.bed_z_water) if (np.isfinite(r.bed_z_water) and fw > 0.005) else mz - msd - 0.25
        bed = min(bed, mz - msd - 0.05)
    rows.append(dict(icm=c + 1, cell=c, A=A, fw=fw, fu=fu, fm=fm, ared=0.1 * fm * A, P=P, bed=bed, mz=mz, msd=msd))
C = pd.DataFrame(rows)

# ---- links: type 1 = channel where the edge is lower than CH_MAX, type 11 = flow over the marsh along the rest of the edge
SHORE = {}
if BAY_LINK_LEN == "shore":
    rp_ = cells.geometry.representative_point()
    for c in range(N):
        if c == BAY: continue
        edge = cells.geometry[c].boundary.intersection(cells.geometry[BAY].buffer(3.0))
        if not edge.is_empty and edge.length >= 5: SHORE[c] = max(50.0, rp_[c].distance(edge.interpolate(0.5, normalized=True)))
L = []; FW = dict(zip(C.cell, C.fw)); dry_drop = []; bay_raised = []
for r in lg.itertuples():
    a, b = int(r.a), int(r.b); ia, ib = a + 1, b + 1; xy = (r.ax, r.ay, r.bx, r.by); zc = float(r.edge_z_p2); med = float(r.edge_z_med)
    if BAY in (a, b) and (b if a == BAY else a) in SHORE: r = r._replace(dist_m=SHORE[b if a == BAY else a])
    # a cell without open water whose ground 20-50 m in stays above CH_MAX has no channel into it
    dry = [c for c, bl in ((a, r.band_low_a), (b, r.band_low_b)) if FW.get(c, 1) < 0.05 and np.isfinite(bl) and bl >= CH_MAX]
    if dry and zc < CH_MAX: dry_drop.append((a, b)); zc = CH_MAX
    if BAY_SILL and BAY in (a, b) and r.water_cross_m < 3:   # no creek mouth: the low edge is the bay's own water
        bl = r.band_low_b if a == BAY else r.band_low_a
        if np.isfinite(bl) and bl > zc: bay_raised.append((a, b)); zc = bl
    if zc < CH_MAX:
        chw = min(max(r.low_width_m, 3.0), 60.0); L.append((ia, ib, *xy, 1, zc, max(med, zc + 0.2), r.dist_m, chw, N_CH, EXY_CH))
    else: chw = 0.0
    if r.shared_m - chw >= 5:
        inv = float(r.marsh_cross_z_p25) if np.isfinite(r.marsh_cross_z_p25) else med
        if chw > 0: inv = max(inv, L[-1][7] + 0.1)
        L.append((ia, ib, *xy, 11, inv, max(med, inv + 0.05), r.dist_m, r.shared_m - chw, N_MA, EXY_MA))
LK = pd.DataFrame(L, columns=["us", "ds", "usx", "usy", "dsx", "dsy", "type", "inv", "bank", "len", "wid", "n", "exy"])
LK.insert(0, "link", range(1, len(LK) + 1)); M = len(LK); BC = [BAY + 1]
print("channel links dropped (dry cell):", dry_drop, "| bay links with a raised sill:", bay_raised)
bl = LK[(LK.us == BAY + 1) | (LK.ds == BAY + 1)].inv
if len(bl) and C.bed[BAY] + 0.05 > bl.min(): print(f"WARNING: bay bed {C.bed[BAY]:.2f} m is not 5 cm below the lowest bay-link sill {bl.min():.3f} m")
if NOPOND_BED:
    low = []; nomouth = {(int(r.a) + 1, int(r.b) + 1) for r in lg.itertuples() if BAY in (r.a, r.b) and r.water_cross_m < 3}
    for c in C.index[(C.fw <= 0.005) & (C.cell != BAY)]:
        inv = LK[((LK.us == c + 1) | (LK.ds == c + 1)) & (LK.type == 1) & ~LK.apply(lambda r: (min(r.us, r.ds), max(r.us, r.ds)) in nomouth, axis=1)].inv
        if len(inv) and C.at[c, "bed"] > inv.min() - 0.05: low.append(int(C.at[c, "cell"])); C.at[c, "bed"] = inv.min() - 0.05
    print("cells without ponds lowered below their lowest channel link:", low)
LK.to_csv(f"{MODEL}/links_table.csv", index=False); C.to_csv(f"{MODEL}/cells_table.csv", index=False)
print(f"links: {M} (channels {int((LK.type == 1).sum())}, over the marsh {int((LK.type == 11).sum())})")
hdr = ("Compartment,TotalArea,AreaWaterPortion,AreaUplandPortion,AreaMarshPortion,MarshEdgeLength,perimeter,bed_elev,bed_depth,"
       "bed_bulk_density,percentForETcalc,sand_adaptionCoeff,Runoff_coeff,RainGage,WindGage,ETGage,ka,cf,sedn,sedcalib,alphased,"
       "KKa,KKdepth,MEE,CSSos,marsh_elev_mean,marsh_elev_stdv,Esho,depo_on_off,marsh_elev_adjust")
cells_csv = [hdr] + [f"{r.icm},{r.A:.1f},{r.fw:.4f},{r.fu:.4f},{r.fm:.4f},{r.ared:.1f},{r.P:.1f},{r.bed:.3f},0.05,0.835,10,0.25,0,"
                     f"1,1,1,0.023,0.001,0.5,1.0E-09,0.008,0.05,0.05,0,-9999,{r.mz:.3f},{r.msd:.3f},0.2,1,{-r.msd:.3f}" for r in C.itertuples()]
links_csv = ["Link,USnode,DSnode,USx,USy,DSx,DSy,type,Latr1,Latr2,Latr3,Latr4,Latr5,Latr6,Latr7,Latr8,Exy,Latr9,Latr10,fa_mult"] + [
    f"{r.link},{r.us},{r.ds},{r.usx:.1f},{r.usy:.1f},{r.dsx:.1f},{r.dsy:.1f},{r.type},{r.inv:.3f},{r.bank:.3f},{r.len:.1f},{r.wid:.1f},{r.n},0,0,0,{r.exy},0,0,1.0"
    for r in LK.itertuples()]

# ---- 500 m grid over the cells (needed by ICM, not used here)
bx0, by0, bx1, by1 = cells.total_bounds
x0, y0 = int(np.floor((bx0 - 250) / 500) * 500 + 250), int(np.floor((by0 - 250) / 500) * 500 + 250)
nc, nr = int(np.ceil((bx1 - x0) / 500)), int(np.ceil((by1 - y0) / 500))
gx = [x0 + 250 + 500 * i for i in range(nc)]; gy = [y0 + 250 + 500 * j for j in range(nr)]
cent = gpd.GeoSeries(gpd.points_from_xy([x for j in range(nr - 1, -1, -1) for x in gx], [gy[j] for j in range(nr - 1, -1, -1) for x in gx]), crs=cells.crs)
rp = cells.geometry.representative_point(); grid = []
for k, p in enumerate(cent, 1):
    ins = cells.contains(p); c = int(cells[ins].cell.iloc[0]) if ins.any() else int(cells.distance(p).idxmin())
    lks = LK[(LK.us == c + 1) | (LK.ds == c + 1)]
    if lks.empty: lks = LK.iloc[:1]
    l1 =int(lks.link.iloc[0]); mx = (lks.usx.iloc[0] + lks.dsx.iloc[0]) / 2; my = (lks.usy.iloc[0] + lks.dsy.iloc[0]) / 2
    grid.append((k, c + 1, l1, p.distance(rp[c]), np.hypot(p.x - mx, p.y - my), C.bed[c], C.mz[c], 100 * (1 - C.fw[c]), 100 * C.fw[c]))
NG = len(grid)

# ---- writers
tpl = open(TEMPLATE).read().splitlines()


def runcontrol(year, sim):
    vals = {1: year, 2: 1, 3: DT, 4: M, 5: N, 6: 1, 7: 1, 8: 1, 9: 0, 10: sim - 1, 11: 1, 12: 0, 15: 1, 16: 1, 17: 1, 18: 1, 19: 1, 20: 1,
            61: 20, 74: M, 75: N, 76: 1, 77: int(MINWATER), 80: 1, 81: 0, 82: 1, 83: 0, 84: 1.0, 88: sim, 89: 1, 100: 0, 101: 0, 102: 0, 103: 0, 104: 0, 105: 0.05, 106: IADV}
    return "\n".join(("%-11s !%s" % (vals[i], l.split("!", 1)[1])) if i in vals else l for i, l in enumerate(tpl, 1)) + "\n"


def write_year(y):
    D = f"{RUNS}/{y}"; shutil.rmtree(D, ignore_errors=True); os.makedirs(D)
    W = lambda n, t: open(f"{D}/{n}", "w").write(t)
    last = min(daily.index.max(), (hourly.index.max() - pd.Timedelta(hours=24)).normalize())   # a partial last year ends with the data
    d0, d1 = pd.Timestamp(f"{y}-01-01"), min(pd.Timestamp(f"{y}-12-31"), last); SIM = len(pd.date_range(d0, d1))
    dd = daily.loc[d0:d1]; hw = hourly.loc[d0:d1 + pd.Timedelta(hours=23)]; ht = hourly.loc[d0:d1 + pd.Timedelta(days=1)]
    assert len(dd) == SIM and len(hw) == SIM * 24 and len(ht) == SIM * 24 + 1
    W("ICM_info_into_EH.txt", "\n".join(["VegWaveAmp.csv", "VegMeanSal.csv", "VegSummerDepth.csv", "VegSummerSal.csv", "VegSummerTemp.csv",
      "VegTreeEstCond.csv", "VegBIHeight.csv", "VegPerLand.csv", str(NG), str(nr), str(nc), str(x0), str(y0), "1", "1", "1", "0", "0", "hydro_run.log"]) + "\n")
    W("RuncontrolR.dat", runcontrol(y, SIM)); W("Cells.csv", "\n".join(cells_csv) + "\n"); W("Links.csv", "\n".join(links_csv) + "\n")
    W("Fetch.csv", "\n".join(f"{r.icm}," + ",".join([f"{max(50.0, np.sqrt(r.fw * r.A)):.0f}"] * 16) for r in C.itertuples()) + "\n")
    W("LinksClosedHours.dat", "Link," + ",".join(f"h{h}" for h in range(23)) + "\n" + "\n".join(f"{i}," + ",".join(["0"] * 23) for i in range(1, M + 1)) + "\n")
    W("LockControlObservedData.csv", "row,date,lock1\n")
    W("TribQ.csv", "1\n" + "\n".join(f"{q:.4f}" for q in dd.trib_Q_cms) + "\n")
    W("TribS.csv", "1\n" + "\n".join(["0.0"] * SIM) + "\n"); W("TribF.csv", "1\n" + "\n".join(["0.0"] * SIM) + "\n")
    W("QMult.csv", "node,t1\n" + "\n".join(f"{j},{1.0 if j == INFLOW_CELL + 1 else 0.0}" for j in range(1, N + 1)) + "\n")
    W("Precip.csv", "day,g1\n" + "\n".join(f"{k + 1},{p:.2f}" for k, p in enumerate(dd.precip_mm)) + "\n")
    W("PET.csv", "day,g1\n" + "\n".join(f"{k + 1},{p:.2f}" for k, p in enumerate(dd.pet_mm)) + "\n")
    W("Meteorology.csv", "air_temp,trib_water_temp\n" + "\n".join(f"{a:.2f},{a:.2f}" for a in dd.tair_C) + "\n")
    W("MissRToC.csv", "\n".join(f"{a:.2f}" for a in dd.tair_C) + "\n")
    W("WindVectorsX.csv", "idx,date,w1\n" + "\n".join(f"{k + 1},{t:%Y-%m-%dT%H},{u:.3f}" for k, (t, u) in enumerate(hw.wind_u_ms.items())) + "\n")
    W("WindVectorsY.csv", "idx,date,w1\n" + "\n".join(f"{k + 1},{t:%Y-%m-%dT%H},{u:.3f}" for k, (t, u) in enumerate(hw.wind_v_ms.items())) + "\n")
    W("TideData.csv", "idx,date,g1\n" + "\n".join(f"{k + 1},{t:%Y-%m-%dT%H},{z:.4f}" for k, (t, z) in enumerate(ht.wl_navd88_m.items())) + "\n")
    W("surge.csv", "idx,date,b1\n" + "\n".join(f"{k + 1},{t:%Y-%m-%dT%H},0.0" for k, t in enumerate(ht.index)) + "\n")
    W("TideTranspose.csv", "gage,compartment,transpose_hr\n" + "".join(f"{i + 1},{c},0\n" for i, c in enumerate(BC))); W("TideWeight.csv", "comp,Ewght,Ecomp,Wwght,Wcomp\n")
    W("KBC.dat", " ".join(str(c) for c in BC) + "\n")
    W("SBC.dat", "".join(f"{c},{dd.sal_bc_ppt.iloc[0]:.3f},0.0,0.1,0.05,0.5,0.1,7.0,5.0,0.02,0.02,0.0\n" for c in BC))
    W("SalBC_daily.csv", "day,bc1\n" + "\n".join(f"{k + 1},{s:.3f}" for k, s in enumerate(dd.sal_bc_ppt)) + "\n")
    W("BCToC2.dat", "\n".join(f"{t:.2f}" for t in dd.twater_C) + "\n")
    W("links_to_write.csv", "link\n" + "\n".join(str(i) for i in range(1, M + 1)) + "\n")
    W("hourly_stage_to_write.csv", "compartment\n" + "\n".join(str(j) for j in range(1, N + 1)) + "\n")
    W("grid_lookup_500m.csv", "grid,compartment" + "".join(f",link{i}" for i in range(1, 21)) + "\n" + "\n".join(f"{g[0]},{g[1]},{g[2]}" + ",-9999" * 19 for g in grid) + "\n")
    W("grid_interp_dist_500m.csv", "grid,dist_comp" + "".join(f",dist_link{i}" for i in range(1, 21)) + "\n" + "\n".join(f"{g[0]},{g[3]:.1f},{g[4]:.1f}" + ",-9999" * 19 for g in grid) + "\n")
    W("grid_data_500m.csv", "grid,bed_elev,land_elev,pct_land,x,pct_water\n" + "\n".join(f"{g[0]},{g[5]:.3f},{g[6]:.3f},{g[7]:.1f},0,{g[8]:.1f}" for g in grid) + "\n")
    W("grid_IDs_Veg_matrix.csv", "\n".join(" ".join(str(r * nc + c + 1) for c in range(nc)) for r in range(nr)) + "\n")
    if y == YEARS[0]:   # cold start: bay level everywhere, salinity from INIT_SAL (or 5 ppt); a spin-up year removes its effect
        z0 = float(hourly.wl_navd88_m.iloc[0]); t0 = float(dd.twater_C.iloc[0]); sb = float(dd.sal_bc_ppt.iloc[0])
        hs = ["comp,stage,sal,css1,css2,css3,css4,tmp," + ",".join(f"chem{i}" for i in range(1, 15)) + ",marsh_stage"]
        for r in C.itertuples():
            es = max(z0, r.bed + 0.06); sal = sb if r.cell == BAY else s0.get(r.cell, np.nan)
            sal = 5.0 if not np.isfinite(sal) else sal
            hs.append(f"{r.icm},{es:.3f},{sal:.2f},0,0,0,0,{t0:.2f}," + ",".join(["0.01"] * 14) + f",{max(es, r.mz - r.msd):.3f}")
        W("hotstart_in.dat", "\n".join(hs) + "\n")
    return SIM


for y in YEARS: print(y, write_year(y), "days ->", f"{RUNS}/{y}")
