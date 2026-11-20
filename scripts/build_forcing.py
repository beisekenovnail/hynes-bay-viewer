"""Clean and gap-fill the open data into daily and hourly model forcing.
Env: T0, T1 = first day and the day after the last (YYYY-MM-DD). Writes FORCING/daily_forcing_TAG.csv,
hourly_waterlevel_wind_TAG.csv, inflow_TAG.csv and QC_forcing_TAG.txt."""
import os
import numpy as np, pandas as pd
from config import E, RAW, FORCING, TAG, CATCHMENT, LAT, TIDE, TIDE_FILL, LEVEL_BACKUP, WIND_BACKUP, SALINITY, WEATHER, RIVER

t0, t1 = pd.Timestamp(E("T0", "2008-12-01")), pd.Timestamp(E("T1", "2026-11-11"))
os.makedirs(FORCING, exist_ok=True)
days = pd.date_range(t0, t1, freq="D"); hours = pd.date_range(t0, t1, freq="h"); qc = [f"period {t0.date()} .. {t1.date()} (raw: {RAW})"]


def gapinfo(s):
    na = s.isna(); runs = (na != na.shift()).cumsum()[na]; L = runs.value_counts().max() if na.any() else 0
    return f"missing {na.sum()} ({na.mean() * 100:.1f}%), longest gap {L}"


def noaa(name):
    f = f"{RAW}/{name}.csv"
    if not os.path.exists(f) or os.path.getsize(f) < 50: return None
    d = pd.read_csv(f, skipinitialspace=True); d.columns = [c.strip() for c in d.columns]; d.index = pd.to_datetime(d["Date Time"]); return d


# water level of the main gauge, hourly
wlh = noaa(f"noaa_{TIDE}_water_level")["Water Level"].astype(float).resample("h").mean().reindex(hours); qc.append("main gauge WL hourly: " + gapinfo(wlh))
# gaps over 6 h: lagged linear regression on nearby gauges, fitted on the hours both report
LAGS = (-6, -3, 0, 2, 4, 6, 9, 12, 18, 24)


def gauge(key):
    st, datum = LEVEL_BACKUP[key]; d = noaa(f"noaa_{st}_water_level_hourly_{datum}")
    return None if d is None else pd.to_numeric(d["Water Level"], errors="coerce").resample("h").mean().reindex(hours)


POC, PL, RP = gauge("poc"), gauge("pl"), gauge("rp")


def design(series):
    cols = {}
    for nm, x in series.items():
        for L in LAGS: cols[f"{nm}{L:+d}"] = x.shift(L)
        cols[f"{nm}_72h"] = x.rolling(72, min_periods=60).mean()
    X = pd.DataFrame(cols); X.insert(0, "c", 1.0); return X


fills = []
for name, ser in (("A: poc + pl", {"poc": POC, "pl": PL}), ("B: poc", {"poc": POC}), ("C: rp", {"rp": RP})):
    if any(v is None for v in ser.values()): continue
    X = design(ser); okx = X.notna().all(axis=1); fit = okx & wlh.notna()
    if fit.sum() < 2000: continue
    tr = fit & (hours < "2022-01-01"); te = fit & (hours >= "2022-01-01")
    if name.startswith("C"): tr = fit & (hours >= "2016-09-01"); te = fit & (hours < "2014-08-01")   # rp: tested on the years before the gap
    b = np.linalg.lstsq(X[tr].values, wlh[tr].values, rcond=None)[0]; e = X[te].values @ b - wlh[te].values
    ed = (pd.Series(X[te].values @ b, index=hours[te]).resample("D").mean() - wlh[te].resample("D").mean()).dropna()
    b = np.linalg.lstsq(X[fit].values, wlh[fit].values, rcond=None)[0]; pred = pd.Series(X.values @ b, index=hours).where(okx)
    fills.append((name, pred))
    qc.append(f"WL model {name}: fit {int(fit.sum())} h; hold-out {int(te.sum())} h: hourly RMSE {np.sqrt(np.mean(e ** 2)) * 100:.1f} cm, "
              f"r {np.corrcoef(X[te].values @ b, wlh[te].values)[0, 1]:.3f}; daily-mean RMSE {np.sqrt(np.mean(ed ** 2)) * 100:.1f} cm, bias {ed.mean() * 100:+.1f} cm")
gap = wlh.isna(); runs = (gap != gap.shift()).cumsum(); glen = gap.groupby(runs).transform("sum"); longgap = gap & (glen > 6)
wl_src = pd.Series(np.where(wlh.notna(), "gauge", ""), index=hours, dtype=object)
for name, pred in fills:
    m = longgap & wlh.isna() & pred.notna(); wlh[m] = pred[m]; wl_src[m] = name.split(":")[0]; qc.append(f"WL: {int(m.sum())} h filled by model {name}")
m = wlh.isna(); wlh = wlh.interpolate(limit_direction="both"); wl_src[m & (wl_src == "")] = "interpolated"; qc.append(f"WL: {int(m.sum())} h interpolated")


# wind: main gauge, gaps from the fill gauge and then the backup, speed scaled over the hours both report
def uv(w):
    spd = pd.to_numeric(w["Speed"], errors="coerce"); dd = pd.to_numeric(w[[c for c in w.columns if c.startswith("Direction")][0]], errors="coerce")
    return (-spd * np.sin(np.radians(dd))).resample("h").mean(), (-spd * np.cos(np.radians(dd))).resample("h").mean()


u, v = [s.reindex(hours) for s in uv(noaa(f"noaa_{TIDE}_wind"))]; qc.append("main gauge wind hourly: " + gapinfo(u))
for st, fname in ((TIDE_FILL, f"noaa_{TIDE_FILL}_wind"), (WIND_BACKUP, f"noaa_{WIND_BACKUP}_wind_hourly")):
    w = noaa(fname)
    if w is None: continue
    up, vp = [s.reindex(hours) for s in uv(w)]; both = u.notna() & up.notna()
    k = float(np.hypot(u[both], v[both]).mean() / np.hypot(up[both], vp[both]).mean()) if both.sum() > 100 else 1.0
    fill = u.isna() & up.notna(); u[fill] = k * up[fill]; v[fill] = k * vp[fill]
    qc.append(f"wind: {int(fill.sum())} h filled from {st} (speed x {k:.2f}, {int(both.sum())} common hours)")
qc.append(f"wind: {int(u.isna().sum())} h still missing, interpolated"); u = u.interpolate(limit_direction="both"); v = v.interpolate(limit_direction="both")


# water temperature: main gauge, then the fill gauge with its median offset, then the salinity station
def twdb(p):
    d = pd.read_csv(f"{RAW}/twdb_{SALINITY}_{p}_daily.csv", comment="#"); d["d"] = pd.to_datetime(d.datetime, utc=True).dt.tz_convert(None).dt.normalize()
    return d.groupby("d").value.mean()


twater = noaa(f"noaa_{TIDE}_water_temperature")["Water Temperature"].astype(float).resample("D").mean().reindex(days); qc.append("main gauge Twater daily: " + gapinfo(twater))
wt = noaa(f"noaa_{TIDE_FILL}_water_temperature")
if wt is not None:
    tp = pd.to_numeric(wt["Water Temperature"], errors="coerce").resample("D").mean().reindex(days); both = twater.notna() & tp.notna()
    off = float((twater[both] - tp[both]).median()) if both.sum() > 30 else 0.0; fill = twater.isna() & tp.notna(); twater[fill] = tp[fill] + off
    qc.append(f"water temperature: {int(fill.sum())} d filled from {TIDE_FILL} (offset {off:+.2f} C, {int(both.sum())} common days)")
tdl = twdb("water_temperature").reindex(days); fill = twater.isna() & tdl.notna(); twater[fill] = tdl[fill]
qc.append(f"water temperature: {int(fill.sum())} d filled from {SALINITY}"); twater = twater.interpolate(limit_direction="both")

# boundary salinity: measured where available (gaps up to 3 d interpolated), else estimated from river flow
Sm = twdb("seawater_salinity")
BAD = ("2022-06-12", "2022-07-07")   # sonde drift, removed
nb = int(Sm[BAD[0]:BAD[1]].notna().sum()); Sm[BAD[0]:BAD[1]] = np.nan; qc.append(f"salinity: {nb} days {BAD[0]}..{BAD[1]} removed (sonde drift)")
qU = pd.read_csv(f"{RAW}/usgs_{RIVER}_discharge_daily.rdb", sep="\t", comment="#").iloc[1:]; qU.index = pd.to_datetime(qU.datetime)
QR = pd.to_numeric(qU[[c for c in qU.columns if c.endswith("_00060_00003")][0]], errors="coerce") * 0.0283168   # m3/s
QU = QR.reindex(pd.date_range(QR.index.min(), QR.index.max())).interpolate(limit=5)


def feats(idx):   # log of the 7, 30 and 90 day mean river flow, and the season
    lq = lambda n: np.log10(QU.rolling(n, min_periods=int(n * 0.8)).mean().clip(lower=1.0)).reindex(idx).values
    doy = 2 * np.pi * idx.dayofyear.values / 365.25
    return pd.DataFrame({"c": 1.0, "q7": lq(7), "q30": lq(30), "q90": lq(90), "sin": np.sin(doy), "cos": np.cos(doy)}, index=idx)


Xs = feats(Sm.index); ok = Xs.notna().all(axis=1) & Sm.notna()
fitb = lambda m: np.linalg.lstsq(Xs[m].values, Sm[m].values, rcond=None)[0]
tr, te = ok & (Sm.index < "2020-01-01"), ok & (Sm.index >= "2020-01-01"); bb = fitb(tr); pt = np.clip(Xs[te].values @ bb, 0, 36); et = pt - Sm[te].values
qc.append(f"salinity-flow relation, hold-out from 2020 ({int(te.sum())} d): r={np.corrcoef(pt, Sm[te])[0, 1]:.2f}, RMSE={np.sqrt(np.mean(et ** 2)):.2f} ppt, bias={et.mean():+.2f} ppt")
bb = fitb(ok); pa = np.clip(Xs[ok].values @ bb, 0, 36); ea = pa - Sm[ok].values
qc.append(f"salinity-flow relation, all {int(ok.sum())} days: r={np.corrcoef(pa, Sm[ok])[0, 1]:.2f}, RMSE={np.sqrt(np.mean(ea ** 2)):.2f} ppt; "
          f"coefficients {dict(zip(Xs.columns, np.round(bb, 2)))}")
meas = Sm.reindex(days); qc.append("measured salinity daily: " + gapinfo(meas))
est = pd.Series(np.clip(feats(days).values @ bb, 0, 36), index=days)
sal = meas.interpolate(limit=3, limit_area="inside"); src = pd.Series(np.where(sal.notna(), "measured", "estimated"), index=days)
sal = sal.fillna(est).interpolate(limit_direction="both").clip(lower=0)
qc.append(f"boundary salinity: {int((src == 'measured').sum())} d measured, {int((src == 'estimated').sum())} d estimated from river flow")


# rain and air temperature: nearest stations first
def ncei(name):
    f = f"{RAW}/ncei_{WEATHER[name]}_{name}_daily.csv"
    if not os.path.exists(f) or os.path.getsize(f) < 50: return pd.DataFrame(index=days, columns=["PRCP", "TMAX", "TMIN"], dtype=float)
    d = pd.read_csv(f); d.index = pd.to_datetime(d.DATE); d = d.reindex(days)
    for k in ("PRCP", "TMAX", "TMIN"):
        if k not in d: d[k] = np.nan
    return d


aus, anw, crn, pc = ncei("austwell"), ncei("aransas_nwr"), ncei("port_aransas_32nne"), ncei("point_comfort")
pr = aus.PRCP.fillna(anw.PRCP).fillna(crn.PRCP).fillna(pc.PRCP); qc.append(f"precipitation: {int(aus.PRCP.isna().sum())} d from other stations, {int(pr.isna().sum())} d set 0"); pr = pr.fillna(0.0)
tmax = crn.TMAX.fillna(anw.TMAX).fillna(aus.TMAX).fillna(pc.TMAX).interpolate(limit_direction="both")
tmin = crn.TMIN.fillna(anw.TMIN).fillna(aus.TMIN).fillna(pc.TMIN).interpolate(limit_direction="both"); tair = (tmax + tmin) / 2
# potential evaporation, Hargreaves-Samani
phi = np.radians(LAT); J = days.dayofyear.values; dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365); de = 0.409 * np.sin(2 * np.pi * J / 365 - 1.39); ws = np.arccos(-np.tan(phi) * np.tan(de))
Ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(de) + np.cos(phi) * np.cos(de) * np.sin(ws))
pet = 0.0023 * 0.408 * Ra * (tair.values + 17.8) * np.sqrt(np.clip(tmax.values - tmin.values, 0, None))
qg = QR.reindex(days)

# upland runoff: SCS curve number with the 5-day antecedent rain, then linear reservoirs in series from the top part down
P = pr.values; trib = np.zeros_like(P)
if os.path.exists(CATCHMENT):
    import geopandas as gpd
    KC = gpd.read_file(CATCHMENT); KC = KC.to_crs(KC.estimate_utm_crs("NAD83")); KC = KC.sort_values("Name", ascending=False)   # top part first
    AREA = KC.area.values; CN2 = float(E("TRIB_CN", "84")); KRES = float(E("TRIB_K", "0.75"))
    P5 = pd.Series(P, index=days).shift(1).rolling(5, min_periods=1).sum().fillna(0).values; grow = (days.month >= 4) & (days.month <= 10)
    amc = np.where(P5 < np.where(grow, 35.6, 12.7), 1, np.where(P5 > np.where(grow, 53.3, 27.9), 3, 2))
    CN = np.select([amc == 1, amc == 3], [CN2 / (2.281 - 0.01281 * CN2), CN2 / (0.427 + 0.00573 * CN2)], CN2); S = 25400 / CN - 254
    Qmm = np.where(P > 0.2 * S, (P - 0.2 * S) ** 2 / (P + 0.8 * S), 0.0)
    a = np.exp(-1 / KRES); store = np.zeros(len(AREA)); out = np.zeros((len(P), len(AREA)))
    for i, qm in enumerate(Qmm):   # exact daily solution of dS/dt = I - S/k (m3)
        up = 0.0
        for j, A in enumerate(AREA):
            I = qm / 1000 * A + up; e = store[j] * a + I * KRES * (1 - a); out[i, j] = store[j] + I - e; store[j] = e; up = out[i, j]
    trib = out[:, -1] / 86400
    inflow = pd.DataFrame({"precip_mm": P.round(2), "amc": amc, "cn": np.round(CN, 1), "runoff_mm": np.round(Qmm, 2), "inflow_cms": np.round(trib, 5)}, index=days)
    inflow.index.name = "date"; inflow.to_csv(f"{FORCING}/inflow_{TAG}.csv")
    qc.append(f"upland inflow: {AREA.sum() / 1e4:.1f} ha, CN {CN2:g}, {len(AREA)} reservoirs k {KRES:g} d; runoff {Qmm.sum() / max(P.sum(), 1e-9) * 100:.1f}% of rain; "
              f"mean {trib.mean():.4f} m3/s, max {trib.max():.3f} m3/s")

daily = pd.DataFrame({"precip_mm": pr.round(2), "pet_mm": np.round(pet, 2), "tair_C": tair.round(2), "tmax_C": tmax.round(2), "tmin_C": tmin.round(2),
                      "twater_C": twater.round(2), "sal_bc_ppt": sal.round(3), "sal_bc_source": src, "river_Q_cms": qg.round(2),
                      "wl_daily_mean_navd88_m": wlh.resample("D").mean().reindex(days).round(4), "trib_Q_cms": np.round(trib, 4)}, index=days)
daily.index.name = "date"; daily.to_csv(f"{FORCING}/daily_forcing_{TAG}.csv")
hourly = pd.DataFrame({"wl_navd88_m": wlh.round(4), "wl_source": wl_src, "wind_u_ms": u.round(3), "wind_v_ms": v.round(3)}, index=hours)
hourly.index.name = "datetime"; hourly.to_csv(f"{FORCING}/hourly_waterlevel_wind_{TAG}.csv")
open(f"{FORCING}/QC_forcing_{TAG}.txt", "w").write("\n".join(qc) + "\n"); print("\n".join(qc))
