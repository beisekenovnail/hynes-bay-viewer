"""Download the open forcing data (NOAA, TWDB, NCEI, USGS) into RAW.
Env: P0, P1 = first and last day (default 2008-12-01 to yesterday), P0_LONG = start of the salinity and river records."""
import io, os, time, urllib.error, urllib.request
import pandas as pd
from config import RAW, TIDE, TIDE_FILL, LEVEL_BACKUP, WIND_BACKUP, SALINITY, WEATHER, RIVER

E = os.environ.get
P0 = pd.Timestamp(E("P0", "2008-12-01"))
P1 = pd.Timestamp(E("P1", f"{pd.Timestamp.today() - pd.Timedelta(days=1):%Y-%m-%d}"))
P0_LONG = pd.Timestamp(E("P0_LONG", "2008-01-01"))
os.makedirs(RAW, exist_ok=True)
NOAA = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?application=icm_marsh_model&time_zone=lst&units=metric&format=csv"


def get(url, tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 400:
                return "Error: no data"   # NOAA answers 400 for a period without data
            if k == tries - 1:
                raise
        except Exception:
            if k == tries - 1:
                raise
        time.sleep(4 * (k + 1))


def noaa(station, products, freq, fname, datum=""):
    """products are tried in order for each period (e.g. verified first, then preliminary)"""
    parts, miss = [], []
    for per in pd.period_range(P0, P1, freq=freq):
        a, b = per.start_time, min(per.end_time, P1 + pd.Timedelta(hours=23))
        for p in products:
            t = get(f"{NOAA}&station={station}&product={p}" + (f"&datum={datum}" if datum else "") + f"&begin_date={a:%Y%m%d}&end_date={b:%Y%m%d}")
            time.sleep(0.8)
            if t.strip() and "Error" not in t[:300] and len(t.strip().splitlines()) > 1:
                d = pd.read_csv(io.StringIO(t), dtype=str, skipinitialspace=True)
                d.columns = [c.strip() for c in d.columns]
                parts.append(d)
                break
        else:
            miss.append(str(per))
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["Date Time"])
    d.to_csv(f"{RAW}/{fname}", index=False)
    print(f"NOAA {station} {products[0]}: {len(d)} rows, periods without data: {len(miss)}", flush=True)


noaa(TIDE, ["water_level"], "M", f"noaa_{TIDE}_water_level.csv", "NAVD")
noaa(TIDE, ["wind&interval=h"], "M", f"noaa_{TIDE}_wind.csv")
noaa(TIDE, ["water_temperature&interval=h"], "M", f"noaa_{TIDE}_water_temperature.csv")
for prod in ("wind", "water_temperature"):
    noaa(TIDE_FILL, [f"{prod}&interval=h"], "Q", f"noaa_{TIDE_FILL}_{prod}.csv")
for st, datum in LEVEL_BACKUP.values():
    noaa(st, ["hourly_height", "water_level&interval=h"], "M", f"noaa_{st}_water_level_hourly_{datum}.csv", datum.upper())
noaa(WIND_BACKUP, ["wind&interval=h"], "M", f"noaa_{WIND_BACKUP}_wind_hourly.csv")

d0, d1, dl = f"{P0:%Y-%m-%d}", f"{P1:%Y-%m-%d}", f"{P0_LONG:%Y-%m-%d}"
for p in ("seawater_salinity", "water_temperature"):
    t = get(f"https://waterdatafortexas.org/coastal/api/stations/{SALINITY}/data/{p}?start_date={dl}&end_date={d1}&binning=day&output_format=csv")
    open(f"{RAW}/twdb_{SALINITY}_{p}_daily.csv", "w").write(t)
    print(f"TWDB {SALINITY} {p}: {sum(1 for l in t.splitlines() if l[:2] in ('19', '20'))} days", flush=True)
for name, st in WEATHER.items():
    t = get(f"https://www.ncei.noaa.gov/access/services/data/v1?dataset=daily-summaries&stations={st}&startDate={d0}&endDate={d1}"
            "&dataTypes=PRCP,TMAX,TMIN&units=metric&format=csv")
    open(f"{RAW}/ncei_{st}_{name}_daily.csv", "w").write(t)
    print(f"NCEI {name}: {len(t.splitlines()) - 1} days", flush=True)
t = get(f"https://waterservices.usgs.gov/nwis/dv/?format=rdb&sites={RIVER}&parameterCd=00060&statCd=00003&startDT={dl}&endDT={d1}")
open(f"{RAW}/usgs_{RIVER}_discharge_daily.rdb", "w").write(t)
print(f"USGS {RIVER}: {sum(1 for l in t.splitlines() if l.startswith('USGS'))} days", flush=True)
print("->", RAW)
