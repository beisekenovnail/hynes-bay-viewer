"""Paths and site settings shared by the scripts. Every path can be set by an environment variable of the same name.
The values below are those of the Hynes Bay example; change them for another marsh."""
import os

E = os.environ.get
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = E("WORK", os.path.join(REPO, "work"))           # one folder per marsh

# inputs
CELLS = E("CELLS", f"{WORK}/gis/cells.shp")            # cell polygons, cell number = row order, projected CRS in metres
DEM = E("DEM", f"{WORK}/gis/dem.tif")                  # ground elevation, m NAVD88
WATER = E("WATER", f"{WORK}/gis/water.gpkg")           # open-water polygons
LANDCOVER = E("LANDCOVER", f"{WORK}/gis/landcover")    # one polygon file per land-cover class
CATCHMENT = E("CATCHMENT", f"{WORK}/gis/catchment.kml")  # optional upland catchment, parts named 1, 2, ... from the marsh upstream
RAW = E("RAW", f"{WORK}/data/raw")                     # downloaded open data
INIT_SAL = E("INIT_SAL", "")                           # optional CSV cell,sal_ppt: salinity at the cold start

# outputs
FORCING = E("FORCING", f"{WORK}/data/forcing")
MODEL = E("MODEL_DIR", f"{WORK}/model")
RUNS = E("RUNS", f"{MODEL}/runs")
RES = E("RES", f"{WORK}/results")
TAG = E("TAG", "2009_2026")                            # name of the forcing files

# model program
BIN = E("BIN", f"{REPO}/model/bin/icm_hydro")
TEMPLATE = f"{REPO}/model/src/RuncontrolR_template.dat"

# site
BAY_CELL = int(E("BAY_CELL", 10))                      # cell with the open-water boundary
INFLOW_CELL = int(E("INFLOW_CELL", 52))                # cell that receives the upland runoff
LAT = float(E("LAT", 28.445))                          # latitude for evaporation
UPLAND = E("UPLAND", "Mixed_forest,Upland_grasses,Cynodon_dactylon").split(",")   # land-cover classes counted as upland

# open-data stations
TIDE = "8773037"                                       # NOAA: water level, wind, water temperature at the marsh
TIDE_FILL = "8773259"                                  # NOAA: fills wind and water temperature gaps
LEVEL_BACKUP = {"poc": ("8773701", "navd"), "pl": ("8773259", "msl"), "rp": ("8774770", "navd")}   # NOAA hourly levels for long gaps
WIND_BACKUP = "8773701"                                # NOAA: last source for wind
SALINITY = "DELT"                                      # TWDB: boundary salinity and temperature
WEATHER = {"austwell": "USC00410436", "port_aransas_32nne": "USW00023906",
           "point_comfort": "USC00417140", "aransas_nwr": "USC00410305"}   # NCEI daily rain and air temperature
RIVER = "08188800"                                     # USGS daily discharge, for salinity when the salinity station is out


def years(default):
    """YEARS env: '2010-2026' or '2017,2018'"""
    s = E("YEARS", default)
    if not s: return []
    if "-" in s: a, b = s.split("-"); return list(range(int(a), int(b) + 1))
    return [int(v) for v in s.split(",")]
