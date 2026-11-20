# ICM coastal marsh model

A modified ICM-Hydro for small coastal marshes. ICM-Hydro is the hydrology part of the Integrated Compartment Model of the Louisiana Coastal Master Plan (Meselhe et al. 2013). The marsh is split into cells. Links join neighbouring cells, either as channels or as flow over the marsh surface. For every cell the model gives daily water level, salinity and temperature. It also gives the daily flow through every link.

## What it is for

Most marshes have few or no gauges. The model fills in the rest of the marsh from open data: the tide gauge, a salinity station, rain and river flow. It answers questions such as:

- which parts of the marsh stay fresh or turn salty in wet and dry years;
- which cells are cut off from the bay, and for how long;
- how far a freshwater delivery would spread.

**Real example.** This repository holds a real example, the Hynes Bay marsh on the Texas coast. The marsh has 55 cells and a bay cell, joined by 235 links. It was run daily from 2009 to 10 November 2026. The interactive map of its results is at https://beisekenovnail.github.io/icm-coastal-marsh-model/.

## Changes to ICM-Hydro

1. The code builds with gfortran. Intel-only format statements were replaced by standard ones.
2. Boundary salinity can change daily (`SalBC_daily.csv`). In the original it is a constant plus a seasonal term fitted for Louisiana.
3. A sign error in evaporation was fixed in `2D_celldQ.f`. With `fpet = 1`, evaporation acted as an inflow of water.

Transport scheme 3 of the original ignores the channel dispersion `Exy`. The example therefore uses scheme 1.

## Quick start: run the example

You need gfortran.

```bash
./model/build.sh
cd example/2025 && ../../model/bin/icm_hydro
```

One year takes about a minute. The results are written into the same folder:

| file | contents |
|---|---|
| `SAL.out`, `STG.out`, `TMP.out` | daily salinity (ppt), water level (m NAVD88) and temperature (°C), one column per cell |
| `STGhr.out` | hourly water level, one column per cell |
| `FLO.out` | daily flow through each link (m³/s), one column per link, in the order of `Links.csv` |

Rows are days and column *k* is cell *k*. In the example, cell 10 is the bay.

```python
import numpy as np
sal = np.loadtxt("example/2025/SAL.out", delimiter=",")   # days x cells
```

## Set up another marsh

You need Python 3 with numpy, pandas, geopandas, rasterio, pyproj and matplotlib.

1. Put the maps into `work/gis/`:
   - `cells.shp`: cell polygons in a metre-based projection, one of them the bay;
   - `dem.tif`: ground elevation in m NAVD88;
   - `water.gpkg`: open water;
   - `landcover/*.shp`: one file per land-cover class;
   - `catchment.kml` (optional): upland area that drains into the marsh.
2. In `scripts/config.py`, set the stations, the bay cell, the inflow cell and the latitude. Any path there can also be set as an environment variable.
3. Run the steps:

```bash
cd scripts
P0=2008-12-01 P1=2026-11-11 python3 fetch_open_data.py   # download open data
T0=2008-12-01 T1=2026-11-11 python3 build_forcing.py     # daily and hourly forcing, gaps filled
python3 cell_stats.py && python3 link_geometry.py          # cells and candidate links
YEARS=2009-2026 python3 make_inputs.py                     # model input files, one folder per year
YEARS=2009-2026 python3 run_model.py                       # years in order, each starts from the end of the last
YEARS=2010-2026 python3 collect_outputs.py                 # tables by cell and link
YEARS=2010-2026 python3 connectivity.py                    # how often and how deep cells are connected
```

The first year is a spin-up and is left out of the results.

The main model settings can be changed through environment variables of `make_inputs.py`:

| setting | example value | meaning |
|---|---|---|
| `N_CH`, `N_MA` | 0.03, 0.3 | Manning n of channel links and of links over the marsh |
| `EXY_CH`, `EXY_MA` | 200, 20 | dispersion, m²/s |
| `IADV` | 1 | transport scheme |
| `DT` | 10 | time step, s |
| `CH_MAX` | 0.6 | a channel link is built where the shared edge is lower than this, m |
| `BAY_BED` | -0.40 | bed of the bay cell, m |

## Folders

| folder | contents |
|---|---|
| `model/` | ICM-Hydro source with the changes, file list and build script |
| `scripts/` | Python steps from open data to results |
| `example/2025/` | input files of the Hynes Bay example for 2025, starting from the state on 1 January |
| `docs/` | the interactive map of the example |

## License and source

The ICM-Hydro source is under GPL-3 (`model/src/LICENSE`). The original description is in `model/src/README.md`.

Meselhe, E., McCorquodale, J.A., Shelden, J., Dortch, M., Brown, T.S., Elkan, P., Rodrigue, M.D., Schindler, J.K., and Wang, Z. 2013. Ecohydrology component of Louisiana's 2012 Coastal Master Plan: mass-balance compartment model. Journal of Coastal Research, Special Issue 67, 16–28.
