# Runnable Examples & Notebooks

The written guides in this repo are reference manuals. This page indexes the
**runnable** companions: Jupyter notebooks that execute end-to-end against
**real public data** (with a committed cached fallback so they also run
offline), plus model post-processing notebooks that use small synthetic-but-
realistic NetCDF files.

## Setup

Always work inside a virtual environment.

```bash
cd meteo-notes
python3 -m venv venv
source venv/bin/activate            # Windows: venv\Scripts\activate
pip install -r requirements.txt     # core stack + Jupyter
pip install -e .                    # makes `import meteo_examples` work anywhere

# register the venv as a Jupyter kernel (optional but recommended)
python -m ipykernel install --user --name meteo-notes --display-name "meteo-notes"

jupyter lab        # or: jupyter notebook
```

Some notebooks need optional extras — each notebook says so at the top:

```bash
pip install -e ".[maps]"     # cartopy (spatial maps)
pip install -e ".[grib]"     # cfgrib + eccodes (raw GRIB2)
pip install -e ".[ocean]"    # argopy, utide, erddapy (ocean obs)
```

## How data works (live + cached fallback)

Helpers in `meteo_examples/dataio.py` try a **live download** from a free,
no-auth source first, then **cache** the result under `meteo_examples/data/`.
If the network is unavailable, they fall back to the committed cached copy, so
every notebook runs reproducibly with or without internet. To force fresh data,
delete the relevant file in `meteo_examples/data/` (or pass `prefer_cache=False`).

Sources: [IEM ASOS](https://mesonet.agron.iastate.edu/request/download.phtml)
(surface obs), [Open-Meteo](https://open-meteo.com/) (archived NWP forecasts &
ensembles), [NDBC](https://www.ndbc.noaa.gov/) (marine buoys).

## Notebook index

### Verification — back [`python-verification-guide.md`](verification/python-verification-guide.md)
| Notebook | What it does | Data |
|---|---|---|
| [`verification/notebooks/01_point_verification.ipynb`](verification/notebooks/01_point_verification.ipynb) | 2-m temperature forecast vs station obs: ME/MAE/RMSE/r, timeseries, scatter, error-by-lead | live |
| [`verification/notebooks/02_categorical_precip.ipynb`](verification/notebooks/02_categorical_precip.ipynb) | Precip contingency tables, POD/FAR/CSI/ETS by threshold, reliability diagram, FSS | live |
| [`verification/notebooks/03_ensemble_crps.ipynb`](verification/notebooks/03_ensemble_crps.ipynb) | Ensemble CRPS, rank histogram, spread–skill | live |

### Observations
| Notebook | What it does | Data |
|---|---|---|
| [`observations/notebooks/surface_obs_qc.ipynb`](observations/notebooks/surface_obs_qc.ipynb) | Pull ASOS, run range/spike/consistency QC, visualize flags — backs [`surface-obs-guide.md`](observations/surface-obs-guide.md) | live |
| [`observations/notebooks/ocean_obs_access.ipynb`](observations/notebooks/ocean_obs_access.ipynb) | NDBC buoy wave/SST access & QC (+ optional Argo/utide) — backs [`ocean-obs-guide.md`](observations/ocean-obs-guide.md) | live |

### Model post-processing (synthetic sample output)
| Notebook | What it does |
|---|---|
| [`models/notebooks/wrf_postprocessing.ipynb`](models/notebooks/wrf_postprocessing.ipynb) | Open a `wrfout`, plot T2/precip/10-m wind, derive period precip |
| [`models/notebooks/roms_croco_output.ipynb`](models/notebooks/roms_croco_output.ipynb) | Open a ROMS/CROCO history file, plot SST, vertical section |
| [`models/notebooks/wave_spectra.ipynb`](models/notebooks/wave_spectra.ipynb) | Plot a 2-D directional wave spectrum, integrate to Hs |

> The model notebooks generate their sample NetCDF on first run via
> `meteo_examples.sampledata`. The variable names, dimensions, and attributes
> match real output, so the same code works on your actual WRF/ROMS/SWAN files —
> just point `xr.open_dataset()` at them instead.
