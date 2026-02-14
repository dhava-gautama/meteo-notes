# Ocean Observations: Complete Guide

> A practical guide to ocean observational data — platforms, data sources,
> quality control, and Python-based access. Covers moored and drifting buoys,
> Argo floats, tide gauges, and satellite ocean products, with a focus on
> preparing observations for ROMS, SWAN, and WW3 verification.
>
> **Companion guides:**
> - [Surface Observations Guide](surface-obs-guide.md) — land-based surface obs
> - [Python Forecast Verification Guide](../verification/python-verification-guide.md) — verification metrics
> - [ROMS Guide](../models/roms-model-guide.md) — ocean model requiring verification data
> - [SWAN Guide](../models/swan-model-guide.md) — wave model verification
> - [WW3 Guide](../models/ww3-model-guide.md) — wave model verification
> - [COAWST Guide](../models/coawst-model-guide.md) — coupled model system

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [In-Situ Platforms](#2-in-situ-platforms)
3. [Satellite Ocean Products](#3-satellite-ocean-products)
4. [Data Sources & Archives](#4-data-sources--archives)
5. [Python Environment Setup](#5-python-environment-setup)
6. [Reading Buoy & Mooring Data](#6-reading-buoy--mooring-data)
7. [Reading Argo Float Data](#7-reading-argo-float-data)
8. [Reading Tide Gauge Data](#8-reading-tide-gauge-data)
9. [Reading Satellite Products](#9-reading-satellite-products)
10. [Quality Control](#10-quality-control)
11. [Preparing Obs for Model Verification](#11-preparing-obs-for-model-verification)
12. [Indonesian & Tropical Waters](#12-indonesian--tropical-waters)
13. [End-to-End Example](#13-end-to-end-example)

**Appendices:**
- [A. Variable Reference](#appendix-a-variable-reference)
- [B. Resources & Links](#appendix-b-resources--links)

---

## 1. Introduction

Ocean observations are essential for validating ocean circulation models (ROMS), wave models (SWAN, WW3), and coupled systems (COAWST). Unlike surface meteorological stations, ocean observing platforms are sparse and expensive to maintain, making data access and quality control critical.

### What this guide covers

| Topic | Description |
|-------|-------------|
| In-situ platforms | Moored buoys, drifters, Argo floats, tide gauges, ship reports |
| Satellite products | Altimetry (SSH, SWH), SST, ocean color, scatterometer winds |
| Data access | Python code for ERDDAP, Argo, CMEMS, NDBC, UHSLC |
| Quality control | Automated QC for each platform type |
| Verification prep | Matching obs to ROMS/SWAN/WW3 output grids |

### What this guide does NOT cover

- Atmospheric observations over the ocean — see [Surface Observations Guide](surface-obs-guide.md)
- Acoustic and biogeochemical instruments (CTDs for research cruises)
- Ocean data assimilation (EnOI, 4D-Var) — see model-specific guides

---

## 2. In-Situ Platforms

### 2.1 Moored Buoys

Fixed buoys anchored to the ocean floor, providing continuous time series.

| Network | Region | Variables | Typical Depth | Access |
|---------|--------|-----------|---------------|--------|
| **NDBC** | US waters, global | Wind, waves, SST, pressure | Surface | `ndbc.noaa.gov` |
| **TAO/TRITON** | Tropical Pacific | T, S, currents, air-sea fluxes | 0–500 m | `pmel.noaa.gov/tao` |
| **PIRATA** | Tropical Atlantic | T, S, currents, fluxes | 0–500 m | `pmel.noaa.gov/pirata` |
| **RAMA** | Tropical Indian Ocean | T, S, currents, fluxes | 0–500 m | `pmel.noaa.gov/rama` |
| **OceanSITES** | Global deep ocean | Full water column | 0–5000 m | `oceansites.org` |
| **BMKG buoys** | Indonesian waters | Waves, SST, wind | Surface | BMKG internal |

**What moored buoys measure:**

Surface (meteorological):
- Wind speed and direction (anemometer, typically 4–5 m height)
- Atmospheric pressure
- Air temperature and humidity
- Sea surface temperature (hull sensor, ~1 m depth)

Waves (accelerometer or wave rider):
- Significant wave height (Hs)
- Peak wave period (Tp) and mean period (Tm)
- Wave direction (directional buoys only)
- Spectral energy density S(f) or S(f,θ)

Subsurface (mooring line):
- Temperature at multiple depths
- Salinity (conductivity)
- Current velocity (ADCP or point current meters)

### 2.2 Drifting Buoys (Drifters)

Free-floating surface buoys tracked by satellite (Argos, GPS).

| Program | Variables | Coverage |
|---------|-----------|----------|
| **GDP** (Global Drifter Program) | SST, SLP, currents (from drift) | Global, ~1500 active |
| **SVP** (Surface Velocity Program) | Surface currents, SST | Global |
| **SVP-B** (with barometer) | + sea level pressure | Selected regions |

- Drifters follow surface currents → Lagrangian velocity estimates
- Drogue at 15 m depth (represents mixed-layer flow)
- Report every 1–6 hours via Argos or Iridium

### 2.3 Argo Floats

Autonomous profiling floats — the backbone of subsurface ocean observing.

| Property | Detail |
|----------|--------|
| Active floats | ~4000 globally |
| Profile depth | Standard: 0–2000 m; Deep Argo: 0–6000 m |
| Cycle time | ~10 days (sink → drift at 1000 m → profile to surface → transmit) |
| Variables | Temperature, salinity, pressure; some: dissolved O₂, pH, chlorophyll |
| Accuracy | T: ±0.002°C, S: ±0.01 PSU, P: ±2.4 dbar |
| Data access | Argo GDAC, `argopy` Python library |

**Argo cycle diagram:**

```
Surface: transmit data (6-12 hours)
    |
    v  Descend (~6 hours)
    |
Parking depth (1000 m): drift (~9 days)
    |
    v  Descend to profile start
    |
Profile depth (2000 m)
    |
    ^  Ascend & profile (~6 hours)
    |     Measures T, S every 2 m
Surface: transmit → repeat
```

### 2.4 Tide Gauges

Fixed instruments measuring sea level at coastal locations.

| Network | Coverage | Sampling | Variables |
|---------|----------|----------|-----------|
| **UHSLC** | Global (~400 stations) | Hourly, daily | Sea level (relative) |
| **IOC/GLOSS** | Global (~290 core) | Various | Sea level, some T/S |
| **BMKG tide gauges** | Indonesia (~150+) | 1-min to hourly | Sea level |
| **INA-CBT** | Indonesia (tsunami) | 1-min | Sea level (real-time) |
| **BIG** | Indonesia (geodetic) | 1-min | Sea level |

Types of tide gauges:
- **Float/stilling well** — traditional, measures water height in a tube
- **Pressure sensor** — bottom-mounted, measures hydrostatic pressure
- **Radar/acoustic** — non-contact, measures distance to water surface
- **GNSS reflectometry** — emerging technique using satellite signals

### 2.5 Voluntary Observing Ships (VOS)

| Property | Detail |
|----------|--------|
| Ships | ~4000 globally (SOOP, VOS) |
| Variables | SST (bucket/intake), wind, pressure, waves (visual), air T/Td |
| Format | BUFR (replacing SHIP FM-13) |
| Archive | ICOADS (International Comprehensive Ocean-Atmosphere Data Set) |

### 2.6 Platform Comparison

| Platform | Spatial | Temporal | Depth | Best For |
|----------|---------|----------|-------|----------|
| Moored buoy | Fixed point | Minutes–hourly | Surface (+subsurface) | Wave verification, air-sea fluxes |
| Drifter | Lagrangian | Hours | Surface only | SST, surface currents |
| Argo | Global profiles | ~10 days | 0–2000 m | T/S profiles, subsurface structure |
| Tide gauge | Coastal fixed | Minutes | Surface | Sea level, tides, storm surge |
| VOS | Ship tracks | Hours | Surface | SST, wind (sparse) |
| HF Radar | Coastal 50–200 km | Hourly | Surface | Surface currents (coastal) |

---

## 3. Satellite Ocean Products

### 3.1 Sea Surface Temperature (SST)

| Product | Resolution | Revisit | Sensor | Provider |
|---------|-----------|---------|--------|----------|
| **GHRSST L4 MUR** | 0.01° (~1 km) | Daily | Multi-sensor | JPL/NASA |
| **OSTIA** | 0.05° (~5 km) | Daily | Multi-sensor | Met Office |
| **GOES-16/17 SST** | 2 km | Hourly | ABI (IR) | NOAA |
| **MODIS SST** | 1–4 km | Daily | IR + microwave | NASA |
| **AMSR2 SST** | 25 km | Daily | Microwave | JAXA |
| **Himawari-9 SST** | 2 km | Hourly | AHI (IR) | JMA |

Key distinction:
- **IR-based SST** = skin temperature (~10 μm depth), cloud-affected
- **Microwave SST** = sub-skin (~1 mm), sees through clouds, coarser
- **L4 analysis** = gap-filled, blended product (best for verification)

### 3.2 Sea Surface Height / Altimetry

| Product | Resolution | Provider | Variables |
|---------|-----------|----------|-----------|
| **CMEMS SEALEVEL** | 0.25° | Copernicus | SLA, ADT, geostrophic currents |
| **AVISO+** | 0.25° | CNES | SLA, MDT, along-track |
| **Jason-3 / Sentinel-6** | Along-track | EUMETSAT | SSH, SWH, wind speed |
| **SWOT** | 2 km swath | NASA/CNES | 2D SSH (wide-swath) |

Altimetry provides:
- **SSH** (sea surface height) → for ROMS verification
- **SWH** (significant wave height) → for SWAN/WW3 verification
- **σ₀** (backscatter) → wind speed estimates

### 3.3 Ocean Surface Waves (Satellite)

| Product | Variable | Resolution | Provider |
|---------|----------|-----------|----------|
| **Altimeter SWH** | Hs (along-track) | 7 km | CMEMS |
| **SAR wave mode** | Spectral (Sentinel-1) | ~5 km | ESA |
| **CFOSAT SWIM** | Directional spectrum | ~70 km | CNES/CNSA |
| **ERA5 waves** | Hs, Tp, dir (reanalysis) | 0.5° | ECMWF |
| **GlobWave** | Multi-mission Hs | Various | ESA |

### 3.4 Scatterometer Winds

Ocean surface wind vectors from radar backscatter:

| Sensor | Platform | Resolution | Swath | Provider |
|--------|----------|-----------|-------|----------|
| **ASCAT** | MetOp-B/C | 12.5–25 km | 2×550 km | EUMETSAT |
| **RapidScat** | ISS (ended 2016) | 25 km | 900 km | NASA |
| **HY-2B/C/D** | HY-2 series | 25 km | 1300 km | NSOAS (China) |
| **CFOSAT** | CFOSAT | 25 km | 800 km | CNES/CNSA |

Scatterometer winds are useful for:
- Verifying atmospheric model surface winds over the ocean
- Driving wave models (SWAN, WW3) in hindcast mode
- Computing wind stress for ROMS forcing

### 3.5 Ocean Color

| Product | Resolution | Variable | Provider |
|---------|-----------|----------|----------|
| **MODIS Aqua** | 1 km | Chl-a, Kd | NASA |
| **Sentinel-3 OLCI** | 300 m | Chl-a, TSM | EUMETSAT |
| **Himawari chl-a** | 5 km | Chl-a (hourly) | JAXA |

Primarily for biogeochemical model validation rather than physical oceanography.

---

## 4. Data Sources & Archives

### 4.1 Primary Archives

| Archive | Content | Format | Access Method |
|---------|---------|--------|---------------|
| **NDBC** | US buoys (met + waves) | Text/NetCDF | Web / ERDDAP |
| **Argo GDAC** | All Argo profiles | NetCDF | FTP / `argopy` |
| **CMEMS** | Satellite + in-situ | NetCDF | API (`copernicusmarine`) |
| **UHSLC** | Tide gauge sea level | NetCDF/CSV | FTP / web |
| **ICOADS** | Ship + buoy (historical) | IMMA | NCAR RDA |
| **ERDDAP servers** | Multi-source (federated) | Various | `erddapy` |
| **OceanObs PMEL** | TAO/TRITON/PIRATA/RAMA | NetCDF/ASCII | Web / ERDDAP |

### 4.2 ERDDAP Servers

ERDDAP is a data server framework that standardizes access to diverse oceanographic datasets. Key servers:

| Server | URL | Focus |
|--------|-----|-------|
| **CoastWatch** | `coastwatch.pfeg.noaa.gov/erddap` | US coastal + global SST |
| **IOOS** | `erddap.sensors.ioos.us/erddap` | US integrated ocean obs |
| **PMEL** | `data.pmel.noaa.gov/pmel/erddap` | Tropical moorings |
| **Copernicus** | `nrt.cmems-du.eu/erddap` | European/global |
| **INCOIS** | `erddap.incois.gov.in/erddap` | Indian Ocean |

### 4.3 Copernicus Marine Service (CMEMS)

CMEMS is the primary source for European and global ocean data products:

| Product ID | Description | Resolution |
|------------|-------------|-----------|
| `SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001` | Global L4 SST (OSTIA) | 0.05° daily |
| `SEALEVEL_GLO_PHY_L4_NRT_008_046` | Global sea level (gridded) | 0.25° daily |
| `WAVE_GLO_PHY_SWH_L4_NRT_014_003` | Global SWH (multi-mission) | 0.2° 3-hourly |
| `INSITU_GLO_PHYBGCWAV_DISCRETE_MYNRT_013_030` | Global in-situ (all platforms) | Point data |
| `GLOBAL_ANALYSISFORECAST_PHY_001_024` | GLORYS analysis/forecast | 0.083° daily |

---

## 5. Python Environment Setup

### Conda environment

```yaml
# ocean-obs-env.yml
name: ocean-obs
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - numpy
  - pandas
  - xarray
  - netcdf4
  - scipy
  - matplotlib
  - cartopy
  - requests
  - erddapy          # ERDDAP data access
  - argopy            # Argo float data
  - gsw               # TEOS-10 seawater properties
  - pip:
    - copernicusmarine  # CMEMS data access
    - wavespectra        # wave spectral analysis
    - utide              # tidal harmonic analysis
```

```bash
conda env create -f ocean-obs-env.yml
conda activate ocean-obs
```

### Quick import check

```python
import numpy as np
import pandas as pd
import xarray as xr
import erddapy
import argopy
import gsw
print("All ocean obs packages loaded successfully.")
```

---

## 6. Reading Buoy & Mooring Data

### 6.1 NDBC Buoy Data (Standard Meteorological)

```python
import pandas as pd
import numpy as np

def read_ndbc_stdmet(station_id, year):
    """Read NDBC standard meteorological data for one station-year.

    Parameters
    ----------
    station_id : str — NDBC station ID (e.g., "41001")
    year       : int — year to download

    Returns
    -------
    pd.DataFrame with datetime index and standard met variables
    """
    url = (f"https://www.ndbc.noaa.gov/view_text_file.php?"
           f"filename={station_id}h{year}.txt.gz&dir=data/historical/stdmet/")

    # NDBC header has two lines (column names + units)
    df = pd.read_csv(url, sep=r"\s+", header=[0, 1], na_values=[99, 999,
                     9999, 99.0, 999.0, 9999.0])

    # Flatten multi-level column names
    df.columns = [c[0] for c in df.columns]

    # Build datetime
    if "mm" in df.columns:
        df["datetime"] = pd.to_datetime(
            df[["YY", "MM", "DD", "hh", "mm"]].rename(
                columns={"YY":"year","MM":"month","DD":"day",
                         "hh":"hour","mm":"minute"})
        )
    else:
        df["datetime"] = pd.to_datetime(
            df[["YY", "MM", "DD", "hh"]].rename(
                columns={"YY":"year","MM":"month","DD":"day","hh":"hour"})
        )

    df = df.set_index("datetime")

    # Rename to standard names
    rename = {
        "WDIR": "wind_dir_deg",
        "WSPD": "wind_speed_ms",
        "GST":  "wind_gust_ms",
        "WVHT": "swh_m",         # significant wave height
        "DPD":  "peak_period_s",  # dominant wave period
        "APD":  "mean_period_s",  # average wave period
        "MWD":  "wave_dir_deg",   # mean wave direction
        "PRES": "slp_hpa",
        "ATMP": "air_temp_c",
        "WTMP": "sst_c",         # sea surface temperature
        "DEWP": "dewpoint_c",
        "VIS":  "visibility_nmi",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    return df
```

### 6.2 NDBC Spectral Wave Data

```python
def read_ndbc_spectral(station_id, year):
    """Read NDBC spectral wave data (energy density vs frequency).

    Returns
    -------
    pd.DataFrame with columns: datetime, freq_hz (array), energy (array)
    """
    # Alpha1 (mean direction) and alpha2 (principal direction)
    # are in separate files — here we read just the energy density
    url = (f"https://www.ndbc.noaa.gov/view_text_file.php?"
           f"filename={station_id}w{year}.txt.gz"
           f"&dir=data/historical/swden/")

    lines = pd.io.common.urlopen(url).read().decode().split("\n")

    # First line: header, second line: frequencies
    freq_line = lines[1].strip().split()
    freqs = np.array([float(f) for f in freq_line[5:]])  # Skip date columns

    records = []
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        try:
            dt = pd.Timestamp(
                year=int(parts[0]), month=int(parts[1]),
                day=int(parts[2]), hour=int(parts[3]),
                minute=int(parts[4])
            )
            energy = np.array([float(x) for x in parts[5:]])
            records.append({"datetime": dt, "freq_hz": freqs, "energy": energy})
        except (ValueError, IndexError):
            continue

    return pd.DataFrame(records)
```

### 6.3 ERDDAP Access (Generic)

```python
from erddapy import ERDDAP

def fetch_erddap(server_url, dataset_id, variables, constraints,
                 response="csv"):
    """Fetch data from any ERDDAP server.

    Parameters
    ----------
    server_url  : str — e.g., "https://coastwatch.pfeg.noaa.gov/erddap"
    dataset_id  : str — ERDDAP dataset ID
    variables   : list of str — variable names to fetch
    constraints : dict — e.g., {"time>=": "2024-01-01", "latitude>=": -10}

    Returns
    -------
    pd.DataFrame or xr.Dataset
    """
    e = ERDDAP(server=server_url, protocol="tabledap")
    e.dataset_id = dataset_id
    e.variables = variables
    e.constraints = constraints

    if response == "csv":
        return e.to_pandas()
    else:
        return e.to_xarray()

# Example: Fetch TAO/TRITON mooring SST
tao_sst = fetch_erddap(
    server_url="https://data.pmel.noaa.gov/pmel/erddap",
    dataset_id="pmel_tao_daily",
    variables=["station", "time", "latitude", "longitude", "T_25"],
    constraints={
        "time>=": "2024-01-01",
        "time<=": "2024-03-01",
        "latitude>=": -10,
        "latitude<=": 10,
        "longitude>=": 90,
        "longitude<=": 140,
    }
)
```

### 6.4 TAO/TRITON/RAMA Moorings

```python
def fetch_tao_mooring(lon, lat, depth_m=1, start="2024-01-01",
                      end="2024-03-01"):
    """Fetch TAO/TRITON/RAMA temperature data near a location.

    Uses PMEL ERDDAP. The array spans:
    - TAO: Tropical Pacific (137°E–95°W, 8°S–8°N)
    - TRITON: Western Pacific (137°E–156°E)
    - RAMA: Indian Ocean (55°E–100°E)
    """
    e = ERDDAP(
        server="https://data.pmel.noaa.gov/pmel/erddap",
        protocol="tabledap"
    )
    e.dataset_id = "pmel_tao_daily"
    e.variables = ["station", "time", "latitude", "longitude",
                   "T_25", "S_41", "U_320", "V_321"]
    e.constraints = {
        "time>=": start,
        "time<=": end,
        "latitude>=": lat - 1,
        "latitude<=": lat + 1,
        "longitude>=": lon - 1,
        "longitude<=": lon + 1,
    }

    try:
        df = e.to_pandas(parse_dates=["time (UTC)"])
        return df
    except Exception as ex:
        print(f"ERDDAP query failed: {ex}")
        return None
```

---

## 7. Reading Argo Float Data

### 7.1 Using argopy

`argopy` is the recommended Python library for Argo data access.

```python
import argopy
from argopy import DataFetcher

def fetch_argo_region(lon_range, lat_range, depth_range=(0, 2000),
                      start="2024-01-01", end="2024-03-01"):
    """Fetch Argo profiles in a region.

    Parameters
    ----------
    lon_range  : tuple — (lon_min, lon_max)
    lat_range  : tuple — (lat_min, lat_max)
    depth_range: tuple — (depth_min, depth_max) in meters
    start, end : str — date range

    Returns
    -------
    xr.Dataset with Argo profiles
    """
    loader = DataFetcher(src="gdac")  # or "erddap", "argovis"
    ds = loader.region(
        [lon_range[0], lon_range[1],
         lat_range[0], lat_range[1],
         depth_range[0], depth_range[1],
         start, end]
    ).to_xarray()

    return ds

# Example: Indonesian throughflow region
argo_itf = fetch_argo_region(
    lon_range=(115, 135),
    lat_range=(-10, 0),
    start="2024-01-01",
    end="2024-06-01"
)
```

### 7.2 Using Argo GDAC Directly

```python
def fetch_argo_profile(wmo_id, cycle, gdac="https://data-argo.ifremer.fr"):
    """Fetch a single Argo profile from the GDAC.

    Parameters
    ----------
    wmo_id : str — float WMO number (e.g., "2902696")
    cycle  : int — cycle number
    gdac   : str — GDAC base URL

    Returns
    -------
    xr.Dataset with the profile
    """
    # Argo file structure: dac/{dac_center}/{wmo_id}/profiles/R{wmo_id}_{cycle:03d}.nc
    # The DAC center varies — use the index file or argopy for reliable access
    url = f"{gdac}/dac/*/{wmo_id}/profiles/*{wmo_id}_{cycle:03d}.nc"
    # Better: use argopy's float accessor
    loader = DataFetcher(src="gdac").float(int(wmo_id))
    return loader.to_xarray()
```

### 7.3 Argo Data Processing

```python
import gsw  # TEOS-10 Gibbs SeaWater

def process_argo_profile(ds):
    """Process raw Argo profile: compute derived variables.

    Uses TEOS-10 (GSW) for accurate seawater properties.
    """
    # Extract core variables
    p = ds["PRES"].values       # pressure (dbar)
    t = ds["TEMP"].values       # in-situ temperature (°C)
    s = ds["PSAL"].values       # practical salinity (PSU)
    lat = float(ds["LATITUDE"].values)
    lon = float(ds["LONGITUDE"].values)

    # Convert to TEOS-10 conservative temperature and absolute salinity
    sa = gsw.SA_from_SP(s, p, lon, lat)          # absolute salinity
    ct = gsw.CT_from_t(sa, t, p)                 # conservative temperature
    pt = gsw.pt0_from_t(sa, t, p)                # potential temperature

    # Derived variables
    rho = gsw.rho(sa, ct, p)                     # in-situ density
    sigma0 = gsw.sigma0(sa, ct)                  # potential density anomaly
    n2, p_mid = gsw.Nsquared(sa, ct, p, lat=lat) # buoyancy frequency
    spice = gsw.spiciness0(sa, ct)               # spiciness

    # Mixed layer depth (density threshold method)
    ref_idx = np.argmin(np.abs(p - 10))  # reference at 10 dbar
    rho_ref = rho[ref_idx]
    mld_mask = rho > rho_ref + 0.03  # 0.03 kg/m³ threshold
    mld = p[mld_mask][0] if mld_mask.any() else np.nan

    return {
        "pressure": p, "temperature": ct, "salinity": sa,
        "potential_temp": pt, "density": rho, "sigma0": sigma0,
        "n2": n2, "mld": mld,
    }
```

---

## 8. Reading Tide Gauge Data

### 8.1 UHSLC (University of Hawaii Sea Level Center)

```python
def fetch_uhslc_hourly(station_id):
    """Fetch hourly sea level data from UHSLC.

    Parameters
    ----------
    station_id : int — UHSLC station number

    Returns
    -------
    xr.Dataset with sea level time series
    """
    url = (f"https://uhslc.soest.hawaii.edu/data/netcdf/"
           f"rqh/rqh{station_id:03d}a.nc")
    ds = xr.open_dataset(url)
    return ds

def fetch_uhslc_daily(station_id):
    """Fetch daily-mean sea level from UHSLC."""
    url = (f"https://uhslc.soest.hawaii.edu/data/netcdf/"
           f"rqd/rqd{station_id:03d}a.nc")
    ds = xr.open_dataset(url)
    return ds
```

### 8.2 IOC Sea Level Monitoring

```python
def fetch_ioc_sealevel(station_code, start, end):
    """Fetch sea level data from IOC/GLOSS via their API.

    Parameters
    ----------
    station_code : str — IOC station code (e.g., "jaka" for Jakarta)
    start, end   : str — "YYYY-MM-DD"

    Returns
    -------
    pd.DataFrame with time and sea_level_m
    """
    url = (f"https://www.ioc-sealevelmonitoring.org/service.php?"
           f"query=data&code={station_code}"
           f"&timestart={start}&timestop={end}&format=json")

    import requests
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    if not data:
        print(f"No data for station {station_code}")
        return pd.DataFrame()

    df = pd.DataFrame(data)
    df["time"] = pd.to_datetime(df["stime"])
    df["sea_level_m"] = pd.to_numeric(df["slevel"], errors="coerce") / 1000
    return df[["time", "sea_level_m"]].dropna()
```

### 8.3 Tidal Analysis with utide

```python
import utide

def tidal_analysis(time, sea_level, lat):
    """Perform tidal harmonic analysis using utide.

    Parameters
    ----------
    time      : array of datetime or float (days)
    sea_level : array of sea level values (m)
    lat       : float — station latitude

    Returns
    -------
    coef : utide Bunch — harmonic coefficients
    """
    coef = utide.solve(
        time, sea_level,
        lat=lat,
        method="ols",        # ordinary least squares
        conf_int="MC",       # Monte Carlo confidence intervals
        Rayleigh_min=0.95,   # minimum Rayleigh criterion
    )

    # Reconstruct tidal signal
    tide = utide.reconstruct(time, coef)

    # Non-tidal residual (surge + other signals)
    residual = sea_level - tide.h

    return coef, tide.h, residual


def print_tidal_constituents(coef, n_top=10):
    """Print the largest tidal constituents."""
    # Sort by amplitude
    idx = np.argsort(coef.A)[::-1]

    print(f"{'Constituent':<6} {'Amplitude (m)':>13} {'Phase (°)':>10}")
    print("-" * 32)
    for i in idx[:n_top]:
        name = coef.name[i]
        amp = coef.A[i]
        phase = coef.g[i]
        print(f"{name:<6} {amp:>13.4f} {phase:>10.1f}")
```

---

## 9. Reading Satellite Products

### 9.1 SST from CMEMS

```python
import copernicusmarine

def fetch_cmems_sst(lon_range, lat_range, start, end):
    """Fetch daily L4 SST from Copernicus Marine Service.

    Requires CMEMS credentials (set via copernicusmarine login).

    Parameters
    ----------
    lon_range : tuple — (lon_min, lon_max)
    lat_range : tuple — (lat_min, lat_max)
    start     : str — "YYYY-MM-DD"
    end       : str — "YYYY-MM-DD"

    Returns
    -------
    xr.Dataset with analysed_sst variable
    """
    ds = copernicusmarine.open_dataset(
        dataset_id="cmems_obs-sst_glo_phy-sst_nrt_diurnal-oi-0.01-deg_PT1H-m",
        variables=["analysed_sst"],
        minimum_longitude=lon_range[0],
        maximum_longitude=lon_range[1],
        minimum_latitude=lat_range[0],
        maximum_latitude=lat_range[1],
        start_datetime=f"{start}T00:00:00",
        end_datetime=f"{end}T23:59:59",
    )
    return ds
```

### 9.2 Satellite Altimetry (SWH + SSH)

```python
def fetch_cmems_altimetry(lon_range, lat_range, start, end):
    """Fetch gridded sea level anomaly from CMEMS."""
    ds = copernicusmarine.open_dataset(
        dataset_id="cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D",
        variables=["sla", "adt", "ugos", "vgos"],
        minimum_longitude=lon_range[0],
        maximum_longitude=lon_range[1],
        minimum_latitude=lat_range[0],
        maximum_latitude=lat_range[1],
        start_datetime=f"{start}T00:00:00",
        end_datetime=f"{end}T23:59:59",
    )
    return ds

def fetch_cmems_wave(lon_range, lat_range, start, end):
    """Fetch satellite-based significant wave height from CMEMS."""
    ds = copernicusmarine.open_dataset(
        dataset_id="cmems_obs-wave_glo_phy-swh_nrt_multi-l4-2deg_PT1H",
        variables=["VHMO"],  # SWH
        minimum_longitude=lon_range[0],
        maximum_longitude=lon_range[1],
        minimum_latitude=lat_range[0],
        maximum_latitude=lat_range[1],
        start_datetime=f"{start}T00:00:00",
        end_datetime=f"{end}T23:59:59",
    )
    return ds
```

### 9.3 Satellite SST from ERDDAP

```python
def fetch_sst_erddap(lon_range, lat_range, start, end):
    """Fetch daily MUR SST from CoastWatch ERDDAP."""
    e = ERDDAP(
        server="https://coastwatch.pfeg.noaa.gov/erddap",
        protocol="griddap"
    )
    e.dataset_id = "jplMURSST41"
    e.variables = ["analysed_sst"]
    e.constraints = {
        "time>=": start,
        "time<=": end,
        "latitude>=": lat_range[0],
        "latitude<=": lat_range[1],
        "longitude>=": lon_range[0],
        "longitude<=": lon_range[1],
    }

    ds = e.to_xarray()
    # Convert from K to °C if needed
    if ds["analysed_sst"].max() > 100:
        ds["analysed_sst"] = ds["analysed_sst"] - 273.15

    return ds
```

### 9.4 Scatterometer Winds

```python
def fetch_ascat_winds(lon_range, lat_range, start, end):
    """Fetch ASCAT wind vectors from CMEMS."""
    ds = copernicusmarine.open_dataset(
        dataset_id="cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H",
        variables=["eastward_wind", "northward_wind", "wind_speed",
                   "wind_stress"],
        minimum_longitude=lon_range[0],
        maximum_longitude=lon_range[1],
        minimum_latitude=lat_range[0],
        maximum_latitude=lat_range[1],
        start_datetime=f"{start}T00:00:00",
        end_datetime=f"{end}T23:59:59",
    )
    return ds
```

---

## 10. Quality Control

### 10.1 SST Quality Control

```python
def qc_sst(df, lat=None):
    """Quality control for SST observations.

    Parameters
    ----------
    df  : DataFrame with 'sst_c' column
    lat : float — station latitude (for setting tropical limits)
    """
    # Range check
    if lat is not None and abs(lat) < 30:
        # Tropical: SST typically 20–35°C
        sst_min, sst_max = 15.0, 36.0
    else:
        sst_min, sst_max = -2.0, 35.0

    df["sst_qc"] = 0  # good
    df.loc[df["sst_c"].isna(), "sst_qc"] = 9  # missing
    df.loc[(df["sst_c"] < sst_min) | (df["sst_c"] > sst_max), "sst_qc"] = 2

    # Spike check (max 2°C/hour for buoys)
    if isinstance(df.index, pd.DatetimeIndex):
        diff = df["sst_c"].diff().abs()
        dt_hours = df.index.to_series().diff().dt.total_seconds() / 3600
        rate = diff / dt_hours.replace(0, np.nan)
        df.loc[rate > 2.0, "sst_qc"] = np.maximum(
            df.loc[rate > 2.0, "sst_qc"], 1  # suspect
        )

    return df
```

### 10.2 Wave Data Quality Control

```python
def qc_waves(df):
    """Quality control for wave observations (buoy or altimeter).

    Expected columns: swh_m, peak_period_s, mean_period_s
    """
    # SWH range check
    if "swh_m" in df.columns:
        df["swh_qc"] = 0
        df.loc[df["swh_m"].isna(), "swh_qc"] = 9
        df.loc[(df["swh_m"] < 0) | (df["swh_m"] > 25), "swh_qc"] = 2

        # Spike check: max 3 m/hour
        if isinstance(df.index, pd.DatetimeIndex):
            diff = df["swh_m"].diff().abs()
            dt_h = df.index.to_series().diff().dt.total_seconds() / 3600
            rate = diff / dt_h.replace(0, np.nan)
            df.loc[rate > 3.0, "swh_qc"] = np.maximum(
                df.loc[rate > 3.0, "swh_qc"], 1
            )

    # Period range check
    for col in ["peak_period_s", "mean_period_s"]:
        if col in df.columns:
            qc_col = f"{col}_qc"
            df[qc_col] = 0
            df.loc[df[col].isna(), qc_col] = 9
            df.loc[(df[col] < 1) | (df[col] > 30), qc_col] = 2

    # Consistency: peak period should be >= mean period
    if "peak_period_s" in df.columns and "mean_period_s" in df.columns:
        mask = df["mean_period_s"] > df["peak_period_s"] * 1.5
        df.loc[mask, "mean_period_s_qc"] = np.maximum(
            df.loc[mask, "mean_period_s_qc"], 1
        )

    return df
```

### 10.3 Sea Level Quality Control

```python
def qc_sea_level(df, expected_range_m=5.0):
    """Quality control for tide gauge sea level data.

    Parameters
    ----------
    df              : DataFrame with 'sea_level_m' column (datetime index)
    expected_range_m : maximum expected tidal range (location-dependent)
    """
    df["sl_qc"] = 0
    df.loc[df["sea_level_m"].isna(), "sl_qc"] = 9

    # Remove mean to work with anomalies
    mean_sl = df["sea_level_m"].mean()
    anom = df["sea_level_m"] - mean_sl

    # Range check: values beyond expected tidal range + margin
    limit = expected_range_m * 1.5
    df.loc[anom.abs() > limit, "sl_qc"] = 2

    # Spike check: max 0.5 m change per minute (tsunami-scale)
    if isinstance(df.index, pd.DatetimeIndex):
        diff = df["sea_level_m"].diff().abs()
        dt_min = df.index.to_series().diff().dt.total_seconds() / 60
        rate = diff / dt_min.replace(0, np.nan)
        # Flag jumps > 0.5 m/min (instrument glitch, not physical)
        df.loc[rate > 0.5, "sl_qc"] = np.maximum(df.loc[rate > 0.5, "sl_qc"], 1)

    # Flat line check (stuck sensor): same value for > 3 hours
    is_flat = df["sea_level_m"].diff().abs() < 0.001
    flat_groups = is_flat.ne(is_flat.shift()).cumsum()
    flat_counts = flat_groups.map(flat_groups.value_counts())
    flat_duration_check = flat_counts > 180  # >3 hrs at 1-min sampling
    df.loc[flat_duration_check & is_flat, "sl_qc"] = 2

    return df
```

### 10.4 Argo Profile Quality Control

Argo data comes with QC flags from the data centers. Use them:

```python
ARGO_QC_FLAGS = {
    "0": "No QC",
    "1": "Good",
    "2": "Probably good",
    "3": "Probably bad (potentially correctable)",
    "4": "Bad",
    "5": "Changed (value adjusted)",
    "8": "Interpolated",
    "9": "Missing",
}

def filter_argo_qc(ds, max_qc="2"):
    """Filter Argo data to keep only good quality.

    Parameters
    ----------
    ds     : xr.Dataset from argopy
    max_qc : str — maximum acceptable QC flag ("1" = good only,
             "2" = good + probably good)

    Returns
    -------
    xr.Dataset with bad data masked
    """
    max_flag = int(max_qc)

    for var in ["TEMP", "PSAL", "PRES"]:
        qc_var = f"{var}_QC"
        if qc_var in ds:
            qc = ds[qc_var].astype(int)
            ds[var] = ds[var].where(qc <= max_flag)

    return ds
```

---

## 11. Preparing Obs for Model Verification

### 11.1 ROMS Verification

```python
from scipy.interpolate import RegularGridInterpolator

def extract_roms_at_obs(roms_ds, obs_lats, obs_lons, obs_depths=None,
                        var_name="temp"):
    """Extract ROMS output at observation locations.

    Parameters
    ----------
    roms_ds    : xr.Dataset — ROMS output (with lat_rho, lon_rho, s_rho)
    obs_lats   : array — observation latitudes
    obs_lons   : array — observation longitudes
    obs_depths : array — observation depths (m, positive down) or None
    var_name   : str — ROMS variable name

    Returns
    -------
    Array of model values at observation locations
    """
    # ROMS uses curvilinear grids — need 2D interpolation
    mlat = roms_ds["lat_rho"].values
    mlon = roms_ds["lon_rho"].values
    data = roms_ds[var_name].values.squeeze()

    if data.ndim == 2:
        # Surface field
        tree = cKDTree(np.column_stack([mlat.ravel(), mlon.ravel()]))
        _, idx = tree.query(np.column_stack([obs_lats, obs_lons]))
        return data.ravel()[idx]

    elif data.ndim == 3 and obs_depths is not None:
        # 3D field — need depth interpolation
        # ROMS s-coordinate → actual depth requires h and zeta
        h = roms_ds["h"].values
        zeta = roms_ds["zeta"].values.squeeze() if "zeta" in roms_ds else 0
        s = roms_ds["s_rho"].values
        Cs = roms_ds["Cs_r"].values
        hc = float(roms_ds["hc"].values)

        # Extract profiles at observation locations
        tree = cKDTree(np.column_stack([mlat.ravel(), mlon.ravel()]))
        _, idx = tree.query(np.column_stack([obs_lats, obs_lons]))

        results = []
        for i, oi in enumerate(idx):
            j, ii = np.unravel_index(oi, mlat.shape)
            profile = data[:, j, ii]
            h_local = h[j, ii]
            z_local = zeta_to_z(s, Cs, hc, h_local,
                                zeta[j, ii] if np.ndim(zeta) > 0 else zeta)
            # Interpolate to observation depth
            from scipy.interpolate import interp1d
            f = interp1d(-z_local, profile, bounds_error=False,
                        fill_value=np.nan)
            results.append(f(obs_depths[i]))

        return np.array(results)


def zeta_to_z(s, Cs, hc, h, zeta=0):
    """Convert ROMS s-coordinate to depth (z, negative down).

    Parameters
    ----------
    s    : array — s-coordinate values [-1, 0]
    Cs   : array — stretching function values
    hc   : float — critical depth
    h    : float — bottom depth at this point
    zeta : float — sea surface height

    Returns
    -------
    z : array — depth values (negative down)
    """
    z = (hc * s + h * Cs) / (hc + h)
    z = zeta + (zeta + h) * z
    return z
```

### 11.2 SWAN / WW3 Wave Verification

```python
def extract_wave_model_at_buoy(model_ds, buoy_lat, buoy_lon,
                                hs_var="hs", tp_var="tp", dir_var="dir"):
    """Extract wave model output at a buoy location.

    Works for both SWAN and WW3 output (variable names may differ).

    SWAN typical: Hsig, TPsmoo, Dir
    WW3 typical:  hs, tp, dp
    """
    mlat = model_ds["latitude"].values
    mlon = model_ds["longitude"].values

    if mlat.ndim == 1:
        # Regular grid — use nearest
        lat_idx = np.argmin(np.abs(mlat - buoy_lat))
        lon_idx = np.argmin(np.abs(mlon - buoy_lon))

        result = pd.DataFrame()
        result["time"] = model_ds["time"].values
        result["model_hs"] = model_ds[hs_var][:, lat_idx, lon_idx].values

        if tp_var in model_ds:
            result["model_tp"] = model_ds[tp_var][:, lat_idx, lon_idx].values
        if dir_var in model_ds:
            result["model_dir"] = model_ds[dir_var][:, lat_idx, lon_idx].values

        return result

    else:
        # Curvilinear/unstructured — use KD-tree
        tree = cKDTree(np.column_stack([mlat.ravel(), mlon.ravel()]))
        _, idx = tree.query([buoy_lat, buoy_lon])

        if mlat.ndim == 2:
            j, i = np.unravel_index(idx, mlat.shape)
            result = pd.DataFrame()
            result["time"] = model_ds["time"].values
            result["model_hs"] = model_ds[hs_var][:, j, i].values
            return result


def wave_verification_metrics(obs_hs, model_hs):
    """Compute standard wave verification metrics.

    See Python Verification Guide Section 12 for wave-specific methods.

    Returns
    -------
    dict with ME, MAE, RMSE, scatter_index, r, skill_score
    """
    mask = ~np.isnan(obs_hs) & ~np.isnan(model_hs)
    o = obs_hs[mask]
    m = model_hs[mask]

    me = np.mean(m - o)
    mae = np.mean(np.abs(m - o))
    rmse = np.sqrt(np.mean((m - o)**2))
    r = np.corrcoef(m, o)[0, 1] if len(o) > 2 else np.nan

    # Scatter index (normalized RMSE, standard for waves)
    si = rmse / np.mean(o) if np.mean(o) > 0 else np.nan

    # Skill score (vs persistence or climatology)
    ss = 1 - np.sum((m - o)**2) / np.sum((o - np.mean(o))**2)

    return {
        "N": len(o), "ME": me, "MAE": mae, "RMSE": rmse,
        "SI": si, "r": r, "skill_score": ss,
    }
```

### 11.3 SSH / Sea Level Verification

```python
def sea_level_verification(obs_sl, model_sl, obs_time, model_time,
                           do_detide=True, lat=None):
    """Compare tide gauge sea level with ROMS output.

    Parameters
    ----------
    obs_sl     : array — observed sea level (m)
    model_sl   : array — modeled sea level / zeta (m)
    obs_time   : array — observation times
    model_time : array — model output times
    do_detide  : bool — remove tidal signal before comparing
    lat        : float — station latitude (for tidal analysis)
    """
    # Interpolate model to observation times
    from scipy.interpolate import interp1d

    # Convert times to float for interpolation
    t0 = min(obs_time[0], model_time[0])
    obs_hours = (pd.DatetimeIndex(obs_time) - t0).total_seconds() / 3600
    mod_hours = (pd.DatetimeIndex(model_time) - t0).total_seconds() / 3600

    f = interp1d(mod_hours, model_sl, bounds_error=False, fill_value=np.nan)
    model_interp = f(obs_hours)

    if do_detide and lat is not None:
        # Remove tides from both signals
        _, obs_tide, obs_resid = tidal_analysis(obs_hours / 24, obs_sl, lat)
        _, mod_tide, mod_resid = tidal_analysis(obs_hours / 24, model_interp, lat)

        # Compare: total, tidal, and non-tidal components
        results = {}
        for label, o, m in [("total", obs_sl, model_interp),
                             ("tidal", obs_tide, mod_tide),
                             ("non_tidal", obs_resid, mod_resid)]:
            mask = ~np.isnan(o) & ~np.isnan(m)
            results[label] = {
                "RMSE": np.sqrt(np.mean((m[mask] - o[mask])**2)),
                "ME":   np.mean(m[mask] - o[mask]),
                "r":    np.corrcoef(m[mask], o[mask])[0, 1],
            }
        return results
    else:
        mask = ~np.isnan(obs_sl) & ~np.isnan(model_interp)
        o, m = obs_sl[mask], model_interp[mask]
        return {
            "total": {
                "RMSE": np.sqrt(np.mean((m - o)**2)),
                "ME":   np.mean(m - o),
                "r":    np.corrcoef(m, o)[0, 1],
            }
        }
```

### 11.4 SST Verification

```python
def sst_verification(obs_sst, sat_sst_ds, obs_lats, obs_lons, obs_times):
    """Compare buoy SST against satellite SST product.

    Useful for:
    1. Validating satellite SST before using it to verify ROMS
    2. Direct buoy vs ROMS SST comparison

    Parameters
    ----------
    obs_sst    : array — buoy SST values (°C)
    sat_sst_ds : xr.Dataset — satellite SST (e.g., MUR, OSTIA)
    obs_lats   : array — buoy latitudes
    obs_lons   : array — buoy longitudes
    obs_times  : array — buoy observation times

    Returns
    -------
    dict of verification metrics
    """
    sat_values = []
    for t, lat, lon in zip(obs_times, obs_lats, obs_lons):
        try:
            val = sat_sst_ds["analysed_sst"].sel(
                time=t, latitude=lat, longitude=lon,
                method="nearest",
                tolerance=pd.Timedelta("1D")
            ).values
            # Convert K → °C if needed
            if val > 100:
                val -= 273.15
            sat_values.append(float(val))
        except (KeyError, IndexError):
            sat_values.append(np.nan)

    sat_values = np.array(sat_values)
    mask = ~np.isnan(obs_sst) & ~np.isnan(sat_values)

    o = obs_sst[mask]
    s = sat_values[mask]

    return {
        "N": len(o),
        "ME": np.mean(s - o),
        "MAE": np.mean(np.abs(s - o)),
        "RMSE": np.sqrt(np.mean((s - o)**2)),
        "r": np.corrcoef(s, o)[0, 1] if len(o) > 2 else np.nan,
    }
```

---

## 12. Indonesian & Tropical Waters

### 12.1 Key Oceanographic Features

Indonesia sits at the crossroads of the Pacific and Indian Oceans, with unique oceanographic features:

| Feature | Description | Relevance |
|---------|-------------|-----------|
| **Indonesian Throughflow (ITF)** | Pacific-to-Indian ocean transport through Indonesian straits | Affects SST, salinity, ROMS validation |
| **Sunda Strait** | Narrow strait between Java and Sumatra, strong tidal currents | COAWST/ROMS application area |
| **Makassar Strait** | Primary ITF pathway, ~80% of transport | Key validation region |
| **Java Sea** | Shallow shelf sea, strong tidal mixing | Wave model validation |
| **Banda Sea** | Deep basin, upwelling | Subsurface T/S validation |
| **Indo-Pacific Warm Pool** | SST consistently >28°C | SST verification baseline |
| **MJO influence** | Intraseasonal wind/wave variability | Affects verification timing |

### 12.2 Available Observation Networks

| Platform | Coverage | Access |
|----------|----------|--------|
| **RAMA moorings** | Eastern Indian Ocean (90°E–100°E, 15°S–15°N) | PMEL ERDDAP |
| **TRITON moorings** | Western Pacific (137°E–156°E) | PMEL ERDDAP |
| **Argo floats** | Indonesian seas (limited in shallow regions) | `argopy` |
| **BMKG tide gauges** | ~150 stations across archipelago | BMKG / IOC |
| **INA-CBT** | Real-time tsunami buoys | BMKG |
| **BIG tide gauges** | Geodetic sea level stations | BIG |
| **Satellite altimetry** | Full coverage (clouds not an issue) | CMEMS |
| **Satellite SST** | Full coverage (frequent cloud gaps in tropics) | CMEMS / CoastWatch |

### 12.3 Indonesian Tide Gauge Stations

| Station | Code | Lat | Lon | Data Source |
|---------|------|-----|-----|-------------|
| Tanjung Priok (Jakarta) | jaka | -6.10 | 106.87 | IOC/UHSLC |
| Surabaya | sura | -7.20 | 112.73 | IOC |
| Benoa (Bali) | beno | -8.74 | 115.21 | IOC/UHSLC |
| Padang | pada | -0.95 | 100.37 | IOC |
| Cilacap | cila | -7.75 | 109.02 | IOC/UHSLC |
| Makassar | maka | -5.11 | 119.42 | IOC |
| Ambon | ambo | -3.69 | 128.18 | IOC |
| Jayapura | jaya | -2.53 | 140.72 | IOC |
| Sabang | saba | 5.89 | 95.32 | IOC |
| Merauke | mera | -8.47 | 140.33 | IOC |

### 12.4 ITF Monitoring

```python
def fetch_itf_data():
    """Fetch data relevant to Indonesian Throughflow monitoring.

    Key moorings and their approximate locations:
    - INSTANT (2004-2006): Makassar, Lombok, Ombai, Timor straits
    - RAMA 0°, 90°E and 5°S, 95°E: near ITF exit
    """
    # RAMA moorings near the ITF exit
    rama_data = fetch_tao_mooring(lon=95, lat=-5)

    # Argo floats in the ITF region
    argo_itf = fetch_argo_region(
        lon_range=(115, 135),
        lat_range=(-12, 2),
        start="2024-01-01",
        end="2024-06-01"
    )

    return {"rama": rama_data, "argo": argo_itf}
```

### 12.5 Sunda Strait Observations

For COAWST validation in the Sunda Strait area (see [COAWST Guide](../models/coawst-model-guide.md)):

```python
def sunda_strait_obs():
    """Collect available observations for Sunda Strait verification.

    The Sunda Strait is a narrow (~24 km) channel between Java and
    Sumatra with strong tidal currents and complex bathymetry.
    """
    # 1. Nearest tide gauges
    stations = {
        "Merak":    {"lat": -5.93, "lon": 106.00, "code": "mera"},
        "Panjang":  {"lat": -5.46, "lon": 105.32, "code": None},  # BMKG
        "Ciwandan": {"lat": -6.02, "lon": 106.05, "code": None},  # BMKG
    }

    # 2. Satellite SST in the strait
    sst = fetch_sst_erddap(
        lon_range=(104.5, 106.5),
        lat_range=(-7.0, -5.5),
        start="2024-01-01",
        end="2024-02-01"
    )

    # 3. Satellite altimetry (limited in narrow straits)
    # Note: Altimetry accuracy degrades near coasts (<20 km)

    return {"tide_stations": stations, "sst": sst}
```

### 12.6 Tropical Ocean QC Adjustments

```python
TROPICAL_OCEAN_LIMITS = {
    "sst_c":          (20.0, 35.0),    # Indo-Pacific warm pool
    "salinity_psu":   (28.0, 36.0),    # Fresh river influence
    "swh_m":          (0.0,  10.0),    # Lower than mid-latitudes
    "peak_period_s":  (2.0,  18.0),    # Mixed swell + wind-sea
    "sea_level_m":    (-3.0,  3.0),    # Relative to mean
    "current_ms":     (0.0,   3.0),    # ITF can be strong
}
```

---

## 13. End-to-End Example

### Workflow: Verify ROMS SST Against Satellite + Buoy Data

```
[1. Download satellite SST (CMEMS)]
        |
[2. Download buoy SST (ERDDAP)]
        |
[3. QC buoy data]
        |
[4. Validate satellite vs buoy]     <-- check satellite product quality
        |
[5. Load ROMS output]
        |
[6. Extract ROMS at buoy locations]
        |
[7. Extract ROMS on satellite grid]
        |
[8. Compute metrics & visualize]
```

### Complete Script

```python
"""
End-to-end ocean verification:
  Satellite SST + buoy SST → ROMS SST verification
  Region: Eastern Indian Ocean / Indonesian waters
"""
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ── Step 1: Define region ──────────────────────────────────────

REGION = {
    "lon_min": 95, "lon_max": 125,
    "lat_min": -15, "lat_max": 5,
}
START = "2024-01-15"
END   = "2024-01-31"

# ── Step 2: Fetch satellite SST ────────────────────────────────

# Option A: CMEMS (requires login)
# sat_sst = fetch_cmems_sst(
#     (REGION["lon_min"], REGION["lon_max"]),
#     (REGION["lat_min"], REGION["lat_max"]),
#     START, END
# )

# Option B: ERDDAP (open access)
sat_sst = fetch_sst_erddap(
    lon_range=(REGION["lon_min"], REGION["lon_max"]),
    lat_range=(REGION["lat_min"], REGION["lat_max"]),
    start=START, end=END
)
print(f"Satellite SST: {sat_sst.dims}")

# ── Step 3: Fetch buoy/mooring SST ────────────────────────────

# RAMA mooring near 0°N, 90°E
rama = fetch_tao_mooring(lon=90, lat=0, start=START, end=END)

# Also try NDBC buoys if any in region
# ndbc = read_ndbc_stdmet("52001", 2024)  # if available

# ── Step 4: QC buoy data ──────────────────────────────────────

if rama is not None and len(rama) > 0:
    # Standardize column names (ERDDAP returns verbose names)
    rama = rama.rename(columns=lambda c: c.strip().lower())
    if "t_25" in rama.columns:
        rama["sst_c"] = pd.to_numeric(rama["t_25"], errors="coerce")
    rama = qc_sst(rama, lat=0)
    # Keep only good data
    rama.loc[rama["sst_qc"] >= 2, "sst_c"] = np.nan
    rama_clean = rama.dropna(subset=["sst_c"])
    print(f"RAMA observations after QC: {len(rama_clean)}")

# ── Step 5: Load ROMS output ──────────────────────────────────

roms = xr.open_dataset("roms_his.nc")
roms_sst = roms["temp"].isel(s_rho=-1)  # surface layer
roms_lat = roms["lat_rho"].values
roms_lon = roms["lon_rho"].values
roms_time = pd.DatetimeIndex(roms["ocean_time"].values)

# ── Step 6: Extract ROMS at buoy location ─────────────────────

if rama is not None and len(rama_clean) > 0:
    buoy_lat = float(rama_clean["latitude"].iloc[0])
    buoy_lon = float(rama_clean["longitude"].iloc[0])

    tree = cKDTree(np.column_stack([roms_lat.ravel(), roms_lon.ravel()]))
    _, idx = tree.query([buoy_lat, buoy_lon])
    j, i = np.unravel_index(idx, roms_lat.shape)

    roms_sst_buoy = roms_sst[:, j, i].values
    print(f"ROMS SST at buoy: {len(roms_sst_buoy)} time steps")

    # Match times
    paired = []
    for _, row in rama_clean.iterrows():
        obs_time = pd.Timestamp(row["time (utc)"])
        diffs = abs(roms_time - obs_time)
        best = diffs.argmin()
        if diffs[best] <= pd.Timedelta("3h"):
            paired.append({
                "time": obs_time,
                "obs_sst": row["sst_c"],
                "model_sst": float(roms_sst_buoy[best]),
            })

    paired_df = pd.DataFrame(paired)
    print(f"Paired buoy-ROMS: {len(paired_df)} points")

# ── Step 7: Compute metrics ───────────────────────────────────

if len(paired_df) > 5:
    obs = paired_df["obs_sst"].values
    mod = paired_df["model_sst"].values

    metrics = {
        "N":    len(obs),
        "ME":   np.mean(mod - obs),
        "MAE":  np.mean(np.abs(mod - obs)),
        "RMSE": np.sqrt(np.mean((mod - obs)**2)),
        "r":    np.corrcoef(mod, obs)[0, 1],
    }

    print("\n=== ROMS SST Verification (vs RAMA buoy) ===")
    for k, v in metrics.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

# ── Step 8: Visualize ─────────────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (a) Spatial SST map (satellite)
ax = fig.add_subplot(1, 3, 1, projection=ccrs.PlateCarree())
sat_mean = sat_sst["analysed_sst"].mean(dim="time")
sat_mean.plot(ax=ax, transform=ccrs.PlateCarree(),
              cmap="RdYlBu_r", vmin=26, vmax=32,
              cbar_kwargs={"label": "SST (°C)"})
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor="lightgray")
if rama is not None:
    ax.plot(buoy_lon, buoy_lat, "k^", markersize=10,
            transform=ccrs.PlateCarree(), label="RAMA buoy")
    ax.legend()
ax.set_title("Mean Satellite SST")
ax.set_extent([REGION["lon_min"], REGION["lon_max"],
               REGION["lat_min"], REGION["lat_max"]])

# (b) Time series
ax = axes[1]
if len(paired_df) > 0:
    ax.plot(paired_df["time"], paired_df["obs_sst"],
            "ko-", markersize=3, label="RAMA buoy")
    ax.plot(paired_df["time"], paired_df["model_sst"],
            "r.-", markersize=3, label="ROMS")
    ax.set_ylabel("SST (°C)")
    ax.set_title("SST Time Series")
    ax.legend()
    ax.tick_params(axis="x", rotation=45)

# (c) Scatter
ax = axes[2]
if len(paired_df) > 0:
    ax.scatter(paired_df["obs_sst"], paired_df["model_sst"],
               s=15, alpha=0.6)
    vmin = min(obs.min(), mod.min()) - 0.5
    vmax = max(obs.max(), mod.max()) + 0.5
    ax.plot([vmin, vmax], [vmin, vmax], "k--", linewidth=0.5)
    ax.set_xlabel("Observed SST (°C)")
    ax.set_ylabel("ROMS SST (°C)")
    ax.set_title(f"RMSE={metrics['RMSE']:.3f}°C, r={metrics['r']:.3f}")
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_aspect("equal")

plt.tight_layout()
plt.savefig("ocean_sst_verification.png", dpi=150)
plt.show()
```

---

## Appendix A: Variable Reference

### Standard Ocean Variables

| Variable | Symbol | Unit | Typical Range | In-Situ Source | Satellite Source |
|----------|--------|------|---------------|----------------|-----------------|
| Sea surface temperature | SST | °C | -2 to 35 | Buoy, drifter, ship | IR/MW radiometer |
| Subsurface temperature | T | °C | -2 to 30 | Argo, mooring | — |
| Salinity | S | PSU | 0 to 40 | Argo, mooring | SMOS, SMAP |
| Significant wave height | Hs | m | 0 to 20+ | Wave buoy | Altimeter, SAR |
| Peak wave period | Tp | s | 1 to 25 | Wave buoy | SAR (limited) |
| Mean wave period | Tm | s | 1 to 15 | Wave buoy | — |
| Wave direction | θ | ° | 0 to 360 | Directional buoy | SAR |
| Sea level | η | m | ±5 | Tide gauge | Altimeter |
| Sea surface height anomaly | SLA | m | ±1 | — | Altimeter |
| Surface current (U) | u | m/s | ±3 | Drifter, ADCP, HF radar | Altimeter (geostrophic) |
| Surface current (V) | v | m/s | ±3 | Drifter, ADCP, HF radar | Altimeter (geostrophic) |
| Ocean surface wind speed | U₁₀ | m/s | 0 to 40+ | Buoy anemometer | Scatterometer |

### Wave Spectral Parameters

| Parameter | Symbol | Definition | Unit |
|-----------|--------|------------|------|
| Spectral energy density | S(f) | Energy per frequency band | m²/Hz |
| Directional spectrum | S(f,θ) | Energy per frequency and direction | m²/Hz/rad |
| Zeroth moment | m₀ | ∫S(f)df | m² |
| Significant wave height | Hs = Hm0 | 4√m₀ | m |
| Mean period (Tm01) | Tm01 | m₀/m₁ | s |
| Mean period (Tm02) | Tm02 | √(m₀/m₂) | s |
| Peak period | Tp | 1/f_peak | s |
| Peak direction | θp | direction at peak frequency | ° |

---

## Appendix B: Resources & Links

### Documentation

| Resource | URL |
|----------|-----|
| Argo Data Management | `argo.ucsd.edu/data` |
| argopy Documentation | `argopy.readthedocs.io` |
| erddapy Documentation | `ioos.github.io/erddapy` |
| CMEMS Documentation | `help.marine.copernicus.eu` |
| GSW (TEOS-10) | `teos-10.github.io/GSW-Python` |
| utide Documentation | `github.com/wesleybowman/UTide` |
| wavespectra | `wavespectra.readthedocs.io` |
| NDBC Data Guide | `ndbc.noaa.gov/docs` |
| UHSLC Data | `uhslc.soest.hawaii.edu` |
| IOC Sea Level | `ioc-sealevelmonitoring.org` |

### Companion Guides in This Repo

| Guide | Relevance |
|-------|-----------|
| [Surface Observations Guide](surface-obs-guide.md) | Land-based met stations, SYNOP, METAR |
| [Python Verification Guide](../verification/python-verification-guide.md) | Metrics implementation (use with ocean obs) |
| [MET Verification Guide](../verification/met-verification-guide.md) | MET tools for point/grid verification |
| [ROMS Guide](../models/roms-model-guide.md) | Ocean model that these obs verify |
| [SWAN Guide](../models/swan-model-guide.md) | Wave model (buoy/altimeter verification) |
| [WW3 Guide](../models/ww3-model-guide.md) | Wave model verification |
| [COAWST Guide](../models/coawst-model-guide.md) | Coupled system verification |

### Quick Reference: Which Data for Which Model

| Model | Variable | Best In-Situ Source | Best Satellite Source |
|-------|----------|--------------------|-----------------------|
| ROMS SST | Temperature | Argo, mooring, drifter | MUR SST, OSTIA |
| ROMS T/S profile | Temperature, salinity | Argo | — |
| ROMS SSH | Sea level | Tide gauge | Altimetry (CMEMS) |
| ROMS currents | U, V | Drifter, ADCP, HF radar | Altimetry (geostrophic) |
| SWAN/WW3 Hs | Significant wave height | Wave buoy (NDBC) | Altimetry (CMEMS) |
| SWAN/WW3 Tp | Peak period | Wave buoy | SAR (limited) |
| SWAN/WW3 dir | Wave direction | Directional buoy | SAR |
| WRF surface wind | U₁₀, V₁₀ | Buoy anemometer | Scatterometer (ASCAT) |
