# Surface Observations: Complete Guide

> A practical guide to surface meteorological observations — formats, data sources,
> quality control, and Python-based access. Covers SYNOP, METAR, and automated
> station data with a focus on preparing observations for NWP verification.
>
> **Companion guides:**
> - [Forecast Verification: Methods and Key Concepts](../verification/forecast-verification-notes.md)
> - [Python Forecast Verification Guide](../verification/python-verification-guide.md) — uses the observation data prepared here
> - [MET Verification Guide](../verification/met-verification-guide.md) — MET ASCII obs format covered in Section 10
> - [Ocean Observations Guide](ocean-obs-guide.md)

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Surface Observation Types](#2-surface-observation-types)
3. [Data Formats](#3-data-formats)
4. [Data Sources](#4-data-sources)
5. [Python Environment Setup](#5-python-environment-setup)
6. [Reading & Decoding Surface Data](#6-reading--decoding-surface-data)
7. [Quality Control](#7-quality-control)
8. [Data Processing](#8-data-processing)
9. [Station Metadata & Network Selection](#9-station-metadata--network-selection)
10. [Preparing Obs for Verification](#10-preparing-obs-for-verification)
11. [Indonesian & Tropical Context](#11-indonesian--tropical-context)
12. [End-to-End Example](#12-end-to-end-example)

**Appendices:**
- [A. Variable Code Tables](#appendix-a-variable-code-tables)
- [B. Resources & Links](#appendix-b-resources--links)

---

## 1. Introduction

Surface observations are the backbone of NWP verification and data assimilation. They provide ground-truth measurements of temperature, humidity, wind, pressure, precipitation, and visibility at specific locations and times.

### What this guide covers

| Topic | Description |
|-------|-------------|
| Observation types | SYNOP, METAR, ASOS/AWS, and their reporting conventions |
| Formats | TAC (alphanumeric), BUFR, ISD, CSV, PrepBUFR |
| Data access | Free global and regional data sources with Python code |
| Quality control | Automated QC checks you can implement in Python |
| Verification prep | Converting obs to formats usable by MET and Python verification |

### What this guide does NOT cover

- Upper-air observations (radiosondes, pilot balloons) — planned for a separate guide
- Radar observations — see precipitation verification in the [Python Verification Guide](../verification/python-verification-guide.md)
- Satellite retrievals — see the [Ocean Observations Guide](ocean-obs-guide.md) for satellite SST/altimetry
- Data assimilation cycling — see the [WRF Guide](../models/wrf-model-guide.md) for WRF-DA

---

## 2. Surface Observation Types

### 2.1 SYNOP (FM-12)

The **SYNOPtic** report is the WMO standard for surface weather observations from land stations.

| Property | Detail |
|----------|--------|
| Full name | FM-12 Surface Synoptic Observation |
| Reporting times | Main: 00, 06, 12, 18 UTC; Intermediate: 03, 09, 15, 21 UTC |
| Station ID | 5-digit WMO number (Block + Station), e.g., `96749` = Block 96 (Indonesia), Station 749 |
| Spatial coverage | ~11,000 stations globally on the GTS |
| Variables | MSLP, station pressure, T, Td, wind (dd/ff), weather, cloud, visibility, precip, etc. |

**TAC (Traditional Alphanumeric Code) example:**

```
AAXX 15124
96749 42998 10502 10268 20231 30056 40117 52012 60001 70200 84500
333 10285 20218 55300 69942 80000
```

Decoding the header:
- `AAXX` = land SYNOP
- `15124` = day 15, 12 UTC, wind in m/s (indicator 4)
- `96749` = WMO station number

Key groups:
- `iihVV` (42998) — cloud height, visibility
- `Nddff` (10502) — total cloud, wind direction 050°, speed 02 m/s
- `1sTTT` (10268) — temperature +26.8°C
- `2sTTT` (20231) — dewpoint +23.1°C
- `3PPPP` (30056) — station pressure 1005.6 hPa
- `4PPPP` (40117) — MSL pressure 1011.7 hPa
- `6RRRt` (60001) — precipitation

### 2.2 METAR (FM-15/16)

**METeorological Aerodrome Report** — aviation weather observations from airports.

| Property | Detail |
|----------|--------|
| Full name | FM-15 (routine) / FM-16 (SPECI — special) |
| Reporting times | Hourly or half-hourly (varies by country) |
| Station ID | 4-letter ICAO code, e.g., `WIII` = Soekarno-Hatta (Jakarta) |
| Spatial coverage | ~10,000 airports globally |
| Variables | Wind, visibility, weather phenomena, cloud, T/Td, QNH, RVR |

**Example:**

```
METAR WIII 150700Z 32005KT 9999 FEW020CB SCT025 31/24 Q1010 NOSIG
```

Decoding:
- `WIII` — ICAO station
- `150700Z` — day 15, 07:00 UTC
- `32005KT` — wind from 320° at 5 knots
- `9999` — visibility ≥ 10 km
- `FEW020CB` — few cumulonimbus at 2000 ft
- `SCT025` — scattered clouds at 2500 ft
- `31/24` — temperature 31°C, dewpoint 24°C
- `Q1010` — QNH 1010 hPa

### 2.3 ASOS / AWS (Automated Stations)

| Property | Detail |
|----------|--------|
| Type | Automated Surface Observing System (ASOS) / Automatic Weather Station (AWS) |
| Reporting interval | 1-minute to hourly (varies) |
| Variables | Same as SYNOP + additional (solar radiation, soil temperature, etc.) |
| Format | Varies by agency — CSV, binary, proprietary |
| Examples | US ASOS (1-min), BMKG AWS, JMA AMeDAS, BoM AWS |

Key difference from SYNOP/METAR: AWS data is typically continuous (sub-hourly) and stored in agency-specific formats, not on the GTS.

### 2.4 Comparison Table

| Feature | SYNOP | METAR | AWS |
|---------|-------|-------|-----|
| Primary use | Climatology, NWP | Aviation | Real-time monitoring |
| Reporting freq | 3–6 hourly | 1 hourly | 1–60 min |
| Station ID | WMO 5-digit | ICAO 4-letter | Agency-specific |
| On GTS | Yes | Yes | Partial |
| Human observer | Often | Sometimes | No |
| Precipitation | 6/12/24-hr accum | No (remarks only) | Tipping bucket (hourly) |
| Cloud obs | Detailed (Nh, Cl, Cm, Ch) | Layer-based (FEW/SCT/BKN/OVC) | Ceilometer only |
| Wind averaging | 10-min mean | 2-min mean (US) / 10-min (ICAO) | Configurable |

---

## 3. Data Formats

### 3.1 TAC (Traditional Alphanumeric Code)

The legacy text-based format. Still widely used but WMO is migrating to BUFR.

```
AAXX 15064
72493 41583 83210 10156 20089 30180 40210 52014 69921 72562 84501
333 10189 20072 55300
```

- Fixed-width groups of 5 digits
- Position-encoded (group 1 = pressure, group 2 = temperature, etc.)
- Must decode using WMO Manual on Codes (WMO-306)

### 3.2 BUFR (Binary Universal Form for the Representation of meteorological data)

WMO's modern binary format. Self-describing, table-driven.

| Property | Detail |
|----------|--------|
| Structure | Sections 0–5 (Indicator, Identification, Data Description, Data, End) |
| Descriptors | Tables B (elements), C (operators), D (sequences) |
| SYNOP template | `TM 307080` — WMO standard SYNOP BUFR sequence |
| METAR template | `TM 307021` — WMO standard METAR BUFR sequence |
| Key advantage | Supports associated quality flags per element |

### 3.3 ISD (Integrated Surface Database)

NOAA's unified archive — merges SYNOP, METAR, and other surface reports into a single format.

| Property | Detail |
|----------|--------|
| Provider | NOAA NCEI |
| Coverage | 1901–present, ~35,000 stations |
| Format | Fixed-width ASCII (control + mandatory + additional data sections) |
| Access | FTP, AWS S3, or via Python (see Section 6) |
| File naming | `{USAF}-{WBAN}-{YYYY}.gz` |

ISD record layout (mandatory section):

| Column | Width | Field |
|--------|-------|-------|
| 1–4 | 4 | Total chars |
| 5–10 | 6 | USAF ID |
| 11–15 | 5 | WBAN ID |
| 16–23 | 8 | Date (YYYYMMDD) |
| 24–27 | 4 | Time (HHMM UTC) |
| 28 | 1 | Data source flag |
| 29–34 | 6 | Latitude (×1000) |
| 35–41 | 7 | Longitude (×1000) |
| 46–51 | 6 | Elevation (m) |
| 60–63 | 4 | Wind direction (degrees) |
| 64 | 1 | Wind direction quality |
| 65 | 1 | Wind observation type |
| 66–69 | 4 | Wind speed (×10 m/s) |
| 70 | 1 | Wind speed quality |
| 88–92 | 5 | Air temperature (×10 °C) |
| 93 | 1 | Temperature quality |
| 94–98 | 5 | Dewpoint (×10 °C) |
| 99 | 1 | Dewpoint quality |
| 100–104 | 5 | Sea level pressure (×10 hPa) |
| 105 | 1 | SLP quality |

### 3.4 CSV / ASCII (Agency-Specific)

Many national met services provide data in CSV or flat-file formats:

```csv
station_id,datetime,temp_c,rh_pct,wind_dir,wind_spd_ms,slp_hpa,precip_mm
96749,2024-01-15 06:00,26.8,82,50,2.0,1011.7,0.0
96749,2024-01-15 12:00,31.2,65,120,3.5,1009.4,0.0
```

### 3.5 PrepBUFR (NCEP Prepared BUFR)

NCEP's quality-controlled BUFR files used for data assimilation.

| Property | Detail |
|----------|--------|
| Provider | NCEP / GDAS |
| Content | Surface, upper-air, satellite — all QC'd and ready for DA |
| Tools | `prepbufr_decode` (Fortran), GSI utilities, or Python `ncepbufr` |
| Use case | Input for WRF-DA, or extracting QC'd obs for verification |

---

## 4. Data Sources

### 4.1 Global Sources

| Source | Coverage | Format | Access | URL |
|--------|----------|--------|--------|-----|
| **ISD (NCEI)** | Global, 1901–now | ISD ASCII | FTP / S3 / API | `ncei.noaa.gov/products/land-based-station/integrated-surface-database` |
| **MADIS** | US + some global | NetCDF, XML | API (registration) | `madis-data.ncep.noaa.gov` |
| **IEM (Iowa)** | US ASOS/AWOS | CSV, JSON | API (free) | `mesonet.agron.iastate.edu` |
| **OGIMET** | Global SYNOP/METAR | TAC text | Web scraping | `ogimet.com` |
| **WMO WIS** | Global GTS | BUFR | API (WIS 2.0) | `wis2.wmo.int` |
| **ECMWF MARS** | Global | BUFR | API (registered) | `apps.ecmwf.int/mars-catalogue` |
| **Copernicus CDS** | Global reanalysis obs | NetCDF | API | `cds.climate.copernicus.eu` |
| **GTS via Unidata IDD** | Real-time GTS | BUFR, TAC | LDM (academic) | `unidata.ucar.edu/data` |

### 4.2 Regional / National Sources

| Source | Country/Region | Format | Access |
|--------|----------------|--------|--------|
| **BMKG Open Data** | Indonesia | CSV, API | `dataonline.bmkg.go.id` |
| **JMA** | Japan | CSV, BUFR | `jma.go.jp/jma/en/Activities/` |
| **BoM** | Australia | JSON, CSV | `bom.gov.au/climate/data` |
| **Met Office MIDAS** | UK | ASCII | CEDA archive |
| **Météo-France** | France | CSV | `donneespubliques.meteofrance.fr` |
| **DWD CDC** | Germany | CSV, BUFR | `opendata.dwd.de/climate_environment/CDC/` |

### 4.3 Choosing the Right Source

```
Need global coverage?
├── YES → Need real-time?
│         ├── YES → WMO WIS 2.0 or GTS via LDM
│         └── NO  → ISD (most complete archive)
└── NO  → Need US data?
          ├── YES → IEM (easiest) or MADIS (most vars)
          └── NO  → Check national met service first
                    Then fall back to OGIMET (SYNOP/METAR scraping)
```

---

## 5. Python Environment Setup

### Conda environment

```yaml
# surface-obs-env.yml
name: surface-obs
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - numpy
  - pandas
  - xarray
  - scipy
  - matplotlib
  - cartopy
  - requests
  - eccodes           # BUFR/GRIB decoding (C library)
  - python-eccodes    # Python bindings for eccodes
  - cfgrib            # GRIB via xarray
  - metpy             # meteorological calculations
  - siphon            # Unidata data access (IDD, TDS, MADIS)
  - beautifulsoup4    # web scraping (OGIMET)
  - pip:
    - metar            # METAR string parsing
```

```bash
conda env create -f surface-obs-env.yml
conda activate surface-obs
```

### Pip alternative

```bash
pip install numpy pandas xarray scipy matplotlib cartopy \
    requests eccodes cfgrib metpy siphon beautifulsoup4 metar
```

### Quick import check

```python
import numpy as np
import pandas as pd
import xarray as xr
import eccodes
import metpy.calc as mpcalc
from metpy.units import units
from siphon.simplewebservice.iastate import IAStateUpperAir
print("All surface obs packages loaded successfully.")
```

---

## 6. Reading & Decoding Surface Data

### 6.1 ISD Data (NOAA NCEI)

ISD is the most accessible global surface archive. Each file covers one station-year.

**Download a station file:**

```python
import requests
import pandas as pd
import numpy as np

def download_isd(usaf, wban, year, out_dir="."):
    """Download one ISD station-year file from NOAA."""
    url = (f"https://www.ncei.noaa.gov/data/"
           f"global-hourly/access/{year}/{usaf}{wban}.csv")
    resp = requests.get(url)
    resp.raise_for_status()
    out_path = f"{out_dir}/{usaf}-{wban}-{year}.csv"
    with open(out_path, "w") as f:
        f.write(resp.text)
    return out_path
```

**Parse ISD CSV (global-hourly format):**

```python
def read_isd_csv(filepath):
    """Read NOAA ISD global-hourly CSV into a clean DataFrame."""
    df = pd.read_csv(filepath, low_memory=False, parse_dates=["DATE"])

    # Extract key fields from the compound columns
    result = pd.DataFrame()
    result["datetime"] = df["DATE"]
    result["station"] = df["STATION"]
    result["latitude"] = df["LATITUDE"]
    result["longitude"] = df["LONGITUDE"]
    result["elevation"] = df["ELEVATION"]

    # Temperature (stored as °C × 10 in TMP field: "value,quality")
    result["temp_c"] = df["TMP"].apply(parse_isd_field)
    # Dewpoint
    result["dewpoint_c"] = df["DEW"].apply(parse_isd_field)
    # Sea level pressure
    result["slp_hpa"] = df["SLP"].apply(parse_isd_field)
    # Wind speed
    result["wind_speed_ms"] = df["WND"].apply(
        lambda x: parse_isd_wind(x, "speed"))
    # Wind direction
    result["wind_dir_deg"] = df["WND"].apply(
        lambda x: parse_isd_wind(x, "direction"))
    # Visibility
    result["visibility_m"] = df["VIS"].apply(parse_isd_field)

    return result

def parse_isd_field(val):
    """Parse ISD compound field 'value,quality' → float or NaN."""
    if pd.isna(val):
        return np.nan
    parts = str(val).split(",")
    try:
        v = float(parts[0])
        if v == 9999 or v == 99999 or v == 999.9:
            return np.nan
        # Temperature and dewpoint are scaled ×10 in some columns
        return v
    except (ValueError, IndexError):
        return np.nan

def parse_isd_wind(val, component):
    """Parse ISD WND field: 'dir,dir_q,obs_type,speed,speed_q'."""
    if pd.isna(val):
        return np.nan
    parts = str(val).split(",")
    try:
        if component == "direction":
            v = float(parts[0])
            return np.nan if v == 999 else v
        else:  # speed
            v = float(parts[3])
            return np.nan if v == 9999 else v / 10.0  # stored as ×10
    except (ValueError, IndexError):
        return np.nan
```

### 6.2 METAR Strings

Use the `metar` Python library to parse METAR strings.

```python
from metar import Metar

def parse_metar(metar_string):
    """Parse a METAR string into a dict of values."""
    obs = Metar.Metar(metar_string)
    return {
        "station":    obs.station_id,
        "time":       obs.time,
        "temp_c":     obs.temp.value("C") if obs.temp else None,
        "dewpoint_c": obs.dewpoint.value("C") if obs.dewpoint else None,
        "wind_dir":   obs.wind_dir.value() if obs.wind_dir else None,
        "wind_speed_kt": obs.wind_speed.value("KT") if obs.wind_speed else None,
        "wind_gust_kt":  obs.wind_gust.value("KT") if obs.wind_gust else None,
        "visibility_m":  obs.vis.value("M") if obs.vis else None,
        "pressure_hpa":  obs.press.value("HPA") if obs.press else None,
        "weather":    obs.present_weather() if obs.weather else None,
        "sky":        str(obs.sky) if obs.sky else None,
    }

# Example
metar_str = "METAR WIII 150700Z 32005KT 9999 FEW020CB SCT025 31/24 Q1010 NOSIG"
parsed = parse_metar(metar_str)
print(f"Temp: {parsed['temp_c']:.1f}°C, Wind: {parsed['wind_dir']}° @ {parsed['wind_speed_kt']} kt")
```

**Bulk METAR from Iowa Environmental Mesonet (IEM):**

```python
def fetch_metar_iem(station, start, end):
    """Fetch METAR data from IEM for one station.

    Parameters
    ----------
    station : str   — ICAO code (e.g., "WIII")
    start   : str   — start date "YYYY-MM-DD"
    end     : str   — end date "YYYY-MM-DD"

    Returns
    -------
    pd.DataFrame with parsed METAR fields
    """
    url = (
        f"https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"
        f"station={station}&data=all"
        f"&year1={start[:4]}&month1={start[5:7]}&day1={start[8:10]}"
        f"&year2={end[:4]}&month2={end[5:7]}&day2={end[8:10]}"
        f"&tz=Etc/UTC&format=onlycomma&latlon=yes&elev=yes"
        f"&missing=M&trace=T&direct=no&report_type=3"
    )
    df = pd.read_csv(url, na_values=["M", "T"])
    df["valid"] = pd.to_datetime(df["valid"])
    return df
```

### 6.3 BUFR via eccodes

For BUFR files (e.g., from ECMWF, WMO WIS, or national services):

```python
import eccodes

def read_synop_bufr(filepath):
    """Read SYNOP reports from a BUFR file using eccodes."""
    records = []

    with open(filepath, "rb") as f:
        while True:
            bufr = eccodes.codes_bufr_new_from_file(f)
            if bufr is None:
                break

            try:
                # Unpack the data section
                eccodes.codes_set(bufr, "unpack", 1)

                rec = {}
                # Station identification
                rec["block"]   = eccodes.codes_get(bufr, "blockNumber")
                rec["station"] = eccodes.codes_get(bufr, "stationNumber")
                rec["wmo_id"]  = f"{rec['block']:02d}{rec['station']:03d}"
                rec["lat"]     = eccodes.codes_get(bufr, "latitude")
                rec["lon"]     = eccodes.codes_get(bufr, "longitude")

                # Time
                rec["year"]   = eccodes.codes_get(bufr, "year")
                rec["month"]  = eccodes.codes_get(bufr, "month")
                rec["day"]    = eccodes.codes_get(bufr, "day")
                rec["hour"]   = eccodes.codes_get(bufr, "hour")
                rec["minute"] = eccodes.codes_get(bufr, "minute")

                # Meteorological fields
                for key in ["airTemperature", "dewpointTemperature",
                            "windDirection", "windSpeed",
                            "pressureReducedToMeanSeaLevel",
                            "stationPressure", "horizontalVisibility",
                            "totalCloudCover"]:
                    try:
                        val = eccodes.codes_get(bufr, key)
                        # eccodes uses 2147483647 or similar for missing
                        rec[key] = val if val < 1e9 else None
                    except eccodes.CodesInternalError:
                        rec[key] = None

                records.append(rec)
            finally:
                eccodes.codes_release(bufr)

    df = pd.DataFrame(records)
    # Convert temperature from K to °C
    if "airTemperature" in df.columns:
        df["temp_c"] = df["airTemperature"] - 273.15
    if "dewpointTemperature" in df.columns:
        df["dewpoint_c"] = df["dewpointTemperature"] - 273.15
    # Pressure from Pa to hPa
    if "pressureReducedToMeanSeaLevel" in df.columns:
        df["slp_hpa"] = df["pressureReducedToMeanSeaLevel"] / 100.0

    return df
```

### 6.4 Siphon (Unidata Data Access)

Siphon provides Python access to Unidata's Thredds Data Server and related services.

```python
from siphon.simplewebservice.iastate import IAStateUpperAir
from siphon.catalog import TDSCatalog
from datetime import datetime

# Access latest surface obs from a Thredds catalog
cat = TDSCatalog(
    "https://thredds.ucar.edu/thredds/catalog/nws/metar/ncdecoded/catalog.xml"
)
latest = list(cat.datasets.values())[-1]
ncss = latest.subset()

# Query a bounding box
query = ncss.query()
query.lonlat_box(west=95, east=141, south=-11, north=6)  # Indonesia
query.time(datetime.utcnow())
query.variables("air_temperature", "dew_point_temperature",
                "wind_speed", "wind_from_direction", "air_pressure_at_sea_level")
data = ncss.get_data(query)
```

### 6.5 CSV from National Services

For agency CSV files (common pattern):

```python
def read_station_csv(filepath, datetime_col="datetime",
                     datetime_format="%Y-%m-%d %H:%M"):
    """Read a generic station CSV with consistent column naming."""
    df = pd.read_csv(filepath)
    df["datetime"] = pd.to_datetime(df[datetime_col], format=datetime_format)
    df = df.set_index("datetime").sort_index()

    # Standardize column names (adapt to your agency's naming)
    rename_map = {
        "t": "temp_c", "temperature": "temp_c", "temp": "temp_c",
        "rh": "rh_pct", "humidity": "rh_pct",
        "ws": "wind_speed_ms", "ff": "wind_speed_ms",
        "wd": "wind_dir_deg", "dd": "wind_dir_deg",
        "slp": "slp_hpa", "mslp": "slp_hpa",
        "rain": "precip_mm", "precipitation": "precip_mm",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    return df
```

---

## 7. Quality Control

Quality control (QC) is essential before using observations for verification. Bad data will produce meaningless verification scores.

### 7.1 QC Flag Framework

Define a consistent flag scheme:

```python
# QC flag values (following NCEP/MADIS convention)
QC_GOOD     = 0  # Passed all checks
QC_SUSPECT  = 1  # Suspicious but not rejected
QC_BAD      = 2  # Failed checks, should not be used
QC_MISSING  = 9  # Missing data

def init_qc_flags(df, variables):
    """Initialize QC flag columns for each variable."""
    for var in variables:
        df[f"{var}_qc"] = QC_GOOD
        df.loc[df[var].isna(), f"{var}_qc"] = QC_MISSING
    return df
```

### 7.2 Range Checks

Physical plausibility limits:

```python
# Climatological limits (adjust for your region)
RANGE_LIMITS = {
    "temp_c":       (-90.0,  60.0),
    "dewpoint_c":   (-90.0,  40.0),
    "slp_hpa":      (870.0,  1084.0),
    "wind_speed_ms": (0.0,   100.0),
    "wind_dir_deg":  (0.0,   360.0),
    "rh_pct":        (0.0,   100.0),
    "precip_mm":     (0.0,   500.0),  # per reporting period
    "visibility_m":  (0.0,   100000.0),
}

def range_check(df, limits=RANGE_LIMITS):
    """Flag values outside physical limits."""
    for var, (vmin, vmax) in limits.items():
        if var in df.columns:
            mask = (df[var] < vmin) | (df[var] > vmax)
            qc_col = f"{var}_qc"
            if qc_col in df.columns:
                df.loc[mask, qc_col] = QC_BAD
            else:
                df[qc_col] = QC_GOOD
                df.loc[mask, qc_col] = QC_BAD
    return df
```

### 7.3 Internal Consistency Checks

```python
def consistency_check(df):
    """Check physical relationships between variables."""
    # Dewpoint must not exceed temperature
    if "temp_c" in df.columns and "dewpoint_c" in df.columns:
        mask = df["dewpoint_c"] > df["temp_c"] + 0.5  # 0.5°C tolerance
        df.loc[mask, "dewpoint_c_qc"] = QC_BAD
        df.loc[mask, "temp_c_qc"] = np.where(
            df.loc[mask, "temp_c_qc"] == QC_GOOD, QC_SUSPECT,
            df.loc[mask, "temp_c_qc"]
        )

    # Wind direction must be missing or 0 when wind speed is 0
    if "wind_speed_ms" in df.columns and "wind_dir_deg" in df.columns:
        calm = df["wind_speed_ms"] == 0
        has_dir = df["wind_dir_deg"] > 0
        mask = calm & has_dir
        df.loc[mask, "wind_dir_deg_qc"] = QC_SUSPECT

    # RH and dewpoint should be consistent (if both present)
    if "rh_pct" in df.columns and "temp_c" in df.columns:
        mask = (df["rh_pct"] > 100.5)
        df.loc[mask, "rh_pct_qc"] = QC_BAD

    return df
```

### 7.4 Temporal Consistency (Spike Check)

```python
def spike_check(df, var, max_change, dt_hours=1):
    """Flag sudden jumps that exceed max_change per dt_hours.

    Parameters
    ----------
    df         : DataFrame with datetime index
    var        : column name to check
    max_change : maximum allowed change per time step
    dt_hours   : expected time step in hours
    """
    if var not in df.columns:
        return df

    diff = df[var].diff().abs()
    time_diff = df.index.to_series().diff().dt.total_seconds() / 3600

    # Normalize change rate to per-hour
    rate = diff / time_diff.replace(0, np.nan)
    max_rate = max_change / dt_hours

    mask = rate > max_rate
    qc_col = f"{var}_qc"
    if qc_col in df.columns:
        df.loc[mask, qc_col] = np.maximum(df.loc[mask, qc_col], QC_SUSPECT)

    return df

# Typical spike thresholds (per hour)
SPIKE_THRESHOLDS = {
    "temp_c":        10.0,   # 10°C/hr
    "dewpoint_c":    10.0,
    "slp_hpa":       6.0,    # 6 hPa/hr (rapid deepening threshold)
    "wind_speed_ms": 20.0,
}
```

### 7.5 Spatial Buddy Check

Compare a station's value against its neighbors:

```python
from scipy.spatial import cKDTree

def buddy_check(df_stations, var, max_deviation=3.0, min_buddies=3,
                search_radius_km=100):
    """Spatial buddy check — flag values deviating from neighbors.

    Parameters
    ----------
    df_stations    : DataFrame with columns [lat, lon, var, ...]
                     One row per station (single time step).
    var            : variable to check
    max_deviation  : max number of std deviations from buddy mean
    min_buddies    : minimum number of neighbors required
    search_radius_km : search radius in km
    """
    coords = np.deg2rad(df_stations[["lat", "lon"]].values)
    tree = cKDTree(coords)

    R_earth = 6371.0  # km
    search_rad = search_radius_km / R_earth  # radians

    qc_flags = np.full(len(df_stations), QC_GOOD)

    for i in range(len(df_stations)):
        neighbors = tree.query_ball_point(coords[i], search_rad)
        neighbors = [j for j in neighbors if j != i]

        if len(neighbors) < min_buddies:
            continue

        buddy_vals = df_stations[var].iloc[neighbors].dropna()
        if len(buddy_vals) < min_buddies:
            continue

        station_val = df_stations[var].iloc[i]
        if np.isnan(station_val):
            continue

        buddy_mean = buddy_vals.mean()
        buddy_std = buddy_vals.std()

        if buddy_std > 0:
            z_score = abs(station_val - buddy_mean) / buddy_std
            if z_score > max_deviation:
                qc_flags[i] = QC_BAD

    df_stations[f"{var}_qc"] = qc_flags
    return df_stations
```

### 7.6 Persistence Check

Flag stations that report the same value for too long:

```python
def persistence_check(df, var, max_repeats=6):
    """Flag values that repeat unchanged for max_repeats consecutive times."""
    if var not in df.columns:
        return df

    # Count consecutive identical values
    changes = df[var].ne(df[var].shift())
    groups = changes.cumsum()
    counts = groups.map(groups.value_counts())

    mask = counts >= max_repeats
    qc_col = f"{var}_qc"
    if qc_col in df.columns:
        df.loc[mask, qc_col] = np.maximum(df.loc[mask, qc_col], QC_SUSPECT)

    return df
```

### 7.7 Complete QC Pipeline

```python
def run_full_qc(df, variables=None):
    """Run all QC checks in sequence."""
    if variables is None:
        variables = ["temp_c", "dewpoint_c", "slp_hpa",
                     "wind_speed_ms", "wind_dir_deg"]

    # 1. Initialize flags
    df = init_qc_flags(df, variables)

    # 2. Range check
    df = range_check(df)

    # 3. Internal consistency
    df = consistency_check(df)

    # 4. Temporal spike check
    for var, threshold in SPIKE_THRESHOLDS.items():
        if var in variables:
            df = spike_check(df, var, threshold)

    # 5. Persistence check
    for var in variables:
        df = persistence_check(df, var)

    # Summary
    for var in variables:
        qc_col = f"{var}_qc"
        if qc_col in df.columns:
            total = (df[qc_col] != QC_MISSING).sum()
            bad = (df[qc_col] == QC_BAD).sum()
            suspect = (df[qc_col] == QC_SUSPECT).sum()
            if total > 0:
                print(f"{var}: {total} obs, "
                      f"{bad} bad ({100*bad/total:.1f}%), "
                      f"{suspect} suspect ({100*suspect/total:.1f}%)")

    return df
```

---

## 8. Data Processing

### 8.1 Unit Conversions

```python
def convert_units(df):
    """Apply standard unit conversions."""
    conversions = {
        # Wind: knots → m/s
        "wind_speed_kt": ("wind_speed_ms", lambda x: x * 0.514444),
        # Temperature: K → °C
        "temp_k": ("temp_c", lambda x: x - 273.15),
        # Temperature: °F → °C
        "temp_f": ("temp_c", lambda x: (x - 32) * 5/9),
        # Pressure: Pa → hPa
        "pressure_pa": ("slp_hpa", lambda x: x / 100.0),
        # Pressure: inHg → hPa
        "pressure_inhg": ("slp_hpa", lambda x: x * 33.8639),
        # Visibility: miles → meters
        "visibility_mi": ("visibility_m", lambda x: x * 1609.34),
        # Precipitation: inches → mm
        "precip_in": ("precip_mm", lambda x: x * 25.4),
    }

    for src_col, (dst_col, func) in conversions.items():
        if src_col in df.columns and dst_col not in df.columns:
            df[dst_col] = func(df[src_col])

    return df
```

### 8.2 Derived Variables

```python
def compute_derived(df):
    """Compute commonly needed derived variables."""
    # Relative humidity from T and Td (Magnus formula)
    if "temp_c" in df.columns and "dewpoint_c" in df.columns:
        if "rh_pct" not in df.columns:
            a, b = 17.625, 243.04
            gamma_t  = (a * df["temp_c"])     / (b + df["temp_c"])
            gamma_td = (a * df["dewpoint_c"]) / (b + df["dewpoint_c"])
            df["rh_pct"] = 100.0 * np.exp(gamma_td - gamma_t)

    # Wind U/V components from speed and direction
    if "wind_speed_ms" in df.columns and "wind_dir_deg" in df.columns:
        wd_rad = np.deg2rad(df["wind_dir_deg"])
        df["wind_u"] = -df["wind_speed_ms"] * np.sin(wd_rad)
        df["wind_v"] = -df["wind_speed_ms"] * np.cos(wd_rad)

    # Wet-bulb temperature (Stull 2011 approximation)
    if "temp_c" in df.columns and "rh_pct" in df.columns:
        T = df["temp_c"]
        RH = df["rh_pct"]
        df["wetbulb_c"] = (
            T * np.arctan(0.151977 * np.sqrt(RH + 8.313659))
            + np.arctan(T + RH)
            - np.arctan(RH - 1.676331)
            + 0.00391838 * RH**1.5 * np.arctan(0.023101 * RH)
            - 4.686035
        )

    # Heat index (for tropical stations)
    if "temp_c" in df.columns and "rh_pct" in df.columns:
        T_f = df["temp_c"] * 9/5 + 32  # Fahrenheit
        RH = df["rh_pct"]
        mask = T_f >= 80
        HI = pd.Series(np.nan, index=df.index)
        HI[mask] = (
            -42.379 + 2.04901523 * T_f[mask]
            + 10.14333127 * RH[mask]
            - 0.22475541 * T_f[mask] * RH[mask]
            - 0.00683783 * T_f[mask]**2
            - 0.05481717 * RH[mask]**2
            + 0.00122874 * T_f[mask]**2 * RH[mask]
            + 0.00085282 * T_f[mask] * RH[mask]**2
            - 0.00000199 * T_f[mask]**2 * RH[mask]**2
        )
        df["heat_index_c"] = (HI - 32) * 5/9

    return df
```

### 8.3 Temporal Aggregation

```python
def temporal_aggregate(df, freq="1h", methods=None):
    """Aggregate sub-hourly data to a target frequency.

    Parameters
    ----------
    df      : DataFrame with datetime index
    freq    : target frequency ('1h', '3h', '6h', '1D')
    methods : dict of {column: aggregation_method}
    """
    if methods is None:
        methods = {
            "temp_c":        "mean",
            "dewpoint_c":    "mean",
            "rh_pct":        "mean",
            "slp_hpa":       "mean",
            "wind_speed_ms": "mean",
            "wind_dir_deg":  "first",      # direction: use vector mean instead
            "wind_u":        "mean",        # average U/V, then recompute dir
            "wind_v":        "mean",
            "precip_mm":     "sum",         # precipitation = accumulation
            "visibility_m":  "min",         # worst visibility in period
            "temp_c_max":    "max",
            "temp_c_min":    "min",
        }

    # Only aggregate columns that exist
    agg_dict = {k: v for k, v in methods.items() if k in df.columns}
    result = df.resample(freq).agg(agg_dict)

    # Recompute wind direction from mean U/V
    if "wind_u" in result.columns and "wind_v" in result.columns:
        result["wind_dir_deg"] = (
            np.rad2deg(np.arctan2(-result["wind_u"], -result["wind_v"])) % 360
        )
        result["wind_speed_ms"] = np.sqrt(
            result["wind_u"]**2 + result["wind_v"]**2
        )

    return result
```

### 8.4 Station Pressure Reduction to MSL

```python
def reduce_to_msl(station_pressure_hpa, temp_c, elevation_m):
    """Reduce station pressure to mean sea level using the hypsometric equation.

    Parameters
    ----------
    station_pressure_hpa : float or array — station pressure in hPa
    temp_c               : float or array — temperature in °C
    elevation_m          : float — station elevation in meters

    Returns
    -------
    MSL pressure in hPa
    """
    # Virtual temperature (simplified, assuming dry)
    Tv = temp_c + 273.15  # K
    # Hypsometric equation
    g = 9.80665    # m/s²
    Rd = 287.05    # J/(kg·K)
    mslp = station_pressure_hpa * np.exp(g * elevation_m / (Rd * Tv))
    return mslp
```

---

## 9. Station Metadata & Network Selection

### 9.1 ISD Station History

The ISD station history file maps station IDs to locations and periods of record.

```python
def load_isd_history():
    """Load the ISD station history file."""
    url = ("https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv")
    df = pd.read_csv(url)
    df.columns = df.columns.str.strip()

    # Parse begin/end dates
    df["BEGIN"] = pd.to_datetime(df["BEGIN"], format="%Y%m%d", errors="coerce")
    df["END"]   = pd.to_datetime(df["END"],   format="%Y%m%d", errors="coerce")

    # Combine USAF + WBAN for unique ID
    df["station_id"] = df["USAF"].astype(str) + "-" + df["WBAN"].astype(str)

    return df

def find_stations_in_box(history, lat_min, lat_max, lon_min, lon_max,
                         min_year=2000):
    """Find ISD stations within a lat/lon bounding box."""
    mask = (
        (history["LAT"] >= lat_min) & (history["LAT"] <= lat_max) &
        (history["LON"] >= lon_min) & (history["LON"] <= lon_max) &
        (history["END"].dt.year >= min_year)
    )
    result = history[mask].copy()
    result = result.sort_values("END", ascending=False)
    return result[["station_id", "STATION NAME", "CTRY", "LAT", "LON",
                   "ELEV(M)", "BEGIN", "END"]]
```

### 9.2 WMO to ICAO Mapping

```python
def wmo_to_icao_mapping():
    """Common WMO block → ICAO prefix mapping for Southeast Asia."""
    return {
        # WMO Block : ICAO prefix (approximate)
        96: "WI",   # Indonesia (western)
        97: "WA",   # Indonesia (eastern)
        48: "VT",   # Thailand
        48: "WB",   # Malaysia / Brunei
        98: "RP",   # Philippines
        48: "WS",   # Singapore
        48: "VV",   # Vietnam
    }
```

### 9.3 Station Selection for Verification

```python
def select_verification_stations(stations_df, model_domain,
                                 min_record_years=5,
                                 max_elevation_diff_m=100):
    """Select stations suitable for NWP verification.

    Criteria:
    1. Within the model domain
    2. Sufficient record length
    3. Elevation compatible with model topography
    4. Not at extreme locations (mountain peaks, deep valleys)
    """
    selected = stations_df.copy()

    # Filter by domain
    selected = selected[
        (selected["lat"] >= model_domain["lat_min"]) &
        (selected["lat"] <= model_domain["lat_max"]) &
        (selected["lon"] >= model_domain["lon_min"]) &
        (selected["lon"] <= model_domain["lon_max"])
    ]

    # Filter by record length
    if "begin" in selected.columns and "end" in selected.columns:
        record_years = (selected["end"] - selected["begin"]).dt.days / 365.25
        selected = selected[record_years >= min_record_years]

    # Note: elevation check against model terrain requires the model
    # topography field — see Section 10 for implementation

    print(f"Selected {len(selected)} stations for verification")
    return selected
```

### 9.4 Station Map

```python
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def plot_station_map(stations_df, lat_col="lat", lon_col="lon",
                     name_col=None, title="Station Network",
                     extent=None):
    """Plot station locations on a map."""
    fig, ax = plt.subplots(
        figsize=(12, 8),
        subplot_kw={"projection": ccrs.PlateCarree()}
    )

    ax.add_feature(cfeature.LAND, facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle="--")

    ax.scatter(
        stations_df[lon_col], stations_df[lat_col],
        c="red", s=20, zorder=5, transform=ccrs.PlateCarree()
    )

    if name_col and len(stations_df) <= 50:
        for _, row in stations_df.iterrows():
            ax.text(row[lon_col] + 0.1, row[lat_col] + 0.1,
                    row[name_col], fontsize=6,
                    transform=ccrs.PlateCarree())

    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    else:
        ax.set_global()

    ax.gridlines(draw_labels=True, linewidth=0.3)
    ax.set_title(title)
    plt.tight_layout()
    return fig, ax
```

---

## 10. Preparing Obs for Verification

### 10.1 MET ASCII Point Observation Format

MET expects point observations in a specific 11-column ASCII format (or NetCDF). Here we generate the ASCII format:

```
Message_Type Station_ID Valid_Time Lat Lon Elev Var_Name Level Height QC_String Obs_Value
```

```python
# GRIB code mapping for MET
MET_VAR_CODES = {
    "temp_c":        ("TMP",  "2",   "Z2"),     # 2m temperature
    "dewpoint_c":    ("DPT",  "2",   "Z2"),     # 2m dewpoint
    "rh_pct":        ("RH",   "2",   "Z2"),     # 2m relative humidity
    "slp_hpa":       ("PRMSL","0",   "Z0"),     # sea level pressure
    "wind_speed_ms": ("WIND", "10",  "Z10"),    # 10m wind speed
    "wind_dir_deg":  ("WDIR", "10",  "Z10"),    # 10m wind direction
    "wind_u":        ("UGRD", "10",  "Z10"),    # 10m U-wind
    "wind_v":        ("VGRD", "10",  "Z10"),    # 10m V-wind
    "precip_mm":     ("APCP", "0",   "Z0"),     # accumulated precip
    "visibility_m":  ("VIS",  "0",   "Z0"),     # visibility
}

def obs_to_met_ascii(df, output_path, msg_type="ADPSFC"):
    """Convert a station DataFrame to MET ASCII point observation format.

    Parameters
    ----------
    df          : DataFrame with columns:
                  station_id, datetime, lat, lon, elevation,
                  and one or more variable columns (temp_c, slp_hpa, etc.)
    output_path : str — output file path
    msg_type    : str — MET message type (ADPSFC for surface)

    The DataFrame should have QC columns (var_qc) where 0=good, 2=bad.
    Only QC_GOOD (0) observations are written.
    """
    lines = []

    for _, row in df.iterrows():
        dt_str = row["datetime"].strftime("%Y%m%d_%H%M%S")
        sid = str(row["station_id"])
        lat = f"{row['lat']:.4f}"
        lon = f"{row['lon']:.4f}"
        elev = f"{row.get('elevation', -9999):.1f}"

        for var_col, (met_name, level, height) in MET_VAR_CODES.items():
            if var_col not in df.columns:
                continue

            val = row[var_col]
            if pd.isna(val):
                continue

            # Check QC flag if available
            qc_col = f"{var_col}_qc"
            if qc_col in df.columns and row[qc_col] == QC_BAD:
                continue

            # MET expects temperature in K, pressure in Pa
            if met_name == "TMP" or met_name == "DPT":
                val = val + 273.15  # °C → K
            elif met_name == "PRMSL":
                val = val * 100.0   # hPa → Pa

            qc_str = "NA"
            lines.append(
                f'{msg_type} {sid} {dt_str} {lat} {lon} {elev} '
                f'{met_name} {level} {height} {qc_str} {val:.4f}'
            )

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Wrote {len(lines)} observations to {output_path}")
```

### 10.2 Matching Observations to Model Grid

For point verification, you need to extract model values at station locations:

```python
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator

def match_obs_to_model(obs_df, model_ds, var_name,
                       method="nearest"):
    """Extract model values at observation station locations.

    Parameters
    ----------
    obs_df    : DataFrame with lat, lon columns
    model_ds  : xarray Dataset with lat/lon coordinates
    var_name  : model variable name to extract
    method    : "nearest" or "bilinear"

    Returns
    -------
    Array of model values at observation locations
    """
    obs_lats = obs_df["lat"].values
    obs_lons = obs_df["lon"].values

    if method == "nearest":
        # Build KD-tree from model grid
        if model_ds[var_name].ndim >= 2:
            mlat = model_ds["lat"].values
            mlon = model_ds["lon"].values

            if mlat.ndim == 1:
                mlon2d, mlat2d = np.meshgrid(mlon, mlat)
            else:
                mlat2d, mlon2d = mlat, mlon

            coords = np.column_stack([mlat2d.ravel(), mlon2d.ravel()])
            tree = cKDTree(coords)

            obs_coords = np.column_stack([obs_lats, obs_lons])
            _, indices = tree.query(obs_coords)

            data = model_ds[var_name].values
            if data.ndim > 2:
                data = data.squeeze()
            return data.ravel()[indices]

    elif method == "bilinear":
        mlat = model_ds["lat"].values
        mlon = model_ds["lon"].values
        data = model_ds[var_name].values.squeeze()

        interp = RegularGridInterpolator(
            (mlat, mlon), data,
            method="linear", bounds_error=False, fill_value=np.nan
        )
        return interp(np.column_stack([obs_lats, obs_lons]))

    else:
        raise ValueError(f"Unknown method: {method}")
```

### 10.3 Elevation Correction

Model and station elevations often differ. Apply a lapse-rate correction:

```python
def elevation_correction(model_temp, model_elev, station_elev,
                         lapse_rate=0.0065):
    """Correct model temperature for elevation difference.

    Parameters
    ----------
    model_temp    : model temperature at grid point (°C)
    model_elev    : model terrain height at grid point (m)
    station_elev  : actual station elevation (m)
    lapse_rate    : temperature lapse rate (°C/m), default 6.5°C/km

    Returns
    -------
    Corrected temperature (°C)
    """
    dz = station_elev - model_elev
    return model_temp - lapse_rate * dz
```

### 10.4 Time Matching

```python
def match_times(obs_df, model_times, tolerance_minutes=30):
    """Match observation times to nearest model output time.

    Parameters
    ----------
    obs_df         : DataFrame with datetime index or column
    model_times    : array of model output datetimes
    tolerance_minutes : max allowed time difference

    Returns
    -------
    DataFrame with added 'matched_model_time' column
    """
    model_times = pd.DatetimeIndex(model_times)
    obs_times = pd.DatetimeIndex(obs_df["datetime"])

    # Find nearest model time for each obs
    matched = []
    for ot in obs_times:
        diffs = abs(model_times - ot)
        idx = diffs.argmin()
        if diffs[idx] <= pd.Timedelta(minutes=tolerance_minutes):
            matched.append(model_times[idx])
        else:
            matched.append(pd.NaT)

    obs_df = obs_df.copy()
    obs_df["matched_model_time"] = matched

    n_matched = obs_df["matched_model_time"].notna().sum()
    print(f"Matched {n_matched}/{len(obs_df)} observations to model times")

    return obs_df
```

---

## 11. Indonesian & Tropical Context

### 11.1 BMKG Station Network

Indonesia's meteorological agency (BMKG) operates stations across the archipelago:

| Category | Count (approx.) | WMO Blocks | Notes |
|----------|-----------------|------------|-------|
| Synoptic stations | ~180 | 96, 97 | Report on GTS |
| AWS | ~300+ | N/A | Agency network, CSV data |
| Maritime stations | ~50 | 96, 97 | Coastal and island |
| Climate stations | ~150 | 96, 97 | Long records, less frequent reporting |

**WMO Block numbers for Indonesia:**
- Block 96: Western Indonesia (Sumatra, Java, Kalimantan, Sulawesi west)
- Block 97: Eastern Indonesia (Sulawesi east, Maluku, Papua, Nusa Tenggara)

### 11.2 Key Indonesian SYNOP Stations

| WMO ID | ICAO | Name | Lat | Lon | Elev (m) |
|--------|------|------|-----|-----|----------|
| 96749 | WIII | Soekarno-Hatta (Jakarta) | -6.12 | 106.66 | 8 |
| 96839 | WIIJ | Juanda (Surabaya) | -7.38 | 112.79 | 3 |
| 96745 | WIMM | Polonia/Kualanamu (Medan) | 3.56 | 98.67 | 25 |
| 96805 | WAAA | Hasanuddin (Makassar) | -5.06 | 119.55 | 14 |
| 97072 | WAJJ | Sentani (Jayapura) | -2.58 | 140.52 | 79 |
| 96933 | WRKK | El Tari (Kupang) | -10.17 | 123.67 | 105 |
| 96753 | WIBB | Sultan Syarif Kasim (Pekanbaru) | 0.46 | 101.45 | 31 |
| 96801 | WADD | Ngurah Rai (Bali) | -8.75 | 115.17 | 1 |

### 11.3 Tropical QC Considerations

Standard QC limits need adjustment for tropical regions:

```python
TROPICAL_RANGE_LIMITS = {
    "temp_c":        (10.0,   45.0),   # No sub-zero at sea level
    "dewpoint_c":    (10.0,   30.0),   # High dewpoints normal
    "slp_hpa":       (990.0,  1025.0), # Narrow range near equator
    "wind_speed_ms": (0.0,    50.0),   # Lower max (no mid-lat storms)
    "rh_pct":        (20.0,   100.0),  # Rarely below 20% in tropics
    "precip_mm":     (0.0,    300.0),  # High hourly rates possible (convective)
}

# Tropical-specific spike thresholds
TROPICAL_SPIKE_THRESHOLDS = {
    "temp_c":        8.0,    # Slightly lower (less diurnal variation)
    "dewpoint_c":    8.0,
    "slp_hpa":       4.0,    # Smaller pressure changes
    "wind_speed_ms": 15.0,
}
```

### 11.4 Accessing BMKG Data

```python
def fetch_bmkg_data(station_wmo, start_date, end_date):
    """Fetch BMKG observation data via the public API.

    Note: BMKG's data portal availability may change. Check
    dataonline.bmkg.go.id for current access methods.
    """
    # BMKG public data — format varies by product
    # The BMKG open data portal provides CSV downloads
    # For GTS data, use ISD as the most reliable source

    # Fallback: use ISD which includes Indonesian stations
    usaf = str(station_wmo)  # ISD USAF often matches WMO ID

    # Search ISD history for the station
    history = load_isd_history()
    matches = history[
        history["USAF"].astype(str).str.contains(usaf[:5])
    ]

    if len(matches) == 0:
        print(f"Station {station_wmo} not found in ISD. "
              "Try OGIMET for SYNOP text.")
        return None

    print(f"Found {len(matches)} ISD record(s) for station {station_wmo}")
    return matches
```

### 11.5 OGIMET SYNOP Retrieval

OGIMET provides web access to SYNOP reports. Use with care (rate limits apply):

```python
import requests
from bs4 import BeautifulSoup

def fetch_synop_ogimet(wmo_id, start, end):
    """Fetch decoded SYNOP data from OGIMET.

    Parameters
    ----------
    wmo_id : str — 5-digit WMO station ID
    start  : str — "YYYYMMDD0000"
    end    : str — "YYYYMMDD2300"

    Returns
    -------
    str — raw text response (decoded SYNOP)

    Note: Be respectful of OGIMET's resources. Add delays between
    requests and do not make bulk automated downloads.
    """
    url = (
        f"https://www.ogimet.com/cgi-bin/getsynop?"
        f"block={wmo_id}&begin={start}&end={end}"
    )

    import time
    time.sleep(2)  # Rate limiting

    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.text
```

### 11.6 Precipitation Verification in the Tropics

Tropical precipitation is dominated by convection, making it especially challenging to verify:

| Challenge | Impact | Mitigation |
|-----------|--------|------------|
| High spatial variability | Point obs ≠ gridbox average | Use multiple stations per gridbox; compute areal reduction factors |
| Diurnal cycle | Model timing errors penalize scores | Verify accumulation periods (6h, 24h) rather than instantaneous |
| Double penalty | Displacement errors count twice | Use FSS or object-based methods (MODE) |
| Gauge under-catch | Wind-induced losses, splash-out | Apply correction factors (e.g., Koschmieder) |
| Heavy tail distribution | Rare events dominate scores | Use ETS, GSS rather than simple accuracy |

---

## 12. End-to-End Example

### Workflow

```
[1. Download ISD data]
        |
[2. Parse & clean]
        |
[3. Run QC]
        |
[4. Load WRF output]
        |
[5. Match stations to grid]
        |
[6. Compute verification metrics]
        |
[7. Visualize results]
```

### Complete Script

```python
"""
End-to-end surface observation workflow:
  ISD download → QC → match to WRF → verification metrics
"""
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

# ── Step 1: Download & parse ISD data ──────────────────────────

# Jakarta Soekarno-Hatta (USAF=967490, WBAN=99999)
isd_file = download_isd("967490", "99999", 2024)
obs = read_isd_csv(isd_file)
print(f"Loaded {len(obs)} observations")

# ── Step 2: Filter to study period ─────────────────────────────

obs = obs.set_index("datetime").sort_index()
obs = obs["2024-01-15":"2024-01-20"]  # 5-day case study
print(f"Study period: {len(obs)} observations")

# ── Step 3: Quality control ────────────────────────────────────

# Use tropical-specific limits for Indonesian stations
RANGE_LIMITS.update(TROPICAL_RANGE_LIMITS)
SPIKE_THRESHOLDS.update(TROPICAL_SPIKE_THRESHOLDS)

obs = run_full_qc(obs, variables=["temp_c", "dewpoint_c", "slp_hpa",
                                   "wind_speed_ms", "wind_dir_deg"])

# Keep only good observations
qc_vars = ["temp_c", "dewpoint_c", "slp_hpa", "wind_speed_ms"]
for var in qc_vars:
    qc_col = f"{var}_qc"
    if qc_col in obs.columns:
        obs.loc[obs[qc_col] == QC_BAD, var] = np.nan

# ── Step 4: Load WRF output ───────────────────────────────────

wrf = xr.open_dataset("wrfout_d02_2024-01-15_00.nc")
t2 = wrf["T2"]          # 2-m temperature (K)
psfc = wrf["PSFC"]       # surface pressure (Pa)
u10 = wrf["U10"]         # 10-m U-wind (m/s)
v10 = wrf["V10"]         # 10-m V-wind (m/s)
wrf_lat = wrf["XLAT"].isel(Time=0).values
wrf_lon = wrf["XLONG"].isel(Time=0).values

# ── Step 5: Match stations to WRF grid ─────────────────────────

# Build KD-tree from WRF grid
wrf_coords = np.column_stack([wrf_lat.ravel(), wrf_lon.ravel()])
tree = cKDTree(wrf_coords)

# Station location
stn_lat, stn_lon = obs["latitude"].iloc[0], obs["longitude"].iloc[0]
_, grid_idx = tree.query([stn_lat, stn_lon])
j, i = np.unravel_index(grid_idx, wrf_lat.shape)

# Extract WRF time series at the station
wrf_t2_stn   = t2[:, j, i].values - 273.15   # K → °C
wrf_u10_stn  = u10[:, j, i].values
wrf_v10_stn  = v10[:, j, i].values
wrf_ws_stn   = np.sqrt(wrf_u10_stn**2 + wrf_v10_stn**2)
wrf_times    = pd.DatetimeIndex(wrf["Times"].values)

# Elevation correction
wrf_hgt = wrf["HGT"].isel(Time=0).values[j, i]
stn_elev = obs["elevation"].iloc[0]
wrf_t2_stn = elevation_correction(wrf_t2_stn, wrf_hgt, stn_elev)

# ── Step 6: Match times ───────────────────────────────────────

obs_reset = obs.reset_index()
obs_matched = match_times(obs_reset, wrf_times, tolerance_minutes=30)
obs_matched = obs_matched.dropna(subset=["matched_model_time"])

# Build paired arrays
pairs = []
for _, row in obs_matched.iterrows():
    t_idx = np.where(wrf_times == row["matched_model_time"])[0]
    if len(t_idx) == 0:
        continue
    t_idx = t_idx[0]
    pairs.append({
        "datetime":   row["datetime"],
        "obs_temp":   row["temp_c"],
        "fcst_temp":  wrf_t2_stn[t_idx],
        "obs_ws":     row["wind_speed_ms"],
        "fcst_ws":    wrf_ws_stn[t_idx],
    })

paired = pd.DataFrame(pairs).dropna()
print(f"Paired observations: {len(paired)}")

# ── Step 7: Verification metrics ──────────────────────────────

# Temperature metrics (using functions from Python Verification Guide)
# See: ../verification/python-verification-guide.md Section 4
obs_t = paired["obs_temp"].values
fcst_t = paired["fcst_temp"].values

me   = np.mean(fcst_t - obs_t)
mae  = np.mean(np.abs(fcst_t - obs_t))
rmse = np.sqrt(np.mean((fcst_t - obs_t)**2))
r    = np.corrcoef(fcst_t, obs_t)[0, 1]

print(f"\n=== Temperature Verification (2m) ===")
print(f"  ME:   {me:+.2f} °C")
print(f"  MAE:  {mae:.2f} °C")
print(f"  RMSE: {rmse:.2f} °C")
print(f"  r:    {r:.3f}")

# Wind speed metrics
obs_w = paired["obs_ws"].values
fcst_w = paired["fcst_ws"].values
mask = ~np.isnan(obs_w) & ~np.isnan(fcst_w)
obs_w, fcst_w = obs_w[mask], fcst_w[mask]

me_w   = np.mean(fcst_w - obs_w)
mae_w  = np.mean(np.abs(fcst_w - obs_w))
rmse_w = np.sqrt(np.mean((fcst_w - obs_w)**2))

print(f"\n=== Wind Speed Verification (10m) ===")
print(f"  ME:   {me_w:+.2f} m/s")
print(f"  MAE:  {mae_w:.2f} m/s")
print(f"  RMSE: {rmse_w:.2f} m/s")

# ── Step 8: Visualization ─────────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# (a) Time series
ax = axes[0]
ax.plot(paired["datetime"], paired["obs_temp"],
        "ko-", markersize=3, label="Observed")
ax.plot(paired["datetime"], paired["fcst_temp"],
        "r.-", markersize=3, label="WRF")
ax.set_ylabel("Temperature (°C)")
ax.set_title("2m Temperature")
ax.legend()
ax.tick_params(axis="x", rotation=45)

# (b) Scatter plot
ax = axes[1]
vmin = min(obs_t.min(), fcst_t.min()) - 1
vmax = max(obs_t.max(), fcst_t.max()) + 1
ax.scatter(obs_t, fcst_t, s=10, alpha=0.5)
ax.plot([vmin, vmax], [vmin, vmax], "k--", linewidth=0.5)
ax.set_xlabel("Observed (°C)")
ax.set_ylabel("Forecast (°C)")
ax.set_title(f"T2m (RMSE={rmse:.2f}°C, r={r:.3f})")
ax.set_xlim(vmin, vmax)
ax.set_ylim(vmin, vmax)
ax.set_aspect("equal")

# (c) Error distribution
ax = axes[2]
errors = fcst_t - obs_t
ax.hist(errors, bins=20, edgecolor="black", alpha=0.7)
ax.axvline(me, color="red", linestyle="--", label=f"ME={me:+.2f}")
ax.set_xlabel("Forecast Error (°C)")
ax.set_ylabel("Count")
ax.set_title("Error Distribution")
ax.legend()

plt.tight_layout()
plt.savefig("surface_verification.png", dpi=150)
plt.show()
```

---

## Appendix A: Variable Code Tables

### Common SYNOP/METAR Variables

| Variable | SYNOP Group | BUFR Descriptor | ISD Column | Unit |
|----------|-------------|-----------------|------------|------|
| Temperature | 1sTTT | 012101 | TMP | °C (×10) |
| Dewpoint | 2sTTT | 012103 | DEW | °C (×10) |
| Station pressure | 3PPPP | 010004 | — | hPa (×10) |
| MSL pressure | 4PPPP | 010051 | SLP | hPa (×10) |
| Wind direction | Nddff (dd) | 011001 | WND | degrees |
| Wind speed | Nddff (ff) | 011002 | WND | m/s (×10 in ISD) |
| Visibility | VV | 020001 | VIS | m |
| Total cloud | N | 020010 | — | oktas |
| Weather (present) | ww | 020003 | — | WMO code 4677 |
| Precipitation | 6RRRt | 013011 | AA1 | mm |
| Max temperature | 1sTTT (333) | 012111 | — | °C |
| Min temperature | 2sTTT (333) | 012112 | — | °C |

### WMO Cloud Codes

| Code | Meaning | METAR Equivalent |
|------|---------|-----------------|
| 0 | 0 oktas | SKC / CLR |
| 1-2 | 1-2 oktas | FEW |
| 3-4 | 3-4 oktas | SCT |
| 5-7 | 5-7 oktas | BKN |
| 8 | 8 oktas | OVC |
| 9 | Sky obscured | VV (vertical visibility) |

### Present Weather (WMO Code 4677) — Common Values

| Code | Weather |
|------|---------|
| 0-3 | No significant weather |
| 4 | Visibility reduced by smoke |
| 5 | Haze |
| 10 | Mist |
| 13 | Lightning |
| 21 | Recent rain |
| 51-55 | Drizzle (slight → heavy) |
| 61-65 | Rain (slight → heavy) |
| 80-82 | Rain showers (slight → violent) |
| 95 | Thunderstorm with rain |
| 97 | Thunderstorm, heavy rain |

---

## Appendix B: Resources & Links

### Documentation

| Resource | URL |
|----------|-----|
| WMO Manual on Codes (WMO-306) | `library.wmo.int/idurl/4/35713` |
| ISD Format Documentation | `ncei.noaa.gov/data/global-hourly/doc/` |
| BUFR Tables | `community.wmo.int/en/activity-areas/wis/latest-version` |
| eccodes Documentation | `confluence.ecmwf.int/display/ECC` |
| MetPy Documentation | `unidata.github.io/MetPy/` |
| Siphon Documentation | `unidata.github.io/siphon/` |
| metar (Python library) | `pypi.org/project/metar/` |

### Companion Guides in This Repo

| Guide | Focus |
|-------|-------|
| [Forecast Verification Notes](../verification/forecast-verification-notes.md) | Theory & formulas |
| [Python Verification Guide](../verification/python-verification-guide.md) | Python metrics implementation |
| [MET Verification Guide](../verification/met-verification-guide.md) | MET tools & workflows |
| [Ocean Observations Guide](ocean-obs-guide.md) | Buoys, Argo, tide gauges, satellite |
| [WRF Model Guide](../models/wrf-model-guide.md) | NWP model (generates forecasts to verify) |

### Data Access Quick Reference

| Need | Best Source | Python Tool |
|------|-----------|-------------|
| Global SYNOP/METAR archive | ISD (NCEI) | `requests` + `pandas` |
| US ASOS real-time | IEM | `requests` (IEM API) |
| METAR string parsing | Any METAR source | `metar` library |
| BUFR decoding | ECMWF / WMO | `eccodes` |
| Real-time GTS feed | Unidata IDD | `siphon` |
| Indonesian stations | ISD or OGIMET | `requests` + `beautifulsoup4` |
