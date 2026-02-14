# MET (Model Evaluation Tools): Practical Verification Guide

> A practical guide to MET v12.0 and METplus v6.0 for forecast verification workflows.
> Complements the theory-focused [Forecast Verification: Methods and Key Concepts](forecast-verification-notes.md).

---

## Table of Contents

1. [Introduction & Overview](#1-introduction--overview)
2. [Installation & Setup](#2-installation--setup)
3. [Input Data Formats](#3-input-data-formats)
4. [Configuration File System](#4-configuration-file-system)
5. [Core Verification Tools](#5-core-verification-tools)
   - [Point-Stat](#51-point-stat)
   - [Grid-Stat](#52-grid-stat)
   - [Ensemble-Stat](#53-ensemble-stat)
   - [Stat-Analysis](#54-stat-analysis)
6. [Spatial & Object-Based Verification](#6-spatial--object-based-verification)
   - [MODE](#61-mode)
   - [MTD](#62-mtd)
   - [Grid-Diag](#63-grid-diag)
7. [Tropical Cyclone Verification](#7-tropical-cyclone-verification)
   - [TC-Pairs](#71-tc-pairs)
   - [TC-Stat](#72-tc-stat)
   - [TC-Gen](#73-tc-gen)
   - [TC-RMW](#74-tc-rmw)
8. [Specialized Tools](#8-specialized-tools)
   - [PCP-Combine](#81-pcp-combine)
   - [Regrid-Data-Plane](#82-regrid-data-plane)
   - [Gen-Vx-Mask](#83-gen-vx-mask)
   - [Observation Converters](#84-observation-converters)
   - [Series-Analysis](#85-series-analysis)
   - [Wavelet-Stat](#86-wavelet-stat)
   - [IODA2NC](#87-ioda2nc)
9. [METplus Wrappers & Automation](#9-metplus-wrappers--automation)
10. [Visualization & Database Tools](#10-visualization--database-tools)
11. [Python Embedding](#11-python-embedding)
12. [Practical Workflows & Examples](#12-practical-workflows--examples)
13. [Tips, Troubleshooting & Best Practices](#13-tips-troubleshooting--best-practices)

---

## 1. Introduction & Overview

### What is MET?

**MET (Model Evaluation Tools)** is a comprehensive, community-supported verification toolkit developed by the **Developmental Testbed Center (DTC)**, a partnership between NCAR, NOAA, and the US Air Force. The core is written in **C++** for performance and released under the **Apache 2.0** license.

MET provides:
- Traditional verification statistics (continuous, categorical, probabilistic)
- Spatial / object-based verification methods (MODE, MTD)
- Tropical cyclone track and intensity verification
- Ensemble verification (rank histograms, CRPS, PIT)
- Aggregation, stratification, and confidence interval estimation

### MET vs METplus

| Component | Description |
|---|---|
| **MET** | Core C++ verification engine (standalone executables) |
| **METplus** | Python wrapper framework that orchestrates MET tools via configuration files |
| **METviewer** | MySQL-backed web application for interactive stat exploration |
| **METplotpy** | Python plotting library (performance diagrams, Taylor diagrams, etc.) |
| **METcalcpy** | Python calculation engine underlying METplotpy |
| **METdataio** | Data I/O utilities |
| **METexpress** | Lightweight web interface for quick stat viewing |

**Relationship:** METplus wraps MET. You can use MET tools directly (command-line), but METplus adds looping over times/fields, standardized configuration, and reproducible workflows. Most operational centers use both together.

### Version Alignment

| METplus Version | MET Version | Python Requirement |
|---|---|---|
| 6.0 | 12.0 | 3.10.4+ |
| 5.1 | 11.1 | 3.8.6+ |
| 5.0 | 11.0 | 3.8.6+ |
| 4.1 | 10.1 | 3.6.3+ |

### Comparison: MET vs WFRT/verif (Lightweight Alternatives)

| Feature | MET/METplus | WFRT `verif` |
|---|---|---|
| **Language** | C++ core + Python wrappers | Pure Python |
| **Scope** | Full operational verification suite | Quick point verification |
| **Spatial methods** | MODE, MTD, Wavelet-Stat | None |
| **TC verification** | TC-Pairs, TC-Stat, TC-Gen, TC-RMW | None |
| **Ensemble** | Full (CRPS, rank hist, PIT, spread-skill) | Limited |
| **Setup complexity** | Moderate–High (dependencies, compilation) | `pip install verif` |
| **Best for** | Operational centers, research, comprehensive studies | Quick exploratory checks, teaching |

**When to use which:** Use MET/METplus for formal verification studies, operational monitoring, and any workflow needing spatial or ensemble methods. Use `verif` for quick sanity checks during development.

### Relationship to the Verification Theory Notes

The companion document [forecast-verification-notes.md](forecast-verification-notes.md) covers the **theory**: what each metric means, why it matters, and how to interpret it. This guide covers the **practice**: how to compute those metrics using MET/METplus tools, configure the software, and build end-to-end verification workflows.

---

## 2. Installation & Setup

### 2.1 Source Compilation

MET compilation requires several external libraries. The DTC provides a helper script `compile_MET_all.sh` that automates the build.

#### Prerequisites

| Library | Minimum Version | Notes |
|---|---|---|
| **NetCDF-C** | 4.6+ | Required. With HDF5 support |
| **NetCDF-CXX** | 4.3+ | C++ bindings |
| **HDF5** | 1.10+ | Required by NetCDF |
| **zlib** | 1.2.5+ | Compression support |
| **GSL** | 2.1+ | GNU Scientific Library for statistics |
| **BUFRLIB** | 11.3+ | NCEP BUFR decoder |
| **G2CLIB** | 1.6+ | GRIB2 support |
| **GRIB2C** | — | GRIB2 C library |

**Optional libraries:**
- **Cairo + FreeType** — Graphics support (needed for MODE/MTD plots)
- **Python 3.10+** — Python embedding support
- **ECKIT / ATLAS** — ECMWF grid support

#### Compilation Steps

```bash
# Download the compile script and source tarball
wget https://dtcenter.org/sites/default/files/community-code/met/compile_scripts/compile_MET_all.sh
wget https://dtcenter.org/sites/default/files/community-code/met/met-12.0.0.tar.gz

# Set environment
export TEST_BASE=/opt/met
export COMPILER=gnu_9.3.0
export MET_SUBDIR=${TEST_BASE}
export MET_TARBALL=met-12.0.0.tar.gz
export USE_MODULES=FALSE

# Optional: enable Python embedding
export MET_PYTHON=/usr/bin/python3
export MET_PYTHON_BIN_EXE=/usr/bin/python3
export MET_PYTHON_CC="-I/usr/include/python3.10"
export MET_PYTHON_LD="-L/usr/lib -lpython3.10"

# Run the build (compiles all external libs + MET)
chmod +x compile_MET_all.sh
./compile_MET_all.sh 2>&1 | tee compile.log
```

#### Verify Installation

```bash
# Check MET version
${TEST_BASE}/bin/met_config_info

# Run a quick test
${TEST_BASE}/bin/point_stat --version
```

### 2.2 Docker

The simplest way to get started — no compilation required.

```bash
# Pull the official MET image
docker pull dtcenter/met:12.0.0

# Run interactively
docker run -it --rm \
    -v /path/to/data:/data \
    dtcenter/met:12.0.0 /bin/bash

# Run a specific tool
docker run --rm \
    -v /path/to/data:/data \
    dtcenter/met:12.0.0 \
    point_stat /data/forecast.grib2 /data/obs.nc /data/PointStatConfig

# METplus Docker (includes MET + METplus + METplotpy)
docker pull dtcenter/metplus:6.0.0
docker run -it --rm \
    -v /path/to/data:/data \
    dtcenter/metplus:6.0.0 /bin/bash
```

### 2.3 Singularity / Apptainer (HPC)

Most HPC clusters don't allow Docker. Use Singularity (now Apptainer) instead.

```bash
# Build from Docker Hub
singularity pull met_12.0.0.sif docker://dtcenter/met:12.0.0
singularity pull metplus_6.0.0.sif docker://dtcenter/metplus:6.0.0

# Run interactively
singularity shell --bind /scratch:/scratch met_12.0.0.sif

# Run a tool
singularity exec --bind /scratch:/scratch met_12.0.0.sif \
    point_stat /scratch/fcst.grb2 /scratch/obs.nc /scratch/PointStatConfig

# In a SLURM job script
#!/bin/bash
#SBATCH --job-name=met_verify
#SBATCH --ntasks=1
#SBATCH --time=02:00:00

module load singularity   # or apptainer
singularity exec --bind /scratch met_12.0.0.sif \
    grid_stat /scratch/fcst.grb2 /scratch/obs.nc /scratch/GridStatConfig
```

### 2.4 METplus Installation (Python)

```bash
# Clone METplus
git clone https://github.com/dtcenter/METplus.git
cd METplus
git checkout main_v6.0

# Install dependencies
pip install -r requirements.txt

# Set environment variables
export METPLUS_PATH=/path/to/METplus
export MET_INSTALL_DIR=/opt/met          # path to compiled MET
export PATH=${METPLUS_PATH}/ush:${PATH}

# Verify
run_metplus.py --version
```

#### METplus Configuration (minimum)

Create a user configuration file `user.conf`:

```ini
[config]
INPUT_BASE = /data/input
OUTPUT_BASE = /data/output
MET_INSTALL_DIR = /opt/met
LOG_LEVEL = DEBUG
```

### 2.5 Environment Verification Checklist

```bash
# Check MET tools are accessible
which point_stat grid_stat ensemble_stat mode tc_pairs

# Check METplus
run_metplus.py --version

# Check Python embedding (optional)
python3 -c "import met_point_obs; print('MET Python embedding OK')"

# Check data I/O support
ncdump --version        # NetCDF
wgrib2 --version        # GRIB2
```

---

## 3. Input Data Formats

### 3.1 Gridded Data (Forecasts)

MET natively reads:

| Format | Extension | Notes |
|---|---|---|
| **GRIB1** | `.grb`, `.grib` | Legacy format, still common for reanalysis |
| **GRIB2** | `.grb2`, `.grib2` | Standard for most operational NWP models (GFS, ECMWF, etc.) |
| **NetCDF** | `.nc` | Must be **CF-compliant** with recognizable coordinate variables |
| **Python embedding** | any | Custom Python script provides data to MET |

#### GRIB Field Specification

Fields in GRIB files are specified by name and level:

```
# MET config syntax for GRIB fields
field = [
   { name = "TMP";  level = "P850"; },     // Temperature at 850 hPa
   { name = "UGRD"; level = "P500"; },     // U-wind at 500 hPa
   { name = "APCP"; level = "A06"; },      // 6-h accumulated precip
   { name = "APCP"; level = "A24"; },      // 24-h accumulated precip
   { name = "TMP";  level = "Z2"; },       // 2-m temperature
   { name = "RH";   level = "P700"; }      // Relative humidity 700 hPa
];
```

#### NetCDF Requirements

MET expects NetCDF files with:
- CF-compliant coordinate variables (`lat`, `lon`, `time`)
- Standard variable names, or explicit field specification in config
- Projection information in global attributes or via config override

```
# For non-standard NetCDF, specify variable directly
field = [
   { name = "temperature"; level = "(*, *)"; }
];
```

### 3.2 Point Observations

MET expects observations in its internal **NetCDF point observation format** (11-column structure). Raw observations must be converted first:

| Source Format | Converter Tool | Description |
|---|---|---|
| **PrepBUFR** | `pb2nc` | NCEP operational obs (radiosondes, METAR, ships, buoys) |
| **ASCII text** | `ascii2nc` | Custom CSV/text observation files |
| **MADIS** | `madis2nc` | NOAA Meteorological Assimilation Data Ingest System |
| **IODA** | `ioda2nc` | JEDI observation format |
| **Python embedding** | `PYTHON_PANDAS` | Custom Python script provides point obs |

#### ASCII Point Observation Format

MET's `ascii2nc` expects a specific 11-column format:

```
# Column layout for ascii2nc input:
# Message_Type Station_ID Valid_Time(YYYYMMDD_HHMMSS) Lat Lon Elev
# GRIB_Code Level Height QC_String Obs_Value
#
# Example:
ADPSFC  KDEN  20240115_120000  39.86  -104.67  1625.0  11  0  2.0  NA  273.15
ADPSFC  KDEN  20240115_120000  39.86  -104.67  1625.0  32  0  10.0  NA  5.2
ADPSFC  KDEN  20240115_120000  39.86  -104.67  1625.0  33  0  10.0  NA  8.1
```

| Column | Description | Example |
|---|---|---|
| Message_Type | Observation type (ADPSFC, ADPUPA, SFCSHP, etc.) | `ADPSFC` |
| Station_ID | Station identifier | `KDEN` |
| Valid_Time | Observation time (YYYYMMDD_HHMMSS) | `20240115_120000` |
| Lat | Latitude (degrees north) | `39.86` |
| Lon | Longitude (degrees east) | `-104.67` |
| Elev | Elevation (meters) | `1625.0` |
| GRIB_Code | Variable code (11=TMP, 32=UGRD, 33=VGRD, 61=APCP...) | `11` |
| Level | Pressure level (hPa) or 0 for surface | `0` |
| Height | Height above ground (m) | `2.0` |
| QC_String | Quality control flag (NA if none) | `NA` |
| Obs_Value | Observation value | `273.15` |

#### Common GRIB Codes for Point Observations

| Code | Variable | Units |
|---|---|---|
| 11 | Temperature (TMP) | K |
| 32 | U-wind (UGRD) | m/s |
| 33 | V-wind (VGRD) | m/s |
| 52 | Relative Humidity (RH) | % |
| 61 | Total Precipitation (APCP) | kg/m² |
| 1 | Pressure (PRES) | Pa |
| 7 | Geopotential Height (HGT) | gpm |
| 17 | Dewpoint Temperature (DPT) | K |
| 20 | Visibility (VIS) | m |
| 31 | Wind Direction (WDIR) | deg |
| 32 | Wind Speed (WIND) | m/s |

### 3.3 Data Format Conversion Tips

```bash
# Convert PrepBUFR to MET NetCDF
pb2nc prepbufr.gdas.20240115.t12z.nr \
      pb2nc_output.nc \
      PB2NCConfig

# Convert ASCII text to MET NetCDF
ascii2nc my_observations.txt ascii2nc_output.nc

# Convert MADIS surface obs
madis2nc /path/to/madis/metar/20240115_1200.nc \
         madis2nc_output.nc \
         -type metar

# Convert IODA obs (JEDI format)
ioda2nc /path/to/ioda_obs.nc ioda2nc_output.nc IODA2NCConfig

# Check what's in a GRIB2 file
wgrib2 -s forecast.grb2 | head -20

# Check what's in a NetCDF file
ncdump -h forecast.nc
```

---

## 4. Configuration File System

MET uses two levels of configuration:

1. **METplus conf files** (`.conf`) — Python-style INI, orchestrates workflows
2. **MET config files** — C++-style syntax, controls individual tool behavior

### 4.1 METplus Configuration (.conf)

#### Basic Syntax

```ini
[config]
# Tool to run
PROCESS_LIST = PcpCombine, GridStat

# Time looping
LOOP_BY = VALID
VALID_TIME_FMT = %Y%m%d%H
VALID_BEG = 2024011500
VALID_END = 2024011600
VALID_INCREMENT = 12H

# Input templates
FCST_GRID_STAT_INPUT_DIR = /data/forecasts
FCST_GRID_STAT_INPUT_TEMPLATE = wrf_d01_{valid?fmt=%Y%m%d%H}.grb2
OBS_GRID_STAT_INPUT_DIR = /data/obs
OBS_GRID_STAT_INPUT_TEMPLATE = obs_{valid?fmt=%Y%m%d%H}.nc

# Output
GRID_STAT_OUTPUT_DIR = {OUTPUT_BASE}/grid_stat

[dir]
INPUT_BASE = /data
OUTPUT_BASE = /data/output

[filename_templates]
# Templates can also go here
```

#### Variable Substitution

METplus supports powerful variable substitution:

| Syntax | Description | Example |
|---|---|---|
| `{valid?fmt=%Y%m%d%H}` | Valid time with format | `2024011512` |
| `{init?fmt=%Y%m%d%H}` | Initialization time | `2024011500` |
| `{lead?fmt=%HHH}` | Lead time (hours, 3-digit) | `012` |
| `{lead?fmt=%3H}` | Lead time (hours, 3-digit) | `012` |
| `{INPUT_BASE}` | Reference to another variable | `/data` |
| `{ENV[HOME]}` | Environment variable | `/home/user` |
| `{valid?fmt=%Y%m%d%H?shift=-6H}` | Time with shift | `2024011506` |

#### Looping Strategies

```ini
# Loop by valid time (most common)
LOOP_BY = VALID
VALID_TIME_FMT = %Y%m%d%H
VALID_BEG = 2024011500
VALID_END = 2024013100
VALID_INCREMENT = 6H

# Loop by init time
LOOP_BY = INIT
INIT_TIME_FMT = %Y%m%d%H
INIT_BEG = 2024011500
INIT_END = 2024013100
INIT_INCREMENT = 12H

# Multiple lead times
LEAD_SEQ = 6, 12, 18, 24, 36, 48, 72

# Multiple lead time groups
LEAD_SEQ_1 = 0, 6, 12, 18, 24
LEAD_SEQ_1_LABEL = Day1
LEAD_SEQ_2 = 30, 36, 42, 48
LEAD_SEQ_2_LABEL = Day2
```

### 4.2 MET Configuration Files

MET config files use a C++-like syntax with named blocks:

#### General Structure

```c
// MET config file (e.g., GridStatConfig)
// Comments use // or /* */

model = "WRF";         // Model name tag
obtype = "ANALYS";     // Observation type tag

////////////////////////////////////////////////////////////////////////////////
// Forecast and observation field specification

fcst = {
   field = [
      { name = "TMP"; level = "P850"; },
      { name = "TMP"; level = "P500"; },
      { name = "UGRD"; level = "P250"; },
      { name = "VGRD"; level = "P250"; }
   ];
}

obs = fcst;   // Obs uses same field definitions as forecast

////////////////////////////////////////////////////////////////////////////////
// Regridding (match forecast to obs grid or vice versa)

regrid = {
   to_grid    = FCST;           // FCST, OBS, or a grid specification
   method     = BILINEAR;       // BILINEAR, NEAREST, BUDGET, etc.
   width      = 2;
   vld_thresh = 0.5;
}

////////////////////////////////////////////////////////////////////////////////
// Masking — where to compute statistics

mask = {
   grid  = [ "FULL" ];          // Full domain
   poly  = [ "/masks/CONUS.poly", "/masks/east.poly" ];
   sid   = [];                  // Station ID list (Point-Stat only)
}
```

### 4.3 Common Configuration Blocks

#### Output Flag — Control Which Line Types Are Written

```c
// Grid-Stat output flags
output_flag = {
   fho    = BOTH;    // Forecast-Hit-Observation counts → .stat
   ctc    = BOTH;    // Contingency Table Counts
   cts    = BOTH;    // Contingency Table Statistics (CSI, GSS, etc.)
   mctc   = NONE;    // Multi-Category Contingency Table Counts
   mcts   = NONE;    // Multi-Category Contingency Table Statistics
   cnt    = BOTH;    // Continuous Statistics (ME, RMSE, correlation, etc.)
   sl1l2  = BOTH;    // Scalar L1L2 partial sums (for aggregation)
   sal1l2 = NONE;    // Scalar Anomaly L1L2 partial sums
   vl1l2  = BOTH;    // Vector L1L2 partial sums
   val1l2 = NONE;    // Vector Anomaly L1L2
   vcnt   = BOTH;    // Vector Continuous Statistics
   pct    = NONE;    // Probability Contingency Table
   pstd   = NONE;    // Probability Statistics (BSS, reliability)
   pjc    = NONE;    // Probability Joint/Conditional
   prc    = NONE;    // Probability ROC
   eclv   = NONE;    // Economic Cost/Loss Value
   nbrctc = NONE;    // Neighborhood Contingency Table Counts
   nbrcts = NONE;    // Neighborhood Contingency Table Statistics
   nbrcnt = NONE;    // Neighborhood Continuous Statistics
   grad   = NONE;    // Gradient Statistics (S1 score)
   dmap   = NONE;    // Distance Map statistics
}
// Values: NONE = don't compute, STAT = write to .stat only, BOTH = .stat + individual files
```

#### Interpolation Methods

```c
interp = {
   field      = BOTH;           // FCST, OBS, or BOTH
   vld_thresh = 1.0;
   shape      = SQUARE;         // SQUARE or CIRCLE
   type = [
      { method = NEAREST; width = 1; },    // Single grid point
      { method = BILINEAR; width = 2; },   // Bilinear interpolation
      { method = UW_MEAN; width = 3; }     // Unweighted mean 3×3
   ];
}
```

#### Categorical Thresholds

```c
// Categorical threshold syntax
// Operators: >, >=, <, <=, ==, !=, between (inclusive)
cat_thresh = [ >0.0, >=1.0, >=5.0, >=10.0, >=25.0, >=50.0 ];

// Threshold shorthand:
//   gt0     = >0
//   ge1     = >=1
//   lt5     = <5
//   le10    = <=10
//   eq0     = ==0
//   ne0     = !=0
//   bt1&5   = between 1 and 5 (inclusive)

// Probability thresholds (for probabilistic verification)
prob_cat_thresh = [ >=0.10, >=0.25, >=0.50, >=0.75, >=0.90 ];
```

#### Climatology

```c
climo_mean = {
   file_name = [ "/data/climo/climo_mean.nc" ];
   field     = [];              // Defaults to match forecast field
   regrid    = { method = BILINEAR; width = 2; vld_thresh = 0.5; }
   time_interp_method = DW_MEAN;   // Day-weighted mean
   day_interval = 31;
   hour_interval = 6;
}

climo_stdev = {
   file_name = [ "/data/climo/climo_stdev.nc" ];
   field     = [];
   regrid    = { method = BILINEAR; width = 2; vld_thresh = 0.5; }
   time_interp_method = DW_MEAN;
   day_interval = 31;
   hour_interval = 6;
}

// With climatology, MET computes anomaly correlation (ACC) and skill scores
// relative to climatology.
```

### 4.4 Field Naming and Level Specification

#### Level Types

| Prefix | Meaning | Example |
|---|---|---|
| `P` | Pressure level (hPa) | `P850`, `P500`, `P250` |
| `Z` | Height above ground (m) | `Z2` (2-m), `Z10` (10-m) |
| `L` | Generic level | `L0` (surface), `L1` |
| `A` | Accumulation period (hours) | `A06`, `A24` |
| `R` | Record number | `R1`, `R2` |
| `(*)` | All levels (NetCDF) | `(*)` |
| `(*,*)` | 2D field (NetCDF) | `(*,*)` |

#### Common Field Specifications

```c
// Pressure-level fields
field = [
   { name = "TMP";  level = [ "P1000", "P850", "P700", "P500", "P250" ]; },
   { name = "HGT";  level = [ "P500" ]; },
   { name = "UGRD"; level = [ "P850", "P250" ]; },
   { name = "VGRD"; level = [ "P850", "P250" ]; },
   { name = "RH";   level = [ "P850", "P700" ]; }
];

// Surface / near-surface fields
field = [
   { name = "TMP";   level = "Z2"; },       // 2-m temperature
   { name = "DPT";   level = "Z2"; },       // 2-m dewpoint
   { name = "UGRD";  level = "Z10"; },      // 10-m U wind
   { name = "VGRD";  level = "Z10"; },      // 10-m V wind
   { name = "PRMSL"; level = "Z0"; },       // Mean sea-level pressure
   { name = "APCP";  level = "A24"; }       // 24-h accumulated precip
];

// NetCDF fields (variable name directly)
field = [
   { name = "T2";    level = "(*,*)"; },                         // WRF 2-m temp
   { name = "RAINNC"; level = "(*,*)"; cat_thresh = [ >0.0 ]; }  // WRF non-convective rain
];

// Derived / combined wind speed
field = [
   { name = "WIND"; level = "Z10"; }       // MET computes from UGRD + VGRD
];
```

---

## 5. Core Verification Tools

### 5.1 Point-Stat

**Purpose:** Verify gridded forecasts against point observations (stations, radiosondes, buoys, etc.).

**Produces:** Continuous stats (ME, RMSE, MAE, correlation), categorical stats (CSI, POD, FAR), partial sums for later aggregation.

#### Command

```bash
point_stat \
   forecast.grb2 \
   obs_preprocessed.nc \
   PointStatConfig \
   -outdir /output/point_stat \
   -v 2
```

#### Key Config Sections

```c
// PointStatConfig

model = "WRF_d01";
obtype = "ADPSFC";

////////////////////////////////////////////////////////////////////////////////
// Fields to verify

fcst = {
   field = [
      { name = "TMP";  level = "Z2";  message_type = [ "ADPSFC" ]; },
      { name = "DPT";  level = "Z2";  message_type = [ "ADPSFC" ]; },
      { name = "UGRD"; level = "Z10"; message_type = [ "ADPSFC" ]; },
      { name = "VGRD"; level = "Z10"; message_type = [ "ADPSFC" ]; },
      { name = "TMP";  level = [ "P850", "P500", "P250" ];
                        message_type = [ "ADPUPA" ]; }
   ];
}

obs = fcst;

////////////////////////////////////////////////////////////////////////////////
// Spatial matching of obs to forecast grid

interp = {
   vld_thresh = 1.0;
   shape      = SQUARE;
   type = [
      { method = BILINEAR; width = 2; }
   ];
}

////////////////////////////////////////////////////////////////////////////////
// Time window for matching obs to forecast valid time

obs_window = {
   beg = -5400;    // -1.5 hours (seconds)
   end =  5400;    // +1.5 hours
}

////////////////////////////////////////////////////////////////////////////////
// Message types to process

message_type = [ "ADPSFC", "ADPUPA", "SFCSHP" ];

////////////////////////////////////////////////////////////////////////////////
// Verification masks

mask = {
   grid  = [ "FULL" ];
   poly  = [];
   sid   = [ "/masks/stations_java.txt" ];   // Specific station list
}

////////////////////////////////////////////////////////////////////////////////
// Duplicate handling and quality marks

duplicate_flag = UNIQUE;   // NONE, UNIQUE
obs_quality_inc = [];      // Include all quality marks
obs_quality_exc = [];

////////////////////////////////////////////////////////////////////////////////
// Output
output_flag = {
   fho    = BOTH;
   ctc    = BOTH;
   cts    = BOTH;
   cnt    = BOTH;
   sl1l2  = BOTH;
   vl1l2  = BOTH;
   vcnt   = BOTH;
   mpr    = NONE;    // Matched Pair Records — set BOTH for detailed pair output
}
```

#### Output Line Types

| Line Type | Content |
|---|---|
| **FHO** | Forecast, Hit, Observation rates |
| **CTC** | Contingency Table Counts (hits, misses, false alarms, correct negatives) |
| **CTS** | Contingency Table Statistics (BASER, CSI, POD, FAR, GSS, HK, HSS, etc.) |
| **CNT** | Continuous statistics (ME, MAE, RMSE, ESTDEV, FBAR, OBAR, PR_CORR, etc.) |
| **SL1L2** | Scalar partial sums (FBAR, OBAR, FFBAR, FOBAR, OOBAR) |
| **VL1L2** | Vector partial sums (U/V components) |
| **VCNT** | Vector continuous stats (speed bias, direction error) |
| **MPR** | Matched pairs (individual forecast-observation pairs) |

### 5.2 Grid-Stat

**Purpose:** Verify gridded forecasts against gridded observations (analyses, satellite, radar composites). Supports neighborhood (fuzzy) verification.

#### Command

```bash
grid_stat \
   forecast.grb2 \
   gridded_obs.nc \
   GridStatConfig \
   -outdir /output/grid_stat \
   -v 2
```

#### Key Config

```c
// GridStatConfig

model = "WRF_d01";
obtype = "STAGE4";      // Stage IV radar-gauge analysis

fcst = {
   field = [
      {
         name       = "APCP";
         level      = "A24";
         cat_thresh = [ >0.0, >=2.54, >=6.35, >=12.7, >=25.4, >=50.8 ];
         // mm thresholds: trace, 0.1", 0.25", 0.5", 1.0", 2.0"
      }
   ];
}

obs = {
   field = [
      {
         name       = "APCP_24";
         level      = "(*,*)";
         cat_thresh = [ >0.0, >=2.54, >=6.35, >=12.7, >=25.4, >=50.8 ];
      }
   ];
}

regrid = {
   to_grid    = FCST;
   method     = BUDGET;         // Budget interpolation for precip
   width      = 2;
   vld_thresh = 0.5;
}

////////////////////////////////////////////////////////////////////////////////
// Neighborhood verification

nbrhd = {
   field      = BOTH;
   shape      = CIRCLE;
   width      = [ 1, 3, 5, 7, 9 ];   // Grid squares (radius in grid points)
   cov_thresh = [ >=0.5 ];
   vld_thresh = 1.0;
}

////////////////////////////////////////////////////////////////////////////////
// Output — enable neighborhood line types

output_flag = {
   fho    = BOTH;
   ctc    = BOTH;
   cts    = BOTH;
   cnt    = BOTH;
   sl1l2  = BOTH;
   nbrctc = BOTH;      // Neighborhood Contingency Table Counts
   nbrcts = BOTH;      // Neighborhood Contingency Table Statistics (FSS)
   nbrcnt = BOTH;      // Neighborhood Continuous Statistics
   grad   = BOTH;      // Gradient statistics (S1 score)
   dmap   = BOTH;      // Distance map (Hausdorff, Baddeley, etc.)
}
```

#### Neighborhood (Fuzzy) Verification

Grid-Stat computes the **Fractions Skill Score (FSS)** using neighborhood windows. FSS answers: "At what spatial scale does the forecast have useful skill?"

```
FSS = 1 - FBS / FBS_worst

where:
  FBS = (1/N) Σ (P_fcst - P_obs)²     # Fractions Brier Score
  P_fcst = fraction of grid points exceeding threshold within neighborhood
  P_obs  = same for observations
```

Key output columns in NBRCTS:
- `NBRHD_WIDTH` — neighborhood size
- `FSS` — Fractions Skill Score (0=no skill, 1=perfect)

### 5.3 Ensemble-Stat

**Purpose:** Verify ensemble forecasts. Computes rank histograms (Talagrand), CRPS, ensemble spread, PIT histograms, and reliability.

#### Command

```bash
ensemble_stat \
   -ens forecast_mem01.grb2 forecast_mem02.grb2 ... forecast_mem20.grb2 \
   -point_obs obs_stations.nc \
   -grid_obs gridded_obs.nc \
   EnsembleStatConfig \
   -outdir /output/ensemble_stat \
   -v 2

# Or with a file list
ls /data/ens/mem*.grb2 > ens_file_list.txt
ensemble_stat \
   -ens ens_file_list.txt \
   -point_obs obs_stations.nc \
   EnsembleStatConfig \
   -outdir /output/ensemble_stat
```

#### Key Config

```c
// EnsembleStatConfig

model = "GEFS";
obtype = "ADPSFC";

////////////////////////////////////////////////////////////////////////////////
// Ensemble fields

ens = {
   ens_thresh = 0.75;    // Require 75% of members present
   vld_thresh = 1.0;

   field = [
      { name = "TMP"; level = "Z2"; },
      { name = "APCP"; level = "A06";
        cat_thresh = [ >0.0, >=5.0, >=10.0, >=25.0 ]; }
   ];
}

////////////////////////////////////////////////////////////////////////////////
// Forecast fields (for ensemble mean / probability verification)

fcst = {
   field = [
      { name = "TMP"; level = "Z2"; },
      { name = "APCP"; level = "A06";
        cat_thresh = [ >0.0, >=5.0, >=10.0, >=25.0 ]; }
   ];
}

obs = fcst;

////////////////////////////////////////////////////////////////////////////////
// Climatology for CRPSS

climo_mean = {
   file_name = [ "/data/climo/climo_mean.nc" ];
   field     = [];
   regrid    = { method = BILINEAR; width = 2; vld_thresh = 0.5; }
}

////////////////////////////////////////////////////////////////////////////////
// Output options

output_flag = {
   ecnt  = BOTH;    // Ensemble Continuous statistics (CRPS, CRPSS, spread, ME)
   rps   = BOTH;    // Ranked Probability Score
   rhist = BOTH;    // Rank Histogram
   phist = BOTH;    // PIT Histogram
   orank = BOTH;    // Observation Rank (individual matched pairs)
   ssvar = BOTH;    // Spread-Skill Variance
   relp  = BOTH;    // Relative Position
}

////////////////////////////////////////////////////////////////////////////////
// Ensemble product output (write mean, spread, probability as NetCDF)

ensemble_flag = {
   latlon    = TRUE;
   mean      = TRUE;     // Ensemble mean
   stdev     = TRUE;     // Ensemble standard deviation (spread)
   minus     = FALSE;
   plus      = FALSE;
   min       = FALSE;
   max       = FALSE;
   range     = FALSE;
   vld_count = TRUE;
   frequency = TRUE;     // Probability of exceeding thresholds
   nep       = FALSE;    // Neighborhood Ensemble Probability
   nmep      = FALSE;    // Neighborhood Maximum Ensemble Probability
}
```

#### Key Output Statistics

| Line Type | Statistics |
|---|---|
| **ECNT** | CRPS, CRPSS, IGN, ME, RMSE, SPREAD, SPREAD_OERR |
| **RPS** | RPS, RPSS (Ranked Probability Score / Skill Score) |
| **RHIST** | Rank histogram counts (should be flat for calibrated ensemble) |
| **PHIST** | PIT histogram (Probability Integral Transform) |
| **SSVAR** | Spread-Skill variance relationship by spread bin |
| **RELP** | Relative position histogram |

### 5.4 Stat-Analysis

**Purpose:** Post-process MET `.stat` output files. Aggregates partial sums, computes derived statistics, filters by time/region/level, and generates confidence intervals.

This is the **key tool** for turning raw MET output into publication-quality statistics.

#### Command

```bash
# Filter and aggregate
stat_analysis \
   -lookin /output/point_stat/ \
   -job aggregate_stat \
   -line_type SL1L2 \
   -out_line_type CNT \
   -by FCST_VAR,FCST_LEV \
   -out_stat /output/agg_cnt.stat \
   -v 2

# Using a config file
stat_analysis \
   -lookin /output/point_stat/ \
   -config StatAnalysisConfig \
   -v 2
```

#### Common Jobs

```c
// StatAnalysisConfig

jobs = [

   // 1. Aggregate continuous partial sums → continuous stats
   "-job aggregate_stat -line_type SL1L2 -out_line_type CNT "
   "-fcst_var TMP -fcst_lev Z2 "
   "-by FCST_LEAD "
   "-out_stat agg_tmp2m_by_lead.stat",

   // 2. Filter and dump specific line types
   "-job filter -line_type CTS "
   "-fcst_var APCP -fcst_lev A24 "
   "-dump_row filtered_precip_cts.stat",

   // 3. Aggregate CTC (Contingency Table Counts) → CTS (Stats)
   "-job aggregate_stat -line_type CTC -out_line_type CTS "
   "-fcst_var APCP -fcst_lev A24 "
   "-fcst_thresh >=25.4 "
   "-by FCST_LEAD",

   // 4. Compute anomaly correlation from SAL1L2
   "-job aggregate_stat -line_type SAL1L2 -out_line_type CNT "
   "-fcst_var HGT -fcst_lev P500 "
   "-by FCST_LEAD",

   // 5. Bootstrap confidence intervals
   "-job aggregate_stat -line_type SL1L2 -out_line_type CNT "
   "-fcst_var TMP -fcst_lev Z2 "
   "-n_boot 1000 -boot_rng mt19937 "
   "-out_stat agg_tmp2m_boot.stat"
];
```

#### Filtering Options

| Flag | Description | Example |
|---|---|---|
| `-model` | Filter by model name | `-model WRF_d01` |
| `-fcst_var` | Forecast variable | `-fcst_var TMP` |
| `-fcst_lev` | Forecast level | `-fcst_lev Z2` |
| `-fcst_lead` | Forecast lead time (HHMMSS) | `-fcst_lead 120000` |
| `-vx_mask` | Verification mask name | `-vx_mask CONUS` |
| `-line_type` | Line type to read | `-line_type CNT` |
| `-fcst_thresh` | Threshold filter | `-fcst_thresh >=25.4` |
| `-valid_beg` / `-valid_end` | Valid time range | `-valid_beg 20240101 -valid_end 20240131` |
| `-init_beg` / `-init_end` | Init time range | `-init_beg 20240101` |
| `-by` | Group results by column | `-by FCST_LEAD,VX_MASK` |
| `-out_stat` | Output file | `-out_stat result.stat` |

#### Aggregation: Why Partial Sums Matter

MET tools write **partial sums** (SL1L2, SAL1L2, VL1L2, CTC) in addition to derived statistics. Partial sums are essential because:

1. They can be **aggregated** across times, regions, or models correctly
2. Direct averaging of RMSE, correlation, or CSI gives incorrect results
3. `stat_analysis` aggregates partial sums first, then derives the final statistics

```
# Wrong: averaging RMSE values directly
RMSE_wrong = mean(RMSE_day1, RMSE_day2, ..., RMSE_dayN)

# Correct: aggregate partial sums, then compute RMSE
# stat_analysis aggregates FBAR, OBAR, FFBAR, FOBAR, OOBAR from SL1L2
# then computes RMSE = sqrt(FFBAR - 2*FOBAR + OOBAR)
```

---

## 6. Spatial & Object-Based Verification

### 6.1 MODE (Method for Object-based Diagnostic Evaluation)

**Purpose:** Identify, match, and compare spatial features (objects) in forecast and observation fields. Addresses the "double penalty" problem of traditional grid-point verification.

Traditional verification penalizes a displaced forecast twice — once for missing at the observed location, once for falsely forecasting at the displaced location. MODE identifies objects and evaluates location error, size error, and intensity error separately.

#### How MODE Works

1. **Threshold / Convolve** — Apply a convolution radius + threshold to identify objects
2. **Identify Objects** — Connected components above threshold become objects
3. **Merge** — Nearby objects are optionally merged
4. **Match** — Forecast and observed objects are matched using an interest function
5. **Compute Attributes** — Area, centroid distance, angle difference, intensity ratio, etc.

#### Command

```bash
mode \
   forecast.grb2 \
   obs_analysis.nc \
   MODEConfig \
   -outdir /output/mode \
   -v 2
```

#### Key Config

```c
// MODEConfig

model = "WRF_d01";
obtype = "STAGE4";

////////////////////////////////////////////////////////////////////////////////
// Field specification

fcst = {
   field = {
      name       = "APCP";
      level      = "A06";
   };

   // Convolution and thresholding for object identification
   conv_radius = 5;            // Grid points — smoothing radius
   conv_thresh = >=5.0;        // Threshold after convolution (mm)
   merge_thresh = >=2.5;       // Threshold for merging nearby objects
   merge_flag  = THRESH;       // NONE, THRESH, ENGINE, BOTH
}

obs = {
   field = {
      name       = "APCP_06";
      level      = "(*,*)";
   };

   conv_radius = 5;
   conv_thresh = >=5.0;
   merge_thresh = >=2.5;
   merge_flag  = THRESH;
}

////////////////////////////////////////////////////////////////////////////////
// Matching and merging

match_flag  = MERGE_BOTH;    // NONE, MERGE_BOTH, MERGE_FCST, NO_MERGE
max_centroid_dist = 500.0;   // km — max distance for potential matches

////////////////////////////////////////////////////////////////////////////////
// Interest function weights for matching

weight = {
   centroid_dist    = 5.0;
   boundary_dist    = 0.0;
   convex_hull_dist = 0.0;
   angle_diff       = 2.0;
   aspect_diff      = 0.0;
   area_ratio       = 3.0;
   int_area_ratio   = 0.0;
   curvature_ratio  = 0.0;
   complexity_ratio = 0.0;
   inten_perc_ratio = 0.0;
   inten_perc_value = 50;
}

total_interest_thresh = 0.7;  // Minimum interest for a match

////////////////////////////////////////////////////////////////////////////////
// Output

ps_plot_flag  = TRUE;    // PostScript plots
ct_stats_flag = TRUE;    // Contingency table statistics on raw field
```

#### MODE Object Attributes

| Attribute | Description |
|---|---|
| CENTROID_X, CENTROID_Y | Object centroid position (grid coordinates) |
| CENTROID_LAT, CENTROID_LON | Centroid in geographic coordinates |
| AREA | Object area (grid squares) |
| INTENSITY_nn | nth percentile of intensity within object |
| CURVATURE | Object boundary curvature |
| COMPLEXITY | Ratio of perimeter² to area (shape complexity) |
| ANGLE | Orientation angle of fitted ellipse |
| ASPECT_RATIO | Ratio of minor to major axis of fitted ellipse |
| CENTROID_DIST | Distance between matched centroids (km) |
| AREA_RATIO | Ratio of forecast to observed area |
| INTEREST | Total interest score for the match |

### 6.2 MTD (MODE Time Domain)

**Purpose:** Extends MODE to 3D (x, y, time) — tracks objects through time and evaluates temporal evolution.

MTD identifies space-time objects across a sequence of 2D fields and matches forecast objects with observed objects.

#### Command

```bash
mtd \
   -fcst fcst_01.grb2 fcst_02.grb2 fcst_03.grb2 fcst_04.grb2 \
   -obs  obs_01.nc obs_02.nc obs_03.nc obs_04.nc \
   MTDConfig \
   -outdir /output/mtd \
   -v 2
```

#### Key Config

```c
// MTDConfig

model = "WRF_d01";

fcst = {
   field = {
      name  = "APCP";
      level = "A01";
   };

   conv_radius = 5;
   conv_thresh = >=2.0;
}

obs = {
   field = {
      name  = "APCP_01";
      level = "(*,*)";
   };

   conv_radius = 5;
   conv_thresh = >=2.0;
}

// Min volume (grid squares × time steps) for an object
min_volume = 2000;

// Time domain specific
total_interest_thresh = 0.7;

weight = {
   space_centroid_dist = 4.0;
   time_centroid_delta = 3.0;
   speed_delta         = 2.0;
   direction_diff      = 2.0;
   volume_ratio        = 1.0;
   axis_diff           = 1.0;
}
```

#### MTD Additional Attributes

| Attribute | Description |
|---|---|
| VOLUME | Space-time volume of object (grid squares × time steps) |
| START_TIME | First time step where object appears |
| END_TIME | Last time step |
| CDIST_TRAVELLED | Distance travelled by centroid |
| SPEED | Mean object speed (km/time step) |
| DIRECTION | Mean object direction |
| TIME_CENTROID_DELTA | Time offset between matched object centroids |

### 6.3 Grid-Diag

**Purpose:** Compute joint frequency distributions between two or more gridded fields. Useful for diagnosing relationships between forecast errors and environmental conditions.

#### Command

```bash
grid_diag \
   -data field1.nc field2.nc \
   GridDiagConfig \
   -outdir /output/grid_diag \
   -v 2
```

#### Key Config

```c
// GridDiagConfig

data = {
   field = [
      { name = "TMP"; level = "Z2"; },
      { name = "RH";  level = "Z2"; }
   ];
}

// Bins for joint distribution
mask = { grid = [ "FULL" ]; }
```

---

## 7. Tropical Cyclone Verification

MET includes a complete suite for tropical cyclone (TC) verification using the **ATCF (Automated Tropical Cyclone Forecasting)** data format.

### ATCF Data Format

The standard TC track format used by most TC forecast centers:

```
# A-deck (forecast): basin, cyclone number, date, tech number/name, forecast period,
#   lat, lon, max wind, min pressure, ...
AL, 09, 2024082912, 03, AVNO, 000, 261N,  806W,  65, 987, ...
AL, 09, 2024082912, 03, AVNO, 012, 268N,  812W,  70, 983, ...
AL, 09, 2024082912, 03, AVNO, 024, 278N,  820W,  80, 975, ...

# B-deck (best track / observations):
AL, 09, 2024082912,   , BEST,   0, 261N,  806W,  65, 987, HU, ...
AL, 09, 2024083000,   , BEST,   0, 268N,  813W,  75, 980, HU, ...
```

### 7.1 TC-Pairs

**Purpose:** Match forecast TC tracks with observed (best track) tracks. Computes track and intensity differences at each lead time.

#### Command

```bash
tc_pairs \
   -adeck /data/tc/adeck/a*.dat \
   -bdeck /data/tc/bdeck/b*.dat \
   -config TCPairsConfig \
   -out /output/tc/tc_pairs \
   -v 2
```

#### Key Config

```c
// TCPairsConfig

model = [ "AVNO", "HWRF", "HMON" ];

// Storm selection
storm_id   = [];           // Empty = all storms
basin      = [ "AL" ];     // Atlantic
cyclone    = [];           // All cyclones
storm_name = [];

// Time matching
init_beg = "20240801";
init_end = "20241130";
valid_beg = "";
valid_end = "";

// Consensus models (optional)
consensus = [];

// Lead times to process (hours)
lead_req = [];     // Empty = all leads

// Matching criteria
match_points = TRUE;

// Diagnostics
diag_info_map = [];
diag_convert_map = [];

// Output the TCMPR line type
check_dup = TRUE;
interp12  = REPLACE;
```

### 7.2 TC-Stat

**Purpose:** Analyze TC-Pairs output. Filter, aggregate, and compute summary statistics for track error, intensity error, along/cross-track error, etc.

#### Command

```bash
tc_stat \
   -lookin /output/tc/tc_pairs*.tcst \
   -config TCStatConfig \
   -out /output/tc/tc_stat_summary.tcst \
   -v 2
```

#### Key Config

```c
// TCStatConfig

jobs = [
   // Summary of track error by lead time for each model
   "-job summary -line_type TCMPR "
   "-column TK_ERR -by AMODEL,LEAD "
   "-out /output/tc/track_error_by_lead.stat",

   // Intensity (max wind) error
   "-job summary -line_type TCMPR "
   "-column AMAX_WIND-BMAX_WIND -by AMODEL,LEAD "
   "-out /output/tc/intensity_error_by_lead.stat",

   // Along-track and cross-track errors
   "-job summary -line_type TCMPR "
   "-column ALTK_ERR,CRTK_ERR -by AMODEL,LEAD "
   "-out /output/tc/along_cross_track.stat",

   // Filter to specific intensity category
   "-job filter -line_type TCMPR "
   "-column BMAX_WIND -column_thresh >=64 "
   "-dump_row /output/tc/hurricane_only.tcst",

   // Rapid intensification events
   "-job rirw -line_type TCMPR "
   "-rirw_track BDECK "
   "-rirw_time 24 "
   "-rirw_exact FALSE "
   "-rirw_thresh >=30 "
   "-by AMODEL "
   "-out /output/tc/ri_stats.stat"
];
```

#### Key TC Statistics

| Statistic | Description |
|---|---|
| **TK_ERR** | Great-circle track error (km) |
| **ALTK_ERR** | Along-track error (km, + = slow bias) |
| **CRTK_ERR** | Cross-track error (km, + = right of track) |
| **AMAX_WIND** | Forecast maximum wind (kt) |
| **BMAX_WIND** | Best track maximum wind (kt) |
| **AMSLP** | Forecast MSLP (hPa) |
| **BMSLP** | Best track MSLP (hPa) |
| **TK_ERR** by lead | Track error growth with lead time |

### 7.3 TC-Gen

**Purpose:** Verify tropical cyclogenesis forecasts. Evaluates whether and when genesis was predicted relative to observations.

#### Command

```bash
tc_gen \
   -genesis /data/tc/genesis_forecasts.dat \
   -track /data/tc/bdeck/*.dat \
   -config TCGenConfig \
   -out /output/tc/tc_gen \
   -v 2
```

#### Key Config

```c
// TCGenConfig

model = [ "AVNO" ];
init_freq = 6;           // hours between init times

// Genesis matching criteria
genesis_match_radius = 500;    // km
genesis_match_window = {
   beg = -24;                  // hours before best track genesis
   end =  24;                  // hours after
}

// Lead time window
fcst_hr_window = {
   beg = 0;
   end = 120;
}

// Output
ops_hit     = TRUE;
ops_miss    = TRUE;
ops_false_alarm = TRUE;

// Scoring
dev_hit_radius  = 500;
dev_hit_window = {
   beg = -24;
   end =  24;
}
```

### 7.4 TC-RMW

**Purpose:** Compute azimuthally averaged radial profiles of TC fields (wind, rain, etc.) relative to the radius of maximum wind (RMW).

#### Command

```bash
tc_rmw \
   -data forecast.grb2 \
   -deck /data/tc/adeck.dat \
   -config TCRMWConfig \
   -out /output/tc/tc_rmw_output.nc \
   -v 2
```

#### Key Config

```c
// TCRMWConfig

model = "HWRF";

// Radial grid
n_range  = 100;        // Number of radial bins
n_azimuth = 180;       // Number of azimuthal bins
delta_range_km = 10.0; // Radial bin width (km)

// Fields to compute profiles
data = {
   field = [
      { name = "UGRD"; level = "P850"; },
      { name = "VGRD"; level = "P850"; },
      { name = "TMP";  level = "P850"; },
      { name = "APCP"; level = "A06"; }
   ];
}

// RMW source
rmw_source = BDECK;    // Use best track for RMW
```

---

## 8. Specialized Tools

### 8.1 PCP-Combine

**Purpose:** Combine, accumulate, subtract, or derive fields from gridded data. Essential preprocessing step for precipitation verification.

#### Accumulation Mode

```bash
# Sum 6-hourly precip to create 24-hour accumulation
pcp_combine \
   -sum 20240115_000000 6 20240116_000000 24 \
   -pcpdir /data/forecasts \
   -pcprx "wrf_d01_.*\.grb2" \
   /output/wrf_apcp24h_20240116.nc

# Explicit file listing
pcp_combine \
   -sum 00000000_000000 6 00000000_000000 24 \
   /data/fcst_f006.grb2 \
   /data/fcst_f012.grb2 \
   /data/fcst_f018.grb2 \
   /data/fcst_f024.grb2 \
   /output/apcp_24h.nc
```

#### Subtraction Mode

```bash
# Derive 6-h precip from running totals (e.g., WRF RAINNC)
# apcp_f012 - apcp_f006 = 6-h precip valid at f012
pcp_combine \
   -subtract \
   /data/fcst_f012.grb2 'name="APCP"; level="A12";' \
   /data/fcst_f006.grb2 'name="APCP"; level="A06";' \
   /output/apcp_06h_f012.nc
```

#### Derive Mode

```bash
# Compute wind speed from U and V
pcp_combine \
   -derive WIND \
   /data/forecast.grb2 \
   'name="UGRD"; level="Z10";' \
   'name="VGRD"; level="Z10";' \
   /output/wind_speed_10m.nc

# Other derive operations: SUM, MIN, MAX, RANGE, MEAN, STDEV, VWIND, WIND
```

#### User-Defined Derivation

```bash
# Custom math using -derive with user-defined names
pcp_combine \
   -add \
   /data/convective_rain.grb2 'name="ACPCP"; level="A06";' \
   /data/nonconv_rain.grb2 'name="NCPCP"; level="A06";' \
   /output/total_precip_06h.nc \
   -name "APCP_06"
```

### 8.2 Regrid-Data-Plane

**Purpose:** Regrid data from one grid projection/resolution to another. Crucial for aligning forecast and observation grids.

#### Command

```bash
# Regrid to match a reference file's grid
regrid_data_plane \
   input.grb2 \
   reference_grid.grb2 \
   output.nc \
   -field 'name="TMP"; level="Z2";' \
   -method BILINEAR \
   -width 2 \
   -v 2

# Regrid to a specific grid definition string
regrid_data_plane \
   input.grb2 \
   "latlon 720 361 -90.0 0.0 0.5 0.5" \
   output.nc \
   -field 'name="TMP"; level="P850";' \
   -method BILINEAR \
   -width 2

# Regrid to a named grid
regrid_data_plane \
   input.grb2 \
   G212 \
   output.nc \
   -field 'name="HGT"; level="P500";' \
   -method BILINEAR \
   -width 2
```

#### Interpolation Methods

| Method | Description | Best For |
|---|---|---|
| `NEAREST` | Nearest neighbor | Categorical fields, land/sea mask |
| `BILINEAR` | Bilinear interpolation | Smooth continuous fields (temperature, height) |
| `BUDGET` | Budget interpolation (area-weighted) | Precipitation, fluxes |
| `UW_MEAN` | Unweighted mean of surrounding points | General smoothing |
| `DW_MEAN` | Distance-weighted mean | General purpose |
| `MIN` | Minimum value | Conservative estimates |
| `MAX` | Maximum value | Severe weather thresholds |
| `MEDIAN` | Median value | Noise-resistant interpolation |

#### Common Grid Specification Strings

```bash
# Lat-Lon grid: "latlon Nx Ny lat_ll lon_ll dlat dlon"
"latlon 720 361 -90.0 0.0 0.5 0.5"           # 0.5-deg global
"latlon 1440 721 -90.0 0.0 0.25 0.25"         # 0.25-deg global
"latlon 100 100 -10.0 95.0 0.1 0.1"           # 0.1-deg Indonesia domain

# Lambert Conformal: "lambert Nx Ny lat_ll lon_ll lon_orient d_km ..."
# Stereographic: "stereo Nx Ny ..."
# Mercator: "mercator Nx Ny ..."

# Named NCEP grids
# G212 = CONUS 40-km Lambert Conformal
# G218 = CONUS 12-km Lambert Conformal
# G004 = Global 0.5-deg
```

### 8.3 Gen-Vx-Mask

**Purpose:** Create verification masks defining where statistics should be computed. Essential for regional verification.

#### Polyline Mask

```bash
# Mask from a polygon file
gen_vx_mask \
   reference_grid.grb2 \
   /masks/indonesia.poly \
   indonesia_mask.nc \
   -type poly \
   -v 2

# Polygon file format (indonesia.poly):
# INDONESIA
#   95.0  -11.0
#  141.0  -11.0
#  141.0    6.0
#   95.0    6.0
#   95.0  -11.0
```

#### Circle Mask

```bash
# Circular mask centered on a point (e.g., 200 km around Jakarta)
gen_vx_mask \
   reference_grid.grb2 \
   "106.85 -6.2" \
   jakarta_200km_mask.nc \
   -type circle \
   -thresh <=200 \
   -v 2
```

#### Grid Mask

```bash
# Mask based on a grid (e.g., land/sea mask from another file)
gen_vx_mask \
   reference_grid.grb2 \
   land_sea_mask.nc \
   land_only_mask.nc \
   -type data \
   -input_field 'name="LAND"; level="L0";' \
   -thresh ==1 \
   -v 2
```

#### Combining Masks

```bash
# Intersection: apply second mask within first
gen_vx_mask \
   indonesia_mask.nc \
   land_sea_mask.nc \
   indonesia_land_mask.nc \
   -type data \
   -input_field 'name="LAND"; level="L0";' \
   -thresh ==1 \
   -intersection \
   -v 2

# Union: combine two masks
gen_vx_mask \
   mask1.nc \
   mask2.nc \
   combined_mask.nc \
   -type data \
   -input_field 'name="mask2_var"; level="(*,*)";' \
   -thresh >0 \
   -union \
   -v 2
```

#### Station Mask

```bash
# Mask using a station ID list (for Point-Stat)
gen_vx_mask \
   reference_grid.grb2 \
   station_list.txt \
   station_mask.nc \
   -type sid \
   -v 2

# station_list.txt format:
# Station_List
# WIII
# WICC
# WIDD
# WIMM
# WIPP
```

### 8.4 Observation Converters

#### pb2nc — PrepBUFR to NetCDF

```bash
pb2nc \
   prepbufr.gdas.20240115.t12z.nr \
   pb2nc_output.nc \
   PB2NCConfig \
   -v 2
```

```c
// PB2NCConfig

// Observation types to extract
message_type = [ "ADPSFC", "ADPUPA", "SFCSHP", "AIRCFT" ];

// Time window around reference time
obs_window = {
   beg = -5400;    // -1.5 hr
   end =  5400;    // +1.5 hr
}

// Station selection (empty = all)
station_id = [];

// Geographic domain
mask = {
   grid = "";
   poly = "";
}

// Level selection
level_range = {
   beg = 1;
   end = 255;
}

// Quality marks (0=best, 9=worst for NCEP)
quality_mark_thresh = 3;    // Include obs with quality mark ≤ 3

// Variables to process
obs_bufr_var = [ "TOB", "UOB", "VOB", "QOB", "POB", "ZOB", "D_RH" ];
```

#### ascii2nc — ASCII Text to NetCDF

```bash
ascii2nc \
   my_surface_obs.txt \
   ascii2nc_output.nc \
   -format met_point \
   -v 2

# Supported formats: met_point, little_r, surfrad, wwsis, aeronet, isd, python
```

#### madis2nc — MADIS to NetCDF

```bash
madis2nc \
   /data/madis/metar/20240115_1200.nc \
   madis2nc_output.nc \
   -type metar \
   -rec_beg 0 \
   -rec_end 100000 \
   -v 2

# Types: metar, raob, profiler, maritime, mesonet, acarsProfiles
```

### 8.5 Series-Analysis

**Purpose:** Compute verification statistics at each grid point over a time series of matched forecast-observation pairs. Produces a spatial map of statistics.

#### Command

```bash
series_analysis \
   -fcst fcst_list.txt \
   -obs obs_list.txt \
   -paired \
   -config SeriesAnalysisConfig \
   -out /output/series_analysis.nc \
   -v 2
```

#### Key Config

```c
// SeriesAnalysisConfig

fcst = {
   field = [
      { name = "TMP"; level = "Z2"; }
   ];
}

obs = fcst;

// Statistics to compute at each grid point
output_stats = {
   cnt  = [ "ME", "RMSE", "PR_CORR" ];
   ctc  = [];
   cts  = [ "CSI", "GSS" ];
   sl1l2 = [];
}

// Block size for processing (memory management)
block_size = 1024;

// Bootstrapping
boot = {
   interval = PCTILE;
   rep_prop = 1.0;
   n_rep    = 0;          // 0 = no bootstrapping
   seed     = "";
}
```

**Output:** A NetCDF file with spatial maps of the requested statistics (e.g., a grid of RMSE values showing where the model performs well or poorly).

### 8.6 Wavelet-Stat

**Purpose:** Scale-decomposed verification using 2D Haar wavelet decomposition. Shows at which spatial scales the forecast has skill.

#### Command

```bash
wavelet_stat \
   forecast.grb2 \
   observation.nc \
   WaveletStatConfig \
   -outdir /output/wavelet_stat \
   -v 2
```

#### Key Config

```c
// WaveletStatConfig

fcst = {
   field = {
      name  = "APCP";
      level = "A24";
      cat_thresh = [ >=5.0 ];
   };
}

obs = {
   field = {
      name  = "APCP_24";
      level = "(*,*)";
      cat_thresh = [ >=5.0 ];
   };
}

// Wavelet decomposition
grid_decomp_flag = AUTO;    // AUTO or TILE
tile = {
   width = 0;
   location = [];
}

wavelet = {
   type   = HAAR;           // Currently only HAAR supported
   member = 2;
}
```

### 8.7 IODA2NC

**Purpose:** Convert observations in JEDI's IODA (Interface for Observation Data Access) format to MET's internal NetCDF point observation format.

```bash
ioda2nc \
   /data/jedi/ioda_obs_surface.nc \
   ioda2nc_output.nc \
   IODA2NCConfig \
   -v 2
```

```c
// IODA2NCConfig

// Message type mapping
message_type_map = [
   { key = "surface_obs"; val = "ADPSFC"; },
   { key = "sonde_obs";   val = "ADPUPA"; }
];

// Observation variables to process
obs_var = [ "air_temperature", "specific_humidity", "wind_speed" ];

// Quality control
quality_mark_thresh = 0;

// Time window
obs_window = {
   beg = -5400;
   end =  5400;
}
```

---

## 9. METplus Wrappers & Automation

### 9.1 Wrapper Architecture

METplus provides **Python wrappers** for each MET tool, plus additional utility wrappers:

| Wrapper | MET Tool | Purpose |
|---|---|---|
| `PointStat` | `point_stat` | Point verification |
| `GridStat` | `grid_stat` | Gridded verification |
| `EnsembleStat` | `ensemble_stat` | Ensemble verification |
| `MODE` | `mode` | Object-based verification |
| `MTD` | `mtd` | MODE Time Domain |
| `SeriesAnalysis` | `series_analysis` | Time series at grid points |
| `StatAnalysis` | `stat_analysis` | Statistical aggregation |
| `TCPairs` | `tc_pairs` | TC track matching |
| `TCStat` | `tc_stat` | TC statistics |
| `TCGen` | `tc_gen` | TC genesis verification |
| `PcpCombine` | `pcp_combine` | Precip accumulation |
| `RegridDataPlane` | `regrid_data_plane` | Regridding |
| `GenVxMask` | `gen_vx_mask` | Mask generation |
| `PB2NC` | `pb2nc` | PrepBUFR conversion |
| `ASCII2NC` | `ascii2nc` | ASCII conversion |
| `WaveletStat` | `wavelet_stat` | Scale-decomposed verification |
| `GridDiag` | `grid_diag` | Joint distributions |
| `IODA2NC` | `ioda2nc` | JEDI obs conversion |
| `UserScript` | — | Run custom scripts |
| `CyclonePlotter` | — | Plot TC tracks |
| `MakePlots` | — | Generate standard plots |

### 9.2 Running Use Cases

METplus ships with extensive use cases organized by category.

```bash
# List available use cases
ls ${METPLUS_PATH}/parm/use_cases/

# Categories:
#   met_tool_wrapper/    — Single-tool examples
#   model_applications/  — Multi-tool real-world workflows
#   s2s/                 — Sub-seasonal to seasonal
#   marine_and_cryosphere/
#   medium_range/
#   precipitation/
#   convection_allowing_models/

# Run a use case
run_metplus.py \
   ${METPLUS_PATH}/parm/use_cases/met_tool_wrapper/GridStat/GridStat.conf \
   user.conf
```

### 9.3 Custom Use Case Creation

Create a configuration file that chains multiple tools:

```ini
# File: wrf_sfc_verification.conf

[config]
# Chain tools in order
PROCESS_LIST = PB2NC, PcpCombine, PointStat, StatAnalysis

# Time configuration
LOOP_BY = VALID
VALID_TIME_FMT = %Y%m%d%H
VALID_BEG = 2024011500
VALID_END = 2024013100
VALID_INCREMENT = 6H

LEAD_SEQ = 6, 12, 18, 24, 36, 48, 72

# Model info
MODEL = WRF_d01
OBTYPE = GDAS

# ===== PB2NC settings =====
PB2NC_INPUT_DIR = /data/prepbufr
PB2NC_INPUT_TEMPLATE = prepbufr.gdas.{valid?fmt=%Y%m%d}.t{valid?fmt=%H}z.nr
PB2NC_OUTPUT_DIR = {OUTPUT_BASE}/pb2nc
PB2NC_OUTPUT_TEMPLATE = pb2nc_{valid?fmt=%Y%m%d%H}.nc

PB2NC_MESSAGE_TYPE = ADPSFC, SFCSHP
PB2NC_OBS_WINDOW_BEGIN = -5400
PB2NC_OBS_WINDOW_END = 5400
PB2NC_QUALITY_MARK_THRESH = 3

# ===== PcpCombine settings (derive 6-h precip) =====
FCST_PCP_COMBINE_RUN = TRUE
FCST_PCP_COMBINE_METHOD = SUBTRACT
FCST_PCP_COMBINE_INPUT_DIR = /data/wrf
FCST_PCP_COMBINE_INPUT_TEMPLATE = wrfout_d01_{init?fmt=%Y-%m-%d_%H}:00:00.grb2
FCST_PCP_COMBINE_OUTPUT_DIR = {OUTPUT_BASE}/pcp_combine
FCST_PCP_COMBINE_OUTPUT_TEMPLATE = wrf_apcp06_{valid?fmt=%Y%m%d%H}.nc

# ===== Point-Stat settings =====
FCST_POINT_STAT_INPUT_DIR = /data/wrf
FCST_POINT_STAT_INPUT_TEMPLATE = wrfout_d01_{init?fmt=%Y-%m-%d_%H}:00:00.grb2
OBS_POINT_STAT_INPUT_DIR = {OUTPUT_BASE}/pb2nc
OBS_POINT_STAT_INPUT_TEMPLATE = pb2nc_{valid?fmt=%Y%m%d%H}.nc

POINT_STAT_OUTPUT_DIR = {OUTPUT_BASE}/point_stat

FCST_VAR1_NAME = TMP
FCST_VAR1_LEVELS = Z2
OBS_VAR1_NAME = TMP
OBS_VAR1_LEVELS = Z2
OBS_VAR1_OPTIONS = message_type = "ADPSFC";

FCST_VAR2_NAME = UGRD
FCST_VAR2_LEVELS = Z10
OBS_VAR2_NAME = UGRD
OBS_VAR2_LEVELS = Z10

FCST_VAR3_NAME = VGRD
FCST_VAR3_LEVELS = Z10
OBS_VAR3_NAME = VGRD
OBS_VAR3_LEVELS = Z10

POINT_STAT_MESSAGE_TYPE = ADPSFC
POINT_STAT_INTERP_TYPE_METHOD = BILINEAR
POINT_STAT_INTERP_TYPE_WIDTH = 2

POINT_STAT_OUTPUT_FLAG_CNT = BOTH
POINT_STAT_OUTPUT_FLAG_CTS = BOTH
POINT_STAT_OUTPUT_FLAG_SL1L2 = BOTH
POINT_STAT_OUTPUT_FLAG_VL1L2 = BOTH

# ===== Stat-Analysis settings =====
STAT_ANALYSIS_INPUT_DIR = {OUTPUT_BASE}/point_stat
STAT_ANALYSIS_OUTPUT_DIR = {OUTPUT_BASE}/stat_analysis

STAT_ANALYSIS_JOB1 = -job aggregate_stat -line_type SL1L2 -out_line_type CNT -by FCST_VAR,FCST_LEV,FCST_LEAD
STAT_ANALYSIS_JOB2 = -job aggregate_stat -line_type VL1L2 -out_line_type VCNT -by FCST_VAR,FCST_LEV,FCST_LEAD

[dir]
INPUT_BASE = /data
OUTPUT_BASE = /data/output

[filename_templates]
```

### 9.4 Running METplus

```bash
# Run with one or more config files (later configs override earlier ones)
run_metplus.py \
   wrf_sfc_verification.conf \
   user.conf

# Override individual settings from command line
run_metplus.py \
   wrf_sfc_verification.conf \
   user.conf \
   config.VALID_BEG=2024020100 \
   config.VALID_END=2024022800

# Dry run (show commands without executing)
run_metplus.py \
   wrf_sfc_verification.conf \
   user.conf \
   config.LOG_LEVEL=DEBUG
```

---

## 10. Visualization & Database Tools

### 10.1 METviewer

**Purpose:** MySQL-backed web application for interactive exploration of MET `.stat` output.

#### Setup

```bash
# Docker-based setup (recommended)
docker pull dtcenter/metviewer:latest

# Start MySQL container
docker run -d --name metviewer-mysql \
   -e MYSQL_ROOT_PASSWORD=metviewer \
   -e MYSQL_DATABASE=mv_database \
   -p 3306:3306 \
   mysql:8.0

# Start METviewer container
docker run -d --name metviewer \
   -p 8080:8080 \
   --link metviewer-mysql:mysql \
   -v /path/to/data:/data \
   dtcenter/metviewer:latest

# Access at http://localhost:8080/metviewer
```

#### Loading Data

```bash
# Use mv_load.sh to load .stat files into MySQL
docker exec metviewer /METviewer/bin/mv_load.sh \
   /data/load_spec.xml

# load_spec.xml example:
# <load_spec>
#   <connection>
#     <host>mysql:3306</host>
#     <database>mv_database</database>
#     <user>root</user>
#     <password>metviewer</password>
#   </connection>
#   <folder_tmpl>/data/output/{model}</folder_tmpl>
#   <load_val>
#     <field name="model"><val>WRF_d01</val></field>
#   </load_val>
# </load_spec>
```

#### METviewer Plot Types

| Plot Type | Description |
|---|---|
| Series | Statistics as a function of lead time, valid time, or other variable |
| Box | Box-and-whisker plots of distributions |
| Bar | Bar charts comparing models/thresholds |
| Scorecard | Red/green/yellow comparison matrix |
| Performance | Performance diagram (SR vs POD with CSI/bias contours) |
| Taylor | Taylor diagram (correlation vs std dev ratio) |
| ROC | Receiver Operating Characteristic curves |
| Reliability | Reliability diagrams for probabilistic forecasts |
| Contour | Contour plots of 2D statistics |
| Histogram | Rank histograms, PIT histograms |
| Ensemble SS | Ensemble spread-skill diagrams |

### 10.2 METplotpy

**Purpose:** Standalone Python plotting library for MET output. Use when you want more control than METviewer or need batch plotting.

#### Installation

```bash
pip install metplotpy
```

#### Example: Performance Diagram

```python
import metplotpy.plots.performance_diagram.performance_diagram as pd
import metcalcpy.util.read_env_vars_in_config as readconfig
import yaml

# Load config
with open('performance_diagram.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Create plot
diagram = pd.PerformanceDiagram(config)
diagram.save_to_file()
```

```yaml
# performance_diagram.yaml
stat_input: /output/stat_analysis/agg_precip_cts.stat
plot_filename: /output/plots/performance_diagram.png

series_val_1:
  model: [WRF_d01, GFS]

indy_var: FCST_THRESH
indy_vals: [">0.0", ">=2.54", ">=6.35", ">=12.7", ">=25.4"]

plot_ci: [boot]
n_boot: 1000

title: "24-h Precipitation Verification"
xaxis: "Success Ratio (1-FAR)"
yaxis: "Probability of Detection"
```

#### Example: Taylor Diagram

```python
import metplotpy.plots.taylor_diagram.taylor_diagram as td
import yaml

with open('taylor_diagram.yaml', 'r') as f:
    config = yaml.safe_load(f)

diagram = td.TaylorDiagram(config)
diagram.save_to_file()
```

#### Example: ROC Curve

```python
import metplotpy.plots.roc_diagram.roc_diagram as roc
import yaml

with open('roc_diagram.yaml', 'r') as f:
    config = yaml.safe_load(f)

diagram = roc.RocDiagram(config)
diagram.save_to_file()
```

### 10.3 METcalcpy

**Purpose:** Calculation engine underlying METplotpy. Can be used independently for custom analysis scripts.

```python
import metcalcpy.util.ctc_statistics as ctc

# Compute contingency table statistics from counts
# fy_oy = hits, fy_on = false alarms, fn_oy = misses, fn_on = correct negatives
stats = ctc.calculate_ctc_stats(
    fy_oy=150, fy_on=30, fn_oy=20, fn_on=800
)
print(f"CSI:  {stats['csi']:.3f}")
print(f"POD:  {stats['pod']:.3f}")
print(f"FAR:  {stats['far']:.3f}")
print(f"GSS:  {stats['gss']:.3f}")
print(f"BIAS: {stats['fbias']:.3f}")
```

### 10.4 METdataio

**Purpose:** Data I/O utilities for reading MET output files into pandas DataFrames.

```python
from metdataio.read_stat import parse_stat

# Read .stat file into a pandas DataFrame
df = parse_stat('/output/point_stat/point_stat_120000L_20240115_120000V.stat')

# Filter to CNT line type
cnt = df[df['LINE_TYPE'] == 'CNT']

# Extract specific statistics
me = cnt[cnt['FCST_VAR'] == 'TMP']['ME'].astype(float)
rmse = cnt[cnt['FCST_VAR'] == 'TMP']['RMSE'].astype(float)
```

### 10.5 METexpress

**Purpose:** Lightweight web interface for quick verification stat viewing. Simpler than METviewer, designed for operational monitoring.

```bash
# Docker setup
docker pull dtcenter/metexpress:latest
docker run -d -p 3000:3000 \
   --link metviewer-mysql:mysql \
   dtcenter/metexpress:latest

# Access at http://localhost:3000
```

METexpress provides pre-configured "apps" for common verification tasks:
- Upper Air (radiosonde verification)
- Surface (surface station verification)
- Anomaly Correlation (ACC by lead time)
- Precipitation
- Ensemble
- TC and more

---

## 11. Python Embedding

Python embedding allows MET to read **any** data format by calling a user-written Python script at runtime. The script passes data to MET via NumPy arrays, xarray Datasets, or pandas DataFrames.

### 11.1 Three Interfaces

| Interface | Keyword | Python Return Type | Use For |
|---|---|---|---|
| **PYTHON_NUMPY** | `file_type = PYTHON_NUMPY;` | dict with `met_data` (NumPy array) + attrs | Gridded data in any format |
| **PYTHON_XARRAY** | `file_type = PYTHON_XARRAY;` | `xarray.Dataset` | CF-compliant-style gridded data |
| **PYTHON_PANDAS** | `file_type = PYTHON_PANDAS;` | `pandas.DataFrame` | Point observations |

### 11.2 PYTHON_NUMPY Example — Custom Model Format

```python
#!/usr/bin/env python3
# File: read_custom_model.py
# Called by MET via Python embedding to read a custom binary format

import numpy as np
import sys

def main():
    # MET passes the filename as the first argument
    input_file = sys.argv[1]

    # Read your custom format
    data = np.fromfile(input_file, dtype=np.float32).reshape(100, 200)

    # Build the attrs dictionary MET expects
    attrs = {
        'valid':  '20240115_120000',
        'init':   '20240115_000000',
        'lead':   '120000',          # HHMMSS format
        'accum':  '060000',          # Accumulation period
        'name':   'APCP',
        'long_name': 'Total Precipitation',
        'level':  'A06',
        'units':  'mm',

        # Grid definition (lat-lon grid)
        'grid': {
            'type': 'LatLon',
            'name': 'custom_grid',
            'lat_ll':  -10.0,
            'lon_ll':   95.0,
            'delta_lat': 0.1,
            'delta_lon': 0.1,
            'Nlat': 100,
            'Nlon': 200,
        }
    }

    return {'met_data': data, 'attrs': attrs}

# MET calls this at import time
met_data = main()
```

#### Using in MET Config

```c
// GridStatConfig with Python embedding

fcst = {
   field = [
      {
         name = "PYTHON_NUMPY";
         // The script path + any arguments
         // MET calls: python3 read_custom_model.py /data/custom_fcst.bin
         level = "read_custom_model.py /data/custom_fcst.bin";
      }
   ];
}
```

### 11.3 PYTHON_XARRAY Example

```python
#!/usr/bin/env python3
# File: read_model_xarray.py

import xarray as xr
import sys

input_file = sys.argv[1]

# Read with xarray (e.g., a non-standard NetCDF)
ds = xr.open_dataset(input_file)

# Select and rename the variable MET expects
met_data = ds.rename({'temp_2m': 'TMP'})

# MET reads this variable
```

```c
// Config
fcst = {
   field = [
      {
         name = "PYTHON_XARRAY";
         level = "read_model_xarray.py /data/model_output.nc";
      }
   ];
}
```

### 11.4 PYTHON_PANDAS Example — Point Observations

```python
#!/usr/bin/env python3
# File: read_custom_obs.py
# Provides point observations to MET via pandas DataFrame

import pandas as pd
import sys

input_file = sys.argv[1]

# Read observations from custom format
raw = pd.read_csv(input_file)

# Build MET point observation DataFrame
# Required columns (in order):
# typ, sid, vld, lat, lon, elv, var, lvl, hgt, qc, obs
obs = pd.DataFrame({
    'typ': raw['station_type'],        # Message type (e.g., ADPSFC)
    'sid': raw['station_id'],          # Station ID
    'vld': raw['datetime'],            # Valid time (YYYYMMDD_HHMMSS)
    'lat': raw['latitude'],
    'lon': raw['longitude'],
    'elv': raw['elevation'],
    'var': raw['variable_name'],       # Variable name (TMP, UGRD, etc.)
    'lvl': raw['pressure_level'],      # Pressure level (hPa) or 0
    'hgt': raw['height_agl'],          # Height above ground (m)
    'qc':  raw['quality_flag'],        # QC string
    'obs': raw['value']                # Observation value
})

met_point_data = obs
```

```c
// PointStatConfig
obs = {
   field = [
      {
         name = "PYTHON_PANDAS";
         level = "read_custom_obs.py /data/surface_obs.csv";
         message_type = "ADPSFC";
      }
   ];
}
```

### 11.5 Environment Setup for Python Embedding

```bash
# MET must be compiled with Python embedding support
# Set these before compilation:
export MET_PYTHON=/usr/bin/python3
export MET_PYTHON_BIN_EXE=/usr/bin/python3
export MET_PYTHON_CC="-I/usr/include/python3.10"
export MET_PYTHON_LD="-L/usr/lib/x86_64-linux-gnu -lpython3.10"

# At runtime, ensure your script's dependencies are available:
export PYTHONPATH=/path/to/your/scripts:${PYTHONPATH}

# For conda environments:
export MET_PYTHON_EXE=$(which python3)
```

---

## 12. Practical Workflows & Examples

### 12.1 WRF Surface Verification Workflow

**Goal:** Verify WRF 2-m temperature and 10-m wind against surface stations.

```
WRF wrfout files (GRIB2 or NetCDF)
        |
        v
PrepBUFR obs  →  [pb2nc]  →  MET NetCDF obs
        |                           |
        v                           v
[point_stat]  ← forecast + obs matching
        |
        v
.stat files (SL1L2, CTC, CTS, CNT per time)
        |
        v
[stat_analysis]  →  aggregated stats by lead time
        |
        v
Plotting (METplotpy / custom matplotlib)
```

#### Step-by-Step Commands

```bash
# 1. Convert PrepBUFR observations
pb2nc \
   /data/prepbufr/gdas.20240115.t12z.prepbufr \
   /output/pb2nc/obs_20240115_12.nc \
   PB2NCConfig -v 2

# 2. Run Point-Stat for each valid time / lead
for lead in 006 012 018 024 036 048 072; do
   point_stat \
      /data/wrf/wrfprs_d01_f${lead}.grb2 \
      /output/pb2nc/obs_${valid_time}.nc \
      PointStatConfig \
      -outdir /output/point_stat \
      -v 2
done

# 3. Aggregate results
stat_analysis \
   -lookin /output/point_stat/ \
   -config StatAnalysisConfig_surface.conf \
   -v 2

# 4. Plot (example with matplotlib)
python3 plot_verification_scores.py \
   --input /output/stat_analysis/agg_surface.stat \
   --vars TMP,UGRD,VGRD \
   --output /output/plots/
```

### 12.2 Precipitation Verification (QPF vs Gauge/Radar)

**Goal:** Verify WRF 24-h precipitation against Stage IV radar-gauge analysis.

```bash
# 1. Accumulate model precipitation to 24-h totals
pcp_combine \
   -subtract \
   /data/wrf/wrfprs_d01_f024.grb2 'name="APCP"; level="A24";' \
   /data/wrf/wrfprs_d01_f000.grb2 'name="APCP"; level="A0";' \
   /output/pcp_combine/wrf_apcp24_20240116.nc

# 2. Create verification mask (e.g., Java Island)
gen_vx_mask \
   /output/pcp_combine/wrf_apcp24_20240116.nc \
   /masks/java.poly \
   /output/masks/java_mask.nc \
   -type poly -v 2

# 3. Run Grid-Stat with neighborhood verification
grid_stat \
   /output/pcp_combine/wrf_apcp24_20240116.nc \
   /data/obs/stage4_24h_20240116.nc \
   GridStatConfig_precip \
   -outdir /output/grid_stat \
   -v 2

# 4. Aggregate and analyze
stat_analysis \
   -lookin /output/grid_stat/ \
   -job aggregate_stat -line_type CTC -out_line_type CTS \
   -fcst_var APCP -by FCST_THRESH \
   -out_stat /output/precip_cts_by_thresh.stat
```

### 12.3 Ensemble Verification (GEFS)

**Goal:** Verify GEFS ensemble against surface observations (rank histogram, CRPS, spread-skill).

```bash
# 1. List ensemble members
ls /data/gefs/gep*.t00z.pgrb2a.0p50.f024 > /tmp/ens_list.txt
# Contains: gep01, gep02, ..., gep30

# 2. Convert observations
pb2nc \
   /data/prepbufr/gdas.20240116.t00z.prepbufr \
   /output/pb2nc/obs_20240116_00.nc \
   PB2NCConfig -v 2

# 3. Run Ensemble-Stat
ensemble_stat \
   -ens /tmp/ens_list.txt \
   -point_obs /output/pb2nc/obs_20240116_00.nc \
   EnsembleStatConfig \
   -outdir /output/ensemble_stat \
   -v 2

# 4. Check output
# /output/ensemble_stat/ensemble_stat_*.stat
# Contains: ECNT (CRPS, spread), RHIST (rank histogram), SSVAR, etc.

# 5. Extract rank histogram data
stat_analysis \
   -lookin /output/ensemble_stat/ \
   -job filter -line_type RHIST \
   -fcst_var TMP -fcst_lev Z2 \
   -dump_row /output/rhist_tmp2m.stat
```

### 12.4 Tropical Cyclone Verification

**Goal:** Verify GFS TC track forecasts against best track for the Atlantic basin.

```bash
# 1. Match forecast tracks with best tracks
tc_pairs \
   -adeck /data/tc/adeck/aal*.dat \
   -bdeck /data/tc/bdeck/bal*.dat \
   -config TCPairsConfig \
   -out /output/tc/tc_pairs_2024 \
   -v 2

# 2. Compute summary statistics
tc_stat \
   -lookin /output/tc/tc_pairs_2024*.tcst \
   -job summary -line_type TCMPR \
   -column TK_ERR \
   -by AMODEL,LEAD \
   -out /output/tc/track_error_summary.stat

# 3. Intensity error
tc_stat \
   -lookin /output/tc/tc_pairs_2024*.tcst \
   -job summary -line_type TCMPR \
   -column AMAX_WIND-BMAX_WIND \
   -by AMODEL,LEAD \
   -out /output/tc/intensity_error_summary.stat

# 4. Homogeneous comparison (only cases where all models have forecasts)
tc_stat \
   -lookin /output/tc/tc_pairs_2024*.tcst \
   -job summary -line_type TCMPR \
   -column TK_ERR \
   -by AMODEL,LEAD \
   -amodel AVNO,HWRF,HMON \
   -column_thresh TOTAL>=20 \
   -out /output/tc/homogeneous_track_error.stat
```

### 12.5 WFRT/verif Quick-Start (Lightweight Alternative)

For quick, lightweight verification without the full MET stack:

```bash
# Install
pip install verif

# Basic usage — verify a NetCDF forecast against observations
verif forecast.nc -m mae rmse corr -x leadtime

# Compare multiple models
verif model1.nc model2.nc -m mae rmse -x leadtime -leg "WRF" "GFS"

# Reliability diagram for probability forecasts
verif prob_forecast.nc -m reliability -r 0

# Performance diagram
verif forecast.nc -m performance -r 0.1

# Map of scores
verif forecast.nc -m mae -type map

# Multiple thresholds
verif forecast.nc -m ets -r 0.1,1,5,10,25
```

`verif` is great for exploratory analysis but lacks MET's spatial methods, TC tools, and aggregation capabilities.

---

## 13. Tips, Troubleshooting & Best Practices

### 13.1 Common Errors and Solutions

| Error | Cause | Solution |
|---|---|---|
| `ERROR: No data found` | File path or field name mismatch | Check `wgrib2 -s` or `ncdump -h` for exact field names |
| `ERROR: Could not open file` | File doesn't exist or wrong permissions | Verify path; check Docker bind mounts |
| `WARNING: No matching observations` | Time window too narrow or wrong message type | Widen `obs_window`; check `message_type` list |
| `ERROR: Gridded data has no coordinate info` | Non-CF-compliant NetCDF | Add grid info via config or use Python embedding |
| `ERROR: Mismatch in grid dimensions` | Forecast and obs on different grids | Use `regrid` block or `regrid_data_plane` first |
| `WARNING: 0 matched pairs` | No co-located forecast-obs pairs | Check valid times match; check spatial overlap |
| `ERROR: Cannot find Python` | Python embedding not compiled in | Recompile MET with Python support |
| `Segmentation fault` | Memory overflow with large grids | Increase memory; reduce domain; use tiles |

### 13.2 Debugging Tips

```bash
# Increase verbosity (1=minimal, 5=maximum)
point_stat forecast.grb2 obs.nc PointStatConfig -v 5

# Check what MET sees in a GRIB2 file
wgrib2 -s forecast.grb2 | grep TMP

# Check what MET sees in a NetCDF file
ncdump -h observation.nc

# Dump matched pairs for manual inspection
# Set in config: output_flag = { mpr = BOTH; }
# Then examine the MPR lines in the .stat output

# METplus: check the MET commands that were generated
cat ${OUTPUT_BASE}/logs/*.log | grep "Running command"
```

### 13.3 Performance Tuning

```bash
# For large datasets, control memory usage:
# 1. Process one field at a time instead of all at once
# 2. Use -compress option for output
# 3. Limit the verification domain with masks

# Parallel processing with METplus (process multiple valid times)
# METplus doesn't natively parallelize, but you can:
# a) Split time ranges and run multiple METplus instances
# b) Use GNU parallel:
cat valid_times.txt | parallel -j 8 \
   "run_metplus.py my_config.conf config.VALID_BEG={} config.VALID_END={}"

# c) Use SLURM job arrays:
#!/bin/bash
#SBATCH --array=0-30
DATES=($(seq -f "%Y%m%d" 20240101 1 20240131))
run_metplus.py my_config.conf config.VALID_BEG=${DATES[$SLURM_ARRAY_TASK_ID]}
```

### 13.4 Bootstrapping Confidence Intervals

MET supports bootstrapping for confidence intervals in most tools:

```c
// In any MET config file
boot = {
   interval = PCTILE;    // PCTILE or BCA (bias-corrected accelerated)
   rep_prop = 1.0;       // Proportion of sample for each replicate
   n_rep    = 1000;      // Number of bootstrap replicates (0=disable)
   seed     = "";        // Random seed (empty=random)
}
```

In `stat_analysis`:

```bash
stat_analysis \
   -lookin /output/point_stat/ \
   -job aggregate_stat -line_type SL1L2 -out_line_type CNT \
   -fcst_var TMP -fcst_lev Z2 \
   -n_boot 1000 \
   -boot_rng mt19937 \
   -out_stat agg_with_ci.stat
```

Output includes `*_BCL` (lower bound) and `*_BCU` (upper bound) columns for each statistic.

### 13.5 Climatology and Reference Forecasts

For skill scores (ACC, CRPSS, BSS, RPSS), you need a reference forecast (usually climatology):

```bash
# Common climatology sources:
# 1. NCEP Climate Forecast System Reanalysis (CFSR) climatology
# 2. ERA5 climatological means
# 3. Custom climatology computed from historical data

# Generate climatology from historical data
# (Using CDO — Climate Data Operators)
cdo ymonmean -selyear,1991/2020 historical_data.nc climo_mean.nc
cdo ymonstd  -selyear,1991/2020 historical_data.nc climo_stdev.nc
```

Then reference in MET config:

```c
climo_mean = {
   file_name = [ "/data/climo/climo_mean.nc" ];
   field     = [];
   regrid    = { method = BILINEAR; width = 2; vld_thresh = 0.5; }
   time_interp_method = DW_MEAN;
   day_interval = 31;
   hour_interval = 6;
}
```

### 13.6 Verification Best Practices

**From the companion theory notes — practical application in MET:**

| Principle | MET Implementation |
|---|---|
| **Use multiple metrics** — no single score tells the whole story | Enable multiple `output_flag` line types; use `stat_analysis` to compute CNT + CTS + NBRCTS together |
| **Stratify results** — score separately by region, season, intensity | Use `mask` in config; use `-by` in stat_analysis; use `-vx_mask`, `-valid_beg/-end` filters |
| **Use skill scores** — compare against a reference (climatology, persistence) | Configure `climo_mean` / `climo_stdev` to get ACC, BSS, CRPSS, RPSS |
| **Bootstrap confidence intervals** — determine if differences are significant | Set `n_rep > 0` in boot config; use `-n_boot` in stat_analysis |
| **Avoid double penalty** — use spatial methods for displaced features | Use MODE/MTD for object-based; use Grid-Stat neighborhood (FSS) for fuzzy |
| **Verify at appropriate scales** — don't expect grid-point accuracy for convection | Use Grid-Stat `nbrhd` with multiple widths; use Wavelet-Stat for scale decomposition |
| **Aggregate correctly** — use partial sums, not averaged statistics | Always aggregate SL1L2/CTC via stat_analysis, never average RMSE/CSI directly |
| **Homogeneous comparison** — compare models on the same set of cases | Use stat_analysis `-by` and `-column_thresh TOTAL>=N` to ensure common samples |

### 13.7 Quick Reference: Key Statistics and Their MET Line Types

| Statistic | Full Name | Line Type | Tool |
|---|---|---|---|
| **ME** | Mean Error (bias) | CNT | Point-Stat, Grid-Stat |
| **MAE** | Mean Absolute Error | CNT | Point-Stat, Grid-Stat |
| **RMSE** | Root Mean Square Error | CNT | Point-Stat, Grid-Stat |
| **PR_CORR** | Pearson Correlation | CNT | Point-Stat, Grid-Stat |
| **ACC** | Anomaly Correlation Coefficient | CNT (from SAL1L2) | Stat-Analysis |
| **CSI** | Critical Success Index (Threat Score) | CTS | Point-Stat, Grid-Stat |
| **GSS** | Gilbert Skill Score (Equitable Threat Score) | CTS | Point-Stat, Grid-Stat |
| **POD** | Probability of Detection (Hit Rate) | CTS | Point-Stat, Grid-Stat |
| **FAR** | False Alarm Ratio | CTS | Point-Stat, Grid-Stat |
| **HK** | Hanssen-Kuipers Discriminant (TSS) | CTS | Point-Stat, Grid-Stat |
| **HSS** | Heidke Skill Score | CTS | Point-Stat, Grid-Stat |
| **FBIAS** | Frequency Bias | CTS | Point-Stat, Grid-Stat |
| **FSS** | Fractions Skill Score | NBRCTS | Grid-Stat |
| **BSS** | Brier Skill Score | PSTD | Point-Stat, Grid-Stat |
| **CRPS** | Continuous Ranked Probability Score | ECNT | Ensemble-Stat |
| **CRPSS** | CRPS Skill Score | ECNT | Ensemble-Stat |
| **RPS** | Ranked Probability Score | RPS | Ensemble-Stat |
| **RPSS** | RPS Skill Score | RPS | Ensemble-Stat |
| **SPREAD** | Ensemble Spread | ECNT | Ensemble-Stat |
| **TK_ERR** | TC Track Error | TCMPR | TC-Stat |
| **S1** | S1 Gradient Score | GRAD | Grid-Stat |

---

## Appendix A: MET Tool Quick Reference

| Tool | Primary Input | Purpose |
|---|---|---|
| `point_stat` | Grid + point obs | Point verification |
| `grid_stat` | Grid + grid obs | Gridded verification + neighborhood |
| `ensemble_stat` | Ensemble members + obs | Ensemble verification |
| `stat_analysis` | `.stat` files | Aggregation, filtering, CI |
| `mode` | Grid + grid obs | Object-based verification |
| `mtd` | Time series of grids | 3D object tracking |
| `grid_diag` | Two+ grids | Joint distributions |
| `tc_pairs` | A-deck + B-deck | TC track matching |
| `tc_stat` | TC-Pairs output | TC statistics |
| `tc_gen` | Genesis + B-deck | Genesis verification |
| `tc_rmw` | Grid + deck | Radial wind profiles |
| `pcp_combine` | Gridded files | Accumulation, subtraction, derivation |
| `regrid_data_plane` | Gridded file | Regridding |
| `gen_vx_mask` | Grid + poly/stations | Create masks |
| `pb2nc` | PrepBUFR | Convert obs |
| `ascii2nc` | ASCII text | Convert obs |
| `madis2nc` | MADIS NetCDF | Convert obs |
| `ioda2nc` | JEDI IODA | Convert obs |
| `series_analysis` | Paired grids | Grid-point time series stats |
| `wavelet_stat` | Grid + grid obs | Scale-decomposed verification |
| `plot_data_plane` | Grid | Quick visualization of a field |
| `plot_point_obs` | Point obs NetCDF | Quick visualization of obs |
| `shift_data_plane` | Grid | Shift / rotate a field |

---

## Appendix B: METplus Environment Variables Reference

| Variable | Description | Example |
|---|---|---|
| `MET_INSTALL_DIR` | Path to compiled MET binaries | `/opt/met` |
| `METPLUS_PATH` | Path to METplus installation | `/opt/METplus` |
| `INPUT_BASE` | Base directory for input data | `/data/input` |
| `OUTPUT_BASE` | Base directory for output | `/data/output` |
| `LOG_LEVEL` | Logging level (DEBUG, INFO, WARNING, ERROR) | `DEBUG` |
| `MET_TMP_DIR` | Temporary file directory | `/tmp/met` |
| `METPLUS_CONF` | Auto-generated final config (for debugging) | `{OUTPUT_BASE}/metplus_final.conf` |
| `STAGING_DIR` | Directory for intermediate files | `{OUTPUT_BASE}/staging` |

---

## Appendix C: Useful External Resources

| Resource | URL |
|---|---|
| MET User's Guide | https://met.readthedocs.io/ |
| METplus User's Guide | https://metplus.readthedocs.io/ |
| DTC MET Homepage | https://dtcenter.org/community-code/model-evaluation-tools-met |
| METplus Use Cases | https://metplus.readthedocs.io/en/latest/Users_Guide/usecases.html |
| MET GitHub | https://github.com/dtcenter/MET |
| METplus GitHub | https://github.com/dtcenter/METplus |
| METviewer GitHub | https://github.com/dtcenter/METviewer |
| METplotpy GitHub | https://github.com/dtcenter/METplotpy |
| DTC Tutorial Videos | https://dtcenter.org/community-code/metplus/online-tutorial |
| WFRT verif | https://github.com/WFRT/verif |
| Verification Theory (companion) | [forecast-verification-notes.md](forecast-verification-notes.md) |
