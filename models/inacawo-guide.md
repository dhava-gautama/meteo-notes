# InaCAWO — Indonesia Coupled Atmosphere-Wave-Ocean Model Guide

InaCAWO is BMKG's operational coupled weather-ocean-wave forecasting system. It was developed by Baron Weather (USA) under the MMS-1 (Maritime Meteorological System) project, with CLS (France) and Atos/Bull (HPC). Operational since mid-2023 at BMKG headquarters in Jakarta.

This is one of the largest COAWST implementations in the world — fully coupled WRF + ROMS + SWAN across all Indonesian waters at 3 km resolution.

## What does InaCAWO actually do?

It runs three models simultaneously and makes them talk to each other every few minutes:

- **WRF** (atmosphere) computes wind, rain, pressure, temperature
- **ROMS** (ocean) computes currents, SST, salinity, water level
- **SWAN** (waves) computes wave height, period, direction, swell

The coupling matters because in reality these systems interact constantly — waves roughen the surface and change how wind transfers momentum, ocean temperature affects atmospheric convection, wind drives currents and waves, currents modify wave propagation. Running them separately misses all of this.

## System Architecture

### The COAWST Framework

InaCAWO is built on **COAWST** (Coupled Ocean-Atmosphere-Wave-Sediment Transport) from USGS. COAWST uses the **Model Coupling Toolkit (MCT)** to exchange data between models. All three models compile into a **single executable** — they run in parallel on separate processor sets and synchronize at defined intervals.

```
┌─────────────────────────────────────────────────────┐
│                    COAWST (MCT)                     │
│                                                     │
│   ┌─────────┐     ┌─────────┐     ┌─────────┐     │
│   │   WRF   │◄───►│  ROMS   │◄───►│  SWAN   │     │
│   │  (atm)  │     │ (ocean) │     │ (waves) │     │
│   └─────────┘     └─────────┘     └─────────┘     │
│        │               ▲               │           │
│        └───────────────┼───────────────┘           │
│              All exchange via MCT                   │
└─────────────────────────────────────────────────────┘
```

### What gets exchanged?

**WRF → ROMS:**
- Surface pressure, 2m temperature, 10m winds (U, V)
- Relative humidity, shortwave/longwave radiation, rain
- These drive the ocean surface forcing

**ROMS → WRF:**
- Sea surface temperature (SST)
- This updates the lower boundary condition for the atmosphere in real time — crucial for tropical convection

**ROMS → SWAN:**
- Depth-averaged currents, water surface elevation, bathymetry
- Currents modify wave propagation (wave-current interaction)

**SWAN → ROMS:**
- Wave height, direction, length, period, breaking percentage, energy dissipation
- Waves affect ocean mixing and bottom stress

**SWAN → WRF:**
- Sea surface roughness (modified by wave state)
- Young, steep waves make the surface rougher → changes wind profile

### Domain and Resolution

| Parameter | Value |
|-----------|-------|
| Horizontal resolution | 3 km (some sources say 2.5 km) |
| Domain extent | 90°E–145°E, 15°S–15°N |
| Coverage | All Indonesian waters + surrounding seas |
| Planned upgrade | 1 km for critical areas (Sunda Strait) |

This domain is enormous — roughly 55° longitude × 30° latitude at 3 km spacing. That's on the order of **2000 × 1000 grid points** for the atmosphere alone, which is why it needs a 1.31 Pflops supercomputer.

Whether WRF uses nesting (e.g., a coarser outer domain feeding a 3 km inner domain) or runs as a single 3 km domain is not publicly documented. Given the domain size, a single 3 km domain would be extremely expensive — nesting is likely but unconfirmed.

### Operational Schedule

| Parameter | Value |
|-----------|-------|
| Update cycles | 4x daily |
| Forecast length | 10 days |
| Output interval | 3-hourly |
| Workflow management | Baron ROME |

## Boundary Conditions

InaCAWO needs initial and boundary conditions from global models to start and to constrain the lateral boundaries throughout the forecast.

### Atmospheric BC (for WRF)

| Source | What it provides | Likely resolution |
|--------|-----------------|-------------------|
| **GFS** (NCEP) | Atmosphere IC/BC: wind, temperature, humidity, pressure on pressure levels | 0.25° / 3-hourly |
| **ECMWF** (IFS/HRES) | Alternative atmosphere IC/BC (higher quality, but needs license) | 0.1° or 0.25° |

How to download GFS (the free option):
```bash
# GFS 0.25-degree, grib2 format
# Available from NCEP NOMADS or AWS
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.YYYYMMDD/HH/atmos/gfs.tHHz.pgrb2.0p25.fFFF

# Or via AWS (faster, no quota)
aws s3 ls --no-sign-request s3://noaa-gfs-bdp-pds/gfs.YYYYMMDD/HH/atmos/
```

How to download ECMWF (requires MARS access or CDS API):
```bash
# Using ECMWF CDS API (Copernicus Climate Data Store)
pip install cdsapi
# Then use cdsapi.Client() to request operational forecasts or ERA5 for research
```

### Ocean BC (for ROMS)

| Source | What it provides | Resolution |
|--------|-----------------|------------|
| **HYCOM** (NRL) | Ocean IC/BC: temperature, salinity, currents, SSH | 1/12° (~9 km) |
| **Mercator/GLORYS** (Copernicus) | Alternative ocean IC/BC (often better in tropics) | 1/12° |

How to download HYCOM:
```bash
# HYCOM GLBy0.08 (1/12 degree) — free, no registration
# Via OPeNDAP / THREDDS
# Example using ncks (from NCO tools):
ncks -d time,0 -d lat,-15.0,15.0 -d lon,90.0,145.0 \
  "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ssh" \
  -o hycom_ssh.nc
```

How to download Mercator (Copernicus Marine):
```bash
pip install copernicusmarine
copernicusmarine subset \
  --dataset-id cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i \
  --variable uo --variable vo --variable thetao --variable so --variable zos \
  --minimum-longitude 90 --maximum-longitude 145 \
  --minimum-latitude -15 --maximum-latitude 15 \
  --start-datetime "2024-01-01" --end-datetime "2024-01-11"
```

### Wave BC (for SWAN)

SWAN can be forced at the open boundaries with spectral wave data from a global wave model:

| Source | What it provides |
|--------|-----------------|
| **GFS-Wave** (NCEP) | Global wave spectra, significant wave height, period, direction |
| **ECMWF WAM** | Alternative global wave model output |
| **WaveWatch III** (NOAA) | Can also be used for boundary spectra |

### Bathymetry and Terrain

For ROMS and SWAN (not confirmed for InaCAWO specifically, but standard practice):

| Data | Source | Resolution |
|------|--------|-----------|
| Bathymetry | **GEBCO 2023** or **SRTM15+** | 15 arc-second (~450 m) |
| Terrain (WRF) | Standard WRF geo data (USGS or MODIS) | 30-arc-second or finer |
| Coastline | GSHHS | High resolution |
| Tidal forcing | **TPXO9** or **FES2014** | For ROMS tidal boundaries |

## Running InaCAWO (COAWST Workflow)

The exact operational scripts are proprietary (Baron ROME), but the COAWST workflow follows a standard pattern. Here's what it looks like if you were to set up something similar.

### Step 1: Build COAWST

COAWST compiles WRF + ROMS + SWAN into one executable. You need all three model source codes plus the COAWST coupler layer.

```bash
# Clone COAWST
git clone https://github.com/USGS-CMG/COAWST.git
cd COAWST

# Edit coawst.bash to set:
# - WRF configuration (compiler, nesting, etc.)
# - ROMS application name and header file
# - SWAN coupling flags
# - MCT library path

# Build
./coawst.bash
```

The key configuration files:
- `coupling_esmf.dat` or `coupling.dat` — defines which models run and exchange intervals
- `namelist.input` — WRF configuration (physics, domains, timestep)
- `ocean_app.h` + `ocean.in` — ROMS configuration
- `INPUT` — SWAN configuration

### Step 2: Preprocess WRF (WPS)

Standard WRF preprocessing, same as standalone WRF:

```bash
# 1. Define domains
./geogrid.exe

# 2. Extract atmospheric data from GFS/ECMWF grib files
./ungrib.exe

# 3. Interpolate to WRF grid
./metgrid.exe

# 4. Generate initial and boundary condition files
./real.exe
# Produces: wrfinput_d01, wrfbdy_d01
```

### Step 3: Preprocess ROMS

ROMS needs its own grid, initial conditions, boundary conditions, and tidal forcing.

```bash
# 1. Create ROMS grid (using pyroms, gridgen, or MATLAB tools)
#    Define bathymetry, land mask, coordinate stretching
python create_roms_grid.py  # custom script

# 2. Interpolate HYCOM/Mercator onto ROMS grid for IC
python create_roms_ini.py   # interpolate T, S, u, v, zeta from HYCOM

# 3. Create boundary condition files
python create_roms_bry.py   # open boundary T, S, u, v, zeta from HYCOM

# 4. Create tidal forcing (if using tides)
python create_roms_tides.py  # from TPXO9 or FES2014
```

Common tools for ROMS preprocessing:
- `pyroms` — Python tools for ROMS grid generation and interpolation
- `ROMS IRG` — MATLAB-based ROMS tools
- `xESMF` — ESMF-based Python regridder (good for interpolation)
- `OTPS` — for extracting tidal constituents from TPXO

### Step 4: Preprocess SWAN

SWAN needs a grid (can share the ROMS grid or have its own), bathymetry, and boundary spectral conditions.

```bash
# SWAN input is controlled by an INPUT file
# Key settings:
# - Computational grid (can match ROMS grid)
# - Bathymetry input
# - Boundary spectra (from global wave model)
# - Physics: wind-wave generation, whitecapping, bottom friction, breaking
```

### Step 5: Configure Coupling

The coupling configuration tells COAWST which models run, on how many processors each, and how often they exchange data.

Example `coupling.dat`:
```
# Coupling parameters
NMODELS: 3
NGRIDS_ATM: 1
NGRIDS_OCN: 1
NGRIDS_WAV: 1
ATM_NPROCS: 512
OCN_NPROCS: 256
WAV_NPROCS: 128

# Synchronization interval (seconds)
ATM2OCN_COUPLING_INTERVAL: 600
OCN2ATM_COUPLING_INTERVAL: 600
ATM2WAV_COUPLING_INTERVAL: 600
WAV2ATM_COUPLING_INTERVAL: 600
OCN2WAV_COUPLING_INTERVAL: 600
WAV2OCN_COUPLING_INTERVAL: 600
```

The coupling interval is a trade-off — shorter intervals mean more accurate coupling but slower performance. Typical values range from 60 seconds to 1800 seconds depending on the application.

### Step 6: Run

```bash
# Run COAWST (single executable, all models)
mpirun -np 896 ./coawstM coupling.dat

# This launches WRF, ROMS, and SWAN simultaneously
# Each on its assigned processor set
# MCT handles the data exchange at synchronization points
```

For a domain this large at 3 km resolution, expect it to need **hundreds of cores** and many hours of wall time per forecast cycle.

### Step 7: Post-processing and Output

Output files from each model:
- **WRF**: `wrfout_d01_YYYY-MM-DD_HH:00:00` (NetCDF) — wind, rain, temperature, pressure, radiation
- **ROMS**: `ocean_his.nc`, `ocean_avg.nc` (NetCDF) — currents, SST, salinity, water level
- **SWAN**: spectral output files — wave height, period, direction, swell

```bash
# Example: extract significant wave height from SWAN output
python -c "
import xarray as xr
ds = xr.open_dataset('ocean_his.nc')
print(ds['Hwave'])  # SWAN wave height stored in ROMS output when coupled
"
```

In COAWST, SWAN output variables are often written into the ROMS history file for convenience — you get wave fields alongside ocean fields in the same NetCDF.

## Output Products (Operational InaCAWO)

BMKG disseminates InaCAWO forecasts through:

| Product | Variables | Platform |
|---------|----------|----------|
| Maritime forecast | Wind, waves, swell, currents | [peta-maritim.bmkg.go.id/ofs](https://peta-maritim.bmkg.go.id/ofs/) |
| MMS portal | Full model output visualization | [mms.bmkg.go.id](https://mms.bmkg.go.id/) |
| ROME web interface | Operations monitoring + output maps | Internal (Baron) |
| Baron Lynx | TV broadcast weather graphics | BMKG TV presentations |

Output variables at 3-hourly intervals for 10 days:

| Variable | Source Model | Use |
|----------|-------------|-----|
| Wind speed & direction | WRF | Shipping, sailing, fishing |
| Precipitation | WRF | Flash flood warnings |
| Significant wave height | SWAN | Maritime safety |
| Wave direction & period | SWAN | Shipping route planning |
| Swell | SWAN | Coastal hazard assessment |
| Ocean current speed & direction | ROMS | Navigation, SAR operations |
| Sea surface temperature | ROMS | Fisheries, coral bleaching |
| Salinity | ROMS | Estuarine forecasting |
| Water level / tides | ROMS | Port operations, flooding |

## HPC Infrastructure

The MMS-1 supercomputer at BMKG Jakarta:

| Spec | Value |
|------|-------|
| Total nodes | 372 (124 liquid-cooled blades × 3 nodes) |
| Peak performance | 1.31 Pflops |
| Storage | 2 PB low-latency |
| Vendor | Atos/Bull |
| Deployed | February 2023 |
| Operational | Mid-2023 |

MMS-2 is under development — targeting 2 Pflops with 6 PB storage.

## Before InaCAWO: Legacy Systems

InaCAWO replaced several standalone (uncoupled) systems:

| System | Model | Resolution | Status |
|--------|-------|-----------|--------|
| InaNWP | WRF (with data assimilation) | Various | Still operational alongside InaCAWO |
| Ina-Waves | WaveWatch III | 3 nested grids | Replaced by SWAN in InaCAWO |
| Ina-Flows | FVCOM | Unstructured | Replaced by ROMS in InaCAWO |

The upgrade from uncoupled to coupled is significant — in the old system, WRF didn't know about waves, FVCOM didn't know about the atmosphere, and WW3 used prescribed winds. Now they all interact dynamically.

## What's Not Publicly Known

The Baron whitepaper (June 2024) likely contains these details but couldn't be parsed:

- **WRF physics parameterizations** — we don't know which microphysics, cumulus, PBL, radiation, or land surface schemes InaCAWO uses. Given the 3 km resolution over a maritime domain, reasonable guesses based on standard COAWST practice would be Thompson/Morrison microphysics, no cumulus (3 km resolves convection), MYNN PBL, RRTMG radiation — but this is speculation.
- **WRF nesting strategy** — single 3 km domain or nested?
- **ROMS vertical levels** and coordinate stretching
- **ROMS tidal forcing** configuration (TPXO vs FES)
- **SWAN spectral settings** (frequency/direction bins)
- **Coupling synchronization interval** (600s? 300s? 900s?)
- **Data assimilation** — whether InaCAWO assimilates observations or runs cold from global IC
- **Exact model versions** (WRF v4.x, ROMS v3.x, SWAN v41.xx)

## Validation

From the AGU Ocean Sciences 2024 presentation by BMKG researchers:
- InaCAWO reports **>80% accuracy** for all output variables
- This is a general claim — detailed verification against buoys, altimetry, or reanalysis has not been published in peer-reviewed literature yet

For the older Ina-Waves (WW3), verification against Sentinel-3 altimetry showed RMSE of 0.5–1.3 m for significant wave height depending on the basin.

## Key References

1. Geonova et al. (2024), "The Development of a Novel High-Resolution Indonesia Coupled Atmospheric-Wave-Ocean (Ina-CAWO) Model" — AGU Ocean Sciences Meeting 2024
2. Baron Weather (2024), "The InaCAWO Metocean Operational Forecast System" — [whitepaper PDF](https://baronweather.com/hubfs/Brochures%20and%20eBooks/Whitepapers/Baron_BMKG-CAWO%20Whitepaper_06-2024.pdf)
3. Warner et al. (2010), "Development of a Coupled Ocean-Atmosphere-Wave-Sediment Transport (COAWST) Modeling System" — Ocean Modelling
4. CLS (2023), "Indonesian Maritime Meteorological System HPC operational at BMKG HQ" — [cls.fr](https://www.cls.fr/en/indonesian-maritime-meteorological-system-high-performance-computing-is-now-operational-at-bmkg-hq/)
5. Baron Weather, "The Ocean's Power Meets Baron Modeling Technology" — [baronweather.com](https://baronweather.com/weather-insights/baron-technology-ocean-modeling)
6. BMKG OFS Portal — [peta-maritim.bmkg.go.id/ofs](https://peta-maritim.bmkg.go.id/ofs/)
