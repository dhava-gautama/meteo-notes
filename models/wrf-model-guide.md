# WRF Model: Complete Guide from Data Download to Visualization

> A practical guide covering the full WRF (Weather Research and Forecasting) workflow.

---

## Table of Contents

1. [Workflow Overview](#1-workflow-overview)
2. [Downloading Initial & Boundary Condition Data](#2-downloading-initial--boundary-condition-data)
3. [WPS — WRF Preprocessing System](#3-wps--wrf-preprocessing-system)
4. [Running WRF](#4-running-wrf)
5. [Physics Options Reference](#5-physics-options-reference)
6. [Post-Processing & Visualization](#6-post-processing--visualization)
7. [System Requirements](#7-system-requirements)
8. [Troubleshooting & Tips](#8-troubleshooting--tips)

---

## 1. Workflow Overview

```
[Download GRIB data from RDA / NCEP]
        |
        v
[WPS: geogrid.exe]  -->  geo_em.d0*.nc           (static terrain / land use)
        |
[WPS: ungrib.exe]   -->  FILE:YYYY-MM-DD_HH      (intermediate met data)
        |
[WPS: metgrid.exe]  -->  met_em.d0*.YYYY-MM-DD_HH:00:00.nc
        |
        v
[WRF: real.exe]     -->  wrfinput_d0* + wrfbdy_d01
        |
        v
[WRF: wrf.exe]      -->  wrfout_d0*_YYYY-MM-DD_HH:MM:SS  (NetCDF output)
        |
        v
[Post-processing: wrf-python / xWRF / NCL / VAPOR / UPP]
        |
        v
[Visualization: matplotlib + cartopy / GrADS / VAPOR 3D]
```

---

## 2. Downloading Initial & Boundary Condition Data

### Primary Source: NCAR Research Data Archive

Register (free) at [gdex.ucar.edu](https://gdex.ucar.edu/) (formerly rda.ucar.edu).

### Common Datasets

| Dataset | ID | Resolution | Frequency | Period | Format |
|---|---|---|---|---|---|
| **NCEP GFS-FNL** | ds083.3 | 0.25 deg | 6-hourly | 2015–present | GRIB2 |
| **NCEP GFS-FNL** | ds083.2 | 1.0 deg | 6-hourly | 1999–present | GRIB |
| **NCEP GFS Forecast** | ds084.1 | 0.25 deg | 3-hr / 12-hr | 2015–present | GRIB2 |
| **ERA5 Reanalysis** | ds633.0 | ~31 km | 1-hourly | 1979–present | GRIB |
| **NCEP CFSR** | ds093.0 | 0.3–2.5 deg | 6-hourly | 1979–2011 | GRIB |
| **NCEP CFSv2** | ds094.0 | 0.2–2.5 deg | 6-hourly | 2011–present | GRIB |
| **NCEP NAM** | ds609.0 | 12 km | 6-hourly | 2012–present | GRIB |

### Which Dataset to Use?

- **Case studies / severe weather:** FNL (ds083.3) at 0.25 deg — incorporates more observational data than GFS (prepared ~1 hour later to allow extra obs ingestion)
- **Operational-style forecasting:** GFS (ds084.1) at 0.25 deg
- **Reanalysis / climate studies:** ERA5 (ds633.0) — offers 1-hourly temporal resolution
- **Real-time forecasts:** NOAA GFS from `ftpprd.ncep.noaa.gov`

### ERA5 Special Note

ERA5 requires **both** pressure-level AND surface-level data files. Download the complete global files — do not subset individual variables; `ungrib.exe` handles variable extraction.

### Download Methods

- **Web interface** — browse, select time range, download
- **wget/cURL scripts** — RDA generates scripts after you submit a data request (limit: 10 concurrent streams)
- **Globus** — for large transfers
- **Python API** — [github.com/NCAR/rda-apps-clients](https://github.com/NCAR/rda-apps-clients)

---

## 3. WPS — WRF Preprocessing System

Three programs that share `namelist.wps`, run in sequence: **geogrid** -> **ungrib** -> **metgrid**.

### 3.1 geogrid.exe

**Purpose:** Define simulation domains and interpolate static geographical data (terrain, land use, soil) onto the model grids.

**Output:** `geo_em.d01.nc`, `geo_em.d02.nc`, etc.

**Prerequisite:** Download static geographical data (~29 GB full resolution) from NCAR. Set `geog_data_path` in namelist.

**Key config (`&geogrid`):**
```
&geogrid
 parent_id         = 1, 1,
 parent_grid_ratio = 1, 3,
 i_parent_start    = 1, 31,
 j_parent_start    = 1, 17,
 e_we              = 100, 112,
 e_sn              = 100, 97,
 geog_data_res     = 'default', 'default',
 dx                = 15000,
 dy                = 15000,
 map_proj          = 'lambert',
 ref_lat           = 33.00,
 ref_lon           = -79.00,
 truelat1          = 30.0,
 truelat2          = 60.0,
 stand_lon         = -79.0,
 geog_data_path    = '/data/WPS_GEOG/',
/
```

### 3.2 ungrib.exe

**Purpose:** Extract meteorological fields from GRIB/GRIB2 files into WPS intermediate format.

**Output:** `FILE:YYYY-MM-DD_HH` (or custom prefix)

**Workflow:**
```bash
# 1. Link the Vtable for your data source
ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable

# 2. Link GRIB files
./link_grib.csh /path/to/grib/files/*

# 3. Run
./ungrib.exe
```

**Common Vtables:** `Vtable.GFS`, `Vtable.ERA-interim.pl`, `Vtable.ECMWF`, `Vtable.NAM`, `Vtable.CFSR` (40+ available)

**For ERA5** — run ungrib twice:
```bash
# Pressure levels
ln -sf Vtable.ERA-interim.pl Vtable
# set prefix = 'PRES' in namelist.wps
./ungrib.exe

# Surface data
ln -sf Vtable.ERA-interim.sfc Vtable
# set prefix = 'SFC' in namelist.wps
./ungrib.exe
```

Then in metgrid: `fg_name = './PRES', './SFC'`

### 3.3 metgrid.exe

**Purpose:** Horizontally interpolate intermediate meteorological data onto the simulation domains.

**Output:** `met_em.d01.YYYY-MM-DD_HH:00:00.nc` (one per domain per time step) — direct input for WRF's `real.exe`.

**Key config (`&metgrid`):**
```
&metgrid
 fg_name         = './FILE',
 io_form_metgrid = 2,
/
```

### Map Projections

| Projection | `map_proj` | Best For | Key Parameters |
|---|---|---|---|
| Lambert Conformal | `'lambert'` | Mid-latitudes | `truelat1`, `truelat2`, `stand_lon` |
| Polar Stereographic | `'polar'` | High latitudes | `truelat1`, `stand_lon` |
| Mercator | `'mercator'` | Low latitudes / E-W extent | `truelat1` |
| Lat-Lon | `'lat-lon'` | Global domains | `pole_lat`, `pole_lon` |

### Domain Nesting Rules

- `parent_grid_ratio` must be an **odd integer** (3:1 is standard, 5:1 is acceptable)
- Nested domain dimensions must satisfy: `(e_we - 1) % parent_grid_ratio == 0`
- `i_parent_start` / `j_parent_start` locate the nest's lower-left corner in the parent grid
- Nests must be at least **5 parent grid cells** from the parent domain boundary
- All domains must use the same number of vertical levels

### Complete namelist.wps Example

```
&share
 wrf_core         = 'ARW',
 max_dom          = 2,
 start_date       = '2024-01-15_00:00:00', '2024-01-15_00:00:00',
 end_date         = '2024-01-17_00:00:00', '2024-01-17_00:00:00',
 interval_seconds = 21600,
 io_form_geogrid  = 2,
/

&geogrid
 parent_id         = 1, 1,
 parent_grid_ratio = 1, 3,
 i_parent_start    = 1, 31,
 j_parent_start    = 1, 17,
 e_we              = 100, 112,
 e_sn              = 100, 97,
 geog_data_res     = 'default', 'default',
 dx                = 15000,
 dy                = 15000,
 map_proj          = 'lambert',
 ref_lat           = 33.00,
 ref_lon           = -79.00,
 truelat1          = 30.0,
 truelat2          = 60.0,
 stand_lon         = -79.0,
 geog_data_path    = '/data/WPS_GEOG/',
/

&ungrib
 out_format = 'WPS',
 prefix     = 'FILE',
/

&metgrid
 fg_name         = './FILE',
 io_form_metgrid = 2,
/
```

---

## 4. Running WRF

### real.exe vs ideal.exe

| | real.exe | ideal.exe |
|---|---|---|
| **Input** | `met_em.d0*.nc` from metgrid | `input_sounding` files |
| **Output** | `wrfinput_d0*` + `wrfbdy_d01` | `wrfinput_d01` only |
| **Use case** | Real weather events | Theoretical studies (supercell, LES) |
| **Parallelism** | Can run parallel | Must run serial |

### Execution

```bash
# 1. Run real.exe
mpirun -np 4 ./real.exe

# Verify success
tail rsl.out.0000
# Should show: "real_em: SUCCESS COMPLETE REAL_EM INIT"

# Check output files
ls -l wrfinput_d0* wrfbdy_d01

# 2. Run wrf.exe
mpirun -np 64 ./wrf.exe

# Verify success
tail rsl.out.0000
# Should show: "wrf: SUCCESS COMPLETE WRF"

# Output
ls wrfout_d0*
```

### namelist.input — Key Sections

**&time_control:**
```
&time_control
 run_days                 = 2,
 run_hours                = 0,
 start_year               = 2024, 2024,
 start_month              = 01,   01,
 start_day                = 15,   15,
 start_hour               = 00,   00,
 end_year                 = 2024, 2024,
 end_month                = 01,   01,
 end_day                  = 17,   17,
 end_hour                 = 00,   00,
 interval_seconds         = 21600,
 input_from_file          = .true., .true.,
 history_interval         = 60,    60,
 frames_per_outfile       = 1,     1,
 restart                  = .false.,
 restart_interval         = 1440,
 io_form_history          = 2,
 io_form_restart          = 2,
 io_form_input            = 2,
 io_form_boundary         = 2,
/
```

**&domains:**
```
&domains
 time_step                = 90,
 max_dom                  = 2,
 e_we                     = 100, 112,
 e_sn                     = 100, 97,
 e_vert                   = 45,  45,
 dzstretch_s              = 1.1,
 p_top_requested          = 5000,
 num_metgrid_levels       = 34,
 num_metgrid_soil_levels  = 4,
 dx                       = 15000, 5000,
 dy                       = 15000, 5000,
 grid_id                  = 1, 2,
 parent_id                = 0, 1,
 i_parent_start           = 1, 31,
 j_parent_start           = 1, 17,
 parent_grid_ratio        = 1, 3,
 parent_time_step_ratio   = 1, 3,
 feedback                 = 1,
/
```

### Time Step Rule of Thumb

```
time_step (seconds) = 6 x dx (km)
```

| Grid Spacing | Max time_step | Conservative |
|---|---|---|
| 36 km | 216 s | 120–144 s |
| 12 km | 72 s | 48–60 s |
| 4 km | 24 s | 16–20 s |
| 1 km | 6 s | 4–5 s |

### Adaptive Time Stepping (recommended)

```
 use_adaptive_time_step = .true.,
 step_to_output_time    = .true.,
 target_cfl             = 1.2,
 max_step_increase_pct  = 5, 51,
```

### Domain Nesting: 1-way vs 2-way

| | 2-way (`feedback = 1`) | 1-way (`feedback = 0` or ndown.exe) |
|---|---|---|
| **Communication** | Bidirectional — child feeds back to parent | Parent -> child only |
| **When to use** | Nest covers significant portion of parent | Large ratio jumps (> 5:1) |
| **Max ratio** | ~5:1 | Any (use intermediate nests) |
| **ndown.exe** | Not used | Required for separate runs |

---

## 5. Physics Options Reference

### Microphysics (`mp_physics`)

| Value | Scheme | Notes |
|---|---|---|
| 1 | Kessler | Simple warm-rain, no ice |
| 2 | Purdue Lin | 6-class with ice, snow, graupel |
| 4 | WSM5 | Mixed-phase |
| 6 | WSM6 | Ice, snow, graupel |
| 8 | **Thompson** | Popular, well-tested |
| 10 | **Morrison 2-moment** | Double-moment, widely used for research |
| 28 | **Thompson Aerosol-Aware** | Thompson with aerosol interactions |
| 50–53 | P3 variants | Predicted particle properties |

### Cumulus Parameterization (`cu_physics`)

| Value | Scheme | Notes |
|---|---|---|
| 0 | **None** | **Use when dx < ~4 km** (convection-permitting) |
| 1 | **Kain-Fritsch** | Most widely used; deep + shallow |
| 2 | Betts-Miller-Janjic | Adjustment scheme |
| 3 | **Grell-Freitas** | Scale-aware, works across resolutions |
| 6 | Tiedtke | Mass-flux scheme |
| 16 | New Tiedtke | Updated mass-flux |

> **Critical:** Turn OFF cumulus parameterization when dx < ~4 km. Use scale-aware schemes (Grell-Freitas) in the "grey zone" (4–10 km).

### Planetary Boundary Layer (`bl_pbl_physics`)

| Value | Scheme | Type | Paired Surface Layer (`sf_sfclay_physics`) |
|---|---|---|---|
| 1 | **YSU** | Non-local | 1 (Revised MM5) |
| 2 | **MYJ** | Local, TKE | 2 (Eta Similarity) |
| 5 | **MYNN2** | Local, TKE | 5 (MYNN) |
| 7 | ACM2 | Hybrid | 7 (Pleim-Xiu) or 1 |
| 8 | BouLac | Local, TKE | 1 or 2 |

> **Important:** PBL scheme MUST be paired with the correct surface layer scheme.

### Radiation

| Value | LW Scheme | SW Scheme | Notes |
|---|---|---|---|
| 1 | RRTM | Dudhia | Standard, fast |
| 4 | **RRTMG** | **RRTMG** | Most recommended |
| 5 | New Goddard | New Goddard | |
| 99 | GFDL | GFDL | Hurricane heritage |

Set `radt` (radiation interval, minutes) approximately equal to `dx` in km. Never exceed 30 min.

### Land Surface Model (`sf_surface_physics`)

| Value | Scheme | Soil Layers | Notes |
|---|---|---|---|
| 1 | 5-layer Thermal Diffusion | 5 | Simple, fast |
| 2 | **Noah** | 4 | Standard, widely used |
| 4 | **Noah-MP** | 4 | Most advanced, multi-physics |
| 5 | CLM4 | 10 | Community Land Model |

### Recommended Combinations

**General-purpose (research):**
```
mp_physics  = 8      (Thompson)
cu_physics  = 1      (KF for dx > 10 km) or 0 (dx < 4 km)
bl_pbl_physics = 1   (YSU)
sf_sfclay_physics = 1 (Revised MM5)
ra_lw_physics = 4    (RRTMG)
ra_sw_physics = 4    (RRTMG)
sf_surface_physics = 4 (Noah-MP)
```

**Convection-permitting (dx < 4 km):**
```
mp_physics  = 10     (Morrison 2-moment)
cu_physics  = 0      (OFF)
bl_pbl_physics = 5   (MYNN2)
sf_sfclay_physics = 5 (MYNN)
ra_lw_physics = 4    (RRTMG)
ra_sw_physics = 4    (RRTMG)
sf_surface_physics = 4 (Noah-MP)
```

---

## 6. Post-Processing & Visualization

### WRF Output Variables

Output files: `wrfout_d<domain>_<YYYY-MM-DD_HH:MM:SS>` in NetCDF format.

**Key native variables:**

| Variable | Description | Units |
|---|---|---|
| T | Perturbation potential temperature (total = T + 300 K) | K |
| P, PB | Perturbation and base pressure (total = P + PB) | Pa |
| PH, PHB | Perturbation and base geopotential (height = (PH+PHB)/9.81) | m^2/s^2 |
| U, V, W | Wind components (staggered grid, grid-relative) | m/s |
| QVAPOR | Water vapor mixing ratio | kg/kg |
| RAINC | Accumulated convective precipitation | mm |
| RAINNC | Accumulated grid-scale precipitation | mm |
| T2 | 2-meter temperature | K |
| Q2 | 2-meter mixing ratio | kg/kg |
| U10, V10 | 10-meter wind (grid-relative) | m/s |
| PSFC | Surface pressure | Pa |
| XLAT, XLONG | Latitude, Longitude | degrees |
| HGT | Terrain height | m |
| SWDOWN | Downward shortwave at surface | W/m^2 |
| OLR | Outgoing longwave at TOA | W/m^2 |
| SST | Sea surface temperature | K |

### Python with wrf-python (Recommended)

```python
from netCDF4 import Dataset
from wrf import getvar, interplevel, to_np, get_cartopy, latlon_coords
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature

# Open WRF output
ncfile = Dataset("wrfout_d01_2024-01-15_00:00:00")

# Extract diagnostic variables (wrf-python handles destaggering + unit conversion)
slp     = getvar(ncfile, "slp")         # Sea level pressure (hPa)
t2      = getvar(ncfile, "T2")          # 2-m temperature (K)
td2     = getvar(ncfile, "td2")         # 2-m dewpoint (C)
ua      = getvar(ncfile, "ua")          # U-wind earth-relative (m/s)
va      = getvar(ncfile, "va")          # V-wind earth-relative (m/s)
z       = getvar(ncfile, "z")           # Geopotential height (m)
p       = getvar(ncfile, "pressure")    # Pressure (hPa)
rh      = getvar(ncfile, "rh")          # Relative humidity (%)
cape_2d = getvar(ncfile, "cape_2d")     # CAPE and CIN
dbz     = getvar(ncfile, "dbz")         # Reflectivity (dBZ)
rain    = getvar(ncfile, "RAINC") + getvar(ncfile, "RAINNC")  # Total precip

# Interpolate to pressure level
ht_500 = interplevel(z, p, 500.)        # 500 hPa geopotential height

# Get projection and coordinates
cart_proj = get_cartopy(slp)
lats, lons = latlon_coords(slp)

# Plot
fig, ax = plt.subplots(subplot_kw={"projection": cart_proj})
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES, linewidth=0.3)

contours = ax.contourf(
    to_np(lons), to_np(lats), to_np(slp),
    levels=20, transform=crs.PlateCarree(), cmap="coolwarm"
)
plt.colorbar(contours, ax=ax, shrink=0.7, label="SLP (hPa)")
plt.title("Sea Level Pressure")
plt.savefig("slp_plot.png", dpi=150)
```

**Available wrf-python diagnostics:** `slp`, `tk`, `td`, `rh`, `z`, `pressure`, `ua`, `va`, `wspd_wdir`, `cape_2d`, `cape_3d`, `dbz`, `mdbz`, `avo`, `pvo`, `eth`, `omega`, `pw`, `helicity`, `updraft_helicity`, `cloudfrac`, and more.

### Python with xWRF (Modern xarray approach)

```python
import xarray as xr
import xwrf

ds = xr.open_dataset("wrfout_d01_2024-01-15_00:00:00").xwrf.postprocess()
# CF-compliant coordinates, destaggered variables,
# computed diagnostics (air_pressure, air_potential_temperature)
# Supports Dask for parallel processing of large datasets
```

### Other Tools

| Tool | Notes |
|---|---|
| **NCL** | Mature, publication-quality. Uses `wrf_user_getvar()`. In maintenance mode — Python recommended going forward. |
| **VAPOR** | Interactive 3D visualization from NCAR. Imports wrfout directly. Great for volumetric rendering, streamlines, isosurfaces. |
| **GrADS** | Via ARWpost converter. Mature 2D plotting, less common now. |
| **UPP** | NCEP's Unified Post Processor. Outputs GRIB2 for operational use. Computes derived fields on standard pressure levels. |

---

## 7. System Requirements

### Operating System

- **Linux** — primary supported platform (Ubuntu, CentOS, RHEL)
- macOS — supported, less common for production
- Windows — NOT natively supported (use WSL2)

### Compilers

| Compiler | Notes |
|---|---|
| **gfortran** (GNU) | Free, v9+ recommended |
| **ifort/ifx** (Intel) | Generally faster executables; Intel oneAPI is free |
| **gcc/g++** | Required for C components regardless of Fortran compiler |

### Required Libraries

| Library | Purpose |
|---|---|
| **NetCDF-C** (4.6.1+) | Primary I/O (mandatory) |
| **NetCDF-Fortran** (4.4+) | Fortran bindings (mandatory) |
| **HDF5** (1.8+) | Compression for NetCDF-4 |
| **zlib** | Compression for GRIB2 |
| **libpng** | Compression for GRIB2 |
| **JasPer** (1.900) | JPEG2000 for GRIB2 (WPS 4.4+ can build internally) |
| **MPI** (OpenMPI / MPICH / Intel MPI) | Parallel execution |

> **NetCDF-C and NetCDF-Fortran must be compiled with the same compiler as WRF.**

### Environment Variables

```bash
export NETCDF=/usr/local/netcdf
export HDF5=/usr/local/hdf5
export JASPERLIB=/usr/local/jasper/lib
export JASPERINC=/usr/local/jasper/include
export WRF_DIR=/path/to/WRF          # Required for WPS compilation
export PATH=$NETCDF/bin:$PATH
export LD_LIBRARY_PATH=$NETCDF/lib:$HDF5/lib:$LD_LIBRARY_PATH
```

### Hardware

| Configuration | CPU | RAM | Storage | Use Case |
|---|---|---|---|---|
| Learning | 4 cores | 8–16 GB | 50 GB | Small domain, coarse res |
| Desktop Research | 8–16 cores | 32–64 GB | 500 GB | Regional (15–30 km) |
| Serious Research | 64–256 cores | 128–512 GB | 2+ TB | High-res (1–3 km), nested |
| HPC / Operational | 256–1000+ cores | 512 GB+ | 10+ TB | Convection-permitting, climate |

### Build Process

```bash
# 1. Configure WRF
cd WRF
./configure
# Select compiler + parallelism (e.g., GNU dmpar)
# Select nesting (1=basic, 2=prescribed moves, 3=vortex-following)

# 2. Compile
./compile em_real >& compile.log &
# Verify: ls -l main/real.exe main/wrf.exe main/ndown.exe

# 3. Configure WPS
cd ../WPS
export WRF_DIR=../WRF
./configure              # or: ./configure --build-grib2-libs

# 4. Compile WPS
./compile >& compile.log &
# Verify: ls -l geogrid.exe ungrib.exe metgrid.exe
```

---

## 8. Troubleshooting & Tips

### CFL Errors

CFL violations = numerical instability. Fixes (in order):

1. **Reduce time step** — from `6*dx_km` down to `4*dx_km` or `3*dx_km`
2. **Enable adaptive time stepping** — `use_adaptive_time_step = .true.`
3. **Enable vertical velocity damping** — `w_damping = 1` in `&dynamics`
4. **Increase sound step damping** — `epssm = 0.2` (up to 0.5 for complex terrain)
5. **Smooth topography** — `smooth_cg_topo = .true.` in `&domains`
6. **Check input data** — inspect `met_em` files for unrealistic values

> Occasional CFL messages in `rsl.error` files that don't crash the model are fine — the model recovered. Only persistent crashes require intervention.

### Spin-up Time

| Simulation Type | Recommended Spin-up |
|---|---|
| Short-term forecast | 6–12 hours |
| General weather | 12–24 hours |
| Precipitation studies | 24–48 hours |
| Climate downscaling | ~1 month |

Start your simulation 12–24 hours before the period of interest, then discard spin-up from analysis. Grid nudging (FDDA) can reduce spin-up effects.

### SIZE MISMATCH Errors

Domain dimensions in `namelist.input` (e_we, e_sn, e_vert) don't match the input files. Adjust namelist to match `met_em` / `wrfinput` dimensions.

### Common Debugging Commands

```bash
# Check for CFL errors
grep -i cfl rsl.error.*

# Check completion
tail rsl.out.0000

# Inspect NetCDF headers
ncdump -h wrfout_d01_2024-01-15_00:00:00

# Check timestamps
ncdump -v Times wrfout_d01_2024-01-15_00:00:00

# Check namelist errors
grep -i error rsl.error.0000

# Verify dimensions match
ncdump -h met_em.d01.2024-01-15_00:00:00.nc | grep DIMENSION

# Clean and rebuild
./clean -a && ./configure && ./compile em_real >& compile.log &
```

### Memory Issues

- Run `ulimit -s unlimited` before execution (NOT for OpenMP/smpar builds)
- Increase MPI processes rather than memory per node for large domains
- Use `io_form_restart = 102` for restart files > 4 GB
- Ensure sufficient disk space — high-res output can be many GB per time step

---

## References

- [WRF Users Guide](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/)
- [WRF GitHub](https://github.com/wrf-model/WRF) | [WPS GitHub](https://github.com/wrf-model/WPS)
- [NCAR Geoscience Data Exchange (RDA)](https://gdex.ucar.edu/)
- [wrf-python Documentation](https://wrf-python.readthedocs.io/)
- [xWRF Documentation](https://xwrf.readthedocs.io/)
- [NCAR WRF Tutorial Presentations](https://www2.mmm.ucar.edu/wrf/users/tutorial/tutorial.html)
- [WRF Forum](https://forum.mmm.ucar.edu/)
