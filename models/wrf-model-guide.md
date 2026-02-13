# WRF Model: Complete Guide from Data Download to Visualization

> A practical guide covering the full WRF (Weather Research and Forecasting) workflow.
> Based on **WRF v4.7.1** (June 2025) — the latest stable release.

---

## Table of Contents

1. [Workflow Overview](#1-workflow-overview)
2. [Downloading Initial & Boundary Condition Data](#2-downloading-initial--boundary-condition-data)
3. [WPS — WRF Preprocessing System](#3-wps--wrf-preprocessing-system)
4. [Running WRF](#4-running-wrf)
5. [Physics Options Reference](#5-physics-options-reference)
6. [Outputting Hidden / Non-Default Variables](#6-outputting-hidden--non-default-variables)
7. [Post-Processing & Visualization](#7-post-processing--visualization)
8. [WRF Variants & Extensions](#8-wrf-variants--extensions)
9. [JMA NWP Settings & Recommendations](#9-jma-nwp-settings--recommendations)
10. [System Requirements](#10-system-requirements)
11. [Troubleshooting & Tips](#11-troubleshooting--tips)
12. [What's New in WRF v4.7](#12-whats-new-in-wrf-v47)

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

- **Case studies / severe weather:** FNL (ds083.3) at 0.25 deg — incorporates more observational data than GFS
- **Operational-style forecasting:** GFS (ds084.1) at 0.25 deg
- **Reanalysis / climate studies:** ERA5 (ds633.0) — offers 1-hourly temporal resolution
- **Real-time forecasts:** NOAA GFS from `ftpprd.ncep.noaa.gov`

### ERA5 Special Note

ERA5 requires **both** pressure-level AND surface-level data files. Download complete global files — `ungrib.exe` handles variable extraction.

### Download Methods

- **Web interface** — browse, select time range, download
- **wget/cURL scripts** — RDA generates scripts after data request (limit: 10 concurrent streams)
- **Globus** — for large transfers
- **Python API** — [github.com/NCAR/rda-apps-clients](https://github.com/NCAR/rda-apps-clients)

---

## 3. WPS — WRF Preprocessing System

Three programs sharing `namelist.wps`, run in sequence: **geogrid** -> **ungrib** -> **metgrid**.

### 3.1 geogrid.exe

**Purpose:** Define simulation domains and interpolate static geographical data (terrain, land use, soil) onto model grids.

**Output:** `geo_em.d01.nc`, `geo_em.d02.nc`, etc.

**Prerequisite:** Download static geographical data (~29 GB full resolution) from NCAR.

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

**Purpose:** Horizontally interpolate intermediate meteorological data onto simulation domains.

**Output:** `met_em.d01.YYYY-MM-DD_HH:00:00.nc` — direct input for `real.exe`.

### Map Projections

| Projection | `map_proj` | Best For | Key Parameters |
|---|---|---|---|
| Lambert Conformal | `'lambert'` | Mid-latitudes | `truelat1`, `truelat2`, `stand_lon` |
| Polar Stereographic | `'polar'` | High latitudes | `truelat1`, `stand_lon` |
| Mercator | `'mercator'` | Low latitudes / E-W extent | `truelat1` |
| Lat-Lon | `'lat-lon'` | Global domains | `pole_lat`, `pole_lon` |

### Domain Nesting Rules

- `parent_grid_ratio` must be an **odd integer** (3:1 standard, 5:1 acceptable)
- Nested dimensions: `(e_we - 1) % parent_grid_ratio == 0`
- Nests must be at least **5 parent grid cells** from parent boundary
- All domains share the same number of vertical levels

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
tail rsl.out.0000   # Should show: "SUCCESS COMPLETE REAL_EM INIT"

# 2. Run wrf.exe
mpirun -np 64 ./wrf.exe
tail rsl.out.0000   # Should show: "SUCCESS COMPLETE WRF"
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

### Dynamics & Damping Options (`&dynamics`)

The `&dynamics` section controls diffusion, damping, and numerical stability. These are critical for preventing model blow-ups, especially over complex terrain or in high-resolution runs.

**Diffusion:**

| Parameter | Values | Description |
|---|---|---|
| `diff_opt` | **1** = simple diffusion (gradient on model levels) | How horizontal diffusion is computed |
| | **2** = full diffusion (evaluated on constant-z surfaces; recommended for real cases) | |
| `km_opt` | **1** = constant coefficient | How diffusion coefficient K is computed |
| | **2** = 1.5-order TKE closure (recommended for LES/dx < 1 km) | |
| | **3** = Smagorinsky first-order (requires `diff_opt=2`) | |
| | **4** = horizontal Smagorinsky, vertical from PBL (recommended for most real cases) | |
| `diff_6th_opt` | **0** = off | 6th-order numerical diffusion — removes 2Δx noise |
| | **1** = on (applied everywhere) | |
| | **2** = on (not applied in vertical; recommended) | |
| `diff_6th_factor` | 0.12 (default) | Strength of 6th-order diffusion. Increase to 0.2–0.25 for noisy runs. |

**Damping (upper boundary):**

| Parameter | Values | Description |
|---|---|---|
| `damp_opt` | **0** = no upper damping | Prevents wave reflection from model top |
| | **1** = increased diffusion in upper layers | |
| | **2** = Rayleigh relaxation (recommended for dx > ~5 km) | |
| | **3** = w-Rayleigh (damps vertical velocity only; recommended for convection-permitting) | |
| `dampcoef` | 0.2 (typical) | Damping coefficient. For `damp_opt=2,3`: inverse e-folding time (1/s); 0.01–0.2. |
| `zdamp` | 5000 (m) | Depth of damping layer below model top |

**Stability controls:**

| Parameter | Default | Description |
|---|---|---|
| `w_damping` | 0 | Vertical velocity damping. Set to **1** to stabilize runs with strong convection or steep terrain. |
| `epssm` | 0.1 | Off-centering for vertical sound waves. Increase to **0.2–0.5** for complex terrain to reduce CFL errors. |
| `smdiv` | 0.1 | Divergence damping. Increase up to 0.2 for stability. |
| `emdiv` | 0.01 | External-mode divergence damping. |
| `time_step_sound` | 0 (auto) | Number of sound steps per time step. Auto-calculated from CFL. Override with 4–6 if unstable. |

**Typical `&dynamics` for real cases:**
```
&dynamics
 hybrid_opt           = 2,          ! Hybrid sigma-pressure vertical coordinate (v4.0+)
 diff_opt             = 2, 2,
 km_opt               = 4, 4,
 diff_6th_opt         = 2, 2,
 diff_6th_factor      = 0.12, 0.12,
 damp_opt             = 3,
 zdamp                = 5000., 5000.,
 dampcoef             = 0.2, 0.2,
 w_damping            = 1,
 epssm                = 0.1, 0.1,
/
```

> **`hybrid_opt = 2`** (available since v4.0): Uses hybrid sigma-pressure vertical coordinate instead of pure terrain-following sigma. Significantly reduces pressure-gradient errors over steep terrain. **Recommended for all real-data runs.**

### Vertical Level Configuration

WRF uses **terrain-following eta (η) levels** where η=1 at the surface and η=0 at the model top. The distribution of levels controls how well the model resolves the boundary layer, jet streams, and the tropopause.

**Default behavior:** If you only set `e_vert` (number of levels), WRF auto-generates levels using `dzstretch_s` (surface stretching) and `dzstretch_d` (deep stretching).

| Parameter | Default | Description |
|---|---|---|
| `e_vert` | 45 | Total number of vertical levels (all domains must match) |
| `p_top_requested` | 5000 (Pa) | Model top pressure. 5000 Pa ≈ 50 km. Use 1000–2000 Pa for stratospheric studies. |
| `dzstretch_s` | 1.1 | Stretching factor near surface. Lower = more levels in PBL. |
| `dzstretch_d` | 1.1 | Stretching factor aloft. |
| `auto_levels_opt` | 2 | Level generation algorithm. 2 = improved (v4.1+). |

**Custom eta levels — manual specification:**

For full control, specify all eta values explicitly in `namelist.input`:
```
&domains
 e_vert           = 51,
 eta_levels       = 1.000, 0.9975, 0.995, 0.990, 0.985,
                    0.980, 0.970, 0.960, 0.945, 0.930,
                    0.910, 0.890, 0.865, 0.840, 0.810,
                    0.780, 0.750, 0.715, 0.680, 0.640,
                    0.600, 0.560, 0.520, 0.480, 0.440,
                    0.400, 0.365, 0.330, 0.300, 0.270,
                    0.240, 0.215, 0.190, 0.170, 0.150,
                    0.130, 0.115, 0.100, 0.085, 0.070,
                    0.060, 0.050, 0.040, 0.030, 0.023,
                    0.016, 0.010, 0.006, 0.003, 0.001,
                    0.000,
/
```

**Guidelines for vertical level design:**

| Purpose | Recommendation |
|---|---|
| PBL resolution | At least **8–10 levels** below 1 km AGL. First level at ~10–20 m. |
| Convection-permitting | 50–80 levels. Dense spacing in PBL and near tropopause. |
| JMA-style (MSM) | 96 levels with ~20 m first level and 37.5 km model top. |
| Climate downscaling | 40–50 levels sufficient. |
| Urban / LES | 80+ levels. First level at ~5 m. Very fine spacing in lowest 500 m. |

> **Common mistake:** Too few levels in the PBL (lowest 1–2 km). This causes poor representation of the nocturnal boundary layer, low-level jets, and surface fluxes. If your 2-m temperature or 10-m wind verification is poor, check your lowest-level spacing first.

### Lateral Boundary Conditions

Lateral boundary conditions (LBCs) feed external data into the outermost domain. The outermost domain is the only one that needs LBCs — nested domains get their boundaries from the parent.

**Types (`specified_bdy`):**

| Type | Setting | Description |
|---|---|---|
| **Specified** | `specified = .true.` (default for real cases) | Values at boundary come from external data (GFS, ERA5, etc.). Relaxation zone blends external data with model solution. Standard for all real-weather runs. |
| **Periodic** | `periodic_x = .true.` and/or `periodic_y = .true.` | Wrap-around boundaries. For idealized simulations (e.g., channel flow, squall lines). |
| **Open** | `open_xs = .true.`, etc. | Radiation/open boundary — waves pass through. For idealized cases. |
| **Symmetric** | `symmetric_xs = .true.`, etc. | Mirror boundary. For idealized cases. |

**Relaxation zone (`spec_bdy_width`):**

```
&bdy_control
 spec_bdy_width  = 5,        ! Number of grid points in relaxation zone (default 5)
 spec_zone       = 1,        ! Points with full specification (default 1)
 relax_zone      = 4,        ! Points with relaxation/nudging blend (default 4)
 specified       = .true.,  .false.,
 nested          = .false., .true.,
/
```

- `spec_bdy_width = spec_zone + relax_zone`
- Larger relaxation zone (8–10 points) helps prevent boundary artifacts for long simulations
- Nested domains always use `nested = .true.` — they get their LBCs from the parent domain

**Boundary update frequency (`interval_seconds`):**

| Source Data | Typical `interval_seconds` | Notes |
|---|---|---|
| GFS FNL (6-hourly) | 21600 | Standard. WRF linearly interpolates between updates. |
| ERA5 (1-hourly) | 3600 | Better for rapid weather changes. |
| GFS forecast (3-hourly) | 10800 | Good temporal resolution. |

> **Important:** Coarse temporal LBC updates (e.g., 6-hourly) can cause boundary artifacts — sudden jumps in wind/temperature every 6 hours that propagate inward. If you see periodic noise near boundaries, increase `interval_seconds` by using higher-frequency input data, or increase `spec_bdy_width`.

### Restart & Cycling Runs

**Restart runs** let you break long simulations into segments, resume after crashes, or branch experiments from a checkpoint.

**Namelist configuration:**
```
&time_control
 ! --- Initial run ---
 restart              = .false.,
 restart_interval     = 720,          ! Write restart file every 720 min (12 hr)
 io_form_restart      = 2,            ! NetCDF format

 ! --- Restart run (change these) ---
 ! restart             = .true.,
 ! start_year/month/day/hour = restart time
 ! (keep end time as original)
/
```

**Restart procedure:**
```bash
# 1. Locate restart file
ls wrfrst_d01_2024-01-16_00:00:00

# 2. Edit namelist.input
#    - Set restart = .true.
#    - Set start time to match restart file timestamp
#    - Keep end time unchanged

# 3. Run wrf.exe (do NOT re-run real.exe)
mpirun -np 64 ./wrf.exe
```

**Key rules:**
- Restart files must match the exact domain configuration (grid size, physics, vertical levels)
- For restart files > 4 GB: use `io_form_restart = 102` (parallel NetCDF)
- Physics options can generally be changed between restart segments, but microphysics changes may cause instabilities
- If you change physics, consider a short spin-up period in the new segment

**Cycling (continuous forecasting):**

For operational-style cycling (e.g., 6-hourly analysis-forecast cycles):
```
Cycle 1: real.exe -> wrf.exe (00Z–06Z) -> wrfout at 06Z
              |
Cycle 2: wrfout(06Z) -> WRF-DA -> updated wrfinput -> wrf.exe (06Z–12Z)
              |
Cycle 3: wrfout(12Z) -> WRF-DA -> updated wrfinput -> wrf.exe (12Z–18Z)
              ...
```

Each cycle uses the previous forecast as the first guess for data assimilation, producing improved initial conditions.

### SST Update for Multi-Day Runs

By default, WRF uses a static SST from the initial condition throughout the entire simulation. For runs longer than ~48 hours over ocean, time-varying SST significantly improves surface fluxes and maritime weather.

**Enable SST update:**
```
&physics
 sst_update = 1,
/

&time_control
 auxinput4_inname     = "wrflowinp_d<domain>"
 auxinput4_interval   = 360, 360,    ! SST update interval (minutes)
 io_form_auxinput4    = 2,
/
```

**Workflow:**
```bash
# 1. Include SST data in WPS processing
#    (GFS already contains SST; ERA5 requires surface-level files)

# 2. real.exe will automatically create wrflowinp_d01 if it detects
#    time-varying SST fields in met_em files

# 3. Verify:
ncdump -h wrflowinp_d01 | grep SST
```

**Fields updated from `wrflowinp`:**

| Field | Description |
|---|---|
| SST | Sea surface temperature |
| VEGFRA | Vegetation fraction (seasonal) |
| ALBEDO | Surface albedo (seasonal) |
| LAI | Leaf area index |

> **Tip:** For tropical cyclone simulations, SST update is essential — storm-induced cooling from ocean mixing can reduce TC intensity by 20–50%. Consider coupling with an ocean model (WRF-Hydro or POM) for full feedback.

---

## 5. Physics Options Reference

### Microphysics (`mp_physics`)

| Value | Scheme | Description |
|---|---|---|
| 1 | Kessler | Simple warm-rain only (cloud, rain). No ice-phase processes. Educational / idealized cases only. |
| 2 | Purdue Lin | Single-moment, 6-class (vapor, cloud, rain, ice, snow, graupel). Computationally cheap. Suitable for coarse-resolution runs where detailed microphysics is not critical. |
| 6 | WSM6 | Single-moment, 6-class. Improved ice nucleation over Lin. Good balance of speed and accuracy for general simulations. Developed by NCEP/KMA (Hong & Lim 2006). |
| 8 | **Thompson** | Hybrid single/double-moment — double-moment for rain and ice, single for others. Excellent for mixed-phase clouds, winter precipitation, and snowfall. One of the most widely validated and popular schemes. |
| 10 | **Morrison 2-moment** | Full double-moment for cloud, rain, ice, snow, graupel (predicts both mass and number concentration). Best for research requiring accurate size distributions. Standard choice for convection-permitting simulations. |
| 18 | **NSSL 2-moment** | Full double-moment with optional hail category and CCN prediction. Strong for severe convective storms — explicitly predicts graupel/hail size. Consolidated in v4.6+ (controls via `nssl_*` flags). |
| 27 | **UFS Double Moment** | 7-class, double-moment with aerosol-aware CCN activation. Developed for NOAA's Unified Forecast System. New in v4.7 (from S.-Y. Hong). |
| 28 | **Thompson Aerosol-Aware** | Extension of Thompson with prognostic aerosol (water- and ice-friendly). Aerosol-cloud interactions affect droplet number and precipitation efficiency. Good for air quality-sensitive forecasts. |
| 29 | **RCON** | Improved representation of warm rain and drizzle processes. Better collision-coalescence and autoconversion. New in v4.7 (Conrick et al. 2023). |
| 50–53 | P3 variants | Predicted Particle Properties — represents ice with a single free category using prognostic properties (mass, number, rime fraction, rime density) instead of fixed categories. Eliminates artificial ice-type conversions. |

> **NSSL consolidation (v4.6+):** `mp_physics=17, 19, 21, 22` are deprecated. Use `mp_physics=18` with `nssl_2moment_on`, `nssl_ccn_on`, `nssl_hail_on` flags instead.

### Cumulus Parameterization (`cu_physics`)

| Value | Scheme | Description |
|---|---|---|
| 0 | **None** | No cumulus parameterization. **Required when dx < ~4 km** — at convection-permitting resolution, convection is explicitly resolved by the model grid. |
| 1 | **Kain-Fritsch** | Mass-flux scheme with CAPE-removal closure. Handles both deep and shallow convection. The most widely used and validated cumulus scheme. Triggers when CAPE exceeds a threshold and removes CAPE over a specified timescale. Best for dx > 10 km. |
| 2 | Betts-Miller-Janjic | Relaxation/adjustment scheme — adjusts thermodynamic profiles toward observed post-convective reference profiles. Tends to produce lighter, more widespread precipitation. Commonly used for tropical and hurricane simulations. |
| 3 | **Grell-Freitas** | Scale-aware ensemble mass-flux scheme. Automatically reduces its contribution as resolution increases toward convection-permitting scales. Ideal for domains spanning the "grey zone" (4–10 km) and multi-resolution nesting without manually toggling cumulus on/off. |
| 6 | Tiedtke | Mass-flux scheme with moisture-convergence closure. Originally from ECMWF. Handles deep, shallow, and mid-level convection separately. Commonly used in global and tropical applications. |
| 16 | New Tiedtke | Updated version of Tiedtke with improved detrainment and CAPE-based closure for deep convection. Used operationally in several global models. Better representation of organized convection and diurnal cycle. |

> **Critical:** Turn OFF cumulus parameterization when dx < ~4 km. Use scale-aware schemes (Grell-Freitas) in the "grey zone" (4–10 km).

### Shallow Cumulus (`shcu_physics`)

Separate parameterization for shallow (non-precipitating) convection. Most PBL and deep cumulus schemes already include some shallow convection treatment, so this is only needed in specific configurations.

| Value | Scheme | Description |
|---|---|---|
| 0 | **None** (default) | Shallow convection handled internally by PBL scheme (MYNN-EDMF) or cumulus scheme (KF). Sufficient for most applications. |
| -1 | Use `cu_physics` | Relies on the deep cumulus scheme to also handle shallow convection. Works with KF and Grell-Freitas. |
| 2 | UW (Park & Bretherton) | University of Washington scheme. Mass-flux approach based on convective inhibition and TKE. Good for marine stratocumulus transition to trade cumulus. Used in CESM/CAM. |
| 3 | GRIMS (Global/Regional) | Korean GRIMS shallow convection. Mass-flux with moisture-convergence closure. |
| 5 | **Deng (MYNN-based)** | Designed specifically for WRF-Solar. Uses MYNN TKE to drive shallow cumulus mass flux. Requires MYNN (`bl_pbl_physics=5`) or MYJ (`bl_pbl_physics=2`) PBL. Improves cloud-radiation interaction for solar irradiance forecasting. |

> **When to use:** Only activate `shcu_physics` if your deep cumulus scheme does not handle shallow convection (e.g., BMJ) or if you need explicit control for solar/marine applications. Using `shcu_physics=5` with MYNN PBL is recommended for WRF-Solar.

### Planetary Boundary Layer (`bl_pbl_physics`)

| Value | Scheme | Type | Paired SfcLay | Description |
|---|---|---|---|---|
| 1 | **YSU** | Non-local | 1 (Revised MM5) | Yonsei University scheme. Uses counter-gradient flux to represent large eddies that transport heat/moisture across the entire PBL depth. Well-mixed boundary layers, good for clear-sky daytime convective conditions. Most popular general-purpose PBL scheme. |
| 2 | **MYJ** | Local, TKE | 2 (Eta Similarity) | Mellor-Yamada-Janjic Level 2.5. Prognostic TKE with local diffusion only — mixing depends on local gradients. Tends to produce shallower, less well-mixed PBLs than YSU. Good for stable/nocturnal conditions. Heritage from NAM/Eta model. |
| 5 | **MYNN-EDMF** | Local, TKE | 5 (MYNN) | Mellor-Yamada-Nakanishi-Niino Level 2.5 with Eddy-Diffusivity Mass-Flux extension for shallow convection. Prognostic TKE and sub-grid clouds. Superior for marine/coastal environments, fog, and low cloud prediction. Increasingly the recommended default. |
| 7 | ACM2 | Hybrid | 7 (Pleim-Xiu) or 1 | Asymmetric Convective Model v2. Combines non-local upward transport (convective plumes) with local downward diffusion. Designed for air quality applications — pairs with Pleim-Xiu land surface for EPA CMAQ modeling chain. |
| 8 | BouLac | Local, TKE | 1 or 2 | Bougeault-Lacarrère scheme. TKE-based local closure with diagnostic mixing length from stability and terrain geometry. Designed for complex mountainous terrain where terrain-influenced turbulence is important. |

> **v4.7 note:** MYNN-EDMF is now a git submodule, refactored to k-only scheme (10–15% faster). Module names changed from `*_mynn_*` to `*_mynnedmf_*`.

> **Important:** PBL scheme MUST be paired with the correct surface layer scheme.

### Surface Layer (`sf_sfclay_physics`)

The surface layer scheme computes friction velocity, exchange coefficients, and surface fluxes (heat, moisture, momentum) that are passed to the PBL scheme. Each PBL scheme requires a specific compatible surface layer.

| Value | Scheme | Compatible PBL | Description |
|---|---|---|---|
| 1 | **Revised MM5** | YSU (1), BouLac (8) | Monin-Obukhov similarity theory with Carlson-Boland viscous sub-layer and standard stability functions (Paulson 1970, Dyer-Hicks, Webb). Computes bulk transfer coefficients for heat, moisture, and momentum over land and water. The most widely used surface layer scheme. |
| 2 | **Eta Similarity** | MYJ (2) | Janjic's Monin-Obukhov implementation from the Eta/NAM model. Iterative solution for surface fluxes with Zilitinkevich thermal roughness. More conservative (lower) surface fluxes than MM5 in strong instability. Required for MYJ PBL. |
| 5 | **MYNN** | MYNN-EDMF (5) | Nakanishi-Niino surface layer. Similar to Eta similarity but with improved roughness length formulation over water (Charnock + viscous) and better handling of very stable conditions. Includes sub-grid cloud fraction coupling. Required for MYNN PBL. |
| 7 | **Pleim-Xiu** | ACM2 (7) | Designed for EPA air quality modeling. Uses Pleim (2006) stability functions optimized for the ACM2 PBL. Indirect soil moisture adjustment from 2-m observations. Required for ACM2 PBL. |
| 91 | **Old MM5** | Any (legacy) | Original MM5 surface layer without viscous sub-layer. Retained for backwards compatibility only. Not recommended for new simulations. |

> **Rule:** Always check the PBL–surface layer pairing. Using mismatched schemes (e.g., MYNN PBL with MM5 surface layer) will produce physically inconsistent results — the model may run but surface fluxes will be wrong.

### Radiation

| Value | LW Scheme | SW Scheme | Description |
|---|---|---|---|
| 1 | RRTM | Dudhia | **LW:** Rapid Radiative Transfer Model — correlated-k method with 16 spectral bands. Accurate longwave fluxes and cooling rates. **SW:** Dudhia — simple broadband scheme, fast but less accurate for clear-sky direct beam. Good default for quick simulations where radiation speed matters. |
| 4 | **RRTMG** | **RRTMG** | RRTM for GCMs — updated RRTM with Monte Carlo Independent Column Approximation (McICA) for sub-grid cloud variability. 14 SW bands and 16 LW bands. Accounts for cloud overlap effects. The most recommended radiation scheme — balances accuracy and speed. Supports aerosol direct effects. |
| 5 | New Goddard | New Goddard | NASA Goddard scheme with 11 SW bands and 10 LW bands. Good for tropical and cloud-radiation interaction studies. Includes ozone and CO2 effects. Slightly more expensive than RRTMG. |
| 99 | GFDL | GFDL | NOAA Geophysical Fluid Dynamics Laboratory scheme. Heritage from HWRF hurricane model. Optimized for tropical cyclone simulations — well-tested for warm-core vortex radiation balance. Simmons broadband approach. |

Set `radt` (radiation interval, minutes) approximately equal to `dx` in km. Never exceed 30 min.

### Land Surface Model (`sf_surface_physics`)

| Value | Scheme | Soil Layers | Description |
|---|---|---|---|
| 1 | 5-layer Thermal Diffusion | 5 | Simplest option — soil temperature diffusion only. No vegetation canopy, no explicit soil moisture evolution, no evapotranspiration. Fast but unrealistic surface fluxes. Only for quick tests. |
| 2 | **Noah** | 4 | NCEP operational standard. Predicts soil temperature, soil moisture, skin temperature, snowpack, and canopy water. Includes frozen soil physics, evapotranspiration via Penman scheme, and 24-category USGS or 20-category MODIS land use. Reliable and well-tested for general use. |
| 4 | **Noah-MP** | 4 | Noah Multi-Physics — most advanced option. Adds separate vegetation canopy, dynamic vegetation, multi-layer snowpack, groundwater table, and multiple options for each process (e.g., Ball-Berry stomatal resistance, runoff schemes). Allows toggling sub-models on/off. Best for research requiring realistic land-atmosphere coupling. Now a git submodule in v4.7. |
| 5 | CLM4 | 10 | Community Land Model v4 (from NCAR CESM). 10 soil layers to 3.4 m depth, up to 5 snow layers, sophisticated biogeophysics (photosynthesis-conductance, urban canopy). Most detailed treatment of sub-grid land heterogeneity (PFTs). Computationally expensive. Best for climate and long-duration simulations. |
| 7 | Pleim-Xiu | 2 | Two-layer soil model with indirect soil moisture nudging from observed 2-m temperature/humidity. Designed for EPA air quality modeling chain (pairs with ACM2 PBL). Supports 61-category MODIS LCZ in v4.7 for urban studies. |

### Urban Canopy Model (`sf_urban_physics`)

Urban physics parameterizations represent the effects of buildings, streets, and urban materials on surface energy balance, wind flow, and temperature. Essential for urban heat island (UHI) studies and city-scale weather forecasting.

| Value | Scheme | Description |
|---|---|---|
| 0 | **None** (default) | Urban grid cells treated as rough, impervious land surface by the LSM. Adequate if urban areas are not the focus. |
| 1 | **Single-layer UCM** | Single-layer Urban Canopy Model. Represents a simplified urban canyon (roof, wall, road) with bulk energy balance. Accounts for shadowing, reflection, and anthropogenic heat. Moderate computational cost. Good starting point for UHI studies. |
| 2 | **BEP** (Multi-layer) | Building Energy Parameterization. Multi-layer scheme where buildings interact with multiple atmospheric levels — not just the surface. Explicitly resolves drag, turbulence, and radiation trapping at each building level. Requires `bl_pbl_physics = 2 or 8` (local TKE schemes only). |
| 3 | **BEP+BEM** | BEP with Building Energy Model. Adds indoor energy budget — HVAC systems, waste heat, indoor temperature. Captures anthropogenic heat release from air conditioning. Most realistic for megacity simulations and energy demand studies. Requires `bl_pbl_physics = 2 or 8`. |

**Key settings:**
```
&physics
 sf_urban_physics   = 1, 1,          ! UCM option
 use_wudapt_lcz     = 1,             ! Use WUDAPT Local Climate Zones (v4.3+)
 num_urban_layers   = 1040,          ! Number of urban layers (BEP/BEM only)
/
```

**Land use requirements:**
- UCM requires urban land-use categories in the WPS static data
- **WUDAPT LCZ** (Local Climate Zones) provides 10 urban sub-categories (compact high-rise, open low-rise, etc.) instead of a single "urban" class — greatly improves UHI representation
- v4.7: Pleim-Xiu LSM now supports 61-category MODIS with LCZ

> **BEP/BEM constraint:** Only works with local TKE-based PBL schemes (MYJ=2 or BouLac=8) because multi-layer urban drag directly modifies the vertical diffusion profile. Non-local schemes (YSU) are incompatible.

### Ocean Physics (`sf_ocean_physics`)

Simple ocean mixed-layer models that provide SST feedback without requiring a full ocean model coupling. Useful for tropical cyclone simulations where storm-induced ocean cooling affects intensity.

| Value | Scheme | Description |
|---|---|---|
| 0 | **None** (default) | SST is prescribed (static or from `sst_update`). No ocean response to atmospheric forcing. Fine for most applications over land or short runs. |
| 1 | **Simple Mixed-Layer** | 1D mixed-layer ocean model. Prognostic SST responds to surface heat/momentum fluxes. Fixed mixed-layer depth (set by `ocean_z_levels` or from initialization). Captures basic storm-induced SST cooling. Low cost. |
| 2 | **3D Price-Weller-Pinkel** | 3D ocean mixed-layer model with horizontal advection. Better representation of ocean current response and upwelling under tropical cyclones. More realistic SST cooling pattern (cold wake). Higher cost but still much cheaper than full ocean coupling. |

**Namelist example:**
```
&physics
 sf_ocean_physics     = 1,
 oml_hml0             = 50,        ! Initial mixed-layer depth (m)
 oml_gamma            = 0.14,      ! Temperature lapse rate below mixed layer (K/m)
/
```

> **When to use:** Primarily for tropical cyclone studies lasting > 48 hours, where SST cooling by 1–3°C under the storm eye significantly reduces intensity. For runs over land-dominated domains, keep `sf_ocean_physics = 0`.

### Gravity Wave Drag (`gwd_opt`)

Parameterizes the drag on the large-scale flow caused by sub-grid terrain-generated gravity waves that propagate vertically and break aloft. Important at coarser resolutions where individual mountain features are not resolved.

| Value | Scheme | Description |
|---|---|---|
| 0 | **None** (default) | No gravity wave drag. Acceptable for high-resolution (dx < ~10 km) where terrain is well resolved, or for flat domains. |
| 1 | **GWD** | Gravity wave drag + low-level blocking. Based on Kim & Arakawa (1995). Represents both mountain wave breaking aloft (decelerates upper-level flow) and flow blocking around sub-grid orography (decelerates low-level flow). |
| 3 | **GWD + Small-scale GWD** | Adds sub-grid scale orographic drag from small terrain features not captured even by the resolved terrain. Developed for the Unified Forecast System (UFS). New in v4.6. |

**When to activate:**
- dx > ~20 km: **Recommended** — sub-grid mountains have significant unresolved drag
- dx 10–20 km: Consider enabling if domain has significant topography
- dx < 10 km: Generally not needed — terrain is well-resolved

```
&physics
 gwd_opt = 1,
/
```

### Land Use / Land Cover Datasets

WRF uses land use categories for surface properties (albedo, roughness, emissivity, soil type). The dataset is chosen during WPS geogrid processing.

**Available datasets:**

| Dataset | Categories | Description |
|---|---|---|
| **MODIS** (default) | 20 (+1 lake) | MODIS-based IGBP land classification at 30-arc-second (~1 km). Default since WRF v3.8. Includes 1 urban category. Modern and globally consistent. |
| **USGS** | 24 (+1 lake) | USGS 24-category land use at 30-arc-second. Legacy dataset from 1990s AVHRR. 3 urban categories (commercial, residential, industrial). Some users prefer it for urban studies, though WUDAPT LCZ is now better. |
| **MODIS + LCZ** | 20 + 10 LCZ | MODIS with WUDAPT Local Climate Zone urban sub-categories. 10 built types (compact high-rise through lightweight low-rise) and 7 land types. Best for urban heat island studies. Requires `use_wudapt_lcz = 1`. |
| **NLCD** (US only) | 40 | National Land Cover Database. 30-m resolution, US-only. Most detailed land cover for CONUS simulations. Requires NLCD-to-WRF conversion. |

**Selection in geogrid:**
```
&geogrid
 geog_data_res = 'modis_lakes+15s+modis_fpar+modis_lai', 'modis_lakes+15s+modis_fpar+modis_lai',
/
```

**Switching between MODIS and USGS:**

In `GEOGRID.TBL`, the `landuse_*` entry controls which dataset is used. MODIS is the default. To use USGS:
```
 rel_path = default:landuse_30s/
```
changes to:
```
 rel_path = usgs:landuse_30s_with_lakes/
```

> **Recommendation:** Use MODIS (default) for most applications. Switch to MODIS+LCZ for urban studies. Only use USGS if you have legacy configurations that depend on the 24-category classification.

### Recommended Combinations

**General-purpose research (dx > 10 km):**
```
mp_physics         = 8      (Thompson)
cu_physics         = 1      (KF) or 3 (Grell-Freitas, scale-aware)
bl_pbl_physics     = 1      (YSU)
sf_sfclay_physics  = 1      (Revised MM5)
ra_lw_physics      = 4      (RRTMG)
ra_sw_physics      = 4      (RRTMG)
sf_surface_physics = 4      (Noah-MP)
sf_urban_physics   = 0      (off, unless studying urban areas)
gwd_opt            = 1      (gravity wave drag ON for coarse grids)
```

**Convection-permitting (dx < 4 km):**
```
mp_physics         = 10     (Morrison 2-moment)
cu_physics         = 0      (OFF)
bl_pbl_physics     = 5      (MYNN-EDMF)
sf_sfclay_physics  = 5      (MYNN)
ra_lw_physics      = 4      (RRTMG)
ra_sw_physics      = 4      (RRTMG)
sf_surface_physics = 4      (Noah-MP)
sf_urban_physics   = 0      (or 1 UCM for urban domains)
gwd_opt            = 0      (terrain well-resolved at this scale)
```

**Tropical cyclone (ocean-dominant):**
```
mp_physics         = 8      (Thompson)
cu_physics         = 3      (Grell-Freitas, scale-aware for multi-nest)
bl_pbl_physics     = 1      (YSU)
sf_sfclay_physics  = 1      (Revised MM5)
ra_lw_physics      = 4      (RRTMG)
ra_sw_physics      = 4      (RRTMG)
sf_surface_physics = 2      (Noah)
sf_ocean_physics   = 1      (mixed-layer for SST feedback)
isftcflx           = 1      (modified surface exchange coefficients for high wind)
sst_update         = 1      (time-varying SST)
```

**Urban heat island study:**
```
mp_physics         = 8      (Thompson)
cu_physics         = 0      (OFF, convection-permitting)
bl_pbl_physics     = 2      (MYJ — required for BEP/BEM)
sf_sfclay_physics  = 2      (Eta Similarity)
ra_lw_physics      = 4      (RRTMG)
ra_sw_physics      = 4      (RRTMG)
sf_surface_physics = 2      (Noah)
sf_urban_physics   = 2      (BEP multi-layer) or 3 (BEP+BEM)
use_wudapt_lcz    = 1      (Local Climate Zones for urban sub-categories)
```

**Solar energy forecasting:**
```
mp_physics         = 28     (Thompson Aerosol-Aware)
cu_physics         = 3      (Grell-Freitas)
bl_pbl_physics     = 5      (MYNN-EDMF)
sf_sfclay_physics  = 5      (MYNN)
shcu_physics       = 5      (Deng shallow cumulus — improves cloud-radiation)
ra_lw_physics      = 4      (RRTMG)
ra_sw_physics      = 4      (RRTMG)
swint_opt          = 2      (FARMS sub-timestep irradiance)
aer_opt            = 1      (aerosol climatology)
sf_surface_physics = 4      (Noah-MP)
solar_diagnostics  = 1      (output GHI, DNI, DHI)
```

---

## 6. Outputting Hidden / Non-Default Variables

Many useful variables are computed internally by WRF but not written to `wrfout` by default. There are several ways to access them — most require **no recompilation**.

### 6.1 Runtime I/O (`iofields_filename`) — No Recompile

The preferred method. Add/remove any state variable from output at runtime.

**In `namelist.input`:**
```
&time_control
 iofields_filename       = "my_iofields_d01.txt", "my_iofields_d02.txt"
 ignore_iofields_warning = .true.
/
```

**Text file syntax:**
```
op:streamtype:streamid:VAR1,VAR2,VAR3
```

| Field | Values | Meaning |
|---|---|---|
| `op` | `+` or `-` | Add or remove variable |
| `streamtype` | `h` or `i` | History (output) or Input |
| `streamid` | `0`–`24` | Stream number (0 = default wrfout) |
| variables | comma-separated | Registry variable names, **NO SPACES** |

**Examples:**

Remove variables to reduce file size:
```
-:h:0:RAINC,RAINNC,GLW,OLR,XLAT_U,XLONG_U
```

Add physics tendencies to wrfout:
```
+:h:0:RTHCUTEN,RTHBLTEN,RTHRATEN,H_DIABATIC
```

Add variables to a custom auxiliary stream:
```
+:h:7:REFL_10CM,REFD_COM,ECHOTOP
```

**Important rules:**
- No spaces between variable names (commas only)
- **256-character line limit** — split long lists across multiple lines
- Variable names must exactly match the Registry `Dname`
- Only `wrf.exe` needs to be rerun (not `real.exe`)

### 6.2 Output Streams (auxhist)

WRF supports **25 output streams** (0–24). Stream 0 = default `wrfout`. Others produce separate auxiliary files.

**Reserved streams:**

| Stream | Reserved For |
|---|---|
| 0 | Default history (`wrfout`) |
| 2 | AFWA diagnostics |
| 3 | Climate extreme diagnostics (`wrfxtrm`) |
| 22 | Height-level diagnostics |
| 23 | Pressure-level diagnostics |

**Configure a custom stream (e.g., stream 7):**
```
&time_control
 iofields_filename       = "my_fields.txt"
 auxhist7_outname        = "wrfrefl_d<domain>_<date>"
 auxhist7_interval       = 15, 15
 frames_per_auxhist7     = 4, 4
 io_form_auxhist7        = 2
/
```

### 6.3 AFWA Diagnostics — No Recompile

Enables a comprehensive suite of operational weather diagnostics.

```
&afwa
 afwa_diag_opt    = 1         ! Master switch
 afwa_buoy_opt    = 1, 1      ! CAPE, CIN, LFC, lifted index
 afwa_severe_opt  = 1, 1      ! Tornado param, LLWS, hail
 afwa_radar_opt   = 1, 1      ! Reflectivity, echo top
 afwa_ptype_opt   = 1, 1      ! Precipitation type
 afwa_vil_opt     = 1, 1      ! Vertically integrated liquid
 afwa_vis_opt     = 1, 1      ! Visibility
 afwa_cloud_opt   = 1, 1      ! Cloud cover and ceiling
 afwa_therm_opt   = 1, 1      ! Heat index, wind chill
 afwa_turb_opt    = 1, 1      ! Turbulence indices
 afwa_icing_opt   = 1, 1      ! Icing severity
/
```

**Key output variables:**

| Package | Variables |
|---|---|
| Buoyancy | `AFWA_CAPE`, `AFWA_CAPE_MU`, `AFWA_CIN`, `AFWA_CIN_MU`, `AFWA_ZLFC`, `AFWA_PLFC`, `AFWA_LIDX` |
| Severe | `AFWA_TORNADO`, `AFWA_LLWS`, `AFWA_HAIL`, `WSPD10MAX` |
| Radar | `REFD` (3D), `REFD_COM` (composite), `ECHOTOP` |
| Precip type | `AFWA_RAIN`, `AFWA_SNOW`, `AFWA_ICE`, `AFWA_FZRA` |
| Cloud | `AFWA_CLOUD`, `AFWA_CLOUD_CEIL` |

> **Constraint:** AFWA cannot compile with OpenMP (smpar). Use distributed-memory (dmpar) only.

### 6.4 Pressure-Level & Height-Level Diagnostics

**Pressure-level output (stream 23):**
```
&diags
 p_lev_diags      = 1
 num_press_levels  = 9
 press_levels      = 92500, 85000, 70000, 60000, 50000, 40000, 30000, 20000, 10000
 use_tot_or_hyd_p  = 1         ! 1=total pressure, 2=hydrostatic
/

&time_control
 auxhist23_outname    = "wrfpress_d<domain>_<date>"
 auxhist23_interval   = 60, 60
 frames_per_auxhist23 = 100, 100
 io_form_auxhist23    = 2
/
```

> Pressure values are in **Pa** (not hPa). 850 hPa = 85000 Pa.

**Height-level output (stream 22):**
```
&diags
 z_lev_diags    = 1
 num_z_levels   = 4
 z_levels       = -30, -80, -120, -150
/
```

> Negative values = above ground level (AGL). `-80` means 80 m AGL. Useful for wind energy (hub heights).

### 6.5 NWP Convective Diagnostics

```
&time_control
 nwp_diagnostics = 1
/
&physics
 do_radar_ref    = 1
/
```

Output (max/mean values between output times, reset each interval):

| Variable | Description |
|---|---|
| `WSPD10MAX` | Max 10 m wind speed |
| `UP_HELI_MAX` | Max updraft helicity (2–5 km) |
| `W_UP_MAX` | Max updraft velocity |
| `W_DN_MAX` | Max downdraft velocity |
| `GRPL_MAX` | Max column-integrated graupel |
| `REFD_MAX` | Max derived reflectivity |

### 6.6 Other Useful Namelist Switches

| What You Want | Setting | Recompile? |
|---|---|---|
| Precipitation accumulation buckets | `prec_acc_dt = 60` in `&physics` | No |
| Physics tendency output (16 vars) | `acc_phy_tend = 1` in `&physics` | No |
| Climate extremes (T2MAX, T2MIN, etc.) | `output_diagnostics = 1` + auxhist3 config | No |
| Solar diagnostics (GHI, DNI, DHI) | `solar_diagnostics = 1` in `&diags` | No |
| Incremental Analysis Update | `iau = 1` in `&time_control` (new in v4.7) | No |

### 6.7 Post-Processing for Variables Not in WRF Output

Some variables (vorticity, PV, storm-relative helicity) are **never** in wrfout. Compute them in post-processing:

```python
from wrf import getvar

avo = getvar(ncfile, "avo")        # Absolute vorticity
pvo = getvar(ncfile, "pvo")        # Potential vorticity
srh = getvar(ncfile, "helicity")   # Storm-relative helicity
cape = getvar(ncfile, "cape_2d")   # MCAPE, MCIN, LCL, LFC
dbz = getvar(ncfile, "dbz")        # Simulated reflectivity
rh  = getvar(ncfile, "rh")         # Relative humidity
```

### 6.8 The Registry System (Advanced — Requires Recompile)

For variables not exposed by any namelist switch, edit the Registry files directly.

WRF controls variable I/O through `Registry/Registry.EM_COMMON`. Each state entry has an IO column:
- `h` = write to history (wrfout)
- `i` = read from input
- `r` = write to restart
- `h0` = stream 0, `h1` = stream 1, etc.

**Example — enable a hidden variable:**

Change:
```
state  real  rthcuten  ikj  misc  1  -  r  "RTHCUTEN"  "..."  "..."
```
To:
```
state  real  rthcuten  ikj  misc  1  -  rh  "RTHCUTEN"  "..."  "..."
```

Then clean and recompile:
```bash
./clean -a && ./configure && ./compile em_real
```

### Quick Decision Guide

| What You Want | Method | Recompile? |
|---|---|---|
| Remove variables from wrfout | `iofields_filename` with `-:h:0:...` | No |
| Add existing state var to wrfout | `iofields_filename` with `+:h:0:...` | No |
| CAPE, CIN, lifted index | `&afwa` buoy/severe options | No |
| Composite reflectivity | `afwa_radar_opt=1` or `do_radar_ref=1` | No |
| Updraft helicity, max winds | `nwp_diagnostics=1` | No |
| Pressure-level interpolation | `p_lev_diags=1` + auxhist23 | No |
| Vorticity, PV, SRH | Post-processing (wrf-python) | No |
| Precipitation buckets | `prec_acc_dt=60` | No |
| Brand new variable not in Registry | Edit Registry + Fortran code | **Yes** |

---

## 7. Post-Processing & Visualization

### WRF Output Variables

Output files: `wrfout_d<domain>_<YYYY-MM-DD_HH:MM:SS>` in NetCDF format.

| Variable | Description | Units |
|---|---|---|
| T | Perturbation potential temperature (total = T + 300 K) | K |
| P, PB | Perturbation and base pressure (total = P + PB) | Pa |
| PH, PHB | Perturbation and base geopotential (height = (PH+PHB)/9.81) | m^2/s^2 |
| U, V, W | Wind components (staggered, grid-relative) | m/s |
| QVAPOR | Water vapor mixing ratio | kg/kg |
| RAINC, RAINNC | Accumulated convective / grid-scale precipitation | mm |
| T2 | 2-meter temperature | K |
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

ncfile = Dataset("wrfout_d01_2024-01-15_00:00:00")

# Diagnostic variables (wrf-python handles destaggering + unit conversion)
slp  = getvar(ncfile, "slp")        # Sea level pressure (hPa)
t2   = getvar(ncfile, "T2")         # 2-m temperature (K)
z    = getvar(ncfile, "z")          # Geopotential height (m)
p    = getvar(ncfile, "pressure")   # Pressure (hPa)

# Interpolate to pressure level
ht_500 = interplevel(z, p, 500.)

# Plot
cart_proj = get_cartopy(slp)
lats, lons = latlon_coords(slp)

fig, ax = plt.subplots(subplot_kw={"projection": cart_proj})
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
contours = ax.contourf(to_np(lons), to_np(lats), to_np(slp),
                        levels=20, transform=crs.PlateCarree(), cmap="coolwarm")
plt.colorbar(contours, ax=ax, shrink=0.7, label="SLP (hPa)")
plt.savefig("slp_plot.png", dpi=150)
```

**wrf-python diagnostics:** `slp`, `tk`, `td`, `rh`, `z`, `pressure`, `ua`, `va`, `wspd_wdir`, `cape_2d`, `cape_3d`, `dbz`, `mdbz`, `avo`, `pvo`, `eth`, `omega`, `pw`, `helicity`, `updraft_helicity`, `cloudfrac`, and more.

### Python with xWRF (Modern xarray approach)

```python
import xarray as xr
import xwrf

ds = xr.open_dataset("wrfout_d01_2024-01-15_00:00:00").xwrf.postprocess()
# CF-compliant coordinates, destaggered, supports Dask for parallel processing
```

### Other Tools

| Tool | Notes |
|---|---|
| **NCL** | Mature, publication-quality. In maintenance mode — Python recommended. |
| **VAPOR** | Interactive 3D visualization from NCAR. Imports wrfout directly. |
| **GrADS** | Via ARWpost converter. |
| **UPP** | NCEP's Unified Post Processor. Outputs GRIB2 for operational use. |

---

## 8. WRF Variants & Extensions

### 8.1 WRF-DA (Data Assimilation)

Optimally blends observations with a prior forecast to produce improved initial conditions.

**Supported methods:**
- **3DVAR** — single-time analysis with static background error covariance (`be.dat`)
- **4DVAR** — time-window analysis using adjoint/tangent-linear models (requires WRFPLUS)
- **Hybrid Ensemble-Variational** — combines variational framework with flow-dependent ensemble covariances

**Workflow:**
```
WPS + real.exe  -->  wrfinput (first guess)
     +
obsproc.exe     -->  ob.ascii (processed observations)
     +
gen_be          -->  be.dat (background error statistics, generated once)
     |
     v
da_wrfvar.exe   -->  wrfvar_output (analysis)
     |
da_update_bc.exe --> updated wrfbdy_d01
     |
     v
wrf.exe          -->  forecast (becomes next cycle's first guess)
```

**Key components:**

| Program | Purpose |
|---|---|
| `obsproc.exe` | Processes observations into WRFDA format; QC, error assignment |
| `gen_be` | Generates static background error covariance (`be.dat`) via NMC method |
| `da_wrfvar.exe` | Performs cost function minimization (the actual DA) |
| `da_update_bc.exe` | Updates lateral boundaries to match new analysis |

**Observation types supported:** SYNOP, METAR, SHIP, BUOY, radiosonde, aircraft (AMDAR/ACARS), wind profilers, GPS (PW, RO), AMVs, scatterometer winds, Doppler radar (radial velocity + reflectivity), satellite radiances (AMSU, AIRS, IASI, CrIS, ATMS, AHI, etc.)

**Compilation:**
```bash
cd WRFDA
./configure wrfda
./compile all_wrfvar    # Builds 44 executables

# For 4DVAR: compile WRFPLUS first
cd WRFPLUS && ./configure wrfplus && ./compile wrfplus
```

### 8.2 WRF-Chem (Chemistry)

Online-coupled atmospheric chemistry and aerosol modeling. Chemistry and meteorology integrated at every time step — enables aerosol-radiation and aerosol-cloud feedbacks.

**Key chemistry mechanisms (`chem_opt`):**

| Value | Gas-Phase | Aerosol | Notes |
|---|---|---|---|
| 1 | RADM2 | None | Gas-phase only |
| 2 | RADM2 | MADE/SORGAM | Modal aerosol |
| 8 | CBMZ | MOSAIC 8-bin | Sectional aerosol |
| 112 | MOZART | GOCART bulk | MOZCART |
| 202 | MOZART | MOSAIC 4-bin + VBS SOA | Full treatment |

**Required emission inputs:**
- `wrfchemi_*` — anthropogenic emissions (from **anthro_emis** preprocessor)
- `wrffirechemi_*` — fire emissions (from **fire_emiss** / FINN / GFED)
- `wrfbiochemi_*` — biogenic VOC emissions (from **MEGAN / bio_emiss**)
- Chemical IC/BC from **mozbc** (maps global model output to WRF-Chem grid)

**Compilation:**
```bash
export WRF_CHEM=1
export WRF_KPP=1          # Optional, requires FLEX library
./configure && ./compile em_real
```

### 8.3 WRF-Fire (SFIRE)

Coupled atmosphere-wildfire model using the **level-set method** with **Rothermel fire spread formula**.

- Fire grid is a refined sub-grid (e.g., 10:1 ratio via `sr_x`, `sr_y`)
- Two-way coupling: winds drive fire spread, fire heat fluxes modify atmosphere
- Uses **Anderson 13 fuel categories** (from LANDFIRE data for US domains)

**Key namelist (`&fire`):**
```
ifire              = 2           ! Fire model ON
fire_fuel_read     = -1          ! Read fuel from WPS (-1=real data)
fire_num_ignitions = 1
fire_ignition_start_x1 = 1000.
fire_ignition_start_y1 = 500.
fire_ignition_radius1  = 50.
fire_ignition_time1    = 2.
fire_wind_height       = 6.096   ! ~20 ft, matches Rothermel
```

**Compilation:** Part of standard WRF — no special flags needed. `ifire` namelist controls activation. Use `./compile em_fire` for idealized test cases.

### 8.4 WRF-Hydro

Couples atmospheric and terrestrial hydrologic processes. Disaggregates land surface states to a high-resolution routing grid (100–250 m), computes lateral water redistribution, aggregates back.

**Routing components:** overland flow (diffusive wave), subsurface lateral flow, channel routing (Muskingum-Cunge or diffusive wave), reservoir routing (level-pool), baseflow bucket model.

**NOAA's National Water Model** is built on WRF-Hydro — forecasts streamflow on 2.7 million river reaches.

**Compilation:**
```bash
# Standalone
cd trunk/NDHMS && ./configure && ./compile_offline_NoahMP.sh

# Coupled with WRF
export WRF_HYDRO=1
cd WRF && ./configure && ./compile em_real
```

### 8.5 WRF-Solar

Solar energy forecasting enhancement — improves surface irradiance prediction.

**Key outputs:** GHI (`SWDOWN`), DNI (`SWDDNI`), DHI (`SWDDIF`)

**Key settings:**
```
swint_opt = 2        ! FARMS sub-timestep irradiance (recommended)
aer_opt   = 1        ! Aerosol climatology (or 3 for full aerosol-cloud coupling)
solar_diagnostics = 1
shcu_physics = 5     ! Deng shallow cumulus (requires MYNN or MYJ PBL)
```

Part of standard WRF since v4.2. No special compile flags — activated by namelist.

### 8.6 HWRF / HAFS (Hurricane Modeling)

**HWRF** (Hurricane WRF) was NOAA's primary TC model (2007–2023). Featured movable vortex-following nests, ocean coupling (MPIPOM-TC/HYCOM), and specialized air-sea exchange at high wind speeds.

**HAFS** (Hurricane Analysis and Forecast System) replaced HWRF operationally in June 2023. Uses the **FV3 cubed-sphere** dynamical core instead of WRF-NMM. Provides ~10–15% track improvement and better rapid intensification prediction.

### 8.7 MPAS (Model for Prediction Across Scales)

Global non-hydrostatic model using **unstructured Voronoi mesh** with **smooth variable resolution** — no abrupt nesting boundaries. Conceptual successor to WRF for global/variable-resolution applications. Shares WRF's physics packages.

**Key advantages over WRF:**
- No lateral boundary conditions needed (global model)
- Smooth mesh transitions instead of sharp nesting interfaces
- No polar singularity
- Single continuous integration from global to convection-permitting scales

**Pre-built meshes:** 240 km, 120 km, 60 km, 30 km, 15 km (uniform); 60-15 km, 60-10 km, 60-3 km (variable)

### Summary Comparison

| Variant | Part of WRF? | Compile Flag / Target | Application |
|---|---|---|---|
| **WRF-DA** | Yes | `compile all_wrfvar` | Data assimilation cycling |
| **WRF-Chem** | Yes | `WRF_CHEM=1` | Air quality, aerosols |
| **WRF-Fire** | Yes (built-in) | None (namelist) | Wildfire behavior |
| **WRF-Hydro** | Separate repo | `WRF_HYDRO=1` | Streamflow, floods |
| **WRF-Solar** | Yes (v4.2+) | None (namelist) | Solar irradiance |
| **HWRF** | Separate (retired) | dmpar only | Tropical cyclones |
| **HAFS** | Separate (FV3-based) | UFS build | Tropical cyclones (current) |
| **MPAS** | Separate model | `CORE=atmosphere` | Global variable-resolution |

---

## 9. JMA NWP Settings & Recommendations

The Japan Meteorological Agency (JMA) operates one of the world's most advanced NWP suites. Their configurations provide excellent reference points for WRF users targeting the Asia-Pacific region.

### JMA Operational NWP Systems

| Model | Resolution | Domain | Vert. Levels | Model Top | Forecast Range | DA Method |
|---|---|---|---|---|---|---|
| **GSM** (Global) | 13 km (TQ959) | Global | 128 | 0.01 hPa | 264 hr | Hybrid 4D-Var + LETKF |
| **MSM** (Mesoscale) | 5 km | Japan + environs | 96 | 37.5 km | 78 hr | 4D-Var (ASUCA-based) |
| **LFM** (Local) | 2 km | Japan | 76 | 21.8 km | 18 hr | Hybrid DA |
| **MEPS** (Meso Ensemble) | 2 km | Japan + environs | 96 | 37.5 km | 39 hr | 4D-Var + ensemble perturbations |
| **GEPS** (Global Ensemble) | 25–37.5 km | Global | 128 | 0.01 hPa | 11+ days | Ensemble perturbations |

### Key JMA Design Principles for WRF Users

1. **High vertical resolution** — 76–128 levels with model top at 0.01 hPa (GSM/GEPS) or 37.5 km (MSM). This is considerably more than the WRF default of 30–45 levels.

2. **Frequent updates** — MSM runs every 3 hours, LFM runs every hour. Frequent cycling with data assimilation greatly improves short-range forecasts.

3. **Nonhydrostatic mesoscale** — The MSM and LFM use JMA's in-house nonhydrostatic model ASUCA. For WRF-based work targeting similar domains, use the ARW (nonhydrostatic) core.

4. **Scale-appropriate convection** — LFM at 2 km runs convection-permitting (no cumulus parameterization), consistent with WRF best practice for dx < 4 km.

### JMA-Inspired WRF Configuration for Japan / Western Pacific

Based on JMA's operational choices, here is a recommended WRF setup for the region:

**Domain setup (3-domain nesting):**
```
&domains
 max_dom          = 3,
 dx               = 15000, 5000, 1667,     ! 15 km -> 5 km -> ~1.7 km
 dy               = 15000, 5000, 1667,
 e_vert           = 76,    76,   76,        ! High vertical resolution (JMA MSM uses 96)
 p_top_requested  = 1000,                   ! ~50 km model top (similar to MSM's 37.5 km)
 parent_grid_ratio = 1, 3, 3,
 parent_time_step_ratio = 1, 3, 3,
/
```

**Physics (matching JMA's operational approach):**
```
&physics
 ! Microphysics: double-moment for convection-permitting
 mp_physics         = 10, 10, 10,           ! Morrison 2-moment (or 8 Thompson)

 ! Cumulus: ON for 15km, scale-aware for 5km, OFF for 1.7km
 cu_physics         = 3, 3, 0,              ! Grell-Freitas scale-aware
                                             ! Turns itself off at convection-permitting scales

 ! PBL: MYNN (good for maritime/coastal Japan)
 bl_pbl_physics     = 5, 5, 5,              ! MYNN-EDMF
 sf_sfclay_physics  = 5, 5, 5,              ! MYNN surface layer

 ! Radiation: RRTMG (JMA uses broadband radiation similar to RRTMG)
 ra_lw_physics      = 4, 4, 4,              ! RRTMG
 ra_sw_physics      = 4, 4, 4,              ! RRTMG
 radt               = 15, 5, 2,             ! ~dx in km

 ! Land surface: Noah-MP for advanced land processes
 sf_surface_physics = 4, 4, 4,              ! Noah-MP

 ! SST update (critical for maritime domains)
 sst_update         = 1,
/
```

**Data assimilation approach (JMA-inspired cycling):**

JMA uses hybrid 4D-Var for their global model and 4D-Var for their mesoscale model, both with ensemble-derived background errors. For WRF-DA equivalent:

```
# Use WRF-DA 3DVAR or hybrid as a practical equivalent
# JMA cycles every 3 hours (MSM) — WRF-DA can replicate this
# Key: assimilate radar, satellite, and surface obs for Japan domain
```

### JMA Modernization Timeline (Reference)

| Year | Change |
|---|---|
| 2023 | GSM resolution upgraded 20 km -> 13 km |
| 2022 | MSM vertical layers increased to 96; GEPS resolution enhanced |
| 2021 | GSM/GEPS vertical layers increased to 128 |
| 2020 | GSM extended to 264 hours; 4D-Var incorporated into MSM |
| 2019 | MEPS introduced; LFM extended |

> **Source:** [JMA NWP Activities](https://www.jma.go.jp/jma/en/Activities/nwp.html)

---

## 10. System Requirements

### Operating System

- **Linux** — primary supported platform (Ubuntu, CentOS, RHEL)
- macOS — supported, less common for production
- Windows — NOT natively supported (use WSL2)

### Compilers

| Compiler | Notes |
|---|---|
| **gfortran** (GNU) | Free, v9+ recommended. v4.7.1 fixes GCC 15 compatibility. |
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

> NetCDF-C and NetCDF-Fortran must be compiled with the **same compiler** as WRF.

### Environment Variables

```bash
export NETCDF=/usr/local/netcdf
export HDF5=/usr/local/hdf5
export JASPERLIB=/usr/local/jasper/lib
export JASPERINC=/usr/local/jasper/include
export WRF_DIR=/path/to/WRF
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

### Build Process (v4.7)

```bash
# 1. Configure WRF
cd WRF
./configure       # Select compiler + parallelism (e.g., GNU dmpar)
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

> **v4.7 download note:** When downloading from GitHub, use the `v4.7.x.tar.gz` file — NOT the auto-generated "Source Code" archive, which misses required submodules (NoahMP, MYNN-EDMF, MMM-physics).

---

## 11. Troubleshooting & Tips

### CFL Errors

CFL violations = numerical instability. Fixes (in order):

1. **Reduce time step** — from `6*dx_km` down to `4*dx_km` or `3*dx_km`
2. **Enable adaptive time stepping** — `use_adaptive_time_step = .true.`
3. **Enable vertical velocity damping** — `w_damping = 1` in `&dynamics`
4. **Increase sound step damping** — `epssm = 0.2` (up to 0.5 for complex terrain)
5. **Smooth topography** — `smooth_cg_topo = .true.` in `&domains`
6. **Check input data** — inspect `met_em` files for unrealistic values

### Spin-up Time

| Simulation Type | Recommended Spin-up |
|---|---|
| Short-term forecast | 6–12 hours |
| General weather | 12–24 hours |
| Precipitation studies | 24–48 hours |
| Climate downscaling | ~1 month |

### SIZE MISMATCH Errors

Domain dimensions in `namelist.input` (e_we, e_sn, e_vert) don't match the input files.

### Common Debugging Commands

```bash
grep -i cfl rsl.error.*                       # Check CFL errors
tail rsl.out.0000                              # Check completion
ncdump -h wrfout_d01_2024-01-15_00:00:00      # Inspect NetCDF headers
ncdump -v Times wrfout_d01_*                   # Check timestamps
grep -i error rsl.error.0000                   # Namelist errors
./clean -a && ./configure && ./compile em_real  # Clean rebuild
```

### Memory Issues

- Run `ulimit -s unlimited` before execution (NOT for OpenMP/smpar builds)
- Increase MPI processes rather than memory per node for large domains
- Use `io_form_restart = 102` for restart files > 4 GB

---

## 12. What's New in WRF v4.7

### Version History

| Version | Date | Type |
|---|---|---|
| **v4.7.1** | June 2, 2025 | Bug-fix (GCC 15 compatibility) |
| **v4.7.0** | April 25, 2025 | Major release |
| v4.6.1 | October 16, 2024 | Bug-fix |
| v4.6.0 | May 9, 2024 | Major release |

### Key v4.7 Additions

| Feature | Details |
|---|---|
| **UFS Double Moment microphysics** | `mp_physics=27`; 7-class, aerosol-aware CCN |
| **RCON microphysics** | `mp_physics=29`; improved warm rain/drizzle |
| **MYNN-EDMF refactor** | Now a git submodule; 10–15% faster; module names changed to `*_mynnedmf_*` |
| **Incremental Analysis Update** | `iau=1` in `&time_control`; for 3DVAR cycling |
| **Pleim-Xiu + LCZ** | Now supports 61-category MODIS LCZ dataset |
| **CMake Chem+KPP** | CMake build system now supports WRF-Chem and WRFPLUS |
| **Submodule architecture** | NoahMP, MYNN-EDMF, and MMM-physics as git submodules |

### Breaking Changes

- `mp_physics=17, 19, 21, 22` deprecated (use `mp_physics=18` with `nssl_*` flags)
- MYNN module names changed: `*_mynn_*` -> `*_mynnedmf_*`
- `mp_zero_out` now only affects `moist` array; use `mp_zero_out_all` for old behavior

---

## References

- [WRF Users Guide](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/)
- [WRF GitHub](https://github.com/wrf-model/WRF) | [WPS GitHub](https://github.com/wrf-model/WPS)
- [WRF v4.7.0 Release Notes](https://github.com/wrf-model/WRF/releases/tag/v4.7.0)
- [NCAR Geoscience Data Exchange (RDA)](https://gdex.ucar.edu/)
- [wrf-python Documentation](https://wrf-python.readthedocs.io/)
- [xWRF Documentation](https://xwrf.readthedocs.io/)
- [NCAR WRF Tutorial](https://www2.mmm.ucar.edu/wrf/users/tutorial/tutorial.html)
- [WRF Forum](https://forum.mmm.ucar.edu/)
- [WRFDA User Guide](https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/v4.3/users_guide_chap6.html)
- [WRF-Chem User Guide](https://ruc.noaa.gov/wrf/wrf-chem/Users_guide.pdf)
- [WRF-Fire Users Guide](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/fire.html)
- [WRF-Hydro Documentation](https://wrf-hydro.readthedocs.io/)
- [AFWA Diagnostics in WRF](https://www2.mmm.ucar.edu/wrf/users/docs/AFWA_Diagnostics_in_WRF.pdf)
- [JMA NWP Activities](https://www.jma.go.jp/jma/en/Activities/nwp.html)
- [MPAS-Atmosphere](https://www2.mmm.ucar.edu/projects/mpas/)
