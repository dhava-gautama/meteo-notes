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

---

## 5. Physics Options Reference

### Microphysics (`mp_physics`)

| Value | Scheme | Notes |
|---|---|---|
| 1 | Kessler | Simple warm-rain, no ice |
| 2 | Purdue Lin | 6-class with ice, snow, graupel |
| 6 | WSM6 | Ice, snow, graupel |
| 8 | **Thompson** | Popular, well-tested |
| 10 | **Morrison 2-moment** | Double-moment, widely used for research |
| 18 | **NSSL 2-moment** | Consolidated in v4.6+; controls via `nssl_*` flags |
| 27 | **UFS Double Moment** | New in v4.7; 7-class, aerosol-aware CCN (from S.-Y. Hong) |
| 28 | **Thompson Aerosol-Aware** | Thompson with aerosol interactions |
| 29 | **RCON** | New in v4.7; improved warm rain / drizzle (Conrick et al. 2023) |
| 50–53 | P3 variants | Predicted particle properties |

> **NSSL consolidation (v4.6+):** `mp_physics=17, 19, 21, 22` are deprecated. Use `mp_physics=18` with `nssl_2moment_on`, `nssl_ccn_on`, `nssl_hail_on` flags instead.

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
| 5 | **MYNN-EDMF** | Local, TKE | 5 (MYNN) |
| 7 | ACM2 | Hybrid | 7 (Pleim-Xiu) or 1 |
| 8 | BouLac | Local, TKE | 1 or 2 |

> **v4.7 note:** MYNN-EDMF is now a git submodule, refactored to k-only scheme (10–15% faster). Module names changed from `*_mynn_*` to `*_mynnedmf_*`.

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
| 4 | **Noah-MP** | 4 | Most advanced, multi-physics (now a git submodule) |
| 5 | CLM4 | 10 | Community Land Model |
| 7 | Pleim-Xiu | 2 | v4.7: now supports 61-category MODIS LCZ |

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
