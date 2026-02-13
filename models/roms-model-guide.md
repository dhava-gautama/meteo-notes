# ROMS: Complete Guide from Grid Setup to Visualization

> A practical guide covering the full ROMS (Regional Ocean Modeling System) workflow for coastal and regional ocean modeling.
> Covers **ROMS-Rutgers**, **ROMS-UCLA**, and **CROCO** forks.

---

## Table of Contents

1. [Overview & Key Concepts](#1-overview--key-concepts)
2. [ROMS Versions & Forks](#2-roms-versions--forks)
3. [Grid Generation](#3-grid-generation)
4. [Forcing Data](#4-forcing-data)
5. [Compilation & System Requirements](#5-compilation--system-requirements)
6. [Configuration Files](#6-configuration-files)
7. [Vertical S-Coordinate Grid](#7-vertical-s-coordinate-grid)
8. [Physics Options Reference](#8-physics-options-reference)
9. [Lateral Boundary Conditions](#9-lateral-boundary-conditions)
10. [Nesting](#10-nesting)
11. [Running the Model](#11-running-the-model)
12. [Output & Post-Processing](#12-output--post-processing)
13. [Indonesian Waters & Sunda Strait Application](#13-indonesian-waters--sunda-strait-application)
14. [Troubleshooting & Tips](#14-troubleshooting--tips)

---

## 1. Overview & Key Concepts

### What Is ROMS?

ROMS (Regional Ocean Modeling System) is a free-surface, terrain-following, primitive-equations ocean model widely used for coastal and regional ocean applications. Developed through collaboration between Rutgers University (Hernan Arango) and UCLA (Alexander Shchepetkin), ROMS simulates ocean response to wind, heating, tides, river discharge, and open-ocean boundary conditions.

### Governing Equations

ROMS solves the **hydrostatic primitive equations** under the **Boussinesq approximation** with a **free surface**:

**Horizontal momentum:**
```
du/dt + (u·∇)u − fv = −(1/ρ₀) ∂P/∂x + Fᵤ + Dᵤ
dv/dt + (u·∇)v + fu = −(1/ρ₀) ∂P/∂y + Fᵥ + Dᵥ
```

**Hydrostatic approximation:**
```
∂P/∂z = −ρg
```

**Free-surface equation:**
```
∂ζ/∂t + ∂(Hu)/∂x + ∂(Hv)/∂y = 0
```

**Tracer equations** (temperature and salinity):
```
∂T/∂t + u·∇T = Fₜ + Dₜ
∂S/∂t + u·∇S = Fₛ + Dₛ
```

**Equation of state:** ρ = f(T, S, P) — nonlinear (UNESCO/Jackett-McDougall) or linear.

### Key Features

| Feature | Description |
|---|---|
| **Vertical coordinate** | Terrain-following S-coordinate (sigma-like). Levels follow bathymetry at bottom, free surface at top. |
| **Horizontal grid** | Arakawa C-grid, orthogonal curvilinear coordinates |
| **Time stepping** | Split-explicit: fast barotropic (2D) mode separated from slow baroclinic (3D) mode |
| **Free surface** | Prognostic free-surface elevation |
| **Boussinesq** | Reference density ρ₀ ≈ 1025 kg/m³; full density retained only in pressure gradient |

### Arakawa C-Grid

```
+-------v-------+-------v-------+
|               |               |
u      rho      u      rho      u
|               |               |
+-------v-------+-------v-------+
|               |               |
u      rho      u      rho      u
|               |               |
+-------v-------+-------v-------+
  psi             psi
```

| Point Type | Location | Variables |
|---|---|---|
| **RHO-points** | Cell centers | ζ, T, S, ρ, all tracers |
| **U-points** | West/East edges | Zonal velocity (u) |
| **V-points** | South/North edges | Meridional velocity (v) |
| **PSI-points** | Cell corners | Vorticity, streamfunction |

### Split-Explicit Time Stepping

- **Barotropic (2D) mode:** Integrated with small time step (`dt / NDTFAST`) to resolve fast surface gravity waves
- **Baroclinic (3D) mode:** Integrated with larger time step (`dt`)
- A predictor-corrector coupling ensures volume conservation and suppresses mode-splitting errors

### Workflow Overview

```
[Define Domain & Generate Grid]  -->  grid.nc
        |
[Download Bathymetry (GEBCO/ETOPO)]
        |
[Smooth Bathymetry (rx0, rx1)]
        |
[Prepare Forcing Data]
   ├── Atmospheric (ERA5/GFS)     -->  forcing.nc
   ├── Tidal (TPXO9/FES2014)     -->  tides.nc
   ├── Boundary (GLORYS/HYCOM)    -->  bry.nc
   ├── Initial Conditions         -->  ini.nc
   └── Rivers (optional)          -->  rivers.nc
        |
[Configure CPP flags + ocean.in]
        |
[Compile ROMS]  -->  romsM / romsO / romsS
        |
[Run Model]  -->  history.nc, avg.nc, rst.nc
        |
[Post-Processing & Visualization]
   ├── Python (xarray, xroms, matplotlib)
   └── MATLAB (ROMS Toolbox)
```

---

## 2. ROMS Versions & Forks

### History

ROMS emerged from collaboration between Rutgers (Arango) and UCLA (Shchepetkin, McWilliams) in the late 1990s. The time-stepping kernels were identical until ~1999, when UCLA switched to predictor-corrector schemes for 3D momentum. The codes have since diverged into distinct forks.

### Comparison

| Feature | ROMS-Rutgers | ROMS-UCLA | CROCO |
|---|---|---|---|
| **Lead** | Hernan Arango (Rutgers) | Shchepetkin & McWilliams (UCLA) | IRD/INRIA/IFREMER/SHOM/CNRS (France) |
| **Repository** | [github.com/myroms/roms](https://github.com/myroms/roms) | [github.com/CESR-lab/ucla-roms](https://github.com/CESR-lab/ucla-roms) | [gitlab.inria.fr/croco-ocean](https://gitlab.inria.fr/croco-ocean) |
| **License** | MIT/X | — | MIT/X + CeCILL-C (AGRIF) |
| **2D time stepping** | LF-AM3 with feedback | FB AB3-AM4 | FB AB3-AM4 |
| **3D momentum** | Adams-Bashforth 3 | LF-AM3 | LF-AM3 |
| **Max 3D Courant** | 0.72 | 1.58 | 1.58 |
| **Data assimilation** | **4D-Var** (IS4DVAR, I4DVAR) | 3D-Var | 3D-Var |
| **Online nesting** | No (offline only) | No | **Yes (AGRIF)** |
| **Non-hydrostatic** | No | Experimental | **Yes (CROCO-NH)** |
| **Community** | myroms.org forum | — | croco-ocean.org |
| **Latest release** | Continuous (GitHub) | — | v2.1.2 (Nov 2025) |

### Which Fork to Choose

- **ROMS-Rutgers:** Best for **data assimilation** (4D-Var), largest community, most comprehensive documentation
- **ROMS-UCLA:** Best for **numerical algorithm development**; larger stable Courant numbers allow larger time steps
- **CROCO:** Best for **online nesting** (AGRIF 2-way), non-hydrostatic dynamics, French/European community

### COAWST

Not a separate fork but a coupled system built on ROMS-Rutgers by John Warner (USGS). Couples ROMS (ocean) + WRF (atmosphere) + SWAN (waves) + CSTMS (sediment) via MCT.

---

## 3. Grid Generation

### Grid File Contents

A ROMS grid is a NetCDF file containing:

| Variable | Description |
|---|---|
| `lat_rho`, `lon_rho` | Coordinates at RHO-points |
| `lat_u`, `lon_u` | Coordinates at U-points |
| `lat_v`, `lon_v` | Coordinates at V-points |
| `lat_psi`, `lon_psi` | Coordinates at PSI-points |
| `h` | Bathymetry depth (m) at RHO-points |
| `mask_rho`, `mask_u`, `mask_v`, `mask_psi` | Land/sea masks (1=water, 0=land) |
| `pm`, `pn` | Inverse grid spacing (1/Δx, 1/Δy) |
| `f` | Coriolis parameter |
| `angle` | Grid rotation angle from east |

### Grid Generation Tools

| Tool | Language | Notes |
|---|---|---|
| **ROMS-Tools** | Python | Modern package from CWorthy-ocean. Creates grids, forcing, tides, IC/BC. [github.com/CWorthy-ocean/roms-tools](https://github.com/CWorthy-ocean/roms-tools) |
| **CROCO_PYTOOLS** | Python | Full pre/post-processing for CROCO. v2.0.3. |
| **CROCO_TOOLS** | MATLAB | Classic MATLAB toolbox for CROCO/ROMS |
| **pyroms** | Python | Grid creation, nesting, pre/post-processing |
| **GridBuilder** | MATLAB | GUI-based interactive grid design |
| **EasyGrid** | MATLAB | Simple script for beginners |
| **gridgen** | C + Python | Pavel Sakov's orthogonal grid generator |

### Domain Definition

1. **Geographic extent:** Define corner points or center + extent
2. **Resolution:** Regional 1–10 km; coastal 100 m – 1 km
3. **Grid dimensions:** `Lm` × `Mm` (interior RHO-points). NetCDF dimensions are `Lm+2` × `Mm+2`
4. **Grid type:** Rectangular, curvilinear (follow coastlines), or rotated

### Bathymetry Sources

| Dataset | Resolution | Notes |
|---|---|---|
| **GEBCO 2025** | 15 arc-sec (~450 m) | Most authoritative global bathymetry |
| **ETOPO 2022** | 15 arc-sec | NOAA product, land + ocean |
| **SRTM30_PLUS** | 30 arc-sec (~900 m) | Smith & Sandwell |
| **Regional surveys** | Varies | Multi-beam, LIDAR for high-res coastal |

### Bathymetry Smoothing

Terrain-following coordinates generate **pressure gradient errors** over steep topography. Two key stiffness metrics:

**rx0 (Beckmann-Haidvogel number):**
```
rx0 = max |h(i) − h(j)| / [h(i) + h(j)]
```
Maximum fractional depth change between adjacent cells. **Target: rx0 ≤ 0.2–0.4**

**rx1 (Haney number):**
```
rx1 = max |(z_k(i) − z_k(j)) + (z_{k-1}(i) − z_{k-1}(j))| /
          |(z_k(i) − z_{k-1}(i)) + (z_k(j) − z_{k-1}(j))|
```
Slope of S-coordinate surfaces. Depends on both bathymetry and vertical stretching. **Target: rx1 ≤ 3–7.** rx1 matters more than rx0.

**Smoothing algorithms:**
- **Negative Adjustment Filter:** Iteratively adjusts cells exceeding rx0 target
- **Shapiro Filter:** Low-pass filter preserving large-scale features
- **Linear Programming (LP):** Minimum modification to achieve target stiffness (Sikirić et al., 2009)
- **Selective smoothing:** Smooth only where stiffness exceeds threshold

> **Warning:** Over-smoothing destroys shelf breaks, canyons, and critical sill depths. Balance smoothing with numerical stability.

---

## 4. Forcing Data

### 4.1 Atmospheric Forcing

**Two modes:**
1. **Direct flux forcing:** Provide pre-computed wind stress, heat flux, freshwater flux
2. **Bulk flux formulation** (recommended, `BULK_FLUXES` CPP flag): Provide meteorological variables; ROMS computes air-sea fluxes internally using COARE algorithm

**Bulk flux input variables:**

| Variable | NetCDF Name | Units | Description |
|---|---|---|---|
| U-wind (10 m) | `Uwind` | m/s | Eastward component |
| V-wind (10 m) | `Vwind` | m/s | Northward component |
| Air temperature | `Tair` | °C | At 2 m or 10 m |
| Relative humidity | `Qair` | fraction | Surface air humidity |
| Air pressure | `Pair` | mbar | Surface pressure |
| Precipitation | `rain` | kg/m²/s | Rain rate |
| Shortwave radiation | `swrad` | W/m² | Downward shortwave |
| Longwave radiation | `lwrad_down` | W/m² | Downwelling longwave |

**Atmospheric products:**

| Product | Provider | Resolution | Period | Notes |
|---|---|---|---|---|
| **ERA5** | ECMWF/Copernicus | 0.25°, hourly | 1940–present | Gold standard reanalysis. Recommended for hindcast. |
| **NCEP-CFSR/CFSv2** | NOAA | ~0.5°, hourly | 1979–present | Coupled reanalysis |
| **GFS** | NOAA/NCEP | 0.25°, 3-hourly | Forecast | For forecast applications |
| **JRA-55** | JMA | ~0.56° | 1958–present | Alternative reanalysis |

**ERA5 unit conversions for ROMS:**
- Temperature: K → °C (subtract 273.15)
- Surface pressure: Pa → mbar (÷ 100)
- Dewpoint → relative humidity
- Accumulated radiation (J/m²) → W/m² (÷ accumulation period)
- Accumulated precipitation (m) → kg/m²/s

### 4.2 Tidal Forcing

**CPP flags:**

| Flag | Purpose |
|---|---|
| `SSH_TIDES` | Tidal elevation at boundaries |
| `UV_TIDES` | Tidal barotropic velocity at boundaries |
| `TIDE_GENERATING_FORCES` | Astronomical tide-generating potential (recommended) |
| `RAMP_TIDES` | Gradual ramp-in during spin-up (recommended) |
| `ADD_FSOBC` | Add tidal elevation to boundary free-surface |
| `ADD_M2OBC` | Add tidal currents to boundary barotropic velocity |

**Global tidal models:**

| Model | Resolution | Constituents | Notes |
|---|---|---|---|
| **TPXO9** | 1/6° (~18 km); 1/30° regional | 15 (M2, S2, N2, K2, K1, O1, P1, Q1, MM, MF, M4, MN4, MS4, 2N2, S1) | Most widely used with ROMS |
| **TPXO10** | 1/30° global | Same | Latest version |
| **FES2014** | 1/16° (~7 km) | 34 constituents | Broader constituent spectrum |
| **FES2022** | 1/16° | Updated | Most recent FES release |

**Tidal forcing file:** Contains amplitude, phase (elevation and currents) for each constituent.

### 4.3 Boundary & Initial Conditions

| Product | Provider | Resolution | Period | Notes |
|---|---|---|---|---|
| **GLORYS12V1** | Copernicus Marine | 1/12° (~8 km), 50 levels | 1993–present | Eddy-resolving, NEMO-based. Most commonly used today. |
| **HYCOM** | US Navy/NOAA | 1/12° (~9 km), hybrid levels | 1994–present | Real-time analysis + reanalysis |
| **SODA** | UMD | 0.25–0.5° | 1980–present | Good for longer-term climate studies |
| **WOA** | NOAA/NCEI | 0.25–1° | Climatology | Monthly climatological T/S |

**Processing workflow:**
1. Download global product covering ROMS domain + buffer
2. Horizontally interpolate to ROMS grid
3. Vertically interpolate from z-levels to S-coordinate levels
4. Write to ROMS-formatted NetCDF

**Tools:** ROMS-Tools (Python), CROCO_PYTOOLS, model2roms, pyroms

### 4.4 River Input

Rivers are **point sources/sinks**. Two methods:

| Method | CPP/Config | Description |
|---|---|---|
| **LuvSrc** (recommended) | `LuvSrc = T` | Horizontal transport through U/V face |
| **LwSrc** | `LwSrc = T` | Vertical influx into water column |

**River forcing file variables:**

| Variable | Description |
|---|---|
| `river_Xposition` | Xi-direction grid index |
| `river_Eposition` | Eta-direction grid index |
| `river_direction` | Flow direction (0=xi, 1=eta) |
| `river_Vshape` | Vertical distribution (fractions summing to 1) |
| `river_transport` | Volume transport (m³/s), time-varying |
| `river_temp` | Water temperature (°C) |
| `river_salt` | Salinity (PSU, typically 0) |

---

## 5. Compilation & System Requirements

### Dependencies

| Library | Purpose | Notes |
|---|---|---|
| **Fortran compiler** | Primary language | Intel `ifort`/`ifx`, GNU `gfortran`, PGI, Cray |
| **C compiler** | Preprocessing, C routines | `gcc`, `icc`/`icx` |
| **NetCDF-Fortran** | I/O | Required. Set `NETCDF_INCDIR`, `NETCDF_LIBDIR` |
| **NetCDF-C** | Underlying NetCDF | Required by NetCDF-Fortran |
| **HDF5** | NetCDF-4 support | For compression and parallel I/O |
| **MPI** | Distributed-memory parallelism | OpenMPI, MPICH, Intel MPI |
| **GNU Make / CMake** | Build system | Both supported |
| **Perl** | Dependency generation | `sfmakedepend` script |

### Build System

The primary build scripts are `build_roms.sh` (bash) or `build_roms.csh` (csh) in `ROMS/Bin/`.

**Key build variables:**

| Variable | Description |
|---|---|
| `ROMS_APPLICATION` | Application name (must match header file) |
| `MY_ROOT_DIR` | Root directory of ROMS source |
| `MY_HEADER_DIR` | Directory with application header file |
| `MY_ANALYTICAL_DIR` | Custom analytical functions |
| `SCRATCH_DIR` | Build scratch directory |
| `USE_MPI` | Enable MPI (`on`) |
| `USE_NETCDF4` | Use NetCDF-4/HDF5 |
| `USE_LARGE` | Large file support (> 2 GB) |
| `FORT` | Compiler identifier (`ifort`, `gfortran`, `pgi`) |

**Build commands:**
```bash
# Standard build
./build_roms.sh

# Parallel compilation
./build_roms.sh -j 8
```

**Executables produced:**
- `romsM` — MPI version
- `romsO` — OpenMP version
- `romsS` — Serial version

### CPP Options

ROMS uses C preprocessor flags extensively, defined in an **application header file** (e.g., `my_app.h`). The master file is `ROMS/Include/cppdefs.h`.

**Example header file for realistic coastal application:**
```c
/* sunda_strait.h */

/* Dynamics */
#define UV_ADV                    /* Momentum advection */
#define UV_COR                    /* Coriolis force */
#define UV_VIS2                   /* Laplacian horizontal viscosity */
#define MIX_S_UV                  /* Mix momentum along S-surfaces */
#define SOLVE3D                   /* 3D primitive equations */
#define SALINITY                  /* Include salinity tracer */
#define NONLIN_EOS                /* Nonlinear equation of state */
#define CURVGRID                  /* Curvilinear grid */
#define SPHERICAL                 /* Spherical coordinates */
#define MASKING                   /* Land/sea masking */
#define SPLINES_VDIFF             /* Spline vertical diffusion */
#define SPLINES_VVISC             /* Spline vertical viscosity */
#define RI_SPLINES                /* Spline Richardson number */

/* Tracer advection & diffusion */
#define TS_DIF2                   /* Laplacian tracer diffusion */
#define MIX_GEO_TS                /* Mix tracers along geopotentials */

/* Surface forcing */
#define BULK_FLUXES               /* Compute fluxes via COARE bulk formulae */
#define LONGWAVE                  /* Use downwelling longwave */
#define SOLAR_SOURCE              /* Solar shortwave penetration */

/* Vertical mixing */
#define GLS_MIXING                /* GLS turbulence closure */
#define CANUTO_A                  /* Canuto A stability function */
#define N2S2_HORAVG               /* Horizontal averaging of N² and S² */

/* Tides */
#define SSH_TIDES                 /* Tidal elevation forcing */
#define UV_TIDES                  /* Tidal current forcing */
#define TIDE_GENERATING_FORCES    /* Tide-generating potential */
#define RAMP_TIDES                /* Gradual tidal ramp-in */

/* Pressure gradient */
#define DJ_GRADPS                 /* Spline density Jacobian (recommended) */

/* Bottom */
#define UV_QDRAG                  /* Quadratic bottom friction */

/* Output */
#define AVERAGES                  /* Time-averaged output */
#define DIAGNOSTICS_TS            /* Tracer budget diagnostics */
#define DIAGNOSTICS_UV            /* Momentum budget diagnostics */

/* Boundary analytical */
#define ANA_BTFLUX                /* Analytical bottom T flux (= 0) */
#define ANA_BSFLUX                /* Analytical bottom S flux (= 0) */
```

---

## 6. Configuration Files

### 6.1 ocean.in (Main Input File)

**Application & I/O:**
```
TITLE = Sunda Strait Simulation
MyAppCPP = SUNDA_STRAIT
VARNAME = /path/to/varinfo.dat
```

**Grid Dimensions:**
```
Lm == 200            ! Interior RHO-points in xi (I)
Mm == 150            ! Interior RHO-points in eta (J)
N  == 40             ! Vertical S-levels
NAT == 2             ! Active tracers (T and S)
```

**Domain Decomposition:**
```
NtileI == 4          ! Tiles in I-direction
NtileJ == 8          ! Tiles in J-direction
                     ! Total MPI tasks = 4 × 8 = 32
```

**Time Stepping:**
```
NTIMES == 86400      ! Total baroclinic time steps
DT == 60.0d0         ! Baroclinic time step (seconds)
NDTFAST == 30        ! Barotropic steps per baroclinic step
                     ! → barotropic dt = 60/30 = 2 seconds
```

**Vertical Coordinate:**
```
Vtransform == 2      ! Improved transformation (recommended)
Vstretching == 4     ! Shchepetkin 2010 (recommended)
THETA_S == 6.0d0     ! Surface stretching (0–10)
THETA_B == 2.0d0     ! Bottom stretching (0–4)
TCLINE == 50.0d0     ! Critical depth (m)
```

**Physical Parameters:**
```
RDRG2 == 3.0d-03     ! Quadratic bottom drag coefficient
Zob == 0.02d0        ! Bottom roughness (m)
VISC2 == 5.0d0       ! Laplacian horizontal viscosity (m²/s)
TNU2 == 5.0d0 5.0d0  ! Laplacian horizontal diffusion for T, S
AKT_BAK == 1.0d-6 1.0d-6  ! Background vertical diffusion for T, S
AKV_BAK == 1.0d-5    ! Background vertical viscosity
```

**Lateral Boundary Conditions** (W, S, E, N):
```
LBC(isFsur) ==   Cha     Cha     Cha     Cha        ! free-surface
LBC(isUbar) ==   Fla     Fla     Fla     Fla        ! 2D U-momentum
LBC(isVbar) ==   Fla     Fla     Fla     Fla        ! 2D V-momentum
LBC(isUvel) ==   RadNud  RadNud  RadNud  RadNud     ! 3D U-momentum
LBC(isVvel) ==   RadNud  RadNud  RadNud  RadNud     ! 3D V-momentum
LBC(isMtke) ==   Rad     Rad     Rad     Rad        ! mixing TKE
LBC(isTvar) ==   RadNud  RadNud  RadNud  RadNud     ! temperature
                  RadNud  RadNud  RadNud  RadNud     ! salinity
```

**Nudging Coefficients:**
```
TNUDG == 10.0d0 10.0d0   ! Nudging time scale (days) for T, S
ZNUDG == 10.0d0           ! Free surface nudging (days)
M2NUDG == 10.0d0          ! 2D momentum nudging (days)
M3NUDG == 10.0d0          ! 3D momentum nudging (days)
OBCFAC == 120.0d0         ! Boundary-to-interior nudging ratio
```

**Output Controls:**
```
NRREC == 0               ! Restart record (0 = new run, -1 = latest)
LcycleRST == T           ! Cycle restart file
NRST == 86400            ! Steps between restart writes
NHIS == 360              ! Steps between history writes
NAVG == 360              ! Steps between average writes
NDIA == 360              ! Steps between diagnostic writes
NDEFHIS == 0             ! Steps between new history files (0 = one file)
```

**Input File Paths:**
```
GRDNAME == /path/to/grid.nc
ININAME == /path/to/ini.nc
BRYNAME == /path/to/bry.nc
TIDENAME == /path/to/tides.nc
FRCNAME == /path/to/forcing.nc
```

**Output File Paths:**
```
HISNAME == /path/to/output_his.nc
AVGNAME == /path/to/output_avg.nc
DIANAME == /path/to/output_dia.nc
RSTNAME == /path/to/output_rst.nc
```

### 6.2 varinfo.dat (Variable Metadata)

Defines metadata for ALL ROMS I/O variables. Path set by `VARNAME` in `ocean.in`. Each entry has 7 fields: variable name, long name, units, field type, time variable, index code, and staggering type.

**Staggering types:** `r2dvar` (2D RHO), `u2dvar` (2D U), `v2dvar` (2D V), `r3dvar` (3D RHO), `u3dvar` (3D U), `v3dvar` (3D V), `w3dvar` (3D W-point).

**Example entries:**
```
'temp'
  'potential temperature'
  'Celsius'
  'temp, scalar, series'
  'ocean_time'
  'idTvar(itemp)'
  'r3dvar'

'zeta'
  'free-surface'
  'meter'
  'free-surface, scalar, series'
  'ocean_time'
  'idFsur'
  'r2dvar'
```

> Variable names in `varinfo.dat` must match names in your NetCDF input files. Rarely needs modification unless adding custom variables.

---

## 7. Vertical S-Coordinate Grid

### Transformation Equation (Vtransform)

| Value | Description |
|---|---|
| **Vtransform = 1** | Original (Song & Haidvogel 1994). In shallow water (h < hc), S-levels become nearly horizontal. Can have issues in very shallow areas. |
| **Vtransform = 2** | Improved (Shchepetkin 2005). **Recommended.** Better behavior in shallow water. z = ζ + (ζ + h)·S, where S = (hc·σ + h·C(σ))/(hc + h). |

### Stretching Function (Vstretching)

| Value | Author | Description |
|---|---|---|
| **1** | Song & Haidvogel (1994) | Original. theta_s: 0–20, theta_b: 0–1 |
| **2** | Shchepetkin (2005) | UCLA-ROMS. Uses cosh. theta_s: 0–20, theta_b: 0–1 |
| **3** | R. Geyer | Designed for high **bottom boundary layer** resolution in shallow applications |
| **4** | Shchepetkin (2010) | **Recommended.** Double stretching with independent surface/bottom refinement. theta_s: 0–10, theta_b: 0–4 |
| **5** | Souza et al. (2015) | Quadratic Legendre polynomial for higher surface resolution |

### Parameters

| Parameter | Description | Typical Range |
|---|---|---|
| `theta_s` | Surface stretching. Higher → more levels near surface. | 0–10 (Vstretching=4) |
| `theta_b` | Bottom stretching. Higher → more levels near bottom. | 0–4 (Vstretching=4) |
| `Tcline` / `hc` | Transition depth (m). Above hc: surface stretching dominates. Below: sigma-like. | 10–300 m |
| `N` | Number of vertical S-levels | 20–60 typically |

### Recommended Configurations

**Deep ocean (h > 2000 m):**
```
Vtransform = 2, Vstretching = 4
theta_s = 6–7, theta_b = 2–4, hc = 250–300 m, N = 40–60
```

**Shallow strait (h ~ 50–500 m, e.g., Sunda Strait):**
```
Vtransform = 2, Vstretching = 4
theta_s = 5–7, theta_b = 1–3, hc = 10–50 m, N = 30–50
```

**Rule of thumb:** Set `hc` close to or smaller than the minimum depth in your domain.

> **Critical:** When you change theta_s, theta_b, or Vstretching, ALL depth-dependent input files (IC, BC, climatology) must be regenerated.

---

## 8. Physics Options Reference

### 8.1 Vertical Mixing / Turbulence Closure

Four mutually exclusive frameworks:

#### GLS (Generic Length Scale) — `GLS_MIXING`

Two-equation turbulence model (Warner et al., 2005) that encompasses classical closures:

| Closure | GLS_P | GLS_M | GLS_N |
|---|---|---|---|
| **k-epsilon** | 3.0 | 1.5 | -1.0 |
| **k-omega** | -1.0 | 0.5 | -1.0 |
| **gen** (Umlauf & Burchard) | 2.0 | 1.0 | -0.67 |

**Sub-option CPP flags:**

| Flag | Description |
|---|---|
| `CANUTO_A` | Canuto A stability function (recommended) |
| `CANUTO_B` | Canuto B stability function |
| `KANTHA_CLAYSON` | Kantha-Clayson stability function |
| `CHARNOK` | Charnock surface roughness from wind |
| `CRAIG_BANNER` | Craig & Banner surface TKE flux |
| `N2S2_HORAVG` | Horizontal averaging of buoyancy and shear |
| `K_C2ADVECTION` | 2nd-order centered TKE/GLS advection |
| `K_C4ADVECTION` | 4th-order centered TKE/GLS advection |

**GLS parameters in ocean.in:** `GLS_Kmin`, `GLS_Pmin`, `GLS_CMU0`, `GLS_C1`, `GLS_C2`, `GLS_C3M`, `GLS_C3P`, `GLS_SIGK`, `GLS_SIGP`

#### KPP (K-Profile Parameterization) — `LMD_MIXING`

Large, McWilliams, Doney (1994). First-order closure matching surface boundary layer and interior parameterizations.

| Flag | Description |
|---|---|
| `LMD_SKPP` | Surface boundary layer KPP |
| `LMD_BKPP` | Bottom boundary layer KPP |
| `LMD_RIMIX` | Shear instability (Richardson number) |
| `LMD_CONVEC` | Convective mixing |
| `LMD_DDMIX` | Double-diffusive mixing |
| `LMD_NONLOCAL` | Nonlocal transport (counter-gradient) |
| `LMD_SHAPIRO` | Shapiro filter for boundary layer depth |

**Typical KPP configuration:**
```
LMD_MIXING + LMD_SKPP + LMD_BKPP + LMD_RIMIX + LMD_CONVEC + LMD_NONLOCAL
```

#### MY2.5 — `MY25_MIXING`

Mellor-Yamada Level 2.5 classic two-equation closure.

#### BVF Mixing — `BVF_MIXING`

Simple Brunt-Väisälä frequency-based mixing. Unstable: K = 0.1 m²/s. Stable: K proportional to 1/√N².

### 8.2 Advection Schemes

**Tracer advection** (since ROMS 3.7, set in ocean.in via `Hadvection`/`Vadvection`; older versions use CPP flags):

| Algorithm | Keyword | CPP Flag | Notes |
|---|---|---|---|
| 4th-order Akima | `A4` | `TS_A4HADVECTION` | Horizontal and vertical |
| 2nd-order centered | `C2` | `TS_C2HADVECTION` | |
| 4th-order centered | `C4` | `TS_C4HADVECTION` | |
| 3rd-order upstream | `U3` | `TS_U3HADVECTION` | Default for most applications |
| 3rd-order upstream split | `SU3` | `TS_U3ADV_SPLIT` | |
| MPDATA | `MPDATA` | `TS_MPDATA` | Recursive flux-corrected |
| HSIMT (TVD) | `HSIMT` | `TS_HSIMT` | Monotonic; same for H & V |

**Momentum advection (CPP flags):**

| Flag | Description |
|---|---|
| `UV_ADV` | Enable momentum advection |
| `UV_C2ADVECTION` | 2nd-order centered |
| `UV_C4ADVECTION` | 4th-order centered |
| `UV_U3ADV_SPLIT` | 3rd-order upstream split |
| `UV_SADVECTION` | Splines vertical advection |

**CROCO additional advection options:**

| Flag | Description |
|---|---|
| `TS_HADV_RSUP3` (default) | Split & rotated 3rd-order upstream — recommended |
| `TS_HADV_WENO5` | 5th-order WENOZ (quasi-monotonic) |
| `UV_HADV_UP3` (default) | 3rd-order upstream momentum |
| `UV_HADV_WENO5` | 5th-order WENOZ momentum |
| `VADV_ADAPT_IMP` | Adaptive implicit vertical advection (unconditionally stable) |

### 8.3 Pressure Gradient Schemes

Critical for terrain-following coordinates over steep topography:

| CPP Flag | Description |
|---|---|
| **`DJ_GRADPS`** | Spline density Jacobian (Shchepetkin & McWilliams 2003). **Recommended, default.** 3rd-order accurate, robust for steep slopes. |
| `PJ_GRADP` | Finite volume Pressure Jacobian |
| `PJ_GRADPQ2` | Quartic-2 Pressure Jacobian |
| `PJ_GRADPQ4` | Quartic-4 Pressure Jacobian |
| `WJ_GRADP` | Weighted density Jacobian |

### 8.4 Equation of State

| CPP Flag | Description |
|---|---|
| `NONLIN_EOS` | Nonlinear EOS (UNESCO/Jackett-McDougall 1995): ρ = f(T, S, P). **Recommended for realistic simulations.** |
| (not defined) | Linear EOS: ρ = ρ₀[1 − α(T−T₀) + β(S−S₀)] |

### 8.5 Bottom Boundary Layer

| CPP Flag | Description |
|---|---|
| `SSW_BBL` | Sherwood-Signell-Warner wave-current BBL |
| `SG_BBL` | Styles & Glenn wave-current BBL |
| `MB_BBL` | Meinte Blaas (Soulsby 1995) wave-current BBL |

**Bottom friction (independent of BBL):**

| Flag | Description |
|---|---|
| `UV_QDRAG` | Quadratic bottom friction (most common) |
| `UV_LDRAG` | Linear bottom friction |
| `UV_LOGDRAG` | Logarithmic bottom friction |

### 8.6 Horizontal Mixing

| Flag | Description |
|---|---|
| `UV_VIS2` | Laplacian horizontal viscosity |
| `UV_VIS4` | Biharmonic horizontal viscosity |
| `TS_DIF2` | Laplacian horizontal tracer diffusion |
| `TS_DIF4` | Biharmonic horizontal tracer diffusion |
| `MIX_S_UV` | Mix momentum along S-surfaces |
| `MIX_GEO_UV` | Mix momentum along geopotentials |
| `MIX_S_TS` | Mix tracers along S-surfaces |
| `MIX_GEO_TS` | Mix tracers along geopotentials |
| `MIX_ISO_TS` | Mix tracers along isopycnals |
| `UV_SMAGORINSKY` | Smagorinsky-type viscosity |
| `TS_SMAGORINSKY` | Smagorinsky-type diffusion |

---

## 9. Lateral Boundary Conditions

Since ROMS 3.6, boundary conditions are specified via the `LBC` keyword in `ocean.in` per variable and per boundary edge (West, South, East, North).

### Available LBC Types

| Keyword | Name | Best For |
|---|---|---|
| `Cha` | Chapman (implicit) | Free-surface. Outgoing signals leave at √(gH). |
| `Fla` | Flather | 2D velocity. Radiates deviations at gravity wave speed. |
| `Cla` | Clamped | All variables. Dirichlet-type. Can cause reflections. |
| `Clo` | Closed | Land boundaries. Wall, no flow. |
| `Gra` | Gradient | Zero-gradient (Neumann). |
| `Rad` | Radiation | Orlanski-type. Allows outgoing waves to pass. |
| `RadNud` | Radiation + Nudging | **Recommended for 3D variables.** Radiation on outflow, nudge to data on inflow. |
| `Per` | Periodic | Must match in pairs (W-E or S-N). |
| `Nes` | Nested | For nesting applications. |
| `Red` | Reduced Physics | 2D momentum alternative. |
| `Shc` | Shchepetkin | 2D momentum alternative. |

### Typical Configuration

```
LBC(isFsur) ==   Cha     Cha     Cha     Cha        ! free-surface
LBC(isUbar) ==   Fla     Fla     Fla     Fla        ! 2D U
LBC(isVbar) ==   Fla     Fla     Fla     Fla        ! 2D V
LBC(isUvel) ==   RadNud  RadNud  RadNud  RadNud     ! 3D U
LBC(isVvel) ==   RadNud  RadNud  RadNud  RadNud     ! 3D V
LBC(isTvar) ==   RadNud  RadNud  RadNud  RadNud     ! tracers
```

---

## 10. Nesting

### Types

| Type | Description |
|---|---|
| **Refinement** | Fine grid embedded within coarse grid (like WRF nesting) |
| **Composition** | Non-overlapping grids connected along edges |
| **Composite** | Overlapping grids, possibly non-aligned |

### One-Way vs Two-Way

- **One-way** (`ONE_WAY` CPP): Coarse → fine only. Simple but can cause boundary inconsistencies.
- **Two-way** (default): Fine feeds back to coarse. ROMS philosophy: refinement should always be two-way.

### Key Points

- Only **odd** refinement ratios: 3, 5, 7, 9 (recommended: 3 or 5)
- Requires `NESTING` + `REFINEMENT` (or `COMPOSITION`) CPP flags
- Contact file generated by `contact.m` (MATLAB) or equivalent Python tools
- Set `Ngrids` in `ocean.in` for number of grids

### Online vs Offline Nesting

| Method | Description |
|---|---|
| **Online** (ROMS built-in) | Parent and child run in same executable simultaneously |
| **Offline** | Parent runs first; output interpolated to child BC. Simpler, no two-way feedback. |
| **AGRIF** (CROCO only) | Online adaptive refinement. **Only available in CROCO.** |

---

## 11. Running the Model

### Parallel Execution

| Mode | Description |
|---|---|
| Serial | Single processor |
| OpenMP | Shared-memory threads on single node |
| **MPI** | Distributed-memory across nodes (most common) |

```bash
# MPI execution (total tasks = NtileI × NtileJ)
mpirun -np 32 romsM ocean.in

# OpenMP
export OMP_NUM_THREADS=8
./romsO ocean.in
```

**Domain decomposition tips:**
- Total MPI processes must equal `NtileI × NtileJ`
- Make tiles approximately square for optimal communication
- Minimum ~10 interior points per tile dimension

### Timestep Selection

**Barotropic CFL:**
```
Cg = dtfast × √(g × h_max) × √(1/dx² + 1/dy²) < 1
```

**Baroclinic CFL:**
```
C = dt × c_internal × √(1/dx² + 1/dy²) < 1
```

where `c_internal ≈ N·H/π` for constant buoyancy frequency N.

**NDTFAST guidelines:**
- Deep ocean: NDTFAST > 30
- Shallow water (< 100 m): NDTFAST ~ 20
- NDTFAST = 10 is generally too small

**Typical dt values:**

| Resolution | dt (seconds) | NDTFAST |
|---|---|---|
| 10 km | 300–600 | 30–60 |
| 3 km | 60–180 | 30 |
| 1 km | 20–60 | 20–30 |
| 100 m | 2–10 | 15–20 |

### Restart

| Parameter | Description |
|---|---|
| `NRST` | Steps between restart file writes |
| `NRREC` | Record to read (0 = new run, -1 = latest) |
| `PERFECT_RESTART` / `EXACT_RESTART` (CPP) | Saves 2 consecutive records for exact LF-AM3 restart |

**Restart procedure:**
1. Set `NRREC = -1` in `ocean.in`
2. Point `ININAME` to the restart file
3. If using `EXACT_RESTART`, set `NRREC` to the second record
4. If originally compiled with `ANA_INITIAL`, recompile with `#undef ANA_INITIAL`

> **On blowup:** ROMS writes a 3rd record to the restart file capturing the state at blowup — useful for debugging.

---

## 12. Output & Post-Processing

### Output File Types

| Type | Frequency Parameter | CPP Flag | Description |
|---|---|---|---|
| **History** | `NHIS` | Always available | Instantaneous snapshots |
| **Averages** | `NAVG`, `NTSAVG` | `AVERAGES` | Time-averaged fields |
| **Diagnostics** | `NDIA`, `NTSDIA` | `DIAGNOSTICS_TS`, `DIAGNOSTICS_UV` | Budget terms |
| **Restart** | `NRST` | Always available | Full state for restart |
| **Stations** | `NSTA` | `STATIONS` | Time series at specified points |
| **Floats** | `NFLT` | `FLOATS` | Lagrangian particle trajectories |

### Key Output Variables

| Variable | Description | Grid |
|---|---|---|
| `zeta` | Free-surface height (m) | 2D, RHO |
| `ubar` / `vbar` | Depth-averaged velocity (m/s) | 2D, U/V |
| `u` / `v` | 3D velocity (m/s) | 3D, U/V |
| `temp` | Temperature (°C) | 3D, RHO |
| `salt` | Salinity (PSU) | 3D, RHO |
| `AKv` | Vertical viscosity (m²/s) | 3D, W |
| `AKt` | Vertical diffusivity (m²/s) | 3D, W |
| `omega` | S-coordinate vertical velocity | 3D, W |
| `w` | True vertical velocity (m/s) | 3D, W |

**Diagnostic terms** (when enabled): time rate, horizontal advection, vertical advection, horizontal diffusion, vertical diffusion, pressure gradient, Coriolis for both tracers and momentum.

### Post-Processing Tools

**Python:**

| Tool | Description |
|---|---|
| **[xroms](https://github.com/xoceanmodel/xroms)** | xarray-based ROMS tools. Auto-computes z-coordinates with correct grid metrics. |
| **[xarray](https://docs.xarray.dev/)** | Native ROMS NetCDF support. Use with xgcm for staggered grids. |
| **[xcroco](https://github.com/jaard/xcroco)** | Python tools for CROCO/ROMS output |
| **matplotlib + cartopy** | Plotting |

**MATLAB:**
- ROMS Toolbox (from myroms.org SVN)
- ROMSTOOLS

**Example Python visualization:**
```python
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Open ROMS history file
ds = xr.open_dataset("output_his.nc")

# Plot SST (surface temperature)
sst = ds["temp"].isel(ocean_time=-1, s_rho=-1)  # Last time, surface level
lon = ds["lon_rho"]
lat = ds["lat_rho"]

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(10, 8))
ax.coastlines(resolution="10m")
ax.add_feature(cfeature.LAND, color="lightgray")
cf = ax.pcolormesh(lon, lat, sst, cmap="RdYlBu_r", transform=ccrs.PlateCarree())
plt.colorbar(cf, ax=ax, label="SST (°C)")
plt.title("ROMS Sea Surface Temperature")
plt.savefig("roms_sst.png", dpi=150, bbox_inches="tight")
```

---

## 13. Indonesian Waters & Sunda Strait Application

### Sunda Strait Characteristics

| Feature | Value |
|---|---|
| Location | Between Java and Sumatra, connecting Java Sea to Indian Ocean |
| Narrowest width | ~24 km |
| Maximum depth | ~80 m (sill) to ~200 m (deepest) |
| Tidal regime | Strong semidiurnal and diurnal; mixed tides |
| Key features | Anak Krakatau volcano, strong tidal currents, internal waves |
| Indonesian Throughflow | Minor pathway (Lombok, Makassar are primary) |

### Published ROMS Studies in Indonesian Seas

1. **Robertson & Ffield (2008)** — Barotropic + baroclinic tides at ~5 km. 4 constituents (M2, S2, K1, O1). Active internal tide generation at straits. (*JGR Oceans*)
2. **Nugroho et al. (2017)** — Explicit tides essential for correct SST simulation. (*Marine Pollution Bulletin*)
3. **Lagrangian transport studies (~3 km)** — Water sources through Lombok, Ombai, Timor passages. (*Frontiers in Marine Science*, 2023)

### Recommended Configuration for Sunda Strait

**Grid:**
- Parent domain: Indonesian Seas, ~3–5 km resolution
- Nested child: Sunda Strait focus, ~500 m – 1 km resolution
- Bathymetry: GEBCO 2025 (15 arc-sec), supplemented with local survey data if available
- Careful smoothing: preserve the ~80 m sill depth while keeping rx0 < 0.3, rx1 < 7

**Vertical grid:**
```
Vtransform  = 2
Vstretching = 4
N           = 40          ! 40 levels for resolving tidal mixing
theta_s     = 6.0         ! Strong surface refinement
theta_b     = 2.0         ! Moderate bottom refinement
Tcline      = 20.0        ! Small, since minimum depth is shallow
```

**Physics:**
```c
#define GLS_MIXING            /* GLS k-epsilon for tidal mixing */
#define CANUTO_A
#define N2S2_HORAVG
#define NONLIN_EOS
#define DJ_GRADPS             /* Essential for steep bathymetry */
#define UV_QDRAG              /* Quadratic bottom friction */
#define BULK_FLUXES
```

**Tides (critical for Sunda Strait):**
```c
#define SSH_TIDES
#define UV_TIDES
#define TIDE_GENERATING_FORCES
#define RAMP_TIDES
```

Use at least **8 constituents**: M2, S2, N2, K2 (semidiurnal) + K1, O1, P1, Q1 (diurnal). Indonesian seas have strong mixed tidal signals.

**Forcing:**
- Atmospheric: ERA5 (hourly, 0.25°)
- Tidal: TPXO10 or FES2022
- IC/BC: GLORYS12V1 (1/12° eddy-resolving)
- Rivers: Ciujung River (Banten), though minor influence in strait

**Boundaries:**
```
LBC(isFsur) ==   Cha     Cha     Cha     Cha
LBC(isUbar) ==   Fla     Fla     Fla     Fla
LBC(isVbar) ==   Fla     Fla     Fla     Fla
LBC(isUvel) ==   RadNud  RadNud  RadNud  RadNud
LBC(isVvel) ==   RadNud  RadNud  RadNud  RadNud
LBC(isTvar) ==   RadNud  RadNud  RadNud  RadNud
```

**Timestep (for ~1 km resolution):**
```
DT      = 30.0           ! 30 seconds baroclinic
NDTFAST = 30             ! 1-second barotropic
```

### Known Challenges

| Challenge | Solution |
|---|---|
| **Steep bathymetry** | Use `DJ_GRADPS`, smooth carefully preserving sill depths, use Vtransform=2/Vstretching=4 |
| **Narrow strait (few grid cells)** | Use nesting to increase resolution to ~500 m at the strait |
| **Strong tidal mixing** | GLS k-epsilon performs well. Need adequate vertical resolution (N ≥ 30). |
| **Internal tides** | Generated at sills and ridges. Need sufficient vertical levels and small enough dt. |
| **Pressure gradient errors** | Monitor with resting ocean test (stratified, no forcing). Spurious currents < 1 cm/s is acceptable. |
| **Wetting/drying** | Enable `WET_DRY` CPP flag if very shallow tidal flats are present near the strait |

---

## 14. Troubleshooting & Tips

### Common Errors

| Problem | Causes | Solutions |
|---|---|---|
| **BLOWUP** | Steep bathymetry, CFL violation, bad forcing, wetting/drying | Smooth bathymetry, reduce dt, check boundary data, increase viscosity |
| **CFL violation** | Timestep too large | Reduce dt; increase NDTFAST if barotropic CFL violated |
| **Pressure gradient error** | Steep topography, poor vertical stretching | Smooth bathymetry (rx0 < 0.4, rx1 < 7), use DJ_GRADPS, use Vtransform=2/Vstretching=4 |
| **Spurious deep currents** | Violated hydrostatic consistency | Smooth bathymetry, check rx1 |
| **NetCDF errors** | Time units mismatch, file permissions | Ensure consistent time reference (e.g., "seconds since 1948-01-01 00:00:00") |
| **Model drift** | Poor boundary conditions | Use RadNud boundaries, check nudging timescales, verify forcing quality |

### Debugging Strategy

1. Check standard output log for grid point with extreme values
2. Examine 3rd restart record (blowup state) to locate problematic region
3. Investigate bathymetry, sigma-level distribution, and BC at that location
4. Run a **resting ocean test** (stratified, no forcing) to check pressure gradient errors

### Performance Tips

- **Tile decomposition:** Make tiles approximately square; avoid < 10 interior points per dimension
- **Output frequency:** Reduce history output for long runs — I/O can dominate
- **NetCDF-4 compression:** Enable `USE_NETCDF4` for significantly smaller output files
- **Single precision output:** Use `OUT_DOUBLE = F` in `ocean.in` to halve output size

---

## References

### Foundational Papers

- Shchepetkin & McWilliams (2005): "The Regional Oceanic Modeling System (ROMS): A split-explicit, free-surface, topography-following-coordinate oceanic model." *Ocean Modelling*, 9, 347–404
- Shchepetkin & McWilliams (2003): "A method for computing horizontal pressure-gradient force in an oceanic model with a nonaligned vertical coordinate." *JGR*, 108(C3), 3090
- Haidvogel et al. (2008): "Ocean forecasting in terrain-following coordinates." *J. Comput. Phys.*, 227, 3595–3624
- Large, McWilliams & Doney (1994): "Oceanic vertical mixing: A review and a model with a nonlocal boundary layer parameterization." *Rev. Geophys.*, 32(4), 363–403
- Warner et al. (2005): "Performance of four turbulence closure models implemented using a generic length scale method." *Ocean Modelling*, 8, 81–113
- Warner et al. (2010): "Development of a COAWST Modeling System." *Ocean Modelling*, 35, 230–244
- Debreu et al. (2012): "Two-way nesting in split-explicit ocean models." *Ocean Modelling*, 49–50, 1–21
- Song & Haidvogel (1994): "A semi-implicit ocean circulation model using a generalized topography-following coordinate system." *J. Comput. Phys.*, 115, 228–244

### Indonesian Seas

- Robertson & Ffield (2008): "Baroclinic tides in the Indonesian seas." *JGR Oceans*, 113, C07031
- Nugroho et al. (2017): "Modelling explicit tides in the Indonesian seas." *Marine Pollution Bulletin*
- Nagai et al. (2015): "Internal tides and vertical mixing in the Indonesian Archipelago." *JGR Oceans*

### Online Resources

- [ROMS-Rutgers Wiki & Forum (myroms.org)](https://www.myroms.org)
- [ROMS-Rutgers GitHub](https://github.com/myroms/roms)
- [ROMS-UCLA GitHub](https://github.com/CESR-lab/ucla-roms)
- [CROCO Website](https://www.croco-ocean.org)
- [CROCO Documentation](https://croco-ocean.gitlabpages.inria.fr/croco_doc/)
- [CROCO GitLab](https://gitlab.inria.fr/croco-ocean)
- [ROMS-Tools (Python)](https://github.com/CWorthy-ocean/roms-tools)
- [model2roms](https://github.com/trondkr/model2roms)
- [xroms](https://github.com/xoceanmodel/xroms)
- [xcroco](https://github.com/jaard/xcroco)
- [Copernicus Marine (GLORYS)](https://data.marine.copernicus.eu)
- [ERA5 Climate Data Store](https://cds.climate.copernicus.eu)
- [GEBCO Bathymetry](https://www.gebco.net)
- [TPXO Tidal Models](https://www.tpxo.net)
- [COAWST](https://code.usgs.gov/coawstmodel/COAWST)
