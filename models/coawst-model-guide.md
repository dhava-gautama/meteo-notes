# COAWST: Complete Guide to Coupled Ocean-Atmosphere-Wave-Sediment Modeling

> A practical guide covering the COAWST (Coupled Ocean-Atmosphere-Wave-Sediment Transport) system — coupling WRF, ROMS, and SWAN for coastal and maritime applications.
> Based on **COAWST v3.8**. Includes a standalone SWAN reference and a Sunda Strait application guide.

---

## Table of Contents

1. [Overview](#1-overview)
2. [SWAN Wave Model Reference](#2-swan-wave-model-reference)
3. [Coupling Mechanisms & Field Exchanges](#3-coupling-mechanisms--field-exchanges)
4. [Compilation](#4-compilation)
5. [SCRIP Weight Generation](#5-scrip-weight-generation)
6. [Configuration Files](#6-configuration-files)
7. [Running COAWST](#7-running-coawst)
8. [Sunda Strait Application: Building Your Own Ina-CAWO](#8-sunda-strait-application-building-your-own-ina-cawo)
9. [Troubleshooting & Tips](#9-troubleshooting--tips)
10. [SWAN Grid Pre-Processing](#10-swan-grid-pre-processing)
11. [Output Files & Post-Processing](#11-output-files--post-processing)
12. [Verifying Coupling Works](#12-verifying-coupling-works)
13. [Spin-Up Strategies](#13-spin-up-strategies)
14. [Restart & Cycling](#14-restart--cycling)
15. [Sediment Transport (CSTMS)](#15-sediment-transport-cstms)
16. [Data Download Sources](#16-data-download-sources)
17. [Project Directory Organization](#17-project-directory-organization)
18. [Running a Test Case (Inlet_test)](#18-running-a-test-case-inlet_test)
19. [SWAN Boundary File Preparation](#19-swan-boundary-file-preparation)
20. [Nesting in COAWST](#20-nesting-in-coawst)
21. [WRF Setup in COAWST (8-Step Workflow)](#21-wrf-setup-in-coawst-8-step-workflow)
22. [WRF Map Projections](#22-wrf-map-projections)
23. [WaveWatch III in COAWST](#23-wavewatch-iii-in-coawst)
24. [ROMS Build System & External Libraries](#24-roms-build-system--external-libraries)
25. [Parallel I/O (PIO)](#25-parallel-io-pio)
26. [Vegetation Module](#26-vegetation-module)
27. [InWave (Infragravity Waves)](#27-inwave-infragravity-waves)
28. [Biology & Biogeochemistry](#28-biology--biogeochemistry)
29. [COAWST MATLAB Tools Reference](#29-coawst-matlab-tools-reference)

---

## 1. Overview

### What Is COAWST?

COAWST (Coupled Ocean-Atmosphere-Wave-Sediment Transport) is an open-source coupled modeling system developed by **Dr. John C. Warner** at the **USGS Woods Hole Coastal and Marine Science Center**. It dynamically couples atmosphere, ocean, wave, and sediment transport models for coastal applications.

COAWST is the open-source framework behind BMKG's operational **Ina-CAWO** system for Indonesian maritime forecasting.

### Component Models

| Component | Model | Purpose |
|---|---|---|
| Atmosphere | **WRF** | Mesoscale weather simulation |
| Ocean | **ROMS** | 3D ocean hydrodynamics, currents, temperature, salinity |
| Waves | **SWAN** | Spectral wave generation and propagation |
| Waves (alt.) | **WAVEWATCH III** | Global/regional wave modeling |
| Sediment | **CSTMS** | Suspended/bedload sediment transport, morphodynamics |
| Hydrology | **WRF-Hydro** | Hydrological routing |
| Sea Ice | **CICE** | Sea ice dynamics |

### Coupling Infrastructure

| Framework | Description |
|---|---|
| **MCT** (Model Coupling Toolkit) | Production standard. Parallel communication, domain decomposition mapping, sparse matrix interpolation. From Argonne National Lab. |
| **ESMF/NUOPC** | Newer alternative. Earth System Modeling Framework with NUOPC layer. Requires ESMF v8+. |

### Repository

| | |
|---|---|
| **GitHub** | [github.com/DOI-USGS/COAWST](https://github.com/DOI-USGS/COAWST) |
| **USGS GitLab** | [code.usgs.gov/coawstmodel/COAWST](https://code.usgs.gov/coawstmodel/COAWST) |
| **Version** | 3.8 (main branch) |
| **License** | Public Domain (USGS) / CC0 1.0 Universal |
| **Users** | 800+ registered worldwide |

### Workflow Overview

```
[WRF Setup]           [ROMS Setup]         [SWAN Setup]
WPS → wrfinput,       Grid → ini.nc,       Grid → coord.grd,
      wrfbdy           bry.nc, frc.nc       bathy.bot
    |                     |                     |
    +----------+----------+----------+----------+
               |                     |
     [SCRIP Weight Generation]
     scrip_coawst → weights.nc
               |
     [Create coupling_*.in]
     (processor allocation, exchange intervals)
               |
     [Create Application Header (.h)]
     (CPP flags: ROMS_MODEL, SWAN_MODEL, WRF_MODEL, MCT_LIB)
               |
     [Compile: build_coawst.sh]
               |
     [Run: mpirun -np N coawstM coupling.in]
               |
     [Output: ocean_his.nc, wrfout_d01, swan output]
```

---

## 2. SWAN Wave Model Reference

### What Is SWAN?

SWAN (Simulating WAves Nearshore) is a third-generation spectral wave model developed at **TU Delft** (Netherlands). It solves the **spectral action balance equation** for wind-generated waves in coastal regions. Current version: **41.51**.

**Action balance equation:**
```
∂N/∂t + ∂(cₓN)/∂x + ∂(cᵧN)/∂y + ∂(c_σN)/∂σ + ∂(c_θN)/∂θ = S_tot/σ
```

Where N = E/σ is action density (energy/frequency), c are propagation velocities, and S_tot is the total source/sink term.

**Source terms:** S_tot = S_in (wind input) + S_nl4 (quadruplet interactions) + S_nl3 (triad interactions) + S_ds,w (whitecapping) + S_ds,br (depth-induced breaking) + S_ds,b (bottom friction)

### SWAN vs WAM vs WAVEWATCH III

| Feature | SWAN | WAM | WW3 |
|---|---|---|---|
| **Numerics** | Implicit (efficient in shallow water) | Explicit | Explicit |
| **Best for** | Coastal, nearshore, estuaries | Deep ocean, global | Deep ocean, global |
| **Grid types** | Regular, curvilinear, unstructured | Regular | Regular, unstructured |
| **Developer** | TU Delft | ECMWF heritage | NOAA/NCEP |
| **License** | GPL | — | Public domain |

### Grid Types

| Type | Description |
|---|---|
| **Regular** | Uniform rectangular. `CGRID REGular ...` |
| **Curvilinear** | Curved quadrilateral. Good for following coastlines. |
| **Unstructured** | Triangular mesh. Optimal resolution adaptation. |

### Physics Options

#### Wind Input (GEN3 command)

| Formulation | Description |
|---|---|
| `KOMEN` | Komen et al. (1984). WAM Cycle 3 heritage. |
| `JANSSEN` | Janssen (1989, 1991). Quasi-linear Miles mechanism. WAM Cycle 4. |
| `WESTH` | Van der Westhuysen et al. (2007). Yan (1987) wind input. **Default since v41.45.** |
| `ST6` | Rogers et al. (2012). Newest, observational-based. |

#### Whitecapping Dissipation

| Formulation | Description |
|---|---|
| `KOMEN` | Steepness-dependent. WAM Cycle 3. |
| `JANSSEN` | WAM Cycle 4 with different coefficients. |
| `AB` (Alves-Banner) | **Default.** Saturation-based. Distinguishes breaking vs non-breaking. Better swell behavior. |
| `ST6` | Local + cumulative dissipation terms. |

#### Bottom Friction

| Formulation | Default Coefficient | Notes |
|---|---|---|
| `JONSWAP` | C = **0.038** m²/s³ (swell over sand) | 0.067 for wind-sea in shallow water. **Default.** |
| `COLLINS` | C_f = **0.015** | Drag-law model |
| `MADSEN` | K_N = **0.05** m | Jonsson friction factor with roughness length |
| `RIPPLES` | — | Movable bed with sediment properties |

#### Depth-Induced Breaking

Battjes & Janssen (1978): `D = −α·Q_b·σ_mean·H_max² / (8π)`

| Parameter | Default | Description |
|---|---|---|
| `alpha` | 1.0 | Dissipation rate coefficient |
| `gamma` | 0.73 | Breaker index (H_max = γ·d). Tune 0.6–0.8 depending on slope. |

#### Triad Wave-Wave Interactions

| Method | Description |
|---|---|
| `DCTA` | Consistent collinear triad approximation. **Default since v41.45.** |
| `LTA` | Lumped Triad Approximation (Eldeberky 1996). Earlier default. |
| `SPB` | Becq-Girard et al. (1999). |
| `FTIM` | Full Triad Interaction Model. Newest. |

#### Quadruplet Interactions

DIA (Discrete Interaction Approximation, Hasselmann et al. 1985). Controls spectral shape in deep water. λ = 0.25, C_nl4 = 3×10⁷. Enabled by `QUADRUPL` command.

### SWAN Input File Structure

```
$ sunda_swan.swn
$
PROJECT 'Sunda' '01'
SET LEVEL 0.0
MODE NONSTATIONARY TWODIMENSIONAL
COORDINATES SPHERICAL
$
$ --- Computational grid ---
CGRID REGular -107.0 -7.5 0.0 4.0 3.0 200 150 &
      CIRCLE 36 0.03 1.0 30
$                 ↑MDC ↑flow ↑fhigh ↑MSC
$
$ --- Bottom ---
INPGRID BOTTOM REG -107.0 -7.5 0.0 200 150 0.02 0.02
READINP BOTTOM 1.0 'sunda_bathy.bot' 3 0 FREE
$
$ --- Wind ---
INPGRID WIND REG -107.0 -7.5 0.0 100 75 0.04 0.04 NONSTAT 20240901.000000 1 HR 20240906.000000
READINP WIND 1.0 SERIES 'wind_files.dat' 3 0 FREE
$
$ --- Boundary conditions ---
BOUNDSPEC SIDE WEST CCW CON FILE 'wave_boundary.sp2' 1
BOUNDSPEC SIDE SOUTH CCW CON FILE 'wave_boundary.sp2' 1
$
$ --- Initial conditions ---
INITIAL DEFAULT
$
$ --- Physics ---
GEN3 WESTH
BREAKING CONSTANT 1.0 0.73
FRICTION JONSWAP 0.038
TRIAD
QUADRUPL
$
$ --- Numerics ---
NUMERIC ACCUR 0.02 0.02 0.02 98 NONSTAT 20
$
$ --- Output ---
BLOCK 'COMPGRID' NOHEADER 'hsig.mat' LAY 3 HSIGN 1. OUTPUT 20240901.000000 1 HR
BLOCK 'COMPGRID' NOHEADER 'tps.mat' LAY 3 TPS 1. OUTPUT 20240901.000000 1 HR
BLOCK 'COMPGRID' NOHEADER 'dir.mat' LAY 3 DIR 1. OUTPUT 20240901.000000 1 HR
$
$ --- Compute ---
COMPUTE NONSTAT 20240901.000000 60 SEC 20240906.000000
STOP
```

### Key SWAN Output Variables

| Variable | Description | Units |
|---|---|---|
| `HSIGN` | Significant wave height | m |
| `TPS` | Smoothed peak period | s |
| `TM01` | Mean period (1st moment) | s |
| `TM02` | Mean period (2nd moment) | s |
| `TMM10` | Energy period (−1st moment) | s |
| `DIR` | Mean wave direction | degrees |
| `DSPR` | Directional spreading | degrees |
| `UBOT` | Near-bottom orbital velocity | m/s |
| `HSWELL` | Swell significant wave height | m |
| `WLEN` | Mean wavelength | m |
| `SETUP` | Wave-induced water level | m |
| `QB` | Fraction of breaking waves | — |
| `FORCE` | Wave-induced force (radiation stress gradient) | N/m² |
| `DISBOT` | Bottom friction dissipation | W/m² |
| `DISSURF` | Depth-induced breaking dissipation | W/m² |
| `DISWCAP` | Whitecapping dissipation | W/m² |

### Running SWAN (Standalone)

```bash
# Compile
make config
make ser            # Serial
make omp            # OpenMP
make mpi            # MPI

# Run
./swanrun -input myrun                    # Serial
./swanrun -input myrun -omp 4            # OpenMP
./swanrun -input myrun -mpi 8            # MPI
```

---

## 3. Coupling Mechanisms & Field Exchanges

Models exchange physical fields via MCT at user-specified intervals. **SCRIP** interpolation weights handle remapping between different grids.

### WRF → ROMS (Atmosphere to Ocean)

| Field | Description | Units |
|---|---|---|
| GSW | Net shortwave radiation | W/m² |
| GLW | Net longwave radiation | W/m² |
| LH | Latent heat flux | W/m² |
| HFX | Sensible heat flux | W/m² |
| USTRESS / VSTRESS | Surface wind stress | Pa |
| MSLP | Mean sea level pressure | Pa |
| T2 | 2-m air temperature | °C |
| U10 / V10 | 10-m wind components | m/s |
| RELH | Relative humidity | % |
| RAIN | Precipitation rate | m/s |
| CLDFRA | Cloud fraction | 0–1 |

**Flux computation** (choose one):
- `ATM2OCN_FLUXES`: WRF computes fluxes, passes them to ROMS
- `BULK_FLUXES`: ROMS receives state variables, computes fluxes internally

### ROMS → WRF (Ocean to Atmosphere)

| Field | Description |
|---|---|
| SST | Sea surface temperature (dynamically evolving) |

### ROMS → SWAN (Ocean to Wave)

| Field | Description |
|---|---|
| DEPTH | Bathymetry |
| WLEV | Free-surface elevation |
| VELX / VELY | Depth-integrated currents |
| ZO | Bottom roughness |

### SWAN → ROMS (Wave to Ocean)

| Field | Description |
|---|---|
| HSIGN | Significant wave height |
| RTP | Peak wave period |
| TMBOT | Bottom mean wave period |
| UBOT | Bottom orbital velocity |
| DIRE / DIRN | Mean wave direction components |
| LWAVE / LWAVEP | Mean/peak wavelength |
| DISBOT | Bottom friction dissipation |
| DISSURF | Breaking dissipation |
| DISWCAP | Whitecapping dissipation |
| QB | Percent wave breaking |

### WRF → SWAN

| Field | Description |
|---|---|
| U10 / V10 | 10-m wind components |

### SWAN → WRF

| Field | Description |
|---|---|
| HSIGN | Significant wave height |
| LWAVEP | Peak wavelength |
| Charnock coefficient | Wave-age-dependent sea surface roughness |

> **Key physics:** The wave-dependent Charnock coefficient is one of the most important coupling pathways. Young, steep waves under strong winds produce higher roughness than mature swell — critical for hurricane/storm simulations.

---

## 4. Compilation

### Dependencies

| Library | Notes |
|---|---|
| Fortran compiler | ifort, gfortran (v10+: add `-fallow-argument-mismatch`) |
| MPI | OpenMPI, Intel MPI, or MPICH |
| NetCDF-4 + Fortran | Required for all models |
| HDF5 | Required by NetCDF-4 |
| MCT | Included in `Lib/MCT/`, compiled as part of build |

### Step 1: Compile MCT (one-time)

```bash
cd Lib/MCT
./configure
make
cd ../..
```

### Step 2: Configure `build_coawst.sh`

Key variables to set:

```bash
COAWST_APPLICATION=SUNDA_STRAIT    # Must match header file name
ROMS_APPLICATION=${COAWST_APPLICATION}

MY_ROOT_DIR=/path/to/COAWST
MY_HEADER_DIR=${MY_ROOT_DIR}/Projects/Sunda_Strait
MY_ANALYTICAL_DIR=${MY_ROMS_SRC}/ROMS/Functionals
SCRATCH_DIR=./Build_roms

USE_MPI=on
USE_MPIF90=on
FORT=ifort                          # or gfortran
USE_LARGE=on
USE_NETCDF4=on
```

### Step 3: Build

```bash
./build_coawst.sh -j 8
```

**Build sequence:** Clean → MCT params → WRF → SWAN → ROMS/COAWST → Link → `coawstM`

**Build options:**

| Flag | Description |
|---|---|
| `-j N` | Parallel compilation with N threads |
| `-noclean` | Skip cleaning all components |
| `-nocleanwrf` | Skip cleaning WRF only |
| `-nocleanswan` | Skip cleaning SWAN only |

### Application Header File

```c
/* sunda_strait.h */

/* Models to couple */
#define ROMS_MODEL
#define SWAN_MODEL
#define WRF_MODEL
#define MCT_LIB

/* Grid interpolation */
#define MCT_INTERP_OC2AT      /* Ocean ↔ Atmosphere */
#define MCT_INTERP_OC2WV      /* Ocean ↔ Wave */
#define MCT_INTERP_WV2AT      /* Wave ↔ Atmosphere */

/* Flux approach */
#define ATM2OCN_FLUXES        /* WRF computes fluxes */

/* Wave-current interaction */
#define WEC_VF                /* Vortex force formalism */
#define WDISS_WAVEMOD         /* Wave dissipation from SWAN */
#define UV_KIRBY              /* Kirby-Chen mean current for Doppler */

/* ROMS dynamics */
#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPHERICAL
#define MASKING
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define RI_SPLINES
#define DJ_GRADPS             /* Spline density Jacobian pressure gradient */

/* Tracers */
#define TS_DIF2
#define MIX_GEO_TS

/* Vertical mixing */
#define GLS_MIXING
#define CANUTO_A
#define N2S2_HORAVG

/* Tides */
#define SSH_TIDES
#define UV_TIDES
#define TIDE_GENERATING_FORCES
#define RAMP_TIDES

/* Bottom */
#define UV_QDRAG

/* Output */
#define AVERAGES
#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV

/* Analytical */
#define ANA_BTFLUX
#define ANA_BSFLUX
```

---

## 5. SCRIP Weight Generation

SCRIP computes interpolation weights between model grids. Required whenever models use different grids (which is almost always).

### Configure SCRIP Input

```
$ scrip_coawst_sunda.in
$INPUTS
  OUTPUT_NCFILE = 'scrip_sunda_weights.nc'

  ROMS_GRIDS  = 1
  SWAN_GRIDS  = 1
  WW3_GRIDS   = 0
  WRF_GRIDS   = 1
  HYD_GRIDS   = 0

  ROMS_GRIDFILE1 = 'Projects/Sunda_Strait/roms_grid.nc'

  SWAN_COORD1   = 'Projects/Sunda_Strait/swan_coord.grd'
  SWAN_BATH1    = 'Projects/Sunda_Strait/swan_bathy.bot'
  SWAN_NUMX1    = 200
  SWAN_NUMY1    = 150
  CARTESIAN1    = 0             ! 0=Spherical

  WRF_GRIDS1    = 'Projects/Sunda_Strait/wrfinput_d01'
$END
```

### Build and Run SCRIP

```bash
cd Lib/SCRIP_COAWST
make clean && make
./scrip_coawst < scrip_coawst_sunda.in
# Produces: scrip_sunda_weights.nc
```

---

## 6. Configuration Files

### 6.1 Coupling Input File (`coupling_sunda.in`)

The **single argument** passed to `coawstM` at runtime.

```
! Processor allocation (total must match mpirun -np)
  NnodesATM =  64              ! WRF MPI ranks
  NnodesWAV =   8              ! SWAN MPI ranks
  NnodesOCN =  48              ! ROMS MPI ranks
  NnodesHYD =   0              ! WRF-Hydro (disabled)

! Coupling exchange intervals (seconds)
  TI_ATM2WAV =  600.0          ! WRF → SWAN every 10 min
  TI_ATM2OCN =  600.0          ! WRF → ROMS every 10 min
  TI_WAV2ATM =  600.0          ! SWAN → WRF every 10 min
  TI_WAV2OCN =  600.0          ! SWAN → ROMS every 10 min
  TI_OCN2WAV =  600.0          ! ROMS → SWAN every 10 min
  TI_OCN2ATM =  600.0          ! ROMS → WRF every 10 min

! Input files for each model
  ATM_name = namelist.input
  WAV_name = Projects/Sunda_Strait/swan_sunda.in
  OCN_name = Projects/Sunda_Strait/ocean_sunda.in

! SCRIP weights
  SCRIP_WEIGHT_OPTION = 1
  SCRIP_COAWST_NAME = Projects/Sunda_Strait/scrip_sunda_weights.nc
```

### 6.2 Individual Model Configs

Each component uses its own standard configuration:

| Model | Config File | Refer To |
|---|---|---|
| **WRF** | `namelist.input` | WRF Model Guide |
| **ROMS** | `ocean_sunda.in` | ROMS Model Guide |
| **SWAN** | `swan_sunda.in` | Section 2 above |

> **Important:** `NtileI × NtileJ` in ROMS `ocean.in` must equal `NnodesOCN` in the coupling file.

---

## 7. Running COAWST

### MPI Domain Decomposition

Processors are partitioned **sequentially** by MPI rank:

```
Example: 120 total MPI processes
  WRF:  ranks 0–63   (NnodesATM = 64)
  SWAN: ranks 64–71  (NnodesWAV = 8)
  ROMS: ranks 72–119 (NnodesOCN = 48)
```

**Critical:** `NnodesATM + NnodesWAV + NnodesOCN + NnodesHYD` must **exactly** equal `mpirun -np N`.

### Execution

```bash
# Run
mpirun -np 120 ./coawstM Projects/Sunda_Strait/coupling_sunda.in > coawst.log

# SLURM example
#!/bin/bash
#SBATCH --job-name=coawst_sunda
#SBATCH --ntasks=120
#SBATCH --time=48:00:00
#SBATCH --output=coawst_%j.out

mpirun ./coawstM Projects/Sunda_Strait/coupling_sunda.in > coawst.log
```

### Execution Flow

1. MPI initializes; `coawstM` reads `coupling_*.in`
2. Processor groups assigned by model "color" via `MPI_COMM_SPLIT`
3. MCTWorld initialized; routers established between model pairs
4. Each model group initializes its component
5. SCRIP weights loaded; interpolation matrices built
6. Time-stepping loop with coupling exchanges at specified intervals
7. Finalization: MCT cleanup + `MPI_FINALIZE`

### Processor Allocation Tips

| Model | Typical Share | Notes |
|---|---|---|
| WRF | 50–60% | Most expensive component |
| ROMS | 30–40% | Second most expensive |
| SWAN | 5–10% | Usually least expensive |

---

## 8. Sunda Strait Application: Building Your Own Ina-CAWO

### What Ina-CAWO Does

| Feature | Ina-CAWO | Your Setup |
|---|---|---|
| Domain | All Indonesian waters (15°S–15°N, 90°E–145°E) | Sunda Strait focus (~4°×3° region) |
| Resolution | 3 km | 1 km (ROMS/SWAN) + 3 km (WRF parent) + 1 km (WRF nested) |
| Components | WRF + ROMS + SWAN | Same |
| Coupling | COAWST/MCT | Same |
| Forecast | 10 days, 4 cycles/day | Configurable |

### Step-by-Step Setup

#### Step 1: WRF Domain Design

Two nested WRF domains:
- **d01** (parent): 3 km, covers wider Indonesian maritime region
- **d02** (nest): 1 km, covers Sunda Strait area

Key `namelist.input` physics for tropical maritime:
```
&physics
 mp_physics        = 8, 8,         ! Thompson microphysics
 cu_physics        = 0, 0,         ! No cumulus at 3 km / 1 km
 ra_lw_physics     = 4, 4,         ! RRTMG longwave
 ra_sw_physics     = 4, 4,         ! RRTMG shortwave
 bl_pbl_physics    = 5, 5,         ! MYNN PBL
 sf_sfclay_physics = 5, 5,         ! MYNN surface layer
 sf_surface_physics= 2, 2,         ! Noah LSM
 sst_update        = 1,            ! SST from ROMS coupling
/
```

#### Step 2: ROMS Grid (Sunda Strait)

- Resolution: ~1 km
- Domain: ~104.5°E–106.5°E, ~7.5°S–5.5°S
- Bathymetry: GEBCO 2025, preserve 80 m sill depth
- Smoothing: rx0 < 0.3, rx1 < 7
- Vertical: Vtransform=2, Vstretching=4, N=40, theta_s=6, theta_b=2, hc=20

Key ROMS physics:
```c
#define GLS_MIXING            /* k-epsilon for tidal mixing */
#define DJ_GRADPS             /* Essential for steep bathymetry */
#define SSH_TIDES             /* 8+ tidal constituents */
#define UV_TIDES
#define TIDE_GENERATING_FORCES
#define RAMP_TIDES
#define WEC_VF                /* Wave-current interaction */
#define WDISS_WAVEMOD
```

#### Step 3: SWAN Grid (Sunda Strait)

- Same domain as ROMS or slightly larger
- Resolution: ~1 km (can differ from ROMS; SCRIP handles interpolation)
- 36 directional bins, 30 frequency bins (0.03–1.0 Hz)
- Physics: GEN3 WESTH, JONSWAP bottom friction (0.038), Battjes-Janssen breaking

#### Step 4: Forcing Data

| Data | Source | For |
|---|---|---|
| Atmospheric IC/BC | GFS 0.25° or ERA5 | WRF (via WPS) |
| Ocean IC/BC | GLORYS12V1 (1/12°) | ROMS |
| Tidal forcing | TPXO10 or FES2022 | ROMS (8 constituents: M2, S2, N2, K2, K1, O1, P1, Q1) |
| Wave boundaries | ERA5 wave or WAVEWATCH III | SWAN (spectral boundary) |

#### Step 5: SCRIP Weights

Generate interpolation weights between WRF d02, ROMS, and SWAN grids.

#### Step 6: Coupling Configuration

```
  NnodesATM = 64               ! WRF (both domains)
  NnodesWAV =  8               ! SWAN
  NnodesOCN = 48               ! ROMS
  NnodesHYD =  0

  TI_ATM2WAV = 600.0           ! 10-minute exchange
  TI_ATM2OCN = 600.0
  TI_WAV2ATM = 600.0
  TI_WAV2OCN = 600.0
  TI_OCN2WAV = 600.0
  TI_OCN2ATM = 600.0
```

### Why Coupling Matters for Sunda Strait

| Process | Without Coupling | With Coupling |
|---|---|---|
| **Sea surface roughness** | Constant Charnock (0.0185) | Wave-age-dependent roughness from SWAN → WRF |
| **SST** | Fixed or slowly updated | Dynamic SST from ROMS → WRF (upwelling, tidal mixing) |
| **Wave-current interaction** | Ignored | Tidal currents modify wave propagation (ROMS → SWAN) |
| **Wave setup** | Ignored | Wave radiation stress modifies water level (SWAN → ROMS) |
| **Bottom orbital velocity** | Ignored | SWAN bottom velocity drives sediment resuspension in ROMS |
| **Wind-wave generation** | May use different wind | Consistent wind fields from WRF to both ROMS and SWAN |

---

## 9. Troubleshooting & Tips

### Compilation Issues

| Problem | Solution |
|---|---|
| gfortran 10+ type mismatch errors | Add `-fallow-argument-mismatch -fallow-invalid-boz` to WRF `configure.wrf` |
| Missing `module_wrf_top.o` | WRF didn't compile fully. Check WRF logs. Don't use `-nocleanwrf` on first build. |
| Missing `ocn_comm_world` | Verify both `ROMS_MODEL` and `MCT_LIB` defined in header |
| `wec_vf.f90` error | Add `#define WDISS_WAVEMOD` when using `#define WEC_VF` |

### Runtime Issues

| Problem | Solution |
|---|---|
| Processor count mismatch | Sum of `Nnodes*` must exactly equal `mpirun -np N` |
| "Unable to open input file" | Paths in `coupling_*.in` are relative to launch directory. Use absolute paths. |
| CFL blowup | Reduce time steps. In coupled mode, SST feedback can alter WRF dynamics → reduce WRF dt too. |
| `MPI_Send` fatal error | Check coupling intervals are consistent. Ensure no model waits for data that won't arrive. |

### Coupling Issues

| Problem | Solution |
|---|---|
| Model doesn't receive data | Check: (1) `TI_*` interval is non-zero, (2) `Nnodes* > 0`, (3) `MCT_INTERP_*` CPP flags defined, (4) SCRIP weights file exists |
| Interpolation artifacts / NaN | Regenerate SCRIP weights with correct grid files and masks |
| SCRIP crashes | Verify grid dimensions match files. Check `CARTESIAN` flag (0=spherical, 1=Cartesian). |

### Performance Tuning

- **Coupling interval:** 600–900 s for mesoscale. Decrease for rapidly evolving events (storms).
- **Processor balance:** Profile each model's wall time separately, then allocate proportionally.
- **I/O:** Enable `USE_PARALLEL_IO` + `USE_HDF5` for large runs.
- **Minimum depth:** Set ROMS minimum bathymetry ≥ 5–10 m to avoid shallow-water CFL issues.

---

## 10. SWAN Grid Pre-Processing

### Grid File Formats

SWAN needs two files for curvilinear grid input: a **coordinate file** (`.grd`) and a **bathymetry file** (`.bot`). For **regular grids** (`CGRID REGular`), no `.grd` file is needed — the grid is defined entirely by the `CGRID` command parameters.

#### Coordinate File (`.grd`)

Plain ASCII with **one value per line** in **Fortran column-major order**. Two consecutive blocks: all x-coordinates (longitudes) first, then all y-coordinates (latitudes).

For a grid with NUMX × NUMY points, the file has `2 × NUMX × NUMY` lines total:

```
 104.500000000000     ← lon(1,1)
 104.510000000000     ← lon(2,1)
 104.520000000000     ← lon(3,1)
 ...                  ← (all longitudes, column-major)
 106.500000000000     ← lon(NUMX,NUMY)
  -7.200000000000     ← lat(1,1)  — second half begins
  -7.200000000000     ← lat(2,1)
 ...
  -5.500000000000     ← lat(NUMX,NUMY)
```

**Column-major order** means the x-index varies fastest: for grid point `(i,j)`, the longitude is at line `NUMX*(j-1) + i`, and the latitude is at line `NUMX*NUMY + NUMX*(j-1) + i`.

#### Bathymetry File (`.bot`)

Plain ASCII, NUMY rows with NUMX space-separated depth values each. Depths are **positive downward**. Land points use the exception value **9999**.

```
   9999.0000   9999.0000     80.2000     79.5000     78.1000  ...
   9999.0000     82.1000     81.3000     80.0000     78.5000  ...
   ...
```

#### The `idla` Parameter

The `idla` parameter in `READGRID`/`READINP` controls spatial ordering. COAWST universally uses **`idla=4`** (data from lower-left corner, line breaks not required), matching ROMS Fortran array ordering.

| idla | Description |
|---|---|
| 1 | From upper-left, row by row, line breaks required |
| 3 | From lower-left, row by row, line breaks required |
| 4 | From lower-left, free format (line breaks not required) — **COAWST standard** |

### Dimension Gotcha: CGRID vs NUMX

This is a critical source of confusion:

- **`CGRID CURVILINEAR mxc myc`**: `mxc` and `myc` are the number of **meshes** (cells), not grid points
- **Grid points** = meshes + 1: `NUMX = mxc + 1`, `NUMY = myc + 1`
- **`scrip_coawst.in`** requires the full number of **grid points** (NUMX, NUMY)

Example: ROMS grid with dimensions `(eta_rho=171, xi_rho=201)`:
- SWAN NUMX = 201, NUMY = 171
- CGRID CURVILINEAR **200 170** (meshes = points − 1)
- scrip_coawst.in: `SWAN_NUMX=201, SWAN_NUMY=171` (grid points)

### COAWST MATLAB Tool: `roms2swan.m`

The primary tool is `Tools/mfiles/mtools/roms2swan.m`. It reads a ROMS grid and writes both SWAN files:

```matlab
% Usage — just one command:
roms2swan('sunda_roms_grid.nc')
% Creates: swan_coord.grd and swan_bathy.bot
```

What it does internally:
1. Reads `lon_rho`, `lat_rho`, `h`, `mask_rho` from the ROMS grid
2. Sets land points (`mask_rho == 0`) to 9999 in bathymetry
3. Writes coordinates with `fprintf(fid,'%18.12f\n', x_rho)` — MATLAB's `fprintf` automatically outputs arrays in Fortran column-major order
4. Writes bathymetry row by row

Other COAWST MATLAB tools:

| Script | Location | Purpose |
|---|---|---|
| `roms2swan.m` | `Tools/mfiles/mtools/` | Create `.grd` and `.bot` from ROMS grid |
| `create_swan_2dspec_from_ww3.m` | `Tools/mfiles/mtools/` | SWAN spectral boundary from WW3 |
| `create_swanTpar_from_WW3.m` | `Tools/mfiles/swan_forc/` | TPAR boundary from WW3 |

### Python Equivalent: `roms2swan`

```python
import numpy as np
import netCDF4 as nc

def roms2swan(roms_grid_file,
              coord_file='swan_coord.grd',
              bath_file='swan_bathy.bot'):
    """
    Python equivalent of COAWST's roms2swan.m.
    Creates SWAN .grd and .bot files from a ROMS grid.
    """
    ds = nc.Dataset(roms_grid_file)
    lon_rho = ds.variables['lon_rho'][:]   # shape: (eta_rho, xi_rho)
    lat_rho = ds.variables['lat_rho'][:]
    h = ds.variables['h'][:]
    mask_rho = ds.variables['mask_rho'][:]
    ds.close()

    # Mask land with 9999
    h_swan = h.copy()
    h_swan[mask_rho == 0] = 9999.0

    ny, nx = lon_rho.shape
    print(f"Grid points: NUMX={nx}, NUMY={ny}")
    print(f"CGRID CURVILINEAR {nx-1} {ny-1}")

    # Write coordinate file — one value per line, Fortran column-major order
    # .T then flatten gives column-major ordering
    with open(coord_file, 'w') as f:
        for val in lon_rho.T.flatten():
            f.write(f'{val:18.12f}\n')
        for val in lat_rho.T.flatten():
            f.write(f'{val:18.12f}\n')

    # Write bathymetry file — NUMY rows, NUMX values each
    with open(bath_file, 'w') as f:
        for j in range(ny):
            line = '   '.join(f'{h_swan[j, i]:12.4f}' for i in range(nx))
            f.write(f'   {line}\n')

    print(f"Created {coord_file} and {bath_file}")
    return nx, ny

# Usage:
numx, numy = roms2swan('sunda_roms_grid.nc')
```

### Creating a Regular Grid from GEBCO (Python)

If you don't have a ROMS grid yet, create a standalone SWAN regular grid:

```python
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import netCDF4 as nc

# Define SWAN domain
lon_min, lon_max = 104.5, 106.5
lat_min, lat_max = -7.5, -5.5
dx, dy = 0.01, 0.01  # ~1 km

lons = np.arange(lon_min, lon_max + dx/2, dx)
lats = np.arange(lat_min, lat_max + dy/2, dy)
numx, numy = len(lons), len(lats)

# Load and interpolate GEBCO
gebco = nc.Dataset('GEBCO_2024.nc')
g_lon = gebco.variables['lon'][:]
g_lat = gebco.variables['lat'][:]
mask_lon = (g_lon >= lon_min - 0.5) & (g_lon <= lon_max + 0.5)
mask_lat = (g_lat >= lat_min - 0.5) & (g_lat <= lat_max + 0.5)
elev = gebco.variables['elevation'][np.ix_(mask_lat, mask_lon)]
gebco.close()

interp = RegularGridInterpolator(
    (g_lat[mask_lat], g_lon[mask_lon]), elev, method='linear')
lon_grid, lat_grid = np.meshgrid(lons, lats)
depth = -interp(np.column_stack([lat_grid.ravel(),
         lon_grid.ravel()])).reshape(lat_grid.shape)

# SWAN conventions: positive down, land = 9999
depth[depth <= 0] = 9999.0
depth = np.maximum(depth, 5.0)  # minimum depth for stability

np.savetxt('swan_bathy.bot', depth, fmt='%10.2f')
print(f"Regular grid: {numx} x {numy}")
```

Corresponding SWAN input (no `.grd` file needed for regular grids):

```
CGRID REGular 104.5 -7.5 0.0 2.0 2.0 200 200 &
      CIRCLE 36 0.03 1.0 30

INPGRID BOTTOM REGular 104.5 -7.5 0.0 200 200 0.01 0.01
READINP BOTTOM 1.0 'swan_bathy.bot' 4 0 FREE
```

### Corresponding SWAN Input for Curvilinear Grid

```
CGRID CURVilinear 200 170 EXC 9.999000e+003 &
      CIRCLE 36 0.04 1.0 24
READGRID COORDINATES 1 'swan_coord.grd' 4 0 0 FREE

INPGRID BOTTOM CURVilinear 0 0 200 170 EXC 9.999000e+003
READINP BOTTOM 1 'swan_bathy.bot' 4 0 FREE
```

### Grid Alignment Tips

- **Same grid as ROMS** (recommended): Use `roms2swan.m` or the Python equivalent. SCRIP weights become trivial (near-identity mapping). Simplifies debugging.
- **Slightly larger than ROMS**: SWAN domain can extend beyond ROMS to avoid wave energy losses at boundaries. SCRIP handles interpolation.
- **Different resolution**: SWAN and ROMS grids can differ in resolution. SWAN often uses coarser resolution since wave physics doesn't require the same resolution as tidal/current dynamics.
- **Spectral resolution**: 36 directional bins and 24–30 frequency bins (0.04–1.0 Hz) is standard. Increase directional bins for narrow straits.

---

## 11. Output Files & Post-Processing

### Output Files by Component

#### ROMS Output

| File | Config Parameter | Description |
|---|---|---|
| `ocean_his.nc` | `NHIS` | Snapshots at specified intervals. Full 3D fields. |
| `ocean_avg.nc` | `NAVG` | Time-averaged fields. Removes tidal signals. |
| `ocean_rst.nc` | `NRST` | Restart file. Ping-pong (2 records). |
| `ocean_dia.nc` | `NDIA` | Momentum/tracer equation term diagnostics. Requires `DIAGNOSTICS_UV`/`DIAGNOSTICS_TS`. |
| `ocean_sta.nc` | `NSTA` | Time series at specified station locations. |
| `ocean_flt.nc` | — | Lagrangian float trajectories. Requires `FLOATS`. |

Key ROMS coupled variables in `ocean_his.nc`:

| Variable | Description | Indicates Coupling With |
|---|---|---|
| `Hwave` | Significant wave height | SWAN → ROMS |
| `Dwave` | Mean wave direction | SWAN → ROMS |
| `Lwave` | Mean wavelength | SWAN → ROMS |
| `Pwave_bot` | Bottom wave period | SWAN → ROMS |
| `sustr`, `svstr` | Surface wind stress | WRF → ROMS (if `ATM2OCN_FLUXES`) |
| `shflux` | Net surface heat flux | WRF → ROMS |

#### WRF Output

| File | Config Parameter | Description |
|---|---|---|
| `wrfout_d0X_*` | `history_interval` (min) | Main output. All atmospheric variables. |
| `wrfrst_d0X_*` | `restart_interval` (min) | Restart file. |

Key WRF coupled variable: `SST` (should change over time if ROMS → WRF coupling active).

#### SWAN Output

| File | Format | Description |
|---|---|---|
| `swan_*_his.nc` | NetCDF | Wave field snapshots (Hsig, Dir, TPS, etc.) |
| `swan_rst.dat` | ASCII | Hotstart/restart (2D spectral energy at all grid points) |

### Python Post-Processing

```bash
pip install xarray netCDF4 matplotlib cartopy cmocean dask
```

#### Surface Temperature Map

```python
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean

ds = xr.open_dataset('ocean_his.nc')
sst = ds['temp'].isel(ocean_time=-1, s_rho=-1)

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(12, 8))
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.COASTLINE)
pc = ax.pcolormesh(ds['lon_rho'], ds['lat_rho'], sst,
                    transform=ccrs.PlateCarree(),
                    cmap=cmocean.cm.thermal, shading='auto')
plt.colorbar(pc, ax=ax, label='Temperature (°C)')
ax.set_title('ROMS Sea Surface Temperature')
ax.gridlines(draw_labels=True)
plt.savefig('sst_map.png', dpi=150, bbox_inches='tight')
```

#### Vertical Section

```python
import numpy as np

ds = xr.open_dataset('ocean_his.nc')
eta_idx = 50  # cross-section index

temp = ds['temp'].isel(ocean_time=-1, eta_rho=eta_idx)
h = ds['h'].isel(eta_rho=eta_idx)
zeta = ds['zeta'].isel(ocean_time=-1, eta_rho=eta_idx)
Cs_r = ds['Cs_r']
hc = float(ds['hc'])

# Compute depth (Vtransform=2)
S = (hc * Cs_r + h * Cs_r) / (hc + h)
z = zeta + (zeta + h) * S

lon = ds['lon_rho'].isel(eta_rho=eta_idx)

fig, ax = plt.subplots(figsize=(14, 6))
pc = ax.pcolormesh(
    lon.values[np.newaxis, :] * np.ones_like(z.values),
    z.values, temp.values,
    cmap=cmocean.cm.thermal, shading='auto')
plt.colorbar(pc, ax=ax, label='Temperature (°C)')
ax.set_xlabel('Longitude')
ax.set_ylabel('Depth (m)')
ax.set_title('Vertical Temperature Section')
plt.savefig('vertical_section.png', dpi=150, bbox_inches='tight')
```

#### Compare Coupled Fields

```python
import numpy as np

# WRF wind speed
wrf = xr.open_dataset('wrfout_d01_2024-09-01_00:00:00')
wspd = np.sqrt(wrf['U10']**2 + wrf['V10']**2).isel(Time=-1)

# ROMS wave height (from SWAN coupling)
roms = xr.open_dataset('ocean_his.nc')
hwave = roms['Hwave'].isel(ocean_time=-1)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
wspd.plot(ax=axes[0], cmap='viridis')
axes[0].set_title('WRF 10m Wind Speed (m/s)')
hwave.plot(ax=axes[1], cmap=cmocean.cm.amp)
axes[1].set_title('ROMS Wave Height from SWAN (m)')
plt.tight_layout()
plt.savefig('coupled_fields.png', dpi=150)
```

---

## 12. Verifying Coupling Works

After your first coupled run, **always verify** that models are actually exchanging data. A run that completes without errors does not guarantee coupling is active.

### Step 1: Check Log Output

COAWST prints MCT initialization messages at startup. Look for:
- `"MCT: Initialize coupler..."` for each model pair
- Periodic exchange messages at your specified `TI_*` interval
- No warnings about zero-valued fields or failed interpolation

### Step 2: SST Exchange (WRF ↔ ROMS)

The easiest coupling check. If ROMS sends SST to WRF, the `SST` field in `wrfout` should **change over time** and match the ROMS surface temperature.

```python
import xarray as xr

wrf = xr.open_dataset('wrfout_d01_2024-09-01_00:00:00')
sst_wrf = wrf['SST']

# If SST is constant in time → coupling NOT working
print("SST std over time:", float(sst_wrf.std('Time')))
# Should be > 0 for active coupling

roms = xr.open_dataset('ocean_his.nc')
sst_roms = roms['temp'].isel(s_rho=-1)
print("ROMS SST range:", float(sst_roms.min()), "to", float(sst_roms.max()))
```

### Step 3: Wave Fields in ROMS (SWAN → ROMS)

If coupling works, ROMS `ocean_his.nc` should contain wave variables with non-zero, physically reasonable values.

```python
roms = xr.open_dataset('ocean_his.nc')

for var in ['Hwave', 'Dwave', 'Lwave', 'Pwave_bot']:
    if var in roms:
        print(f"{var}: min={float(roms[var].min()):.3f}, "
              f"max={float(roms[var].max()):.3f}")
    else:
        print(f"WARNING: {var} not found — SWAN→ROMS coupling inactive")
```

Expected values:
- `Hwave`: 0.1–5 m (coastal), should vary spatially
- `Lwave`: 10–300 m
- `Pwave_bot`: 3–15 s

### Step 4: Coupling Interval Consistency

The coupling interval must divide evenly into each model's time step:

```
TI_ATM2OCN = 600 s, WRF dt = 20 s  → 600/20 = 30 ✓
TI_WAV2OCN = 600 s, ROMS dt = 60 s → 600/60 = 10 ✓
TI_WAV2OCN = 600 s, SWAN dt = 3 s  → 600/3 = 200 ✓
```

If the division is not exact, coupling exchanges will be skipped or misaligned.

### Step 5: Coastal Masking

A known issue: WRF can send zero-valued forcing at coastal strips due to land/sea mask mismatch. Check wind stress near coastlines in ROMS output for suspicious zero bands.

### Quick Checklist

| Check | How | Pass Criteria |
|---|---|---|
| SST varies in `wrfout` | Plot SST at multiple times | SST changes between timesteps |
| `Hwave` exists in ROMS | `ncdump -h ocean_his.nc` | Variable present and non-zero |
| Wind stress non-zero | Plot `sustr`/`svstr` | No zero bands at coast |
| Log shows exchanges | `grep "MCT" coawst.log` | Exchange messages at `TI_*` intervals |
| Energy conservation | Compare domain-total KE over time | No sudden jumps at exchange times |

---

## 13. Spin-Up Strategies

### How Long to Spin Up

| Component | Typical Spin-Up | Monitor |
|---|---|---|
| **ROMS (barotropic/tidal)** | 1–7 days | Sea surface height, tidal amplitudes reaching equilibrium |
| **ROMS (baroclinic)** | Weeks to months | Domain-averaged kinetic energy, temperature/salinity drift |
| **WRF** | 6–12 hours | Precipitation patterns, boundary layer development |
| **SWAN** | Hours to ~1 day | Wave energy reaching steady state (depends on fetch) |

**Key diagnostic:** Monitor domain-averaged **kinetic energy** over time. When KE stabilizes (oscillates around a quasi-steady value), spin-up is complete.

### Cold Start vs Warm Start

| | Cold Start | Warm Start |
|---|---|---|
| **ROMS** | Initialize from GLORYS12V1 or climatology | Use restart file from previous run |
| **WRF** | Initialize from GFS/ERA5 analysis | Use `wrfrst_d0X_*` restart file |
| **SWAN** | Flat sea (no wave energy) | Hotstart file (`swan_rst.dat`) |
| **Pro** | Simple, no prior runs needed | Short spin-up, immediate realism |
| **Con** | Long spin-up, first days unrealistic | Requires managing restart files |

### Recommended Approach: Sequential Spin-Up

For first-time users, spin up components separately before coupling:

```
Phase 1 (days -14 to -7):  Run ROMS uncoupled → develop ocean circulation
Phase 2 (days -7 to -1):   Run SWAN standalone → develop wave field
Phase 3 (day 0 onward):    Start coupled COAWST using ROMS restart +
                            SWAN hotstart as initial conditions
```

This avoids the shock of sudden coupling feedback during adjustment.

### Alternative: Simultaneous Start with Ramp

Start all models together but enable `RAMP_TIDES` to gradually introduce tidal forcing:

```c
#define RAMP_TIDES    /* ramp tides from zero over ~1 day */
```

Allow at least 1–2 days after the ramp completes before using the output as "production" data.

---

## 14. Restart & Cycling

### Restart Files by Component

| Component | Restart File | Config Parameter |
|---|---|---|
| ROMS | `ocean_rst.nc` | `NRST` (timesteps between writes) |
| WRF | `wrfrst_d0X_YYYY-MM-DD_HH:MM:SS` | `restart_interval` (minutes) |
| SWAN | `swan_rst.dat` | `RESTART` command in SWAN input |

**Critical:** Set all three restart intervals to the **same wall-clock time** so you have a consistent set of restart files.

### ROMS Restart Configuration

In `ocean.in`:

```
! Write restart every 6 hours (if DT=2.5s: 6*3600/2.5 = 8640)
NRST == 8640

! Cold start:
NRREC == 0
ININAME == ocean_ini.nc

! For restart:
! NRREC == -1              ! read last record automatically
! ININAME == ocean_rst.nc  ! use restart file as initial condition

! Create new output files on restart (T) or append (F)
LDEFOUT == T
```

### WRF Restart Configuration

In `namelist.input`:

```
&time_control
 restart                  = .false.      ! .true. for restart
 restart_interval         = 360          ! minutes between restart writes
 override_restart_timers  = .true.       ! recommended for coupled restarts
/
```

### SWAN Restart Configuration

In the SWAN input file:

```
$ Cold start (omit INIT HOTSTART):
INIT DEFAULT

$ Write restart every 6 hours:
RESTART 'swan_rst.dat' FREE 6 HR

$ For warm start (add INIT HOTSTART):
$ INIT HOTSTART 'swan_rst.dat'
```

### Coupled Restart Procedure

1. **Verify all restart files exist at the same time:**
   - `ocean_rst.nc`
   - `wrfrst_d01_YYYY-MM-DD_HH:MM:SS`
   - `swan_rst.dat`

2. **Modify configuration files:**

   ROMS (`ocean.in`):
   ```
   NRREC == -1
   ININAME == ocean_rst.nc
   LDEFOUT == T
   ```

   WRF (`namelist.input`):
   ```
   restart = .true.
   ```

   SWAN input:
   ```
   INIT HOTSTART 'swan_rst.dat'
   ```

3. **Run as normal:**
   ```bash
   mpirun -np 120 ./coawstM coupling_sunda.in > coawst_restart.log
   ```

### Operational Cycling (Daily Forecast)

```
Day 1 (initialization):
  1. Run ROMS uncoupled for 7–14 days spin-up
  2. Run SWAN standalone to develop wave field

Day 2+ (daily cycle):
  1. Download latest GFS atmospheric forcing
  2. Download latest ocean BC (GLORYS near-real-time, if available)
  3. Use previous cycle's restart files:
     ROMS:  NRREC=-1, ININAME=ocean_rst.nc
     WRF:   restart=.true.
     SWAN:  INIT HOTSTART 'swan_rst.dat'
  4. Run coupled COAWST for 24–48 hours
  5. Save restart files; archive output
  6. Repeat next day
```

---

## 15. Sediment Transport (CSTMS)

COAWST integrates the **Community Sediment Transport Modeling System (CSTMS)** into ROMS for suspended load, bedload, and morphodynamic evolution.

### Enabling Sediment Transport

Minimum CPP flags in your application header:

```c
/* Core sediment */
#define SEDIMENT              /* master switch */
#define SUSPLOAD              /* suspended load transport */

/* Bottom boundary layer (choose one) */
#define SSW_BBL               /* Sherwood-Signell-Warner (recommended with waves) */
/* #define MB_BBL */          /* Meinte Blaas */
/* #define SG_BBL */          /* Styles-Glenn */

/* Bedload transport (choose one or more) */
#define BEDLOAD_SOULSBY       /* waves + currents (continental shelf) */
/* #define BEDLOAD_MPM */     /* currents only (rivers, estuaries) */
/* #define BEDLOAD_VANDERA */ /* nearshore, wave nonlinearity */

/* Morphodynamics */
#define SED_MORPH             /* allow bed elevation to evolve */

/* Optional */
#define SED_DENS              /* sediment effects on seawater density */
#define ANA_SEDIMENT          /* analytical initial sediment fields */
```

### CPP Flag Reference

| Flag | Description |
|---|---|
| `SEDIMENT` | Master switch |
| `SUSPLOAD` | Suspended load transport |
| `BEDLOAD_MPM` | Meyer-Peter-Mueller bedload (currents only) |
| `BEDLOAD_SOULSBY` | Soulsby-Damgaard bedload (waves + currents) |
| `BEDLOAD_VANDERA` | van der A bedload (nearshore, wave asymmetry) |
| `SED_MORPH` | Bed elevation evolves with erosion/deposition |
| `SED_DENS` | Sediment affects seawater density |
| `SED_BIODIFF` | Biodiffusion of sediment in bed |
| `SSW_BBL` | Bottom boundary layer with wave-current interaction |
| `MB_BBL` | Meinte Blaas BBL |
| `SG_BBL` | Styles-Glenn BBL |

### The `sediment.in` Configuration File

#### Sediment Class Properties

Number of sediment classes is set at compile time: `NCS` (cohesive/mud) + `NNS` (non-cohesive/sand) = `NST` (total).

**Mud (cohesive) parameters:**

```
MUD_SD50   == 0.015d0      ! Median grain diameter (mm) — 15 µm
MUD_SRHO   == 2650.0d0     ! Grain density (kg/m³) — quartz
MUD_CSED   == 0.0d0        ! Initial concentration in water (kg/m³)
MUD_WSED   == 0.1d0        ! Settling velocity (mm/s)
MUD_TAU_CE == 0.05d0       ! Critical shear stress for erosion (N/m²)
MUD_TAU_CD == 0.05d0       ! Critical shear stress for deposition (N/m²)
MUD_ERATE  == 5.0d-4       ! Surface erosion rate (kg/m²/s)
MUD_BFRAC  == 0.5d0        ! Bed fraction (0–1)
```

**Sand (non-cohesive) parameters:**

```
SAND_SD50   == 0.25d0      ! Median grain diameter (mm) — fine sand
SAND_SRHO   == 2650.0d0    ! Grain density (kg/m³)
SAND_CSED   == 0.0d0       ! Initial concentration in water (kg/m³)
SAND_WSED   == 36.0d0      ! Settling velocity (mm/s)
SAND_TAU_CE == 0.16d0      ! Critical shear stress for erosion (N/m²)
SAND_TAU_CD == 0.16d0      ! Critical shear stress for deposition (N/m²)
SAND_ERATE  == 5.0d-4      ! Surface erosion rate (kg/m²/s)
SAND_BFRAC  == 0.5d0       ! Bed fraction (0–1)
```

#### Bed Layer Parameters

```
Nbed == 3                  ! Number of bed layers

! Initial thickness of each layer (m): active, sub-surface, deep
BTHK == 0.01d0  0.10d0  1.00d0

! Porosity of each layer
BPOR == 0.5d0  0.5d0  0.5d0

! New layer forms when deposition exceeds this thickness (m)
NEWLAYER_THICK == 0.001d0

! Bedload active layer thickness factor
BEDLOAD_COEFF == 0.05d0
```

### Settling Velocity Reference

| Sediment Type | Grain Size (mm) | Settling Velocity (mm/s) |
|---|---|---|
| Fine clay | 0.001–0.004 | 0.001–0.01 |
| Silt | 0.004–0.063 | 0.01–1.0 |
| Very fine sand | 0.063–0.125 | 1–10 |
| Fine sand | 0.125–0.25 | 10–25 |
| Medium sand | 0.25–0.50 | 25–60 |
| Coarse sand | 0.50–2.0 | 60–200 |

### Bedload Formulations

| Method | Best For | Accounts for Waves? |
|---|---|---|
| `BEDLOAD_MPM` | Rivers, estuaries, tidal channels | No |
| `BEDLOAD_SOULSBY` | Continental shelf, moderate wave influence | Yes |
| `BEDLOAD_VANDERA` | Surf zone, nearshore bars, wave-dominated | Yes (including wave asymmetry) |

### Wave-Sediment Coupling Feedback

When SWAN is coupled to ROMS with sediment enabled, a powerful feedback loop occurs:

```
SWAN provides wave orbital velocity → ROMS BBL computes enhanced bed stress →
Sediment is resuspended/transported → Bed elevation changes (SED_MORPH) →
Updated bathymetry sent back to SWAN → Wave refraction/shoaling changes →
Modified wave field affects stress → ...
```

Key fields from SWAN that drive sediment processes:
- `Hwave` — controls orbital velocities
- `Pwave_bot` — determines orbital velocity magnitude
- `Dwave` — determines wave-driven transport direction

Required CPP flags for wave-sediment interaction:

```c
#define SEDIMENT
#define SUSPLOAD
#define BEDLOAD_SOULSBY       /* or BEDLOAD_VANDERA */
#define SED_MORPH             /* if you want bed evolution */
#define SSW_BBL               /* wave-current BBL */
#define WEC_VF                /* vortex force */
#define WDISS_WAVEMOD
```

---

## 16. Data Download Sources

### Atmospheric Forcing (for WRF)

#### GFS (Global Forecast System)

| Resource | URL |
|---|---|
| NOMADS main page | https://nomads.ncep.noaa.gov/ |
| GFS 0.25° GRIB filter | https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl |
| GFS products | https://www.nco.ncep.noaa.gov/pmb/products/gfs/ |

```bash
# Download GFS 0.25° (example for a single forecast hour)
wget "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20240901/00/atmos/gfs.t00z.pgrb2.0p25.f012"
```

#### ERA5 Reanalysis

| Resource | URL |
|---|---|
| CDS portal | https://cds.climate.copernicus.eu/ |
| ERA5 single levels | https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels |
| ERA5 pressure levels | https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels |
| CDS API setup | https://cds.climate.copernicus.eu/how-to-api |

```bash
pip install cdsapi
# Create ~/.cdsapirc with your API key from the CDS portal
```

### Ocean IC/BC (for ROMS)

#### GLORYS12V1 (Copernicus Marine)

| Resource | URL |
|---|---|
| Product page | https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030 |
| Copernicus Marine Toolbox | https://pypi.org/project/copernicusmarine/ |

Resolution: 1/12° (~8 km), 50 vertical levels, 1993–present.

```bash
pip install copernicusmarine
copernicusmarine login
```

```python
import copernicusmarine

copernicusmarine.subset(
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    variables=["thetao", "so", "uo", "vo", "zos"],
    minimum_longitude=104, maximum_longitude=107,
    minimum_latitude=-8, maximum_latitude=-5,
    start_datetime="2024-09-01", end_datetime="2024-09-30",
    output_filename="glorys12v1_sunda.nc",
)
```

### Tidal Forcing (for ROMS)

#### TPXO9

| Resource | URL |
|---|---|
| TPXO main page | https://www.tpxo.net/ |
| Registration | https://www.tpxo.net/tpxo-products-and-registration |

Free for academic use. Register to download TPXO9-atlas-v5 NetCDF files. Use `otps2roms` scripts to interpolate tidal constituents onto your ROMS grid.

#### FES2022 (Alternative)

| Resource | URL |
|---|---|
| FES2022 release | https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes/release-fes22.html |
| PyFES Python package | https://github.com/CNES/aviso-fes |

Resolution: 1/30° (~3.7 km), 34 tidal constituents. Register at AVISO.

### Bathymetry

#### GEBCO

| Resource | URL |
|---|---|
| Download page | https://download.gebco.net/ |
| GEBCO 2025 info | https://www.gebco.net/data-products/gridded-bathymetry-data/ |

Resolution: 15 arc-seconds (~450 m). Download full global NetCDF (~8 GB) or use the interactive tool to subset your region.

### Wave Boundaries (for SWAN)

#### ERA5 Waves

Same CDS portal as atmospheric ERA5. Key variables: `significant_height_of_combined_wind_waves_and_swell`, `mean_wave_direction`, `mean_wave_period`. Resolution: 0.5°.

#### WAVEWATCH III Hindcast

| Resource | URL |
|---|---|
| WW3 hindcasts | https://polar.ncep.noaa.gov/waves/hindcasts/nopp-phase2.php |
| WW3 production data | https://polar.ncep.noaa.gov/waves/hindcasts/prod-multi_1.php |

Extract spectral or bulk wave parameters along your SWAN domain boundaries and convert to SWAN TPAR format (Hs, Tp, Dir, Dspr).

### Quick Reference

| Data | Source | URL |
|---|---|---|
| GFS 0.25° | NCEP NOMADS | https://nomads.ncep.noaa.gov/ |
| ERA5 | Copernicus CDS | https://cds.climate.copernicus.eu/ |
| GLORYS12V1 | Copernicus Marine | https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030 |
| TPXO9 tides | Oregon State | https://www.tpxo.net/ |
| FES2022 tides | AVISO/CNES | https://www.aviso.altimetry.fr/ |
| GEBCO bathymetry | GEBCO | https://download.gebco.net/ |
| ERA5 waves | Copernicus CDS | https://cds.climate.copernicus.eu/ |
| WW3 hindcast | NOAA NCEP | https://polar.ncep.noaa.gov/waves/ |

---

## 17. Project Directory Organization

Recommended directory structure for a COAWST project:

```
COAWST/                              # COAWST source code (git clone)
├── Projects/
│   └── Sunda_Strait/                # Your project directory
│       ├── sunda_strait.h           # Application header (CPP flags)
│       ├── coupling_sunda.in        # Coupling configuration
│       ├── ocean_sunda.in           # ROMS configuration
│       ├── swan_sunda.in            # SWAN input file
│       ├── namelist.input           # WRF configuration
│       ├── sediment_sunda.in        # Sediment configuration (if used)
│       ├── varinfo.dat              # ROMS variable metadata
│       ├── stations.in              # ROMS station locations (if used)
│       │
│       ├── Grid/                    # Grid files
│       │   ├── roms_grid.nc         # ROMS grid + bathymetry
│       │   ├── swan_coord.grd       # SWAN coordinates (curvilinear)
│       │   ├── swan_bathy.bot       # SWAN bathymetry
│       │   └── scrip_sunda.nc       # SCRIP interpolation weights
│       │
│       ├── Forcing/                 # Forcing data
│       │   ├── wrfinput_d01         # WRF initial conditions
│       │   ├── wrfbdy_d01           # WRF boundary conditions
│       │   ├── wrflowinp_d01        # WRF lower boundary (SST update)
│       │   ├── ocean_ini.nc         # ROMS initial conditions
│       │   ├── ocean_bry.nc         # ROMS boundary conditions
│       │   ├── ocean_frc.nc         # ROMS surface forcing (if BULK_FLUXES)
│       │   ├── ocean_tides.nc       # ROMS tidal forcing
│       │   └── swan_boundary.sp2    # SWAN spectral boundary
│       │
│       ├── Output/                  # Model output (gitignored)
│       │   ├── ocean_his.nc
│       │   ├── ocean_avg.nc
│       │   ├── ocean_rst.nc
│       │   ├── wrfout_d01_*
│       │   ├── wrfrst_d01_*
│       │   ├── swan_his.nc
│       │   └── swan_rst.dat
│       │
│       └── Scripts/                 # Pre/post-processing scripts
│           ├── make_roms_grid.py
│           ├── make_swan_grid.py
│           ├── make_forcing.py
│           ├── make_tides.py
│           ├── make_scrip_weights.sh
│           ├── run_coawst.sh        # Job submission script
│           └── plot_results.py
```

### Tips

- Keep source code (`COAWST/`) separate from output. Never modify source files for project-specific settings — use the header file and `.in` files.
- Add `Output/` to `.gitignore` — output files are large and regenerable.
- Use relative paths in `coupling_sunda.in` (relative to the directory where you launch `coawstM`).
- Store scripts that generate grids and forcing — reproducibility matters.
- For multiple experiments, create separate subdirectories under `Projects/` (e.g., `Sunda_Strait_tides_only/`, `Sunda_Strait_coupled/`).

---

## 18. Running a Test Case (Inlet_test)

Before attempting your own project, run the built-in **Inlet_test** to verify your build works. It's a small, fast test case that exercises ROMS+SWAN coupling with sediment transport.

### What It Simulates

An **idealized tidal inlet** — a narrow channel connecting the open ocean to a back-barrier bay:
- Sinusoidal ±1 m tidal signal at the north boundary (12-hour period)
- Constant JONSWAP waves from the north (Hs = 1.0 m, Tp = 10.0 s)
- One sand class (d50 = 1.0 mm) with bedload and morphodynamic updating
- Cartesian coordinates, no Coriolis, 12-hour simulation

| Parameter | Value |
|---|---|
| Grid | 75 × 70 × 8 levels (curvilinear) |
| DT (baroclinic) | 60 s |
| NDTFAST | 20 |
| Total time steps | 720 (= 12 hours) |
| Coupling interval | 600 s (10 minutes) |
| MPI processes | 2 (1 ROMS + 1 SWAN) |
| Runtime | 1–5 minutes on a modern machine |

### Available Variants

| Variant | Directory | Description |
|---|---|---|
| **Coupled** | `Projects/Inlet_test/Coupled/` | Default. ROMS + SWAN on same grid. |
| **DiffGrid** | `Projects/Inlet_test/DiffGrid/` | ROMS + SWAN on different grids (tests SCRIP interpolation). |
| **Refined** | `Projects/Inlet_test/Refined/` | Nested/refined grids. |
| **InWave** | `Projects/Inlet_test/InWave/` | Uses built-in InWave instead of SWAN. |
| **Swanonly** | `Projects/Inlet_test/Swanonly/` | SWAN-only with nested grids, no ocean. |
| **WW3** | `Projects/Inlet_test/WW3/` | ROMS + WAVEWATCH III coupling. |

### Step-by-Step: Build and Run (Coupled variant)

#### Step 1: Clone

```bash
git clone https://github.com/DOI-USGS/COAWST.git
cd COAWST
```

#### Step 2: Edit `build_coawst.sh`

Set these variables:

```bash
COAWST_APPLICATION=INLET_TEST          # must be uppercase
MY_ROOT_DIR=/path/to/COAWST
MY_HEADER_DIR=${MY_ROOT_DIR}/Projects/Inlet_test/Coupled
MY_ANALYTICAL_DIR=${MY_ROOT_DIR}/Projects/Inlet_test/Coupled

FORT=gfortran                           # or ifort
USE_MPI=on
USE_MPIF90=on
which_MPI=openmpi                       # or mpich, intel
USE_NETCDF4=on
```

#### Step 3: Compile

```bash
./build_coawst.sh -j 4
```

This produces the executable `coawstM` in the root directory.

#### Step 4: Run

```bash
mpirun -np 2 ./coawstM Projects/Inlet_test/Coupled/coupling_inlet_test.in > inlet_test.log
```

**Critical:** Run from the COAWST root directory. All file paths in the input files are relative to this location.

#### Step 5: Verify

Check the log for:
- ROMS initialization header with grid dimensions
- SWAN initialization messages
- Time-stepping output (KINETIC_ENRG, POTEN_ENRG columns)
- Coupling exchange messages every 600 s (10 exchanges total)
- Clean termination

Check output files:
- `ocean_his_coupled.nc` — history snapshots
- `ocean_avg_coupled.nc` — time-averaged fields
- `*.mat` — SWAN block outputs (Hsig, depth, etc.)

### Key Files in Coupled/

| File | Purpose |
|---|---|
| `inlet_test.h` | CPP flags (ROMS_MODEL, SWAN_MODEL, MCT_LIB, SEDIMENT, WEC_VF, etc.) |
| `ocean_inlet_test.in` | ROMS config (grid, time stepping, physics, output) |
| `swan_inlet_test.in` | SWAN config (grid, physics, boundaries, output) |
| `coupling_inlet_test.in` | Coupling config (NnodesWAV=1, NnodesOCN=1, TI=600s) |
| `sediment_inlet_test.in` | Sediment class properties |
| `inlet_test_grid.nc` | ROMS grid NetCDF |
| `inlet_test_grid_coord.grd` | SWAN coordinate file |
| `inlet_test_bathy.bot` | SWAN bathymetry file |

### Common Beginner Issues

| Issue | Solution |
|---|---|
| NetCDF not found | Verify `nf-config --all` works. Set `USE_NETCDF4=on`. |
| `mpif90` not found | Install MPI and ensure it's on PATH. |
| SCRIP weight file missing | Pre-generated file should be in the repo. If not, build and run `Lib/SCRIP_COAWST/scrip_coawst`. |
| Wrong MPI count | `mpirun -np` must equal `NnodesWAV + NnodesOCN` (default 2). |
| `NtileI × NtileJ ≠ NnodesOCN` | Both must be consistent — default is 1×1 = 1. |
| File not found errors | Run from the COAWST root, not a subdirectory. |
| `COAWST_APPLICATION` case mismatch | Must be `INLET_TEST` (uppercase) in build script. |

---

## 19. SWAN Boundary File Preparation

### TPAR Format (Parametric)

The simplest boundary format — bulk wave parameters at each time step. One file per boundary point.

```
TPAR
20240901.000000  1.500   8.50  225.00   25.0
20240901.030000  1.620   8.80  228.00   24.0
20240901.060000  1.800   9.20  230.00   22.0
20240901.090000  2.100  10.00  235.00   20.0
```

| Column | Parameter | Unit |
|---|---|---|
| 1 | Time | `YYYYMMDD.HHMMSS` |
| 2 | Hs | m (significant wave height) |
| 3 | Tp | s (peak wave period) |
| 4 | Dir | degrees (direction waves come FROM, clockwise from N) |
| 5 | Dspr | degrees (directional spreading width) |

First line must be `TPAR`. SWAN interpolates linearly between time steps.

### Spectral Format (.sp2)

Full 2D directional-frequency spectrum. More accurate for multi-modal seas (combined swell + wind sea):

```
SWAN   1                                Swan standard spectral file
TIME                                    time-dependent data
     1                                  time coding option
LONLAT                                  locations in lon/lat
     3                                  number of locations
  105.0000   -7.0000
  105.0000   -6.5000
  105.0000   -6.0000
APTS                                    absolute spectral densities
    30                                  number of frequencies
  0.0400  0.0440  0.0484  ...           (frequency bins in Hz)
    36                                  number of directions
  0.0000  10.0000  20.0000  ...         (direction bins in degrees)
QUANT
     1
VaDens                                  variance densities
m2/Hz/degr                              unit
  -0.9900E+02                           exception value
20240901.000000                         date and time
LOCATION     1
  0.0012  0.0015  0.0008  ...           (nfreq × ndir spectral values)
LOCATION     2
  ...
```

### BOUNDSPEC Command

Place `BOUND SHAPE` before `BOUNDSPEC`:

```
BOUND SHAPE JONSWAP 3.3 PEAK DSPR DEGREES
```

#### Constant Spectrum Along One Side

```
BOUNDSPEC SIDE WEST CONSTANT FILE 'tpar_west.txt' 1
```

#### Variable Spectrum (Multiple Points Along a Side)

```
BOUNDSPEC SIDE WEST CCW VARIABLE FILE &
    'tpar_west/tpar_bnd001.txt' 1 &
    'tpar_west/tpar_bnd002.txt' 2 &
    'tpar_west/tpar_bnd003.txt' 3 &
    'tpar_west/tpar_bnd004.txt' 4 &
    'tpar_west/tpar_bnd005.txt' 5
```

`CCW` = counter-clockwise ordering. For `SIDE WEST CCW`, points go from SW corner to NW corner.

#### Multiple Boundaries

```
! Western boundary (5 TPAR files, S to N)
BOUNDSPEC SIDE WEST CCW VARIABLE FILE &
    'tpar_west/tpar_001.txt' 1 &
    'tpar_west/tpar_002.txt' 2 &
    'tpar_west/tpar_003.txt' 3 &
    'tpar_west/tpar_004.txt' 4 &
    'tpar_west/tpar_005.txt' 5

! Southern boundary (5 TPAR files, W to E)
BOUNDSPEC SIDE SOUTH CCW VARIABLE FILE &
    'tpar_south/tpar_001.txt' 1 &
    'tpar_south/tpar_002.txt' 2 &
    'tpar_south/tpar_003.txt' 3 &
    'tpar_south/tpar_004.txt' 4 &
    'tpar_south/tpar_005.txt' 5

! Northern boundary - constant small waves (sheltered)
BOUNDSPEC SIDE NORTH CONSTANT PAR 0.5 6.0 0.0 30.0

! Eastern boundary - land, no boundary needed (default = no incoming energy)
```

#### Using a Spectral File

```
BOUNDSPEC SIDE WEST CCW VARIABLE FILE 'swan_west_boundary.sp2' 1
```

### Python: ERA5 to SWAN TPAR

```python
import numpy as np
from netCDF4 import Dataset, num2date
from scipy.interpolate import interpn
import os

def era5_to_tpar(nc_file, boundary_points, output_dir='./tpar_files',
                  default_dspr=25.0):
    """
    Convert ERA5 wave data to SWAN TPAR files.

    Parameters
    ----------
    nc_file : str
        ERA5 netCDF file with swh, pp1d (or mwp), mwd.
    boundary_points : list of (lon, lat) tuples
        Points along the SWAN boundary.
    """
    os.makedirs(output_dir, exist_ok=True)

    ds = Dataset(nc_file)
    lons = np.array(ds.variables['longitude'][:])
    lats = np.array(ds.variables['latitude'][:])
    time_var = ds.variables['time']
    times = num2date(time_var[:], units=time_var.units,
                      calendar='standard')

    hs = np.array(ds.variables['swh'][:])
    # Use peak period if available, else mean period
    if 'pp1d' in ds.variables:
        tp = np.array(ds.variables['pp1d'][:])
    else:
        tp = np.array(ds.variables['mwp'][:]) * 1.1  # Tm01 → Tp approx
    mwd = np.array(ds.variables['mwd'][:])
    ds.close()

    # ERA5 lats are often descending — fix for interpolation
    if lats[0] > lats[-1]:
        lats = lats[::-1]
        hs = hs[:, ::-1, :]
        tp = tp[:, ::-1, :]
        mwd = mwd[:, ::-1, :]

    pts = np.array([(lat, lon) for lon, lat in boundary_points])

    tpar_files = []
    for ip, (blon, blat) in enumerate(boundary_points):
        fname = os.path.join(output_dir, f'tpar_bnd{ip+1:03d}.txt')
        tpar_files.append(fname)

        with open(fname, 'w') as f:
            f.write('TPAR\n')
            for it, t in enumerate(times):
                hs_val = max(float(interpn((lats, lons), hs[it],
                         [[blat, blon]])[0]), 0.01)
                tp_val = max(float(interpn((lats, lons), tp[it],
                         [[blat, blon]])[0]), 1.0)

                # Circular interpolation for direction (avoid 0/360 wrap)
                dir_rad = np.deg2rad(mwd[it])
                sin_val = float(interpn((lats, lons), np.sin(dir_rad),
                          [[blat, blon]])[0])
                cos_val = float(interpn((lats, lons), np.cos(dir_rad),
                          [[blat, blon]])[0])
                dir_val = np.rad2deg(np.arctan2(sin_val, cos_val)) % 360

                tstr = t.strftime('%Y%m%d.%H%M%S')
                f.write(f'{tstr} {hs_val:7.3f} {tp_val:7.2f} '
                        f'{dir_val:7.2f} {default_dspr:6.1f}\n')

        print(f'  Wrote {fname} ({len(times)} time steps)')
    return tpar_files

# === Usage for Sunda Strait ===
# Western boundary points (S to N along 104.5°E)
west_pts = [(104.5, lat) for lat in np.arange(-7.5, -5.4, 0.5)]
# Southern boundary points (W to E along 7.5°S)
south_pts = [(lon, -7.5) for lon in np.arange(104.5, 106.6, 0.5)]

era5_to_tpar('era5_waves_202409.nc', west_pts, './tpar_west')
era5_to_tpar('era5_waves_202409.nc', south_pts, './tpar_south')
```

> **Direction convention:** ERA5 uses meteorological convention (direction waves come FROM, clockwise from N) — same as SWAN's `NAUTICAL`. No conversion needed.

### COAWST MATLAB Tools for WW3 Boundaries

COAWST includes tools in `Tools/mfiles/` for creating SWAN boundaries from WAVEWATCH III:

| Script | Purpose |
|---|---|
| `create_swanTpar_from_WW3.m` | WW3 bulk parameters → SWAN TPAR files |
| `create_swan_2dspec_from_ww3.m` | WW3 2D spectra → SWAN `.sp2` files |

**Unit conversion for WW3 spectra:** WW3 outputs spectral density in m²/Hz/rad. SWAN expects m²/Hz/deg. Multiply by π/180.

**Direction conversion for WW3:** WW3 typically uses "going to" convention. SWAN uses "coming from". Add 180° to WW3 directions.

### SWAN-to-SWAN Nesting (NESTout / BOUNDNEST)

For one-way nesting from a coarse SWAN run to a fine run:

**Coarse run** — output nested boundary spectra:

```
NGRID 'fine' 104.5 -7.5 0.0 2.0 2.0 200 200
NESTOUT 'fine' 'nest_spectra.sp2' OUTPUT 20240901.000000 1 HR
```

**Fine run** — read nested boundary spectra:

```
BOUNDNEST1 NEST 'nest_spectra.sp2' CLOSED
```

| Command | Purpose |
|---|---|
| `BOUNDNEST1` | Nesting from SWAN (reads SWAN `.sp2` format) |
| `BOUNDNEST2` | Nesting from WAM |
| `BOUNDNEST3` | Nesting from WAVEWATCH III |
| `CLOSED` | Apply boundary on all 4 sides |
| `OPEN` | Leave shoreward side open |

---

## 20. Nesting in COAWST

COAWST supports multiple grids for each component model. This enables higher resolution in areas of interest while keeping the outer domain coarser.

### Multiple SWAN Grids

The most common nesting in COAWST. Each SWAN grid runs as a separate instance.

#### Configuration

In `coupling.in`, list multiple SWAN input files:

```
NnodesWAV = 4              ! total SWAN processors (split among grids)

WAV_name = Projects/Sunda/swan_parent.in \
           Projects/Sunda/swan_child.in
```

In `scrip_coawst.in`, specify multiple SWAN grids:

```
SWAN_GRIDS = 2

SWAN_COORD1 = 'swan_parent_coord.grd'
SWAN_BATH1  = 'swan_parent_bathy.bot'
SWAN_NUMX1  = 201
SWAN_NUMY1  = 171
CARTESIAN1  = 0

SWAN_COORD2 = 'swan_child_coord.grd'
SWAN_BATH2  = 'swan_child_bathy.bot'
SWAN_NUMX2  = 301
SWAN_NUMY2  = 251
CARTESIAN2  = 0
```

SCRIP generates interpolation weights for all grid pairs automatically.

#### How It Works

- The parent SWAN grid couples with WRF/ROMS at the specified interval
- The child SWAN grid also couples with ROMS
- SWAN grids can have different resolutions; the child is typically finer
- Each SWAN grid receives wind from WRF and currents/water levels from ROMS
- Each SWAN grid sends wave fields back to ROMS

### Multiple ROMS Grids (Nesting)

ROMS supports two nesting modes:

| Mode | CPP Flag | Description |
|---|---|---|
| **Refinement** | `NESTING` | Parent and child share the same domain region; child has finer resolution. Two-way feedback optional. |
| **Composition** | `NESTING` + `COMPOSED` | Child covers a different region (e.g., adjacent domains). |

#### Configuration for Nested ROMS in COAWST

In the header file:

```c
#define NESTING           /* enable nesting */
/* #define COMPOSED */    /* for composition mode */
```

In `ocean.in`, define multiple grids:

```
! Parent grid: Lm=150, Mm=120, N=30
! Child grid:  Lm=200, Mm=160, N=30 (finer resolution)
Ngrids = 2

NtileI == 4  2             ! I-tiles for each grid
NtileJ == 4  2             ! J-tiles for each grid
```

In `coupling.in`:

```
NnodesOCN = 20             ! total ROMS processors
                           ! parent: NtileI(1)*NtileJ(1) = 16
                           ! child:  NtileI(2)*NtileJ(2) = 4
```

#### Which Grid Couples to SWAN/WRF?

By default, all ROMS grids exchange data with SWAN and WRF. The SCRIP weight file must contain weights for each ROMS grid paired with each SWAN/WRF grid.

### WRF Nested Domains

WRF nesting is handled internally by WRF (not by MCT). The coupling to ROMS/SWAN typically uses the **finest WRF domain** that overlaps with the ocean/wave grids.

In `namelist.input`:

```
&domains
 max_dom     = 2,
 dx          = 9000, 3000,
 dy          = 9000, 3000,
 parent_id   = 1, 1,
 i_parent_start = 1, 40,
 j_parent_start = 1, 30,
 parent_grid_ratio = 1, 3,
/
```

The SCRIP weight file must match the WRF domain that couples to ROMS/SWAN. In `scrip_coawst.in`:

```
WRF_GRIDS = 1
WRF_GRIDS1 = 'wrfinput_d02'     ! couple using the fine domain
```

### Processor Allocation with Nesting

The total `mpirun -np N` must equal the sum of all model processors:

```
NnodesATM = 64    ! WRF (handles all WRF domains internally)
NnodesWAV = 8     ! SWAN (split among SWAN grids)
NnodesOCN = 48    ! ROMS (split among ROMS grids via NtileI/NtileJ)
NnodesHYD = 0

! Total: 64 + 8 + 48 = 120
! mpirun -np 120 ./coawstM coupling.in
```

**Processor split for nested ROMS:**
- ROMS distributes processors among grids based on `NtileI × NtileJ` per grid
- The sum across all grids must equal `NnodesOCN`

**Processor split for multiple SWAN grids:**
- SWAN distributes `NnodesWAV` processors among grids proportionally
- Cannot control per-grid allocation directly in `coupling.in`

### Example: Sunda Strait with Nesting

A practical nested setup:

```
WRF:  d01 (9 km) + d02 (3 km) — WRF handles nesting internally
ROMS: parent (3 km, wider Indonesian waters) + child (1 km, Sunda Strait)
SWAN: parent (3 km, matching ROMS parent) + child (1 km, matching ROMS child)

Processor allocation:
  NnodesATM = 64
  NnodesWAV = 8     (4 per SWAN grid)
  NnodesOCN = 48    (36 parent + 12 child)
  Total: 120

Coupling: WRF d02 ↔ both ROMS grids ↔ both SWAN grids
```

### Tips for Nesting

- **Start simple**: Get the unnested single-grid case working first, then add nesting.
- **Grid alignment**: The child grid boundaries should align with parent grid cells for clean interpolation.
- **Time step**: The child grid may need a smaller time step than the parent. Set `DT` per grid in `ocean.in`.
- **Coupling interval**: Use the same coupling interval for all grids, or match each to its time step.
- **SCRIP weights**: Must include all grid pairs. The weight file grows with the number of grids.
- **Debugging**: Test each model's nesting independently (ROMS-only, SWAN-only) before coupling.

---

## 21. WRF Setup in COAWST (8-Step Workflow)

> Based on COAWST 2021 Training Day 1 — WRF Setup and Configuration (John Warner, Maitane Olabarrieta)

### Overview

WRF (Weather Research and Forecasting) is the atmospheric component of COAWST. Setting up WRF follows an 8-step workflow using the WRF Preprocessing System (WPS) and the WRF model itself.

### The 8-Step Workflow

```
Step 1: Define domain (namelist.wps)
Step 2: geogrid.exe → Static terrain/land-use data
Step 3: Download met data (GFS/ERA5/NARR)
Step 4: ungrib.exe → Extract met fields from GRIB
Step 5: metgrid.exe → Interpolate to WRF grid
Step 6: real.exe → Create wrfinput_d01, wrfbdy_d01
Step 7: wrf.exe → Run WRF standalone (test)
Step 8: Couple with COAWST
```

### Step 1: Define Domain (`namelist.wps`)

Key parameters in `&share` and `&geogrid` sections:

```
&share
  max_dom = 2,                    ! number of domains
  start_date = '2019-09-01_00:00:00', '2019-09-01_00:00:00',
  end_date   = '2019-09-05_00:00:00', '2019-09-05_00:00:00',
  interval_seconds = 21600        ! input data interval (6h = GFS)
/

&geogrid
  parent_id         = 1, 1,
  parent_grid_ratio = 1, 3,       ! refinement ratio
  i_parent_start    = 1, 31,      ! child start position in parent (i)
  j_parent_start    = 1, 17,      ! child start position in parent (j)
  e_we          = 74, 112,        ! grid points in x (west-east)
  e_sn          = 61, 97,         ! grid points in y (south-north)
  dx = 27000,                     ! parent grid spacing (meters)
  dy = 27000,
  map_proj = 'mercator',          ! projection
  ref_lat  = 30.0,                ! center latitude
  ref_lon  = -80.0,               ! center longitude
  truelat1 = 30.0,                ! true latitude 1
  geog_data_res = 'default', 'default',
  geog_data_path = '/path/to/WPS_GEOG/'
/
```

**Grid spacing rules:**
- Parent `dx`/`dy` in meters
- Child grid spacing = parent / `parent_grid_ratio`
- Child grid dimensions must satisfy: `(e_we_child - 1)` must be divisible by `parent_grid_ratio`

### Step 2: Run `geogrid.exe`

Generates static geographic data files (`geo_em.d01.nc`, `geo_em.d02.nc`):
- Terrain height (HGT_M)
- Land use categories (LU_INDEX)
- Soil types, vegetation fraction, albedo
- Land mask (LANDMASK)

```bash
./geogrid.exe >& geogrid.log
# Verify: check for "Successful completion of geogrid"
ncdump -h geo_em.d01.nc   # inspect output
```

### Step 3: Download Meteorological Data

Common sources:
| Source | Resolution | Interval | Use Case |
|--------|-----------|----------|----------|
| GFS | 0.25° | 6h (3h available) | Operational forecasting |
| ERA5 | 0.25° | 1h | Reanalysis, hindcast |
| NARR | 32 km | 3h | North America regional |
| FNL | 1° | 6h | Historical cases |

Download in GRIB2 format. GFS from NCEP NOMADS, ERA5 from Copernicus CDS.

### Step 4: Run `ungrib.exe`

Extracts meteorological fields from GRIB files:

```bash
# Link GRIB files
./link_grib.csh /path/to/gfs_files/gfs*

# Link the appropriate Vtable
ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable

# Run ungrib
./ungrib.exe >& ungrib.log
# Output: FILE:YYYY-MM-DD_HH intermediate files
```

**Vtable** maps GRIB fields to WPS intermediate format. Different Vtables for GFS, ERA5, NARR, etc.

### Step 5: Run `metgrid.exe`

Horizontally interpolates meteorological data onto WRF grid:

```bash
./metgrid.exe >& metgrid.log
# Output: met_em.d01.YYYY-MM-DD_HH:00:00.nc (one per time step per domain)
```

### Step 6: Run `real.exe`

Vertical interpolation + boundary/initial condition generation:

```bash
# Set namelist.input for real.exe
# Key: num_metgrid_levels must match met_em files
mpirun -np 1 ./real.exe
# Output: wrfinput_d01, wrfinput_d02, wrfbdy_d01
```

Key `namelist.input` settings:
```
&domains
  time_step = 120,               ! main time step (seconds), ~6*dx(km)
  e_vert = 45, 45,               ! vertical levels (same for all domains)
  num_metgrid_levels = 34,       ! from met_em files
  eta_levels = 1.000, 0.9975, ... 0.000   ! optional explicit levels
/
```

### Step 7: Run WRF Standalone (Test)

```bash
mpirun -np 4 ./wrf.exe
# Output: wrfout_d01_YYYY-MM-DD_HH:00:00
# Verify: check rsl.error.0000 for completion message
```

### Step 8: Couple with COAWST

When running in COAWST:
- WRF is compiled as part of the COAWST build system
- `wrfinput_d0X` and `wrfbdy_d01` are still needed
- WRF exchanges fields with ROMS/SWAN via MCT/ESMF coupler
- Set `NnodesATM` in `coupling_coawst.in` for processor allocation

### Moving Nests in WRF

WRF supports **vortex-following moving nests** for hurricane tracking:

```
&domains
  num_moves = -99,               ! automatic vortex tracking
  move_id   = 2,                 ! which domain moves
  move_interval = 15,            ! check interval (minutes)
  move_cd_x = 0, 1,              ! movement flags
  move_cd_y = 0, 1,
/
```

- The child domain automatically follows the storm center
- Useful for hurricane/typhoon simulations
- Moving nest in COAWST: ROMS/SWAN grids remain fixed, only WRF nest moves

---

## 22. WRF Map Projections

> Based on COAWST 2021 Training Day 1 — Map Projections for WRF

### Available Projections

WRF supports four map projections set via `map_proj` in `namelist.wps`:

| Projection | `map_proj` | Best For | Key Parameters |
|------------|-----------|----------|----------------|
| **Lambert Conformal Conic** | `'lambert'` | Mid-latitudes | `truelat1`, `truelat2` |
| **Mercator** | `'mercator'` | Tropics, equatorial regions | `truelat1` |
| **Polar Stereographic** | `'polar'` | Polar regions | `truelat1` |
| **Lat-Lon (Cylindrical Equidistant)** | `'lat-lon'` | Global or large tropical domains | None special |

### Lambert Conformal Conic

- **Standard for mid-latitude applications** (30°-60° N/S)
- Conformal: preserves angles and local shapes
- Two true latitudes where scale is exact
- Distortion increases away from true latitudes

```
map_proj = 'lambert',
truelat1 = 30.0,        ! first true latitude
truelat2 = 60.0,        ! second true latitude
stand_lon = -98.0,       ! standard longitude (center meridian)
ref_lat = 40.0,
ref_lon = -98.0,
```

### Mercator

- **Best for tropical/equatorial domains** (within ±30°)
- Conformal projection
- Constant grid spacing along parallels
- Scale distortion increases toward poles
- **Recommended for Indonesian/tropical COAWST applications**

```
map_proj = 'mercator',
truelat1 = 0.0,         ! true latitude (equator for minimal distortion)
ref_lat = -6.0,
ref_lon = 106.0,
```

### Polar Stereographic

- **Best for high-latitude/polar domains**
- Conformal at the pole
- True scale at one latitude

```
map_proj = 'polar',
truelat1 = 60.0,        ! true latitude
ref_lat = 90.0,
ref_lon = 0.0,
```

### Lat-Lon (Cylindrical Equidistant)

- **For global or very large domains**
- Grid spacing in degrees, not meters
- `dx`/`dy` specified in degrees
- Not conformal — significant shape distortion away from equator

```
map_proj = 'lat-lon',
dx = 0.25,              ! grid spacing in degrees
dy = 0.25,
ref_lat = 0.0,
ref_lon = 180.0,
```

### Choosing the Right Projection

```
Latitude of domain center:
  |
  ├── > 60° → Polar Stereographic
  ├── 30°-60° → Lambert Conformal
  ├── < 30° → Mercator
  └── Global/very large → Lat-Lon
```

**For Indonesia/Sunda Strait:** Use **Mercator** with `truelat1` near the equator.

---

## 23. WaveWatch III in COAWST

> Based on COAWST 2021 Training Day 3 — WaveWatch III (Maitane Olabarrieta)

### Overview

WaveWatch III (WW3) is an alternative wave model to SWAN in COAWST. Developed by NOAA/NCEP, it's designed for **regional to global scale** wave modeling.

**SWAN vs WW3:**
| Feature | SWAN | WW3 |
|---------|------|-----|
| Scale | Nearshore/coastal | Regional/global |
| Grid | Structured | Structured + unstructured |
| Physics | Depth-limited, triads | Deep water, global |
| Coupling | Well-tested in COAWST | Available since COAWST v3.4+ |
| Best for | Small coastal domains | Large open-ocean domains |

### WW3 Switch File

WW3 uses a **switch file** to configure compile-time options (similar to ROMS CPP flags):

```
F90 NOGRB DIST MPI PR3 UQ FLX0 LN1 ST4 STAB0 NL1 BT1 DB1 MLIM
TR0 BS0 IC0 IS0 REF0 XX0 WNT1 WNX1 RWND CRT1 CRN1 O0 O1 O2 O3
O4 O5 O6 O7 O14 O15
```

Key switches:
| Switch | Description |
|--------|-------------|
| `ST4` | Source term package (Ardhuin et al. 2010) |
| `NL1` | DIA nonlinear interactions |
| `BT1` | JONSWAP bottom friction |
| `DIST MPI` | Distributed memory parallel |
| `PR3` | 3rd-order propagation |
| `FLX0` | Flux computation (default) |
| `LN1` | Linear input (Cavaleri & Malanotte-Rizzoli) |
| `MLIM` | Miche-style depth limiting |

### WW3 Grid Creation

```bash
# 1. Define grid in ww3_grid.inp
#    - Spatial grid (resolution, extent)
#    - Spectral grid (frequencies, directions)
#    - Bathymetry
#    - Time steps

# 2. Run grid preprocessor
ww3_grid
# Output: mod_def.ww3 (model definition file)
```

Spectral grid settings:
```
# Frequency range and number of bins
  0.0345  1.10  32        ! f_low, f_ratio, n_freq
# Direction: 36 bins = 10° resolution
  36   0.0                 ! n_dir, first_dir
```

### WW3 Preprocessors (4 Required)

```
1. ww3_grid  → mod_def.ww3      (grid/spectral definition)
2. ww3_strt  → restart.ww3      (initial conditions)
3. ww3_bound → nest.ww3         (boundary spectra, if nested)
4. ww3_prnc  → wind.ww3         (wind forcing input)
```

**`ww3_prnc`** converts NetCDF wind fields to WW3 format:
```
# ww3_prnc.inp example
'WND' 'T' 'T'                    ! field type, time interp, space interp
 0. 0. 0. 0.                     ! time offset
'wind_forcing.nc'                 ! input file
'longitude' 'latitude' 'uwnd' 'vwnd'  ! variable names
```

### Coupling WW3 with COAWST

In the COAWST header file (`.h`):
```c
#define WW3_MODEL          /* Use WW3 instead of SWAN */
#undef  SWAN_MODEL         /* Disable SWAN */
#define MCT_LIB            /* MCT coupling */
```

In `coupling_coawst.in`:
```
NnodesWAV = 16             ! processors for WW3
WAV_model = WW3            ! specify WW3
```

**Key difference from SWAN coupling:**
- WW3 provides wave parameters (Hs, Tp, Dir, Tm) to ROMS
- ROMS provides currents and water levels to WW3
- WRF provides 10m winds to WW3
- WW3 can provide wave-dependent roughness to WRF

### Running WW3 in COAWST

```bash
# 1. Prepare WW3 inputs (grid, restart, boundary, wind)
# 2. Set switches and compile COAWST with WW3_MODEL defined
# 3. Configure coupling_coawst.in
# 4. Run as usual:
mpirun -np 64 coawstM coupling_coawst.in
```

---

## 24. ROMS Build System & External Libraries

> Based on COAWST 2021 Training Day 2 — ROMS Build System (Hernan Arango)

### COAWST/ROMS Build System

The COAWST build uses a shell script (`build_coawst.sh`) that configures compilers, paths, and flags before invoking `make`.

### `build_coawst.sh` Key Variables

```bash
# Application
export   COAWST_APPLICATION=INLET_TEST    # Your application name
export   MY_HEADER_DIR=${MY_ROOT_DIR}/Projects/Inlet_test

# Compiler
export   FORT=ifort          # ifort, gfortran, pgi
export   USE_MPI=on          # Enable MPI
export   USE_MPIF90=on       # Use mpif90 wrapper

# Parallel options
export   USE_OpenMP=          # Leave blank to disable
export   USE_LARGE=on         # >2GB files (64-bit offsets)
export   USE_NETCDF4=on       # NetCDF-4/HDF5 support

# Paths
export   MY_ROOT_DIR=/home/user/COAWST
export   NETCDF_INCDIR=/usr/local/include
export   NETCDF_LIBDIR=/usr/local/lib
```

### External Libraries

COAWST depends on several external libraries, managed in `Lib/` directory:

| Library | Purpose | Required? |
|---------|---------|-----------|
| **NetCDF-C + NetCDF-Fortran** | I/O for all model components | Yes |
| **MCT** | Model Coupling Toolkit (coupler) | Yes (if coupling) |
| **ARPACK** | Eigenvalue solver for 4D-Var, adjoint, tangent linear | Only for data assimilation |
| **ESMF/NUOPC** | Alternative coupling framework | Optional (instead of MCT) |
| **PIO/SCORPIO** | Parallel I/O library | Optional (for large runs) |
| **HDF5** | Hierarchical Data Format (used by NetCDF-4) | Yes (if USE_NETCDF4=on) |

### MCT Build

MCT is built as part of the COAWST compilation:

```bash
cd Lib/MCT
./configure FC=ifort CC=icc --prefix=/path/to/install
make
make install
```

### ARPACK Build

Only needed for ROMS data assimilation modes (4D-Var, adjoint):

```bash
cd Lib/ARPACK
# Edit ARmake.inc: set FC, FFLAGS, ARPACKLIB path
make lib       # builds libarpack.a
```

### Compilation Flow

```
build_coawst.sh
  ├── Set environment variables
  ├── Source compiler-specific makefile (Linux-ifort.mk, etc.)
  ├── Build MCT (if coupled)
  ├── Build ARPACK (if 4D-Var)
  ├── Build PIO (if enabled)
  ├── Preprocess CPP flags from header file (.h)
  └── make -j N coawstM
       ├── Compile ROMS modules
       ├── Compile SWAN/WW3 (if enabled)
       ├── Compile WRF (if enabled)
       └── Link → coawstM executable
```

### Compiler-Specific Notes

| Compiler | Makefile | Notes |
|----------|----------|-------|
| Intel `ifort` | `Linux-ifort.mk` | Recommended, best performance |
| GNU `gfortran` | `Linux-gfortran.mk` | Free, widely available |
| PGI/NVIDIA `pgf90` | `Linux-pgi.mk` | GPU support possible |

### Common Build Issues

- **NetCDF not found**: Ensure `NETCDF_INCDIR` and `NETCDF_LIBDIR` point to correct paths. Check `nc-config --includedir` and `nc-config --libdir`.
- **MPI errors**: Make sure `mpif90 --version` works and matches the Fortran compiler.
- **MCT build fails**: Often a compiler mismatch. MCT must be built with the same Fortran compiler as COAWST.
- **Undefined symbols at link**: Usually missing libraries. Check `-L` paths and `-l` flags in the makefile.

---

## 25. Parallel I/O (PIO)

> Based on COAWST 2021 Training Day 2 — Parallel I/O (Hernan Arango)

### Why Parallel I/O?

Standard ROMS I/O has a bottleneck: **only the master processor** reads/writes NetCDF files. All other processors send/receive data to/from the master via MPI. For large grids and many processors, this becomes a severe performance limitation.

**PIO** (Parallel I/O) allows multiple processors to read/write simultaneously, dramatically improving I/O performance for large simulations.

### I/O Strategies in ROMS

ROMS supports 4 I/O strategies:

| Strategy | CPP Flag | Description |
|----------|----------|-------------|
| **Serial** | (default) | Master PE does all I/O |
| **PIO (Parallel I/O)** | `PIO_LIB` | Multiple I/O PEs via PIO library |
| **SCORPIO** | `SCORPIO_LIB` | Successor to PIO, better async |
| **Parallel NetCDF** | `PARALLEL_IO` | Direct parallel NetCDF-4 (HDF5) |

### PIO Configuration in `roms.in`

```
! I/O strategy
   Nway = 4          ! Number of I/O processors
   IOrec = 0         ! I/O record (0 = auto)

! PIO-specific
   PIO_NUMTASKS = 4           ! I/O processor count
   PIO_STRIDE = 16            ! Stride between I/O PEs
   PIO_BASE = 0               ! First I/O PE rank
   PIO_TYPENAME = "pnetcdf"   ! pnetcdf or netcdf4p
   PIO_REARRANGER = "box"     ! box or subset
```

### Synchronous vs Asynchronous I/O

- **Synchronous**: I/O PEs are compute PEs that pause to do I/O → simpler but computation stalls during writes
- **Asynchronous**: Dedicated I/O PEs that only handle I/O → better overlap but uses extra processors

```
Synchronous:   [Compute + I/O] [Compute + I/O] [Compute + I/O] ...
Asynchronous:  [Compute PEs] → hand off → [Dedicated I/O PEs]
```

### When to Use PIO

| Scenario | Recommendation |
|----------|---------------|
| Small grid (<500×500), few PEs (<32) | Standard serial I/O |
| Medium grid, 32-128 PEs | Consider PIO with 4 I/O PEs |
| Large grid (>1000×1000), >128 PEs | PIO or SCORPIO strongly recommended |
| High-frequency output | PIO + async for minimal compute impact |

### Build with PIO

```bash
# In build_coawst.sh:
export USE_PIO=on
export PIO_INCDIR=/path/to/pio/include
export PIO_LIBDIR=/path/to/pio/lib

# In header file:
#define PIO_LIB
```

---

## 26. Vegetation Module

> Based on COAWST 2021 Training Day 4 — Vegetation (Tarandeep Kalra)

### Overview

The COAWST vegetation module simulates the effects of **submerged aquatic vegetation (SAV)** and **marsh vegetation** on hydrodynamics, waves, and sediment transport.

### Vegetation Effects

```
Vegetation influences:
  ├── Hydrodynamics
  │   ├── Drag force on currents (reduces flow velocity)
  │   ├── Turbulence modification (TKE production + dissipation)
  │   └── Flexible plant bending (effective height reduction)
  ├── Waves
  │   ├── Wave energy dissipation (Mendez & Losada 2004)
  │   └── Wave thrust on plants
  └── Sediment
      ├── Reduced bed shear stress under canopy
      ├── Enhanced settling in vegetation patches
      └── Modified erosion/deposition patterns
```

### CPP Options

```c
/* Vegetation master switches */
#define VEG_DRAG          /* Vegetation drag on currents */
#define VEG_FLEX          /* Flexible vegetation (bending) */
#define VEG_TURB          /* Vegetation-enhanced turbulence */
#define VEG_SWAN_COUPLING /* Wave-vegetation coupling via SWAN */
#define VEG_STREAMING     /* Wave streaming effect */

/* Marsh dynamics */
#define MARSH_DYNAMICS    /* Lateral marsh erosion/accretion */
#define MARSH_SED_EROSION /* Sediment-marsh interaction */
#define MARSH_WAVE_THRUST /* Wave thrust on marsh edge */
#define MARSH_TIDAL_RANGE /* Tidal range control on marsh */
```

### Vegetation Input (`vegetation.in`)

```
! Number of vegetation types
  NVEG = 1

! Vegetation properties per type (can vary spatially)
  VEG_DENS = 500.0       ! stem density (stems/m²)
  VEG_HGHT = 0.3         ! stem height (m)
  VEG_DIAM = 0.006       ! stem diameter (m)
  VEG_CD   = 1.0         ! drag coefficient (dimensionless)
  VEG_FLEX_MODULUS = 0.5  ! Young's modulus for flexibility (GPa)
  VEG_THICK = 0.001      ! blade thickness (m) — for flex
```

### Drag Force Formulation

Vegetation drag is modeled as a body force in the momentum equations:

```
F_drag = -0.5 * Cd * a * |u| * u

Where:
  Cd = drag coefficient (typically 1.0-1.5)
  a  = frontal area per unit volume = N * d * h_eff / h
  N  = stem density (stems/m²)
  d  = stem diameter (m)
  h_eff = effective vegetation height (after bending)
  |u| = current speed magnitude
```

### Flexible Vegetation

When `VEG_FLEX` is defined, plants bend under combined current and wave forcing:

```
Effective height: h_eff = h * f(Ca)
  Where Ca = Cauchy number = ρ * Cd * a * u² * h³ / (E * I)
  E = Young's modulus
  I = second moment of area

Large Ca → plant bends flat → small effective height
Small Ca → plant stands upright → full height
```

### Wave-Vegetation Coupling (SWAN)

With `VEG_SWAN_COUPLING`, ROMS sends vegetation properties to SWAN:
- SWAN calculates wave dissipation by vegetation using Mendez & Losada (2004)
- Dissipation depends on stem density, diameter, height, and drag coefficient
- SWAN returns reduced wave energy to ROMS

### Marsh Dynamics

`MARSH_DYNAMICS` enables lateral marsh erosion and accretion:

```
Marsh erosion rate:
  dB/dt = -Be * (P_wave / L_marsh)

Where:
  Be = marsh erodibility coefficient
  P_wave = wave power at marsh edge
  L_marsh = marsh width

Marsh accretion:
  Vertical accretion from sediment deposition during inundation
  Controlled by tidal range and sediment supply
```

### Example: SAV in an Estuary

```c
/* Header file for SAV application */
#define VEG_DRAG
#define VEG_FLEX
#define VEG_TURB
#define VEG_SWAN_COUPLING
#define SEDIMENT
```

```
! vegetation.in
NVEG = 2
! Type 1: Seagrass (Zostera marina)
  VEG_DENS(1) = 800.0    ! dense meadow
  VEG_HGHT(1) = 0.4
  VEG_DIAM(1) = 0.005
  VEG_CD(1)   = 1.2
! Type 2: Kelp (Macrocystis)
  VEG_DENS(2) = 10.0     ! sparse canopy
  VEG_HGHT(2) = 5.0
  VEG_DIAM(2) = 0.03
  VEG_CD(2)   = 1.5
```

---

## 27. InWave (Infragravity Waves)

> Based on COAWST 2021 Training Day 4 — InWave (Maitane Olabarrieta)

### Overview

**InWave** is an infragravity (IG) wave driver within COAWST that resolves low-frequency wave motions (periods 25s–300s) generated by the nonlinear interaction of short-wave groups.

### Why InWave?

Standard wave models (SWAN, WW3) are **phase-averaged** — they don't resolve individual waves. This misses:
- **Infragravity waves** (25-300s period) from wave group forcing
- **Wave runup** on beaches
- **Harbor resonance** at IG frequencies
- **Surf beat** and edge waves

InWave fills this gap using a **phase-averaged approach for IG waves** (not phase-resolving like Boussinesq models), making it computationally efficient.

### InWave Physics

```
InWave solves:
  ├── Wave action conservation (for IG band)
  ├── Radiation stress from short waves (via SWAN)
  │   └── Short-wave group envelope → IG wave forcing
  ├── IG wave propagation (refraction, shoaling)
  └── IG wave dissipation (breaking, bottom friction)
```

The short-wave envelope (from SWAN) forces IG wave generation through:
1. **Bound IG waves**: Phase-locked to wave groups (shoaling amplification)
2. **Free IG waves**: Released at the breakpoint, propagate independently
3. **Edge waves**: Trapped along the coast by refraction

### CPP Options

```c
#define INWAVE_MODEL       /* Master InWave switch */
#define INWAVE_SWAN_COUPLING  /* Get short-wave forcing from SWAN */

/* InWave physics options */
#define INWAVE_ROLLER      /* Surface roller model */
#define INWAVE_BOTTOM_FRICTION  /* Bottom friction for IG waves */
#define INWAVE_BREAKING    /* IG wave breaking */
```

### InWave Input Configuration

In `inwave.in`:
```
! Spectral discretization for IG band
  IG_NFREQ = 10           ! number of IG frequency bins
  IG_FMIN  = 0.004        ! min IG frequency (Hz) = 250s period
  IG_FMAX  = 0.04         ! max IG frequency (Hz) = 25s period
  IG_NDIR  = 36           ! directional bins

! Breaking parameters
  IG_GAMMA = 0.5          ! breaking ratio for IG waves
  IG_ALPHA = 1.0          ! breaking dissipation coefficient
```

### Coupling with SWAN and ROMS

```
SWAN → short-wave spectra → InWave (IG generation)
                                 ↓
                            IG wave spectra
                                 ↓
                         IG radiation stresses → ROMS (currents, setup)
                                 ↓
                         ROMS → water levels, currents → InWave (refraction)
```

### When to Use InWave

| Application | Use InWave? |
|-------------|------------|
| Open ocean wave forecasting | No — use SWAN/WW3 |
| Coastal engineering (breakwaters, harbors) | Yes — IG resonance matters |
| Beach morphology & runup | Yes — IG drives swash |
| Storm surge with wave setup | Maybe — if IG contribution significant |
| Coral reef hydrodynamics | Yes — IG dominates inside lagoon |

---

## 28. Biology & Biogeochemistry

> Based on COAWST 2021 Training Day 4 — Biology (Neil Ganju, Tarandeep Kalra)

### Overview

ROMS includes multiple biogeochemical (BGC) modules that can be activated through CPP flags. These range from simple nutrient-phytoplankton models to complex multi-element cycles.

### Available Biology Modules

| Module | CPP Flag | Tracers | Description |
|--------|----------|---------|-------------|
| **Fennel** | `BIO_FENNEL` | 7+ | N-based, 1 phyto, 1 zoo (most popular) |
| **NPZD (Powell)** | `NPZD_POWELL` | 4 | Simple N-P-Z-D |
| **NPZD (Franks)** | `NPZD_FRANKS` | 4 | Alternative NPZD |
| **EcoSim** | `ECOSIM` | 30+ | Multi-species ecosystem |
| **Bio-Sediment** | `BEST_NPZ` | 15+ | Bering Sea specific |
| **Estuarine** | `ESTUARINE_BGC` | varies | Estuarine biogeochemistry |
| **Red Tide** | `RED_TIDE` | varies | Harmful algal bloom |
| **Hypoxia** | `HYPOXIA_SRM` | varies | Dissolved oxygen |
| **Carbon** | `CARBON` | +2 | Add-on: DIC + alkalinity |
| **Oxygen** | `OXYGEN` | +1 | Add-on: dissolved O₂ |
| **PCB** | `PCB_BATHY` | varies | Contaminant tracking |
| **Reef Ecosys.** | `REEF_ECOSYS` | varies | Coral reef ecosystem |

### Fennel Model (Most Common)

The Fennel et al. (2006) model is the most widely used ROMS biology module:

```
State Variables:
  NO₃  - Nitrate
  NH₄  - Ammonium
  Phyt - Phytoplankton (chlorophyll-based)
  Zoo  - Zooplankton
  LDet - Large detritus (sinking)
  SDet - Small detritus (suspended)
  Oxy  - Dissolved oxygen (if OXYGEN defined)
  DIC  - Dissolved inorganic carbon (if CARBON defined)
  Alk  - Alkalinity (if CARBON defined)
```

```
Flow diagram:
  NO₃ → Phyto → Zoo → LDet → NH₄ → NO₃ (nitrification)
         ↓        ↓      ↓
        SDet    SDet   sinking
         ↓        ↓
        NH₄     NH₄ (remineralization)
```

### CPP Options for Biology

```c
/* Master biology switch */
#define BIO_FENNEL         /* Fennel nitrogen-based model */

/* Additional components */
#define CARBON             /* Add DIC + alkalinity tracers */
#define OXYGEN             /* Add dissolved oxygen tracer */
#define DENITRIFICATION    /* Include denitrification process */
#define BIO_SEDIMENT       /* Biology-sediment coupling */
#define DIAGNOSTICS_BIO    /* Output biological diagnostics */

/* Light */
#define SHORTWAVE          /* Shortwave radiation penetration */
#define ANA_SPFLUX         /* Analytical surface tracer flux */
#define ANA_BPFLUX         /* Analytical bottom tracer flux */
```

### Biology Input (`bio_Fennel.in`)

```
! Light parameters
  AttSW = 0.04        ! seawater light attenuation (1/m)
  AttChl = 0.024      ! chlorophyll light attenuation (1/(mg Chl m²))
  PARfrac = 0.43      ! fraction of shortwave that is PAR

! Phytoplankton
  Vp0 = 1.5           ! max growth rate at 0°C (1/day)
  I_thNH4 = 1.5       ! NH4 inhibition of NO3 uptake (mmol/m³)
  D_p5NH4 = 0.1       ! half-saturation for NH4 uptake (mmol/m³)
  D_p5NO3 = 0.5       ! half-saturation for NO3 uptake (mmol/m³)

! Zooplankton
  ZooGR = 0.6         ! grazing rate (1/day)
  ZooEEN = 0.3        ! excretion efficiency
  ZooEED = 0.2        ! egestion efficiency

! Remineralization
  LDeRRN = 0.01       ! large detritus remineralization (1/day)
  SDeRRN = 0.03       ! small detritus remineralization (1/day)
  CoagR = 0.005       ! coagulation rate (1/day m³/mmol N)
```

### HydroBioSed Coupling

When biology and sediment are both active, **HydroBioSed** coupling allows:
- Organic matter settling and burial in sediment
- Nutrient release from sediment to water column
- Sediment-mediated denitrification
- Resuspension of buried organic material

```c
#define BIO_FENNEL
#define SEDIMENT
#define BIO_SEDIMENT        /* Enable coupling */
```

### Tips for Running Biology

1. **Start simple**: Get the physics right first (no biology), then add biology
2. **Initial conditions**: Biology needs realistic initial nutrient/chlorophyll profiles (from WOA, GLORYS-BGC, etc.)
3. **Spin-up**: Biology needs longer spin-up than physics (months to years for deep ocean)
4. **Diagnostics**: Always enable `DIAGNOSTICS_BIO` to track biological rates
5. **Time step**: Biology uses the same time step as ROMS — if biology is unstable, reduce `DT`

---

## 29. COAWST MATLAB Tools Reference

> Based on COAWST 2021 Training Day 4 — Tools and Cloud Access

### Overview

COAWST includes MATLAB tools in the `Tools/` directory for pre-processing, grid generation, and post-processing. Additionally, the community is moving toward **Python (Pangeo)** workflows.

### MATLAB Grid Tools (`Tools/mfiles/`)

| Tool | Purpose |
|------|---------|
| `create_roms_grid.m` | Interactive ROMS grid generation |
| `create_swan_grid.m` | Generate SWAN grid from ROMS grid |
| `create_scrip_weights.m` | Generate SCRIP interpolation weights |
| `add_bathy.m` | Add bathymetry to ROMS grid |
| `add_mask.m` | Edit land/sea mask interactively |
| `coawst_explorer.m` | GUI for exploring COAWST output |

### COAWST Explorer

A MATLAB GUI for visualizing COAWST output:

```matlab
% Launch the explorer
coawst_explorer

% Features:
% - Load multiple NetCDF output files (ROMS, SWAN, WRF)
% - 2D/3D visualization of any variable
% - Time animation
% - Cross-section plots
% - Vector overlay (currents, winds)
% - Difference plots (two runs)
```

### Python/Pangeo Alternatives

The COAWST community increasingly uses Python tools:

| Python Tool | Replaces MATLAB... | Purpose |
|-------------|-------------------|---------|
| `xarray` | NetCDF MATLAB tools | Read/analyze NetCDF output |
| `xgcm` | Grid-aware analysis | Grid-aware computations (ROMS C-grid) |
| `xroms` | ROMS MATLAB tools | ROMS-specific xarray extensions |
| `cartopy` | m_map | Map projections and coastlines |
| `holoviews`/`hvplot` | Visualization | Interactive visualization |
| `dask` | N/A | Parallel/out-of-core computation |
| `pyroms` | ROMS MATLAB tools | ROMS grid gen + pre/post-processing |
| `gridgen-c` | create_roms_grid | Grid generation (C library) |

### Cloud Access for COAWST

Running COAWST on cloud platforms:

```
Cloud workflow:
  1. Launch compute instance (AWS EC2, Google Cloud, Azure)
  2. Install dependencies (compilers, MPI, NetCDF, HDF5)
  3. Clone COAWST
  4. Compile and run
  5. Analyze output via Jupyter
```

**Jupyter Remote Access:**
```bash
# On remote server:
jupyter lab --no-browser --port=8888

# On local machine (SSH tunnel):
ssh -N -L 8888:localhost:8888 user@remote-server

# Open in browser: http://localhost:8888
```

### Key Data Analysis Patterns (Python)

```python
import xarray as xr

# Open ROMS output
ds = xr.open_dataset('ocean_his.nc')

# Plot sea surface temperature
ds.temp.isel(s_rho=-1, ocean_time=0).plot()

# Time series at a point
ds.temp.isel(s_rho=-1, eta_rho=50, xi_rho=80).plot()

# Cross-section
ds.temp.isel(ocean_time=0, xi_rho=80).plot(x='eta_rho', y='s_rho')

# Compute depth-averaged current
import xroms
ds_roms = xroms.open_netcdf('ocean_his.nc')
u_avg = ds_roms.u.mean(dim='s_rho')
```

### Pre-Processing Checklist

```
□ Grid: create_roms_grid.m or pyroms
  □ Bathymetry smoothed (rx0 < 0.35)
  □ Land mask edited
  □ Minimum depth set (h_min = 5m typical)
□ SCRIP weights generated (for coupling)
□ Initial conditions: from GLORYS, HYCOM, or analytical
□ Boundary conditions: from same source, time-varying
□ Forcing files: atmospheric (ERA5/GFS), tidal (TPXO/FES)
□ SWAN grid + bottom + boundary
□ WRF inputs (wrfinput, wrfbdy) via WPS
□ All grids registered in coupling_coawst.in
```

---

## References

### Foundational Papers

- Warner, J.C., Armstrong, B., He, R., and Zambon, J.B. (2010): "Development of a Coupled Ocean-Atmosphere-Wave-Sediment Transport (COAWST) modeling system." *Ocean Modelling*, 35(3), 230–244
- Kumar, N., Voulgaris, G., Warner, J.C., and Olabarrieta, M. (2012): "Implementation of the vortex force formalism in COAWST." *Ocean Modelling*, 47, 65–95
- Booij, N., Ris, R.C., and Holthuijsen, L.H. (1999): "A third-generation wave model for coastal regions." *J. Geophys. Res.*, 104(C4), 7649–7666

### Sediment Transport

- Warner, J.C., Sherwood, C.R., Signell, R.P., Harris, C.K., and Arango, H.G. (2008): "Development of a three-dimensional, regional, coupled wave, current, and sediment-transport model." *Computers & Geosciences*, 34(10), 1284–1306

### Online Resources

- [COAWST GitHub](https://github.com/DOI-USGS/COAWST)
- [COAWST USGS Page](https://www.usgs.gov/programs/coastal-and-marine-hazards-and-resources-program/science/coupled-ocean-atmosphere-wave)
- [COAWST Publications (112+ papers)](https://github.com/DOI-USGS/COAWST/wiki/Publications)
- [COAWST User's Guide (3rd Edition)](https://www.researchgate.net/publication/345813288_COAWST_User's_Guide_-_Third_Edition)
- [SWAN Official Website](https://swanmodel.sourceforge.io/)
- [SWAN User Manual (PDF)](https://swanmodel.sourceforge.io/download/zip/swanuse.pdf)
- [SWAN Technical Documentation (PDF)](https://swanmodel.sourceforge.io/download/zip/swantech.pdf)
- [WikiROMS Sediment Model](https://www.myroms.org/wiki/Sediment_Model)
- [WikiROMS sediment.in](https://www.myroms.org/wiki/sediment.in)
- [Baron Ina-CAWO Whitepaper (PDF)](https://21175713.fs1.hubspotusercontent-na1.net/hubfs/21175713/Brochures%20and%20eBooks/Whitepapers/Baron_BMKG-CAWO%20Whitepaper_06-2024.pdf)
- [BMKG Ina-CAWO Web Interface](https://web-meteo.bmkg.go.id/id/model-prediksi-cuaca/inacawo)

### Data Sources

- [NCEP NOMADS (GFS)](https://nomads.ncep.noaa.gov/)
- [Copernicus CDS (ERA5)](https://cds.climate.copernicus.eu/)
- [Copernicus Marine (GLORYS12V1)](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030)
- [TPXO Tidal Models](https://www.tpxo.net/)
- [AVISO FES2022 Tides](https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html)
- [GEBCO Bathymetry](https://download.gebco.net/)
- [WAVEWATCH III Hindcasts](https://polar.ncep.noaa.gov/waves/hindcasts/nopp-phase2.php)
