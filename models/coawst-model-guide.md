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

## References

### Foundational Papers

- Warner, J.C., Armstrong, B., He, R., and Zambon, J.B. (2010): "Development of a Coupled Ocean-Atmosphere-Wave-Sediment Transport (COAWST) modeling system." *Ocean Modelling*, 35(3), 230–244
- Kumar, N., Voulgaris, G., Warner, J.C., and Olabarrieta, M. (2012): "Implementation of the vortex force formalism in COAWST." *Ocean Modelling*, 47, 65–95
- Booij, N., Ris, R.C., and Holthuijsen, L.H. (1999): "A third-generation wave model for coastal regions." *J. Geophys. Res.*, 104(C4), 7649–7666

### Online Resources

- [COAWST GitHub](https://github.com/DOI-USGS/COAWST)
- [COAWST USGS Page](https://www.usgs.gov/programs/coastal-and-marine-hazards-and-resources-program/science/coupled-ocean-atmosphere-wave)
- [COAWST Publications (112+ papers)](https://github.com/DOI-USGS/COAWST/wiki/Publications)
- [SWAN Official Website](https://swanmodel.sourceforge.io/)
- [SWAN User Manual (PDF)](https://swanmodel.sourceforge.io/download/zip/swanuse.pdf)
- [SWAN Technical Documentation (PDF)](https://swanmodel.sourceforge.io/download/zip/swantech.pdf)
- [Baron Ina-CAWO Whitepaper (PDF)](https://21175713.fs1.hubspotusercontent-na1.net/hubfs/21175713/Brochures%20and%20eBooks/Whitepapers/Baron_BMKG-CAWO%20Whitepaper_06-2024.pdf)
- [BMKG Ina-CAWO Web Interface](https://web-meteo.bmkg.go.id/id/model-prediksi-cuaca/inacawo)
