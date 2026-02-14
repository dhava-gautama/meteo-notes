# WAVEWATCH III: Complete Guide from Grid Setup to Visualization

> A practical guide covering the full WAVEWATCH III (WW3) workflow for global, regional, and coastal wave modeling.
> Covers governing equations, switch file system, all source term packages, grid types, multi-grid mosaic, input/output, nesting, operational use, coupling, and calibration.

---

## Table of Contents

1. [Overview & Key Concepts](#1-overview--key-concepts)
2. [Installation & Compilation](#2-installation--compilation)
3. [Grid Types](#3-grid-types)
4. [Physics Options (All Source Terms)](#4-physics-options-all-source-terms)
5. [Input Files](#5-input-files)
6. [Multi-Grid (Mosaic) System](#6-multi-grid-mosaic-system)
7. [Running WW3](#7-running-ww3)
8. [Output & Post-Processing](#8-output--post-processing)
9. [Nesting & Boundary Conditions](#9-nesting--boundary-conditions)
10. [Operational Use](#10-operational-use)
11. [Coupling with Other Models](#11-coupling-with-other-models)
12. [Calibration & Validation](#12-calibration--validation)

---

## 1. Overview & Key Concepts

### What Is WAVEWATCH III?

WAVEWATCH III (WW3) is a third-generation, phase-averaged, spectral wind-wave model developed and maintained at NOAA/NCEP Environmental Modeling Center (EMC). The model framework was originally conceived by Hendrik Tolman and has evolved into a community-driven wave modeling framework with 50+ contributing scientists spanning most continents.

**Current stable release:** Version 6.07.1 (April 2024). The development branch (v7.14) includes CMake builds, C preprocessor directives, and `.F90` source format with NetCDF4 as default I/O.

**License:** GNU Lesser General Public License v3 (LGPL v3)

**GitHub repository:** https://github.com/NOAA-EMC/WW3

**Code language:** ANSI standard Fortran 90, fully modular and fully allocatable.

**Documentation:**
- Official NOAA model page: https://polar.ncep.noaa.gov/waves/wavewatch/
- GitHub Wiki: https://github.com/NOAA-EMC/WW3/wiki/
- Online Tutorial (v6.07): https://ww3-docs.readthedocs.io/en/latest/
- User Manual v6.07 PDF: https://raw.githubusercontent.com/wiki/NOAA-EMC/WW3/files/manual.pdf

### Governing Equation

WW3 solves the random phase spectral action density balance equation for wavenumber-direction spectra. The action density `N(k, theta)` is defined as:

```
N(k, theta; x, t) = F(k, theta; x, t) / sigma
```

where `F` is the energy density spectrum and `sigma` is the intrinsic (relative) frequency. Action density is conserved in the presence of currents.

**General form:**

```
dN/dt + nabla_x . [(c_g + U) N] + d/dk [k_dot N] + d/dtheta [theta_dot N] = S_tot / sigma
```

| Term | Meaning |
|------|---------|
| `dN/dt` | Local rate of change of action density |
| `nabla_x . [(c_g + U) N]` | Propagation in geographic space (c_g = group velocity, U = current) |
| `d/dk [k_dot N]` | Shifting in wavenumber space (shoaling, frequency shifting) |
| `d/dtheta [theta_dot N]` | Refraction (depth and current-induced directional turning) |
| `S_tot / sigma` | All source and sink terms |

**Total source function:**

```
S_tot = S_ln + S_in + S_nl + S_ds + S_bot + S_db + S_tr + S_bs + S_ice + S_ref + S_xx
```

| Symbol | Process |
|--------|---------|
| `S_ln` | Linear (initial) wind input |
| `S_in` | Exponential wind input (wind-wave growth) |
| `S_nl` | Nonlinear quadruplet wave-wave interactions |
| `S_ds` | Whitecapping dissipation |
| `S_bot` | Bottom friction |
| `S_db` | Depth-induced breaking |
| `S_tr` | Triad wave-wave interactions |
| `S_bs` | Bottom scattering |
| `S_ice` | Wave-ice interactions |
| `S_ref` | Wave reflection |
| `S_xx` | User-defined / additional source terms |

### Phase-Averaged Spectral Model

WW3 is a phase-averaged model. The sea state is described statistically by a 2D spectrum (frequency vs. direction), discretized using:

- **Constant directional increment** covering all directions. Typically 24 directions (15-degree spacing) or 36 directions (10-degree spacing).
- **Logarithmic intrinsic frequency grid**: `f_n = f_1 * r^(n-1)` where `r` is typically 1.10 (10% increment). Typically 25-36 frequency bins with lowest frequency ~0.0345-0.0418 Hz.

### Differences from SWAN and WAM

| Feature | WW3 | SWAN | WAM |
|---------|-----|------|-----|
| Primary scale | Global to regional | Nearshore, coastal | Deep ocean, global |
| Propagation scheme | 3rd-order ULTIMATE QUICKEST (PR3/UQ) | Implicit (unconditionally stable) | 1st-order explicit upwind |
| Grid types | Regular, curvilinear, unstructured, SMC, tripole, mosaic | Regular, curvilinear, unstructured | Regular |
| Multi-grid | Full two-way mosaic (ww3_multi) | One-way nesting | One-way nesting |
| Parallelization | MPI (spectral decomp), OpenMP, hybrid | MPI, OpenMP | MPI |
| Source term flexibility | Multiple packages (ST1-ST6) | Multiple (KOMEN, JANSSEN, WESTH, ST6) | Single set |
| Time stepping | Dynamic, explicit (structured); implicit (unstructured) | Implicit | Explicit |
| Developer | NOAA/NCEP EMC + community | TU Delft | ECMWF |
| Operational users | NOAA, Met Office, BoM, Meteo-France | Coastal agencies, Deltares, USACE | ECMWF, DWD |

### Version History (Major Releases)

| Version | Year | Key Features |
|---------|------|-------------|
| v1.18 | 1999 | Public domain; regular grids only |
| v2.22 | 2002 | Improved source terms; GSE alleviation |
| v3.14 | 2008 | Mosaic grids (ww3_multi); moving grid for hurricanes |
| v4.18 | 2013 | Curvilinear grids; unstructured grids; ST4 source terms; SMC grids |
| v5.16 | 2016 | Community-driven; hybrid MPI-OpenMP; tripole grids; OASIS coupling |
| v6.07 | 2024 | Open GitHub development; PDLIB domain decomposition; ESMF interface; IC5 ice |
| v7.14 | dev | CMake build; CPP directives; `.F90` format; NetCDF4 default |

---

## 2. Installation & Compilation

### Obtaining the Source Code

```bash
# Fork then clone:
git clone --recurse-submodules https://github.com/<your-username>/WW3.git
cd WW3

# Download required binary/NetCDF data for regression tests:
sh model/bin/ww3_from_ftp.sh
```

### Dependencies

**Required:**
- Fortran 90 compiler: gfortran (GNU >= 9), ifort/ifx (Intel oneAPI), pgfortran (NVIDIA HPC SDK)
- C compiler: gcc, icc, or compatible
- Make: GNU make (legacy) or CMake >= 3.19 (modern)

**Required for NetCDF output (strongly recommended):**

Build order: zlib --> HDF5 --> NetCDF-C --> NetCDF-Fortran. The `NC4` switch must be included in the switch file.

**Required for parallel execution:**
- MPI library: MPICH, OpenMPI, Intel MPI, or Cray MPI
- Switch file must include `DIST MPI` (instead of `SHRD`)

**Optional:**
- ParMETIS or SCOTCH: Domain decomposition on unstructured grids (`PDLIB` switch)
- ESMF/NUOPC: Coupled model applications (UFS, CESM)
- OASIS-MCT: OASIS coupling framework

### Legacy Build System (w3_setup / w3_make)

**Step 1 -- Environment initialization:**

```bash
cd WW3/model
./bin/w3_setup -c <compiler> -s <switch_file> <install_path>
```

**Step 2 -- Edit compiler/linker scripts:**

Located in `WW3/model/bin/`:
- `comp.gfortran`, `comp.Intel`, etc. -- set compiler, optimization flags, include paths
- `link.gfortran`, `link.Intel`, etc. -- set library paths (`-lnetcdf -lnetcdff`)

Users MUST edit these to match their local system.

**Step 3 -- Compilation:**

```bash
# Compile all programs:
./bin/w3_make

# Compile specific programs:
./bin/w3_make ww3_grid ww3_strt ww3_shel ww3_ounf
```

Executables are placed in `WW3/model/exe/`.

### Modern Build System (CMake)

```bash
cd WW3
mkdir build && cd build
cmake .. -DSWITCH=<switch_name>
make -j8
```

Executables placed in `WW3/build/bin/`. Requires CMake >= 3.19.

### Switch File System

The switch file is a single line of space-separated keywords controlling compile-time options. Located in `WW3/model/bin/` with names like `switch_<config>`.

**Example switch file:**
```
F90 DIST MPI LRB4 SEED ST4 STAB0 FLX0 LN1 NL1 BT4 DB1 TR0 BS0 IC0 IS0 REF0 WNT1 WNX1 CRT1 CRX1 PR3 UQ NC4 NOGRB O0 O1 O2
```

**Switch groups (complete reference):**

| Group | Options | Description |
|-------|---------|-------------|
| Fortran standard | `F90` | ANSI Fortran 90 |
| Hardware/messaging | `SHRD`, `DIST`, `MPI` | Serial/shared or distributed memory |
| Propagation | `PR0`, `PR1`, `PR2/UNO`, `PR3/UQ` | None, 1st, 2nd, 3rd order |
| Flux computation | `FLX0`-`FLX4`, `FLXX` | Friction velocity methods |
| Linear input | `LN0`, `SEED`, `LN1`, `LNX` | Initial wind input |
| Wind/whitecapping | `ST0`-`ST4`, `ST6`, `STX` | Source term packages |
| Nonlinear interactions | `NL0`-`NL5`, `NLX` | DIA, WRT, GMD, TSA, GKE |
| Bottom friction | `BT0`, `BT1`, `BT4`, `BT8`, `BT9`, `BTX` | Various formulations |
| Sea ice damping | `IC0`-`IC5` | None to generalized polynomial |
| Ice scattering | `IS0`-`IS2` | None to floe-size dependent |
| Reflection | `REF0`, `REF1` | None or shoreline/iceberg |
| Depth breaking | `DB0`, `DB1`, `DBX` | None or Battjes-Janssen |
| Triad interactions | `TR0`, `TR1`, `TRX` | None or LTA |
| Bottom scattering | `BS0`, `BS1`, `BSX` | None or Magne-Ardhuin |
| Wind interpolation | `WNT0`-`WNT2`, `WNX0`-`WNX2` | Time and space |
| Current interpolation | `CRT0`-`CRT2`, `CRX0`-`CRX2` | Time and space |
| GRIB output | `NOGRB`, `NCEP1`, `NCEP2` | None, GRIB1, GRIB2 |
| Additional | `NC4`, `LRB4`/`LRB8`, `PDLIB`, `SCRIP`, `TIDE`, `SETUP` | Various |

**Common switch file examples:**

Deep ocean operational (ST4, MPI):
```
F90 DIST MPI LRB4 SEED ST4 STAB0 FLX0 LN1 NL1 BT1 DB0 TR0 BS0 IC0 IS0 REF0 WNT1 WNX1 CRT1 CRX1 PR3 UQ NC4 NOGRB
```

Coastal application (with depth breaking and triads):
```
F90 DIST MPI LRB4 SEED ST4 STAB0 FLX0 LN1 NL1 BT4 DB1 TR1 BS1 IC0 IS0 REF1 WNT1 WNX1 CRT1 CRX1 PR3 UQ NC4 NOGRB
```

Arctic with ice:
```
F90 DIST MPI LRB4 SEED ST4 STAB0 FLX0 LN1 NL1 BT1 DB0 TR0 BS0 IC3 IS2 REF1 WNT1 WNX1 CRT1 CRX1 PR3 UQ NC4 NOGRB
```

Unstructured grid:
```
F90 DIST MPI LRB4 SEED ST4 STAB0 FLX0 LN1 NL1 BT4 DB1 TR1 BS0 IC0 IS0 REF0 PDLIB WNT1 WNX1 CRT1 CRX1 PR3 UQ NC4 NOGRB
```

### Main WW3 Programs

| Program | Purpose |
|---------|---------|
| `ww3_grid` | Grid preprocessor: reads `ww3_grid.inp`, produces `mod_def.ww3` |
| `ww3_strt` | Initial conditions: reads `ww3_strt.inp`, produces `restart.ww3` |
| `ww3_bound` | Boundary conditions: creates `nest.ww3` from coarser grid |
| `ww3_bounc` | Boundary conditions from NetCDF input |
| `ww3_prnc` | NetCDF forcing preprocessor (wind, current, ice, level) |
| `ww3_prep` | General ASCII forcing preprocessor |
| `ww3_shel` | Main model driver (single grid) |
| `ww3_multi` | Main model driver (multi-grid / mosaic) |
| `ww3_outf` | Post-processor: gridded fields to ASCII |
| `ww3_ounf` | Post-processor: gridded fields to NetCDF |
| `ww3_outp` | Post-processor: point output to ASCII |
| `ww3_ounp` | Post-processor: point output to NetCDF |
| `ww3_grib` | Post-processor: gridded fields to GRIB/GRIB2 |
| `ww3_trck` | Post-processor: track output |

---

## 3. Grid Types

### 3.1 Regular Latitude-Longitude Grid

The standard grid type. Specified in `ww3_grid.inp`:

```
'RECT' T 'NONE'
NX NY
SX SY SFAC
X0 Y0 SFAC
```

| Parameter | Description |
|-----------|-------------|
| `NX`, `NY` | Grid dimensions |
| `SX`, `SY` | Grid spacing (degrees or meters) |
| `X0`, `Y0` | Origin coordinates |
| `SFAC` | Scaling factor |
| Closure: `NONE` | Open edges |
| Closure: `SMPL` | Periodic in longitude (for global grids) |

### 3.2 Curvilinear Grid

Added in v4.01. Logically rectangular (i,j indexing) but physical coordinates of each grid point specified explicitly.

```
'CURV' T 'NONE'
NX NY
<coordinate arrays from separate files>
```

Uses: matching coastline geometry, concentrating resolution, aligning with ocean model grids (MOM6, ROMS).

### 3.3 Unstructured Grid (Triangular, PDLIB)

Added in v4.02. Provides highly variable resolution (50 km deep ocean to 100 m nearshore).

```
'UNST' T 'NONE'
```

- Requires `PDLIB` switch for MPI parallelization via ParMETIS/SCOTCH
- Supports implicit numerical scheme (no CFL constraint at fine resolutions)
- Compatible with ESMF/NUOPC coupling
- Used in Great Lakes operational model and E3SM

### 3.4 Spherical Multi-Cell (SMC) Grid

Added in v4.13. Variable resolution while retaining quadrilateral cells.

- Quad-tree-like mesh refinement: cell sizes vary by powers of 2
- Static refinement (not adaptive during runtime)
- Sub-timesteps for refined cells
- Avoids polar singularity through rotated grid options
- Operational use: UK Met Office since 2016

### 3.5 Tripole Grid

Added in v5.03. Avoids North Pole singularity for global ocean configurations. Supports standard tripole definitions as used by GFDL/MOM6.

---

## 4. Physics Options (All Source Terms)

### 4.1 Wind Input & Whitecapping (STn Packages)

Wind input and whitecapping are treated together because they are physically coupled through the energy balance.

#### ST1 -- WAM Cycle 3

- References: Snyder et al. (1981), Komen et al. (1984)
- Simple and computationally efficient
- Known issue: insufficient swell dissipation
- Namelists: `&SIN1`, `&SDS1`

#### ST2 -- Tolman and Chalikov (1996)

- Direct air-sea momentum flux calculation
- Includes atmospheric boundary layer stability correction (`STAB2`)
- Used in Phase 1 of NOAA 30-year hindcast
- Namelists: `&SIN2`, `&SDS2`

#### ST3 -- WAM Cycle 4+ (Bidlot et al.)

- Quasi-linear theory following Janssen (1991)
- Feedback of wave-supported stress on atmospheric boundary layer
- Basis of ECMWF operational WAM physics
- More sensitive to domain characteristics than ST4/ST6
- Namelists: `&SIN3`, `&SDS3`, `&FLX3`

#### ST4 -- Ardhuin et al. (2010) -- RECOMMENDED

- Saturation-based whitecapping with threshold behavior
- Includes swell dissipation (linear viscous + nonlinear turbulent decay)
- **Key tuning parameters:**

| Parameter | Description | Default | Notes |
|-----------|-------------|---------|-------|
| `BETAMAX` | Overall wind input growth rate | 1.52 | 1.33 for CFSR winds, 1.39 for Met Office UM |
| `SDSBR` | Breaking dissipation coefficient | - | Controls breaking intensity |
| `SDSC1` | Saturation-based whitecapping coefficient | - | |
| `SDSBT` | Saturation threshold for breaking | - | |

- Performance: Correlation ~0.94 for Hs, ~0.83 for wave period
- Used in Phase 2 of NOAA 30-year hindcast and most modern applications
- Namelists: `&SIN4`, `&SDS4`

#### ST6 -- BYDRZ / Observation-Based

- References: Donelan et al. (2006), Babanin et al. (2007), Rogers et al. (2012)
- Wind input: Non-linear dependence on wave steepness; accounts for airflow separation at high winds
- Whitecapping: Two mechanisms -- inherent (threshold) breaking + cumulative dissipation
- Recommended for extreme wave simulations (hurricanes, severe storms)
- Namelists: `&SIN6`, `&SDS6`, `&SWL6`

### 4.2 Nonlinear Wave-Wave Interactions (NLn)

| Switch | Method | Cost | Accuracy | Usage |
|--------|--------|------|----------|-------|
| `NL1` | DIA (Hasselmann et al., 1985) | Fast | Adequate | **All operational forecasting** |
| `NL2` | WRT exact interaction | 10-100x DIA | Most accurate | Research only |
| `NL3` | GMD (Tolman, 2013) | Moderate | Better than DIA, especially shallow water | Research/advanced |
| `NL4` | TSA (Resio et al., 2013) | Moderate | Improved over DIA | Research/advanced |
| `NL5` | GKE | Moderate | Recent development | Research |

Key NL1 parameters: `LAMBDA` (coupling coefficient, default 0.25), `NLPROP` (proportionality constant).

### 4.3 Bottom Friction (BTn)

| Switch | Formulation | Key Parameter | Application |
|--------|-------------|---------------|-------------|
| `BT1` | JONSWAP (Hasselmann et al., 1973) | `GAMMA_BT` = 0.038 m^2/s^3 | Sandy bottoms, general |
| `BT4` | SHOWEX (Ardhuin et al., 2003) | Moveable bed roughness | Modern applications with ST4 |
| `BT8` | Dalrymple and Liu | Fluid mud | Muddy coasts (Gulf of Mexico) |
| `BT9` | Ng (2000) | Viscous mud layer | Muddy environments |

### 4.4 Depth-Induced Breaking (DBn)

- `DB0`: No breaking (deep-water-only applications)
- `DB1`: Battjes-Janssen (1978). Key parameter: `GAMMA` (breaking threshold, default 0.73). Namelist: `&SDB1`

### 4.5 Triad Interactions (TRn)

- `TR0`: None (standard for deep water)
- `TR1`: Lumped Triad Approximation (LTA). Important in shallow water for energy transfer to harmonics.

### 4.6 Bottom Scattering (BSn)

- `BS0`: None (standard for most applications)
- `BS1`: Magne and Ardhuin. Scattering by bottom topography variations.

### 4.7 Wave-Ice Interactions

**Ice dissipation:**

| Switch | Method | Description |
|--------|--------|-------------|
| `IC0` | None | No ice effects |
| `IC1` | Simple uniform | Exponential attenuation, frequency-independent |
| `IC2` | Liu et al. | Under-ice friction, frequency-dependent |
| `IC3` | Wang and Shen | Viscoelastic layer, most physically comprehensive |
| `IC4` | Collins and Rogers (2017) | Frequency-dependent empirical, most flexible |
| `IC5` | Mosig et al. | Generalized polynomial, most general |

**Ice scattering:**

| Switch | Method | Description |
|--------|--------|-------------|
| `IS0` | None | No scattering |
| `IS1` | Simple diffusive | Isotropic redistribution |
| `IS2` | Floe-size dependent | Boltzmann equation + creep dissipation |

### 4.8 Reflection (REF1)

- Reference: Ardhuin and Roland (2012)
- Reflection from shorelines and icebergs based on Miche number
- Important for enclosed basins, harbors, steep beaches

### 4.9 Linear Input

- `LN0`: No linear input
- `SEED`: Spectral seeding with f^{-5} tail for numerical stability
- `LN1`: Cavaleri and Malanotte-Rizzoli (1981) with filter

---

## 5. Input Files

### 5.1 ww3_grid.inp -- Grid Definition

Processed by `ww3_grid` to produce `mod_def.ww3`. Structure:

**1. Header:**
```
'GULF OF MEXICO OPERATIONAL RUN'
```

**2. Spectral definition:**
```
freq_factor  f_min  n_freq  n_dir  dir_offset
1.1          0.0418 25      24     0.
```

**3. Model flags** (six T/F values):
```
F T T T F T
```
Controls: dry run, X-propagation, Y-propagation, refraction/k-shift, source terms only, source terms active.

**4. Time steps** (four values in seconds):
```
DT_max  DT_CFL_xy  DT_CFL_ktheta  DT_min_src
900.    950.        900.            300.
```

**5. Namelists** (physics parameter overrides):
```fortran
&SIN4 BETAMAX = 1.33 /
&SDS4 SDSBR = 0.00085 /
&SNL1 LAMBDA = 0.25, NLPROP = 2.5E7 /
&SBT4 ... /
&SDB1 BJALFA = 1.0, BJGAM = 0.73 /
&MISC CICE0 = 0.25, CICEN = 0.75, FLAGTR = 4 /
```

**6. Grid definition:**
```
'RECT' T 'NONE'
NX NY
SX SY SFAC
X0 Y0 SFAC
```

**7. Bathymetry and mask data.**

Depths are positive downward. Land points have negative depth.

**Status map values:**

| Value | Meaning |
|-------|---------|
| 0 | Excluded point |
| 1 | Sea point (active) |
| 2 | Active boundary point (for nesting) |
| 3 | Land point |

### 5.2 ww3_strt.inp -- Initial Conditions

| Type | Description |
|------|-------------|
| 0 | Calm (flat) sea -- all spectral components zero |
| 1 | Parametric JONSWAP spectrum (Hs, peak period, direction, spread) |
| 2 | User-specified spectrum from file |
| 3 | Fetch-limited growth initialization |
| 4 | User-defined full spectral input |
| 5 | Two-component spectrum (wind sea + swell) |

For cold starts, Type 0 is common. Requires 1-7 day spin-up.

### 5.3 ww3_prnc.nml -- NetCDF Forcing Preprocessor

Example for wind forcing:

```fortran
&FORCING_NML
  FORCING%TIMESTART   = '20100101 000000'
  FORCING%TIMESTOP    = '20101231 000000'
  FORCING%FIELD%WINDS = T
  FORCING%GRID%LATLON = T
/

&FILE_NML
  FILE%FILENAME  = 'wind.nc'
  FILE%LONGITUDE = 'longitude'
  FILE%LATITUDE  = 'latitude'
  FILE%VAR(1)    = 'u10'
  FILE%VAR(2)    = 'v10'
/
```

**Common issue:** Latitude values in input NetCDF must be in ascending order. If reversed, fix with:
```bash
ncpdq -h -O -a -lat input.nc output.nc
```

### 5.4 ww3_shel.inp -- Main Run Configuration

**Section 1 -- Forcing flags** (T/F/C for each forcing type):
```
T T   Winds
T F   Ice concentrations
F F   Currents
...
```

**Section 2 -- Time frame:**
```
19680606 000000          $ Start time
19680607 000000          $ End time
```

**Section 3 -- Output server mode** (0=no, 1=shared, 2=dedicated)

**Section 4 -- Gridded field output:**
```
19680606 000000 3600 19680608 000000    $ Start, stride(s), stop
```

**Section 5 -- Point output:**
```
0.0   0.0   'Point_1  '
2.0   1.0   'Point_2  '
0.0   0.0   'STOPSTRING'
```

Modern practice uses `ww3_shel.nml` namelist format (takes priority over `.inp`).

### 5.5 Binary Files Reference

| File | Creator | Consumer | Content |
|------|---------|----------|---------|
| `mod_def.ww3` | `ww3_grid` | All programs | Grid definition, config |
| `restart.ww3` | `ww3_strt`/`ww3_shel` | `ww3_shel` | Spectral data (restart) |
| `out_grd.ww3` | `ww3_shel` | `ww3_outf`/`ww3_ounf`/`ww3_grib` | Gridded field output |
| `out_pnt.ww3` | `ww3_shel` | `ww3_outp`/`ww3_ounp` | Point output |
| `nest.ww3` | `ww3_shel` | `ww3_bound` | Boundary data for nesting |
| `wind.ww3` | `ww3_prnc`/`ww3_prep` | `ww3_shel` | Wind forcing |
| `current.ww3` | `ww3_prnc`/`ww3_prep` | `ww3_shel` | Current forcing |
| `ice.ww3` | `ww3_prnc`/`ww3_prep` | `ww3_shel` | Ice concentration |
| `level.ww3` | `ww3_prnc`/`ww3_prep` | `ww3_shel` | Water level |

---

## 6. Multi-Grid (Mosaic) System

### 6.1 Overview

The `ww3_multi` program (Tolman, 2008) is a major distinguishing feature of WW3. Multiple grids of different resolutions, types, and physics packages run concurrently with full two-way interaction. This is the basis for NOAA's operational multi-grid wave forecasting system.

### 6.2 Grid Classification

**Input grids:** Provide forcing fields (wind, ice, current, water levels). Multiple input grids of different coverage can be defined.

**Computational (wave) grids:** Where the wave action equation is solved. Multiple grids can overlap with different resolutions and time steps.

**Output grids:** For interpolating and writing results. Can have different extent and resolution from computational grids.

### 6.3 How Grids Communicate

**Two-way nesting is always active:**
- Coarse --> Fine: Fine grid boundary conditions interpolated from coarse grid
- Fine --> Coarse: Coarse grid values updated by averaging fine-grid results in overlap areas

**Key features:**
- Communication through interpolation of full 2D spectra
- Each grid advances with its own time step
- Synchronization at defined communication intervals
- Grid ranks determine priority/hierarchy
- Algorithm ensures conservation of wave action across boundaries

### 6.4 ww3_multi.inp Configuration

Key sections:
- Number of wave grids, input grids, unified point output flag
- Input grid definitions (which forcing fields each provides)
- Wave model grid definitions (input grid associations, rank, group, MPI weight)
- Time parameters (start, end, communication interval)
- Output configuration per grid
- Point output locations

**MPI decomposition:** Weight values control processor distribution. Example: With 64 MPI tasks and weights 0.50/0.30/0.20 for three grids, approximately 32/19/13 tasks allocated.

---

## 7. Running WW3

### 7.1 Complete Workflow (Step by Step)

**Step 1 -- Grid preprocessing:**
```bash
# Input:  ww3_grid.inp, bathymetry file, mask file
./ww3_grid
# Output: mod_def.ww3
```

**Step 2 -- Initial conditions:**
```bash
# Input:  mod_def.ww3, ww3_strt.inp
./ww3_strt
# Output: restart.ww3
```

**Step 3 -- Force field preprocessing:**
```bash
# Input:  mod_def.ww3, ww3_prnc.nml, forcing NetCDF files
./ww3_prnc    # Run once per forcing type (wind, current, ice, level)
# Output: wind.ww3, current.ww3, ice.ww3, level.ww3
```

**Step 4 -- Boundary conditions (if nested):**
```bash
# Input:  mod_def.ww3 (child grid), spectral files from parent
./ww3_bound
# Output: nest.ww3
```

**Step 5 -- Model execution:**
```bash
# Single grid:
mpirun -np 32 ./ww3_shel

# Multi-grid:
mpirun -np 64 ./ww3_multi
```

**Step 6 -- Post-processing:**
```bash
# NetCDF field output:
./ww3_ounf
# Output: ww3.YYYYMM.nc

# NetCDF point output:
./ww3_ounp
# Output: ww3.YYYYMM_spec.nc

# GRIB output:
./ww3_grib
```

### 7.2 MPI Execution Details

**Spectral component decomposition (CD):** Default MPI mode. Spectral space distributed across processors. Maximum parallelism limited by total spectral components (`n_freq * n_dir`).

**Domain decomposition (DD):** Available for unstructured grids via `PDLIB`. Uses ParMETIS or SCOTCH for mesh partitioning.

**Hybrid MPI-OpenMP:** Combines MPI across nodes with OpenMP within nodes.

**SLURM example:**

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=36
#SBATCH --time=02:00:00

module load intel netcdf mpi
srun -n 144 ./ww3_shel
```

### 7.3 Complete Workflow Example (Single Grid)

```bash
# === Environment Setup ===
export WW3_DIR=/path/to/WW3/model
export PATH=$WW3_DIR/exe:$PATH

# === Step 1: Grid ===
cd /path/to/run_directory
cp bathymetry.dat bottom.inp
cp landmask.dat mapsta.inp
./ww3_grid                     # --> mod_def.ww3

# === Step 2: Initial Conditions ===
./ww3_strt                     # --> restart.ww3

# === Step 3: Wind Forcing ===
./ww3_prnc                     # --> wind.ww3

# === Step 4: Run Model ===
mpirun -np 32 ./ww3_shel       # --> out_grd.ww3, out_pnt.ww3

# === Step 5: Post-Processing ===
./ww3_ounf                     # --> ww3.YYYYMM.nc
./ww3_ounp                     # --> ww3.YYYYMM_spec.nc
```

---

## 8. Output & Post-Processing

### 8.1 Binary Output Format

WW3 produces internal binary files (`out_grd.ww3`, `out_pnt.ww3`) in model-specific Fortran unformatted format. These are NOT directly readable and must be converted using post-processing programs. Record length controlled by `LRB4` (4-byte) or `LRB8` (8-byte) switch.

### 8.2 Field Output Variables (Complete Reference)

**Group 1 -- Forcing/Input Fields:**

| Code | Description | Units |
|------|-------------|-------|
| `DPT` | Water depth | m |
| `CUR` | Current speed | m/s |
| `WND` | Wind speed at 10m | m/s |
| `AST` | Air-sea temperature difference | K |
| `WLV` | Water level | m |
| `ICE` | Ice concentration | fraction 0-1 |

**Group 2 -- Standard Mean Wave Parameters:**

| Code | Description | Units |
|------|-------------|-------|
| `HS` | Significant wave height | m |
| `LM` | Mean wavelength | m |
| `T02` | Mean zero-crossing period | s |
| `T0M1` | Mean energy period | s |
| `T01` | Mean wave period (1st moment) | s |
| `FP` | Peak frequency | Hz |
| `DIR` | Mean wave direction | degrees |
| `SPR` | Directional spread | degrees |
| `DP` | Peak wave direction | degrees |
| `MXH` | Maximum wave height | m |

**Group 3 -- Spectral Parameters:**

| Code | Description | Units |
|------|-------------|-------|
| `EF` | 1D frequency energy spectrum | m^2/Hz |
| `TH1M` | Mean direction of 1D spectrum | degrees |
| `STH1M` | Directional spread of 1D spectrum | degrees |
| `WN` | Wavenumber | rad/m |

**Group 4 -- Spectral Partition Parameters:**

| Code | Description | Units |
|------|-------------|-------|
| `PHS` | Partitioned Hs | m |
| `PTP` | Partitioned peak period | s |
| `PDIR` | Partitioned mean direction | degrees |
| `PWS` | Partitioned wind sea fraction | - |
| `PNR` | Number of partitions | - |

**Group 5 -- Atmosphere-Waves Coupling:**

| Code | Description | Units |
|------|-------------|-------|
| `UST` | Friction velocity | m/s |
| `CHA` | Charnock coefficient | - |
| `CGE` | Energy flux | W/m |
| `TAW` | Net air-sea stress | N/m^2 |
| `WCC` | Whitecap coverage | - |

**Group 6 -- Wave-Ocean Coupling:**

| Code | Description | Units |
|------|-------------|-------|
| `SXY` | Radiation stress tensor | N/m |
| `TWO` | Wave-to-ocean stress | N/m^2 |
| `FOC` | Wave-to-ocean energy flux | W/m^2 |
| `USS` | Surface Stokes drift | m/s |
| `TUS` | Stokes transport | m^2/s |

**Group 7 -- Wave-Bottom Interaction:**

| Code | Description | Units |
|------|-------------|-------|
| `ABR` | Near-bottom rms excursion amplitude | m |
| `UBR` | Near-bottom rms orbital velocity | m/s |

### 8.3 Post-Processing Tools

**NetCDF conversion (`ww3_ounf`):**

```bash
./ww3_ounf    # out_grd.ww3 --> ww3.YYYYMM.nc
```

Readable by: xarray (Python), CDO, NCO, ncview, Panoply, GrADS.

**GRIB conversion (`ww3_grib`):** Requires `NCEP1` or `NCEP2` switch. Readable by wgrib2, eccodes, cfgrib.

**Point output:** Full 2D spectra, 1D spectra, bulk parameters, source term diagnostics, partition data.

### 8.4 Python Post-Processing

**xarray (standard approach):**

```python
import xarray as xr
ds = xr.open_dataset("ww3_output.nc")
hs = ds['hs']
hs.isel(time=0).plot()
```

**wavespectra library** (`pip install wavespectra`):

```python
import xarray as xr
dset = xr.open_dataset("ww3_spec.nc", engine="ww3")
dset.spec.hs()      # Significant wave height
dset.spec.tp()      # Peak period
dset.spec.dm()      # Mean direction
```

Features: 15+ input formats, 60+ spectral methods, built on xarray, CF-compliant. https://wavespectra.readthedocs.io/

**Other tools:**
- WW3-tools (NOAA official): https://github.com/NOAA-EMC/WW3-tools
- pyww3: https://github.com/caiostringari/pyww3
- APPMAR (GUI for WW3 hindcast analysis): https://github.com/cemanetwork/appmar

---

## 9. Nesting & Boundary Conditions

### 9.1 One-Way Nesting (Traditional)

**Step 1 -- Parent run:**

Run coarse-resolution model. In `ww3_shel.inp`, activate Output Type 5 (boundary data). Produces `nest.ww3`.

**Step 2 -- Boundary extraction:**

On child grid, run `ww3_bound` to read parent's spectral output and interpolate to child boundary points (status=2 in grid mask).

**Step 3 -- Child run:**

Run fine-resolution model with `nest.ww3` providing boundary conditions.

This is one-way only. Up to 9 nested grids can receive boundary data from a single parent.

### 9.2 ww3_bound Configuration

```
Line 1: Interpolation method (1=nearest, 2=linear interpolation)
Line 2: Verbose output flag (T/F)
Lines 3+: Boundary spectra file names
```

`ww3_bounc` provides the same functionality but reads NetCDF-format boundary spectra.

### 9.3 Multi-Grid Approach as Alternative

The `ww3_multi` mosaic provides two-way nesting:

**Advantages:**
- Automatic two-way feedback between all grids
- No sequential separate runs needed
- More physically consistent energy exchange
- Simpler workflow (single command)

**Disadvantages:**
- Higher memory requirements
- More complex configuration
- All grids must use same compiled physics (switch file)

### 9.4 Typical Nesting Strategy

```
Global grid (0.5 degrees)
  |
  |-- Regional grid (0.25-0.1 degrees)
  |     |
  |     |-- Coastal grid (0.01-0.005 degrees, or unstructured)
  |
  |-- Another regional grid
        |
        |-- Local grid (very fine resolution)
```

Boundary data intervals: typically 1-3 hours to capture swell propagation.

---

## 10. Operational Use

### 10.1 NOAA Operational Wave Model Suite

**GFS-Wave (deterministic):**
- Four daily cycles: 00Z, 06Z, 12Z, 18Z
- Forecast: Hourly to 120h, 3-hourly to 384h
- Grid: Global multi-grid with regional nested domains
- Wind forcing: GFS 10m winds at 0.25-degree resolution
- Regional models: 126-hour forecasts

**GEFS-Wave (ensemble):**
- 30 ensemble members
- Resolution: 0.25 x 0.25 degree global
- Wind forcing: 1-hour intervals
- Forecast range: Extended to 16 days

**Great Lakes Waves:**
- Single unstructured grid
- Forecast range: 84 hours at 3-hour intervals

**Validation:** https://polar.ncep.noaa.gov/waves/validation/

### 10.2 Hindcast Products

**30-Year Hindcast (NOPP):**

| Phase | Period | Physics | Grids | Wind |
|-------|--------|---------|-------|------|
| Phase 1 | 1979-2009 | ST2 | 17 regular | CFSR 0.5deg |
| Phase 2 | 1979-2009 | ST4 | 15 regular | CFSR with bias correction |

URL: https://polar.ncep.noaa.gov/waves/hindcasts/

**IFREMER Hindcasts:**
- Databases: HOMERE (European waters), IOWAGA (global)
- URL: https://data-ww3.ifremer.fr/

---

## 11. Coupling with Other Models

### 11.1 ESMF/NUOPC Coupling

WW3 is wrapped as a NUOPC component within ESMF for standardized coupling.

**WW3 imports (from atmosphere/ocean/ice):**
- 10m wind components (U10, V10)
- Surface currents
- Ice concentration, water levels, air density

**WW3 exports (to atmosphere/ocean/ice):**
- Significant wave height, mean period, mean direction
- Surface Stokes drift
- Wave-to-ocean stress, wave-supported stress
- Radiation stress tensor, Charnock coefficient, whitecap coverage

### 11.2 UFS (Unified Forecast System) Integration

WW3 is a core component of NOAA's Unified Forecast System:
- Coupled to FV3 atmosphere and MOM6 ocean via CMEPS mediator
- Supports S2S, medium range, hurricane, and coastal applications
- UFS Coastal includes WW3 + ADCIRC/ROMS/FVCOM/SCHISM

URL: https://github.com/oceanmodeling/ufs-coastal-app

### 11.3 COAWST Coupling

COAWST (USGS) couples WRF + ROMS + SWAN/WW3 + sediment via MCT:
- 800+ registered users worldwide
- Applications: Coastal storms, tropical cyclones, sediment transport
- URL: https://www.usgs.gov/centers/whcmsc/science/coawst

### 11.4 CESM Coupling

WW3 v3.14 in CESM2:
- Coarse global grid (4 x 3.25 degrees)
- Provides Stokes drift for Langmuir turbulence parameterization
- Helps resolve shallow mixed layer bias in Southern Ocean
- Repository: https://github.com/ESCOMP/WW3-CESM

### 11.5 E3SM Integration

WW3 v6.07 with unstructured grids in E3SM:
- Global-to-coastal resolution (50 km to 100 m) in a single mesh
- Leverages PDLIB parallelization

---

## 12. Calibration & Validation

### 12.1 Key Tuning Parameters

**ST4 (most commonly tuned):**

| Parameter | Description | Default | Notes |
|-----------|-------------|---------|-------|
| `BETAMAX` | Wind input growth rate | 1.52 | 1.33 for CFSR, 1.39 for Met Office UM |
| `SDSBR` | Breaking dissipation coefficient | - | |
| `SDSC1` | Saturation-based whitecapping | - | |
| `ZALP` | Wave age tuning | - | |

**General parameters (all packages):**

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `LAMBDA` (NL1) | DIA coupling coefficient | 0.25 |
| `GAMMA_BT` (BT1) | Bottom friction coefficient | 0.038 or 0.067 m^2/s^3 |
| `GAMMA` (DB1) | Breaking threshold Hs/d | 0.73 |

**IMPORTANT:** `BETAMAX` must be tuned for the specific wind forcing product. Different reanalyses (CFSR, ERA5, GFS) have different wind biases and require different values.

### 12.2 Calibration Methodology

1. Define error function (bias, RMSE, scatter index, correlation) across variables (Hs, Tp, direction)
2. Systematically vary key parameters (especially `BETAMAX` for ST4)
3. Run model for calibration period
4. Compare against buoy and altimetry observations
5. Cross-validate against independent period/dataset
6. Test under extreme conditions

### 12.3 Common Validation Datasets

**In-situ observations:**
- NDBC buoys (~200+ moored buoys): https://www.ndbc.noaa.gov/
- CDIP (Scripps, California/Pacific): https://cdip.ucsd.edu/
- OWS, PIRATA/TAO/TRITON arrays

**Satellite altimetry:**
- Jason-3, Sentinel-3A/3B, Sentinel-6, CFOSAT, SARAL/AltiKa

### 12.4 Typical Model Performance

| Metric | Open Ocean | Coastal |
|--------|-----------|---------|
| Hs bias | +/- 5-15 cm | Variable |
| Hs scatter index | 10-20% | 15-30% |
| Hs correlation | 0.93-0.97 | Variable |
| Tp correlation | 0.80-0.90 | Variable |

**Package comparison:**
- ST4: Best overall performance (R=0.94 for Hs)
- ST6: Excellent for extreme events (SI within 19%)
- ST3: Good performance, more sensitive to domain tuning

### 12.5 Recommendations for New Users

1. **START** with ST4 + NL1 + BT4 + DB1 for general applications
2. **USE** ST6 for extreme event studies (hurricanes, severe storms)
3. **ALWAYS** tune `BETAMAX` for your specific wind forcing product
4. **VALIDATE** against local buoy data when available
5. **USE** the 30-year hindcast as a benchmark
6. **CONSIDER** computational cost: NL1 (DIA) for operational; NL2 (WRT) for research only
7. **FOR COASTAL** applications: add TR1, DB1, and consider REF1 and BT4
8. **FOR ARCTIC/ANTARCTIC**: select appropriate IC (IC3 or IC5) and IS (IS2) switches
9. **FOR FIRST RUNS**: Type 0 (calm) initial conditions with 3-7 day spin-up
10. **FOR GLOBAL GRIDS**: use `SMPL` closure for periodic longitude wrapping

---

## References & Resources

### Official Documentation
- NOAA model page: https://polar.ncep.noaa.gov/waves/wavewatch/
- GitHub Wiki: https://github.com/NOAA-EMC/WW3/wiki/
- Online tutorial: https://ww3-docs.readthedocs.io/en/latest/
- User Manual v6.07: https://raw.githubusercontent.com/wiki/NOAA-EMC/WW3/files/manual.pdf
- Operational page: https://polar.ncep.noaa.gov/waves/index2.shtml
- Validation archive: https://polar.ncep.noaa.gov/waves/validation/

### Source Code & Tools
- GitHub (source): https://github.com/NOAA-EMC/WW3
- WW3-tools: https://github.com/NOAA-EMC/WW3-tools
- UFS Coastal: https://github.com/oceanmodeling/ufs-coastal-app
- WW3-CESM: https://github.com/ESCOMP/WW3-CESM

### Python Tools
- wavespectra: https://wavespectra.readthedocs.io/ (`pip install wavespectra`)
- pyww3: https://github.com/caiostringari/pyww3
- APPMAR: https://github.com/cemanetwork/appmar

### Hindcast Data
- NOAA hindcasts: https://polar.ncep.noaa.gov/waves/hindcasts/
- IFREMER: https://data-ww3.ifremer.fr/
- NDBC buoys: https://www.ndbc.noaa.gov/
- CDIP: https://cdip.ucsd.edu/

### Key References

- Ardhuin, F., et al. (2010). Semiempirical dissipation source functions for ocean waves. *J. Phys. Oceanogr.*, 40, 917-941.
- Battjes, J. A. and Janssen, J. P. F. M. (1978). Energy loss and set-up due to breaking of random waves. *Proc. 16th ICCE*, 569-587.
- Hasselmann, S., et al. (1985). Computations and parameterizations of the nonlinear energy transfer. *J. Phys. Oceanogr.*, 15, 1378-1391.
- Rogers, W. E., et al. (2012). Observation-consistent input and whitecapping dissipation. *J. Atmos. Ocean. Tech.*, 29, 1329-1346.
- Tolman, H. L. (2008). A mosaic approach to wind wave modeling. *Ocean Modelling*, 25, 35-47.
- Tolman, H. L. and Chalikov, D. (1996). Source terms in a third-generation wind wave model. *J. Phys. Oceanogr.*, 26, 2497-2518.
