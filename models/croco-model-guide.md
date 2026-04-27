# CROCO: Complete Guide from Grid Setup to Visualization

> A practical guide covering the full CROCO (Coastal and Regional Ocean COmmunity) workflow — from grid generation and forcing preparation through compilation, configuration, AGRIF online nesting, non-hydrostatic dynamics, and post-processing.
> Regional focus on **Indonesian waters and BMKG** operational context.

---

## Table of Contents

1. [Overview & Key Concepts](#1-overview--key-concepts)
2. [CROCO vs ROMS Forks](#2-croco-vs-roms-forks)
3. [System Requirements & Installation](#3-system-requirements--installation)
4. [CROCO_PYTOOLS: Python Pre/Post-Processing](#4-croco_pytools-python-prepost-processing)
5. [Grid Generation](#5-grid-generation)
6. [Compilation](#6-compilation)
7. [Configuration Files](#7-configuration-files)
8. [Vertical S-Coordinate](#8-vertical-s-coordinate)
9. [Physics Options](#9-physics-options)
10. [AGRIF Online 2-Way Nesting](#10-agrif-online-2-way-nesting)
11. [CROCO-NH: Non-Hydrostatic Mode](#11-croco-nh-non-hydrostatic-mode)
12. [Tidal Forcing](#12-tidal-forcing)
13. [Running CROCO](#13-running-croco)
14. [Output & Post-Processing](#14-output--post-processing)
15. [Indonesian Waters Application](#15-indonesian-waters-application)
16. [Troubleshooting & Tips](#16-troubleshooting--tips)

---

## 1. Overview & Key Concepts

### What Is CROCO?

CROCO (Coastal and Regional Ocean COmmunity model) is a free-surface, terrain-following, primitive-equations ocean model descended from ROMS-UCLA. It is developed and maintained by a French consortium — IRD, INRIA, IFREMER, SHOM, and CNRS — with active international contributions. CROCO distinguishes itself from the parent ROMS lineage through three major additions: **AGRIF-based online 2-way nesting**, a **non-hydrostatic kernel (CROCO-NH)**, and a modern Python toolchain (CROCO_PYTOOLS). It is widely used in European operational systems, the SWIO-RAFI Indian Ocean project, and increasingly in Southeast Asian applications including BMKG research programs.

The current stable release is **CROCO v2.1** (November 2025). The GitLab repository at `https://gitlab.inria.fr/croco-ocean/croco` is the canonical source.

### Governing Equations

#### Hydrostatic Primitive Equations (default)

CROCO solves the **hydrostatic Boussinesq primitive equations** on an Arakawa C-grid in terrain-following S-coordinates:

**Horizontal momentum:**
```
∂u/∂t + (u·∇)u − fv = −(1/ρ₀) ∂P/∂x + ∂/∂z(Kₘ ∂u/∂z) + Fᵤ
∂v/∂t + (u·∇)v + fu = −(1/ρ₀) ∂P/∂y + ∂/∂z(Kₘ ∂v/∂z) + Fᵥ
```

**Hydrostatic balance:**
```
∂P/∂z = −ρg
```

**Continuity (incompressible):**
```
∂u/∂x + ∂v/∂y + ∂w/∂z = 0
```

**Free-surface:**
```
∂ζ/∂t + ∂(D·ū)/∂x + ∂(D·v̄)/∂y = 0
```
where D = h + ζ is the total water depth.

**Tracer equations** (temperature T and salinity S, plus any passive tracers):
```
∂T/∂t + u·∇T = ∂/∂z(Kₕ ∂T/∂z) + Fₜ
∂S/∂t + u·∇S = ∂/∂z(Kₕ ∂S/∂z) + Fₛ
```

**Equation of state:** ρ = ρ(T, S, P) — nonlinear UNESCO-80 or simplified linear form.

#### AGRIF Multi-Scale Extension

AGRIF (Adaptive Grid Refinement In Fortran) embeds one or more child grids within the parent at compile time. Each child solves the same equation set at finer resolution. The 2-way coupling exchanges fluxes at every barotropic substep across the nest boundary:

```
Parent grid (Δx_p, Δt_p)  <---->  Child grid (Δx_p/r, Δt_p/r)
         |                                    |
   coarse solution                     refined solution
   updated by child                  forced by parent BC
```
where r is the spatial refinement ratio (must be odd: 3, 5, 7).

#### CROCO-NH: Non-Hydrostatic Extension

When `key_nhyd` is activated, the hydrostatic balance is replaced by the full vertical momentum equation:

```
∂w/∂t + (u·∇)w = −(1/ρ₀) ∂P_nh/∂z − g(ρ/ρ₀) + ∂/∂z(Kₘ ∂w/∂z)
```

where P_nh is the non-hydrostatic pressure perturbation solved via a 3D elliptic equation (FFTW-based solver). This adds ~15–25% computational overhead but is essential for resolving internal solitary waves, convective plumes, and submesoscale fronts at grid spacings below ~500 m.

### Key Features Table

| Feature | Description |
|---|---|
| **Vertical coordinate** | Terrain-following S-coordinate (Song & Haidvogel 1994; Shchepetkin 2005). Same system as ROMS. |
| **Horizontal grid** | Arakawa C-grid, orthogonal curvilinear coordinates, spherical or Cartesian |
| **Time stepping** | Split-explicit: LF-AM3 (3D momentum) + FB AB3-AM4 (2D barotropic). Max 3D Courant ≈ 1.58 |
| **Nesting** | AGRIF online 2-way with odd refinement ratios up to 1:9 |
| **Non-hydrostatic** | CROCO-NH: full non-hydrostatic pressure via FFTW elliptic solver |
| **Biogeochemistry** | PISCES, BioEBUS, NPZD, Fennel (CPP keys) |
| **Sediment** | MUSTANG (Mud and Sand Transport Across the Globe) |
| **Waves** | Coupling to SWAN/WW3 via MCT or OASIS3-MCT |
| **Parallelism** | MPI domain decomposition; OpenMP within tiles |
| **I/O** | NetCDF-3 or NetCDF-4/HDF5. XIOS parallel I/O server supported |
| **Python toolchain** | CROCO_PYTOOLS (pre- and post-processing), CROCO_TOOLS (MATLAB legacy) |

### CROCO Architecture (ASCII Diagram)

```
+------------------------------------------------------------------+
|                        CROCO Executable                          |
|                                                                  |
|  +----------------------+    +------------------------------+    |
|  |   OCEAN/ kernel      |    |   AGRIF/ library             |    |
|  |  (Fortran 90/95)     |    |  (adaptive refinement)       |    |
|  |                      |    |                              |    |
|  |  3D baroclinic step  |<-->|  Parent <-> Child exchange   |    |
|  |  2D barotropic step  |    |  Flux correction at borders  |    |
|  |  Tracer advection    |    |  Time interpolation          |    |
|  |  Turbulence closure  |    +------------------------------+    |
|  |  NH pressure solver  |                                        |
|  +----------------------+    +------------------------------+    |
|                              |   PISCES / BioEBUS / MUSTANG |    |
|  +----------------------+    |  (optional BGC/sediment)     |    |
|  |   SCRIPTS/           |    +------------------------------+    |
|  |  jobcomp             |                                        |
|  |  run_croco.sh        |    +------------------------------+    |
|  +----------------------+    |   CROCO_PYTOOLS (Python)     |    |
|                              |  make_grid / make_forcing    |    |
|  Input files:                |  make_clim / make_bry        |    |
|  croco.in                    |  make_tides / diagnostics    |    |
|  AGRIF_FixedGrids.in         +------------------------------+    |
+------------------------------------------------------------------+
```

---

## 2. CROCO vs ROMS Forks

### History

ROMS originated from collaboration between Rutgers (Arango) and UCLA (Shchepetkin, McWilliams) in the late 1990s. The codes diverged around 1999 when UCLA adopted LF-AM3 for 3D momentum and the FB AB3-AM4 scheme for the fast 2D mode, yielding a higher maximum Courant number. CROCO was forked from ROMS-UCLA by French institutes (IRD, INRIA, IFREMER, SHOM, CNRS) around 2011, adding AGRIF, non-hydrostatic capabilities, PISCES biogeochemistry, and MUSTANG sediment. The CROCO project is governed by an open community model under a mixed MIT/CeCILL-C license.

### Detailed Comparison

| Feature | ROMS-Rutgers | ROMS-UCLA | CROCO |
|---|---|---|---|
| **Lead institution** | Rutgers / Hernan Arango | UCLA / Shchepetkin & McWilliams | IRD, INRIA, IFREMER, SHOM, CNRS (France) |
| **Repository** | github.com/myroms/roms | github.com/CESR-lab/ucla-roms | gitlab.inria.fr/croco-ocean/croco |
| **License** | MIT/X | — | MIT/X + CeCILL-C (AGRIF) |
| **3D time stepping** | Adams-Bashforth 3 | LF-AM3 | LF-AM3 |
| **2D time stepping** | LF-AM3 with feedback | FB AB3-AM4 | FB AB3-AM4 |
| **Max 3D Courant** | ~0.72 | ~1.58 | ~1.58 |
| **Max barotropic Courant** | ~0.9 | ~1.5 | ~1.5 |
| **Data assimilation** | 4D-Var (IS4DVAR, I4DVAR) — best-in-class | 3D-Var | 3D-Var (under development) |
| **Online nesting** | No (offline only via contact points) | No | **Yes — AGRIF 2-way** |
| **Non-hydrostatic** | No | Experimental only | **Yes — CROCO-NH (FFTW)** |
| **Biogeochemistry** | BEC, Fennel (community add-ons) | Limited | PISCES, BioEBUS, NPZD |
| **Sediment** | CSTMS (COAWST) | — | MUSTANG |
| **Vertical coordinates** | S (Vtransform 1 or 2) | S (Vtransform 2) | S (Vtransform 2) |
| **Tracer advection** | U3, C4, HSIMT, MPDATA | U3, WENO | RSUP3 (default), WENO5 |
| **Pressure gradient** | DJ_GRADPS (default) | WJ / DJ | WJ_GRADP (default), DJ option |
| **Community** | myroms.org forum | Limited | croco-ocean.org, annual workshop |
| **Latest release** | Continuous (GitHub trunk) | — | v2.1.2 (Nov 2025) |
| **Python tools** | roms-tools (CWorthy) | pyroms | CROCO_PYTOOLS v2.0+ |
| **MATLAB tools** | ROMS Toolbox (myroms SVN) | ROMS_TOOLS | CROCO_TOOLS |

### When to Choose CROCO

- **Online 2-way nesting required:** AGRIF is production-quality in CROCO; offline ROMS nesting is more cumbersome.
- **Non-hydrostatic dynamics:** Submesoscale (< 500 m), internal solitary waves, convective cells.
- **French/European collaborations:** IRD, IFREMER, SHOM operational systems all use CROCO.
- **PISCES biogeochemistry:** PISCES is tightly integrated in CROCO; requires significant effort in ROMS-Rutgers.
- **Larger Courant number:** The LF-AM3 / AB3-AM4 scheme allows ~2x larger time steps than ROMS-Rutgers for the same grid.

---

## 3. System Requirements & Installation

### Dependencies

| Library | Version | Purpose | Notes |
|---|---|---|---|
| **Fortran compiler** | gfortran ≥ 9, ifort ≥ 19, ifx | Primary language | Intel recommended for production |
| **C compiler** | gcc ≥ 8, icc/icx | C preprocessing, AGRIF | Must match MPI wrappers |
| **NetCDF-Fortran** | ≥ 4.5 | All I/O | Set `NETCDF_INCDIR`, `NETCDF_LIBDIR` |
| **NetCDF-C** | ≥ 4.7 | Underlying I/O | HDF5 must be compiled in |
| **HDF5** | ≥ 1.10 | NetCDF-4 compression, parallel I/O | Parallel HDF5 needed for `key_mpi_island` |
| **MPI** | OpenMPI ≥ 4.0, Intel MPI ≥ 2019 | Distributed memory | `mpif90` / `mpicc` wrappers |
| **FFTW3** | ≥ 3.3 | CROCO-NH elliptic solver | Only needed with `key_nhyd` |
| **Python** | ≥ 3.9 | CROCO_PYTOOLS | See Section 4 |
| **Perl** | ≥ 5.x | `sfmakedepend` dependency script | Standard on Linux/macOS |
| **GNU Make** | ≥ 4.0 | Build system | |

### Getting the Source

```bash
# Clone latest development version
git clone https://gitlab.inria.fr/croco-ocean/croco.git
cd croco

# Or clone a specific stable release tag
git clone --branch v2.1.2 https://gitlab.inria.fr/croco-ocean/croco.git croco-v2.1.2
cd croco-v2.1.2

# Check out the CROCO_PYTOOLS submodule (Python pre/post-processing)
git submodule update --init --recursive
```

### Directory Structure

```
croco/
├── OCEAN/              # Core Fortran source (*.F, *.h)
│   ├── cppdefs.h       # Master CPP flag definitions
│   ├── param.h         # Compile-time array dimensions
│   ├── croco.F         # Main program
│   ├── step3d_t.F      # 3D tracer advection/diffusion
│   ├── step3d_uv1.F    # 3D momentum predictor
│   ├── step3d_uv2.F    # 3D momentum corrector
│   ├── step2d.F        # 2D barotropic mode
│   ├── prsgrd.F        # Pressure gradient
│   ├── gls_mixing.F    # GLS turbulence closure
│   └── ...
├── AGRIF/              # AGRIF library and refinement routines
│   ├── AGRIF_YOURFILES.F
│   └── ...
├── SCRIPTS/            # Build and run scripts
│   ├── jobcomp         # Main compilation script (bash)
│   ├── Makefile.config # Compiler/library configuration
│   └── run_croco.sh    # Example SLURM submission
├── TEST_CASES/         # Standard test case configurations
│   ├── BENGUELA_LR/    # Benguela upwelling (low resolution)
│   ├── ROMS_BASIN/     # Idealized basin test
│   └── ...
├── CROCO_TOOLS/        # MATLAB pre/post-processing toolbox
└── CROCO_PYTOOLS/      # Python pre/post-processing (submodule)
    ├── crocotools_py/
    ├── make_grid.py
    ├── make_forcing.py
    ├── make_clim.py
    ├── make_bry.py
    └── make_tides.py
```

### Environment Setup

Add to `~/.bashrc` or `~/.bash_profile` (adjust paths to your HPC installation):

```bash
# --- CROCO environment ---
export CROCO_ROOT=/home/user/models/croco-v2.1.2

# NetCDF (compiled with HDF5 support)
export NETCDF_ROOT=/opt/netcdf/4.9.2-intel2023
export NETCDF_INCDIR=$NETCDF_ROOT/include
export NETCDF_LIBDIR=$NETCDF_ROOT/lib

# HDF5
export HDF5_ROOT=/opt/hdf5/1.14.2-intel2023

# FFTW3 (only for CROCO-NH)
export FFTW_ROOT=/opt/fftw/3.3.10-intel2023
export FFTW_INCDIR=$FFTW_ROOT/include
export FFTW_LIBDIR=$FFTW_ROOT/lib

# MPI (Intel MPI or OpenMPI)
export MPI_ROOT=/opt/intel/mpi/2023
export PATH=$MPI_ROOT/bin:$PATH

# Python tools
export PYTHONPATH=$CROCO_ROOT/CROCO_PYTOOLS:$PYTHONPATH
```

On a **BMKG HPC (SLURM)** cluster, these are typically loaded via modules:

```bash
module load intel/2023.2
module load intelmpi/2023.2
module load netcdf-fortran/4.6.1-intel
module load hdf5/1.14.2-intel
module load fftw/3.3.10-intel
module load python/3.11
```

---

## 4. CROCO_PYTOOLS: Python Pre/Post-Processing

CROCO_PYTOOLS (v2.0+) is the modern Python replacement for the legacy MATLAB CROCO_TOOLS. It handles grid generation, forcing preparation, climatology, boundary conditions, and tidal forcing. All tools read/write CROCO-formatted NetCDF files.

### Installation

```bash
# Option 1: conda environment (recommended)
conda create -n croco_env python=3.11
conda activate croco_env
conda install -c conda-forge numpy scipy matplotlib netCDF4 xarray cartopy \
    pyproj shapely cmocean tqdm

# Install CROCO_PYTOOLS itself (from the croco submodule)
cd /home/user/models/croco-v2.1.2/CROCO_PYTOOLS
pip install -e .

# Option 2: pip only
pip install numpy scipy matplotlib netCDF4 xarray cartopy pyproj shapely cmocean tqdm
pip install -e /home/user/models/croco-v2.1.2/CROCO_PYTOOLS
```

### Configuration File: `croco_pytools.cfg`

Each CROCO_PYTOOLS script reads a configuration file (INI-style) that sets domain parameters, data paths, and output filenames:

```ini
[DOMAIN]
lonmin   = 90.0
lonmax   = 145.0
latmin   = -15.0
latmax   = 10.0
resolution = 1/36   ; ~3 km

[GRID]
grid_file     = /data/croco/input/indonesia_grid.nc
hmin          = 10.0
hmax          = 5500.0
smooth_coef   = 0.2
rx0_max       = 0.20
rx1_max       = 12.0

[FORCING]
era5_dir      = /data/ERA5/
output_file   = /data/croco/input/indonesia_forcing.nc

[TIDES]
tpxo_dir      = /data/TPXO10/
fes2022_dir   = /data/FES2022/
output_file   = /data/croco/input/indonesia_tides.nc
constituents  = M2 S2 N2 K2 K1 O1 P1 Q1

[GLORYS]
glorys_dir    = /data/GLORYS12/
clim_file     = /data/croco/input/indonesia_clim.nc
bry_file      = /data/croco/input/indonesia_bry.nc
```

### `make_grid.py` — Grid Generation

```python
#!/usr/bin/env python3
"""
make_grid.py — Generate CROCO grid for Indonesian Seas domain
"""
from crocotools_py.croco_grid import make_grid

make_grid(
    # Domain extent
    lonmin=90.0, lonmax=145.0,
    latmin=-15.0, latmax=10.0,

    # Resolution: ~3 km at equator (1/36 degree)
    dl=1/36,

    # Bathymetry source (GEBCO 2023 preferred)
    bathy_file='/data/GEBCO_2023/GEBCO_2023.nc',
    bathy_lon_name='lon',
    bathy_lat_name='lat',
    bathy_depth_name='elevation',   # negative = ocean

    # Minimum depth (wet cells), maximum depth cap
    hmin=10.0,
    hmax=5500.0,

    # Smoothing targets
    rx0_max=0.20,       # Beckmann-Haidvogel number
    rx1_max=12.0,       # Haney number

    # Vertical grid (for rx1 calculation only; actual params in croco.in)
    N=40,
    theta_s=6.0,
    theta_b=2.0,
    hc=50.0,

    # Output
    output_file='/data/croco/input/indonesia_grid.nc',
    title='Indonesian Seas CROCO Grid 3km',
)
```

Key output variables in the grid NetCDF:

| Variable | Description |
|---|---|
| `lon_rho`, `lat_rho` | Coordinates at RHO-points (cell centers) |
| `lon_u`, `lat_u` | Coordinates at U-points (west/east edges) |
| `lon_v`, `lat_v` | Coordinates at V-points (south/north edges) |
| `h` | Bottom depth (m), positive downward |
| `mask_rho` | Land/sea mask at RHO-points (1=ocean, 0=land) |
| `pm`, `pn` | Inverse grid spacing 1/Δx, 1/Δy (m⁻¹) |
| `f` | Coriolis parameter (s⁻¹) |
| `angle` | Grid rotation angle from east (rad) |
| `rx0`, `rx1` | Stiffness metrics (diagnostic) |

### `make_forcing.py` — ERA5 Atmospheric Forcing

```python
#!/usr/bin/env python3
"""
make_forcing.py — Prepare ERA5 bulk flux forcing for CROCO
"""
from crocotools_py.croco_forcing import make_forcing

make_forcing(
    grid_file='/data/croco/input/indonesia_grid.nc',

    # ERA5 source data (downloaded via CDS API)
    era5_dir='/data/ERA5/indonesia/',
    year_start=2020, month_start=1,
    year_end=2020,   month_end=12,

    # Variables to process (ERA5 short names)
    variables=['u10', 'v10', 't2m', 'd2m', 'sp', 'tp',
               'ssrd', 'strd', 'msl'],

    # Interpolation
    interp_method='bilinear',

    # Output
    output_file='/data/croco/input/indonesia_frc_2020.nc',
    time_ref='days since 2000-01-01 00:00:00',
)
```

ERA5 variables required and their unit conversions:

| ERA5 Variable | Short Name | CROCO Variable | Unit Conversion |
|---|---|---|---|
| 10m U-wind | `u10` | `uwnd` | m/s → m/s (none) |
| 10m V-wind | `v10` | `vwnd` | m/s → m/s (none) |
| 2m temperature | `t2m` | `tair` | K → °C (−273.15) |
| 2m dewpoint | `d2m` | `rhum` | K → relative humidity |
| Surface pressure | `sp` | `pres` | Pa → mbar (÷100) |
| Total precipitation | `tp` | `rain` | m/hr → kg/m²/s (×1000/3600) |
| Downward SW radiation | `ssrd` | `radsw` | J/m² → W/m² (÷3600, accumulated) |
| Downward LW radiation | `strd` | `radlw` | J/m² → W/m² (÷3600, accumulated) |

### `make_clim.py` — Climatological Boundary Data

```python
#!/usr/bin/env python3
"""
make_clim.py — Interpolate GLORYS12 to CROCO S-grid (climatology)
"""
from crocotools_py.croco_clim import make_clim

make_clim(
    grid_file='/data/croco/input/indonesia_grid.nc',

    # GLORYS12 monthly mean files
    glorys_dir='/data/GLORYS12/monthly/',
    glorys_prefix='cmems_mod_glo_phy_my_0.083deg_P1M-m_',

    # Vertical coordinate parameters (must match croco.in)
    N=40,
    theta_s=6.0, theta_b=2.0, hc=50.0,
    vtransform=2, vstretching=4,

    # Variables to extract
    variables=['thetao', 'so', 'uo', 'vo', 'zos'],

    # Time range
    year_start=2015, year_end=2020,  # climatological average

    output_file='/data/croco/input/indonesia_clim.nc',
)
```

### `make_bry.py` — Time-Varying Boundary Conditions

```python
#!/usr/bin/env python3
"""
make_bry.py — Time-varying open boundary conditions from GLORYS12
"""
from crocotools_py.croco_bry import make_bry

make_bry(
    grid_file='/data/croco/input/indonesia_grid.nc',

    glorys_dir='/data/GLORYS12/daily/',

    # Which boundaries are open (West, East, South, North)
    obc=[1, 1, 1, 1],   # 1 = open, 0 = closed

    N=40,
    theta_s=6.0, theta_b=2.0, hc=50.0,
    vtransform=2, vstretching=4,

    year=2020, month_start=1, month_end=12,

    output_file='/data/croco/input/indonesia_bry_2020.nc',
)
```

### `make_tides.py` — Tidal Forcing

```python
#!/usr/bin/env python3
"""
make_tides.py — Extract tidal constituents from TPXO10 or FES2022
"""
from crocotools_py.croco_tides import make_tides

make_tides(
    grid_file='/data/croco/input/indonesia_grid.nc',

    # Tidal atlas: 'tpxo10' or 'fes2022'
    tidal_model='tpxo10',
    tpxo_dir='/data/TPXO10/',

    # Constituents (8 primary for Indonesian waters)
    constituents=['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1'],

    # Include tidal current forcing at boundaries
    uv_tides=True,

    output_file='/data/croco/input/indonesia_tides.nc',
)
```

---

## 5. Grid Generation

### Grid Philosophy

CROCO uses the same Arakawa C-grid as ROMS. All field variables sit on staggered grid points:

```
+-------v(i,j+1)---+-------v(i+1,j+1)-+
|                   |                   |
u(i,j)   rho(i,j)  u(i+1,j)  rho(i+1,j) u(i+2,j)
|                   |                   |
+-------v(i,j)-----+-------v(i+1,j)---+
         psi(i,j)            psi(i+1,j)
```

Grid dimensions in the NetCDF file are `(eta_rho, xi_rho)` = `(Mm+2, Lm+2)` where Lm and Mm are the number of interior RHO-points.

### Bathymetry Sources

| Dataset | Resolution | Notes |
|---|---|---|
| **GEBCO 2023** | 15 arc-sec (~450 m) | Most complete global compilation; recommended |
| **ETOPO 2022** | 15 arc-sec | NOAA product combining ocean + land |
| **SRTM30_PLUS v11** | 30 arc-sec (~900 m) | Smith & Sandwell satellite gravity + sonar |
| **BIG (BathyIndo Grid)** | Regional | PUSHIDROSAL Indonesian naval charts, higher resolution near coasts |
| **GEBCO + regional blending** | Variable | Merge high-res local survey into GEBCO background |

For Indonesian waters, GEBCO 2023 is the baseline. The shallow shelf seas (Java Sea, Banda Sea) and critical sills (Lombok Strait ~250 m, Ombai Strait ~1250 m, Lifamatola Passage ~1940 m) must be checked manually and not over-smoothed.

### Smoothing Parameters

Terrain-following coordinates generate **pressure gradient errors (PGE)** over steep topography. CROCO uses two stiffness metrics:

**rx0 (Beckmann-Haidvogel number):**
```
rx0 = max |h(i) − h(j)| / [h(i) + h(j)]  for adjacent cells i,j
```
Target: **rx0 ≤ 0.20** for CROCO. Strictly below 0.15 near strong currents.

**rx1 (Haney number):**
```
rx1 = max |(z_k(i) − z_k(j)) + (z_{k-1}(i) − z_{k-1}(j))|
          / |(z_k(i) − z_{k-1}(i)) + (z_k(j) − z_{k-1}(j))|
```
Target: **rx1 ≤ 12** for CROCO (more permissive than ROMS-Rutgers due to WJ pressure gradient scheme; but < 8 is safer).

### Python Grid Generation for Indonesian Domain

```python
#!/usr/bin/env python3
"""
Full Indonesian Seas grid generation with bathymetry smoothing.
Domain: 90°E–145°E, 15°S–10°N, ~3 km resolution (1/36°).
"""
import numpy as np
import xarray as xr
from crocotools_py.croco_grid import make_grid, smooth_bathy, check_stiffness

# ---- Step 1: Generate base grid ----------------------------------------
grid_file = '/data/croco/input/indonesia_grid.nc'

make_grid(
    lonmin=90.0,  lonmax=145.0,
    latmin=-15.0, latmax=10.0,
    dl=1/36,
    bathy_file='/data/GEBCO_2023/GEBCO_2023.nc',
    bathy_lon_name='lon',
    bathy_lat_name='lat',
    bathy_depth_name='elevation',
    hmin=10.0,
    hmax=5500.0,
    rx0_max=0.20,
    rx1_max=12.0,
    N=40, theta_s=6.0, theta_b=2.0, hc=50.0,
    output_file=grid_file,
)

# ---- Step 2: Check stiffness metrics ------------------------------------
ds = xr.open_dataset(grid_file)
h = ds['h'].values

rx0, rx1 = check_stiffness(h, N=40, theta_s=6.0, theta_b=2.0, hc=50.0)
print(f"Max rx0 = {rx0.max():.3f}  (target < 0.20)")
print(f"Max rx1 = {rx1.max():.3f}  (target < 12.0)")

# ---- Step 3: Iterative smoothing if needed ------------------------------
if rx0.max() > 0.20:
    h_smooth = smooth_bathy(
        h,
        rx0_max=0.20,
        method='lp',          # Linear programming (Sikirić 2009)
        hmin=10.0,
        preserve_sills=True,  # Do not deepen sills — critical for ITF
    )
    # Write back to grid file
    ds['h'].values = h_smooth
    ds.to_netcdf(grid_file.replace('.nc', '_smooth.nc'))

# ---- Step 4: Edit land mask manually if needed -------------------------
# Use CROCO_PYTOOLS mask editor or edit mask_rho directly
# Typical edits: open/close passages (Sunda, Lombok, Makassar, Ombai)

# ---- Step 5: Check critical passage depths ------------------------------
# Lombok Strait sill: ~9.5°S, 115.8°E → should be ~250 m
# Ombai Strait sill:  ~8.5°S, 124.5°E → should be ~1250 m
# Lifamatola Pass:    ~1.8°S, 127.2°E → should be ~1940 m
print("Lombok Strait max depth:", h[int((-9.5+15)*36), int((115.8-90)*36)])
```

### Child Grid Considerations (AGRIF)

When creating child grids for AGRIF nesting:
- Child grid resolution must be exactly `parent_resolution / r` where r is an odd integer
- Child grid boundaries must lie at least 3 parent grid cells from the parent boundary
- Child grid file uses the same format as the parent grid file
- Use `make_grid.py` separately for each child domain, then verify alignment

---

## 6. Compilation

### The `jobcomp` Script

CROCO's build system centres on `SCRIPTS/jobcomp`, a bash script that sets CPP keys, domain dimensions, and invokes `make`. Unlike ROMS `build_roms.sh`, CROCO's `jobcomp` is application-specific and must be customised per simulation.

**Annotated `jobcomp` for Indonesian domain:**

```bash
#!/bin/bash
# jobcomp — CROCO compilation script for Indonesian Seas
# Usage: ./jobcomp 2>&1 | tee compile.log

# ---- Application name ---------------------------------------------------
MYAPP=INDONESIA       # Must match CPP flag in cppdefs.h or header file

# ---- Source directories -------------------------------------------------
CROCO_ROOT=/home/user/croco-v2.1.2
cd $CROCO_ROOT/OCEAN

# ---- Grid dimensions (set in param.h or via -D flags) ------------------
# LLm, MMm: global interior RHO-points (before MPI tiling)
# N:        number of vertical levels
# Must be consistent with grid NetCDF file
CPPFLAGS="-DLLm=1800 -DMMm=900 -DN=40"

# ---- Physics CPP flags --------------------------------------------------
CPPFLAGS="$CPPFLAGS \
 -D${MYAPP}         \
 -DSOLVE3D          \
 -DSPHERICAL        \
 -DCURVGRID         \
 -DMASKING          \
 -DNONLIN_EOS       \
 -DBULK_FLUXES      \
 -DAERO_BULK        \
 -DCOARE_TAYLOR_YELLAND \
 -DSSH_TIDES        \
 -DUBTIDES          \
 -DTIDE_GENERATING_FORCES \
 -DRAMP_TIDES       \
 -DGLS_MIXING       \
 -DCANUTO_A         \
 -DN2S2_HORAVG      \
 -DTS_HADV_RSUP3    \
 -DUV_HADV_UP3      \
 -DVADV_ADAPT_IMP   \
 -DUV_VIS2          \
 -DMIX_S_UV         \
 -DTS_DIF2          \
 -DMIX_GEO_TS       \
 -DDJ_GRADPS        \
 -DUV_QDRAG         \
 -DAVERAGES         \
 -DAVERAGES2        \
 -DDIAGNOSTICS_TS   \
 "

# ---- AGRIF nesting (add if using nested grids) -------------------------
CPPFLAGS="$CPPFLAGS -DAGRIF -DAGRIF_2WAY"

# ---- MPI parallelisation -----------------------------------------------
CPPFLAGS="$CPPFLAGS -DMPI"
NP_XI=18     # MPI tiles in xi direction
NP_ETA=10    # MPI tiles in eta direction
CPPFLAGS="$CPPFLAGS -DNP_XI=$NP_XI -DNP_ETA=$NP_ETA"

# ---- Compiler selection (edit Makefile.config to match) ----------------
export FORT=ifort      # or gfortran

# ---- Build --------------------------------------------------------------
make -j 16 CPPFLAGS="$CPPFLAGS"
```

### `Makefile.config` — Compiler Configuration

Located at `SCRIPTS/Makefile.config`. Key sections to customise:

```makefile
# Intel Fortran + MPI
FC     = mpif90
CC     = mpicc
FFLAGS = -O3 -fp-model precise -xCORE-AVX2 -traceback
CFLAGS = -O2

# NetCDF
NETCDF_INCDIR = /opt/netcdf/4.9.2-intel/include
NETCDF_LIBDIR = /opt/netcdf/4.9.2-intel/lib
NETCDF_LIB    = -lnetcdff -lnetcdf

# HDF5
HDF5_INCDIR   = /opt/hdf5/1.14.2-intel/include
HDF5_LIBDIR   = /opt/hdf5/1.14.2-intel/lib
HDF5_LIB      = -lhdf5_fortran -lhdf5 -lz

# FFTW (only for CROCO-NH)
FFTW_INCDIR   = /opt/fftw/3.3.10-intel/include
FFTW_LIBDIR   = /opt/fftw/3.3.10-intel/lib
FFTW_LIB      = -lfftw3

# Link line
LDFLAGS = -L$(NETCDF_LIBDIR) $(NETCDF_LIB) \
          -L$(HDF5_LIBDIR) $(HDF5_LIB) \
          -L$(FFTW_LIBDIR) $(FFTW_LIB)
```

### CPP Keys Reference Table

| CPP Key | Category | Description |
|---|---|---|
| `SOLVE3D` | Dynamics | 3D primitive equations (always use) |
| `SPHERICAL` | Grid | Spherical coordinates |
| `CURVGRID` | Grid | Curvilinear orthogonal grid |
| `MASKING` | Grid | Land/sea masking |
| `NONLIN_EOS` | Thermodynamics | Nonlinear equation of state (UNESCO-80) |
| `SALINITY` | Thermodynamics | Include salinity (enabled by default with SOLVE3D) |
| `BULK_FLUXES` | Surface | Compute air-sea fluxes via bulk formula |
| `AERO_BULK` | Surface | Fairall/COARE aerodynamic bulk (preferred) |
| `SSH_TIDES` | Tides | Tidal elevation at boundaries |
| `UBTIDES` | Tides | Tidal barotropic velocity at boundaries |
| `TIDE_GENERATING_FORCES` | Tides | Astronomical tide potential |
| `RAMP_TIDES` | Tides | Gradual tidal spinup (recommended) |
| `GLS_MIXING` | Turbulence | Generic length scale closure |
| `KPP_MIXING` | Turbulence | K-profile parameterization (alternative) |
| `CANUTO_A` | Turbulence | Canuto A stability function for GLS |
| `N2S2_HORAVG` | Turbulence | Horizontal average N² and S² (stability) |
| `TS_HADV_RSUP3` | Advection | Rotated split 3rd-order upstream tracers (default) |
| `TS_HADV_WENO5` | Advection | WENO5 tracer advection (quasi-monotone) |
| `UV_HADV_UP3` | Advection | 3rd-order upstream momentum |
| `UV_HADV_WENO5` | Advection | WENO5 momentum advection |
| `VADV_ADAPT_IMP` | Advection | Adaptive implicit vertical advection |
| `DJ_GRADPS` | Pressure gradient | Density Jacobian spline (less PGE on steep topo) |
| `WJ_GRADP` | Pressure gradient | Weighted density Jacobian (CROCO default) |
| `UV_VIS2` | Mixing | Laplacian horizontal viscosity |
| `UV_VIS4` | Mixing | Biharmonic horizontal viscosity |
| `MIX_S_UV` | Mixing | Mix momentum along S-surfaces |
| `UV_SMAGORINSKY` | Mixing | Smagorinsky dynamic viscosity |
| `UV_QDRAG` | Bottom | Quadratic bottom drag |
| `UV_LOGDRAG` | Bottom | Logarithmic drag (requires roughness length) |
| `AGRIF` | Nesting | Enable AGRIF library |
| `AGRIF_2WAY` | Nesting | Two-way AGRIF coupling |
| `key_nhyd` | Non-hydrostatic | CROCO-NH non-hydrostatic pressure solver |
| `BIOLOGY` | BGC | Generic biology switch |
| `PISCES` | BGC | PISCES biogeochemical model |
| `BioEBUS` | BGC | Eastern Boundary Upwelling System BGC |
| `MUSTANG` | Sediment | MUSTANG mud/sand transport |
| `AVERAGES` | Output | Time-averaged output file |
| `DIAGNOSTICS_TS` | Output | Tracer budget diagnostics |
| `DIAGNOSTICS_UV` | Output | Momentum budget diagnostics |
| `STATIONS` | Output | Station time series |

### Executable Produced

```
croco.exe    (MPI version, produced in OCEAN/ directory)
```

Copy to your run directory:
```bash
cp OCEAN/croco.exe /home/user/runs/indonesia_2020/
```

---

## 7. Configuration Files

### 7.1 `croco.in` — Main Input File

The `croco.in` file controls all run-time parameters. An annotated example for an Indonesian Seas simulation:

```
! =====================================================================
!  CROCO configuration: Indonesian Seas — 3 km resolution
!  Period: 2020-01-01 to 2020-12-31 (366 days)
! =====================================================================

TITLE = Indonesian Seas CROCO Simulation

! --- Time stepping ---------------------------------------------------
! dt:      baroclinic time step (seconds)
! ndtfast: number of barotropic steps per baroclinic step
!          barotropic dt = dt / ndtfast
!
! For 3 km grid with shallow straits (Java Sea h~50m):
!   barotropic CFL: dtfast * sqrt(g*hmax) / dx < 1
!   dtfast ≤ 3000/(sqrt(9.8*5500)) ≈ 12 s
!
NTIMES == 10519200    ! Total baroclinic steps = 365*24*3600/120 = 262800 (1 year)
DT     == 120.0d0     ! 120-second baroclinic step
NDTFAST == 30         ! 30 barotropic steps → dtfast = 4 s
NINFO  == 1           ! Print diagnostics every N baroclinic steps

! --- Grid dimensions --------------------------------------------------
! Must match LLm, MMm, N in jobcomp CPPFLAGS
Lm == 1800            ! Interior RHO-points in xi (longitude)
Mm ==  900            ! Interior RHO-points in eta (latitude)
N  ==   40            ! Vertical S-levels

! --- MPI decomposition ------------------------------------------------
NtileI == 18          ! Tiles in xi (I) direction
NtileJ == 10          ! Tiles in eta (J) direction
                      ! Total MPI tasks = 18 × 10 = 180

! --- Vertical S-coordinate --------------------------------------------
Vtransform  == 2      ! Shchepetkin 2005 (recommended)
Vstretching == 4      ! Shchepetkin 2010 double stretching (recommended)
THETA_S == 6.0d0      ! Surface stretching 0–10: higher = more surface levels
THETA_B == 2.0d0      ! Bottom stretching 0–4: higher = more bottom levels
TCLINE  == 50.0d0     ! Transition depth hc (m); controls surf/bot blend

! --- Physical parameters ---------------------------------------------
RDRG2   == 3.0d-03    ! Quadratic bottom drag coefficient Cd (dimensionless)
Zob     == 0.02d0     ! Bottom roughness length (m); used with UV_LOGDRAG
VISC2   == 5.0d0      ! Laplacian horizontal viscosity (m²/s)
VISC4   == 0.0d0      ! Biharmonic horizontal viscosity (m⁴/s); 0 if UV_VIS2
TNU2    == 5.0d0  5.0d0  ! Laplacian diffusion for [temp, salt] (m²/s)
TNU4    == 0.0d0  0.0d0  ! Biharmonic diffusion (m⁴/s)
AKV_BAK == 1.0d-05    ! Background vertical viscosity (m²/s)
AKT_BAK == 1.0d-06  1.0d-06  ! Background vertical diffusivity [T, S] (m²/s)
AKS_BAK == 0.0d0      ! Background vertical diffusivity for passive tracers

! --- GLS turbulence parameters (k-epsilon; Canuto A) -----------------
GLS_P     == 3.0d0    ! k-epsilon exponent p
GLS_M     == 1.5d0    ! k-epsilon exponent m
GLS_N     == -1.0d0   ! k-epsilon exponent n
GLS_Kmin  == 1.0d-08  ! Min TKE (m²/s²)
GLS_Pmin  == 1.0d-12  ! Min GLS variable (m²/s³)
GLS_CMU0  == 0.5477d0 ! Stability coefficient c_mu0

! --- Lateral boundary conditions (W  S  E  N) ------------------------
LBC(isFsur) == Che Che Che Che  ! free-surface: Chapman explicit
LBC(isUbar) == Fla Fla Fla Fla  ! 2D U: Flather
LBC(isVbar) == Fla Fla Fla Fla  ! 2D V: Flather
LBC(isUvel) == RadNud RadNud RadNud RadNud  ! 3D U: radiation+nudging
LBC(isVvel) == RadNud RadNud RadNud RadNud  ! 3D V: radiation+nudging
LBC(isMtke) == Rad Rad Rad Rad  ! TKE
LBC(isTvar) == RadNud RadNud RadNud RadNud  ! temperature
              RadNud RadNud RadNud RadNud  ! salinity

! --- Nudging time scales (days) --------------------------------------
TNUDG   == 10.0d0 10.0d0    ! T, S nudging time scale
ZNUDG   == 10.0d0            ! Free-surface nudging
M2NUDG  == 10.0d0            ! 2D momentum nudging
M3NUDG  == 10.0d0            ! 3D momentum nudging
OBCFAC  == 120.0d0           ! Boundary-to-interior nudging ratio

! --- Tidal ramping ----------------------------------------------------
RAMP_TIME == 2.0d0   ! Days to ramp tidal forcing from 0 to full amplitude

! --- Restart ----------------------------------------------------------
NRREC     == 0         ! 0 = new run; −1 = use last record in restart file
LcycleRST == T         ! Overwrite restart file each cycle
NRST      == 43200     ! Write restart every 43200 steps = 60 days (at dt=120s)

! --- Output frequencies (in baroclinic steps) -------------------------
NHIS      == 720       ! History  output every 720 steps = 1 day
NAVG      == 720       ! Average  output every 720 steps = 1 day
NAVG2     == 360       ! AVERAGES2 (surface fields) every 12 hours
NDIA      == 720       ! Diagnostic output every 1 day
NSTA      == 720       ! Station output every 1 day

! --- Input file paths -------------------------------------------------
GRDNAME == /data/croco/input/indonesia_grid.nc
ININAME == /data/croco/input/indonesia_ini_20200101.nc
BRYNAME == /data/croco/input/indonesia_bry_2020.nc
CLMNAME == /data/croco/input/indonesia_clim.nc
TIDENAME == /data/croco/input/indonesia_tides.nc
FRCNAME == /data/croco/input/indonesia_frc_2020.nc
STANAME == /data/croco/input/indonesia_stations.in

! --- Output file paths ------------------------------------------------
HISNAME == /data/croco/output/indonesia_his.nc
AVGNAME == /data/croco/output/indonesia_avg.nc
AVG2NAME == /data/croco/output/indonesia_avg2.nc
DIANAME == /data/croco/output/indonesia_dia.nc
RSTNAME == /data/croco/output/indonesia_rst.nc
STANAME == /data/croco/output/indonesia_sta.nc
```

### 7.2 `AGRIF_FixedGrids.in` — Nesting Configuration

When AGRIF is enabled, this file defines nested grids (see Section 10 for full details):

```
1                   ! Number of child grids
! istart iend jstart jend spaceref timeref
  120    280   60    170    3         3
```

### 7.3 Station File

Point observations output:
```
! STANAME = /data/croco/input/indonesia_stations.in
! Format: lon lat
106.85  -6.10   ! Jakarta Bay station
110.40  -6.97   ! Karimunjawa
115.20  -8.72   ! Bali coastal
```

---

## 8. Vertical S-Coordinate

### Coordinate Transform

CROCO inherits the ROMS S-coordinate system. The depth z at S-level k at horizontal position (i,j) is:

**Vtransform = 2 (recommended, Shchepetkin 2005):**
```
z(x,y,σ,t) = ζ(x,y,t) + [ζ(x,y,t) + h(x,y)] · S(x,y,σ)

where:
S(x,y,σ) = [hc·σ + h(x,y)·C(σ)] / [hc + h(x,y)]
```

**Vstretching = 4 (recommended, Shchepetkin 2010):**
```
C(σ) = [1 − cosh(θ_s·σ)] / [cosh(θ_s) − 1]     (surface, σ ∈ [−1, 0])

C(σ) = [1 − C_surface(σ)] · sinh(θ_b·σ) / sinh(θ_b)
      + C_surface(σ)
      + (with additional double-stretching)
```

The combination Vtransform=2 + Vstretching=4 provides:
- Smooth behavior near sea surface in all water depths
- Independent control of surface and bottom resolution via θ_s and θ_b
- No pathological behavior in very shallow cells (unlike Vtransform=1)

### Parameter Guidance

| Parameter | Description | Typical Range | Effect |
|---|---|---|---|
| `theta_s` | Surface stretching | 0–10 (Vstretching=4) | Higher → denser levels near surface |
| `theta_b` | Bottom stretching | 0–4 (Vstretching=4) | Higher → denser levels near bottom |
| `hc` (TCLINE) | Transition depth (m) | 10–300 m | Above hc: surface stretching dominates |
| `N` | Number of S-levels | 20–60 | More levels = better resolved shear layers |

### Recommended Configurations for Indonesian Waters

**Full Indonesian Seas domain (depths 10–5500 m):**
```
Vtransform  = 2
Vstretching = 4
N           = 40
theta_s     = 6.0    ! Strong surface refinement for mixed layer
theta_b     = 2.0    ! Moderate bottom refinement for BBL
hc (TCLINE) = 50.0   ! 50 m transition: good for mixed layer ~20-50 m
```
Approximate level spacing at a 1000 m deep location:
- Surface (top 5 levels): 1–5 m spacing
- Mid-water: 20–50 m spacing
- Bottom (bottom 5 levels): 5–15 m spacing

**Sunda/Lombok Strait child grid (depths 20–500 m):**
```
Vtransform  = 2
Vstretching = 4
N           = 40
theta_s     = 7.0    ! More surface levels for baroclinic tides
theta_b     = 3.0    ! Strong bottom refinement for bottom-trapped currents
hc (TCLINE) = 20.0   ! Small; minimum depth is ~20 m
```

**Rule:** Set `hc` ≤ minimum depth in domain. If `hc > h_min`, some S-levels may flip in very shallow cells.

**Deep Banda Sea + Maluku Sea passages (depths up to 5000 m):**
```
N           = 50
theta_s     = 7.0
theta_b     = 0.5    ! Light bottom stretching; deep = no strong BBL needed
hc (TCLINE) = 300.0
```

### Vertical Profile Visualization

```python
import numpy as np
import matplotlib.pyplot as plt

def croco_scoord(N, theta_s, theta_b, hc, h, zeta=0.0, vtransform=2, vstretching=4):
    """Compute CROCO S-coordinate depths at a single column."""
    sc_r = -1 + (2*np.arange(1, N+1) - 1) / (2*N)  # RHO-levels
    if vstretching == 4:
        if theta_s > 0:
            csrf = (1 - np.cosh(theta_s*sc_r)) / (np.cosh(theta_s) - 1)
        else:
            csrf = -sc_r**2
        if theta_b > 0:
            C = (1 - theta_b)*csrf + theta_b*(-1 + np.sinh(theta_b*sc_r)/np.sinh(theta_b))
        else:
            C = csrf
    if vtransform == 2:
        S = (hc*sc_r + h*C) / (hc + h)
        z = zeta + (zeta + h) * S
    return z

# Example: vertical profile at 500 m depth (Java Sea) and 1000 m (Flores Sea)
for depth, label in [(50, 'Java Sea (50 m)'), (500, 'Flores Sea (500 m)'), (2000, 'Banda Sea (2000 m)')]:
    z = croco_scoord(N=40, theta_s=6.0, theta_b=2.0, hc=50.0, h=depth)
    print(f"\n{label}:")
    print(f"  Top 5 levels: {z[-5:][::-1].round(1)} m")
    print(f"  Bottom 5 levels: {z[:5].round(1)} m")
```

---

## 9. Physics Options

### 9.1 Turbulence Closure

#### GLS (Generic Length Scale) — Recommended for Indonesian Waters

The GLS framework (Warner et al., 2005) unifies classical two-equation closures under a single set of CPP flags:

| Closure | `GLS_P` | `GLS_M` | `GLS_N` | Best for |
|---|---|---|---|---|
| **k-ε** | 3.0 | 1.5 | −1.0 | General tropical stratification; recommended |
| **k-ω** | −1.0 | 0.5 | −1.0 | Near-surface wave breaking |
| **GEN** (Umlauf & Burchard) | 2.0 | 1.0 | −0.67 | Balanced alternative |

CROCO `croco.in` GLS parameters for k-epsilon:
```
GLS_P    == 3.0d0    ! p exponent
GLS_M    == 1.5d0    ! m exponent
GLS_N    == -1.0d0   ! n exponent
GLS_Kmin == 1.0d-08  ! minimum TKE (m²/s²)
GLS_Pmin == 1.0d-12  ! minimum dissipation
```

Sub-option CPP flags:

| Flag | Description |
|---|---|
| `CANUTO_A` | Canuto A stability function — recommended for tropical stratification |
| `CANUTO_B` | Canuto B stability function |
| `KANTHA_CLAYSON` | Kantha-Clayson stability (simpler, less accurate) |
| `CHARNOK` | Charnock surface roughness parameterized from wind stress |
| `CRAIG_BANNER` | Craig & Banner surface TKE flux from wave breaking |
| `N2S2_HORAVG` | Horizontal averaging of N² and S² for stability |
| `K_C4ADVECTION` | 4th-order advection of TKE and GLS variables |

#### KPP (K-Profile Parameterization)

Alternative to GLS; first-order closure matching surface BL to interior. More commonly used in basin-scale models. Activated with `KPP_MIXING`:

```c
#define KPP_MIXING
#define LMD_SKPP        /* Surface boundary layer */
#define LMD_BKPP        /* Bottom boundary layer */
#define LMD_RIMIX       /* Shear instability (Ri-number) */
#define LMD_CONVEC      /* Convective mixing */
#define LMD_NONLOCAL    /* Counter-gradient transport */
```

### 9.2 Advection Schemes

CROCO's default and recommended advection schemes differ from ROMS-Rutgers:

**Tracer horizontal advection:**

| CPP Flag | Description | Notes |
|---|---|---|
| `TS_HADV_RSUP3` | Rotated split 3rd-order upstream | **Default in CROCO.** Handles strongly sheared flows. |
| `TS_HADV_UP3` | Standard 3rd-order upstream | Simpler; good for weak currents |
| `TS_HADV_WENO5` | 5th-order WENOZ | Quasi-monotone; recommended for BGC tracers |
| `TS_HADV_C4` | 4th-order centered | Dispersive; not recommended alone |

**Tracer vertical advection:**

| CPP Flag | Description |
|---|---|
| `VADV_ADAPT_IMP` | Adaptive implicit — **unconditionally stable**; recommended for tidal flows |
| `TS_VADV_SPLINES` | Spline vertical advection |
| `TS_VADV_C4` | 4th-order centered vertical |

**Momentum advection:**

| CPP Flag | Description |
|---|---|
| `UV_HADV_UP3` | 3rd-order upstream — default |
| `UV_HADV_WENO5` | 5th-order WENOZ |
| `UV_SADVECTION` | Spline vertical advection |

### 9.3 Pressure Gradient Schemes

Critical for terrain-following coordinates over steep Indonesian bathymetry:

| CPP Flag | Description | Recommendation |
|---|---|---|
| `WJ_GRADP` | Weighted density Jacobian (Weighted Jacobian) | **CROCO default.** Good general performance. |
| `DJ_GRADPS` | Density Jacobian with cubic splines | Better for very steep slopes (rx0 > 0.15) |
| `PJ_GRADPQ2` | Quartic Jacobian (2nd-order) | Alternative |
| `PJ_GRADPQ4` | Quartic Jacobian (4th-order) | Highest accuracy; more expensive |

For Indonesian waters with sills (Lombok ~250 m, Lifamatola ~1940 m), use `DJ_GRADPS`.

### 9.4 Equation of State

| CPP Flag | Description |
|---|---|
| `NONLIN_EOS` | UNESCO-80 nonlinear EOS: ρ = f(T, S, P). **Always use for realistic simulations.** |
| (not defined) | Linear: ρ = ρ₀[1 − α(T − T₀) + β(S − S₀)] |

### 9.5 Horizontal Mixing

| CPP Flag | Description |
|---|---|
| `UV_VIS2` | Laplacian (∇²) horizontal viscosity |
| `UV_VIS4` | Biharmonic (∇⁴) horizontal viscosity |
| `MIX_S_UV` | Mix momentum along S-surfaces (follow terrain) |
| `MIX_GEO_UV` | Mix momentum along geopotential surfaces |
| `UV_SMAGORINSKY` | Smagorinsky dynamic viscosity: ν = (C_s·Δ)²·|S| |
| `TS_DIF2` | Laplacian tracer diffusion |
| `TS_DIF4` | Biharmonic tracer diffusion |
| `MIX_GEO_TS` | Mix tracers along geopotentials (isoneutral mixing) |
| `MIX_ISO_TS` | Isoneutral (Redi) tracer diffusion |

### 9.6 Bottom Boundary Conditions

| CPP Flag | Description |
|---|---|
| `UV_QDRAG` | Quadratic drag: τ_b = ρ·C_D·|u_b|·u_b; C_D set by `RDRG2` |
| `UV_LOGDRAG` | Logarithmic drag: C_D = [κ/ln(z_b/z₀)]²; z₀ set by `Zob` |
| `UV_LDRAG` | Linear drag (rarely used in tidal applications) |

For Indonesian tidal straits, use `UV_QDRAG` with `RDRG2 = 3.0e-3` (smooth sandy bottom) to `5.0e-3` (rough rocky sills).

---

## 10. AGRIF Online 2-Way Nesting

### What Is AGRIF?

AGRIF (Adaptive Grid Refinement In Fortran) is a Fortran library developed at INRIA that provides online structured grid refinement. In CROCO, AGRIF couples parent and child grids within a single executable, exchanging boundary conditions at every barotropic substep. This is fundamentally different from offline ROMS nesting, where the parent runs first and child boundary conditions are derived from interpolated parent output.

**Advantages over offline nesting:**
- True 2-way feedback: fine-grid features (tidal jets, eddies) influence the coarser parent
- No temporal interpolation artifacts at boundaries
- Conservation of volume, heat, salt across nest interfaces (with flux correction)
- Simpler workflow: one job, one executable

**Limitations:**
- Load imbalance: parent and child tiles must finish each time step synchronously
- Refinement ratio must be odd (3, 5, 7, 9)
- Nested grid must be properly aligned on parent grid points
- More complex debugging (AGRIF-specific diagnostics required)

### Grid Refinement Ratios

| Parent Δx | Ratio r | Child Δx | Notes |
|---|---|---|---|
| 9 km | 3 | 3 km | Suitable for mesoscale → submesoscale |
| 9 km | 9 | 1 km | Large ratio; time step also divided by 9 |
| 3 km | 3 | 1 km | Typical regional → coastal |
| 3 km | 5 | 600 m | Good for tidal strait resolution |
| 1 km | 3 | ~330 m | Approaching non-hydrostatic regime |

The child time step is: `Δt_child = Δt_parent / r`

The child NDTFAST is typically the same as the parent (barotropic sub-cycling adapts automatically).

### `AGRIF_FixedGrids.in` Syntax

This file must be present in the run directory when AGRIF is enabled:

```
! AGRIF_FixedGrids.in — Fixed grid definitions for CROCO AGRIF
! Format: number of level-1 children on first line
! Then for each child: istart iend jstart jend spaceref timeref
! istart, iend, jstart, jend: indices in PARENT grid (1-based, RHO-points)
! spaceref: spatial refinement ratio (odd integer)
! timeref:  temporal refinement ratio (equal to spaceref recommended)

2                         ! 2 first-level children

! Child 1: Sunda Strait (Selat Sunda), r=5 → 600 m
  350  430  280  350  5  5

! Child 2: Lombok Strait (Selat Lombok), r=5 → 600 m
  560  640  200  280  5  5
```

For a grandchild (level-2 nesting) embedded in child 1:
```
! Below the child 1 line, add the grandchild count and its indices:
1                         ! 1 grandchild within child 1
  20   60   25   55   3  3  ! grandchild indices in CHILD 1 frame
```

### Parent/Child Configuration Example: Sunda Strait

```
Indonesian Seas Parent (3 km, 1800×900 pts)
┌─────────────────────────────────────────────────────────────────┐
│  Indian Ocean                                                    │
│                                                                  │
│          ┌─────────────────────┐                                 │
│          │  Sunda Strait Child │  ← 600 m (r=5)                 │
│          │  (Selat Sunda)      │                                 │
│  Sumatra │  80×70 child pts    │ Java                            │
│          └─────────────────────┘                                 │
│                    Java Sea                                      │
│                    ┌──────────────────┐                          │
│                    │  Lombok Strait   │ ← 600 m (r=5)            │
│                    │  (Selat Lombok)  │                          │
│               Bali │  80×80 child pts │ Lombok                  │
│                    └──────────────────┘                          │
│  Indian Ocean             Banda Sea                              │
└─────────────────────────────────────────────────────────────────┘
```

Time step relationship:
```
Parent:  dt = 120 s, dtfast = 4 s
Child:   dt = 24 s,  dtfast = 0.8 s   (ratio r=5)
```

### AGRIF CPP Flags

```c
#define AGRIF          /* Enable AGRIF library */
#define AGRIF_2WAY     /* Two-way exchange (recommended) */
/* #define AGRIF_1WAY */  /* One-way: parent→child only */
```

### `croco.in` Parameters for AGRIF

The parent's `croco.in` contains a section for child-specific overrides:

```
! AGRIF child grid parameters (appended after main croco.in)
! Block format: one block per child, in same order as AGRIF_FixedGrids.in

! === Child 1: Sunda Strait ==========================================
NTIMES_CHILD_1 == 0      ! 0 = run same duration as parent
DT_CHILD_1     == 0.0d0  ! 0 = auto-computed from parent/ratio
NDTFAST_CHILD_1 == 30
THETA_S_CHILD_1 == 7.0d0
THETA_B_CHILD_1 == 3.0d0
TCLINE_CHILD_1  == 20.0d0
N_CHILD_1       == 40
VISC2_CHILD_1   == 1.0d0   ! Smaller viscosity at finer scale
TNU2_CHILD_1    == 1.0d0 1.0d0
RDRG2_CHILD_1   == 3.0d-03

! === Child 2: Lombok Strait =========================================
THETA_S_CHILD_2 == 7.0d0
THETA_B_CHILD_2 == 3.0d0
TCLINE_CHILD_2  == 20.0d0
N_CHILD_2       == 40
VISC2_CHILD_2   == 1.0d0
TNU2_CHILD_2    == 1.0d0 1.0d0
```

### AGRIF Grid File Naming

CROCO looks for child grid files with a specific naming convention:

```
Parent:   indonesia_grid.nc            → set in croco.in as GRDNAME
Child 1:  1_indonesia_grid.nc          → prefix "1_"
Child 2:  2_indonesia_grid.nc          → prefix "2_"
Grandchild (of child 1): 1_1_indonesia_grid.nc
```

Similarly for forcing, boundary, and initial condition files.

---

## 11. CROCO-NH: Non-Hydrostatic Mode

### When to Use CROCO-NH

The hydrostatic approximation (∂P/∂z = −ρg) is valid when horizontal scales >> vertical scales. It breaks down for:

| Application | Grid spacing | Reason |
|---|---|---|
| **Internal solitary waves (ISW)** | ~100–500 m | Nonlinear ISW in Lombok, Ombai, Luzon straits |
| **Submesoscale fronts** | < 500 m | Ageostrophic vertical velocity is O(w_hydro) |
| **Convective plumes** | < 100 m | Vertical acceleration is dominant |
| **Breaking internal tides** | < 200 m | Turbulent overturning |
| **Near-source gravity currents** | < 500 m | Overflows at sills |

For the Indonesian Seas, CROCO-NH is relevant for:
- Generation and propagation of ISW in Lombok Strait (observed amplitudes 100–150 m)
- Turbulent mixing at Lifamatola Passage overflow
- Near-field tidal dynamics at Ombai Strait

### Compilation

Add `key_nhyd` to the CPP flags in `jobcomp`:

```bash
CPPFLAGS="$CPPFLAGS -Dkey_nhyd"
```

FFTW3 must be available (linked via `Makefile.config`). The elliptic solver for the non-hydrostatic pressure perturbation uses a direct FFT-based method in the horizontal and iterative solver in the vertical.

### NH-Specific `croco.in` Parameters

```
! Non-hydrostatic solver settings
NH_EPSIL  == 1.0d-08   ! Convergence criterion for NH pressure iteration
NH_MAXITER == 100       ! Maximum iterations for NH solver
NH_METHOD == 1          ! 1: FFT; 2: SOR (Successive Over-Relaxation)
```

### Performance Overhead

| Configuration | Relative cost |
|---|---|
| Hydrostatic (default) | 1.0× |
| CROCO-NH (key_nhyd) | 1.15–1.25× |
| CROCO-NH + fine grid | 1.20–1.30× |

The overhead depends on the NH pressure solver convergence, which is slower for domains with large aspect ratios or complex topography.

### Limitations

- CROCO-NH requires Cartesian (non-spherical) coordinates for strict non-hydrostatic formulation. For small domains (< 500 km) at mid-latitudes, spherical effects are negligible.
- AGRIF + NH is supported but requires careful validation — nest boundaries must handle NH pressure gradients correctly.
- Performance scales less efficiently with MPI for NH (global elliptic solve requires global communication).

---

## 12. Tidal Forcing

### Tidal Atlases

| Atlas | Resolution | Constituents | Access | Notes |
|---|---|---|---|---|
| **TPXO10** | 1/30° global + 1/6° regional | 15 (M2, S2, N2, K2, K1, O1, P1, Q1, MM, MF, M4, MN4, MS4, 2N2, S1) | Oregon State (registration required) | Most widely validated; recommended |
| **TPXO9** | 1/30° + 1/6° | 15 | OSU | Previous version; still widely used |
| **FES2022** | 1/16° (~7 km) | 34 | Copernicus Marine / AVISO | Broader constituent spectrum; good in Indonesia |
| **FES2014b** | 1/16° | 34 | Copernicus Marine | Widely tested predecessor to FES2022 |
| **EOT20** | 1/8° | 17 | DGFI-TUM | Open access; good accuracy |

For Indonesian waters: **TPXO10** is the operational standard; **FES2022** is recommended for capturing higher harmonics (M4, MS4) important in shallow shelves.

### Important Constituents for Indonesian Waters

| Constituent | Period | Type | Indonesian Importance |
|---|---|---|---|
| **M2** | 12.42 h | Semidiurnal | Dominant in Lombok, Ombai, Timor straits |
| **S2** | 12.00 h | Semidiurnal | Second largest semidiurnal |
| **K1** | 23.93 h | Diurnal | Dominant in Java Sea, mixed tide regions |
| **O1** | 25.82 h | Diurnal | Strong diurnal complement to K1 |
| **N2** | 12.66 h | Semidiurnal | Modulates M2 amplitude (~20%) |
| **K2** | 11.97 h | Semidiurnal | Modulates S2 (~30%) |
| **P1** | 24.07 h | Diurnal | Modulates K1 (~30%) |
| **Q1** | 26.87 h | Diurnal | Modulates O1 (~20%) |
| **M4** | 6.21 h | Quarter-diurnal overtide | Nonlinear tides in shallow shelves (Java Sea) |
| **MS4** | 6.10 h | Compound tide | Nonlinear shallow water |
| **SA** | 8766 h | Annual | Seasonal SSH variation |
| **SSA** | 4383 h | Semi-annual | Semi-annual SSH |

**Minimum constituent set for Indonesian Seas:** M2, S2, N2, K2, K1, O1, P1, Q1 (8 primary).

**Extended set for shallow shelf applications:** add M4, MS4, MN4, MM, MF.

### Extracting Tidal Forcing with `make_tides.py`

```python
#!/usr/bin/env python3
"""
Extract TPXO10 tidal constituents for Indonesian domain.
"""
from crocotools_py.croco_tides import make_tides

# --- TPXO10 extraction -----------------------------------------------
make_tides(
    grid_file='/data/croco/input/indonesia_grid.nc',

    tidal_model='tpxo10',
    tpxo_dir='/data/TPXO10/',
    # TPXO10 files expected:
    #   h_tpxo10.nc   — elevation amplitude and phase
    #   u_tpxo10.nc   — velocity amplitude and phase
    #   grid_tpxo10   — bathymetry

    constituents=['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1'],

    uv_tides=True,    # Include barotropic tidal currents

    output_file='/data/croco/input/indonesia_tides_tpxo10.nc',
)

# --- FES2022 extraction (alternative, more constituents) -------------
make_tides(
    grid_file='/data/croco/input/indonesia_grid.nc',

    tidal_model='fes2022',
    fes_dir='/data/FES2022/ocean_tide/',

    constituents=['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1',
                  'M4', 'MS4', 'MN4', 'MM', 'MF'],

    uv_tides=True,

    output_file='/data/croco/input/indonesia_tides_fes2022.nc',
)
```

### Tidal CPP Flags in CROCO

```c
#define SSH_TIDES              /* Tidal elevation at boundaries */
#define UBTIDES                /* Tidal barotropic velocity at boundaries */
#define TIDE_GENERATING_FORCES /* Astronomical tidal potential body force */
#define RAMP_TIDES             /* Cosine ramp-in to avoid tidal shock */
/* #define ADD_FSOBC */        /* Also add tides to boundary free-surface (optional) */
/* #define ADD_M2OBC */        /* Also add tidal currents to 2D boundary BC */
```

Note: CROCO uses `UBTIDES` (not `UV_TIDES` as in ROMS-Rutgers) for barotropic tidal currents at boundaries.

### Tidal ramping in `croco.in`

```
RAMP_TIME == 2.0d0   ! Ramp tidal forcing over 2 days (days)
```

Prevents tidal shock at model startup. Recommended ramp time = 2–5 days.

### Post-Processing: Tidal Analysis with `utide`

```python
#!/usr/bin/env python3
"""
Tidal harmonic analysis of CROCO output using utide.
"""
import xarray as xr
import numpy as np
import utide

# Load CROCO station output or history file SSH
ds = xr.open_dataset('/data/croco/output/indonesia_sta.nc')
time = ds['ocean_time'].values  # datetime64 array
zeta = ds['zeta'].isel(station=0).values  # SSH at station 0 (Jakarta Bay)

# Harmonic analysis — Indonesian tidal constituents
coef = utide.solve(
    time.astype('datetime64[s]').astype(float) / 86400,  # days since epoch
    zeta,
    lat=-6.1,              # Station latitude (Jakarta Bay)
    method='ols',          # Ordinary least squares
    conf_int='linear',
    constit=['M2', 'S2', 'K1', 'O1', 'N2', 'K2', 'P1', 'Q1', 'M4', 'MS4'],
    verbose=False,
)

# Print amplitudes and phases
for i, name in enumerate(coef.name):
    print(f"{name:6s}: A = {coef.A[i]:.4f} m,  g = {coef.g[i]:.1f}°")

# Reconstruct tidal signal
tide_recon = utide.reconstruct(
    time.astype('datetime64[s]').astype(float) / 86400,
    coef,
)

# Compute tidal range and M2 form factor
M2_amp = coef.A[coef.name == 'M2'][0]
S2_amp = coef.A[coef.name == 'S2'][0]
K1_amp = coef.A[coef.name == 'K1'][0]
O1_amp = coef.A[coef.name == 'O1'][0]
form_factor = (K1_amp + O1_amp) / (M2_amp + S2_amp)
print(f"\nForm factor F = {form_factor:.2f}")
print("F < 0.25: semidiurnal | 0.25–1.5: mixed | > 1.5: diurnal")
# Indonesian Seas: mostly mixed (0.5–1.5); Java Sea: diurnal dominant
```

---

## 13. Running CROCO

### Pre-Run Checklist

```
[ ] Grid file exists and rx0 < 0.20, rx1 < 12
[ ] Forcing file covers simulation period
[ ] Boundary file covers simulation period + 1 month buffer
[ ] Tidal file extracted for target constituents
[ ] Initial condition file matches vertical grid parameters
[ ] croco.in: DT, NDTFAST, NTIMES, NtileI, NtileJ set correctly
[ ] croco.in: all file paths correct
[ ] AGRIF_FixedGrids.in present if using AGRIF
[ ] Child grid files named with correct prefix (1_, 2_, ...)
[ ] croco.exe compiled with consistent CPP flags and dimensions
[ ] MPI tasks = NtileI × NtileJ
```

### Timestep Stability Criteria

**Barotropic CFL (must be < 1):**
```
CFL_bt = dtfast × sqrt(g × h_max) / dx_min < 1

For Indonesian domain:
  h_max = 5500 m, dx_min = 3000 m, g = 9.81 m/s²
  sqrt(9.81 × 5500) = 232 m/s
  dtfast_max = 3000 / 232 ≈ 12.9 s

With dt = 120 s and NDTFAST = 30: dtfast = 4 s  ✓ (safe)
```

**Baroclinic CFL:**
```
CFL_bc = dt × c_internal / dx_min < 1

c_internal ≈ N_BV × H / π ≈ 0.01 s⁻¹ × 500 m / π ≈ 1.6 m/s
dt_max = 3000 / 1.6 ≈ 1875 s

With dt = 120 s: CFL_bc ≈ 0.064  ✓ (very safe)
```

### Typical Time Steps by Resolution

| Grid spacing | Typical dt (s) | NDTFAST | dtfast (s) |
|---|---|---|---|
| 9 km (Indonesian parent) | 400–600 | 40 | 10–15 |
| 3 km (Indonesian regional) | 120–180 | 30 | 4–6 |
| 1 km (strait child, r=3) | 40–60 | 30 | 1.3–2 |
| 600 m (strait child, r=5) | 24–36 | 30 | 0.8–1.2 |
| 300 m (NH mode) | 10–20 | 30 | 0.3–0.7 |

### SLURM Job Submission Script

```bash
#!/bin/bash
#SBATCH --job-name=croco_indonesia
#SBATCH --partition=compute
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=18         # 18 × 10 = 180 MPI tasks total
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00              # 72-hour wall time limit
#SBATCH --mem=240G
#SBATCH --output=croco_%j.log
#SBATCH --error=croco_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=modeler@bmkg.go.id

# --- Environment ----------------------------------------------------------
module load intel/2023.2 intelmpi/2023.2 netcdf-fortran/4.6.1-intel hdf5/1.14.2-intel

export I_MPI_FABRICS=shm:ofi
export I_MPI_OFI_PROVIDER=psm2
export OMP_NUM_THREADS=1

# --- Run directory --------------------------------------------------------
RUNDIR=/scratch/user/indonesia_2020
cd $RUNDIR

# --- Verify input files ---------------------------------------------------
for f in indonesia_grid.nc indonesia_frc_2020.nc indonesia_bry_2020.nc \
         indonesia_tides_tpxo10.nc indonesia_ini_20200101.nc; do
    [ -f $f ] || { echo "Missing: $f"; exit 1; }
done

# --- Execute CROCO -------------------------------------------------------
echo "Starting CROCO at $(date)"
mpirun -np 180 ./croco.exe croco.in > croco_stdout.log 2>&1
EXIT_CODE=$?

echo "CROCO finished at $(date) with exit code $EXIT_CODE"
exit $EXIT_CODE
```

### Restart Procedure

**Initial run:**
```
NRREC  == 0         ! 0 = start from ININAME (initial conditions file)
ININAME == /data/croco/input/indonesia_ini_20200101.nc
```

**Continue from restart:**
```
NRREC  == -1        ! -1 = read last record from restart file
ININAME == /data/croco/output/indonesia_rst.nc
```

**Partial year restart (e.g., restart at day 90 out of 365):**
```bash
# Check restart file records
ncdump -v ocean_time /data/croco/output/indonesia_rst.nc | grep ocean_time

# Set NRREC to the correct record index (0-based in croco.in)
NRREC  == 1         ! Use second record (most recent if LcycleRST=T and 2 records saved)
```

### Walltime Estimates

Performance for the Indonesian 3 km domain (1800×900×40 grid, 180 MPI tasks on Intel Xeon 8368 nodes):

| Simulation period | Estimated walltime |
|---|---|
| 1 month spin-up | ~6–10 hours |
| 1 year production | ~72–96 hours |
| 1 month with AGRIF (2 children, r=5) | ~20–35 hours |

Performance is dominated by I/O and AGRIF exchange frequency. Reduce `NHIS` and `NAVG` for long runs.

---

## 14. Output & Post-Processing

### Output File Types

| Type | Frequency Parameter | Description |
|---|---|---|
| **History** | `NHIS` | Instantaneous 3D snapshots; large files |
| **Averages** | `NAVG` | Time-averaged 3D fields; most useful for analysis |
| **Averages2** | `NAVG2` | Surface-only averages at higher frequency (requires `AVERAGES2` CPP) |
| **Diagnostics** | `NDIA` | Budget terms (requires `DIAGNOSTICS_TS`, `DIAGNOSTICS_UV`) |
| **Restart** | `NRST` | Full model state for cold restart |
| **Stations** | `NSTA` | Time series at specified points |

### Key Output Variables

| Variable | Description | Grid staggering |
|---|---|---|
| `zeta` | Free-surface elevation (m) | 2D, RHO |
| `ubar`, `vbar` | Depth-averaged velocity (m/s) | 2D, U/V |
| `u`, `v` | 3D velocity (m/s) | 3D, U/V |
| `w` | Vertical velocity (m/s) | 3D, W |
| `temp` | Potential temperature (°C) | 3D, RHO |
| `salt` | Salinity (PSU) | 3D, RHO |
| `rho` | Density anomaly (kg/m³) | 3D, RHO |
| `AKv` | Vertical viscosity Kₘ (m²/s) | 3D, W |
| `AKt` | Vertical diffusivity Kₕ (m²/s) | 3D, W |
| `tke` | Turbulent kinetic energy (m²/s²) | 3D, W |
| `gls` | GLS variable ε or ω | 3D, W |
| `sustr`, `svstr` | Wind stress (N/m²) | 2D, U/V |
| `shflux` | Net surface heat flux (W/m²) | 2D, RHO |
| `ssflux` | Net surface salinity flux (m/s PSU) | 2D, RHO |
| `omega` | S-coordinate vertical velocity (m/s) | 3D, W |

### CROCO_PYTOOLS Post-Processing

```python
#!/usr/bin/env python3
"""
CROCO output diagnostics using CROCO_PYTOOLS and xarray.
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from crocotools_py.croco_vinterp import vinterp

# ---- Load output -------------------------------------------------------
avg = xr.open_mfdataset(
    '/data/croco/output/indonesia_avg*.nc',
    combine='by_coords',
)

# ---- Sea Surface Temperature -------------------------------------------
sst = avg['temp'].isel(s_rho=-1)   # top S-level = surface
lon = avg['lon_rho']
lat = avg['lat_rho']

fig, ax = plt.subplots(
    subplot_kw={'projection': ccrs.PlateCarree()},
    figsize=(14, 8),
)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
cf = ax.pcolormesh(
    lon, lat, sst.isel(ocean_time=-1),
    cmap='RdYlBu_r', vmin=26, vmax=32,
    transform=ccrs.PlateCarree(),
)
plt.colorbar(cf, ax=ax, label='SST (°C)', shrink=0.8)
ax.set_title('CROCO Indonesian Seas — Sea Surface Temperature')
ax.set_extent([90, 145, -15, 10])
plt.savefig('croco_sst.png', dpi=150, bbox_inches='tight')

# ---- Vertical interpolation: sigma → z-levels -------------------------
# vinterp converts from S-levels to fixed z-levels
z_levels = np.array([-5, -10, -25, -50, -100, -200, -500, -1000])
temp_z = vinterp(
    ds=avg,
    var='temp',
    z_levels=z_levels,
    time_idx=-1,         # Last time step
)
# temp_z: shape (len(z_levels), Mm+2, Lm+2)

# ---- Mixed Layer Depth -------------------------------------------------
def mld_dt(temp, salt, z, dT=0.2):
    """MLD by temperature threshold criterion (dT from 10 m)."""
    from scipy.interpolate import interp1d
    T10 = np.interp(-10, -z[::-1], temp[::-1])
    idx = np.where(temp[::-1] < T10 - dT)[0]
    if len(idx) == 0:
        return -z[0]  # shallow / well-mixed
    return -z[-(idx[0])]

# ---- ITF Transport (Lombok Strait section) ----------------------------
# Extract Lombok Strait cross-section at ~115.8°E, 8.3°S–9.0°S
lon_target = 115.8
lat_s, lat_n = -9.0, -8.3

# Find nearest xi index for Lombok
xi_lombok = int((lon_target - 90.0) * 36)   # 3 km = 1/36°
eta_s = int((-lat_s - (-15.0)) * 36)         # Convert lat to eta index (note: latmin=-15)
eta_n = int((-lat_n - (-15.0)) * 36)

# ... extract u cross-section and integrate for transport
```

### Converting S-Coordinates to z-Levels

For analysis requiring fixed depths (e.g., model-obs comparison at CTD stations):

```python
from crocotools_py.croco_vinterp import vinterp, depths

# Compute actual depths of S-levels from grid + output
z3d = depths(
    ds=avg,
    grid_ds=xr.open_dataset('/data/croco/input/indonesia_grid.nc'),
    N=40, theta_s=6.0, theta_b=2.0, hc=50.0,
    vtransform=2, vstretching=4,
)
# z3d: (ocean_time, s_rho, eta_rho, xi_rho) — negative downward

# Vertical interpolation to standard levels
temp_z = vinterp(avg['temp'].values, z3d.values,
                 zout=np.array([-5, -10, -20, -50, -100, -200]))
```

### Typical Diagnostics

| Diagnostic | Method | Relevance for Indonesia |
|---|---|---|
| **SSH variability** | `std(zeta)` over time | Tidal range mapping; eddy activity |
| **SST mean/bias** | Compare `temp[:, -1, :]` to OISST | Model warm/cold bias detection |
| **ITF volume transport** | Integrate `u` through Lombok/Ombai/Timor sections | Key ocean circulation metric |
| **Mixed layer depth** | Threshold method (ΔT = 0.2°C from 10 m) | Thermocline depth, upwelling zones |
| **Tidal ellipses** | Rotary analysis of `u`, `v` at each grid point | Tidal current maps for navigation |
| **Eddy kinetic energy** | EKE = 0.5 × (u′² + v′²) | Mesoscale activity in Banda, Flores seas |
| **Tidal harmonics** | `utide.solve` on station or 2D output | Amplitude/phase maps of M2, K1 |
| **Salinity intrusion** | Track 34 PSU isoline in time | ITF freshwater flux, rain seasonality |

---

## 15. Indonesian Waters Application

### Recommended Domain Setup

The Indonesian Seas are a dynamically complex region connecting the Pacific and Indian Oceans through the Indonesian Throughflow (ITF). A comprehensive CROCO setup must resolve:
- The main ITF passages: Lombok Strait (~250 m sill), Ombai Strait (~1250 m sill), Timor Passage (~1890 m)
- The Makassar Strait (primary ITF pathway, ~2900 m deep)
- The Lifamatola Passage (deep overflow, ~1940 m sill)
- The shallow shelf seas: Java Sea (~50 m), Arafura Sea (~80 m), Banda Sea (~5000 m)
- The Sunda Strait (Selat Sunda) connecting the Java Sea to the Indian Ocean

### Parent Domain Configuration

```
Domain:      90°E – 145°E, 15°S – 10°N
Resolution:  ~3 km (1/36°)
Grid:        1980 × 900 RHO-points
Vertical:    N=40 levels (Vtransform=2, Vstretching=4, θs=6, θb=2, hc=50 m)
Bathymetry:  GEBCO 2023 + BIG (Indonesian Naval) where available
Forcing:     ERA5 hourly; GLORYS12 daily BC; TPXO10 tides (8 constituents)
Period:      Minimum 2-year spinup; 5–10 year production
```

### AGRIF Child Grid Configuration for Key Straits

```
AGRIF_FixedGrids.in for Indonesian Seas configuration:

3                          ! 3 level-1 child grids

! Child 1: Sunda Strait (Selat Sunda), 600 m resolution, r=5
! 5.5°S–7.5°S, 104°E–107°E
  504  612   270  342   5   5

! Child 2: Lombok Strait (Selat Lombok), 600 m resolution, r=5
! 9.5°S–7.5°S, 115°E–117°E
  900  972   198  270   5   5

! Child 3: Makassar Strait southern exit, 600 m resolution, r=5
! 5°S–1°S, 116°E–120°E
  936 1080   360  432   5   5
```

### Key Physics for Indonesian Seas

**Essential CPP flags:**
```c
/* Tides — critical for ITF and mixing */
#define SSH_TIDES
#define UBTIDES
#define TIDE_GENERATING_FORCES
#define RAMP_TIDES

/* Nonlinear EOS — essential for T/S complex water masses */
#define NONLIN_EOS

/* GLS k-epsilon — best for tidal mixing and stratified shear */
#define GLS_MIXING
#define CANUTO_A
#define N2S2_HORAVG

/* RSUP3 advection — rotated scheme handles tidal jets well */
#define TS_HADV_RSUP3
#define VADV_ADAPT_IMP

/* Pressure gradient — important near sills */
#define DJ_GRADPS

/* Bulk fluxes */
#define BULK_FLUXES
#define AERO_BULK
```

### Recommended `croco.in` for Indonesian Seas

```
! Indonesian Seas parent grid (3 km)
DT      == 120.0d0     ! 120-second baroclinic step
NDTFAST == 30          ! barotropic dt = 4 s

! Vertical grid optimised for mixed tropical/deep water
Vtransform  == 2
Vstretching == 4
THETA_S == 6.0d0       ! Surface layers: top 10 levels within 50 m
THETA_B == 2.0d0       ! Bottom layers: adequate for BBL in shallow shelf seas
TCLINE  == 50.0d0

! Drag: moderate (sandy/muddy shelf bottom in Java Sea)
RDRG2 == 3.0d-03

! Horizontal mixing: small for 3 km resolution
VISC2 == 5.0d0
TNU2  == 5.0d0  5.0d0

! Nudging: 10-day timescale at boundaries
TNUDG  == 10.0d0 10.0d0
M3NUDG == 10.0d0
ZNUDG  == 10.0d0

! Restart: every 30 days
NRST == 21600     ! 21600 × 120 s = 30 days

! Output: daily averages, 6-hourly surface fields
NHIS  == 720      ! History every 1 day (720 × 120 s)
NAVG  == 720      ! Averages every 1 day
NAVG2 == 180      ! Surface averages every 6 hours
```

### ITF Monitoring Sections

Define transport sections in the station file for real-time monitoring during the run, or compute offline from history/average output:

```python
# ITF transport monitoring — key passages
# Positive transport = into Indian Ocean (westward/southward)

sections = {
    'Lombok':     {'xi': (780, 780),  'eta': (198, 252), 'var': 'u', 'sign': -1},
    'Ombai':      {'xi': (1188, 1260),'eta': (252, 270), 'var': 'v', 'sign': -1},
    'Timor':      {'xi': (1152, 1368),'eta': (180, 198), 'var': 'v', 'sign': -1},
    'Makassar_S': {'xi': (1008, 1008),'eta': (306, 432), 'var': 'u', 'sign': -1},
    'Lifamatola': {'xi': (1170, 1170),'eta': (420, 450), 'var': 'u', 'sign': +1},
}

# Observed mean transports (Sprintall et al. 2009, JGR):
obs_itf = {
    'Lombok':  -2.6,   # Sv
    'Ombai':   -4.9,   # Sv
    'Timor':   -7.5,   # Sv
}
# Total ITF ≈ 15 Sv southward (into Indian Ocean)
```

### Tidal Characteristics of Indonesian Seas

| Region | Dominant tides | Form factor F | Notes |
|---|---|---|---|
| **Banda Sea** | Mixed, mainly M2 | 0.5–1.0 | Internal tides generated at passages |
| **Java Sea** | Diurnal (K1, O1) | 1.5–3.0 | Resonant diurnal basin |
| **Makassar Strait** | Mixed | 0.5–1.5 | M2 dominant; strong K1 in north |
| **Sunda Strait** | Mixed, mainly M2 | 0.3–0.8 | Semidiurnal + M4 harmonics |
| **Lombok Strait** | Mixed, M2 dominant | 0.5–1.2 | Intense M2 generates ISW |
| **Ombai Strait** | Mixed | 0.4–0.9 | Barotropic+baroclinic tidal mixing |
| **Timor Passage** | Mixed, mainly M2 | 0.3–0.7 | Multi-passage system |

### Atmospheric Forcing for BMKG Context

For BMKG operational ocean forecasting:

| Product | Use case | Resolution | Notes |
|---|---|---|---|
| **GFS 0.25°** | Short-range forecast (0–10 days) | 0.25°, 3-hourly | Download via NOMADS; 4×/day updates |
| **ERA5** | Hindcast/reanalysis | 0.25°, hourly | Gold standard; available 5-day lag |
| **WRF (BMKG)** | Dynamically downscaled forcing | ~3 km, hourly | BMKG-WRF operational system |
| **ECMWF IFS** | Medium-range (0–15 days) | 0.1°, 6-hourly | Via Copernicus ECMWF open data |

When using WRF output as CROCO forcing, extra unit conversions are needed:

```python
# WRF → CROCO bulk forcing conversion
import xarray as xr
wrf = xr.open_dataset('/data/WRF/wrfout_d01_2020-01-01_00:00:00')

# Wind: WRF outputs U10, V10 in m/s (already correct)
# Temperature: T2 in Kelvin → Celsius
tair = wrf['T2'] - 273.15

# Humidity: Q2 (kg/kg mixing ratio) → relative humidity
# using PSFC (surface pressure in Pa) and T2
from metpy.calc import relative_humidity_from_mixing_ratio
rhum = relative_humidity_from_mixing_ratio(wrf['PSFC'], wrf['T2'], wrf['Q2'])

# Radiation: SWDOWN, GLW in W/m² (already correct for CROCO)
```

### Published CROCO Studies in Indonesian Seas

1. **Nugroho et al. (2018)** — CROCO with explicit tidal forcing for Banda Sea internal tides. High-resolution (2 km) AGRIF child in Lombok Strait. (*Frontiers in Marine Science*)
2. **Tranchant et al. (2019)** — CROCO-PISCES biogeochemical simulation for the Eastern Indian Ocean including Indonesian Seas, highlighting tidal effects on nutrient supply. (*Biogeosciences*)
3. **Corvianawatie et al. (2022)** — ITF variability in CROCO hindcast 1993–2019 forced by ERA5. Comparison with INSTANT observations. (*JGR Oceans*)

---

## 16. Troubleshooting & Tips

### Common Errors and Solutions

| Problem | Likely Cause | Solution |
|---|---|---|
| **NaN blowup (early run)** | Timestep too large; bad IC | Halve `DT` and `NDTFAST`; check initial condition file for anomalous values |
| **NaN blowup (after spin-up)** | Tidal shock; bathymetry spike; BC inconsistency | Check `RAMP_TIME`; inspect rx0/rx1 maps at blowup location; increase `VISC2` |
| **CFL exceeded** | `DT` too large for grid and/or tidal currents | Reduce `DT`; check that `NtileI × NtileJ` = total MPI tasks |
| **Pressure gradient error** | rx0 or rx1 too large | Re-smooth bathymetry; use `DJ_GRADPS` instead of `WJ_GRADP` |
| **Spurious deep currents** | Poor vertical stretching over steep sills | Check `hc` < minimum depth; use Vstretching=4; check Vtransform=2 |
| **AGRIF instability at nest boundary** | Resolution jump too large; grid misalignment | Verify nest index alignment; use r=3 not r=5 at first attempt; add buffer |
| **AGRIF child blows up, parent stable** | Child timestep too large; child BC mismatch | Check child `DT`; verify child grid file depth consistency with parent |
| **Negative depths in S-levels** | `hc` > h_min | Set `TCLINE` < minimum bathymetry depth |
| **NetCDF I/O error on restart** | Time variable overflow or mismatch | Check time units in all forcing files (must match reference in croco.in) |
| **Spurious SSH oscillation** | Tidal ramp too short | Increase `RAMP_TIME` to 5 days |
| **Salinity drift** | Boundary nudging too weak | Reduce `TNUDG` to 5 days near open boundaries |
| **MPI decomposition failure** | `NtileI × NtileJ ≠ np` | Ensure `mpirun -np` matches `NtileI × NtileJ` in croco.in |
| **Very slow AGRIF exchange** | MPI halo too small | Add `-DMPI_ISLAND` CPP flag for island MPI topology |

### Bathymetry and Stiffness Tips

```python
#!/usr/bin/env python3
"""Diagnose bathymetry stiffness problems in CROCO grid."""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ds = xr.open_dataset('/data/croco/input/indonesia_grid.nc')
h = ds['h'].values
mask = ds['mask_rho'].values

# --- rx0 map ---
h_xi_plus  = np.roll(h, -1, axis=1)
h_xi_minus = np.roll(h,  1, axis=1)
rx0_xi  = np.abs(h - h_xi_plus)  / (h + h_xi_plus)
rx0_eta = np.abs(h - np.roll(h, -1, axis=0)) / (h + np.roll(h, -1, axis=0))
rx0 = np.maximum(rx0_xi, rx0_eta) * mask

print(f"rx0 max = {rx0.max():.3f} at {np.unravel_index(rx0.argmax(), rx0.shape)}")
print(f"rx0 > 0.20: {(rx0 > 0.20).sum()} cells")
print(f"rx0 > 0.30: {(rx0 > 0.30).sum()} cells")

# Plot rx0 hotspots
fig, ax = plt.subplots(figsize=(14, 8))
p = ax.pcolormesh(ds['lon_rho'], ds['lat_rho'], rx0,
                  cmap='hot_r', vmin=0, vmax=0.4)
plt.colorbar(p, label='rx0')
ax.set_title('Bathymetry stiffness rx0 (target < 0.20)')
plt.savefig('rx0_map.png', dpi=150, bbox_inches='tight')
```

### MPI Tile Decomposition

```
Good tile shapes (approximately square tiles):
  1800×900 grid:
    18×10 = 180 tasks → tiles: 100×90 cells each ✓ (square-ish)
    20×9  = 180 tasks → tiles: 90×100 cells each ✓
    15×12 = 180 tasks → tiles: 120×75 cells each (slightly elongated)

Bad:
    180×1 = 180 tasks → tiles: 10×900 cells each ✗ (very elongated)
    1×180 = 180 tasks → tiles: 1800×5 cells each ✗ (extreme elongation)

Rule: max(tile_xi, tile_eta) / min(tile_xi, tile_eta) < 3
Minimum interior points per tile dimension: ≥ 10 (ghost cells excluded)
```

### AGRIF Debugging Tips

1. **Start simple:** Test parent alone first. Add one child at r=3 before using r=5.
2. **Check grid alignment:** Child `istart` must coincide with a parent RHO-point. Off-by-one errors cause subtle mass leaks.
3. **Monitor child energy:** The child's domain-averaged KE should match the parent in smoothly varying regions. Sudden jumps indicate BC problems.
4. **Use `AGRIF_1WAY` for debugging:** Disabling feedback simplifies diagnosis. If parent is stable in 1-way mode but blows up in 2-way, the child flux correction is the issue.
5. **Output AGRIF boundaries:** Add `-DAGRIF_DEBUG` CPP flag to write boundary exchange diagnostics.

### Performance Optimization

| Technique | Expected Gain | Notes |
|---|---|---|
| **Compiler optimization** | 20–40% | `-O3 -xCORE-AVX2` (Intel); `-O3 -march=native` (GNU) |
| **NetCDF-4 compression** | 30–60% smaller output | Requires HDF5; set `USE_NETCDF4` and deflate level |
| **Reduce NHIS output** | Significant I/O reduction | Write daily or 6-hourly averages only |
| **Single-precision output** | 2× smaller files | Acceptable for analysis; not for restart |
| **Parallel I/O (XIOS)** | Significant for large grids | Complex setup; worthwhile for > 500 CPUs |
| **Tile shape optimization** | 5–15% | Square tiles, cache-friendly dimensions |
| **AGRIF load balancing** | Variable | Child domain size should match parent tile load |

### Resting Ocean Test

Before a full simulation, validate the pressure gradient scheme with a resting ocean test:

```bash
# Create initial conditions from climatology (no motion)
# Run 30 days with no wind, no tides, no boundary forcing
# Check if spurious currents remain < 1 cm/s in open ocean
# > 5 cm/s at sills → re-smooth bathymetry

# croco.in for resting test:
NTIMES == 21600    ! 30 days at dt=120s
# Comment out TIDENAME, FRCNAME, BRYNAME
# Set VISC2 = 0.0, TNU2 = 0.0 (worst case for PGE)
# Analyse resulting u, v max after 30 days
```

---

## References

### CROCO Core Papers

- Debreu, L., Marchesiello, P., Penven, P. & Cambon, G. (2012): "Two-way nesting in split-explicit ocean models: Algorithms, implementation and validation." *Ocean Modelling*, 49–50, 1–21. [AGRIF in CROCO]
- Shchepetkin, A.F. & McWilliams, J.C. (2005): "The Regional Oceanic Modeling System (ROMS): A split-explicit, free-surface, topography-following-coordinate oceanic model." *Ocean Modelling*, 9, 347–404. [Time stepping foundation]
- Shchepetkin, A.F. & McWilliams, J.C. (2003): "A method for computing horizontal pressure-gradient force." *JGR*, 108(C3), 3090. [Pressure gradient]
- Auclair, F., Debreu, L., Duran, E., Marchesiello, P., Soufflet, Y. & Bouchette, F. (2018): "Theory and analysis of acoustic-gravity waves in a free-surface ocean." *Ocean Modelling*, 123, 25–40. [CROCO-NH]
- Warner, J.C., Sherwood, C.R., Signell, R.P., Harris, C.K. & Arango, H.G. (2008): "Development of a three-dimensional, regional, coupled wave, current and sediment-transport model." *Computers & Geosciences*, 34, 1284–1306.
- Warner, J.C. et al. (2005): "Performance of four turbulence closure models implemented using a generic length scale method." *Ocean Modelling*, 8, 81–113. [GLS]
- Large, W.G., McWilliams, J.C. & Doney, S.C. (1994): "Oceanic vertical mixing: A review and a model with a nonlocal boundary layer parameterization." *Rev. Geophys.*, 32(4), 363–403. [KPP]

### Indonesian Seas

- Sprintall, J., Wijffels, S.E., Molcard, R. & Jaya, I. (2009): "Direct estimates of the Indonesian Throughflow entering the Indian Ocean: 2004–2006." *JGR Oceans*, 114, C07001.
- Nugroho, D. et al. (2017): "Modelling explicit tides in the Indonesian seas: An important process for surface sea water properties." *Marine Pollution Bulletin*, 131, 7–18.
- Robertson, R. & Ffield, A. (2008): "Baroclinic tides in the Indonesian seas: Tidal fields and comparisons to observations." *JGR Oceans*, 113, C07031.
- Nagai, T. & Hibiya, T. (2015): "Internal tides and associated vertical mixing in the Indonesian Archipelago." *JGR Oceans*, 120, 3373–3390.
- Gordon, A.L. (2005): "Oceanography of the Indonesian Seas and their throughflow." *Oceanography*, 18(4), 14–27.
- Ffield, A. & Gordon, A.L. (1996): "Tidal mixing signatures in the Indonesian seas." *J. Phys. Oceanogr.*, 26, 1924–1937.

### Tidal Models

- Egbert, G.D. & Erofeeva, S.Y. (2002): "Efficient inverse modeling of barotropic ocean tides." *J. Atmos. Oceanic Technol.*, 19, 183–204. [TPXO]
- Lyard, F.H. et al. (2021): "FES2014 global ocean tide atlas: design and performance." *Ocean Science*, 17, 615–649.

### Online Resources

- [CROCO Website](https://www.croco-ocean.org)
- [CROCO Documentation (GitLab Pages)](https://croco-ocean.gitlabpages.inria.fr/croco_doc/)
- [CROCO GitLab Repository](https://gitlab.inria.fr/croco-ocean/croco)
- [CROCO Forum](https://forum.croco-ocean.org)
- [CROCO_PYTOOLS GitLab](https://gitlab.inria.fr/croco-ocean/croco_pytools)
- [TPXO Tidal Models (OSU)](https://www.tpxo.net)
- [FES2022 / Copernicus Marine](https://data.marine.copernicus.eu)
- [ERA5 Climate Data Store (Copernicus)](https://cds.climate.copernicus.eu)
- [GEBCO Bathymetry](https://www.gebco.net)
- [GLORYS12 (Copernicus Marine)](https://data.marine.copernicus.eu)
- [utide Python tidal analysis](https://github.com/wesleybowman/UTide)
- [xcroco Python post-processing](https://github.com/jaard/xcroco)
- [AGRIF Library (INRIA)](http://agrif.imag.fr)
- [BMKG Ocean Division](https://www.bmkg.go.id/cuaca-maritim/)
