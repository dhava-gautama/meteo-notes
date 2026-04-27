# SCHISM: Complete Guide from Grid Setup to Visualization

> A practical guide covering the full SCHISM (Semi-implicit Cross-scale Hydroscience Integrated System Model) workflow for coastal, estuarine, and regional ocean modeling.
> Covers **SCHISM v5.12**, with emphasis on Indonesian waters and BMKG operational applications.

---

## Table of Contents

1. [Overview & Key Concepts](#1-overview--key-concepts)
2. [SCHISM vs Other Coastal Models](#2-schism-vs-other-coastal-models)
3. [System Requirements & Installation](#3-system-requirements--installation)
4. [Grid Generation](#4-grid-generation)
5. [Forcing Data & Preprocessing](#5-forcing-data--preprocessing)
6. [Tidal Forcing (bctides.in)](#6-tidal-forcing-bctidesin)
7. [param.nml: Master Configuration](#7-paramnml-master-configuration)
8. [Vertical Coordinate (LSC²)](#8-vertical-coordinate-lsc)
9. [Physics Options](#9-physics-options)
10. [Non-Hydrostatic Mode](#10-non-hydrostatic-mode)
11. [Wave Coupling (WWM)](#11-wave-coupling-wwm)
12. [Sediment Transport (SED3D)](#12-sediment-transport-sed3d)
13. [Running SCHISM](#13-running-schism)
14. [Output & Post-Processing](#14-output--post-processing)
15. [Indonesian Waters Application](#15-indonesian-waters-application)
16. [Troubleshooting & Tips](#16-troubleshooting--tips)

---

## 1. Overview & Key Concepts

### What Is SCHISM?

SCHISM (Semi-implicit Cross-scale Hydroscience Integrated System Model) is a 3D baroclinic circulation model built on **unstructured triangular/quadrilateral grids**, capable of seamlessly resolving scales from continental shelves down to estuarine channels and tidal flats in a single domain. It is the direct successor to SELFE (Estuarine and Coastal Ocean Model), developed primarily by Joseph Zhang and colleagues at the Virginia Institute of Marine Science (VIMS).

The "semi-implicit" in the name refers to its **time-stepping philosophy**: the surface gravity wave terms and Coriolis are treated implicitly, removing the CFL stability constraint on the barotropic mode. This allows SCHISM to run with a single, large time step across the whole unstructured grid — even where resolution is locally fine — without splitting into fast/slow modes like structured-grid models (ROMS, MOM6).

SCHISM is open-source (Apache 2.0), actively developed at VIMS and by an international community, and used operationally by CMEMS (Copernicus), NOAA, and regional agencies worldwide.

### Evolution from SELFE

| Generation | Name | Key Changes |
|---|---|---|
| 2000–2008 | SELFE v1–v3 | Original semi-implicit unstructured solver, Zhang & Baptista |
| 2008–2014 | SELFE v3.1 | GLS turbulence, 3D baroclinic, initial GOTM coupling |
| 2014 | SCHISM v5.0 | Complete rewrite in Fortran 90/95, modular architecture, LSC² coordinates, CMake build |
| 2016 | SCHISM v5.3 | WWM (wave model) integration, SED3D sediment |
| 2019 | SCHISM v5.8 | Non-hydrostatic extension, ICM biogeochemistry |
| 2022 | SCHISM v5.10 | Quadrilateral elements, GPU acceleration experiments |
| 2024 | SCHISM v5.12 | Enhanced ParMETIS partitioning, improved GOTM coupling, updated pyschism |

**Reference:** Zhang, Y. J., et al. (2016). Seamless cross-scale modeling with SCHISM. *Ocean Modelling*, 102, 64–81.

### Governing Equations

SCHISM solves the **3D baroclinic Navier-Stokes equations** under the **hydrostatic** and **Boussinesq** approximations (with optional non-hydrostatic extension):

**Horizontal momentum:**
```
∂u/∂t + u·∇u − fv = −g ∂η/∂x − (1/ρ₀) ∫_z^η ∂ρ'/∂x dz' + ∂/∂z(κ ∂u/∂z) + F_u
∂v/∂t + u·∇v + fu = −g ∂η/∂y − (1/ρ₀) ∫_z^η ∂ρ'/∂y dz' + ∂/∂z(κ ∂v/∂z) + F_v
```

**Continuity (free surface):**
```
∂η/∂t + ∇·∫_{-h}^{η} u dz = Q_source
```

**Tracer (T, S and passive):**
```
∂C/∂t + u·∇C = ∂/∂z(κ_h ∂C/∂z) + F_C
```

**Equation of state (JMFWG06):**
```
ρ = ρ(T, S, P) — UNESCO-compatible, nonlinear
```

**Semi-implicit time discretization** (θ-method, θ = 0.5 → Crank-Nicolson):
```
[I − θ Δt L] u^{n+1} = u^n + (1-θ) Δt L u^n + Δt RHS_explicit
```

where L contains gravity wave and Coriolis terms. This renders the scheme **unconditionally stable** for surface gravity waves.

### Unstructured Grid Philosophy

Unlike ROMS (structured curvilinear) or MOM6 (structured with unstructured-like nesting), SCHISM uses a **single unstructured triangular (or mixed tri/quad) mesh**:

```
Offshore (coarse, ~10 km):          Coastal (medium, ~1 km):       Strait (fine, ~200 m):

  *           *           *              *       *       *          *  *  *  *
    \       /   \       /                 \     / \     /            \ | / \ |
      \   /       \   /                   \   /   \   /              \|/   \|
        *           *                       * --- *                   *-----*
      /   \       /   \                   /   \   /                  /|\   /|
    /       \   /       \               /     \ /                   / | \ / |
  *           *           *            *       *                   *  *  *  *
```

**Resolution specification** is flexible: smooth gradients from offshore to nearshore driven by depth and coastline curvature. Elements are triangles with optional quadrilaterals in channels.

### Key Features

| Feature | Description |
|---|---|
| **Grid type** | Unstructured triangular/quadrilateral; single mesh from shelf to tidal flat |
| **Time stepping** | Semi-implicit (θ-scheme); no barotropic/baroclinic split; single Δt |
| **Vertical coordinate** | LSC² (Localized Sigma Coordinates with Shaved Cells) — avoids sigma PGE |
| **Tidal flats** | Robust wetting/drying algorithm |
| **Tracer advection** | TVD schemes with mass conservation |
| **Turbulence** | GLS k-ω, k-ε, MY2.5 via native or GOTM interface |
| **Wave coupling** | Built-in WWM (Wind Wave Model) — no external coupler needed |
| **Sediment** | SED3D (cohesive + non-cohesive, bed morphodynamics) |
| **Non-hydrostatic** | Optional NHYD extension for sub-kilometre dynamics |
| **Scalability** | Tested to >10,000 MPI cores; ParMETIS domain decomposition |
| **Biogeochemistry** | ICM (Integrated Compartment Model), CoSINE, fabm |
| **Data assimilation** | Offline nudging; EnKF via external frameworks |

### Workflow Overview

```
[Define Domain & Resolution Targets]
          |
[Generate Unstructured Mesh]
   ├── SMS / OceanMesh2D / QMESH
   ├── hgrid.gr3  (node/element connectivity)
   └── hgrid.ll   (lon/lat version if projected)
          |
[Assign Bathymetry & Boundaries]
   ├── GEBCO 2025 → hgrid.gr3 (depth column)
   └── open boundary node lists
          |
[Generate Vertical Grid]  →  vgrid.in
          |
[Prepare Forcing Files]
   ├── Atmospheric sflux/  (sflux_air_1.nc, sflux_prc_1.nc, sflux_rad_1.nc)
   ├── Tidal  →  bctides.in
   ├── Ocean BC  →  elev2D.th.nc / uv3D.th.nc / TEM_3D.th.nc / SAL_3D.th.nc
   ├── Initial conditions  →  hotstart.nc (or cold start)
   └── Rivers  →  source_sink.in / msource.th
          |
[Edit param.nml]
          |
[Compile pschism.exe]
          |
[Run: mpirun -np N pschism.exe]
          |
[Combine Outputs]
   ├── combine_output11.f90  or  pyschism
   └── schout_*.nc  →  xarray analysis / pyschism plots
```

---

## 2. SCHISM vs Other Coastal Models

| Feature | **SCHISM v5.12** | **FVCOM v4.4** | **ROMS-Rutgers** | **Delft3D-FM** |
|---|---|---|---|---|
| **Grid type** | Unstructured tri/quad | Unstructured tri | Structured curvilinear | Unstructured quad/tri |
| **Time stepping** | Semi-implicit, single Δt | Explicit (mode-split) | Split-explicit (2D/3D) | Semi-implicit |
| **Non-hydrostatic** | Yes (`-DUSE_NHYD`) | No | No (CROCO fork only) | Yes (Delft3D-FLOW) |
| **Sediment transport** | SED3D (coh + non-coh) | SED (limited) | CSTMS (COAWST) | Full morphology |
| **Wave coupling** | Built-in WWM or ext. SWAN | Ext. SWAVE/SWAN | Ext. SWAN (COAWST) | Built-in SWAN |
| **Biogeochemistry** | ICM, CoSINE, FABM | FABM, UG-RCA | FENNEL, ECOSIM, BIO | Delft3D-WAQ |
| **Data assimilation** | Nudging; ext. EnKF | Nudging | 4D-Var (IS4DVAR) | Limited |
| **Wetting/drying** | Robust, built-in | Yes | Limited | Robust |
| **Scalability (cores)** | >10,000 (ParMETIS) | ~1,000–2,000 | ~2,000–4,000 | ~2,000 |
| **Language** | Fortran 90/CMake | Fortran 90/Makefile | Fortran 90/Make | Fortran/C++ |
| **License** | Apache 2.0 | LGPL-type | MIT/X | Commercial/open |
| **Community** | VIMS + international | UMASSD/WHOI | myroms.org | Deltares |
| **Latest version** | 5.12 (2024) | 4.4 (2023) | Continuous | 2024 |

### When to Choose SCHISM

- **Complex coastlines:** Indonesian archipelago, deltas, estuaries — unstructured mesh handles irregular geometry natively
- **Multi-scale in one domain:** Shelf + coastal + tidal flat without nesting overhead
- **Tidal flat dynamics:** Wetting/drying is a first-class feature
- **Wave-current interaction:** WWM coupling is compiled-in, no MCT/OASIS needed
- **Internal tides / solitons:** Non-hydrostatic mode for sub-kilometre resolution in straits
- **Open-source end-to-end:** pyschism provides Python preprocessing/postprocessing

---

## 3. System Requirements & Installation

### Dependencies

| Library / Tool | Version | Purpose | Notes |
|---|---|---|---|
| **NetCDF-C** | ≥ 4.7 | I/O | With HDF5 backend for compression |
| **NetCDF-Fortran** | ≥ 4.5 | Fortran I/O bindings | Must match NetCDF-C |
| **MPI** | OpenMPI ≥ 4.0 / Intel MPI | Parallelism | Tested on OpenMPI 4.1, Intel MPI 2021 |
| **CMake** | ≥ 3.12 | Build system | Replaces older Makefile approach |
| **ParMETIS** | ≥ 4.0 | Domain decomposition | Optional but highly recommended for >32 cores |
| **GOTM** | 5.x | Turbulence library | Optional; needed for k-ε/MY2.5 |
| **WWM** | Built-in | Wave model | Compiled as part of SCHISM when `USE_WWM=ON` |
| **Lapack/BLAS** | Any | Linear algebra | System or MKL |
| **Python** | ≥ 3.9 | Pre/post-processing | pyschism, xarray, cartopy |

### Getting the Source

```bash
# Clone SCHISM v5.12 from GitHub
git clone https://github.com/schism-dev/schism.git
cd schism
git checkout v5.12.0          # pin to stable release
git submodule update --init   # WWM and other submodules

# Directory structure:
# schism/
#   src/          — Fortran source (core, modules)
#   src/Core/     — main solver (schism_init.F90, schism_step.F90)
#   src/Hydro/    — hydro routines
#   src/Modules/  — physics modules (SED3D, ICM, WWM interface, NHYD)
#   src/GOTM/     — GOTM turbulence library (if enabled)
#   sample_inputs/ — example param.nml, bctides.in, vgrid.in
#   utility/       — post-processing tools (combine_output11.f90, etc.)
#   pyschism/      — Python package for pre/post-processing
```

### CMake Build Process

```bash
# On a typical HPC (Intel compiler + Intel MPI + MKL):
module load intel/2023.1
module load intel-mpi/2021.9
module load netcdf-fortran/4.6.1
module load hdf5/1.14.0
module load cmake/3.26.4
module load parmetis/4.0.3

cd schism
mkdir build && cd build

cmake ../src \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DCMAKE_BUILD_TYPE=Release \
  -DNETCDF_FORTRAN_DIR=$NETCDF_FORTRAN_HOME \
  -DNETCDF_C_DIR=$NETCDF_HOME \
  -DPARMETIS_DIR=$PARMETIS_HOME \
  -DUSE_GOTM=ON \
  -DUSE_WWM=ON \
  -DUSE_SED=ON \
  -DUSE_ICM=OFF \
  -DUSE_NHYD=OFF \
  2>&1 | tee cmake_config.log

make -j 16 pschism 2>&1 | tee build.log

# Executable produced:
ls -lh bin/pschism.exe
```

**CMake flags reference:**

| Flag | Values | Description |
|---|---|---|
| `USE_WWM` | ON/OFF | Include Wind Wave Model |
| `USE_SED` | ON/OFF | Include SED3D sediment transport |
| `USE_GOTM` | ON/OFF | Use GOTM turbulence library |
| `USE_NHYD` | ON/OFF | Non-hydrostatic pressure extension |
| `USE_ICM` | ON/OFF | Integrated Compartment Model (biogeochemistry) |
| `USE_COSINE` | ON/OFF | CoSINE marine biogeochemistry |
| `USE_FABM` | ON/OFF | FABM biogeochemistry framework |
| `USE_MARSH` | ON/OFF | Marsh/wetland dynamics |
| `CMAKE_BUILD_TYPE` | Release/Debug | Release: full optimization; Debug: no optimization + checks |

### Environment Module Setup on HPC (BMKG Meteorologi HPC)

```bash
# ~/.bashrc or module file for SCHISM
export SCHISM_DIR=/scratch/operational/schism/v5.12
export SCHISM_BIN=$SCHISM_DIR/build/bin
export SCHISM_UTIL=$SCHISM_DIR/utility
export SCHISM_SAMPLE=$SCHISM_DIR/sample_inputs
export PATH=$SCHISM_BIN:$PATH

# Python environment for pyschism
module load python/3.10
conda activate schism_env

# Verify
pschism.exe --version   # Should print SCHISM v5.12.0
python -c "import pyschism; print(pyschism.__version__)"
```

### Building with GNU Compiler (Alternative)

```bash
module load gcc/12.3.0
module load openmpi/4.1.5-gcc12
module load netcdf-fortran/4.6.1-gcc12

cmake ../src \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops" \
  -DNETCDF_FORTRAN_DIR=$NETCDF_FORTRAN_HOME \
  -DNETCDF_C_DIR=$NETCDF_HOME \
  -DUSE_GOTM=ON \
  -DUSE_WWM=ON

make -j 16 pschism
```

---

## 4. Grid Generation

### Philosophy: Unstructured Mesh

SCHISM's grid is stored in `hgrid.gr3` — a plain-text file listing node coordinates (lon, lat, depth) and element connectivity (which nodes form each triangle or quad). Resolution can vary smoothly from ~50 km offshore to ~50 m in narrow straits, all in a single mesh.

**hgrid.gr3 format:**
```
[number_of_elements]  [number_of_nodes]
[node_id]  [x(lon)]  [y(lat)]  [depth(m)]
...
[element_id]  [3 or 4]  [node1]  [node2]  [node3]  [node4]
...
[boundary_block]
```

### Grid Generation Tools

| Tool | Language | License | Notes |
|---|---|---|---|
| **SMS** (Surface Water Modeling System) | GUI/Windows | Commercial | Industry standard; full SCHISM workflow support; v13.1 |
| **OceanMesh2D** | MATLAB/Python | MIT | Open source; depth-adaptive mesh; active development |
| **QMESH** | Python | LGPL | QGIS plugin; polygon-based resolution specification |
| **pyschism** | Python | Apache 2.0 | SCHISM's own Python package; uses gmsh backend |
| **gmsh** | C++/GUI | GPL | General-purpose unstructured mesher; scriptable |

### Resolution Strategy for Indonesian Waters

```
Region                    Target Resolution    Rationale
─────────────────────────────────────────────────────────
Deep Indian/Pacific       5–10 km             Baroclinic Rossby radius ~50 km
Java Sea shelf            2–5 km              Shelf-wave dynamics
Continental shelf edge    1–2 km              Frontal dynamics
Coastal/bay              500 m – 1 km        Wind-driven circulation
Makassar Strait          500 m – 1 km        ITF throughflow
Lombok/Ombai Strait      200–500 m           Internal tide generation
Sunda/Gaspar Strait      200–500 m           Strong tidal currents
Estuary mouths           50–200 m            River plumes, salinity fronts
Tidal flats              50–100 m            Wetting/drying
```

### Minimum Depth (hmin)

Wet nodes must have a minimum depth to avoid numerical singularities in thin water columns:

```
hmin = 0.1 m   (tidal flat / marsh applications)
hmin = 1.0 m   (coastal applications, recommended default)
hmin = 2.0 m   (conservative, less wetting/drying activity)
```

Set in `param.nml` as `h0` (minimum depth for wet/dry). All depths in `hgrid.gr3` shallower than `h0` are still read but treated as dry at initialization.

### Typical Node Counts

| Domain | Nodes | Elements | Typical Use |
|---|---|---|---|
| Indonesian Archipelago (parent, 5–10 km) | 80,000–150,000 | 150,000–280,000 | Basin-scale ITF |
| Indonesian + refinement (to 500 m straits) | 300,000–600,000 | 600,000–1,200,000 | Operational forecast |
| Single strait (Lombok, 200 m) | 50,000–100,000 | 100,000–200,000 | Process study |
| Full domain with fine estuaries | >1,000,000 | >2,000,000 | Research |

### Python Example: Domain Setup with pyschism

```python
from pyschism.mesh import Hgrid
from pyschism.mesh.base import Gr3
import numpy as np

# --- Load an existing hgrid ---
hgrid = Hgrid.open('hgrid.gr3', crs='EPSG:4326')

# --- Basic mesh statistics ---
print(f"Number of nodes:    {hgrid.coords.shape[0]}")
print(f"Number of elements: {len(hgrid.elements.triangles)}")
print(f"Min depth:          {hgrid.values.min():.2f} m")
print(f"Max depth:          {hgrid.values.max():.2f} m")

# --- Plot mesh ---
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

fig, ax = plt.subplots(figsize=(14, 10))
hgrid.make_plot(
    ax=ax,
    show_mesh=True,
    show_bathymetry=True,
    vmin=-5000, vmax=0,
    cmap='Blues_r',
)
ax.set_title('SCHISM Unstructured Grid — Indonesian Waters')
ax.set_xlabel('Longitude (°E)')
ax.set_ylabel('Latitude (°N)')
plt.tight_layout()
plt.savefig('schism_mesh_indonesia.png', dpi=150)

# --- Check element quality ---
from pyschism.mesh.quality import mesh_quality_metrics
quality = mesh_quality_metrics(hgrid)
print(f"Min interior angle: {quality['min_angle']:.1f}°")
print(f"Max skewness:       {quality['max_skewness']:.3f}")
print(f"Elements < 30° min angle: {(quality['angles'] < 30).sum()}")
```

### OceanMesh2D Example (MATLAB/Python)

```matlab
% OceanMesh2D — create mesh for Indonesian Archipelago
addpath('/opt/OceanMesh2D/');

% Define geographic extent
bbox = [95, 145; -15, 10];    % [lon_min lon_max; lat_min lat_max]

% Minimum resolution: 200 m near straits, 10 km offshore
min_el  = 200;    % meters
max_el  = 10000;  % meters

% Load GEBCO 2025 bathymetry
dem = 'gebco_2025_indonesia.nc';

% Define resolution function based on depth gradient
gdat = geodata('dem', dem, 'bbox', bbox, 'h0', min_el);

% Edge-length function: finer where bathymetric gradient is large
fh = edgefx('geodata', gdat, ...
            'fs', 1/100, ...  % 1% of length scale per unit
            'max_el', max_el, ...
            'min_el', min_el);

% Generate mesh
mshopts = meshgen('ef', fh, 'bou', gdat, 'plot_on', 1);
mshopts = mshopts.build;

% Export to SCHISM hgrid.gr3
m = mshopts.grd;
write_schism_hgrid(m, 'hgrid.gr3');
```

---

## 5. Forcing Data & Preprocessing

### 5.1 Atmospheric Forcing (sflux/)

SCHISM reads atmospheric forcing from NetCDF files placed in the `sflux/` subdirectory. Three file types are expected:

| File | Variables | Description |
|---|---|---|
| `sflux_air_1.nc` | `uwind`, `vwind`, `prmsl`, `stmp`, `spfh` | Wind (m/s), pressure (Pa), temp (K), specific humidity |
| `sflux_prc_1.nc` | `prate` | Precipitation rate (kg/m²/s) |
| `sflux_rad_1.nc` | `dlwrf`, `dswrf` | Downward longwave / shortwave radiation (W/m²) |

Filenames are indexed: `sflux_air_1.0001.nc`, `sflux_air_1.0002.nc`, ... with each file covering a fixed time window.

**sflux_inputs.txt** controls the stack:
```
&sflux_inputs
  air_1_relative_weight = 1.
  air_2_relative_weight = 1.
  prc_1_relative_weight = 1.
  rad_1_relative_weight = 1.
  step_air = 10800   ! ERA5 3-hourly = 10800 s
  step_prc = 10800
  step_rad = 10800
/
```

**Preparing ERA5 sflux files with Python:**

```python
import xarray as xr
import numpy as np
from datetime import datetime, timedelta

def era5_to_sflux(era5_u10_file, era5_v10_file, era5_sp_file,
                  era5_t2m_file, era5_d2m_file,
                  output_dir='sflux', chunk_days=10):
    """Convert ERA5 to SCHISM sflux_air format."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    ds_u   = xr.open_dataset(era5_u10_file)['u10']
    ds_v   = xr.open_dataset(era5_v10_file)['v10']
    ds_sp  = xr.open_dataset(era5_sp_file)['sp']      # Pa
    ds_t   = xr.open_dataset(era5_t2m_file)['t2m']    # K
    ds_d   = xr.open_dataset(era5_d2m_file)['d2m']    # K dewpoint

    # Specific humidity from dewpoint (Bolton 1980)
    e_s = 6.112 * np.exp(17.67 * (ds_d - 273.15) / (ds_d - 29.65))  # hPa
    spfh = 0.622 * e_s / (ds_sp / 100.0 - 0.378 * e_s)

    times = ds_u['time'].values
    n_chunks = int(np.ceil(len(times) / (chunk_days * 8)))  # 3-hourly

    for i in range(n_chunks):
        sl = slice(i * chunk_days * 8, (i + 1) * chunk_days * 8)
        t_chunk = times[sl]

        # SCHISM time encoding: days since 2000-01-01 00:00:00
        ref = np.datetime64('2000-01-01T00:00:00')
        t_days = (t_chunk - ref) / np.timedelta64(1, 'D')

        out = xr.Dataset({
            'time':  xr.DataArray(t_days, dims=['time'],
                       attrs={'units': 'days since 2000-01-01 00:00:00 UTC',
                              'calendar': 'standard'}),
            'lon':   ds_u['longitude'].rename({'longitude': 'lon'}),
            'lat':   ds_u['latitude'].rename({'latitude': 'lat'}),
            'uwind': ds_u.isel(time=sl).rename({'longitude': 'nx_grid',
                                                'latitude': 'ny_grid'}),
            'vwind': ds_v.isel(time=sl).rename({'longitude': 'nx_grid',
                                                'latitude': 'ny_grid'}),
            'prmsl': ds_sp.isel(time=sl).rename({'longitude': 'nx_grid',
                                                  'latitude': 'ny_grid'}),
            'stmp':  ds_t.isel(time=sl).rename({'longitude': 'nx_grid',
                                                 'latitude': 'ny_grid'}),
            'spfh':  spfh.isel(time=sl).rename({'longitude': 'nx_grid',
                                                 'latitude': 'ny_grid'}),
        })
        fname = f"{output_dir}/sflux_air_1.{i+1:04d}.nc"
        out.to_netcdf(fname, unlimited_dims=['time'])
        print(f"Written: {fname}")
```

### 5.2 Tidal Forcing

See Section 6 (`bctides.in`) for full detail. The recommended tidal databases for Indonesian waters:

| Database | Resolution | Constituents | Access |
|---|---|---|---|
| **TPXO10** | 1/30° global | 15 | OSU (registration required) |
| **FES2022** | 1/16° | 34 | Copernicus Marine (free) |
| **EOT20** | 1/8° | 17 | DGFI-TUM (free) |
| **NAO.99Jb** | 1/12° | 16 | Useful in East Asian Seas |

### 5.3 Ocean Boundary Conditions

SCHISM requires time-varying open boundary conditions for:
- **Elevation** → `elev2D.th.nc`
- **Barotropic velocity** (optional nudging) → `uv3D.th.nc`
- **Temperature** → `TEM_3D.th.nc`
- **Salinity** → `SAL_3D.th.nc`

**Source products:**

| Product | Resolution | Levels | Period | Notes |
|---|---|---|---|---|
| **GLORYS12V1** | 1/12°, 50 z-levels | 50 | 1993–present | CMEMS; best for hindcast |
| **HYCOM GOFS 3.1** | 1/12°, 40 hybrid | 40 | 1994–present | Near-real-time; forecast |
| **CMEMS Global Forecast** | 1/12°, 50 | 50 | 10-day forecast | Operational boundary |
| **WOA23 climatology** | 0.25° | 57 | Climatology | Cold-start initialization |

**HYCOM/GLORYS to SCHISM boundary conversion:**

```python
from pyschism.forcing.hycom import Hycom
from pyschism.mesh import Hgrid

hgrid = Hgrid.open('hgrid.gr3', crs='EPSG:4326')
vgrid_path = 'vgrid.in'

# GLORYS boundary interpolation via CMEMS
hycom = Hycom(
    hgrid=hgrid,
    vgrid=vgrid_path,
    start_date=datetime(2023, 1, 1),
    end_date=datetime(2023, 3, 31),
    product='GLORYS12V1',       # or 'HYCOM'
    nudge_hgrid=hgrid,          # nudging zone = entire domain
)
hycom.fetch_data()
hycom.write(output_directory='./ocean_bc/', overwrite=True)
# Produces: elev2D.th.nc, TEM_3D.th.nc, SAL_3D.th.nc, uv3D.th.nc
```

### 5.4 Initial Conditions

**Cold start** (from climatology/analysis):

```python
from pyschism.forcing.hycom import DownloadHycom
from pyschism.forcing.ic import Ic

# Download a single GLORYS snapshot for initialization
ic = Ic(
    hgrid='hgrid.gr3',
    vgrid='vgrid.in',
    date=datetime(2023, 1, 1, 0, 0, 0),
    product='GLORYS12V1',
)
ic.write(output_path='./hotstart.nc', overwrite=True)
```

**Hotstart (warm restart)** from a previous run:

```bash
# In param.nml:
ihot = 1            ! 0=cold start, 1=hotstart.nc, 2=hotstart from time 0
hot_nsteps = 0      ! or specific step to restart from
```

### 5.5 River Discharge

Rivers are specified via point source/sink files:

**source_sink.in** — source node list:
```
2                      ! number of sources (rivers)
12345                  ! node ID — Cimanuk River mouth (Java)
67890                  ! node ID — Kapuas River mouth (Kalimantan)
0                      ! number of sinks (0 if none)
```

**msource.th** — time series (time in seconds, then T and S for each source):
```
0.0     28.5  0.0    27.8  0.0
3600.0  28.5  0.0    27.8  0.0
86400.0 28.2  0.0    27.5  0.0
```

**vsource.th** — volume flux (m³/s) per source per time step:
```
0.0      500.0   2500.0
3600.0   510.0   2480.0
86400.0  520.0   2600.0
```

---

## 6. Tidal Forcing (bctides.in)

### Overview

`bctides.in` is a plain-text input file that defines tidal constituents and their amplitudes/phases at each open boundary node. It replaces the need for separate tidal NetCDF files — everything is baked in analytically at runtime.

### Generating bctides.in

The utility `gen_bctides.in` (in `utility/Pre-Processing/`) interpolates from a global tidal model database:

```bash
# Using TPXO10 (recommended):
cd utility/Pre-Processing/
./gen_bctides \
  --grid     ../../run/hgrid.gr3 \
  --tidal-db /data/tpxo10/hf.tpxo10_atlas_30_v2 \
  --outfile  ../../run/bctides.in \
  --ntides   8 \
  --constituents M2 S2 N2 K2 K1 O1 P1 Q1 \
  --cutoff-depth 10.0

# Alternative: pyschism wrapper
python -c "
from pyschism.forcing.tides import Tides
tides = Tides(
    tidal_database='tpxo',
    tpxo_database='/data/tpxo10/hf.tpxo10_atlas_30_v2',
)
tides.use_constituents(['M2','S2','N2','K2','K1','O1','P1','Q1','SA','SSA'])
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3', crs='EPSG:4326')
tides.fetch_data(hgrid)
tides.write('bctides.in')
"
```

### bctides.in File Structure

```
0.1 40             ! tidal potential cutoff depth (m), ntip (number of tidal pot. constituents)
0                  ! tip_dp (not used if ntip > 0 in old format)
8                  ! nbfr — number of tidal boundary frequencies
M2
0.0001405189     0.0   1   1.0   0.0   ! angular freq (rad/s), nodal factor, eq. arg, ...
S2
0.0001454441     0.0   1   1.0   0.0
N2
0.0001378797     0.0   1   1.0   0.0
K2
0.0001458423     0.0   1   1.0   0.0
K1
0.0000729212     0.0   1   1.0   0.0
O1
0.0000675977     0.0   1   1.0   0.0
P1
0.0000725229     0.0   1   1.0   0.0
Q1
0.0000649585     0.0   1   1.0   0.0
1                  ! nope — number of open boundary segments
[NETA]             ! number of nodes on this open boundary
3  0  3  0         ! ibtype flags: elev (3=tidal), flow (0=none), temp (3=relax), salt (3=relax)
[node_id]  [eta_M2_amp]  [eta_M2_pha]  [eta_S2_amp]  ... (for all 8 constituents)
...
```

**ibtype flags** (column 1–4 for each boundary segment):

| Flag | Meaning | Type |
|---|---|---|
| 0 | No forcing | — |
| 1 | Time history from file | elev/flow |
| 2 | Constant value | elev/flow |
| 3 | Tidal harmonic from bctides.in | elev |
| 4 | Tidal + subtidal (from elev2D.th.nc) | elev |
| 5 | Flather radiation condition | flow |
| -1 | Nudging to external 3D data | T/S |

### Tidal Constituents — Indonesian Waters

Indonesian waters exhibit **mixed semidiurnal** tides (both diurnal and semidiurnal significant), driven by resonance in the Java Sea and the geometry of the archipelago.

| Constituent | Period (h) | Physics | Importance in Indonesia |
|---|---|---|---|
| **M2** | 12.42 | Principal lunar semidiurnal | Dominant in Makassar Strait, Flores Sea |
| **S2** | 12.00 | Principal solar semidiurnal | Spring-neap modulation |
| **N2** | 12.66 | Larger lunar elliptic | ~19% of M2 |
| **K2** | 11.97 | Lunisolar semidiurnal | ~12% of M2 |
| **K1** | 23.93 | Lunisolar diurnal | Dominant in Java Sea, Gulf of Thailand |
| **O1** | 25.82 | Principal lunar diurnal | Strong in Makassar, Banda |
| **P1** | 24.07 | Principal solar diurnal | ~30% of K1 |
| **Q1** | 26.87 | Larger lunar elliptic diurnal | ~7% of K1 |
| **SA** | 8766 | Solar annual (seasonal) | Sea level annual cycle |
| **SSA** | 4383 | Solar semi-annual | Semi-annual sea level |
| **M4** | 6.21 | Overtide | Shallow-water nonlinearity in Java Sea |

> **Note:** The mixed semidiurnal character (Form Factor F = (K1+O1)/(M2+S2) > 1.5) across the Java Sea means SCHISM simulations in this region must include both diurnal and semidiurnal constituents. Neglecting K1 underestimates tidal range by ~30–50% in the Java Sea interior.

### Nudging / Relaxation Zones

For open-ocean boundaries where tidal forcing alone is insufficient, a relaxation zone prevents spurious reflection:

```
! In bctides.in, nudging parameters per boundary:
! Format: [relax_in] [relax_out] (days)
! relax_in  = nudging timescale within the sponge zone
! relax_out = nudging timescale at domain interior (usually 0 = no nudging)

! Typical values for Indonesian open boundaries:
! relax_in = 1.0 day (strong nudging at boundary)
! relax_out = 0.0 (no interior nudging)
```

---

## 7. param.nml: Master Configuration

`param.nml` is SCHISM's main configuration file — a Fortran namelist controlling all simulation parameters. It replaces the individual `ocean.in` / `roms.in` approach used in structured models.

### Annotated param.nml

```fortran
&CORE
!--- Time control ---
  rnday       = 90              ! Run duration (days)
  dt          = 120.0           ! Time step (seconds) — semi-implicit; CFL ≈ 1-2 OK
  nspool      = 360             ! Output every nspool steps = 720 min = 12 hours
  ihfskip     = 0               ! Steps to skip in hotstart (0 for fresh hotstart)
  start_year  = 2023            ! Start date
  start_month = 1
  start_day   = 1
  start_hour  = 0.0             ! UTC

!--- Hotstart ---
  ihot        = 0               ! 0=cold start; 1=hotstart.nc; 2=hotstart at t=0
  ibc         = 0               ! 0=baroclinic (3D T/S); 1=barotropic

!--- Domain ---
  ics         = 2               ! 1=Cartesian (m); 2=spherical (lon/lat degrees)
  h0          = 0.1             ! Minimum wet depth (m) for wetting/drying
  h1_bcc      = 50.0            ! Depth above which baroclinic bc applied (m)
  h2_bcc      = 100.0           ! Depth below which baroclinic bc full strength (m)
/

&OPT
!--- Baroclinic physics ---
  ibcc_mean   = 0               ! 0=mean density for PGF; 1=use full rho
  rmaxvel     = 10.0            ! Max velocity (m/s) — capping to detect blow-ups
  velmin_bnd  = 0.0             ! Min normal velocity at boundary (m/s)

!--- Equation of state ---
  ieos_type   = 1               ! 1=linear; 10=nonlinear (JMFWG06, recommended)
  ieos_pres   = 0               ! 0=pressure from z; 1=from p_atm + rho*g*z

!--- Turbulence ---
  itur        = 3               ! 0=const Kz; 1=Pacanowski-Philander; 3=GLS k-kl; 4=GOTM
  dfv0        = 1.e-6           ! Background vertical viscosity (m²/s)
  dfh0        = 1.e-6           ! Background vertical diffusivity (m²/s)
  dfv1        = 1.e-4           ! Max vertical viscosity (m²/s)
  dfh1        = 1.e-4           ! Max vertical diffusivity (m²/s)
  mid         = 'MY25'          ! GOTM turbulence model: 'KPP','MY25','keps','komega'

!--- Horizontal diffusion ---
  msc2        = 0               ! 0=no harmonic diffusion; 1=Smagorinsky
  mdc2        = 0               ! Tracer horizontal diffusion (0=none)
  huu         = 0.025           ! Smagorinsky coefficient (if msc2=1)
  hur         = 0.025           ! Smagorinsky coeff for tracers

!--- Wetting/drying ---
  iwet_dry    = 1               ! 1=enable wetting/drying
  do_wwm      = 0               ! 1=couple with WWM wave model

!--- Advection ---
  nadv        = 1               ! 0=no advection; 1=Eulerian upwind; 2=ELM
  dtb_max     = 120.0           ! Max sub-timestep for TVD tracer advection
  thetai      = 0.9             ! Implicitness for vertical advection (0=explicit; 1=implicit)

!--- Bottom friction ---
  nchi        = -1              ! -1=drag law via z0b; 0=const Cd; >0=Manning
  dzb_min     = 0.5             ! Min bottom boundary layer thickness (m)
  z0b_min     = 0.01e-2         ! Min roughness height (m) = 0.0001 m

!--- Coriolis ---
  ncor        = 1               ! 1=Coriolis ON (always for ocean)

!--- Atmospheric forcing ---
  nws         = 2               ! 0=no wind; 1=constant wind; 2=sflux files; 3=user-defined
  wtiminc     = 10800.0         ! sflux time increment (s) — must match sflux files
  nrampwind   = 1               ! Days to ramp wind from 0 to full strength
  iwindoff    = 0               ! 0=use wind stress from sflux; 1=turn off wind

!--- Precipitation/evaporation ---
  iprecip     = 1               ! 1=include precipitation from sflux_prc

!--- Heat flux ---
  ihconsv     = 1               ! 1=heat conservation (solar penetration)
  isconsv     = 1               ! 1=salinity conservation (evap-precip)

!--- Solar radiation penetration (Jerlov water types) ---
  nramp_elev  = 1               ! Days to ramp tidal elevation from 0
  nu_sum_min  = 1.e-4           ! Min sum of eddy viscosities
/

&SCHOUT
!--- Output control ---
  iout_sta    = 0               ! Station output (0=off; 1=on, reads station.in)
  nspool_sta  = 60              ! Station output every N steps

!--- 2D scalar output ---
  iof_hydro(1)  = 1             ! elevation (eta)
  iof_hydro(2)  = 0             ! air pressure
  iof_hydro(3)  = 0             ! wind speed x
  iof_hydro(4)  = 0             ! wind speed y
  iof_hydro(5)  = 1             ! depth-averaged u velocity
  iof_hydro(6)  = 1             ! depth-averaged v velocity

!--- 3D vector/scalar output ---
  iof_hydro(14) = 1             ! u velocity (3D)
  iof_hydro(15) = 1             ! v velocity (3D)
  iof_hydro(16) = 0             ! w velocity (3D)
  iof_hydro(17) = 1             ! temperature (3D)
  iof_hydro(18) = 1             ! salinity (3D)
  iof_hydro(19) = 0             ! density (3D)
  iof_hydro(20) = 1             ! vertical diffusivity Kh
  iof_hydro(21) = 0             ! vertical viscosity Km
  iof_hydro(22) = 0             ! Richardson number
  iof_hydro(23) = 0             ! bottom stress x
  iof_hydro(24) = 0             ! bottom stress y

!--- Output format ---
  !  Outputs go to outputs/ directory
  !  File naming: schout_[var]_[process].nc.0000.nc, then combined
/
```

### Key Parameter Reference Table

| Parameter | Typical Value | Description |
|---|---|---|
| `dt` | 60–300 s | Time step. Larger allowed than structured models because gravity waves are implicit |
| `rnday` | 30–365 | Simulation duration in days |
| `h0` | 0.1–1.0 m | Minimum wet depth (wetting/drying threshold) |
| `ibc` | 0 | 0=baroclinic 3D run; 1=barotropic (2D only) |
| `ics` | 2 | 2=spherical coordinates (always for ocean) |
| `ieos_type` | 10 | 10=nonlinear JMFWG06 equation of state |
| `itur` | 4 | 4=GOTM turbulence; 3=native GLS |
| `nws` | 2 | 2=sflux files for atmospheric forcing |
| `msc2` | 0 or 1 | 1=Smagorinsky horizontal viscosity |
| `nchi` | -1 | -1=Manning/roughness-based bottom friction |
| `rmaxvel` | 10.0 m/s | Max velocity cap (blow-up detector) |
| `ihot` | 0 or 1 | 0=cold start; 1=restart from hotstart.nc |
| `nspool` | 360 | Output frequency in time steps |
| `iof_hydro` | array | Which variables to output (0=off; 1=on) |

---

## 8. Vertical Coordinate (LSC²)

### Why LSC² — Localized Sigma Coordinates with Shaved Cells

Standard sigma (terrain-following) coordinates generate **pressure gradient errors (PGE)** over steep bathymetry because the coordinate surfaces are not horizontal. SCHISM v5+ uses **LSC²**:

1. Each column of nodes has its own set of sigma-like layers, but the **transition points** between layers are localized to the water column depth.
2. **Shaved cells** at the bottom allow the bottom layer to be trimmed to the actual bathymetric depth, rather than forcing all columns to the same fractional depth.
3. Near-surface layers remain uniform (z-like) for all depths > a transition depth `h_s`, avoiding artificial mixing in deep waters.

```
LSC² Schematic (vertical cross-section):

     Surface: z = η(x,y,t)
    ─────────────────────────────────────────────────
    |  σ-type layers (near surface, fine resolution) |
    |  (follow free surface)                         |
    ─────────────────────────────────────────────────
    |                                                 |
    |  Z-type layers (intermediate depths)            |
    |  (horizontal, z-level like)                     |
    |                                                 |
    ─────────────────────────────────────────────────
    |  σ-type layers near bottom (follow bathymetry) |
    |  [SHAVED CELL: last layer trimmed to seabed]   |
     ╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲
         Bathymetry
```

### vgrid.in Format

```
ivcor                  ! 1 = LSC²; 2 = sigma (not recommended for ocean)
nvrt                   ! Total number of vertical levels

! For ivcor=1 (LSC²):
kz                     ! Number of z-level layers (above transition depth)
h_s                    ! Transition depth (m) between sigma and z portions

! Z-level depths (for kz levels):
[level_index]  [depth_m]
...

! S-type sigma parameters (for sigma layers below/above h_s):
h_c                    ! Surface sigma transition depth (m)
theta_b                ! Bottom sigma stretching (0–1)
theta_f                ! Surface sigma stretching (0–10)
```

### Example vgrid.in — 20 Hybrid Levels for Indonesian Waters

```
ivcor = 1           ! LSC²
nvrt  = 20          ! Total levels
kz    = 2           ! 2 pure z-levels at depth
h_s   = 100.0       ! Sigma transition at 100 m depth

! Z-levels (only activated for nodes deeper than h_s):
1     -1000.0       ! Level 1: 1000 m depth
2     -100.0        ! Level 2: 100 m depth (= h_s transition)

! S-type parameters for sigma portion:
h_c     = 10.0      ! Min depth for surface sigma layers (m)
theta_b = 0.5       ! Bottom stretching
theta_f = 5.0       ! Surface stretching

! Sigma level fractions (20 levels from surface to bottom):
! [level_id]  [sigma_value from 0(surface) to 1(bottom)]
1    0.000
2    0.010
3    0.025
4    0.050
5    0.075
6    0.100
7    0.150
8    0.200
9    0.280
10   0.370
11   0.460
12   0.540
13   0.620
14   0.700
15   0.770
16   0.840
17   0.900
18   0.940
19   0.970
20   1.000
```

**Interpretation:** Levels 1–2 are z-levels (activated for deep water), levels 3–20 are sigma-type following bathymetry. The shaved cell ensures the bottom sigma layer always exactly touches the seafloor.

### Generating vgrid.in with pyschism

```python
from pyschism.mesh.vgrid import Vgrid

# Create 30-level LSC² grid for Indonesian domain
vgrid = Vgrid.lsc2(
    hgrid='hgrid.gr3',
    nvrt=30,           # Total vertical levels
    a_vqs0=-0.3,       # Surface sigma layer fraction
    theta=0.7,         # Sigma stretching
    b_vqs0=0.0,        # Bottom sigma fraction
    h_c=5.0,           # Min surface sigma depth
    h_s=80.0,          # Z/sigma transition depth
)
vgrid.write('vgrid.in')

# Visualize vertical grid at a specific node
import matplotlib.pyplot as plt
node_idx = 1234   # example node in Lombok Strait
z_levels = vgrid.get_zlevels(node_idx, eta=0.0)
plt.figure(figsize=(3, 8))
plt.plot(np.zeros(len(z_levels)), z_levels, 'o-')
plt.ylabel('Depth (m)')
plt.title(f'Vertical levels at node {node_idx}')
plt.grid(True)
plt.tight_layout()
plt.savefig('vgrid_node_lombok.png', dpi=100)
```

---

## 9. Physics Options

### 9.1 Turbulence Closure

SCHISM supports multiple turbulence models via the `itur` switch in `param.nml`:

| `itur` | Model | Description | Recommended |
|---|---|---|---|
| 0 | Constant Kz | Fixed vertical diffusivity | Only for testing |
| 1 | Pacanowski-Philander | Richardson-number-based | Equatorial applications |
| 2 | MY2.5 (native) | Mellor-Yamada 2.5-level | Shallow coastal |
| 3 | GLS k-kl (native) | Generic Length Scale | Good default |
| 4 | GOTM | External library; k-ε, k-ω, MY2.5, KPP | Most comprehensive |

**GOTM configuration** (`gotmturb.nml`, placed in run directory):

```fortran
&turbulence
  turb_method = 3          ! 3=second-order turbulence model
  tke_method  = 2          ! 2=dynamic TKE equation
  len_scale_method = 8     ! 8=generic length scale (GLS)
  stab_method = 3          ! 3=Canuto A stability functions
/
&generic
  ! GLS k-omega settings (Umlauf & Burchard 2003)
  gen_m    = 1.5
  gen_n    = -1.0
  gen_p    = 3.0
  cpsi1    = 1.44
  cpsi2    = 1.92
  cpsi3minus = -0.4
  sig_kpsi = 1.0
  sig_psi  = 1.3
/
```

### 9.2 Horizontal Diffusion

| Option | param.nml | Description |
|---|---|---|
| No diffusion | `msc2=0, mdc2=0` | Default; rely on numerical diffusion |
| Smagorinsky (momentum) | `msc2=1, huu=0.025` | Dynamic viscosity based on deformation tensor |
| Smagorinsky (tracers) | `mdc2=1, hur=0.025` | Dynamic diffusivity |
| Constant (biharmonic) | `msc2=2` | Fixed biharmonic viscosity |

**Smagorinsky coefficient guidance:**
- `huu = 0.01–0.05`: Smaller for higher resolution; larger for coarser grids
- For 500 m resolution Indonesian straits: `huu = 0.025` recommended

### 9.3 Equation of State

SCHISM supports two EOS options via `ieos_type`:

```
ieos_type = 1    ! Linear: ρ = ρ₀(1 - α(T-T₀) + β(S-S₀))
                 !   α = 1.7e-4 K⁻¹; β = 7.6e-4 PSU⁻¹
                 !   Only for idealized or testing runs

ieos_type = 10   ! UNESCO JMFWG06 nonlinear EOS (default for ocean)
                 !   Full polynomial in T, S, p
                 !   Required for accurate thermohaline dynamics
```

### 9.4 Wetting/Drying

Wetting/drying is activated by `iwet_dry = 1` in `param.nml`. Key parameters:

```fortran
h0        = 0.10    ! Minimum wet depth (m) — cells shallower go dry
h1_bcc    = 1.0     ! Depth threshold for tracer BC modification
```

**Algorithm:** At each time step, SCHISM checks if water depth at nodes exceeds `h0`. Dry nodes are excluded from the momentum solver; they rejoin when neighboring cell water rises above `h0 + ε`. Mass is conserved through this process.

> **Warning:** Very small `h0` (< 0.01 m) can cause numerical instability in tidal flat applications. For Indonesian mangrove-fronted coasts, `h0 = 0.05–0.10 m` is typically stable with `dt ≤ 60 s`.

### 9.5 Wind Stress Bulk Formulas

When `nws = 2` (sflux files), SCHISM internally computes wind stress using:

```
τ = ρ_air × Cd × |W| × W

Cd = drag coefficient (Large & Yeager 2004 or Kara et al. 2007):
  Cd = (0.142 + 0.0764 × |W| - 0.0023 × |W|²) × 10⁻³  for |W| < 33 m/s
  Cd = 2.33 × 10⁻³                                       for |W| ≥ 33 m/s
```

The COARE 3.5 bulk algorithm is available when compiled with `USE_BULK_FLUX`:

```fortran
! In param.nml:
ihconsv = 1     ! Sensible + latent heat from sflux variables
isconsv = 1     ! Salinity flux (E-P) from sflux
```

### 9.6 Bottom Boundary Condition

```fortran
! Roughness-based drag (recommended: nchi = -1):
nchi    = -1      ! Use z0b roughness height
z0b_min = 1.e-4   ! Minimum roughness length (m) = 0.1 mm (sand)
                  ! Typical values: 1e-4 (sand), 5e-3 (gravel), 1e-2 (biogenic reef)

! Manning n (nchi > 0 = node-varying Manning, read from manning.gr3):
nchi    = 1       ! Read manning.gr3
                  ! Manning n: 0.02 (sand), 0.03 (gravel), 0.05 (mangrove)

! Constant drag coefficient (nchi = 0):
nchi    = 0
cd_a   = 0.0025   ! Cd = 0.0025 (typical)
```

---

## 10. Non-Hydrostatic Mode

### When to Activate

The hydrostatic approximation (dP/dz = -ρg) breaks down when **vertical accelerations are non-negligible** — i.e., when the horizontal scale of motion approaches the water depth:

| Situation | Resolution Threshold | Example |
|---|---|---|
| Internal solitary waves | < 1 km | Lombok Strait, Sulu Sea |
| Near-field tidal jets | < 500 m | Narrow straits, sills |
| Convective overturning | < 100 m | Deep convection (rare coastal) |
| Short surface gravity waves | < 100 m | Normally handled by WWM |

For Indonesian straits where internal solitary waves (ISWs) propagate (e.g., Lombok, Ombai, Lifamatola), non-hydrostatic mode is scientifically important when resolution drops below ~500 m.

### Compilation with `-DUSE_NHYD`

```bash
cmake ../src \
  -DUSE_NHYD=ON \
  -DUSE_WWM=ON \
  -DCMAKE_BUILD_TYPE=Release \
  ...

make -j 16 pschism
```

### Non-Hydrostatic param.nml Keys

```fortran
&NHYD
  inunfl   = 0          ! 0=off; 1=enable non-hydrostatic pressure
  pnc_iter = 2          ! Number of NH pressure iterations (2–5 typical)
  bpbc     = 0.0        ! Coefficient for NH bottom pressure BC
  alpha_nh = 0.5        ! Time implicitness for NH pressure (0.5=Crank-Nicolson)
/
```

**Performance overhead:** Non-hydrostatic adds ~20–50% computational cost (extra Poisson solve per step). Only activate for domains/regions where it matters.

### Internal Solitary Waves — Lombok Strait Example

```
Lombok Strait (~35 km wide, sill depth ~250 m):

  Indian Ocean          Lombok Strait         Flores Sea
    ─────────────────────────────────────────────────────
                    ╔═══════════╗
    Thermocline ════╣  ISW Pkg  ╠════ propagation →
                    ╚═══════════╝
    Sill (~250 m)   ╲╲╲╲╲╲╲╲╲╲╲╲╲
    ─────────────────────────────────────────────────────

ISW characteristics:
  Wavelength:    2–10 km
  Amplitude:     20–80 m vertical displacement
  Phase speed:   ~2–3 m/s
  Needed resolution: < 500 m (ideally 200–300 m)
  SCHISM dt:     30–60 s (with NHYD)
```

**Validation reference:** Susanto, R.D., et al. (2005). Tidal mixing signatures in the Indonesian Seas from high-resolution sea surface temperature data. *GRL*.

---

## 11. Wave Coupling (WWM)

### Overview

WWM (Wind Wave Model) is SCHISM's built-in spectral wave model, derived from WWMII (Roland et al., 2012). It solves the wave action balance equation on the same unstructured mesh as SCHISM — no grid interpolation needed at the coupling interface.

**Compilation:**
```bash
cmake ../src -DUSE_WWM=ON ...
make -j 16 pschism    # Single executable includes both SCHISM + WWM
```

### wwminput.nml — Key Parameters

```fortran
&PROC
  PROCNAME    = 'Indonesian_WWM'
  DIMMODE     = 2               ! 2=2D spectral (Hs, Tm, Dir) output
  LSTEA       = F               ! Steady-state solution (F=transient)
  LQSTEA      = F               ! Quasi-steady-state
  LSPHE       = T               ! Spherical coordinates
  LNAUTIN     = T               ! Nautical convention for wave directions
  BEGTC       = '20230101.000000'
  DELTC       = 120.0           ! WWM time step (s) — match or sub-multiple of dt
  UNITC       = 'SEC'
  ENDTC       = '20230401.000000'
  DMIN        = 0.1             ! Minimum depth (m)
/

&GRID
  IGRIDTYPE   = 3               ! 3=SCHISM unstructured grid
  FILEGRID    = 'hgrid.gr3'
  IWBNDRY     = 1               ! Open boundary wave input
/

&BOUC
  LBCSE       = T               ! Spectral energy BC at open boundaries
  LBCWA       = T               ! JONSWAP parametric BC
  ! Boundary wave parameters (for open boundary nodes):
  BSPBAG      = 3.3             ! Peak enhancement factor γ (JONSWAP)
  BSPCDIR     = 225.0           ! Mean wave direction at boundary (°, nautical)
  BSPHM       = 2.0             ! Hs at boundary (m) — or from file
  BSPTP       = 10.0            ! Peak period (s) at boundary
/

&NUMS
  ICOMP       = 3               ! 3=fully implicit + iterative
  AMETHOD     = 7               ! 7=FVM (finite volume) propagation
  SMETHOD     = 1               ! Source term integration
  DMETHOD     = 1               ! Diffraction approximation
  RTHETA      = 0.5             ! Implicitness for θ-advection
  LSICE       = F               ! Sea ice (off)
  MESNL       = 1               ! 1=DIA nonlinear quad interactions
/

&PHYSICS
  WINDSCALING = 28.0            ! Wind input scaling
  WHITECAP    = 1               ! 1=Komen whitecapping
  BREAKING    = 1               ! 1=depth-induced breaking (Battjes-Janssen)
  ALPHA_BJ    = 1.0             ! Breaking coefficient
  GAMMA_BJ    = 0.6             ! Breaker index
  FRICTION    = 1               ! 1=JONSWAP bottom friction
  JONSWAP_CF  = 0.067           ! JONSWAP friction coefficient
/

&OUTPUT
  IOUTS(1)    = 1               ! Significant wave height (Hs)
  IOUTS(2)    = 1               ! Mean wave period (Tm01)
  IOUTS(3)    = 1               ! Peak wave period (Tp)
  IOUTS(4)    = 1               ! Mean wave direction
  IOUTS(5)    = 1               ! Directional spread
  IOUTS(6)    = 0               ! 2D wave spectra (large output!)
  IOUTS(7)    = 1               ! Bottom orbital velocity
  IOUTS(8)    = 1               ! Wave-induced surface stress
/
```

### SCHISM–WWM Coupling Fields

| Direction | Field | SCHISM → WWM | WWM → SCHISM |
|---|---|---|---|
| **Surface** | Wind speed | Pass U10, V10 | — |
| **Surface** | Current | Pass depth-avg u,v | — |
| **Surface** | Water level | Pass η | — |
| **Bottom** | Depth | Pass total depth | — |
| **Return** | Radiation stress | — | S_xx, S_xy, S_yy (added to momentum) |
| **Return** | Wave mixing | — | TKE production (added to turbulence) |
| **Return** | Stokes drift | — | Us, Vs (Langmuir turbulence) |
| **Return** | Bottom stress | — | τ_wave (bottom orbital enhancement) |

**Coupling in param.nml:**
```fortran
do_wwm       = 1           ! 1=enable WWM coupling
wwminput_dir = './'        ! Directory for wwminput.nml
iwave_form   = 1           ! 1=radiation stress; 2=vortex force
nrampwwm     = 1           ! Ramp-in period (days)
```

### Comparison: WWM vs Standalone SWAN/WW3

| Feature | **WWM (SCHISM built-in)** | **SWAN standalone** | **WW3 standalone** |
|---|---|---|---|
| Grid | Same unstructured as SCHISM | Structured (nested) | Structured or unstructured |
| Coupling | Built-in (zero interpolation) | OASIS/MCT or file-based | OASIS or file-based |
| Physics | WWMII (similar to SWAN) | SWAN physics | WW3 physics (more complete) |
| Ocean coupling | Seamless | Requires coupler | Requires coupler |
| Setup complexity | Low (single executable) | Medium | High |
| Suitable for | Regional coastal + waves | Pure wave climate studies | Offshore + global waves |

---

## 12. Sediment Transport (SED3D)

### Overview

SED3D is SCHISM's built-in sediment transport module supporting:
- Multiple grain-size classes (cohesive and non-cohesive, mixed)
- Bed composition evolution
- Bed morphodynamics (bathymetry update)
- Wave-current combined bed shear stress (when WWM is active)
- Consolidation for cohesive sediments

**Compilation:**
```bash
cmake ../src -DUSE_SED=ON ...
```

### sediment.nml — Key Parameters

```fortran
&SED_OPT
  sed_morph     = .true.    ! Enable bed morphodynamics (bathymetry update)
  sed_morph_fac = 1.0       ! Morphological acceleration factor (1=real time)
  sand_frac     = 0.5       ! Initial sand fraction (0–1)
  mud_frac      = 0.5       ! Initial mud fraction (0–1)

  ! Grain size classes (up to 10 classes):
  ntracer_sed   = 3         ! 3 sediment classes

  ! Class 1: Fine sand (non-cohesive)
  d50(1)        = 200.e-6   ! Median grain size (m) = 200 μm
  rhos(1)       = 2650.0    ! Grain density (kg/m³)
  wsed(1)       = -1.0      ! Settling velocity (m/s; -1=auto-compute from d50)
  tau_ce(1)     = 0.15      ! Critical shear stress for erosion (Pa)
  tau_cd(1)     = 0.08      ! Critical shear stress for deposition (Pa)
  sed_type(1)   = 1         ! 1=non-cohesive; 2=cohesive

  ! Class 2: Coarse silt (non-cohesive)
  d50(2)        = 50.e-6    ! 50 μm
  rhos(2)       = 2650.0
  wsed(2)       = -1.0
  tau_ce(2)     = 0.08
  tau_cd(2)     = 0.05
  sed_type(2)   = 1

  ! Class 3: Fine mud (cohesive)
  d50(3)        = 10.e-6    ! 10 μm
  rhos(3)       = 2650.0
  wsed(3)       = 1.e-4     ! Settling velocity (m/s) — specify for cohesive
  tau_ce(3)     = 0.05      ! Pa
  tau_cd(3)     = 0.02      ! Pa
  sed_type(3)   = 2         ! cohesive

  ! Bed properties:
  bed_thick_init = 5.0      ! Initial bed thickness (m)
  bed_por        = 0.4      ! Bed porosity
  n_bed_layers   = 3        ! Number of bed layers for stratigraphy

  ! Flocculation (cohesive):
  ised_flocmod   = 0        ! 0=no flocculation; 1=Winterwerp floc model
/
```

### Indonesian Delta/Estuary Application

```
Mahakam Delta, East Kalimantan:
  - River input: ~700 m³/s average, 2000 m³/s peak wet season
  - Fine sediment load: ~3 Mt/year (mostly mud class 3)
  - Delta progradation: ~30–50 m/year at active distributaries
  - SCHISM SED3D: 3 sediment classes + morphodynamics
  - Typical mud concentration: 50–500 mg/L in turbidity maximum zone

Ciliwung/Cisadane Estuary, Jakarta Bay:
  - Anthropogenic loading: high suspended sediment, land-based pollution
  - Sedimentation rate: 2–5 cm/year in Jakarta Bay
  - Critical issue: tidal flat accretion vs. land subsidence
  - SCHISM grid: 50–200 m in estuary, 1 km offshore Jakarta Bay
```

### Output Variables (SED3D)

| Variable | Description |
|---|---|
| `sed_dp` | Active layer thickness (m) |
| `sed_tr` | Depth-integrated suspended concentration (kg/m²) per class |
| `bed_thick` | Total bed thickness (m) |
| `bed_frac` | Grain-size fraction per bed layer |
| `bstress` | Total bed shear stress (Pa) |
| `erosion_rate` | Erosion flux (kg/m²/s) |
| `deposition_rate` | Deposition flux (kg/m²/s) |

---

## 13. Running SCHISM

### Required Input Files Checklist

```
run_directory/
├── param.nml           ! Master configuration
├── hgrid.gr3           ! Horizontal grid (nodes, elements, depths)
├── hgrid.ll            ! Same as hgrid.gr3 but in lon/lat (if ICS=1)
├── vgrid.in            ! Vertical coordinate specification
├── bctides.in          ! Tidal boundary conditions
├── elev2D.th.nc        ! Sub-tidal elevation boundary conditions
├── TEM_3D.th.nc        ! Temperature boundary conditions
├── SAL_3D.th.nc        ! Salinity boundary conditions
├── uv3D.th.nc          ! Velocity boundary conditions (optional)
├── hotstart.nc         ! Initial conditions (if ihot=1)
├── sflux/              ! Atmospheric forcing directory
│   ├── sflux_inputs.txt
│   ├── sflux_air_1.0001.nc  ... sflux_air_1.NNNN.nc
│   ├── sflux_prc_1.0001.nc  ... sflux_prc_1.NNNN.nc
│   └── sflux_rad_1.0001.nc  ... sflux_rad_1.NNNN.nc
├── outputs/            ! Directory for model output (must exist)
└── [optional files]
    ├── source_sink.in  ! River locations
    ├── msource.th      ! River T/S time series
    ├── vsource.th      ! River discharge time series
    ├── manning.gr3     ! Spatially varying Manning n
    ├── diffmin.gr3     ! Minimum horizontal diffusivity field
    ├── wwminput.nml    ! WWM wave model config (if USE_WWM)
    ├── sediment.nml    ! SED3D config (if USE_SED)
    └── station.in      ! Time-series station output locations
```

### SLURM Job Script

```bash
#!/bin/bash
#SBATCH --job-name=schism_indonesia
#SBATCH --partition=compute
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32          # 16 × 32 = 512 MPI tasks
#SBATCH --cpus-per-task=1
#SBATCH --mem=0                       # Use all available memory
#SBATCH --time=48:00:00
#SBATCH --output=schism_%j.out
#SBATCH --error=schism_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=operator@bmkg.go.id

module load intel/2023.1
module load intel-mpi/2021.9
module load netcdf-fortran/4.6.1

export I_MPI_HYDRA_BOOTSTRAP=slurm
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export OMP_NUM_THREADS=1

# Run directory
RUNDIR=/scratch/bmkg/schism/runs/indonesia_2023q1
cd $RUNDIR

# Verify input files exist
for f in param.nml hgrid.gr3 vgrid.in bctides.in; do
    [[ -f $f ]] || { echo "Missing: $f"; exit 1; }
done

mkdir -p outputs

echo "Starting SCHISM at $(date)"
echo "Nodes: $SLURM_NNODES × $SLURM_NTASKS_PER_NODE = $SLURM_NTASKS tasks"

mpirun -np $SLURM_NTASKS \
    /scratch/bmkg/schism/v5.12/build/bin/pschism.exe \
    2>&1 | tee schism_run.log

echo "SCHISM completed at $(date)"

# Quick check for completion
if grep -q "Run completed successfully" schism_run.log; then
    echo "SUCCESS: simulation completed normally"
else
    echo "WARNING: check schism_run.log for errors"
    exit 1
fi
```

### Domain Decomposition (ParMETIS)

SCHISM uses ParMETIS for graph-based domain decomposition. It is called **automatically** at model startup when ParMETIS is linked:

```
[Startup] SCHISM reads hgrid.gr3 → builds element adjacency graph
       ↓
[ParMETIS] Partitions graph into N_tasks subdomains
       ↓
[Output] local_to_global_*.in files (one per MPI task)
       ↓
[Each task] reads only its subdomain elements and nodes
       ↓
[Exchange] MPI halo exchange at subdomain boundaries every time step
```

**ParMETIS memory scaling:**
- ~500,000 elements: ~2 GB total partition memory → 32 tasks fine
- ~5,000,000 elements: ~15 GB total → need ≥ 128 tasks with ≥ 128 GB/node

### Restart from Hotstart

```bash
# 1. Convert last output to hotstart format:
cd outputs/
/path/to/schism/utility/Post-Processing/combine_hotstart7 \
    -i schout_0001.nc \      # or use the specific restart step file
    -o ../hotstart.nc \
    -t 77760000              # time in seconds to extract (= 900 days into run)

# 2. Modify param.nml:
#    ihot = 1
#    start_year = 2023, start_month = 4, start_day = 1  (adjusted to restart time)

# 3. Resubmit SLURM script
```

### Typical Walltimes (Indonesian Domain, 512 cores)

| Domain Size | Resolution | Cores | 30-day walltime |
|---|---|---|---|
| 80k nodes (shelf only) | 5–10 km | 64 | ~2 hours |
| 300k nodes (+ fine straits) | 200–5000 m | 256 | ~6 hours |
| 600k nodes (operational) | 200–5000 m | 512 | ~8 hours |
| 1.2M nodes (research) | 100–5000 m | 1024 | ~10 hours |

---

## 14. Output & Post-Processing

### Output File Structure

SCHISM writes output to the `outputs/` directory. In parallel runs, each MPI task writes its own files; these must be combined afterward.

```
outputs/
├── schout_0000_1.nc       ! Process 0, output stack 1 (2D variables)
├── schout_0001_1.nc       ! Process 1, output stack 1
├── ...
├── schout_[N-1]_1.nc      ! Process N-1, output stack 1
├── schout_0000_2.nc       ! Stack 2 (next time block)
├── ...
├── local_to_global_0000   ! ParMETIS partition maps
├── ...
└── combined/              ! After combination
    ├── schout_2023010100.nc   ! Combined: elevation, velocities, T, S
    └── ...
```

**Output variable naming convention:**

| Variable | SCHISM name | Dimensions | Description |
|---|---|---|---|
| Surface elevation | `elev` | (time, node) | Free surface η (m) |
| Depth-avg u | `dahv` | (time, node) | Depth-averaged eastward velocity (m/s) |
| 3D u velocity | `hvel` | (time, node, sigma) | Horizontal velocity vector (m/s) |
| Temperature | `temp` | (time, node, sigma) | Potential temperature (°C) |
| Salinity | `salt` | (time, node, sigma) | Salinity (PSU) |
| Vertical velocity | `w` | (time, node, sigma) | Vertical velocity (m/s) |
| Diffusivity | `diffusivity` | (time, node, sigma) | Vertical eddy diffusivity (m²/s) |
| Viscosity | `viscosity` | (time, node, sigma) | Vertical eddy viscosity (m²/s) |
| Bottom stress | `bottom_stress` | (time, node) | Bottom shear stress (Pa) |
| Wave height (WWM) | `WWM_1` | (time, node) | Significant wave height Hs (m) |

### Combining Outputs

**Tool 1: combine_output11.f90 (Fortran utility)**

```bash
# Compile the combiner (do once):
cd /opt/schism/v5.12/utility/Post-Processing/
gfortran -O2 -o combine_output11 combine_output11.f90 \
    -I$NETCDF_INC -L$NETCDF_LIB -lnetcdff -lnetcdf

# Combine for a specific output block:
cd /scratch/runs/indonesia_2023q1/outputs/
/opt/schism/v5.12/utility/Post-Processing/combine_output11 \
    -i schout_0000_1.nc \   # any one file from stack 1
    -n 512 \                # number of MPI tasks
    -o combined/schout_combined_1.nc

# Loop over all stacks:
for stack in $(seq 1 90); do
    /opt/schism/.../combine_output11 \
        -i schout_0000_${stack}.nc \
        -n 512 \
        -o combined/schout_combined_${stack}.nc
done
```

**Tool 2: pyschism (Python)**

```python
from pyschism.outputs import SchismOutputCombiner

combiner = SchismOutputCombiner(
    run_directory='/scratch/runs/indonesia_2023q1',
    output_directory='outputs/',
    num_processors=512,
)
combiner.combine(
    stacks=list(range(1, 91)),   # stacks 1–90
    out_dir='outputs/combined/',
    overwrite=True,
)
```

### xarray Analysis

```python
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# --- Load combined output ---
ds = xr.open_mfdataset(
    'outputs/combined/schout_combined_*.nc',
    combine='by_coords',
    decode_times=True,
)

# --- Load grid ---
import pandas as pd
hgrid_data = pd.read_csv('hgrid.gr3', skiprows=1, sep=r'\s+',
                          names=['id','x','y','depth'],
                          nrows=int(open('hgrid.gr3').readline().split()[1]))
x = hgrid_data['x'].values
y = hgrid_data['y'].values

# Read element connectivity for triangulation
with open('hgrid.gr3') as f:
    ne, nn = map(int, f.readline().split())
    for _ in range(nn):
        f.readline()  # skip nodes
    triangles = []
    for _ in range(ne):
        parts = f.readline().split()
        if int(parts[1]) == 3:  # triangle
            triangles.append([int(parts[2])-1, int(parts[3])-1, int(parts[4])-1])

triang = tri.Triangulation(x, y, triangles)

# --- Plot sea surface elevation ---
fig, ax = plt.subplots(figsize=(14, 8))
eta = ds['elev'].isel(time=-1).values
tcf = ax.tricontourf(triang, eta, levels=50, cmap='RdBu_r', vmin=-1.5, vmax=1.5)
ax.tricontour(triang, eta, levels=[0], colors='k', linewidths=0.5)
cbar = plt.colorbar(tcf, ax=ax, label='Sea Surface Elevation (m)')
ax.set_xlabel('Longitude (°E)')
ax.set_ylabel('Latitude (°N)')
ax.set_title('SCHISM — Sea Surface Elevation, Indonesian Waters')
plt.tight_layout()
plt.savefig('schism_eta_indonesia.png', dpi=150, bbox_inches='tight')

# --- Time series at a station ---
# Find nearest node to Lombok Strait mooring
from scipy.spatial import KDTree
tree = KDTree(np.column_stack([x, y]))
lombok_lon, lombok_lat = 116.0, -8.5
dist, node_idx = tree.query([lombok_lon, lombok_lat])

eta_ts = ds['elev'].isel(node=node_idx).values
temp_ts = ds['temp'].isel(node=node_idx, nSCHISM_vgrid_layers=-1).values  # surface

time = ds['time'].values
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
ax1.plot(time, eta_ts, 'b-', lw=0.8)
ax1.set_ylabel('Elevation (m)')
ax1.set_title(f'Lombok Strait ({lombok_lon}°E, {lombok_lat}°N)')
ax2.plot(time, temp_ts, 'r-', lw=0.8)
ax2.set_ylabel('SST (°C)')
ax2.set_xlabel('Time (UTC)')
plt.tight_layout()
plt.savefig('schism_lombok_timeseries.png', dpi=100)

# --- T/S Cross-section ---
# Extract a transect across Lombok Strait
transect_nodes = [node_idx - 50 + i for i in range(100)]  # simplified
temp_section = ds['temp'].isel(node=transect_nodes, time=0).values
salt_section = ds['salt'].isel(node=transect_nodes, time=0).values
depth_section = hgrid_data['depth'].values[transect_nodes]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, var, label, cmap in zip(axes,
                                 [temp_section, salt_section],
                                 ['Temperature (°C)', 'Salinity (PSU)'],
                                 ['RdYlBu_r', 'viridis']):
    im = ax.contourf(range(100), np.linspace(0, -1, var.shape[1]),
                     var.T, levels=20, cmap=cmap)
    plt.colorbar(im, ax=ax, label=label)
    ax.set_xlabel('Node index (W→E across Lombok)')
    ax.set_ylabel('Sigma level (surface→bottom)')
plt.suptitle('SCHISM T/S Cross-section: Lombok Strait')
plt.tight_layout()
plt.savefig('schism_lombok_section.png', dpi=100)
```

### Typical Diagnostics to Check

| Diagnostic | Method | Healthy Sign |
|---|---|---|
| Volume conservation | `grep "Total volume" schism_run.log` | Changes < 0.01% per day |
| Energy | `grep "KE" schism_run.log` | No monotonic growth |
| Min depth | Check `elev` doesn't exceed limits | No node permanently above MWL |
| Salinity range | `ds['salt'].min(), .max()` | 0–40 PSU; no negatives |
| Temperature range | `ds['temp'].min(), .max()` | 15–35 °C for Indonesian waters |
| Velocity extremes | `ds['hvel'].max()` | < 5 m/s in straits; < 10 m/s limit |
| Tidal amplitude | Harmonic analysis at coastal stations | Agree with tide gauge within ~10% |

---

## 15. Indonesian Waters Application

### Domain Design for Indonesian Archipelago

The Indonesian Maritime Continent sits at the confluence of the Pacific and Indian Oceans, connected by a complex system of straits that regulate the **Indonesian Throughflow (ITF)** — the only low-latitude exchange between ocean basins. SCHISM's unstructured mesh is ideally suited for this geometry.

```
SCHISM Indonesian Domain — Resolution Map:

  95°E            115°E           135°E          145°E
   |                 |               |               |
10°N ─────────────────────────────────────────────────── 10°N
     |                                                 |
     |    Andaman Sea       South China Sea            |
 5°N |    [5 km]              [5 km]                   | 5°N
     |                      ┌──────────────────────┐   |
     |  Malacca Strait      │ Natuna Sea [2–5 km]  │   |
 0°  │  [200–500 m]         └──────┬───────────────┘   | 0°
     │                      Karimata [500 m]            |
     │         Sunda Strait Java Sea [2 km]  Makassar  │
 5°S │         [200–300 m]  ─────────────── [500 m]    │ 5°S
     │                       Bali Lombok Sape           │
     │  Indian Ocean         [200 m] Straits [200 m]   │
10°S │  [5–10 km]                         Banda Sea    │ 10°S
     │                                    [2–5 km]     │
     │                       Flores Sea    Ombai [200m]│
15°S │                       [1–2 km]     Timor [500m] │ 15°S
     |                                                 |
     ─────────────────────────────────────────────────
  95°E            115°E           135°E          145°E
```

### ITF Pathways and Key Straits

| Strait | Width | Sill Depth | SCHISM Resolution | Key Process |
|---|---|---|---|---|
| **Makassar** | 120 km | 680 m | 500 m – 2 km | Primary ITF conduit (~80% transport) |
| **Lombok** | 35 km | 300 m | 200–500 m | Internal tide generation; ISWs |
| **Ombai** | 35 km | 3200 m | 200–500 m | Deep ITF outflow to Indian Ocean |
| **Timor** | 100 km | 1890 m | 500 m | Large ITF outflow passage |
| **Lifamatola** | 40 km | 1940 m | 500 m | Deep overflow into Banda Sea |
| **Sunda** | 24 km | 20 m | 100–200 m | Wind-driven seasonal reversal |
| **Karimata** | 100 km | 40 m | 500 m | South China Sea ↔ Java Sea exchange |

### Monsoon-Driven Coastal Currents

```
DJF (Dec–Jan–Feb): NW Monsoon
  ─ Java Sea currents: Westward (wind-driven)
  ─ Banda Sea: Cyclonic gyre
  ─ Upwelling: Suppressed on south Java coast
  ─ SST: Warm pool in eastern Indonesia

JJA (Jun–Jul–Aug): SE Monsoon
  ─ Java Sea currents: Eastward
  ─ South Java: Strong upwelling (Ekman divergence)
  ─ SST: SST anomaly −2°C to −4°C off south Jawa, Nusa Tenggara
  ─ ITF: Maximum transport (stronger westward pressure gradient)
```

### Sample bctides.in for Indonesian Tidal Forcing

```
! bctides.in for Indonesian Archipelago
! 10 tidal constituents: M2 S2 N2 K2 K1 O1 P1 Q1 SA SSA
! Open boundaries: West (Indian Ocean), North (Pacific via South China Sea)
!
0.1  5          ! cutoff depth (m), ntip (tidal potential constituents)
5               ! ntip forcing terms
! Each line: period (days), amp (m), phase (°), love_h, love_k
M2  0.517525  0.242  136.8  0.611  0.302
S2  0.500000  0.113  158.4  0.611  0.302
N2  0.527428  0.046  118.6  0.611  0.302
K1  0.997269  0.141   45.2  0.736  0.415
O1  1.075803  0.100   15.8  0.695  0.348
10             ! nbfr — boundary frequencies
M2
0.0001405189  0.0  1  1.0  0.693
S2
0.0001454441  0.0  1  1.0  0.693
N2
0.0001378797  0.0  1  1.0  0.693
K2
0.0001458423  0.0  1  1.0  0.693
K1
0.0000729212  0.0  1  1.0  0.736
O1
0.0000675977  0.0  1  1.0  0.695
P1
0.0000725229  0.0  1  1.0  0.706
Q1
0.0000649585  0.0  1  1.0  0.695
SA
0.0000001140  0.0  1  1.0  1.0
SSA
0.0000002280  0.0  1  1.0  1.0
2              ! nope — 2 open boundary segments
! Boundary 1: West open boundary (Indian Ocean, 1200 nodes)
1200
4  4  4  4     ! ibtype: elev=tidal+subtidal, flow=Flather, temp=relax, salt=relax
[node_1]  [M2_amp] [M2_pha] [S2_amp] [S2_pha] ... [SSA_amp] [SSA_pha]
...
! Boundary 2: North open boundary (South China Sea, 800 nodes)
800
4  4  4  4
[node_1]  [M2_amp] [M2_pha] ...
...
```

### pyschism Domain Setup — Full Example

```python
"""
SCHISM Indonesian Archipelago Setup
Using pyschism v0.10+ and GLORYS12 boundary conditions
"""

from pyschism.mesh import Hgrid
from pyschism.mesh.vgrid import Vgrid
from pyschism.forcing.tides import Tides
from pyschism.forcing.hycom import Hycom
from pyschism.forcing.ic import Ic
from pyschism.forcing.sflux import ERA5
from pyschism.param import Param
from datetime import datetime
import os

# ---- Configuration ----
RUN_DIR    = '/scratch/bmkg/schism/runs/indonesia_Q1_2023'
MESH_FILE  = '/data/schism/grids/indonesia_v2/hgrid.gr3'
TPXO_DIR   = '/data/tpxo10'
ERA5_DIR   = '/data/era5/indonesia_2023'
START_DATE = datetime(2023,  1,  1)
END_DATE   = datetime(2023,  4,  1)

os.makedirs(RUN_DIR, exist_ok=True)
os.makedirs(f'{RUN_DIR}/sflux', exist_ok=True)
os.makedirs(f'{RUN_DIR}/outputs', exist_ok=True)

# ---- 1. Load grid ----
hgrid = Hgrid.open(MESH_FILE, crs='EPSG:4326')
print(f"Grid loaded: {hgrid.coords.shape[0]} nodes")

# ---- 2. Vertical grid (30 LSC² levels) ----
vgrid = Vgrid.lsc2(hgrid=hgrid, nvrt=30, h_s=80.0, theta_f=5.0, theta_b=0.5)
vgrid.write(f'{RUN_DIR}/vgrid.in')

# ---- 3. Tidal forcing (TPXO10) ----
tides = Tides(tidal_database='tpxo10', resource=TPXO_DIR)
tides.use_constituents(['M2','S2','N2','K2','K1','O1','P1','Q1','SA','SSA'])
tides.fetch_data(hgrid)
tides.write(f'{RUN_DIR}/bctides.in')

# ---- 4. Atmospheric forcing (ERA5 sflux) ----
era5 = ERA5(
    hgrid=hgrid,
    start_date=START_DATE,
    end_date=END_DATE,
    pscr=ERA5_DIR,
)
era5.write(output_directory=f'{RUN_DIR}/sflux', overwrite=True)

# ---- 5. Ocean boundary (GLORYS12) ----
hycom = Hycom(
    hgrid=hgrid,
    vgrid=f'{RUN_DIR}/vgrid.in',
    start_date=START_DATE,
    end_date=END_DATE,
    product='GLORYS12V1',
)
hycom.fetch_data()
hycom.write(output_directory=RUN_DIR, overwrite=True)

# ---- 6. Initial conditions ----
ic = Ic(
    hgrid=MESH_FILE,
    vgrid=f'{RUN_DIR}/vgrid.in',
    date=START_DATE,
    product='GLORYS12V1',
)
ic.write(output_path=f'{RUN_DIR}/hotstart.nc', overwrite=True)

# ---- 7. param.nml ----
param = Param()
param.core.rnday   = (END_DATE - START_DATE).days
param.core.dt      = 120.0
param.core.nspool  = 450         # Output every 450 steps = 15 hours
param.core.ics     = 2           # Spherical
param.core.ihot    = 1           # Use hotstart.nc
param.opt.ieos_type = 10         # Nonlinear EOS
param.opt.itur     = 4           # GOTM turbulence
param.opt.nws      = 2           # sflux atmospheric forcing
param.opt.wtiminc  = 10800.0     # ERA5 3-hourly
param.opt.iwet_dry = 1           # Wetting/drying ON
param.opt.do_wwm   = 0           # No waves this run
param.schout.iof_hydro[1]  = 1   # eta
param.schout.iof_hydro[5]  = 1   # depth-avg u
param.schout.iof_hydro[6]  = 1   # depth-avg v
param.schout.iof_hydro[14] = 1   # 3D u
param.schout.iof_hydro[15] = 1   # 3D v
param.schout.iof_hydro[17] = 1   # temperature
param.schout.iof_hydro[18] = 1   # salinity
param.write(f'{RUN_DIR}/param.nml')

print(f"Setup complete. Run directory: {RUN_DIR}")
print(f"Submit with: sbatch {RUN_DIR}/run_schism.slurm")
```

### BMKG Operational Potential

SCHISM is well-suited for BMKG operational marine forecasting:

| Application | SCHISM Capability | Domain | Update Cycle |
|---|---|---|---|
| Storm surge warning | Wetting/drying + tidal + wind | Indonesian coasts | 6-hourly |
| Tidal flooding (rob) | Tidal + SLR + land subsidence | Java north coast | Daily |
| Marine currents (Pelayaran) | 3D currents, SST | National waters | 12-hourly |
| Oil spill trajectory | Passive tracer + Lagrangian | Regional | On-demand |
| Coral bleaching | SST + thermal stress | Coral Triangle | Weekly |
| HAB (harmful algal bloom) | ICM biogeochemistry | Coastal waters | Seasonal |
| ITF monitoring | 3D baroclinic + straits | Full archipelago | Monthly |

**Operational chain example (BMKG):**

```
00 UTC: GFS/ECMWF forecast arrives
  ↓ 01 UTC: ERA5 → sflux conversion (Python, 10 min)
  ↓ 02 UTC: GLORYS/HYCOM boundary update (Python, 20 min)
  ↓ 03 UTC: SCHISM 5-day forecast run (512 cores, 6 hours walltime)
  ↓ 09 UTC: pyschism post-processing + visualization (30 min)
  ↓ 10 UTC: Products pushed to BMKG public portal
```

---

## 16. Troubleshooting & Tips

### Common Errors and Fixes

#### 1. CFL Violation Despite Semi-Implicit Scheme

Although SCHISM's gravity wave terms are implicit (unconditionally stable), the **advective CFL** for momentum and tracers is still explicit:

```
CFL_advective = |u| × dt / Δx  must be < 1

Symptom: "MAXVEL exceeded" in log; velocities blow up in narrow channel
Fix:
  - Reduce dt (try dt/2)
  - Increase Δx in the problematic region (coarsen the mesh there)
  - Reduce rmaxvel temporarily to isolate blow-up location:
      rmaxvel = 5.0  (instead of 10.0)
  - Check bctides.in — phasing error can induce large spurious currents
  - Inspect mesh quality at the blow-up node (often a sliver triangle)
```

#### 2. Wetting/Drying Instability

```
Symptom: Salt/temp oscillation in intertidal zone; negative depths in log
Fix:
  - Increase h0 (e.g., 0.1 → 0.5 m)
  - Reduce dt in areas with frequent wetting/drying
  - Check hgrid.gr3 for extremely shallow (< h0) nodes that should be land
  - Use manning.gr3 with higher roughness in tidal flat areas (n = 0.05–0.10)
  - Verify sflux precipitation doesn't add unrealistic freshwater on tidal flats
```

#### 3. Negative Salinity

```
Symptom: ds['salt'].min() < 0; usually appears near river mouths
Fix:
  - Check vsource.th: is freshwater flux positive (inflow) or negative?
  - Verify msource.th salinity is 0.0 for river sources
  - Reduce dt_advection or use more conservative tracer scheme
  - Enable tracer clipping (not native in SCHISM — post-process):
      salt = np.maximum(salt, 0.0)
  - Check if river node is on a dry cell at low tide
```

#### 4. Instability at Open Boundaries

```
Symptom: Waves propagating inward; SST/SSS spike at boundaries
Fix:
  - Increase nudging strength (reduce relaxation timescale)
  - In bctides.in, use ibtype=4 (tidal + subtidal) not ibtype=3 (tidal only)
  - Extend sponge zone (nudge interior nodes, not just boundary)
  - Check elev2D.th.nc temporal coverage — gaps cause step-changes
  - Verify GLORYS/HYCOM boundary data covers full run period + 1 day
```

#### 5. ParMETIS Memory Error

```
Symptom: "ParMETIS memory allocation failed" at startup
Fix:
  - Increase memory per task: --mem-per-cpu=8G in SLURM
  - Reduce number of tasks (fewer, larger subdomains)
  - Check if mesh has disconnected components (islands with no bridge nodes)
  - Use serial METIS first to verify graph:
      mpmetis hgrid.gr3 N_tasks
```

#### 6. Missing sflux Files

```
Symptom: "sflux_air_1.0000N.nc not found" early termination
Fix:
  - sflux indexing starts at 0001, not 0000
  - Verify sflux_inputs.txt step_air matches file time interval
  - Check file naming: must be sflux_air_1.NNNN.nc (4 digits, zero-padded)
  - Confirm sflux/ directory is in the same directory as param.nml
```

### Mesh Quality Checks

```python
from pyschism.mesh.quality import check_mesh_quality
from pyschism.mesh import Hgrid

hgrid = Hgrid.open('hgrid.gr3', crs='EPSG:4326')
quality = check_mesh_quality(hgrid)

# Recommended thresholds:
print(f"Minimum angle: {quality['min_angle']:.1f}° (target > 30°)")
print(f"Max skewness: {quality['max_skewness']:.3f} (target < 0.85)")
print(f"Max aspect ratio: {quality['max_aspect_ratio']:.2f} (target < 5)")
print(f"Elements < 20° angle: {quality['n_bad_elements']}")

# Find and plot bad elements:
import matplotlib.pyplot as plt
bad_mask = quality['angles'] < 25.0
bad_elements = quality['element_indices'][bad_mask]
print(f"Bad elements (< 25°): {len(bad_elements)}")
# These should be manually refined in SMS or gmsh
```

**Mesh quality rules of thumb:**

| Metric | Minimum Acceptable | Ideal | Action if Violated |
|---|---|---|---|
| Min interior angle | 20° | > 30° | Refine or smooth in mesh editor |
| Max skewness | < 0.90 | < 0.70 | Laplacian smoothing |
| Aspect ratio | < 10 | < 5 | Avoid elongated triangles in flow direction |
| Size gradient | < 30% per element | < 15% | Smooth transition zone |
| Depth gradient | — | < 100% change/element | Avoid cliff bathymetry at single element |

### hmin Tuning

```
hmin (h0) too small → wetting/drying instability, negative volumes
hmin (h0) too large → artificial damping of tidal flat dynamics

Recommended workflow:
  1. Start with h0 = 1.0 m (conservative, stable)
  2. Run 7-day test with observed tidal signal at coastal gauges
  3. Gradually reduce h0 → 0.5 m → 0.2 m → 0.1 m
  4. Check tidal elevation RMSE at tide gauges at each step
  5. Use smallest h0 that keeps simulation stable for full simulation period

For Indonesian mangrove coasts: h0 = 0.05–0.10 m with dt ≤ 60 s
For open shelf/straits: h0 = 0.5–1.0 m with dt ≤ 120–180 s
```

### MPI Scaling Guidance

```
Elements per core:
  Optimal:   3,000–8,000 elements/core (good cache use)
  Acceptable: 1,000–15,000 elements/core
  Too few:   < 500 elements/core (MPI communication overhead dominates)
  Too many:  > 20,000 elements/core (memory pressure)

Examples:
  300k elements: 50–300 cores (optimal: ~100 cores)
  600k elements: 75–600 cores (optimal: ~200 cores)
  2M elements:   250–2000 cores (optimal: ~500–800 cores)

Strong scaling efficiency:
  SCHISM typically shows 70–85% efficiency doubling cores up to ~2,000 cores
  Beyond ~4,000 cores, efficiency drops unless mesh is very large (>5M elements)

Practical BMKG recommendation:
  Operational 600k-element Indonesian domain → 256–512 cores
  Research 2M-element domain → 512–1024 cores
```

### Performance Optimization Tips

```bash
# 1. Use MPI-IO for collective output (faster parallel writes):
#    In param.nml:
#    nc_out  = 1   ! NetCDF4 with HDF5 parallel I/O

# 2. Disable unused modules at compile time:
cmake ../src -DUSE_ICM=OFF -DUSE_SED=OFF -DUSE_WWM=OFF  # if not needed

# 3. Pin MPI tasks to cores (avoid NUMA effects):
mpirun -np 512 --map-by core --bind-to core pschism.exe

# 4. Use large page memory (Linux):
export MALLOC_MMAP_THRESHOLD_=134217728
export MALLOC_TRIM_THRESHOLD_=134217728

# 5. NetCDF output compression (reduce disk I/O):
#    In param.nml:
#    nc_compress = 1    ! Enable lossless compression
#    nc_dlevel = 4      ! Deflate level (1=fast, 9=max compression)

# 6. Reduce output frequency during spin-up:
#    nspool = 1440 during first 30 days (every 48 h)
#    nspool = 360 for analysis period (every 12 h)
```

### Verification Against BMKG Tide Gauges

```python
import pandas as pd
import numpy as np
from scipy import signal

def tidal_skill_score(obs, model, times_sec):
    """
    Compute standard tidal verification metrics.
    """
    # Tidal harmonic analysis (M2 + K1 for Indonesia)
    from uptide import Tides
    tp = Tides(['M2', 'S2', 'N2', 'K1', 'O1', 'P1'])
    tp.set_initial_time(times_sec[0])

    # RMSE
    rmse = np.sqrt(np.mean((obs - model)**2))

    # Bias
    bias = np.mean(model - obs)

    # Correlation
    r = np.corrcoef(obs, model)[0, 1]

    # Skill (Murphy 1988)
    ss = 1.0 - np.sum((obs - model)**2) / np.sum((obs - np.mean(obs))**2)

    # Willmott index of agreement
    d = 1.0 - (np.sum((obs - model)**2) /
                np.sum((np.abs(model - np.mean(obs)) +
                         np.abs(obs - np.mean(obs)))**2))

    return {'RMSE': rmse, 'bias': bias, 'R': r,
            'skill': ss, 'Willmott_d': d}

# BMKG tide gauges for validation:
gauges = {
    'Tanjung_Priok':   (106.883, -6.100),   # Jakarta
    'Semarang':        (110.433, -6.967),   # Central Java
    'Surabaya':        (112.733, -7.200),   # East Java
    'Benoa':           (115.217, -8.750),   # Bali
    'Lembar':          (116.083, -8.717),   # Lombok
    'Benete':          (116.833, -8.683),   # Sumbawa
    'Makassar':        (119.417, -5.133),   # South Sulawesi
    'Bitung':          (125.200,  1.433),   # North Sulawesi
}

# Target metrics (SCHISM Indonesian waters):
# RMSE < 0.15 m (straits), < 0.25 m (open coast)
# R > 0.95
# Willmott d > 0.97
```

---

*Guide prepared for SCHISM v5.12 — Indonesian Waters focus. For BMKG operational implementation, cross-reference with BMKG Marine HPC environment specifications and CMEMS data access credentials. Community support: https://github.com/schism-dev/schism/discussions*
