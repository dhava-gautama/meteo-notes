# MPAS-Atmosphere: Complete Guide from Mesh Setup to Visualization

> A practical guide covering the full MPAS (Model for Prediction Across Scales) atmosphere workflow.
> Based on **MPAS v8.3.1** (June 2025) — the latest stable release.

---

## Table of Contents

1. [Overview & Key Concepts](#1-overview--key-concepts)
2. [Meshes](#2-meshes)
3. [Downloading Initial & Boundary Condition Data](#3-downloading-initial--boundary-condition-data)
4. [Compilation & System Requirements](#4-compilation--system-requirements)
5. [Initialization (init_atmosphere)](#5-initialization-init_atmosphere)
6. [Running the Model](#6-running-the-model)
7. [Physics Options Reference](#7-physics-options-reference)
8. [Dynamics & Damping Configuration](#8-dynamics--damping-configuration)
9. [Configuring I/O (streams.atmosphere)](#9-configuring-io-streamsatmosphere)
10. [Regional (Limited-Area) Simulations](#10-regional-limited-area-simulations)
11. [Post-Processing & Visualization](#11-post-processing--visualization)
12. [Idealized Simulations](#12-idealized-simulations)
13. [Restart & Cycling](#13-restart--cycling)
14. [Troubleshooting & Tips](#14-troubleshooting--tips)
15. [Version History](#15-version-history)

---

## 1. Overview & Key Concepts

### What is MPAS?

MPAS (Model for Prediction Across Scales) is a global nonhydrostatic atmospheric model developed at NCAR. It uses an **unstructured centroidal Voronoi tessellation (CVT)** mesh instead of a regular lat-lon or rectangular grid. This allows **smooth variable resolution** — a single continuous mesh can transition from coarse (e.g., 60 km) to fine (e.g., 3 km) resolution without the abrupt boundaries of traditional grid nesting.

### MPAS vs WRF

| Feature | WRF | MPAS |
|---|---|---|
| **Grid type** | Structured rectangular (Arakawa-C) | Unstructured Voronoi (C-grid) |
| **Domain** | Regional (requires lateral BCs) | Global (no lateral BCs needed) or regional |
| **Nesting** | Abrupt step changes (3:1, 5:1) | Smooth mesh transitions |
| **Polar treatment** | Polar singularity (lat-lon) | No singularity (Voronoi cells) |
| **Resolution range** | Practical from ~1–100 km | Pre-built from 3–480 km |
| **Physics** | Extensive (50+ options per category) | WRF-derived subset (fewer but tested) |
| **Maturity** | 25+ years, massive user base | Newer (2013+), growing adoption |
| **Use case** | Regional NWP, research, operational | Global/variable-res NWP, climate downscaling |

### Key Concepts

**Voronoi mesh:** The domain is divided into polygons (mostly hexagons, with occasional pentagons/heptagons) where each cell center is equidistant from its cell boundary. This creates a nearly isotropic grid.

**SCVT (Spherical Centroidal Voronoi Tessellation):** The mesh generation process optimizes cell positions on the sphere so that each cell's generating point coincides with its centroid — maximizing mesh quality and isotropy.

**C-grid staggering:** Normal wind components are defined on cell edges (faces); scalar quantities (temperature, moisture, pressure) are at cell centers. This is the same staggering philosophy as WRF's Arakawa-C grid.

**Geometric-height vertical coordinate:** Unlike WRF's terrain-following sigma or hybrid sigma-pressure coordinate, MPAS uses a geometric-height-based coordinate (terrain-following at the surface, transitioning to flat levels aloft). Available as terrain-following or hybrid.

### Workflow Overview

```
[Download/Generate Mesh]
        |
        v
[init_atmosphere: Static fields]  -->  static.nc
        |
[WPS ungrib: GRIB -> Intermediate]
        |
[init_atmosphere: Vertical grid + Field interpolation]  -->  init.nc
        |
(Optional: Surface update)  -->  surface.nc
(Optional: LBCs for regional)  -->  lbc.*.nc
        |
        v
[atmosphere_model]  -->  history.*.nc, diag.*.nc
        |
        v
[convert_mpas / Python]  -->  lat-lon interpolation
        |
        v
[Visualization: matplotlib + cartopy / NCL / ParaView]
```

---

## 2. Meshes

### Pre-Built Meshes

NCAR provides pre-built SCVT meshes for download at [mpas-dev.github.io/atmosphere/atmosphere_meshes.html](https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html).

**Quasi-uniform (global, single resolution):**

| Resolution | Grid Cells | Mesh File Size | Static File Size |
|---|---|---|---|
| 480 km | 2,562 | 1.5 MB | 1.0 MB |
| 240 km | 10,242 | 6.3 MB | 4.0 MB |
| 120 km | 40,962 | 26 MB | 16 MB |
| 60 km | 163,842 | 106 MB | 70 MB |
| 30 km | 655,362 | 436 MB | 296 MB |
| 15 km | 2,621,442 | 1.7 GB | 1.4 GB |
| 10 km | 5,898,242 | 3.9 GB | 2.9 GB |
| 7.5 km | 10,485,762 | 6.9 GB | 6.5 GB |
| 4 km | 36,864,002 | 23.7 GB | 21.5 GB |
| 3 km | 65,536,002 | 42 GB | 38.6 GB |

**Variable-resolution (refinement centered at 0°N, 0°E by default):**

| Name | Grid Cells | Mesh Size | Description |
|---|---|---|---|
| 92–25 km | 163,842 | 135 MB | Coarse global with 25-km refined region |
| 60–15 km | 535,554 | 474 MB | Good for regional climate |
| 60–10 km | 999,426 | 913 MB | High-res regional within global |
| 60–3 km | 835,586 | 628 MB | Convection-permitting in small region |
| 15–3 km (circular) | 6,488,066 | 4.6 GB | Large convection-permitting region |
| 15–3 km (elliptical) | 8,060,930 | 6.0 GB | Extended convection-permitting region |

> **Large-mesh I/O:** Meshes finer than ~10 km require `io_type = "pnetcdf,cdf5"` or `"netcdf4"` in streams configuration to avoid the classic NetCDF 4 GB variable-size limit.

### Mesh Files

Each mesh download provides:

| File | Description |
|---|---|
| `x1.NCELLS.grid.nc` or `x5.NCELLS.grid.nc` | Mesh connectivity, coordinates, geometry |
| `x1.NCELLS.graph.info` | Graph connectivity for METIS partitioning |
| `x1.NCELLS.graph.info.part.N` | Pre-computed decompositions for N MPI tasks |

Naming convention: `x1` = quasi-uniform, `x5` = 5:1 variable resolution.

### Relocating Variable-Resolution Meshes (grid_rotate)

All variable-resolution meshes are centered at 0°N, 0°E by default. Use `grid_rotate` to move the refinement region to your area of interest.

```bash
# Clone MPAS-Tools
git clone https://github.com/MPAS-Dev/MPAS-Tools.git
cd MPAS-Tools/mesh_tools/grid_rotate

# Compile
make

# Configure namelist.input
cat > namelist.input << 'EOF'
&input
   config_original_latitude_degrees = 0.0
   config_original_longitude_degrees = 0.0
   config_new_latitude_degrees = 35.0       ! Move refinement to 35°N
   config_new_longitude_degrees = 135.0     ! 135°E (e.g., Japan)
   config_birdseye_rotation_counter_clockwise_degrees = 0.0
/
EOF

# Run (needs write access to grid file — make a copy first)
cp x5.535554.grid.nc x5.535554.grid.rotated.nc
./grid_rotate x5.535554.grid.rotated.nc
```

### Custom Mesh Generation (JIGSAW)

For meshes with custom refinement regions, use the **JIGSAW** mesh generator with `MPAS-Tools`:

```bash
# Install MPAS-Tools (includes JIGSAW bindings)
pip install mpas_tools

# Workflow:
# 1. Define a cell-width function (Python) specifying resolution across the globe
# 2. JIGSAW generates a triangulation based on the cell-width function
# 3. MPAS-Tools converts the triangulation to a Voronoi mesh (SCVT)
# 4. Optionally optimize with Lloyd iterations for centroidal property
```

**Example cell-width function (refine over Southeast Asia):**
```python
import numpy as np
from mpas_tools.mesh.creation import build_spherical_mesh

def cellWidthVsLatLon():
    lat = np.linspace(-90, 90, 721)
    lon = np.linspace(-180, 180, 1441)
    lon2d, lat2d = np.meshgrid(lon, lat)

    # Base resolution: 60 km
    cellWidth = 60.0 * np.ones_like(lat2d)

    # Refine to 15 km over Southeast Asia
    dist = np.sqrt((lat2d - 5.0)**2 + (lon2d - 110.0)**2)
    cellWidth = np.where(dist < 20.0, 15.0, cellWidth)

    # Smooth transition
    cellWidth = np.where((dist >= 20.0) & (dist < 35.0),
                          15.0 + (60.0 - 15.0) * (dist - 20.0) / 15.0,
                          cellWidth)

    return cellWidth, lon, lat

build_spherical_mesh(cellWidthVsLatLon, outFileName="custom_mesh.nc")
```

### Graph Partitioning (METIS)

To run MPAS in parallel, you need a mesh decomposition file matching your MPI task count.

```bash
# Install METIS (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
# Pre-computed partitions are included in mesh downloads for common task counts

# Generate partition for N MPI tasks
gpmetis -minconn -contig -niter=200 x1.163842.graph.info 64
# Produces: x1.163842.graph.info.part.64

# Copy to run directory
cp x1.163842.graph.info.part.64 graph.info.part.64
```

> **Naming convention in run directory:** MPAS expects `graph.info.part.N` (without the mesh prefix). Copy/rename accordingly.

---

## 3. Downloading Initial & Boundary Condition Data

MPAS uses the same meteorological input data as WRF — GRIB files processed through WPS `ungrib.exe` into the "intermediate" format.

### Supported Data Sources

| Dataset | Resolution | Frequency | Vtable | Notes |
|---|---|---|---|---|
| **NCEP GFS** | 0.25° | 3–6 hourly | `Vtable.GFS` | Most common for NWP |
| **NCEP GFS-FNL** | 0.25° | 6-hourly | `Vtable.GFS` | Analysis, more obs ingested |
| **ERA5** (pressure-level) | 0.25° | 1-hourly | `Vtable.ERA-interim.pl` | Best temporal resolution |
| **ERA5** (surface) | 0.25° | 1-hourly | `Vtable.ERA-interim.sfc` | Needed alongside pressure-level |
| **NCEP CFSv2** | 0.5° | 6-hourly | `Vtable.CFSR` | Reanalysis / climate |

### Ungrib Workflow

MPAS does **not** use `metgrid.exe` — it interpolates directly from WPS intermediate files using `init_atmosphere`. You only need `ungrib.exe` from WPS.

```bash
# In WPS directory
# 1. Link Vtable
ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable

# 2. Link GRIB files
./link_grib.csh /path/to/gfs/grib/files/*

# 3. Configure namelist.wps (only &share and &ungrib matter)
cat > namelist.wps << 'EOF'
&share
 wrf_core       = 'ARW',
 max_dom        = 1,
 start_date     = '2024-09-01_00:00:00',
 end_date       = '2024-09-06_00:00:00',
 interval_seconds = 21600,
/
&ungrib
 out_format = 'WPS',
 prefix     = 'FILE',
/
EOF

# 4. Run ungrib
./ungrib.exe

# Output: FILE:2024-09-01_00, FILE:2024-09-01_06, ..., FILE:2024-09-06_00
```

**For ERA5:** Run ungrib twice (pressure-level + surface), then concatenate:
```bash
# Pressure levels
ln -sf Vtable.ERA-interim.pl Vtable
# set prefix = 'PRES' in namelist.wps
./ungrib.exe

# Surface
ln -sf Vtable.ERA-interim.sfc Vtable
# set prefix = 'SFC' in namelist.wps
./ungrib.exe

# Concatenate for each time step
cat PRES:2024-09-01_00 SFC:2024-09-01_00 > FILE:2024-09-01_00
# Repeat for all time steps...
```

### Static Geographical Data

Download the same WPS static data used by WRF:
[www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html](https://www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html)

| Dataset | Size | Description |
|---|---|---|
| Default | ~11 GB | Standard resolutions, sufficient for most |
| Full resolution | ~29 GB | Required for high-resolution meshes (< 10 km) |

Extract to a directory (e.g., `/data/WPS_GEOG/`) and set `config_geog_data_path` in the init namelist.

**MPAS-specific static data** (alternative to full WPS GEOG):
- Minimal static dataset: [www2.mmm.ucar.edu/projects/mpas/mpas_static.tar.bz2](https://www2.mmm.ucar.edu/projects/mpas/mpas_static.tar.bz2)

**Aerosol climatology** (needed for Thompson aerosol-aware microphysics):
- Download: [www2.mmm.ucar.edu/projects/mpas/QNWFA_QNIFA_SIGMA_MONTHLY.dat](https://www2.mmm.ucar.edu/projects/mpas/QNWFA_QNIFA_SIGMA_MONTHLY.dat)
- Place in run directory when using `mp_thompson_aerosols`

---

## 4. Compilation & System Requirements

### Required Libraries

| Library | Purpose |
|---|---|
| **Fortran compiler** (gfortran 9+, ifort/ifx, nvfortran) | Model compilation |
| **C compiler** (gcc, icc) | Support code |
| **MPI** (OpenMPI, MPICH, Intel MPI) | Parallel execution (MPAS has no serial mode) |
| **NetCDF-C** (4.6+) | I/O (mandatory) |
| **NetCDF-Fortran** (4.4+) | Fortran bindings (mandatory) |
| **HDF5** (1.8+) | NetCDF-4 compression |
| **Parallel-NetCDF** (pnetcdf, 1.11+) | Parallel I/O (recommended, especially for large meshes) |
| **PIO** (Parallel I/O, 1.9+) | MPAS I/O abstraction layer |
| **METIS** (5.1+) | Mesh partitioning |

> **Key difference from WRF:** MPAS requires PIO (Parallel I/O library) and cannot run in serial mode — MPI is always required.

### Environment Variables

```bash
export NETCDF=/usr/local/netcdf
export PNETCDF=/usr/local/pnetcdf
export PIO=/usr/local/pio
export METIS_ROOT=/usr/local/metis
export PATH=$NETCDF/bin:$PATH
export LD_LIBRARY_PATH=$NETCDF/lib:$PNETCDF/lib:$PIO/lib:$LD_LIBRARY_PATH
```

### Compilation

```bash
# Clone MPAS
git clone https://github.com/MPAS-Dev/MPAS-Model.git
cd MPAS-Model

# Compile init_atmosphere core
make -j4 gfortran CORE=init_atmosphere PRECISION=single USE_PIO2=true

# Clean between cores (required)
make clean CORE=atmosphere

# Compile atmosphere core
make -j4 gfortran CORE=atmosphere PRECISION=single USE_PIO2=true
```

**Compiler targets:**

| Target | Compiler |
|---|---|
| `gfortran` | GNU gfortran + gcc |
| `ifort` | Intel Classic (ifort + icc) |
| `ifx` | Intel oneAPI (ifx + icx) |
| `nvfortran` | NVIDIA HPC SDK |
| `llvm` | LLVM flang + clang |

**Build options:**

| Flag | Values | Description |
|---|---|---|
| `CORE=` | `init_atmosphere`, `atmosphere` | Which core to build |
| `PRECISION=` | `single`, `double` | Floating-point precision. Single is recommended for efficiency (~30% faster). |
| `USE_PIO2=` | `true`, `false` | Use PIO 2.x (recommended if using PIO) |
| `AUTOCLEAN=` | `true` | Auto-clean when switching between cores (avoids manual `make clean`) |
| `OPENMP=` | `true` | Enable hybrid MPI+OpenMP (experimental) |
| `DEBUG=` | `true` | Enable debug flags and bounds checking |

**Verify compilation:**
```bash
ls -l init_atmosphere_model    # From init_atmosphere build
ls -l atmosphere_model         # From atmosphere build
```

> **v8.0+ SMIOL:** MPAS 8.0 introduced the **Simple MPAS I/O Layer (SMIOL)**, a built-in alternative to PIO. If the `$PIO` environment variable is **not set**, MPAS defaults to SMIOL, which requires only Parallel-netCDF. This significantly reduces external dependency complexity — PIO compilation is no longer needed.

**Simplified build with SMIOL (no PIO needed):**
```bash
# Only NETCDF and PNETCDF needed
unset PIO
make -j4 gfortran CORE=atmosphere PRECISION=single
```

### Hardware Requirements

| Configuration | MPI Ranks | RAM | Mesh | Use Case |
|---|---|---|---|---|
| Learning | 4–8 | 8–16 GB | 240 km (10K cells) | Tutorial, testing |
| Desktop Research | 16–64 | 32–128 GB | 60 km (164K cells) | Global weather, ~5 day forecasts |
| Serious Research | 64–512 | 128–512 GB | 15 km (2.6M cells) | High-res global or variable-res |
| HPC / Operational | 512–10,000+ | 1+ TB | 3–4 km (37–66M cells) | Convection-permitting global |

---

## 5. Initialization (init_atmosphere)

The `init_atmosphere_model` executable creates initial condition files in two steps: static field interpolation, then meteorological field interpolation.

### Step 1: Static Field Interpolation

Interpolates time-invariant geographical data (terrain, land use, soil type, vegetation) onto the MPAS mesh.

**Required files in run directory:**
- `init_atmosphere_model` (executable)
- Mesh file (e.g., `x1.163842.grid.nc`)
- `namelist.init_atmosphere`
- `streams.init_atmosphere`
- WPS geographic data (pointed to by namelist)

**namelist.init_atmosphere:**
```
&nhyd_model
    config_init_case = 7
    config_start_time = '2024-09-01_00:00:00'
/
&dimensions
    config_nvertlevels = 1
    config_nsoillevels = 1
    config_nfglevels = 1
    config_nfgsoillevels = 1
/
&data_sources
    config_geog_data_path = '/data/WPS_GEOG/'
    config_landuse_data = 'MODIFIED_IGBP_MODIS_NOAH'
    config_topo_data = 'GMTED2010'
    config_vegfrac_data = 'MODIS'
    config_albedo_data = 'MODIS'
    config_maxsnowalbedo_data = 'MODIS'
    config_supersample_factor = 3
/
&preproc_stages
    config_static_interp = true
    config_native_gwd_static = true
    config_vertical_grid = false
    config_met_interp = false
    config_input_sst = false
    config_frac_seaice = false
/
```

**streams.init_atmosphere (for static):**
```xml
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="x1.163842.grid.nc"
                  input_interval="initial_only" />

<immutable_stream name="output"
                  type="output"
                  filename_template="x1.163842.static.nc"
                  output_interval="initial_only" />
</streams>
```

**Run:**
```bash
# Static interpolation MUST run on a single MPI task
mpiexec -n 1 ./init_atmosphere_model
```

> **Critical:** The static interpolation step **must run serially** (1 MPI task). This can take up to an hour for large meshes.

### Step 2: Meteorological Field Interpolation

Interpolates atmospheric fields from WPS intermediate files onto the MPAS mesh, creates the vertical grid, and generates `init.nc`.

**Required files:**
- `init_atmosphere_model`
- `x1.163842.static.nc` (from Step 1)
- WPS intermediate files (`FILE:YYYY-MM-DD_HH`)
- `graph.info.part.N` (if running parallel)

**namelist.init_atmosphere:**
```
&nhyd_model
    config_init_case = 7
    config_start_time = '2024-09-01_00:00:00'
    config_theta_adv_order = 3
    config_coef_3rd_order = 0.25
/
&dimensions
    config_nvertlevels = 55
    config_nsoillevels = 4
    config_nfglevels = 38
    config_nfgsoillevels = 4
/
&data_sources
    config_geog_data_path = '/data/WPS_GEOG/'
    config_met_prefix = 'FILE'
    config_use_spechumd = false
/
&vertical_grid
    config_ztop = 30000.0
    config_nsmterrain = 1
    config_smooth_surfaces = true
    config_blend_bdy_terrain = false
/
&preproc_stages
    config_static_interp = false
    config_native_gwd_static = false
    config_vertical_grid = true
    config_met_interp = true
    config_input_sst = false
    config_frac_seaice = true
/
```

**Key dimension parameters:**

| Parameter | Description | How to Determine |
|---|---|---|
| `config_nvertlevels` | Number of model vertical levels | User choice: 55 (default since v6.0), 41 (older default), or custom |
| `config_nsoillevels` | Soil levels in the model | Typically 4 (Noah LSM) |
| `config_nfglevels` | Number of levels in intermediate file | Check with `rd_intermediate.exe` or count levels in GRIB |
| `config_nfgsoillevels` | Soil levels in intermediate file | GFS=4, ERA5=4, CFSR=4 |

**Finding `config_nfglevels`:**
```bash
# Use ncdump on an intermediate file, or count pressure levels in the GRIB data
# GFS 0.25° (2024): typically 34 pressure levels → config_nfglevels = 34
# ERA5 (137 model levels): config_nfglevels = 137 (if using model levels)
# ERA5 (37 pressure levels): config_nfglevels = 38 (37 pressure + surface)
```

**streams.init_atmosphere (for met interp):**
```xml
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="x1.163842.static.nc"
                  input_interval="initial_only" />

<immutable_stream name="output"
                  type="output"
                  filename_template="x1.163842.init.nc"
                  output_interval="initial_only" />
</streams>
```

**Run:**
```bash
# Can run in parallel for this step
mpiexec -n 8 ./init_atmosphere_model
```

### Step 3 (Optional): Surface Update File (Case 8)

For multi-day simulations, create time-varying SST/sea-ice files:

**namelist.init_atmosphere (surface update):**
```
&nhyd_model
    config_init_case = 8
    config_start_time = '2024-09-01_00:00:00'
    config_stop_time = '2024-09-06_00:00:00'
/
&dimensions
    config_nfglevels = 38
    config_nfgsoillevels = 4
/
&data_sources
    config_met_prefix = 'FILE'
    config_sfc_prefix = 'FILE'
    config_fg_interval = 86400           ! SST update interval (seconds); 86400 = daily
/
&preproc_stages
    config_static_interp = false
    config_native_gwd_static = false
    config_vertical_grid = false
    config_met_interp = false
    config_input_sst = true
    config_frac_seaice = true
/
```

**streams.init_atmosphere (surface):**
```xml
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="x1.163842.static.nc"
                  input_interval="initial_only" />

<immutable_stream name="surface"
                  type="output"
                  filename_template="x1.163842.surface.nc"
                  output_interval="86400"
                  packages="initial_conds" />
</streams>
```

---

## 6. Running the Model

### Required Files in Run Directory

| File | Description |
|---|---|
| `atmosphere_model` | Compiled executable |
| `init.nc` (or `x1.163842.init.nc`) | Initial conditions from init_atmosphere |
| `graph.info.part.N` | METIS decomposition for N MPI tasks |
| `namelist.atmosphere` | Model configuration |
| `streams.atmosphere` | I/O configuration |
| `stream_list.atmosphere.*` | Variable lists for output streams |
| Physics lookup tables | `*.TBL`, `*.DBL`, `RRTMG_*` files (in source directory) |
| `surface.nc` (optional) | Time-varying SST/sea-ice |
| `lbc.*.nc` (optional) | Lateral boundary conditions (regional only) |

> **Lookup tables:** Copy physics lookup tables from the MPAS source directory: `cp MPAS-Model/src/core_atmosphere/physics/physics_wrf/files/* ./`

### namelist.atmosphere

```
&nhyd_model
    config_time_integration_order = 2
    config_dt = 720.0                     ! Timestep (seconds)
    config_start_time = '2024-09-01_00:00:00'
    config_run_duration = '5_00:00:00'    ! 5 days
    config_split_dynamics_transport = true
    config_number_of_sub_steps = 2
    config_dynamics_split_steps = 3
    config_h_mom_eddy_visc2 = 0.0
    config_h_mom_eddy_visc4 = 0.0
    config_h_theta_eddy_visc2 = 0.0
    config_h_theta_eddy_visc4 = 0.0
    config_len_disp = 120000.0            ! Smallest cell distance (m)
    config_visc4_2dsmag = 0.05
    config_theta_adv_order = 3
    config_scalar_adv_order = 3
    config_coef_3rd_order = 0.25
    config_positive_definite = false
    config_monotonic = true
/
&damping
    config_zd = 22000.0                   ! Damping layer depth from model top (m)
    config_xnutr = 0.2                    ! Rayleigh damping coefficient
/
&limited_area
    config_apply_lbcs = false             ! Set true for regional runs
/
&physics
    config_physics_suite = 'mesoscale_reference'
    config_sst_update = true
    config_microp_scheme = 'mp_wsm6'
    config_convection_scheme = 'cu_ntiedtke'
    config_pbl_scheme = 'bl_ysu'
    config_gwdo_scheme = 'bl_ysu_gwdo'
    config_radt_lw_scheme = 'rrtmg_lw'
    config_radt_sw_scheme = 'rrtmg_sw'
    config_radt_cld_scheme = 'cld_fraction'
    config_sfclayer_scheme = 'sf_monin_obukhov_rev'
    config_lsm_scheme = 'sf_noah'
    config_convection_scheme = 'cu_ntiedtke'
    config_radt_lw_interval = '01:00:00'
    config_radt_sw_interval = '01:00:00'
/
&soundings
    config_sounding_interval = 'none'
/
&decomposition
    config_block_decomp_file_prefix = 'graph.info.part.'
/
&restart
    config_do_restart = false
/
&iau
    config_IAU_option = 'off'
/
```

### Timestep and Resolution

The same rule of thumb as WRF: **`config_dt ≈ 6 × dx_km`**

| Minimum Cell Spacing | `config_dt` | `config_len_disp` |
|---|---|---|
| 240 km | 1200–1440 s | 240000 |
| 120 km | 720 s | 120000 |
| 60 km | 360 s | 60000 |
| 30 km | 180 s | 30000 |
| 15 km | 90 s | 15000 |
| 10 km | 60 s | 10000 |
| 3 km | 18 s | 3000 |

> **Variable-resolution meshes:** Use the **finest** (smallest) cell spacing to determine `config_dt` and `config_len_disp`. A 60–15 km mesh uses the 15-km values.

### Execution

```bash
# Run with N MPI tasks (must match graph.info.part.N)
mpiexec -n 64 ./atmosphere_model

# Monitor progress
tail -f log.atmosphere.0000.out

# Successful completion shows:
# "Finished running the atmosphere core"
```

---

## 7. Physics Options Reference

### Physics Suites

MPAS provides pre-tested combinations of physics schemes called **suites**. Setting `config_physics_suite` auto-configures all physics options.

| Suite | Description |
|---|---|
| `'mesoscale_reference'` | Standard suite for dx > ~10 km. WSM6 microphysics, New Tiedtke convection, YSU PBL, RRTMG radiation, Noah LSM. |
| `'convection_permitting'` | For dx < ~10 km. Thompson microphysics, Grell-Freitas (scale-aware) convection, MYNN PBL, RRTMG radiation, Noah LSM. |
| `'none'` | No physics — all schemes set to `'off'`. For idealized dry dynamics tests. |

> **Overriding individual schemes:** After setting `config_physics_suite`, you can override any individual scheme. The suite sets defaults; individual settings take precedence.

### Microphysics (`config_microp_scheme`)

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No microphysics. |
| `'mp_kessler'` | Kessler | Simple warm-rain only (cloud, rain). No ice phase. For idealized cases. |
| `'mp_wsm6'` | WSM6 | Single-moment, 6-class (vapor, cloud, rain, ice, snow, graupel). Good balance of speed and accuracy. Default in `mesoscale_reference` suite. From WRF. |
| `'mp_thompson'` | Thompson | Hybrid single/double-moment — double-moment for rain and ice. Excellent for mixed-phase clouds, winter precipitation. Default in `convection_permitting` suite. From WRF. |
| `'mp_thompson_aerosols'` | Thompson Aerosol-Aware | Thompson with prognostic water- and ice-friendly aerosol. Accounts for aerosol-cloud interactions affecting droplet number and precipitation efficiency. |

### Cumulus Convection (`config_convection_scheme`)

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No cumulus parameterization. Use for convection-permitting (dx < ~4 km). |
| `'cu_kain_fritsch'` | Kain-Fritsch | Mass-flux with CAPE-removal closure. Deep + shallow convection. Most widely used. Best for dx > 10 km. From WRF. |
| `'cu_tiedtke'` | Tiedtke | Mass-flux with moisture-convergence closure. Handles deep, shallow, and mid-level convection. Originally from ECMWF. |
| `'cu_ntiedtke'` | New Tiedtke | Updated Tiedtke with CAPE-based closure for deep convection. Improved diurnal cycle representation. Default in `mesoscale_reference`. From WRF. |
| `'cu_grell_freitas'` | Grell-Freitas | Scale-aware ensemble mass-flux. Automatically reduces contribution as resolution approaches convection-permitting. Ideal for variable-resolution meshes. Default in `convection_permitting`. |

> **Variable-resolution advantage:** Grell-Freitas is particularly well-suited for MPAS because its scale-awareness handles the smooth resolution transition within a single mesh — no need to toggle convection on/off between domains.

### Planetary Boundary Layer (`config_pbl_scheme`)

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No PBL parameterization. |
| `'bl_ysu'` | YSU | Yonsei University non-local scheme. Counter-gradient flux for well-mixed daytime boundary layers. Most popular general-purpose PBL scheme. Default in `mesoscale_reference`. Must pair with `sf_monin_obukhov` or `sf_monin_obukhov_rev`. |
| `'bl_mynn'` | MYNN | Mellor-Yamada-Nakanishi-Niino Level 2.5. Prognostic TKE, local closure with EDMF extension for shallow convection. Superior for marine/coastal, fog, and low cloud. Default in `convection_permitting`. Must pair with `sf_mynn`. |

### Surface Layer (`config_sfclayer_scheme`)

| Value | Scheme | Compatible PBL | Description |
|---|---|---|---|
| `'off'` | None | — | No surface layer. |
| `'sf_monin_obukhov'` | Original M-O | YSU | Monin-Obukhov similarity with standard stability functions. |
| `'sf_monin_obukhov_rev'` | Revised MM5 | YSU | Revised Monin-Obukhov with Carlson-Boland viscous sub-layer. Improved over original. Default in `mesoscale_reference`. |
| `'sf_mynn'` | MYNN | MYNN | Nakanishi-Niino surface layer. Required for MYNN PBL. Default in `convection_permitting`. |

### Radiation

**Longwave (`config_radt_lw_scheme`):**

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No longwave radiation. |
| `'rrtmg_lw'` | RRTMG | Rapid Radiative Transfer Model for GCMs. 16 spectral bands, McICA for sub-grid cloud variability. Default and recommended. |
| `'cam_lw'` | CAM | Community Atmosphere Model radiation. Used in CESM. |

**Shortwave (`config_radt_sw_scheme`):**

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No shortwave radiation. |
| `'rrtmg_sw'` | RRTMG | 14 spectral bands, McICA. Default and recommended. |
| `'cam_sw'` | CAM | Community Atmosphere Model shortwave. |

**Radiation calling interval:**
```
config_radtlw_interval = '01:00:00'     ! Call LW radiation every 1 hour
config_radtsw_interval = '01:00:00'     ! Call SW radiation every 1 hour
```

For convection-permitting (dx < 10 km), reduce to `'00:30:00'` or `'00:15:00'`.

### Physics Timing & Control Parameters

| Parameter | Default | Description |
|---|---|---|
| `config_radtlw_interval` | `'01:00:00'` | Longwave radiation call interval |
| `config_radtsw_interval` | `'01:00:00'` | Shortwave radiation call interval |
| `config_conv_interval` | `'none'` | Convection call interval (`'none'` = every timestep) |
| `config_pbl_interval` | `'none'` | PBL call interval |
| `config_n_microp` | 1 | Microphysics sub-steps per timestep |
| `config_sst_update` | false | Enable SST/sea-ice updates from `surface.nc` |
| `config_frac_seaice` | false | Fractional sea-ice treatment |
| `config_bucket_update` | `'none'` | Precipitation bucket reset interval |
| `config_bucket_rainc` | 100.0 | Convective rain bucket size (mm) |
| `config_bucket_rainnc` | 100.0 | Grid-scale rain bucket size (mm) |

### Cloud Fraction (`config_radt_cld_scheme`)

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No cloud fraction. |
| `'cld_incidence'` | Incidence | Binary: cloud present or not based on condensate threshold. Simplest. |
| `'cld_fraction'` | Xu-Randall | Diagnostic cloud fraction based on relative humidity and condensate. Default. |
| `'cld_fraction_thompson'` | Thompson-specific | Cloud fraction tuned for Thompson microphysics. Use with `mp_thompson`. |

### Land Surface Model (`config_lsm_scheme`)

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No land surface model. |
| `'sf_noah'` | Noah | 4-layer soil model. Predicts soil temperature, moisture, snowpack, canopy water. NCEP operational standard. Default. Updated to WRF 4.5 in MPAS v8.0. |
| `'sf_noahmp'` | Noah-MP | Noah Multi-Physics. Separate vegetation canopy, dynamic vegetation, multi-layer snowpack, groundwater. Most advanced option. |

### Gravity Wave Drag (`config_gwdo_scheme`)

| Value | Scheme | Description |
|---|---|---|
| `'off'` | None | No gravity wave drag. Acceptable for dx < ~10 km. |
| `'bl_ysu_gwdo'` | YSU GWDO | Orographic gravity wave drag + low-level blocking. Default for both suites. Computed on native MPAS mesh (not sub-grid). |
| `'bl_ugwp_gwdo'` | Unified GWP | Unified Gravity Wave Physics. Newer, includes both orographic and non-orographic gravity waves. |

### Recommended Configurations

**Global NWP (60 km quasi-uniform):**
```
config_physics_suite = 'mesoscale_reference'
! Uses: WSM6, New Tiedtke, YSU, RRTMG, Noah, Revised MM5 surface layer
config_dt = 360.0
config_len_disp = 60000.0
```

**Variable-resolution regional focus (60–15 km):**
```
config_physics_suite = 'convection_permitting'
! Uses: Thompson, Grell-Freitas (scale-aware), MYNN, RRTMG, Noah
config_dt = 90.0
config_len_disp = 15000.0
```

**Convection-permitting global (3–4 km):**
```
config_physics_suite = 'convection_permitting'
config_convection_scheme = 'off'     ! Override: fully explicit convection
config_dt = 18.0
config_len_disp = 3000.0
```

---

## 8. Dynamics & Damping Configuration

### Dynamical Core

MPAS solves the fully compressible nonhydrostatic equations using:
- **Split-explicit time integration:** 3rd-order Runge-Kutta for meteorological modes, forward-backward for acoustic modes
- **Vector-invariant horizontal momentum:** Avoids metric terms on the unstructured mesh
- **Scalar transport:** Flux-form, monotonic option with 3rd-order accurate reconstruction

### Key Dynamics Parameters (`&nhyd_model`)

| Parameter | Default | Description |
|---|---|---|
| `config_time_integration_order` | 2 | Runge-Kutta order (2 or 3). 2 is standard. |
| `config_split_dynamics_transport` | true | Split dynamics and transport integration. Recommended. |
| `config_number_of_sub_steps` | 2 | Acoustic sub-steps per RK step. |
| `config_dynamics_split_steps` | 3 | Dynamics steps per transport step. |
| `config_theta_adv_order` | 3 | Potential temperature advection order (2, 3, or 4). |
| `config_scalar_adv_order` | 3 | Scalar advection order (2, 3, or 4). |
| `config_coef_3rd_order` | 0.25 | Blending coefficient for 3rd-order advection (0.0–1.0). Higher = more upwinding. |
| `config_monotonic` | true | Monotonic scalar transport (prevents negative moisture). |
| `config_positive_definite` | false | Positive-definite transport (alternative to monotonic). |

### Diffusion

MPAS uses Smagorinsky-type 4th-order hyperdiffusion scaled by the local mesh spacing.

| Parameter | Default | Description |
|---|---|---|
| `config_visc4_2dsmag` | 0.05 | 2D Smagorinsky coefficient for 4th-order horizontal diffusion. Controls numerical diffusion based on local deformation rate. |
| `config_h_mom_eddy_visc2` | 0.0 | 2nd-order explicit momentum diffusion (m²/s). Usually 0 — let Smagorinsky handle it. |
| `config_h_mom_eddy_visc4` | 0.0 | 4th-order explicit momentum diffusion (m⁴/s). Usually 0. |
| `config_h_theta_eddy_visc2` | 0.0 | 2nd-order theta diffusion. |
| `config_h_theta_eddy_visc4` | 0.0 | 4th-order theta diffusion. |
| `config_len_disp` | — | Smallest cell-to-cell distance (m). Scales the Smagorinsky diffusion. Must be set correctly for your mesh. |

> **Variable-resolution note:** `config_len_disp` should be set to the **smallest** cell spacing in your mesh. The Smagorinsky diffusion automatically scales with local cell size, so the model adapts diffusion to the local resolution.

### Upper Damping (`&damping`)

Prevents spurious wave reflection from the model top.

| Parameter | Default | Description |
|---|---|---|
| `config_zd` | 22000.0 | Depth of damping layer below model top (m). |
| `config_xnutr` | 0.2 | Rayleigh damping coefficient. Higher = stronger damping. |

### Vertical Grid

| Parameter | Default | Description |
|---|---|---|
| `config_ztop` | 30000.0 | Model top height (m). 30 km standard. Use 40–60 km for stratospheric studies. |
| `config_nvertlevels` | 55 | Number of vertical levels. |
| `config_nsmterrain` | 1 | Number of terrain smoothing passes. |
| `config_smooth_surfaces` | true | Apply Laplacian smoothing to terrain. |

---

## 9. Configuring I/O (streams.atmosphere)

MPAS uses XML-based stream configuration instead of Fortran namelists for I/O. This is more flexible than WRF's approach.

### Stream Types

| Type | Description |
|---|---|
| `"input"` | Read-only |
| `"output"` | Write-only |
| `"input;output"` | Read and write (restart files) |
| `"none"` | Field grouping only (not read/written) |

### Default Streams

```xml
<streams>

<!-- Initial conditions input -->
<immutable_stream name="input"
                  type="input"
                  filename_template="init.nc"
                  input_interval="initial_only" />

<!-- Restart files -->
<immutable_stream name="restart"
                  type="input;output"
                  filename_template="restart.$Y-$M-$D_$h.$m.$s.nc"
                  input_interval="initial_only"
                  output_interval="1_00:00:00" />

<!-- Model history output -->
<stream name="output"
        type="output"
        filename_template="history.$Y-$M-$D_$h.$m.$s.nc"
        output_interval="6:00:00" >
    <file name="stream_list.atmosphere.output"/>
</stream>

<!-- 2D diagnostics -->
<stream name="diagnostics"
        type="output"
        filename_template="diag.$Y-$M-$D_$h.$m.$s.nc"
        output_interval="3:00:00" >
    <file name="stream_list.atmosphere.diagnostics"/>
</stream>

<!-- SST / sea-ice updates -->
<stream name="surface"
        type="input"
        filename_template="surface.nc"
        input_interval="24:00:00" >
    <file name="stream_list.atmosphere.surface"/>
</stream>

</streams>
```

### Filename Template Variables

| Variable | Meaning |
|---|---|
| `$Y` | Year (4-digit) |
| `$M` | Month (2-digit) |
| `$D` | Day (2-digit) |
| `$h` | Hour (2-digit) |
| `$m` | Minute (2-digit) |
| `$s` | Second (2-digit) |

### Adding/Removing Variables

**Method 1: Using `<var>` elements:**
```xml
<stream name="custom_output"
        type="output"
        filename_template="custom.$Y-$M-$D_$h.$m.$s.nc"
        output_interval="1:00:00">
    <var name="temperature"/>
    <var name="uReconstructZonal"/>
    <var name="uReconstructMeridional"/>
    <var name="surface_pressure"/>
    <var name="precipw"/>
</stream>
```

**Method 2: Using external file lists:**
```xml
<stream name="output"
        type="output"
        filename_template="history.$Y-$M-$D_$h.$m.$s.nc"
        output_interval="6:00:00">
    <file name="stream_list.atmosphere.output"/>
</stream>
```

Where `stream_list.atmosphere.output` is a text file with one variable name per line:
```
temperature
rho
uReconstructZonal
uReconstructMeridional
pressure_p
theta
```

**Method 3: Using `<var_struct>` to include all variables in a structure:**
```xml
<var_struct name="diag"/>      <!-- All diagnostic variables -->
<var_struct name="state"/>     <!-- All state variables -->
```

### Advanced Options

| Attribute | Description | Example |
|---|---|---|
| `filename_interval` | How often to create a new file | `"01-00_00:00:00"` (monthly) |
| `reference_time` | Anchor for filename intervals | `"2024-01-01_00:00:00"` |
| `clobber_mode` | Overwrite behavior | `"overwrite"`, `"append"`, `"never_modify"` |
| `precision` | Float precision in output | `"single"`, `"double"`, `"native"` |
| `io_type` | I/O backend | `"pnetcdf"`, `"pnetcdf,cdf5"`, `"netcdf4"` |

> **Reducing output size:** Use `precision="single"` for history output — halves file size with negligible impact on post-processing. Use `"pnetcdf,cdf5"` or `"netcdf4"` for meshes with > 4M cells.

---

## 10. Regional (Limited-Area) Simulations

Since v7.0, MPAS supports regional (limited-area) simulations on a subset of a global mesh.

### Creating a Regional Mesh

```bash
# Clone the Limited-Area tool
git clone https://github.com/MiCurry/MPAS-Limited-Area.git
cd MPAS-Limited-Area
```

**Define a region (points file):**
```
Name: Mediterranean
Type: ellipse
Point: 37.9, 18.0
Semi-major-axis: 3200000.
Semi-minor-axis: 1700000.
Orientation-angle: 100.
```

**Region types available:**

| Type | Parameters | Description |
|---|---|---|
| `circle` | Center point + radius | Circular region |
| `ellipse` | Center + semi-axes + orientation | Elliptical region |
| `channel` | Upper/lower latitude bounds | Equatorial/latitudinal band |
| `custom` | List of lat/lon vertices | Arbitrary convex polygon |

**Extract regional mesh:**
```bash
python create_region.py \
    --grid x1.163842.grid.nc \
    --graph x1.163842.graph.info \
    --region Mediterranean.region \
    --output Mediterranean.grid.nc
```

### Regional Initialization

Same as global but with key differences:

**In `namelist.init_atmosphere` (Step 2 — met interp):**
```
&vertical_grid
    config_blend_bdy_terrain = true      ! MUST be true for regional
/
```

**For lateral boundary conditions (init_atmosphere Case 9):**
```
&nhyd_model
    config_init_case = 9
    config_start_time = '2024-09-01_00:00:00'
    config_stop_time = '2024-09-04_00:00:00'
/
&data_sources
    config_met_prefix = 'FILE'
    config_fg_interval = 21600           ! LBC update interval (seconds)
/
```

This produces `lbc.2024-09-01_00.00.00.nc`, `lbc.2024-09-01_06.00.00.nc`, etc.

### Running Regional MPAS

**In `namelist.atmosphere`:**
```
&limited_area
    config_apply_lbcs = true             ! Enable lateral boundary conditions
/
```

**In `streams.atmosphere`:**
```xml
<stream name="lbc_in"
        type="input"
        filename_template="lbc.$Y-$M-$D_$h.$m.$s.nc"
        input_interval="6:00:00" />
```

> **LBC update frequency:** Must match `config_fg_interval` used when generating LBC files. 6-hourly (GFS) or 1-hourly (ERA5) depending on your input data.

---

## 11. Post-Processing & Visualization

### The Unstructured Grid Challenge

MPAS output is on an unstructured Voronoi mesh — most plotting tools expect regular lat-lon grids. Three approaches:

### Approach 1: convert_mpas (Interpolate to Lat-Lon)

The `convert_mpas` utility interpolates MPAS output to a regular lat-lon grid, producing standard NetCDF files compatible with any plotting tool.

```bash
# Clone and compile
git clone https://github.com/mgduda/convert_mpas.git
cd convert_mpas
make

# Basic usage: interpolate to 0.5° grid
./convert_mpas history.2024-09-01_00.00.00.nc

# With target domain (subset region):
cat > target_domain << 'EOF'
startlat=-10.0
startlon=90.0
endlat=30.0
endlon=160.0
EOF
./convert_mpas history.2024-09-01_00.00.00.nc

# With field subsetting:
cat > include_fields << 'EOF'
temperature
uReconstructZonal
uReconstructMeridional
surface_pressure
precipw
EOF
./convert_mpas history.2024-09-01_00.00.00.nc
```

Output: `latlon.nc` — a regular-grid NetCDF file viewable with `ncview`, Python, GrADS, etc.

### Approach 2: Direct Python Plotting (Native Mesh)

Plot directly on the unstructured mesh using `matplotlib` triangulation or cell patches.

```python
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import cartopy.crs as crs
import cartopy.feature as cfeature

# Read MPAS output
nc = Dataset("history.2024-09-01_06.00.00.nc")
latCell = np.degrees(nc.variables["latCell"][0] if nc.variables["latCell"].ndim > 1
                     else nc.variables["latCell"][:])
lonCell = np.degrees(nc.variables["lonCell"][0] if nc.variables["lonCell"].ndim > 1
                     else nc.variables["lonCell"][:])
t2m = nc.variables["t2m"][0, :]  # 2-m temperature

# Create triangulation from cell centers
triang = tri.Triangulation(lonCell, latCell)

# Plot
fig, ax = plt.subplots(subplot_kw={"projection": crs.PlateCarree()}, figsize=(12, 6))
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
cf = ax.tricontourf(triang, t2m - 273.15, levels=20, cmap="RdYlBu_r",
                     transform=crs.PlateCarree())
plt.colorbar(cf, ax=ax, shrink=0.7, label="2-m Temperature (°C)")
plt.title("MPAS 2-m Temperature")
plt.savefig("mpas_t2m.png", dpi=150, bbox_inches="tight")
```

### Approach 3: ParaView (3D Interactive)

ParaView has a built-in **NetCDFMPASreader** that directly reads MPAS unstructured output. For more control, use the MPAS-Tools VTK extractor:

```bash
# Using mpas_tools Python package
python -m mpas_tools.viz.paraview_extractor \
    --filename_pattern "history.*.nc" \
    --mesh_filename init.nc \
    --variable_list temperature,uReconstructZonal \
    --dimension_list nCells,nVertLevels
```

### Approach 4: VAPOR (3D Interactive)

NCAR's VAPOR supports direct import of MPAS unstructured output for interactive 3D visualization.

### Key MPAS Output Variables

| Variable | Description | Units |
|---|---|---|
| `temperature` | Temperature on model levels | K |
| `theta` | Potential temperature | K |
| `rho` | Dry air density | kg/m³ |
| `uReconstructZonal` | Earth-relative zonal wind (reconstructed at cell centers) | m/s |
| `uReconstructMeridional` | Earth-relative meridional wind (reconstructed at cell centers) | m/s |
| `w` | Vertical velocity | m/s |
| `pressure_p` | Perturbation pressure | Pa |
| `surface_pressure` | Surface pressure | Pa |
| `t2m` | 2-m temperature | K |
| `q2` | 2-m specific humidity | kg/kg |
| `u10` | 10-m zonal wind | m/s |
| `v10` | 10-m meridional wind | m/s |
| `precipw` | Precipitable water | kg/m² |
| `rainnc` | Accumulated grid-scale precipitation | mm |
| `rainc` | Accumulated convective precipitation | mm |
| `mslp` | Mean sea level pressure | Pa |
| `olrtoa` | Outgoing longwave at TOA | W/m² |
| `qv` | Water vapor mixing ratio (3D) | kg/kg |
| `qc` | Cloud water mixing ratio (3D) | kg/kg |

> **Earth-relative winds:** Unlike WRF, MPAS provides pre-rotated earth-relative wind components (`uReconstructZonal`, `uReconstructMeridional`) directly in the output. No manual rotation needed.

### Sounding Output

MPAS can output time series at specified locations (similar to WRF's tslist):

```
&soundings
    config_sounding_interval = '01:00:00'
/
```

Define locations in a `sounding_locations` stream or configuration.

---

## 12. Idealized Simulations

MPAS includes several pre-built idealized test cases for model verification and theoretical studies.

### Available Cases

| Case Number | Name | Description |
|---|---|---|
| 1 | Jablonowski-Williamson (no pert) | Steady-state baroclinic wave on sphere — tests model equilibrium. |
| 2 | Jablonowski-Williamson (perturbed) | Baroclinic wave with perturbation — tests growth of baroclinic instability. Standard dynamical core test. |
| 3 | Jablonowski-Williamson (normal mode) | Normal-mode perturbation variant. |
| 4 | Squall line | 2D squall line on a doubly-periodic Cartesian plane. Tests cold pool dynamics and convective organization. |
| 5 | Supercell | 3D supercell thunderstorm on a doubly-periodic Cartesian plane. Tests microphysics and convection. |
| 6 | Mountain wave | 2D flow over a mountain in the xz-plane. Tests terrain-forced gravity waves. |

### Running an Idealized Case

```bash
# 1. Compile init_atmosphere
make -j4 gfortran CORE=init_atmosphere PRECISION=single

# 2. Configure for idealized case (e.g., case 2: baroclinic wave)
# namelist.init_atmosphere:
# &nhyd_model
#     config_init_case = 2
#     config_start_time = '0000-01-01_00:00:00'
# /
# &dimensions
#     config_nvertlevels = 26
# /

# 3. Run init to create initial conditions
mpiexec -n 1 ./init_atmosphere_model

# 4. Compile and run atmosphere
make clean CORE=atmosphere
make -j4 gfortran CORE=atmosphere PRECISION=single
mpiexec -n 8 ./atmosphere_model
```

---

## 13. Restart & Cycling

### Restart Runs

MPAS writes restart files at configurable intervals. To continue from a restart:

**In `streams.atmosphere`:**
```xml
<immutable_stream name="restart"
                  type="input;output"
                  filename_template="restart.$Y-$M-$D_$h.$m.$s.nc"
                  input_interval="initial_only"
                  output_interval="1_00:00:00" />
```

**In `namelist.atmosphere`:**
```
&restart
    config_do_restart = true
/
&nhyd_model
    config_start_time = '2024-09-02_00:00:00'   ! Must match restart file time
    config_run_duration = '3_00:00:00'           ! Remaining duration
/
```

**When restarting:**
- The `input` stream reads the restart file instead of `init.nc`
- Set `config_start_time` to match the restart file timestamp exactly
- Restart files contain the full model state — no need to re-run initialization

### Data Assimilation Cycling (IAU)

MPAS supports Incremental Analysis Update for data assimilation cycling:

```
&iau
    config_IAU_option = 'on'
    config_IAU_window_length_s = 21600.0   ! 6-hour window
/
```

**IAU workflow:**
```
Cycle 1: init.nc -> atmosphere_model (0-6h) -> restart at 6h
              |
         Analysis (e.g., MPAS-JEDI) -> analysis increment
              |
Cycle 2: restart(6h) + IAU increment -> atmosphere_model (6-12h)
              ...
```

MPAS-JEDI (Joint Effort for Data assimilation Integration) provides the variational and ensemble data assimilation framework for MPAS, built on JCSDA's JEDI system.

---

## 14. Troubleshooting & Tips

### Common Errors

| Error | Cause | Fix |
|---|---|---|
| `graph.info.part.N not found` | Missing decomposition file for N MPI tasks | Run `gpmetis` to generate, or use a pre-computed partition |
| `config_nfglevels mismatch` | Wrong number of first-guess levels | Check intermediate file with `rd_intermediate.exe` and update namelist |
| `Model blowup / NaN` | Timestep too large or bad input data | Reduce `config_dt`, check intermediate files for NaN/missing data |
| Static interp crashes with > 1 task | Static interp must be serial | Use `mpiexec -n 1` for static field interpolation |
| `io_type` error for large meshes | Classic NetCDF 4 GB limit | Set `io_type = "pnetcdf,cdf5"` or `"netcdf4"` in streams |
| `Negative moisture` | Advection not monotonic | Set `config_monotonic = true` |
| Physics lookup tables not found | Missing `*.TBL`, `RRTMG_*` files | Copy from `MPAS-Model/src/core_atmosphere/physics/physics_wrf/files/` |

### Performance Tips

- **Timestep:** Start conservative (4 × dx_km), increase toward (6 × dx_km) if stable
- **Partition quality:** `gpmetis -minconn -contig` gives better load balance
- **Single precision:** `PRECISION=single` is ~30% faster with negligible accuracy loss for NWP
- **Output frequency:** Reduce history output frequency for long runs — I/O can dominate runtime
- **Split transport:** `config_split_dynamics_transport = true` allows fewer transport calls (saves ~15%)
- **PIO tuning:** Set `config_pio_num_iotasks` and `config_pio_stride` for optimal parallel I/O

### Monitoring a Run

```bash
# Watch the log file
tail -f log.atmosphere.0000.out

# Check for errors
grep -i "error\|abort\|nan" log.atmosphere.0000.err

# The log shows timing information:
# "Timing for integration step: X.XX s"
```

---

## 15. Version History

| Version | Date | Key Changes |
|---|---|---|
| **v8.3.1** | June 17, 2025 | Bug-fix release (latest stable) |
| **v8.3.0** | June 2, 2025 | Major release |
| **v8.2.3** | May 22, 2025 | Bug-fix |
| **v8.2.2** | Sep 20, 2024 | Bug-fix, portability, MPAS-JEDI compatibility |
| **v8.2.0** | June 27, 2024 | Major release |
| **v8.1.0** | April 18, 2024 | Major release |
| **v8.0.0** | June 16, 2023 | Noah updated to WRF 4.5, CCPP physics, Simple MPAS I/O Layer (SMIOL), parallel static field remapping, CAM-MPAS integration |
| **v7.0** | June 8, 2019 | Limited-area (regional) simulations, physics updated to WRF 4.0.3, GMTED2010 terrain default, MODIS land use default |
| **v6.0** | April 17, 2018 | Default levels increased 41→55, GMTED2010 terrain support, parallel init builds, new logging module |
| **v5.0** | January 7, 2017 | Convection-permitting suite, diagnostics framework, IAU module, sounding output, CAPE/CIN diagnostics |
| **v4.0** | May 22, 2015 | Split dynamics-transport, physics suites concept, ERA-Interim support, io_type per stream |
| **v3.0** | November 18, 2014 | XML-based I/O configuration (streams), multi-core compilation |
| **v1.0** | June 14, 2013 | Initial public release |

---

## References

- [MPAS Home Page](https://mpas-dev.github.io/)
- [MPAS-Atmosphere](https://mpas-dev.github.io/atmosphere/atmosphere.html)
- [MPAS-Atmosphere Mesh Downloads](https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html)
- [MPAS-Atmosphere User's Guide v8.3.0](https://www2.mmm.ucar.edu/projects/mpas/mpas_atmosphere_users_guide_8.3.0.pdf)
- [MPAS GitHub Repository](https://github.com/MPAS-Dev/MPAS-Model)
- [MPAS-Tools](https://github.com/MPAS-Dev/MPAS-Tools)
- [MPAS Limited-Area Tool](https://github.com/MiCurry/MPAS-Limited-Area)
- [convert_mpas](https://github.com/mgduda/convert_mpas)
- [MPAS-Atmosphere Release Notes](https://mpas-dev.github.io/atmosphere/mpas-a_release_notes.html)
- [MPAS Tutorial (Howard University 2024)](https://www2.mmm.ucar.edu/projects/mpas/tutorial/Howard2024/index.html)
- [MPAS Tutorial (Boulder 2019)](https://www2.mmm.ucar.edu/projects/mpas/tutorial/Boulder2019/index.html)
- [WRF & MPAS-A Support Forum](https://forum.mmm.ucar.edu/)
- [NCAR MPAS Page](https://www.mmm.ucar.edu/models/mpas)
- [MPAS Visualization Scripts](https://mpas-dev.github.io/atmosphere/visualization.html)
- [MPAS-Atmosphere Test Cases (v7.0)](https://www2.mmm.ucar.edu/projects/mpas/test_cases/v7.0/)
- [MPAS Static Data](https://www2.mmm.ucar.edu/projects/mpas/mpas_static.tar.bz2)
- Skamarock et al. (2012): "A Multiscale Nonhydrostatic Atmospheric Model Using Centroidal Voronoi Tesselations and C-Grid Staggering", MWR
