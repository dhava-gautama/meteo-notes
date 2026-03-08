# OpenIFS: Complete Guide

A comprehensive guide to ECMWF's **OpenIFS** — the open-source version of the Integrated Forecasting System (IFS) for research and education.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Access and Licensing](#2-access-and-licensing)
3. [System Requirements](#3-system-requirements)
4. [Installation](#4-installation)
5. [Configuration](#5-configuration)
6. [Running 3D NWP Experiments](#6-running-3d-nwp-experiments)
7. [Single-Column Model (SCM)](#7-single-column-model-scm)
8. [Static Data and Initial Conditions](#8-static-data-and-initial-conditions)
9. [Resolution and Grid](#9-resolution-and-grid)
10. [Physics Overview](#10-physics-overview)
11. [Docker Deployment](#11-docker-deployment)
12. [Testing](#12-testing)
13. [Tips and Troubleshooting](#13-tips-and-troubleshooting)
14. [References](#14-references)

---

## 1. Overview

OpenIFS provides the **same forecast capability as the operational ECMWF IFS**, but **without the data assimilation system**. It is designed for:

- Academic teaching and training
- Research on atmospheric dynamics and physics
- Process studies (via the Single-Column Model)
- Resolution sensitivity experiments
- Physics parameterization development

### What OpenIFS Is and Is Not

| Included | Not included |
|----------|--------------|
| Full 3D global spectral atmospheric model | Data assimilation (4D-Var, EDA) |
| IFS physics parameterizations | Ocean model coupling (NEMO) |
| Single-Column Model (SCM, since 48r1) | Wave model (ECWAM) |
| Optional atmospheric chemistry | Ensemble prediction system |
| Low-resolution test cases | Operational workflow scripts |

### Current Versions

| Version | IFS Cycle | Notes |
|---------|-----------|-------|
| OpenIFS 48r1 | CY48R1 | Current release, includes SCM |
| OpenIFS 43r3 | CY43R3 | Earlier release |

---

## 2. Access and Licensing

### License

OpenIFS is released under the **Apache License 2.0**.

### Source Code

```bash
# GitHub repository
https://github.com/ecmwf-ifs/openifs
```

### Getting Initial Conditions

- **OpenIFS Data Hub** — registered users can download prepared initial conditions
- **Copernicus Climate Data Store (CDS)** — create custom initial conditions from ECMWF reanalyses (ERA5, etc.)
- **ECMWF HPC** — member state users can request computing access through their National Meteorological Services

---

## 3. System Requirements

### Supported Platforms

- **Linux** — officially supported and tested
- **macOS** — may work with proper dependencies (Docker recommended)
- Other UNIX-like systems may work but are not guaranteed

### Dependencies

| Category | Packages |
|----------|----------|
| Build tools | `git`, `cmake`, `bison`, `flex` |
| Compilers | GNU (gcc/gfortran), Intel (icc/ifort), or Cray |
| MPI | OpenMPI or Intel MPI |
| Python | Python 3 with `ruamel.yaml` and `yaml` |
| Numerical | LAPACK, ATLAS, Eigen3 |
| Data I/O | NetCDF (C and Fortran), HDF5 |
| Boost | date_time, filesystem, serialization, program_options |

Install on Ubuntu/Debian:

```bash
sudo apt-get install git cmake bison flex openmpi-bin libopenmpi-dev \
    python3 python3-yaml gfortran liblapack-dev libatlas-base-dev \
    libeigen3-dev libnetcdf-dev libnetcdff-dev libboost-all-dev
```

Install Python dependencies:

```bash
pip3 install ruamel.yaml pyyaml
```

---

## 4. Installation

### Clone the Repository

```bash
# Shallow clone (recommended — smaller download)
git clone --depth 1 https://github.com/ecmwf-ifs/openifs.git
cd openifs

# Or clone a specific release branch
git clone --depth 1 --branch release/48r1 https://github.com/ecmwf-ifs/openifs.git

# Or clone a specific tag (create branch to enable commits)
git clone --depth 1 --branch v48r1.0 https://github.com/ecmwf-ifs/openifs.git
cd openifs
git checkout -b my-branch
```

### Directory Structure

```
openifs/
├── ifs-source/          # Model source code
├── ifs-test/            # Low-resolution test cases (T21)
├── scripts/             # Build, test, and execution utilities
├── arch/ecmwf/          # HPC-specific configuration files
├── oifs-config.edit_me.sh   # User configuration (edit this)
└── openifs-test.sh      # Main build/test driver script
```

### Configure

Edit `oifs-config.edit_me.sh` — the central configuration file:

```bash
cp oifs-config.edit_me.sh oifs-config.sh
vi oifs-config.sh
```

Key variables to set:

```bash
# Root directory of OpenIFS source
export OIFS_HOME=/path/to/openifs

# Where static data files are stored
export OIFS_DATA_DIR=/path/to/oifs-data

# Compiler selection
export OIFS_COMPILER=gnu    # or intel, cray

# MPI configuration
export OIFS_NPROC=4         # Number of MPI tasks
export OIFS_NTHREAD=1       # Number of OpenMP threads
```

Source the configuration before building or running:

```bash
source oifs-config.sh
```

### Build

The build process uses the `openifs-test.sh` driver script:

```bash
# Create source directory from bundle
./openifs-test.sh -c

# Compile
./openifs-test.sh -b

# Run test suite (T21 low-resolution cases)
./openifs-test.sh -t

# All steps together
./openifs-test.sh -c -b -t
```

A successful build produces:

```
Good news - ctest has passed, openifs is ready for experiment and SCM testing.
```

The test suite runs **21 3D NWP experiments** and **1 SCM test** at T21 resolution.

---

## 5. Configuration

### Experiment Configuration (exp-config.h)

Each experiment directory contains an `exp-config.h` file with key run parameters:

```bash
# Experiment identifier (4 characters)
EXPERIMENT_ID="ab7z"

# Spectral resolution
OIFS_RES="T159"

# Parallel decomposition
OIFS_NPROC=16          # MPI tasks
OIFS_NTHREAD=2         # OpenMP threads per task

# Forecast length
OIFS_FCLEN=240         # Hours

# Postprocessing
OIFS_PP=1              # Enable (1) or disable (0)

# Custom namelist
OIFS_NAMELIST="fort.4"
```

### Namelist (fort.4)

The IFS uses Fortran namelists for detailed configuration. Key namelist groups:

| Namelist | Controls |
|----------|----------|
| `NAMDYN` | Dynamics (timestep, diffusion) |
| `NAMPHY` | Physics switches |
| `NAERAD` | Radiation |
| `NAMCUMF` | Convection |
| `NAMPPC` | Postprocessing / output fields |
| `NAMGFL` | Prognostic GFL (General Field Layout) variables |

---

## 6. Running 3D NWP Experiments

### Setting Up an Experiment

Download an example experiment (e.g., Hurricane Karl 2016 at T159):

```bash
# Download and extract experiment data
cd $OIFS_DATA_DIR
tar xvf ab7z_T159.tar.gz

# Copy run scripts to experiment directory
cp $OIFS_HOME/scripts/oifs-run $OIFS_DATA_DIR/ab7z/
cp $OIFS_HOME/scripts/exp-config.h $OIFS_DATA_DIR/ab7z/
```

Edit `exp-config.h` for your experiment settings (resolution, MPI tasks, forecast length).

### Interactive Run

```bash
source oifs-config.sh
cd $OIFS_DATA_DIR/ab7z
./oifs-run
```

### Batch (SLURM) Run

```bash
# Copy the batch wrapper script
cp $OIFS_HOME/scripts/run-oifs.ecmwf-hpc2020.job $OIFS_DATA_DIR/ab7z/

# Edit SLURM directives (nodes, walltime, account)
vi run-oifs.ecmwf-hpc2020.job

# Submit
sbatch ./run-oifs.ecmwf-hpc2020.job
```

The batch script reads `exp-config.h` but overrides `LAUNCH`, `OIFS_NPROC`, and `OIFS_NTHREAD` with batch-specific values.

### Output

OpenIFS produces GRIB output files. Postprocessing converts spectral fields to gridpoint fields at the desired output resolution.

---

## 7. Single-Column Model (SCM)

Available since **OpenIFS 48r1**, the SCM isolates the IFS physics column for process studies. It is built automatically alongside the 3D model.

### SCM Environment Variables

Set in `oifs-config.sh`:

```bash
# SCM executable (single or double precision)
export SCM_EXE=$OIFS_HOME/build/bin/master_scm_dp   # double precision

# SCM test directory
export SCM_TEST=$OIFS_HOME/ifs-test/scm

# SCM project directory (where case data lives)
export SCM_PROJ=/path/to/scm-project

# SCM run directory (working directory for output)
export SCM_RUN=/path/to/scm-runs
```

### Standard Test Cases

Download the SCM test package (~45 MB):

```bash
wget https://openifs.ecmwf.int/data/scm_openifs_48r1.tar.gz
tar xvf scm_openifs_48r1.tar.gz -C $SCM_PROJ
```

| Case | Description | Science Focus |
|------|-------------|---------------|
| **DYCOMS** | Marine stratocumulus off California | Boundary layer, cloud-top entrainment |
| **BOMEX** | Trade-wind cumulus (Barbados) | Shallow convection, cloud base mass flux |
| **TWPICE** | Tropical Warm Pool (Darwin, Australia) | Deep convection, multi-day tropical forcing |

### Running the SCM

The `callscm` wrapper script provides a convenient interface:

```bash
# Run all three cases with default settings (dt=450s)
$SCM_TEST/callscm

# Run a single case
$SCM_TEST/callscm -c BOMEX

# Run with specific timesteps
$SCM_TEST/callscm -c BOMEX -t "1800 900 450"

# Custom experiment name
$SCM_TEST/callscm -c BOMEX -t "1800 900" -x "bomex_test"
```

**Command-line options:**

| Flag | Description | Default |
|------|-------------|---------|
| `-c <case>` | Case name (DYCOMS, BOMEX, TWPICE) | All cases |
| `-t <dt>` | Timestep(s) in seconds (space-separated) | 450 |
| `-x <name>` | Experiment name | ref-oifs-scm |

### SCM Output

Output directories follow the naming convention:

```
scmout_<CASE>_<EXPT>_<TIMESTEP>/
```

For example: `scmout_BOMEX_bomex_test_900/`

### HPC Notes for SCM

On ECMWF systems, MPI modules must be loaded before running the SCM:

```bash
# Intel compiler stack
module load intel-mpi

# GNU/OpenMPI
module load openmpi
```

---

## 8. Static Data and Initial Conditions

### Static Data

Download required static files to `$OIFS_DATA_DIR`:

```bash
# Radiation lookup tables (required)
wget https://openifs.ecmwf.int/data/rtables.tar.gz
tar xvf rtables.tar.gz -C $OIFS_DATA_DIR

# Resolution-independent IFS data (required)
wget https://openifs.ecmwf.int/data/ifsdata.tar.gz
tar xvf ifsdata.tar.gz -C $OIFS_DATA_DIR

# Resolution-dependent climate files (per resolution)
wget https://openifs.ecmwf.int/data/48r1_climate.v020_159.tar.gz  # T159
tar xvf 48r1_climate.v020_159.tar.gz -C $OIFS_DATA_DIR
```

### Initial Conditions from ERA5

For custom experiments, create initial conditions from ERA5 via the Copernicus Climate Data Store:

```python
import cdsapi

c = cdsapi.Client()

# Download pressure-level fields
c.retrieve('reanalysis-era5-pressure-levels', {
    'product_type': 'reanalysis',
    'variable': ['geopotential', 'temperature', 'u_component_of_wind',
                 'v_component_of_wind', 'specific_humidity'],
    'pressure_level': ['1', '2', '3', '5', '7', '10', '20', '30', '50',
                       '70', '100', '125', '150', '175', '200', '225',
                       '250', '300', '350', '400', '450', '500', '550',
                       '600', '650', '700', '750', '775', '800', '825',
                       '850', '875', '900', '925', '950', '975', '1000'],
    'year': '2024', 'month': '01', 'day': '01',
    'time': '00:00',
    'format': 'grib',
}, 'era5_pl.grib')

# Download surface fields
c.retrieve('reanalysis-era5-single-levels', {
    'product_type': 'reanalysis',
    'variable': ['skin_temperature', 'sea_surface_temperature',
                 'soil_temperature_level_1', 'soil_temperature_level_2',
                 'volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_2',
                 'snow_depth', '10m_u_component_of_wind',
                 '10m_v_component_of_wind', '2m_temperature'],
    'year': '2024', 'month': '01', 'day': '01',
    'time': '00:00',
    'format': 'grib',
}, 'era5_sfc.grib')
```

---

## 9. Resolution and Grid

OpenIFS uses a **spectral transform method** with a reduced Gaussian grid.

### Common Resolutions

| Spectral | Grid Spacing (approx.) | Grid Points | Typical Use |
|----------|----------------------|-------------|-------------|
| T21 | ~500 km | ~8,000 | Testing only |
| T159 | ~125 km | ~35,000 | Teaching, quick experiments |
| T255 | ~80 km | ~88,000 | Research |
| T511 | ~40 km | ~350,000 | High-resolution research |
| T639 | ~30 km | ~540,000 | Near-operational |
| T1279 | ~16 km | ~2,100,000 | Operational-equivalent |
| TCo1279 | ~9 km (cubic octahedral) | ~6,600,000 | Current ECMWF operational |

### Vertical Levels

OpenIFS typically uses **91 levels** (L91) up to 0.01 hPa (~80 km), or **137 levels** (L137) in newer cycles.

---

## 10. Physics Overview

OpenIFS inherits the full IFS physics suite. Key parameterizations:

| Component | Scheme | Notes |
|-----------|--------|-------|
| **Radiation** | ecRad | Shortwave (McICA-RRTMGP) and longwave |
| **Convection** | Bechtold (2008, 2014) | Mass-flux, deep + shallow |
| **Cloud microphysics** | Single-moment bulk | Prognostic cloud condensate |
| **Cloud scheme** | Tiedtke-Bechtold | Statistical cloud fraction |
| **Boundary layer** | EDMF (Eddy-Diffusivity Mass-Flux) | Unified dry/moist turbulence |
| **Surface** | HTESSEL / ECLand | Tiled land surface, multi-layer soil |
| **Gravity wave drag** | Orographic + non-orographic | Lott & Miller (1997), Scinocca (2003) |
| **Ozone** | Prognostic or climatology | Optional linearized chemistry |
| **Aerosols** | Climatological (Tegen) or prognostic | Optional CAMS integration |

---

## 11. Docker Deployment

For systems where dependency management is difficult, OpenIFS provides a Docker build:

```bash
# Generate a Dockerfile
python3 create-oifs-docker.py

# Build the Docker image
docker build -t openifs .

# Run interactively
docker run -it openifs /bin/bash
```

This automates the entire installation, compilation, and test process. Tested on macOS and Linux.

---

## 12. Testing

### Built-in Tests

```bash
# Run the full test suite after building
./openifs-test.sh -t
```

The test suite:
- Runs **21 low-resolution (T21) 3D NWP experiments**
- Runs **1 SCM test case**
- Confirms the model compiles and runs to completion
- Does **not** guarantee bit-reproducibility with reference outputs (relevant for code modification validation)

### Verifying a Custom Run

Check for successful completion in the model output log:

```bash
# Look for "CNT0 - FORECAST IS OVER" in the output
grep "FORECAST IS OVER" NODE.001_01
```

---

## 13. Tips and Troubleshooting

### Common Issues

| Problem | Solution |
|---------|----------|
| Missing MPI at runtime | Load MPI module: `module load openmpi` or `module load intel-mpi` |
| CMake can't find Eigen3 | Install: `apt-get install libeigen3-dev` or set `Eigen3_DIR` |
| NetCDF not found | Ensure both C and Fortran NetCDF are installed |
| Detached HEAD after tag clone | `git checkout -b my-branch` |
| SCM fails on HPC | Ensure MPI environment is loaded before execution |
| Build fails with GNU compiler | Ensure gfortran >= 9.0; check Boost version compatibility |

### Performance Tips

- **MPI vs OpenMP**: Start with pure MPI (NTHREAD=1) and increase threads only if memory-limited
- **T159 on a laptop**: Feasible with 4 MPI tasks, ~10 min for 10-day forecast
- **T511+**: Requires HPC; typical scaling: 32-128 MPI tasks
- **I/O**: Output frequency is the main bottleneck at high resolution — reduce output fields or frequency

### Indonesian / Tropical Applications

OpenIFS is a **global** model, so Indonesian domain coverage is automatic at any resolution. For tropical Maritime Continent studies:

- T159 (~125 km) is too coarse for island-scale features
- T511 (~40 km) begins to resolve major Indonesian islands
- T1279 (~16 km) or TCo1279 (~9 km) needed for strait-scale flows (Makassar, Malacca)
- Compare against BMKG's operational GFS-based forecasts or ECMWF HRES
- SCM cases (TWPICE) are directly relevant to Maritime Continent deep convection

---

## 14. References

- **Source code**: <https://github.com/ecmwf-ifs/openifs>
- **OpenIFS home**: <https://openifs.ecmwf.int/>
- **Documentation (Confluence)**: <https://confluence.ecmwf.int/display/OIFS>
- **User forum**: <https://forum.ecmwf.int/c/ifs-aifs-and-openifs/openifs/45>
- **Static data downloads**: <https://openifs.ecmwf.int/data/>
- **Copernicus CDS (ERA5)**: <https://cds.climate.copernicus.eu/>
- **IFS Documentation (CY48R1)**: <https://www.ecmwf.int/en/publications/ifs-documentation>

### Key Publications

- ECMWF, "IFS Documentation CY48R1" — full physics and dynamics reference
- Bechtold et al. (2008, 2014) — convection parameterization
- Lott & Miller (1997) — orographic gravity wave drag
- Pincus & Stevens (2013) — ecRad radiation scheme
