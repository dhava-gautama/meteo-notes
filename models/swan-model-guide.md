# SWAN: Complete Guide from Grid Setup to Visualization

> A practical guide covering the full SWAN (Simulating WAves Nearshore) workflow for coastal and nearshore wave modeling.
> Covers spectral wave theory, all physics options, grid types, input/output, nesting, calibration, and running strategies.

---

## Table of Contents

1. [Overview & Key Concepts](#1-overview--key-concepts)
2. [Installation & Compilation](#2-installation--compilation)
3. [Grid Types](#3-grid-types)
4. [Physics Options Reference](#4-physics-options-reference)
5. [Numerical Methods](#5-numerical-methods)
6. [Input File Structure](#6-input-file-structure)
7. [Output Types & Variables](#7-output-types--variables)
8. [Nesting](#8-nesting)
9. [Unstructured Grids](#9-unstructured-grids)
10. [Pre/Post-Processing & Visualization](#10-prepost-processing--visualization)
11. [Running SWAN](#11-running-swan)
12. [Calibration & Validation](#12-calibration--validation)

---

## 1. Overview & Key Concepts

### What Is SWAN?

SWAN (Simulating WAves Nearshore) is a third-generation spectral wave model developed at Delft University of Technology (TU Delft), Netherlands. It computes random, short-crested, wind-generated waves in coastal regions, lakes, and estuaries. SWAN is based on the discrete spectral action balance equation and accounts for all relevant physical processes in nearshore wave evolution.

**Current version:** SWAN Cycle III version 41.51

**Developer:** Environmental Fluid Mechanics section, Faculty of Civil Engineering and Geosciences, TU Delft

**License:** GNU General Public License (GPL). Source code available from SourceForge and TU Delft GitLab.

**Key references:**
- Booij, N., Ris, R.C. and Holthuijsen, L.H. (1999). A third-generation wave model for coastal regions. Part I: Model description and validation. *J. Geophys. Res.*, 104(C4), 7649-7666.
- Ris, R.C., Holthuijsen, L.H. and Booij, N. (1999). A third-generation wave model for coastal regions. Part II: Verification. *J. Geophys. Res.*, 104(C4), 7667-7681.

### Spectral Action Balance Equation

The fundamental equation governing SWAN is the spectral action balance equation. Action density `N(sigma, theta)` is defined as:

```
N(sigma, theta) = E(sigma, theta) / sigma
```

where `E` is the energy density spectrum and `sigma` is the relative (intrinsic) radian frequency. Action density is conserved during propagation in the presence of ambient currents, whereas energy density is not.

**General form:**

```
dN/dt + d(cx*N)/dx + d(cy*N)/dy + d(c_sigma*N)/d_sigma + d(c_theta*N)/d_theta = S_tot / sigma
```

| Term | Meaning |
|------|---------|
| `dN/dt` | Local rate of change of action density in time |
| `d(cx*N)/dx + d(cy*N)/dy` | Propagation in 2D geographic space (includes shoaling) |
| `d(c_sigma*N)/d_sigma` | Frequency shifting due to depth and current variations |
| `d(c_theta*N)/d_theta` | Depth-induced and current-induced refraction |
| `S_tot/sigma` | Source/sink terms: wind input, dissipation, nonlinear interactions |

**Dispersion relation:**

```
sigma^2 = g * |k| * tanh(|k| * d)
```

**Source terms (`S_tot`) include:**

```
S_tot = S_in + S_nl3 + S_nl4 + S_ds,w + S_ds,b + S_ds,br + S_ds,veg + S_ds,mud + S_ds,turb + S_ds,ice
```

| Symbol | Process |
|--------|---------|
| `S_in` | Wind input |
| `S_nl3` | Triad nonlinear interactions |
| `S_nl4` | Quadruplet nonlinear interactions |
| `S_ds,w` | Whitecapping dissipation |
| `S_ds,b` | Bottom friction dissipation |
| `S_ds,br` | Depth-induced breaking dissipation |
| `S_ds,veg` | Vegetation dissipation |
| `S_ds,mud` | Mud-induced dissipation |
| `S_ds,turb` | Turbulent viscosity dissipation |
| `S_ds,ice` | Sea ice dissipation |

### Phase-Averaging vs Phase-Resolving

SWAN is a **phase-averaged** (spectral) model. It does not resolve individual wave crests and troughs. Instead, it describes the wave field statistically through the energy/action density spectrum. This approach:

- Is computationally efficient for large domains
- Cannot resolve individual wave interference, standing wave patterns, or diffraction shadow zones in full detail
- Handles diffraction only in a phase-decoupled approximate manner
- Is suitable for domains from tens of meters to hundreds of kilometers

Phase-resolving models (e.g., BOSZ, SWASH, Boussinesq models) track individual wave phases but are computationally far more expensive.

### Comparison with WAM and WAVEWATCH III

| Feature | SWAN | WAM | WAVEWATCH III |
|---------|------|-----|---------------|
| Primary application | Nearshore, coastal, estuaries, lakes | Deep ocean, global scale | Ocean basins, global to regional |
| Numerical scheme | Implicit (unconditionally stable) | Explicit | Explicit (various options) |
| Shallow-water physics | Comprehensive (triads, breaking, friction, vegetation, mud) | Limited | Added in later versions |
| Grid types | Regular, curvilinear, unstructured | Regular | Regular, curvilinear, unstructured, SMC |
| Parallelization | MPI, OpenMP | MPI | MPI, OpenMP, hybrid |
| Source term packages | Multiple (KOMEN, JANSSEN, WESTH, ST6) | Single set | Multiple packages (ST1-ST6) |
| Efficiency at ocean scale | Not recommended | Excellent | Excellent |
| Multi-grid nesting | One-way (consecutive runs) | One-way | Full two-way mosaic |
| Developer | TU Delft | ECMWF consortium | NOAA/NCEP |

**Key differences:**
- WAM and WW3 were designed for large-scale (ocean/global) applications; SWAN was built specifically for nearshore
- SWAN's implicit numerical scheme allows larger time steps in shallow water without Courant number restrictions
- SWAN contains formulations specifically for shallow water that WAM and WW3 historically lacked
- SWAN can be nested within WAM or WW3 for downscaling from ocean to coast

---

## 2. Installation & Compilation

### Source Code Download

SWAN source code is available from:
- **SourceForge:** https://sourceforge.net/projects/swanmodel/files/swan/
- **TU Delft GitLab:** https://gitlab.tudelft.nl/citg/wavemodels/swan (since version 41.41)

```bash
tar xzf swan4151.tar.gz
```

### Dependencies

| Dependency | Required? | Notes |
|------------|-----------|-------|
| Fortran compiler (ANSI F90+) | Yes | gfortran and Intel Fortran (ifort/ifx) supported |
| Perl | Yes | Required for the `switch.pl` preprocessor script |
| MPI library | For parallel distributed | OpenMPI, MPICH, Intel MPI |
| OpenMP support | For parallel shared memory | Built into modern Fortran compilers |
| NetCDF library (>=4.5.x) | Optional | For NetCDF input/output; Fortran interface required |
| METIS library | Optional | For domain decomposition with unstructured grids in MPI mode |
| GNU Make or CMake | Yes | CMake supported since version 41.41 |

### Compilation Steps (GNU Make)

**Step 1: Generate the configuration file**

```bash
make config
```

This creates `macros.inc` containing machine-dependent settings. Edit to set:
- Compiler name (`F90` variable)
- Compiler flags (`F90_FLAGS`)
- NetCDF root directory (`NETCDFROOT`) if desired
- MPI compiler wrapper (e.g., `mpif90`)

**Step 2: Build SWAN in desired mode**

| Command | Mode | Description |
|---------|------|-------------|
| `make ser` | Serial | Single-processor execution |
| `make omp` | OpenMP | Shared-memory parallel (multi-threaded) |
| `make mpi` | MPI | Distributed-memory parallel (multi-process) |

**Cleanup commands:**
- `make clean` -- removes object files and modules
- `make clobber` -- removes all generated files from previous builds

**Step 3 (alternative): Build with CMake**

```bash
mkdir build && cd build
cmake ..
make
```

### Building with NetCDF Support

Set the `NETCDFROOT` variable in `macros.inc`:

```
NETCDFROOT = /usr/local/netcdf
```

Ensure the NetCDF Fortran interface libraries (`libnetcdff`) are accessible.

### The switch.pl Preprocessor

The `switch.pl` Perl script converts `.ftn` and `.ftn90` source files to standard Fortran based on platform and mode flags. It is invoked automatically by the Makefile. Manual usage:

```bash
perl switch.pl -unix -f95 -mpi *.ftn *.ftn90
```

Flags: `-unix`/`-win32`, `-f95`, `-mpi`/`-omp`/`-ser`

### Platform-Specific Notes

| Platform | Notes |
|----------|-------|
| Linux (gfortran) | Most common; use `make config`, edit `macros.inc` for gfortran |
| Linux (Intel) | Use ifort or ifx; supports OpenMP and Intel MPI |
| Windows | Requires Nmake (from Visual Studio); Intel Fortran recommended |
| macOS | gfortran via Homebrew; OpenMP may require additional flags |
| HPC clusters | Use module system to load compiler, MPI, NetCDF; submit via job scheduler |

---

## 3. Grid Types

### 3.1 Regular Rectangular Grid

The simplest grid type. Defined by origin, orientation, length, and number of cells.

```
CGRID REGULAR [xpc] [ypc] [alpc] [xlenc] [ylenc] [mxc] [myc] ...
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `xpc`, `ypc` | Origin coordinates (m or degrees) | 0.0 |
| `alpc` | Grid orientation (degrees from x-axis) | 0.0 |
| `xlenc`, `ylenc` | Grid extent in x and y directions | (required) |
| `mxc`, `myc` | Number of mesh cells in x and y | (required) |

No special grid generation tools needed; defined directly in the input file.

### 3.2 Curvilinear Grid

Grid lines follow curved coordinate directions; useful for domains with curved coastlines, rivers, and estuaries. Grid cells are quadrilaterals.

```
CGRID CURVILINEAR [mxc] [myc] (EXCEPTION [xexc] [yexc]) ...
```

Coordinates supplied via a separate file read with `READGRID COORDINATES`.

**Grid generation tools:**
- **RGFGRID** (Deltares): Dedicated orthogonal curvilinear grid generator within Delft3D
- **Delft Dashboard**: GUI tool for creating Delft3D/SWAN grids
- **SMS (Surface-water Modeling System)**: Commercial mesh generation by Aquaveo

### 3.3 Unstructured Triangular Grid

Composed of triangular elements with variable resolution. Ideal for complex coastlines, islands, and varying bathymetry.

```
CGRID UNSTRUCTURED ...
```

Mesh supplied via external mesh files (nodes, elements, boundary markers).

**Mesh generation tools:**
- **Triangle** (Jonathan Shewchuk): Free, robust Delaunay triangulator
- **Gmsh**: Open-source (GPL) finite element mesh generator with GUI
- **BatTri**: MATLAB-based mesh generator for coastal applications
- **SMS**: Commercial, supports SWAN-format output

**Key considerations:**
- Only triangular elements supported in current SWAN versions
- Vertex-based, fully implicit finite difference method
- Stability is unconditional for any time step
- METIS library recommended for domain decomposition in MPI mode

### 3.4 Spectral Grid (All Grid Types)

The spectral domain is defined as part of the `CGRID` command:

```
... CIRCLE [mdc] [flow] [fhigh] [msc]
... SECTOR [dir1] [dir2] [mdc] [flow] [fhigh] [msc]
```

| Parameter | Description | Recommendation |
|-----------|-------------|----------------|
| `CIRCLE` | Full 360-degree directional coverage | Default; use for most cases |
| `SECTOR [dir1] [dir2]` | Limited directional range | Pad 30 degrees beyond expected spectrum |
| `mdc` | Number of directional bins | Minimum 3 per quadrant; 36 (10-degree bins) is common |
| `flow` | Lowest frequency (Hz) | >= 0.03 Hz recommended |
| `fhigh` | Highest frequency (Hz) | Typically 1.0 Hz |
| `msc` | Number of frequencies minus 1 | Minimum 4; ~30-40 is typical |

Frequencies are distributed logarithmically with ratio gamma = 1.1 (optimal for quadruplet interactions).

---

## 4. Physics Options Reference

### 4.1 Generation Mode Selection

SWAN supports first-, second-, and third-generation modes. The generation mode command activates wind input and whitecapping automatically.

| Command | Mode | Wind Input | Whitecapping |
|---------|------|-----------|--------------|
| `GEN1` | 1st generation | Linear + exponential | Energy ceiling |
| `GEN2` | 2nd generation | Linear + exponential | Spectral limiter |
| `GEN3` | 3rd generation | Physics-based | Physics-based |

Third generation (`GEN3`) is standard for modern applications:

```
GEN3 JANSSEN [cds1] [delta]
GEN3 KOMEN [cds2] [stpm]
GEN3 WESTH                     $ default since v41.45
GEN3 ST6 [a1sds] [a2sds] [p1sds] [p2sds] ...
```

### 4.2 Wind Input Formulations

**KOMEN (Komen et al., 1984):**
- Exponential growth based on WAM Cycle 3 formulation
- Linear growth via Cavaleri & Malanotte-Rizzoli (1981), activated with `AGROW` keyword
- Key coefficient: `cds2` (whitecapping tied to wind input)

**JANSSEN (Janssen, 1989, 1991):**
- Quasi-linear wind-wave interaction theory
- Accounts for feedback of waves on the wind profile through sea surface roughness
- Key coefficients: `cds1 = 4.5`, `delta = 0.5`

**WESTH (Van der Westhuysen et al., 2007):**
- Nonlinear saturation-based approach for whitecapping
- Uses Yan (1987) wind input expression
- **Default since SWAN version 41.45**
- Exponent transitions between weak wind (`p0=2`) and strong wind (`p0=4`) based on `u*/c` ratio

**ST6 (Rogers et al., 2012):**
- Observation-based source term package
- Wind input combines linear and quadratic growth regimes
- Wind scaling: `U = S_ws * u_star`, default `S_ws = 32`
- Drag formulas: `HWANG` (default), `FAN`, `ECMWF`
- Stress calculation: `VECTAU` (vector, default) or `SCATAU` (scalar)

**DRAG options (wind drag coefficient):**
- `WU`: Wu (1982) formula (default since v41.45)
- `FIT`: Second-order polynomial fit (Zijlema et al., 2012)

### 4.3 Whitecapping Dissipation

**KOMEN formulation (`WCAPPING KOMEN`):**

```
S_ds,w(sigma,theta) = -Gamma * sigma~ * (k/k~) * E(sigma,theta)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `cds2` | Dissipation rate coefficient | 2.36e-5 |
| `stpm` | PM spectrum steepness | 3.02e-3 |
| `powst` | Power of normalized steepness | 2 |
| `delta` | Wave number dependency | 1 |
| `powk` | Power of normalized wave number | 1 |

**Alves & Banner (AB) formulation (`WCAPPING AB`):** (default since v41.45)

Based on spectral saturation:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `cds2` | Proportionality coefficient | 5.0e-5 |
| `br` | Threshold saturation level | 1.75e-3 |
| `CURRENT` | Enable current-induced dissipation | off |
| `cds3` | Current-induced coefficient | 0.8 |

### 4.4 Swell Dissipation

```
SSWELL ARDHUIN [cdsv]    $ cdsv default: 1.2
SSWELL ZIEGER [b1]
```

- `ARDHUIN`: Non-breaking dissipation per Ardhuin et al. (2010)
- `ZIEGER`: Per Young et al. (2013) and Zieger et al. (2015)

### 4.5 Bottom Friction

All formulations share a common structure:

```
S_ds,b = -C_b * (sigma^2 / (g^2 * sinh^2(kd))) * E(sigma,theta)
```

**JONSWAP (Hasselmann et al., 1973):**

```
FRICTION JONSWAP CONSTANT [cfjon]
```

| `cfjon` value | Application |
|---------------|-------------|
| 0.038 | Sandy bottoms (default) |
| 0.019 | Smooth seafloors |
| 0.067 | Depth-limited swell conditions |

**COLLINS (Collins, 1972):**

```
FRICTION COLLINS [cfw]     $ cfw default: 0.015
```

The Collins coefficient can vary spatially via `INPGRID`/`READINP FRICTION`.

**MADSEN (Madsen et al., 1988):**

```
FRICTION MADSEN [kn]       $ kn default: 0.05 m
```

Uses equivalent roughness length. The roughness `kn` can vary spatially.

**RIPPLES (Smith et al., 2011):**

```
FRICTION RIPPLES [S] [D]   $ S=2.65, D=0.0001 m
```

### 4.6 Depth-Induced Breaking

**Battjes & Janssen (1978) -- CONSTANT method (default):**

```
BREAKING CONSTANT [alpha] [gamma]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `alpha` | Dissipation rate proportionality | 1.0 |
| `gamma` | Breaker index (H_max/d ratio) | 0.73 |

Field observations show gamma varies between 0.6 and 0.83. Depth-induced breaking is activated by default and should generally not be turned off.

**BKD method (slope-dependent):**

```
BREAKING BKD [alpha] [gamma0] [a1] [a2] [a3]
```

### 4.7 Triad Wave-Wave Interactions

Near-shore nonlinear triad interactions transfer energy between harmonics.

**DCTA (default since v41.45):**

```
TRIAD DCTA [trfac] [p] [COLL|NONC]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `trfac` | Interaction intensity scaling | 4.4 |
| `p` | High-frequency tail shape | 4/3 |
| `COLL` | Restrict to collinear interactions | default |

**Other methods:**
- `LTA` (Lumped Triad Approximation): `TRIAD LTA [trfac]` (trfac default: 1.0)
- `FTIM` (Full Triad Interaction Model, new in v41.51): `TRIAD FTIM [trfac]`
- `SPB` (Spectral Partitioning Bispectrum): `TRIAD SPB [trfac] [a]`

**Biphase parameterization:**
- `ELDEBERKY [urcrit]`: Based on Ursell number; urcrit default = 0.63 (field), 0.2 (lab)
- `DEWIT [lpar]`: De Wit (2022) method with spatial averaging

**Transfer coefficient:**
- `FG`: Freilich & Guza (1984)
- `MS`: Madsen & Sorensen (1993)
- `BREDMOSE`: Exact second-order Stokes
- `QUADWAVE`: Akrish et al. (2024), default

### 4.8 Quadruplet Wave-Wave Interactions

Deep-water nonlinear four-wave interactions; dominant mechanism for spectral evolution in deep water.

```
QUADRUPL [iquad] [lambda] [Cnl4] [Csh1] [Csh2] [Csh3]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `iquad` | Integration method | 2 |
| `lambda` | DIA configuration coefficient | 0.25 |
| `Cnl4` | Proportionality coefficient | 3e7 |

**`iquad` options:**

| Value | Method | Description |
|-------|--------|-------------|
| 1 | DIA (semi-implicit, sweep) | Semi-implicit per sweep |
| 2 | DIA (explicit, sweep) | Default; explicit per sweep |
| 3 | DIA (explicit, iteration) | Per iteration, full explicit |
| 4 | MDIA (Multiple DIA) | More accurate than single DIA |
| 51 | XNL (deep water) | Exact nonlinear deep water |
| 53 | XNL (finite depth) | Exact nonlinear finite depth |

### 4.9 Wave-Current Interaction

SWAN accounts for wave-current interaction through frequency shifting, refraction, and energy exchange. Currents are input via `INPGRID CURRENT` and `READINP CURRENT`.

Current-induced whitecapping can be activated via:

```
WCAPPING AB [cds2] [br] CURRENT [cds3]
```

### 4.10 Diffraction

SWAN uses a phase-decoupled approximation to diffraction:

```
DIFFRACTION [idiffr] [smpar] [smnum] [cgmod]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `idiffr` | 0=none, 1=included | 1 |
| `smpar` | Smoothing parameter | 0 |
| `smnum` | Number of smoothing steps | 0 |
| `cgmod` | Group velocity adaptation (0 or 1) | 1 |

**Recommendations:**
- Use under-relaxation in NUMERIC command (`alpha = 0.01`) when diffraction is active
- Spatial resolution should be 1/5 to 1/10 of the dominant wavelength near obstacles

### 4.11 Vegetation Dissipation

```
VEGETATION [iveg] <[height] [diamtr] [nstems] [drag]>
```

Multiple vertical layers can be specified by repeating the height/diamtr/nstems/drag group.

| Parameter | Description |
|-----------|-------------|
| `iveg` | Method: 1=Suzuki et al. (2011), 2=Jacobsen et al. (2019) |
| `height` | Plant height per layer (m) |
| `diamtr` | Stem diameter per layer (m) |
| `nstems` | Stems per m^2 per layer |
| `drag` | Drag coefficient per layer |

Stem density can vary spatially using `INPGRID NPLANTS` / `READINP NPLANTS`.

### 4.12 Mud Dissipation

```
MUD [layer] [rhom] [viscm]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `layer` | Mud layer thickness (m) | required |
| `rhom` | Mud density (kg/m^3) | 1300 |
| `viscm` | Kinematic viscosity (m^2/s) | 0.0076 |

Mud thickness can vary spatially via `INPGRID MUDLAY` / `READINP MUDLAY`.

### 4.13 Sea Ice Dissipation

```
SICE { R19 [c0]-[c6] | D15 [Chf] | M18 [Chf] | R21B [Chf] [npf] }
```

- `R19`: Empirical polynomial (default)
- `D15`: Doble et al. (2015)
- `M18`: Meylan et al. (2018)
- `R21B`: Rogers et al. (2021)

Ice concentration and thickness input via `INPGRID AICE, HICE` / `READINP AICE, HICE`.

### 4.14 Default Physics Activation Summary

| Process | Activated by Default? | Command to Toggle |
|---------|----------------------|-------------------|
| Depth-induced breaking | YES | `OFF BREAKING` to disable |
| Wind input + whitecapping | With GEN1/2/3 | `GEN3 WESTH` (default) |
| Quadruplet interactions | With GEN3 | `OFF QUAD` to disable |
| Triad interactions | NO | `TRIAD` to enable |
| Bottom friction | NO | `FRICTION` to enable |
| Diffraction | NO | `DIFFRACTION` to enable |
| Vegetation | NO | `VEGETATION` to enable |
| Mud dissipation | NO | `MUD` to enable |
| Turbulence | NO | `TURBULENCE` to enable |
| Sea ice | NO | `SICE` to enable |
| Wave-induced setup | NO | `SETUP` to enable |

---

## 5. Numerical Methods

### 5.1 Implicit Scheme

SWAN employs fully implicit numerical schemes for propagation in geographic and spectral space. The key advantage: **unconditional stability**, meaning the time step and spatial resolution are not restricted by the CFL condition. This is critical for shallow water applications where the group velocity is high and explicit schemes would require prohibitively small time steps.

### 5.2 Propagation Schemes

**Geographic space (x,y) propagation:**

| Scheme | Order | Notes |
|--------|-------|-------|
| BSBT (backward space, backward time) | 1st | Simple, robust, diffusive |
| SORDUP (default for stationary) | 2nd | Less diffusive, better accuracy |
| S&L / Stelling & Leendertse (default for nonstationary) | 3rd | Best accuracy, least diffusion |

BSBT works with all grid types (regular, curvilinear, unstructured). Higher-order schemes are generally preferred.

**Directional space (theta):**

```
NUMERIC DIRIMPL [cdd]
```

`cdd` controls the scheme: 0 = central differences (accurate), 1 = upwind (diffusive). Default: 0.5 (blended).

### 5.3 Garden-Sprinkler Effect (GSE)

The GSE produces unphysical patterns in wave propagation due to discrete directional bins:

```
NUMERIC GSE [waveage]
```

`waveage` = 0: no diffusion (default). Positive values add directional diffusion.

### 5.4 Sweep Technique

SWAN uses a four-direction Gauss-Seidel iteration (sweep technique) in geographic space. In each iteration, four sweeps are performed, one for each quadrant of propagation directions.

### 5.5 Source Term Integration

Source terms are integrated using a semi-implicit or explicit scheme. The limiter restricts the rate of change per time step:

```
NUMERIC NONSTAT [mxitns] [limiter]
```

Default limiter = 0.1 (maximum fractional change per time step = 10%).

### 5.6 Convergence Criteria (Stationary Computations)

```
NUMERIC STOPC [dabs] [drel] [curvat] [npnts]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `dabs` | Absolute change in Hs between iterations (m) | 0.005 |
| `drel` | Relative change in Hs between iterations | 0.01 |
| `curvat` | Curvature of iteration curve | 0.005 |
| `npnts` | Percentage of grid points meeting criteria | 99.5% |

**Stationary iteration control:**

```
NUMERIC STAT [mxitst] [alfa]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `mxitst` | Maximum number of iterations | 50 |
| `alfa` | Under-relaxation factor | 0.00 |

Use `alfa = 0.01` when diffraction is activated.

### 5.7 Directional and Frequency Limiters

```
NUMERIC CTHETA [cfl]    $ directional CFL limiter, default: 0.9
NUMERIC CSIGMA [cfl]    $ frequency CFL limiter, default: 0.9
```

These limiters prevent excessive refraction or frequency shifting per iteration/time step.

---

## 6. Input File Structure

### 6.1 File Structure Overview

A SWAN input file (extension `.swn`) contains commands in sequential order. Comments start with `$` or `!`. The general order is:

```
1. PROJECT
2. SET (optional)
3. MODE (optional)
4. COORDINATES (optional)
5. CGRID
6. READGRID (if curvilinear/unstructured)
7. INPGRID + READINP (for each input field)
8. BOUNDSPEC / BOUNDNEST / INITIAL
9. Physics commands (GEN3, FRICTION, TRIAD, etc.)
10. NUMERIC (optional)
11. Output location definitions (FRAME, GROUP, CURVE, POINTS, etc.)
12. Output requests (BLOCK, TABLE, SPECOUT, NESTOUT)
13. COMPUTE
14. HOTFILE (optional)
15. STOP
```

### 6.2 PROJECT Command

```
PROJECT 'name' 'nr'
        'title1'
        'title2'
        'title3'
```

### 6.3 SET Command

```
SET [level] [nor] [depmin] [maxmes] [maxerr] [grav] [rho] [cdcap]
    [inrhog] [hsrerr] NAUTICAL|CARTESIAN [pwtail] [froudmax]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `level` | Water level offset (m) | 0 |
| `nor` | North direction relative to x-axis (degrees) | 90 |
| `depmin` | Minimum depth threshold (m) | 0.05 |
| `grav` | Gravitational acceleration (m/s^2) | 9.81 |
| `rho` | Water density (kg/m^3) | 1025 |
| `cdcap` | Maximum wind drag coefficient | 99999 |
| `inrhog` | Output units: 0=variance, 1=true energy | 0 |
| `pwtail` | High-frequency spectral tail power | 4-5 |

### 6.4 MODE & COORDINATES

```
MODE STATIONARY|NONSTATIONARY TWODIMENSIONAL|ONEDIMENSIONAL
COORDINATES CARTESIAN|SPHERICAL [CCM|QC]
```

Defaults: `STATIONARY TWODIMENSIONAL`, `CARTESIAN`. For large domains use `SPHERICAL`.

### 6.5 INPGRID & READINP Commands

```
INPGRID <BOTTOM|WLEV|CURRENT|WIND|FRICTION|...>
        <REGULAR [xpinp] [ypinp] [alpinp] [mxinp] [myinp] [dxinp] [dyinp]
        |CURVILINEAR [stagrx] [stagry]
        |UNSTRUCTURED>
        (NONSTATIONARY [tbeginp] [deltinp] SEC|MIN|HR|DAY [tendinp])

READINP <BOTTOM|WLEV|CURRENT|...>
        [fac] 'fname' [idla] [nhedf] [nhedt] [nhedvec]
        FREE|FORMAT 'form'|UNFORMATTED
```

**`idla` options (data ordering):**

| Value | Order |
|-------|-------|
| 1 | Left to right, top to bottom |
| 2 | Left to right, bottom to top |
| 3 | Right to left, top to bottom |
| 4 | Right to left, bottom to top |
| 5 | Top to bottom, left to right |
| 6 | Top to bottom, right to left |

**Constant wind (alternative to grid):**

```
WIND [vel] [dir]
```

### 6.6 BOUNDSPEC Command

**Boundary location:**

```
BOUNDSPEC SIDE NORTH|SOUTH|EAST|WEST CONSTANT|VARIABLE ...
BOUNDSPEC SEGMENT XY [x1] [y1] [x2] [y2] ... CONSTANT|VARIABLE ...
BOUNDSPEC DEFAULT ...
```

**Spectral parameters:**

```
... PAR [hs] [per] [dir] [dd]
... FILE 'fname' [seq]
```

**Spectral shape (use before BOUNDSPEC):**

```
BOUND SHAPESPEC JONSWAP|PM|GAUSS|BIN|TMA [gamma] PEAK|MEAN DSPR POWER|DEGREES [spread]
```

**TPAR files** (nonstationary parametric boundary conditions):

```
BOUNDSPEC ... FILE 'tpar_file.txt'
```

TPAR file format (one line per time step):
```
YYYYMMDD.HHMMSS  Hs  Per  Dir  Dd
```

### 6.7 INITIAL Conditions

```
INITIAL DEFAULT                    $ from local wind (Kahma & Calkoen)
INITIAL ZERO                      $ zero spectral energy everywhere
INITIAL PAR [hs] [per] [dir] [dd] $ uniform parametric spectrum
INITIAL HOTSTART SINGLE|MULTIPLE 'fname' FREE|UNFORMATTED
```

For MPI runs, `MULTIPLE` reads separate hotfiles per processor. `SINGLE` reads one concatenated file.

### 6.8 COMPUTE Command

**Stationary:**

```
COMPUTE STATIONARY [time]
```

**Nonstationary:**

```
COMPUTE NONSTATIONARY [tbegc] [deltc] SEC|MIN|HR|DAY [tendc]
```

Time step recommendation: at most 10 minutes for nonstationary runs.

### 6.9 HOTFILE & STOP

```
HOTFILE 'fname' FREE|UNFORMATTED
STOP
```

`HOTFILE` must appear immediately after a `COMPUTE` command. `STOP` is required and marks the end of all commands.

---

## 7. Output Types & Variables

### 7.1 Output Location Commands

**FRAME** (rectangular output grid):

```
FRAME 'sname' [xpfr] [ypfr] [alpfr] [xlenfr] [ylenfr] [mxfr] [myfr]
```

**POINTS** (individual output locations):

```
POINTS 'sname' [xp] [yp]         $ repeat for each point
POINTS 'sname' FILE 'fname'      $ read from file
```

**CURVE** (output along a line):

```
CURVE 'sname' [xp1] [yp1] [npts] [xp2] [yp2]
```

**NGRID** (nest output boundary):

```
NGRID 'sname' [xpn] [ypn] [alpn] [xlenn] [ylenn] [mxn] [myn]
```

### 7.2 Output Request Commands

**BLOCK** (gridded field output):

```
BLOCK 'sname' HEADER|NOHEADER 'fname' LAYOUT [idla] <var1> <var2> ...
      OUTPUT [tbeg] [delt] SEC|MIN|HR|DAY
```

Supported formats: ASCII, `.mat` (MATLAB), `.nc` (NetCDF), `.vtk` (ParaView).

**TABLE** (point/location output):

```
TABLE 'sname' HEADER|NOHEADER|INDEXED 'fname' <var1> <var2> ...
      OUTPUT [tbeg] [delt] SEC|MIN|HR|DAY
```

**SPECOUT** (spectral output):

```
SPECOUT 'sname' SPEC1D|SPEC2D ABS|REL 'fname'
        OUTPUT [tbeg] [delt] SEC|MIN|HR|DAY
```

**NESTOUT** (nest boundary spectra):

```
NESTOUT 'sname' 'fname' OUTPUT [tbeg] [delt] SEC|MIN|HR|DAY
```

### 7.3 Complete Output Variable Reference

| Abbreviation | Full Name | Units | Description |
|---|---|---|---|
| `HSIGN` | Significant wave height | m | Hs = 4*sqrt(integral of E) |
| `HSWELL` | Swell wave height | m | Hs below swell frequency |
| `TMM10` | Mean absolute period (m-1,0) | s | Inverse frequency-weighted |
| `TM01` | Mean absolute period (m0,1) | s | First spectral moment period |
| `TM02` | Mean absolute period (m0,2) | s | Second spectral moment period |
| `RTP` | Relative peak period | s | From discrete spectral bin |
| `TPS` | Smoothed peak period | s | Parabolic fit around peak |
| `DIR` | Mean wave direction | degrees | Normal to wave crests |
| `PDIR` | Peak direction | degrees | Direction at spectral peak |
| `TDIR` | Energy transport direction | degrees | Direction of energy flux |
| `DSPR` | Directional spreading | degrees | One-sided width |
| `QP` | Peakedness parameter | - | Spectrum narrowness |
| `DEPTH` | Water depth | m | Including setup and water level |
| `BOTLEV` | Bottom level | m | Raw bottom elevation |
| `WATLEV` | Water level | m | Imposed water level change |
| `VEL` | Current velocity | m/s | x and y components |
| `WIND` | Wind velocity | m/s | x and y components |
| `FORCE` | Wave-induced force | N/m^2 | Radiation stress gradient |
| `UBOT` | Bottom orbital velocity | m/s | RMS value near bottom |
| `STEEPNESS` | Wave steepness | - | Hs / wavelength |
| `WLEN` | Mean wavelength | m | Spectral mean |
| `BFI` | Benjamin-Feir Index | - | Freak wave indicator |
| `QB` | Fraction of breaking waves | - | Battjes-Janssen fraction |
| `TRANSP` | Energy transport | W/m | Wave energy flux |
| `SETUP` | Wave-induced setup | m | Mean water level rise |
| `DISSIP` | Total dissipation | W/m^2 | All dissipation combined |
| `GENERAT` | Wind generation | W/m^2 | Wind input source term |
| `AICE` | Ice concentration | - | Fraction 0-1 |
| `TIME` | Date-time | timestamp | Full date/time string |

### 7.4 QUANTITY Command (Output Configuration)

```
QUANTITY 'output_var' 'short' 'long' [lexp] [hexp] [excv] [fswell] [fmin] [fmax]
         PROBLEMCOORD|FRAME
```

Controls naming, precision, exception values, and frequency ranges for output.

---

## 8. Nesting

### 8.1 One-Way Nesting (Coarse to Fine)

**Step 1: Coarse grid run**

Define a nest output boundary in the coarse-grid input file:

```
NGRID 'nestgrid' [xpn] [ypn] [alpn] [xlenn] [ylenn] [mxn] [myn]
NESTOUT 'nestgrid' 'nest_spectra.swn' OUTPUT [tbeg] [delt] HR
```

**Step 2: Fine grid run**

Read the nest boundary conditions:

```
BOUNDNEST1 NEST 'nest_spectra.swn' CLOSED
```

The fine grid automatically interpolates the coarse-grid spectra to its own frequency and directional resolution.

### 8.2 Cross-Model Nesting

**WAM to SWAN:**

```
BOUNDNEST2 WAM 'wam_output.dat' UNFORMATTED|FREE [xgc] [ygc] [lwdate]
```

**WAVEWATCH III to SWAN:**

```
BOUNDNEST3 WW3 'ww3_output.dat' UNFORMATTED|FREE CLOSED|OPEN [xgc] [ygc]
```

Cross-model boundary conditions may not be fully model-consistent due to different numerical implementations.

### 8.3 Guidelines for Nesting

- Use the same coordinate system (Cartesian or spherical) for both grids
- Place the coarse grid's deep-water boundary where shallow-water effects do not dominate
- Spatial and spectral resolution ratios should be within a factor of 2-3
- For larger resolution jumps, use intermediate nesting levels
- Curvilinear grids can be used, but boundaries must remain rectangular

### 8.4 Multiple Nested Grids

Multiple nesting levels can be chained:

```
Global (WW3/WAM) --> Regional (SWAN coarse) --> Local (SWAN fine) --> Harbor (SWAN finest)
```

Each level generates `NESTOUT` for the next level's `BOUNDNEST1` input. There is no limit on nesting depth, but each level introduces interpolation artifacts.

---

## 9. Unstructured Grids

### 9.1 Triangle-Based Mesh

SWAN supports triangular unstructured meshes with the vertex-based, fully implicit finite difference method.

**Input files required:**
- Node file: coordinates and depth values for each vertex
- Element file: triangle connectivity (3 vertices per element)
- Boundary marker file: identifies boundary vertices

### 9.2 Mesh Generation

| Tool | Description | License |
|------|-------------|---------|
| **Triangle** | Free Delaunay triangulator via `.poly` files | Free |
| **Gmsh** | Open-source FEM mesh generator with GUI | GPL |
| **BatTri** | MATLAB-based mesh generator for coastal apps | Free |
| **SMS** | Commercial mesh generation by Aquaveo | Commercial |

### 9.3 Mesh Quality Considerations

- Avoid highly elongated triangles (aspect ratio > 5:1)
- Gradual transition in element size (ratio of adjacent elements < 1.5:1)
- Resolution should be finest where bathymetric gradients are steepest
- Typical range: 10-50 m nearshore, 1-5 km offshore
- Total node count is usually much less than equivalent structured grid

### 9.4 Performance Considerations

- Per-node computation is slightly more expensive than structured grids due to indirect addressing
- Total computation time is typically less because unstructured grids need fewer nodes
- METIS library recommended for domain decomposition in MPI mode
- Unconditional stability is maintained regardless of time step

---

## 10. Pre/Post-Processing & Visualization

### 10.1 Grid Generation Tools

| Tool | Grid Type | Free/Commercial | Platform |
|------|-----------|-----------------|----------|
| RGFGRID | Curvilinear | Free (Deltares) | Windows/Linux |
| SMS | All types | Commercial | Windows |
| Gmsh | Unstructured | Free (GPL) | All |
| Triangle | Unstructured | Free | All |
| BatTri | Unstructured | Free (MATLAB) | MATLAB |
| Delft Dashboard | Curvilinear | Free | MATLAB |

### 10.2 Boundary Condition Preparation

**From buoy data:**
- Convert measured spectra to SWAN `SPEC2D` format
- Or use parametric input (`PAR`) with Hs, Tp, Dir, Dd

**From global models:**
- ERA5/ECMWF: Extract spectra at boundary points, convert to `SPEC2D`
- WAVEWATCH III: Use `NESTOUT`/`BOUNDNEST3` workflow
- WAM: Use `BOUNDNEST2`

**TPAR files (nonstationary parametric):**

Format -- one line per time step:
```
20240101.000000  2.5  10.0  270.0  30.0
```

Fields: `date.time  Hs(m)  Period(s)  Direction(deg)  Spread(deg)`

### 10.3 Python Tools

**swantools** (`pip install swantools`):
- Read TABLE output to pandas DataFrame
- Read BLOCK output
- Read SPEC2D output
- Write TPAR and SPEC files
- GitHub: https://github.com/caiostringari/swantools

**MHKiT:**
- Read SWAN block and table outputs (ASCII and `.mat`)
- Returns pandas DataFrames

**Custom script patterns:**

Reading TABLE output:
```python
import pandas as pd
df = pd.read_csv('output.tab', comment='%', delim_whitespace=True,
                 header=None, names=['Xp', 'Yp', 'Hsig', 'Dir', 'Tm01'])
```

Reading BLOCK output:
```python
import numpy as np
data = np.loadtxt('output.blk', comments='%')
field = data.reshape((ny+1, nx+1))
```

### 10.4 Visualization

**Python (matplotlib):**

```python
import matplotlib.pyplot as plt
import numpy as np

hs = np.loadtxt('hsig.blk', comments='%').reshape((ny+1, nx+1))
plt.contourf(x, y, hs, levels=20, cmap='jet')
plt.colorbar(label='Hs (m)')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Significant Wave Height')
plt.show()
```

**MATLAB:** Native support for `.mat` block output files with `load('output.mat')`.

**ParaView:** Use `.vtk` output format from BLOCK command for 3D visualization.

**GIS integration:** Export block output to ASCII grid format, import into QGIS or ArcGIS.

---

## 11. Running SWAN

### 11.1 Serial Mode

```bash
./swanrun -input mycase
```

Or manually:

```bash
swan.exe < INPUT
```

The input file must be named `INPUT` (or specified via `swanrun`). Output goes to `PRINT` file and user-specified output files.

### 11.2 OpenMP Mode

```bash
./swanrun -input mycase -omp 4
```

Or set the environment variable:

```bash
export OMP_NUM_THREADS=4
./swan.exe
```

- Good efficiency for up to ~8 threads
- Diminishing returns beyond 8 threads on most architectures
- No domain decomposition; shared memory parallelism within sweeps

### 11.3 MPI Mode

```bash
./swanrun -input mycase -mpi 16
```

Or manually:

```bash
mpirun -np 16 swan.exe
```

- Better scaling than OpenMP for large core counts
- Domain decomposition across processors
- Each processor writes its own output; results automatically concatenated
- For unstructured grids, METIS provides optimal decomposition

**MPI with job scheduler (SLURM example):**

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00

module load swan/41.51
srun swan.exe
```

### 11.4 Hybrid Mode (OpenMP + MPI)

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 swan.exe
```

This uses 8 MPI processes, each with 4 OpenMP threads = 32 total cores.

### 11.5 Hot Start / Restart

**Writing a hotfile:**

```
COMPUTE NONSTATIONARY 20240101.000000 10 MIN 20240102.000000
HOTFILE 'hotfile.swn' FREE
```

**Restarting from a hotfile:**

```
INITIAL HOTSTART 'hotfile.swn' FREE
```

**MPI restart:**
- `MULTIPLE`: Reads separate hotfiles per processor -- requires same processor count as original run
- `SINGLE`: Reads one concatenated hotfile -- more flexible; allows changing processor count

### 11.6 Performance Tips

1. **Spectral resolution:** Reduce `mdc` and `msc` to minimum acceptable values. Each doubling roughly doubles computation time.
2. **Time step:** Use the largest stable time step. 5-10 minutes is typical for nonstationary runs.
3. **Grid resolution:** Use unstructured grids to concentrate resolution where needed rather than uniformly fine structured grids.
4. **Convergence criteria:** For stationary runs, relax `STOPC` criteria if high precision is not needed.
5. **Number of iterations:** For nonstationary runs, 1 iteration per time step (`mxitns=1`) is usually sufficient.
6. **Output frequency:** Reduce output frequency for large nonstationary runs to limit file sizes.
7. **Physics selection:** Only activate physics processes relevant to your domain. Triads are expensive.
8. **Propagation scheme:** BSBT is fastest but most diffusive. Default schemes provide good balance.
9. **MPI scaling:** Typical good scaling up to ~100-200 processors for large domains.
10. **Hot starting:** Break long simulations into segments with hotfiles for checkpointing.

### 11.7 Output Files

| File | Description |
|------|-------------|
| `PRINT` | Main log file: input echo, iteration info, errors, timings |
| `Errfile` | Error messages (only created if errors occur) |
| `ERRPTS` | Grid points with computation errors |
| User-named files | BLOCK, TABLE, SPECOUT, NESTOUT output |
| Hotfiles | Restart files (one per processor in MPI mode) |

### 11.8 Example Input File

```
$ SWAN Example: Stationary nearshore simulation
$
PROJECT 'example' '001'
'Nearshore wave simulation'
'Sandy beach, wind from west'

SET level=0.0 nor=90 depmin=0.05
MODE STATIONARY TWODIMENSIONAL
COORDINATES CARTESIAN

$ Computational grid
CGRID REGULAR 0.0 0.0 0.0 10000.0 5000.0 100 50 &
      CIRCLE 36 0.03 1.0 31

$ Bottom grid and data
INPGRID BOTTOM REGULAR 0.0 0.0 0.0 100 50 100.0 100.0
READINP BOTTOM 1.0 'bottom.bot' 1 0 FREE

$ Wind forcing
WIND 15.0 270.0

$ Boundary conditions
BOUND SHAPESPEC JONSWAP 3.3 PEAK DSPR DEGREES
BOUNDSPEC SIDE WEST CONSTANT PAR 2.0 10.0 270.0 30.0

$ Physics
GEN3 WESTH
FRICTION JONSWAP 0.038
BREAKING CONSTANT 1.0 0.73
TRIAD

$ Numerics
NUMERIC STOPC 0.005 0.01 0.005 99.5 STAT 50

$ Output locations
FRAME 'field' 0.0 0.0 0.0 10000.0 5000.0 100 50
POINTS 'buoys' 5000.0 2500.0 &
                8000.0 2500.0 &
                9500.0 2500.0

$ Output requests
BLOCK 'field' NOHEADER 'hsig.blk' LAYOUT 1 HSIGN
BLOCK 'field' NOHEADER 'dir.blk' LAYOUT 1 DIR
BLOCK 'field' NOHEADER 'tm01.blk' LAYOUT 1 TM01
TABLE 'buoys' HEADER 'buoys.tab' HSIGN DIR TM01 DSPR UBOT
SPECOUT 'buoys' SPEC2D ABS 'buoys.sp2'

$ Compute
COMPUTE STATIONARY

STOP
```

---

## 12. Calibration & Validation

### 12.1 Key Calibration Parameters

**Primary calibration knobs (in order of typical importance):**

| Parameter | Command | Range | Impact |
|-----------|---------|-------|--------|
| Whitecapping coefficient | `GEN3 KOMEN cds2` / `JANSSEN cds1` | 2e-5 to 5e-5 / 3.0 to 6.0 | Overall wave height, especially deep water |
| Breaking index (gamma) | `BREAKING CONSTANT [alpha] [gamma]` | 0.6 - 0.83 | Wave height in surf zone |
| Bottom friction (cfjon) | `FRICTION JONSWAP [cfjon]` | 0.019 - 0.067 | Wave height decay over shallow shelves |
| Triad interaction factor | `TRIAD [method] [trfac]` | 0.5 - 5.0 | Nearshore spectral shape |
| Wind drag | `SET cdcap` | 0.002 - 0.003 | Extreme wind conditions |

### 12.2 Calibration Strategies

**Deep water (wave generation):**
- Whitecapping coefficient is the primary tuning parameter
- ST6 or WESTH generally perform better than KOMEN out-of-the-box
- Calibrate whitecapping first; under high-quality wind forcing, optimal coefficient falls in a narrow range

**Intermediate depth (wave propagation):**
- Bottom friction controls energy decay over continental shelves
- `cfjon = 0.038` for sandy bottoms (default)
- `cfjon = 0.019` for smooth seafloors (e.g., Gulf of Mexico)
- `cfjon = 0.067` for rough bottoms or depth-limited conditions
- MADSEN formulation with spatially varying roughness for heterogeneous seabeds

**Nearshore (wave transformation):**
- Breaking index `gamma` dominates the surf zone (range 0.6-0.83, slope-dependent)
- BKD method for slope-dependent breaking
- Triad interactions affect spectral shape and energy redistribution

**Multi-parameter approach (calibrate from deep to shallow):**
1. Whitecapping first (deep water)
2. Bottom friction (intermediate depth, 5-10 m)
3. Breaking (inside the surf zone)

### 12.3 Validation Metrics

| Metric | Formula | Target |
|--------|---------|--------|
| Bias | mean(model - obs) | < 0.1 m for Hs |
| RMSE | sqrt(mean((model-obs)^2)) | < 15% of mean Hs |
| Scatter Index | RMSE / mean(obs) | < 0.15 |
| Correlation (R) | Pearson correlation | > 0.9 |

### 12.4 Common Pitfalls

- Wind field quality is often the largest source of error; calibrate wind first
- Bathymetry resolution and accuracy are critical for nearshore results
- Numerical diffusion (BSBT scheme) can artificially smooth wave fields
- Insufficient spectral resolution (too few frequencies/directions) degrades results
- Ignoring triads in shallow water can underestimate infragravity energy
- Default physics settings changed significantly between versions (especially v41.45); verify which defaults are active

---

## References & Resources

### Official Documentation
- User Manual: https://swanmodel.sourceforge.io/download/zip/swanuse.pdf
- Technical Documentation: https://swanmodel.sourceforge.io/download/zip/swantech.pdf
- Implementation Manual: https://swanmodel.sourceforge.io/download/zip/swanimp.pdf
- Online Documentation: https://swanmodel.sourceforge.io/online_doc/
- Release Notes: https://swanmodel.sourceforge.io/modifications/modifications.htm

### Source Code
- SourceForge: https://sourceforge.net/projects/swanmodel/files/swan/
- TU Delft GitLab: https://gitlab.tudelft.nl/citg/wavemodels/swan

### Python Tools
- swantools: https://github.com/caiostringari/swantools (`pip install swantools`)
- python-swan: https://pypi.org/project/python-swan/
- MHKiT: https://mhkit-software.github.io/MHKiT/

### Mesh Generation
- Triangle: https://www.cs.cmu.edu/~quake/triangle.html
- Gmsh: https://gmsh.info/
- SMS: https://www.aquaveo.com/software/sms-surface-water-modeling-system
