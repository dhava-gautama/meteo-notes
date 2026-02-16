# WRF-Chem Tutorial Guide: From Dust to Global Emissions

> Practical tutorial notes based on the [WRF-Chem Tutorial Exercises](https://hathewaywill.github.io/wrf-chem-site/tutorialexercises.htm) (Hatheway).
> Covers WRF-Chem v3.9/4.x — dust simulation, volcanic ash, US anthropogenic emissions (NEI), and global multi-source emissions with GOCART.

---

## Table of Contents

1. [WRF-Chem Workflow Overview](#1-wrf-chem-workflow-overview)
2. [Exercise 1: Dust Simulation](#2-exercise-1-dust-simulation)
3. [Exercise 2: Volcanic Ash Simulation](#3-exercise-2-volcanic-ash-simulation)
4. [Exercise 3: US Anthropogenic Emissions (NEI)](#4-exercise-3-us-anthropogenic-emissions-nei)
5. [Exercise 4: Global Emissions (RETRO/EDGAR/GOCART)](#5-exercise-4-global-emissions-retroedgargocart)
6. [Auxiliary Input Channel Reference](#6-auxiliary-input-channel-reference)
7. [Key Chemistry Options Quick Reference](#7-key-chemistry-options-quick-reference)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. WRF-Chem Workflow Overview

WRF-Chem adds an emissions preparation and conversion stage between WPS and the WRF run. The general flow is:

```
[WPS: geogrid → ungrib → metgrid]
        |
        v
[Emissions Preparation]
  ├── Dust:           geogrid GEOGRID.TBL.ARW_CHEM → erodibility fields
  ├── Volcanic ash:   prep_chem_sources → convert_emiss → wrfchemv_d01
  ├── Anthropogenic:  NEI processor / prep_chem_sources → convert_emiss → wrfchemi_*
  ├── Biomass burn:   prep_chem_sources → convert_emiss → wrffirechemi_*
  └── GOCART bkgd:    prep_chem_sources → convert_emiss → wrfchemi_gocart_bg_d01
        |
        v
[real.exe]  →  wrfinput_d01 + wrfbdy_d01  (now includes chem arrays)
        |
        v
[wrf.exe]   →  wrfout_d0* with chemistry species (DUST_*, vash_*, SO2, PM, etc.)
        |
        v
[Analysis]  →  ncview, NCL, ncdiff for comparing experiments
```

**Compilation** (applies to all exercises):
```bash
export WRF_CHEM=1
export WRF_KPP=1          # optional, needs FLEX
./configure && ./compile em_real
```

---

## 2. Exercise 1: Dust Simulation

**Objective:** Simulate Saharan dust transport over the Mediterranean using GOCART dust scheme.

### 2.1 Domain Setup

| Parameter | Value |
|---|---|
| Projection | Lambert Conformal |
| Grid dimensions | 41 × 41 |
| Grid spacing | 100 km |
| Center | 35°N, 25°E |
| Forcing data | GFS 6-hourly analysis |
| Period | 14–19 July 2010 |

### 2.2 Geogrid with Chemistry Tables

Standard WRF geogrid uses `GEOGRID.TBL.ARW`. For chemistry with dust, you must switch to the chemistry table — it includes the erodibility map and soil texture fields:

```bash
cd ${WPS}/geogrid
ln -svf GEOGRID.TBL.ARW_CHEM GEOGRID.TBL
cd ..
```

Also download the **erodibility dataset** from the WPS static geography archive:
- Source: http://www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html

Then run geogrid normally:
```bash
./geogrid.exe
```

**Verify** with `ncview geo_em.d01.nc` — check for:
- `HGT_M` — terrain height
- `LANDMASK` — land/sea mask
- `EROD` — erodibility fraction (must have non-zero values over desert)
- `CLAYFRAC`, `SANDFRAC` — soil texture fractions

### 2.3 Metgrid and Real

```bash
./metgrid.exe

cd ${WRF_CHEM}
ln -svf ${WPS}/met_em.* .
mpirun -np 8 ./real.exe
```

Verify `wrfinput_d01` contains `DUST_1` through `DUST_5`, `EROD`, `CLAYFRAC` fields.

### 2.4 Key Namelist Settings

```
&chem
 chem_opt         = 401,       ! GOCART simple (bulk aerosol)
 dust_opt         = 1,         ! GOCART dust (AFWA scheme)
/
```

### 2.5 Run and Analyze

```bash
mpirun -np 8 ./wrf.exe
```

**Post-processing — compute total dust loading:**
```bash
# Sum all 5 dust size bins into a single field
ncap -v -s "total_dust=DUST_1+DUST_2+DUST_3+DUST_4+DUST_5" \
  wrfout_d01_2010-07-14_00:00:00 dust_total.nc
```

### 2.6 Sensitivity: Comparing Dust Schemes

The tutorial compares different dust emission parameterizations:

| `dust_opt` | Scheme | Description |
|---|---|---|
| 1 | **GOCART (AFWA)** | Original GOCART dust with AFWA erodibility |
| 3 | **GOCART-Thompson** | Modified threshold velocities |
| 4 | **GOCART-UoC** | University of Cologne formulation |

**Differencing two experiments:**
```bash
# Run with dust_opt=3, then difference against dust_opt=1
ncap -v -s "total_dust=DUST_1+DUST_2+DUST_3+DUST_4+DUST_5" \
  wrfout_d01_2010-07-14_00:00:00 dust3_all.nc

ncdiff dust3_all.nc dust1_all.nc diff_dust3_1.nc
ncview diff_dust3_1.nc    # visualize differences
```

This reveals sensitivity of dust emissions to the parameterization choice — useful for calibration studies.

---

## 3. Exercise 2: Volcanic Ash Simulation

**Objective:** Simulate volcanic ash dispersal from an eruption using prep_chem_sources and WRF-Chem.

### 3.1 General Approach (4-Step Method)

1. Run WRF **without** chemistry to produce meteorological fields (`wrfinput_d01`)
2. Run `prep_chem_sources` to generate emissions in intermediate binary format
3. Run `convert_emiss.exe` to convert intermediates into WRF-Chem NetCDF inputs
4. Re-run WRF **with** chemistry enabled

### 3.2 Prep-Chem-Sources Configuration

Edit `prep_chem_sources.inp`:

```fortran
ihour               = 00
iday                = 14
imon                = 07
iyear               = 2010
use_volcanoes       = 1
volcano_index       = 9          ! Mount Vesuvius, Italy
begin_eruption      = '201007140000'   ! UTC format: YYYYMMDDhhmm
```

**Volcano selection:** The file `volc_emissions.f90` contains a table of **1,535 volcanoes** worldwide. Look up the volcano you need and set `volcano_index` to its integer index.

Run the preprocessor:
```bash
./prep_chem_sources_RADM_WRF_FIM.exe
```

**Expected output files:**

| File | Size | Content |
|---|---|---|
| `*-g1-volc.bin` | ~34 KB | Volcanic ash emissions |
| `*-g1-gocartBG.bin` | ~1.1 MB | GOCART background (DMS, EROD, H2O2, OH, NO3) |
| `*-g1-bb.bin` | ~236 KB | Biomass burning |
| `*-g1-ab.bin` | ~182 KB | Anthropogenic background |

### 3.3 Convert Emissions

Create symbolic links to the binary output:
```bash
cd ${Run_convert_emiss}
ln -sf ../Run_prep/*-g1-gocartBG.bin    wrf_gocart_backg
ln -sf ../Run_prep/*-g1-ab.bin          emissopt3_d01
ln -sf ../Run_prep/*-g1-bb.bin          emissfire_d01
ln -sf ../Run_prep/*-g1-volc.bin        volc_d01
```

**Note:** The link names (`wrf_gocart_backg`, `emissopt3_d01`, `emissfire_d01`, `volc_d01`) are **hardcoded** — `convert_emiss.exe` looks for exactly these filenames.

Copy `wrfinput_d01` (from the no-chemistry run) and `namelist.input` into the convert directory.

**Namelist additions for volcanic ash:**
```
io_form_auxinput13 = 2       ! NetCDF format for volcanic ash channel
emiss_opt_vol      = 1       ! Enable volcanic emissions
```

Run conversion:
```bash
./convert_emiss.exe
```

Output: `wrfchemv_d01` — volcanic ash emissions in NetCDF. Verify with `ncview wrfchemv_d01` and check the `E_VASH9` field.

### 3.4 WRF-Chem Run with Volcanic Ash

```
&chem
 chem_opt          = 400,      ! GOCART full (needed for volcanic ash)
 emiss_opt_vol     = 1,        ! Volcanic emissions enabled
/
```

```bash
mpirun -np 8 ./real.exe       # wrfinput now has vash_1..vash_10 (initially zero)
mpirun -np 8 ./wrf.exe        # ~5 min runtime
```

**Verify runtime:** Look in `rsl.out.0000` for:
```
mediation_integrate: med_read_wrf_volc_emiss: Open file wrfchemv_d01
mediation_integrate: med_read_wrf_volc_emiss: Read volcanic ash emissions
```

Output contains 10 volcanic ash size bins: `vash_1` through `vash_10`.

### 3.5 Plume Height Sensitivity

Control the eruption injection height via namelist:
```
emiss_ash_hgt = 20000.      ! Plume height in meters (default varies)
```

Higher plume heights inject ash into upper troposphere/stratosphere — compare experiments with different values to study dispersal patterns. Note that deeper eruptions produce substantially larger emissions, so adjust color scales when comparing.

---

## 4. Exercise 3: US Anthropogenic Emissions (NEI)

**Objective:** Process EPA National Emissions Inventory (NEI-2011) for a CONUS domain using the dedicated Fortran emission processor.

### 4.1 Domain Configuration

| Parameter | Value |
|---|---|
| Projection | Polar Stereographic |
| Grid dimensions | 81 × 51 |
| Grid spacing | 60 km |
| Center | 45°N, 97°W |
| Standard parallels | 45°N, 33°N |

### 4.2 Required Data

- **NEI-2011 inventory:** `em11v1_file1.tar`, `em11v1_file2.tar` (point + area sources)
- **Fortran processor:** `emiss_v04_CONUS60km.F`
- **Pre-processed binaries:** `NEI11_processed/` directory (saves reprocessing time)

### 4.3 Compilation

```bash
ifort -FR -convert big_endian emiss_v04_CONUS60km.F
```

Update paths in the Fortran source before compiling:
```fortran
DATA POINDIR /{your_path}/NEI_input/point/
DATA AREADIR /{your_path}/NEI_input/
```

### 4.4 Convert Emissions (Two 12-Hour Periods)

NEI emissions are split into two 12-hour files. Each period requires a separate `convert_emiss.exe` run.

**First period (00Z–12Z):**
```bash
ln -sf ${path}/NEI11_processed/wrfem_00to12Z wrfem_00to12z_d01
ln -sf ${path}/NEI11_processed/wrfem_12to24Z wrfem_12to24z_d01
ln -sf ${path}/Input_files/wrfinput_d01_00Z  wrfinput_d01
```

Namelist:
```
start_hour  = 00
end_hour    = 11
run_hours   = 12
emiss_opt   = 3         ! NEI style
kemit       = 20        ! 20 vertical emission levels
```

```bash
mpirun -np 8 ./convert_emiss.exe
```

**Second period (12Z–23Z):**
```bash
ln -sf ${path}/Input_files/wrfinput_d01_12Z wrfinput_d01
# Update namelist: start_hour=12, end_hour=23
mpirun -np 8 ./convert_emiss.exe
```

**Verify:**
```bash
ncview wrfchemi_00z_d01
ncdump -v Times wrfchemi_*    # check time stamps
```

### 4.5 Things to Examine

- **Spatial patterns:** Urban emission hotspots, power plant plumes
- **Diurnal variation:** Compare 00Z vs 12Z fields — traffic-related emissions peak during daytime hours
- **Species breakdown:** SO2, NOx, CO, VOC — each has distinct spatial signatures

---

## 5. Exercise 4: Global Emissions (RETRO/EDGAR/GOCART)

**Objective:** Process global anthropogenic, sea salt, and biomass burning emissions using prep_chem_sources with the RACM-KPP-GOCART mechanism.

### 5.1 Prep-Chem-Sources Configuration

```fortran
grid_type           = 'lambert'
ihour=00, iday=14, imon=07, iyear=2010

! Domain (same as Exercise 1)
NGRIDS    = 1
NNXP      = 41,   NNYP    = 41
DELTAX    = 100000.,  DELTAY = 100000.    ! 100 km
POLELAT   = 35.,   POLELON = 25.
STDLAT1   = 30.,   STDLAT2 = 40.

! Emission inventories
use_retro           = 1        ! RETRO global anthropogenic VOC
use_edgar           = 3        ! EDGAR v4.2 global anthropogenic
use_gocart          = 1        ! GOCART aerosol module
use_bioge           = 1        ! Biogenic emissions
use_bbem            = 1        ! Biomass burning emissions
use_bbem_plumerise  = 1        ! Plume rise for biomass burning
use_gocart_bg       = 1        ! GOCART background fields (OH, H2O2, NO3, DMS)
use_volcanoes       = 0        ! Disabled for this exercise
```

### 5.2 Run Prep-Chem-Sources

```bash
./prep_chem_sources_RADM_WRF_FIM_.exe < prep_chem_sources.inp
```

**Output files:**

| File | Size | Content |
|---|---|---|
| `*-g1-ab.bin` | ~181 KB | Anthropogenic (RETRO + EDGAR) |
| `*-g1-bb.bin` | ~235 KB | Biomass burning |
| `*-g1-gocartBG.bin` | ~1.1 MB | GOCART background (DMS, EROD, H2O2, OH, NO3) |

### 5.3 Link and Convert

```bash
# Link intermediate binaries with required names
ln -s *-g1-ab.bin          emissopt3_d01
ln -s *-g1-bb.bin          emissfire_d01
ln -s *-g1-gocartBG.bin    wrf_gocart_backg
```

### 5.4 Key Namelist Settings

```
&chem
 chem_opt             = 301,    ! RACM-KPP gas + GOCART aerosol
 emiss_opt            = 5,      ! RACM-KPP GOCART emissions
 emiss_inpt_opt       = 1,
 kemit                = 1,      ! Global inventories = single surface level

 ! Aerosol processes
 biomass_burn_opt     = 1,
 dust_opt             = 1,
 dmsemis_opt          = 1,      ! DMS sea-air emissions
 seas_opt             = 1,      ! Sea salt emissions

 ! Radiative feedback (start OFF for baseline)
 aer_ra_feedback      = 0,
 aer_op_opt           = 0,
 opt_pars_out         = 0,
/
```

**Emissions I/O:**
```
io_style_emissions      = 2         ! Date/time specific files
auxinput5_interval_m    = 100000    ! Anthropogenic (static — read once)
auxinput7_interval_m    = 100000    ! Biomass burning (static)
auxinput8_interval_m    = 100000    ! GOCART background (static)
io_form_auxinput5       = 2         ! NetCDF
io_form_auxinput7       = 2
io_form_auxinput8       = 2
```

Setting `auxinput*_interval_m` to a very large number (100000 min ≈ 69 days) means the emissions file is read **once** and held constant throughout the simulation. This is appropriate for global inventories that don't have sub-daily temporal resolution.

### 5.5 Run Convert, Real, and WRF

```bash
# Convert emissions → NetCDF
mpirun -np 8 ./convert_emiss.exe

# Output files:
#   wrfchemi_d01             → anthropogenic
#   wrffirechemi_d01         → biomass burning
#   wrfchemi_gocart_bg_d01   → GOCART background

# Create time-stamped link for wrfchemi
ln -s wrfchemi_d01 wrfchemi_d01_2010-07-14_00:00:00

# Run real and wrf
mpirun -np 4 ./real.exe
mpirun -np 8 ./wrf.exe
```

**Verify `wrfinput_d01`** after real.exe:
- `backg_oh` — OH background (from GOCART)
- `ebu_in_co` — biomass burning CO
- These fields should have non-zero values if conversion was successful

### 5.6 Output Analysis

Check output for:
- **Sea salt aerosols** — should appear over ocean and correlate with PM10 coarse
- **SO2** — industrial/urban signatures
- **PM, CO from fires** — elevated above surface in fire-active regions

### 5.7 Aerosol Direct Radiative Feedback Experiment

After a baseline run (feedback OFF), enable aerosol-radiation interaction:

```
aer_ra_feedback = 1      ! Aerosol direct effect on radiation
aer_op_opt      = 1      ! Compute aerosol optical properties
opt_pars_out    = 1      ! Output extinction coefficients
```

Rerun `wrf.exe` and compare:

```bash
# Difference key fields between feedback-ON and feedback-OFF runs
ncdiff -v SWDOWN,RAINC,T2 wrfout_feedback_on wrfout_feedback_off diff.nc
ncview diff.nc
```

**What to look for:**
- `SWDOWN` decrease where aerosol loading is high (dust, urban plumes)
- `T2` cooling under heavy aerosol — reduced surface shortwave
- `EXTCOF55` (extinction coefficient at 550 nm) — spatial correlation with dust and anthropogenic aerosol fields
- `RAINC` — changes in convective precipitation from modified thermodynamic profiles

---

## 6. Auxiliary Input Channel Reference

WRF-Chem uses numbered auxiliary input channels for different emission types. These are mapped in the WRF-Chem registry (column 8 of `registry.chem`):

| Channel | Purpose | Typical File |
|---|---|---|
| `auxinput5` | Anthropogenic emissions | `wrfchemi_d01` / `wrfchemi_00z_d01` |
| `auxinput6` | Biogenic emissions (BEIS, MEGAN) | `wrfbiochemi_d01` |
| `auxinput7` | Surface biomass burning | `wrffirechemi_d01` |
| `auxinput8` | GOCART background fields | `wrfchemi_gocart_bg_d01` |
| `auxinput12` | Chemistry initial conditions | from mozbc / global model |
| `auxinput13` | Volcanic ash emissions | `wrfchemv_d01` |
| `auxinput14` | Aircraft emissions | — |
| `auxinput15` | Greenhouse gas emissions | — |

Each channel needs corresponding `io_form_auxinputN = 2` (NetCDF) and `auxinputN_interval_m` in `namelist.input`.

---

## 7. Key Chemistry Options Quick Reference

### `chem_opt` (Chemistry Mechanism)

| Value | Gas-Phase | Aerosol | Use Case |
|---|---|---|---|
| 0 | None | None | Meteorology only |
| 1 | RADM2 | None | Gas-phase only |
| 2 | RADM2 | MADE/SORGAM | Modal aerosol |
| 8 | CBMZ | MOSAIC 8-bin | Sectional aerosol |
| 112 | MOZART | GOCART bulk | MOZCART |
| 202 | MOZART | MOSAIC 4-bin + VBS SOA | Full treatment |
| 301 | **RACM-KPP** | **GOCART** | Tutorial Ex. 4 |
| 400 | — | **GOCART full** | Volcanic ash (Ex. 2) |
| 401 | — | **GOCART simple** | Dust only (Ex. 1) |

### `dust_opt` (Dust Emission Scheme)

| Value | Scheme | Notes |
|---|---|---|
| 1 | GOCART-AFWA | Original, uses erodibility map |
| 3 | GOCART-Thompson | Modified threshold friction velocities |
| 4 | GOCART-UoC | University of Cologne formulation |

### `emiss_opt` (Emission Input Style)

| Value | Style | Notes |
|---|---|---|
| 3 | NEI (US national inventory) | 12-hourly split, kemit=20 vertical levels |
| 5 | RACM-KPP-GOCART global | Single surface level (kemit=1) |

### Other Key Switches

| Namelist Variable | Values | Purpose |
|---|---|---|
| `emiss_opt_vol` | 0/1 | Volcanic ash emissions on/off |
| `emiss_ash_hgt` | meters | Eruption plume height |
| `biomass_burn_opt` | 0/1 | Biomass burning emissions |
| `dmsemis_opt` | 0/1 | DMS ocean-atmosphere emissions |
| `seas_opt` | 0/1 | Sea salt emissions |
| `aer_ra_feedback` | 0/1 | Aerosol direct radiative effect |
| `aer_op_opt` | 0/1 | Compute aerosol optical properties |
| `opt_pars_out` | 0/1 | Output extinction coefficients |
| `kemit` | int | Vertical levels in emission files |

---

## 8. Troubleshooting

### `convert_emiss.exe` Read Errors

**Symptom:** `PGFIO-F-219/unformatted read/write past end of record` or `forrtl: severe (67): input statement requires too much data`

**Causes and fixes:**
- Wrong symbolic link names — `convert_emiss.exe` expects **exact** filenames (`wrf_gocart_backg`, `emissopt3_d01`, `emissfire_d01`, `volc_d01`)
- Zero-sized files — check `ls -la` on all linked files
- Missing GOCART background file — most chemistry options require it
- Version mismatch — the tutorial notes that WRF-Chem v3.8 `convert_emiss.exe` had bugs; v3.6 or v3.9+ recommended

### Geogrid Missing Erodibility

**Symptom:** `EROD` field is all zeros in `geo_em.d01.nc`

**Fix:** Ensure `GEOGRID.TBL.ARW_CHEM` is linked (not `GEOGRID.TBL.ARW`), and the erodibility static dataset is downloaded and referenced in `GEOGRID.TBL`.

### Volcanic Ash Not Appearing in Output

**Symptom:** `vash_*` fields remain zero throughout simulation

**Check:**
- `wrfchemv_d01` exists and is linked in the run directory
- `io_form_auxinput13 = 2` is set in namelist
- `emiss_opt_vol = 1` is set
- Look for read confirmation in `rsl.out.0000`: `med_read_wrf_volc_emiss: Read volcanic ash emissions`

### Chemistry Fields Zero in wrfinput

**Symptom:** After `real.exe`, all chemistry arrays are zero

**This is normal** for most cases — chemical species are initialized to zero. They build up during the `wrf.exe` spin-up period. For non-zero initial conditions, use **mozbc** to map global chemistry model output (e.g., MOZART, WACCM) onto the WRF-Chem grid before running `real.exe`.

---

## References

- [WRF-Chem Tutorial Exercises (Hatheway)](https://hathewaywill.github.io/wrf-chem-site/tutorialexercises.htm)
- [WRF-Chem User Guide (NOAA)](https://ruc.noaa.gov/wrf/wrf-chem/Users_guide.pdf)
- [Prep-Chem-Sources Documentation](https://ruc.noaa.gov/wrf/wrf-chem/Emission_guide.pdf)
- [WRF-Chem Emissions Guide](https://ruc.noaa.gov/wrf/wrf-chem/)
