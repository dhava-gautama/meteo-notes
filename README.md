# meteo-notes

Personal notes on meteorological modeling, observations, and forecast verification.

## Structure

```
meteo-notes/
├── verification/       # Forecast verification methods and metrics
├── models/             # Numerical weather prediction and modeling
└── observations/       # Observational data, instruments, and QC
```

## Contents

### Verification
- [Forecast Verification: Methods and Key Concepts](verification/forecast-verification-notes.md) — Comprehensive notes on verification metrics for dichotomous, multi-category, continuous, and probabilistic forecasts. Covers spatial verification, ensemble diagnostics, and rare event methods.

### Models
- [WRF Model: Complete Guide](models/wrf-model-guide.md) — End-to-end WRF v4.7 workflow: data download, WPS, physics options, hidden variable output, post-processing, WRF variants (DA, Chem, Fire, Hydro, Solar, MPAS), and JMA NWP reference settings.
- [MPAS-Atmosphere: Complete Guide](models/mpas-model-guide.md) — End-to-end MPAS v8.3 workflow: Voronoi mesh setup, initialization, physics suites, dynamics, streams I/O, regional (limited-area) simulations, post-processing, and visualization.
- [ROMS: Complete Guide](models/roms-model-guide.md) — End-to-end ROMS workflow: grid generation, bathymetry smoothing, atmospheric/tidal/boundary forcing, S-coordinate vertical grid, physics options (GLS/KPP mixing, advection, pressure gradient), nesting, Indonesian waters & Sunda Strait application guide.
- [COAWST: Complete Guide](models/coawst-model-guide.md) — Coupled Ocean-Atmosphere-Wave-Sediment Transport (WRF+ROMS+SWAN) guide: SWAN wave model reference, MCT coupling mechanisms, field exchanges, SCRIP weight generation, compilation, configuration, and Sunda Strait coupled application (Ina-CAWO replication).
- [SWAN: Complete Guide](models/swan-model-guide.md) — End-to-end SWAN v41.51 workflow: spectral action balance equation, grid types (regular/curvilinear/unstructured), all physics options (wind input, whitecapping, breaking, friction, triads, vegetation, ice), numerical methods, input/output, nesting (SWAN-to-SWAN, WAM/WW3-to-SWAN), calibration, and running strategies.
- [WAVEWATCH III: Complete Guide](models/ww3-model-guide.md) — End-to-end WW3 v6.07 workflow: switch file system, all source term packages (ST1-ST6, NL1-NL5, IC0-IC5), grid types (regular/curvilinear/unstructured/SMC), multi-grid mosaic (ww3_multi), nesting, operational use (GFS-Wave, GEFS-Wave), ESMF/NUOPC coupling, and calibration.

### Observations
*Coming soon*
