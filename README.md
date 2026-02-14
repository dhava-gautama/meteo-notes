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

### Observations
*Coming soon*
