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

### Observations
*Coming soon*
