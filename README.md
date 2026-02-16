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
- [MET (Model Evaluation Tools): Practical Verification Guide](verification/met-verification-guide.md) — Practical guide to MET v12.0 and METplus v6.0: installation, configuration, core tools (Point-Stat, Grid-Stat, Ensemble-Stat), spatial/object-based methods (MODE, MTD), tropical cyclone verification, Python embedding, METviewer/METplotpy visualization, and end-to-end workflows for WRF, precipitation, ensemble, and TC verification.
- [Python Forecast Verification Guide](verification/python-verification-guide.md) — Lightweight verification using the scientific Python stack (numpy, pandas, xarray, scipy, scikit-learn, matplotlib): continuous, categorical, probabilistic, spatial (FSS), and ensemble (CRPS, rank histogram) metrics with complete code examples and end-to-end WRF verification workflow. No MET/METplus dependency.

### Models
- [WRF Model: Complete Guide](models/wrf-model-guide.md) — End-to-end WRF v4.7 workflow: data download, WPS, physics options, hidden variable output, post-processing, WRF variants (DA, Chem, Fire, Hydro, Solar, MPAS), and JMA NWP reference settings.
- [WRF-Chem Tutorial Guide](models/wrf-chem-tutorial-guide.md) — Hands-on WRF-Chem tutorials: dust simulation (GOCART schemes), volcanic ash dispersal (prep_chem_sources), US anthropogenic emissions (NEI-2011), global multi-source emissions (RETRO/EDGAR/GOCART with RACM-KPP), aerosol radiative feedback experiments, auxiliary input channel reference, and troubleshooting.
- [MPAS-Atmosphere: Complete Guide](models/mpas-model-guide.md) — End-to-end MPAS v8.3 workflow: Voronoi mesh setup, initialization, physics suites, dynamics, streams I/O, regional (limited-area) simulations, post-processing, and visualization.
- [ROMS: Complete Guide](models/roms-model-guide.md) — End-to-end ROMS workflow: grid generation, bathymetry smoothing, atmospheric/tidal/boundary forcing, S-coordinate vertical grid, physics options (GLS/KPP mixing, advection, pressure gradient), nesting, Indonesian waters & Sunda Strait application guide.
- [COAWST: Complete Guide](models/coawst-model-guide.md) — Coupled Ocean-Atmosphere-Wave-Sediment Transport (WRF+ROMS+SWAN) guide: SWAN wave model reference, MCT coupling mechanisms, field exchanges, SCRIP weight generation, compilation, configuration, and Sunda Strait coupled application (Ina-CAWO replication).
- [SWAN: Complete Guide](models/swan-model-guide.md) — End-to-end SWAN v41.51 workflow: spectral action balance equation, grid types (regular/curvilinear/unstructured), all physics options (wind input, whitecapping, breaking, friction, triads, vegetation, ice), numerical methods, input/output, nesting (SWAN-to-SWAN, WAM/WW3-to-SWAN), calibration, and running strategies.
- [WAVEWATCH III: Complete Guide](models/ww3-model-guide.md) — End-to-end WW3 v6.07 workflow: switch file system, all source term packages (ST1-ST6, NL1-NL5, IC0-IC5), grid types (regular/curvilinear/unstructured/SMC), multi-grid mosaic (ww3_multi), nesting, operational use (GFS-Wave, GEFS-Wave), ESMF/NUOPC coupling, and calibration.

### Observations
- [Surface Observations: Complete Guide](observations/surface-obs-guide.md) — SYNOP, METAR, and AWS formats, data sources (ISD, IEM, OGIMET, BMKG), Python decoding (eccodes, metar, siphon), quality control pipeline (range, consistency, spike, buddy, persistence checks), data processing, and preparation for MET/Python verification. Indonesian station network reference.
- [Ocean Observations: Complete Guide](observations/ocean-obs-guide.md) — Moored buoys (NDBC, TAO/TRITON/RAMA), Argo floats, tide gauges, satellite products (SST, altimetry, scatterometer winds), Python access (erddapy, argopy, copernicusmarine, utide), quality control, and preparation for ROMS/SWAN/WW3 verification. Indonesian waters and ITF monitoring reference.
