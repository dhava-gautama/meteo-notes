# Python Forecast Verification Guide

> Lightweight forecast verification using the scientific Python stack — no MET/METplus required.
> Complements the theory-focused [Forecast Verification: Methods and Key Concepts](forecast-verification-notes.md)
> and the MET-based [MET Practical Verification Guide](met-verification-guide.md).

---

## Table of Contents

1. [Introduction & Scope](#1-introduction--scope)
2. [Environment Setup](#2-environment-setup)
3. [Loading NWP Data](#3-loading-nwp-data)
   - [GRIB2](#31-grib2-via-cfgrib--xarray)
   - [NetCDF](#32-netcdf-via-xarray)
   - [CSV Station Observations](#33-csv-station-observations)
   - [Matching Forecasts to Observations](#34-matching-forecasts-to-observations)
4. [Continuous Verification Metrics](#4-continuous-verification-metrics)
5. [Categorical (Dichotomous) Verification Metrics](#5-categorical-dichotomous-verification-metrics)
6. [Probabilistic Verification Metrics](#6-probabilistic-verification-metrics)
7. [Spatial Verification](#7-spatial-verification)
8. [Ensemble Verification](#8-ensemble-verification)
9. [Visualization](#9-visualization)
   - [Taylor Diagrams](#91-taylor-diagrams)
   - [Performance Diagrams](#92-performance-diagrams)
   - [Reliability Diagrams](#93-reliability-diagrams)
   - [ROC Curves](#94-roc-curves)
   - [Rank Histograms](#95-rank-histograms)
   - [Spatial Error Maps](#96-spatial-error-maps)
10. [End-to-End Example: WRF vs. Station Observations](#10-end-to-end-example-wrf-vs-station-observations)
11. [Tips & Best Practices](#11-tips--best-practices)
- [Appendix A: Metric Quick Reference](#appendix-a-metric-quick-reference)
- [Appendix B: Useful Resources](#appendix-b-useful-resources)

---

## 1. Introduction & Scope

This guide shows how to perform forecast verification **entirely in Python** using standard scientific libraries. No compilation, no Docker, no MET installation — just `pip install` and go.

For formula derivations and metric interpretation, see the companion [Forecast Verification: Methods and Key Concepts](forecast-verification-notes.md). This guide focuses purely on **implementation**.

### When to Use Python-Only vs. MET/METplus

| | **Python-Only** | **MET/METplus** |
|---|---|---|
| **Setup** | `pip install` (~2 min) | Compile C++ or Docker (~30 min) |
| **Best for** | Quick studies, prototyping, teaching, custom workflows | Operational monitoring, formal studies, large-scale verification |
| **Spatial methods** | FSS (manual), basic SAL | MODE, MTD, Wavelet-Stat (full suite) |
| **TC verification** | Not practical | TC-Pairs, TC-Stat, TC-Gen, TC-RMW |
| **Ensemble** | Rank hist, CRPS, spread-skill | Full suite + aggregation + PIT |
| **Aggregation** | Manual (pandas groupby) | stat_analysis (built-in, partial sums) |
| **Scalability** | Good with xarray/dask | Excellent (optimized C++) |
| **Customization** | Total control | Config-driven, less flexible |

**Rule of thumb:** Start with Python for exploration. Graduate to MET when you need operational reliability, MODE/MTD, or TC tools.

---

## 2. Environment Setup

### Conda (Recommended)

```yaml
# environment.yml
name: verification
channels:
  - conda-forge
dependencies:
  - python=3.11
  - numpy
  - pandas
  - xarray
  - scipy
  - scikit-learn
  - matplotlib
  - cartopy
  - cfgrib
  - eccodes
  - netcdf4
  - pip
  - pip:
    - properscoring
```

```bash
conda env create -f environment.yml
conda activate verification
```

### Pip Alternative

```bash
pip install numpy pandas xarray scipy scikit-learn matplotlib cartopy \
            cfgrib eccodes netCDF4 properscoring
```

### Package Purposes

| Package | Purpose |
|---|---|
| **numpy** | Array operations, core metric computation |
| **pandas** | Station data, time series, groupby aggregation |
| **xarray** | Gridded data I/O (GRIB2, NetCDF) and manipulation |
| **cfgrib** | GRIB2 backend for xarray |
| **scipy** | Interpolation, KD-tree, statistical tests, bootstrapping |
| **scikit-learn** | ROC curves, AUC, confusion matrix |
| **matplotlib** | All plotting |
| **cartopy** | Map projections for spatial plots |
| **properscoring** | CRPS and other probabilistic scores |

### Verify Installation

```python
import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats
import sklearn.metrics
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import properscoring as ps
print("All packages OK")
```

---

## 3. Loading NWP Data

### 3.1 GRIB2 via cfgrib + xarray

```python
import xarray as xr

# Load 2-m temperature from a GFS GRIB2 file
ds = xr.open_dataset(
    'gfs.t00z.pgrb2.0p25.f024',
    engine='cfgrib',
    backend_kwargs={
        'filter_by_keys': {
            'shortName': '2t',
            'typeOfLevel': 'heightAboveGround',
            'level': 2
        }
    }
)
t2m = ds['t2m']  # (latitude, longitude), units: K

# Load precipitation (total precip)
ds_precip = xr.open_dataset(
    'gfs.t00z.pgrb2.0p25.f024',
    engine='cfgrib',
    backend_kwargs={
        'filter_by_keys': {'shortName': 'tp'}
    }
)
precip = ds_precip['tp']  # kg/m² (= mm)

# Load multiple pressure levels
ds_upper = xr.open_dataset(
    'gfs.t00z.pgrb2.0p25.f024',
    engine='cfgrib',
    backend_kwargs={
        'filter_by_keys': {
            'shortName': 't',
            'typeOfLevel': 'isobaricInhPa'
        }
    }
)
t850 = ds_upper['t'].sel(isobaricInhPa=850)  # 850 hPa temperature
```

**Tip:** Use `wgrib2 -s file.grb2 | head` or `cfgrib.open_datasets('file.grb2')` to discover available fields and their shortNames.

### 3.2 NetCDF via xarray

```python
# Standard CF-compliant NetCDF
ds = xr.open_dataset('analysis.nc')
t2m = ds['TMP_2mAboveGround']

# WRF output (non-CF naming)
wrf = xr.open_dataset('wrfout_d01_2024-01-15_12:00:00')
t2m_wrf = wrf['T2']              # 2-m temperature (K)
rain = wrf['RAINNC'] + wrf['RAINC']  # Total accumulated precip (mm)
u10 = wrf['U10']                 # 10-m U wind (m/s)
v10 = wrf['V10']                 # 10-m V wind (m/s)
wspd10 = np.sqrt(u10**2 + v10**2)    # 10-m wind speed

lats = wrf['XLAT'].isel(Time=0)
lons = wrf['XLONG'].isel(Time=0)

# Derive 6-h precipitation from accumulated field (two time steps)
wrf_t1 = xr.open_dataset('wrfout_d01_2024-01-15_06:00:00')
wrf_t2 = xr.open_dataset('wrfout_d01_2024-01-15_12:00:00')
precip_6h = (wrf_t2['RAINNC'] + wrf_t2['RAINC']) - \
            (wrf_t1['RAINNC'] + wrf_t1['RAINC'])
```

### 3.3 CSV Station Observations

```python
import pandas as pd

# Typical station observation CSV
# station_id, datetime, lat, lon, elev, temperature_K, rh_pct, wspd_ms, precip_mm
obs = pd.read_csv('surface_obs.csv', parse_dates=['datetime'])
obs = obs.dropna(subset=['temperature_K'])

# Filter to a specific valid time
valid_time = pd.Timestamp('2024-01-15 12:00:00')
obs_valid = obs[obs['datetime'] == valid_time].copy()

print(f"Stations: {len(obs_valid)}")
print(obs_valid[['station_id', 'lat', 'lon', 'temperature_K']].head())
```

### 3.4 Matching Forecasts to Observations

This is the critical step that MET handles internally. You must extract forecast values at each observation location.

#### Nearest-Neighbor (fast, recommended for most cases)

```python
from scipy.spatial import cKDTree

def extract_nearest(fcst_field, fcst_lats, fcst_lons, obs_lats, obs_lons):
    """Extract forecast values at observation locations using nearest neighbor.

    Parameters:
        fcst_field: 2D array (ny, nx) of forecast values
        fcst_lats:  2D array (ny, nx) of latitudes
        fcst_lons:  2D array (ny, nx) of longitudes
        obs_lats:   1D array of observation latitudes
        obs_lons:   1D array of observation longitudes

    Returns:
        1D array of forecast values at observation locations
    """
    grid_points = np.column_stack([fcst_lats.ravel(), fcst_lons.ravel()])
    tree = cKDTree(grid_points)

    obs_points = np.column_stack([obs_lats, obs_lons])
    dist, idx = tree.query(obs_points)

    return fcst_field.values.ravel()[idx]

# Usage
fcst_at_obs = extract_nearest(
    t2m_wrf.isel(Time=0), lats, lons,
    obs_valid['lat'].values, obs_valid['lon'].values
)
```

#### Bilinear Interpolation (smoother, better for continuous fields)

```python
from scipy.interpolate import RegularGridInterpolator

def extract_bilinear(fcst_field, lat_1d, lon_1d, obs_lats, obs_lons):
    """Bilinear interpolation for regular lat-lon grids.

    Parameters:
        fcst_field: 2D array (nlat, nlon)
        lat_1d:     1D array of latitudes (must be monotonic)
        lon_1d:     1D array of longitudes (must be monotonic)
        obs_lats:   1D array of observation latitudes
        obs_lons:   1D array of observation longitudes

    Returns:
        1D array of interpolated forecast values
    """
    interp = RegularGridInterpolator(
        (lat_1d, lon_1d), fcst_field,
        method='linear', bounds_error=False, fill_value=np.nan
    )
    return interp(np.column_stack([obs_lats, obs_lons]))

# Usage (for regular grids like GFS)
fcst_at_obs = extract_bilinear(
    t2m.values, t2m.latitude.values, t2m.longitude.values,
    obs_valid['lat'].values, obs_valid['lon'].values
)
```

---

## 4. Continuous Verification Metrics

> For metric definitions and interpretation, see [Section 6.3 of the verification theory notes](forecast-verification-notes.md#63-methods-for-continuous-forecasts).

### Core Metrics

```python
def continuous_metrics(fcst, obs):
    """Compute standard continuous verification metrics.

    Parameters:
        fcst: array-like of forecast values
        obs:  array-like of observation values

    Returns:
        dict with ME, MAE, RMSE, r, multiplicative bias, skill score, N
    """
    fcst, obs = np.asarray(fcst, dtype=float), np.asarray(obs, dtype=float)
    mask = ~(np.isnan(fcst) | np.isnan(obs))
    f, o = fcst[mask], obs[mask]

    err = f - o
    n = len(f)
    me = np.mean(err)
    mae = np.mean(np.abs(err))
    mse = np.mean(err**2)
    rmse = np.sqrt(mse)
    r = np.corrcoef(f, o)[0, 1] if n > 2 else np.nan
    mult_bias = np.mean(f) / np.mean(o) if np.mean(o) != 0 else np.nan

    # Skill score vs. climatology (using obs mean as reference)
    mse_ref = np.var(o)
    ss = 1 - mse / mse_ref if mse_ref > 0 else np.nan

    return {
        'N': n, 'ME': me, 'MAE': mae, 'RMSE': rmse,
        'r': r, 'Mult_Bias': mult_bias, 'Skill_Score': ss,
        'FBAR': np.mean(f), 'OBAR': np.mean(o)
    }
```

### Example Usage

```python
results = continuous_metrics(fcst_at_obs, obs_valid['temperature_K'].values)
df = pd.DataFrame([results])
print(df.to_string(index=False))

#   N      ME     MAE    RMSE      r  Mult_Bias  Skill_Score    FBAR    OBAR
# 150  -0.32    1.45    1.87  0.952      0.999        0.812  298.12  298.44
```

### MSE Decomposition

Separates error into **unconditional bias**, **conditional bias (amplitude)**, and **phase error** components — see [theory notes](forecast-verification-notes.md#63-methods-for-continuous-forecasts).

```python
def mse_decomposition(fcst, obs):
    """Decompose MSE into bias, amplitude, and phase components.

    MSE = bias² + (std_f - std_o)² + 2·std_f·std_o·(1 - r)
    """
    f, o = np.asarray(fcst, dtype=float), np.asarray(obs, dtype=float)
    r = np.corrcoef(f, o)[0, 1]

    bias_sq = (np.mean(f) - np.mean(o))**2
    amplitude = (np.std(f) - np.std(o))**2
    phase = 2 * np.std(f) * np.std(o) * (1 - r)

    return {
        'MSE': bias_sq + amplitude + phase,
        'bias_squared': bias_sq,
        'amplitude_error': amplitude,
        'phase_error': phase
    }
```

### Bootstrap Confidence Intervals

```python
from scipy.stats import bootstrap

def bootstrap_metric(fcst, obs, metric_func, n_resamples=1000, confidence=0.95):
    """Bootstrap confidence interval for any metric.

    Parameters:
        fcst, obs: arrays of forecast and observation values
        metric_func: function(fcst, obs) → scalar
        n_resamples: number of bootstrap replicates
        confidence: confidence level (e.g., 0.95)

    Returns:
        (estimate, lower_bound, upper_bound)
    """
    def statistic(*args, axis):
        idx = args[0].astype(int)
        return np.array([metric_func(fcst[idx], obs[idx])])

    indices = (np.arange(len(fcst)),)
    res = bootstrap(indices, statistic, n_resamples=n_resamples,
                    confidence_level=confidence, method='percentile')

    estimate = metric_func(fcst, obs)
    return estimate, res.confidence_interval.low, res.confidence_interval.high

# Example: bootstrap RMSE
rmse_func = lambda f, o: np.sqrt(np.mean((f - o)**2))
est, lo, hi = bootstrap_metric(fcst_at_obs, obs_valid['temperature_K'].values, rmse_func)
print(f"RMSE = {est:.2f} [{lo:.2f}, {hi:.2f}]")
```

### Partial Sums (for Correct Aggregation)

Never average RMSE directly — aggregate partial sums first.

```python
def partial_sums(fcst, obs):
    """Compute SL1L2 partial sums for later aggregation (matches MET SL1L2)."""
    f, o = np.asarray(fcst, dtype=float), np.asarray(obs, dtype=float)
    return {
        'N': len(f),
        'FBAR': np.mean(f),
        'OBAR': np.mean(o),
        'FFBAR': np.mean(f**2),
        'FOBAR': np.mean(f * o),
        'OOBAR': np.mean(o**2)
    }

def aggregate_partial_sums(sums_list):
    """Aggregate a list of partial sum dicts, then derive metrics."""
    total_n = sum(s['N'] for s in sums_list)
    fbar  = sum(s['FBAR'] * s['N'] for s in sums_list) / total_n
    obar  = sum(s['OBAR'] * s['N'] for s in sums_list) / total_n
    ffbar = sum(s['FFBAR'] * s['N'] for s in sums_list) / total_n
    fobar = sum(s['FOBAR'] * s['N'] for s in sums_list) / total_n
    oobar = sum(s['OOBAR'] * s['N'] for s in sums_list) / total_n

    me = fbar - obar
    mse = ffbar - 2 * fobar + oobar
    rmse = np.sqrt(mse) if mse >= 0 else np.nan

    return {'N': total_n, 'ME': me, 'RMSE': rmse, 'FBAR': fbar, 'OBAR': obar}
```

---

## 5. Categorical (Dichotomous) Verification Metrics

> For the contingency table framework, see [Section 6.1 of the verification theory notes](forecast-verification-notes.md#61-methods-for-dichotomous-yesno-forecasts).

### Contingency Table

```python
def contingency_table(fcst, obs, threshold):
    """Build 2x2 contingency table for a given threshold.

    Returns:
        a (hits), b (false alarms), c (misses), d (correct negatives)
    """
    f = np.asarray(fcst) >= threshold
    o = np.asarray(obs) >= threshold
    a = int(np.sum(f & o))          # hits
    b = int(np.sum(f & ~o))         # false alarms
    c = int(np.sum(~f & o))         # misses
    d = int(np.sum(~f & ~o))        # correct negatives
    return a, b, c, d
```

### All Categorical Metrics

```python
def categorical_metrics(a, b, c, d):
    """Compute categorical verification metrics from contingency table counts.

    Parameters:
        a: hits, b: false alarms, c: misses, d: correct negatives
    """
    n = a + b + c + d
    a_ref = (a + b) * (a + c) / n if n > 0 else 0  # expected random hits

    pod  = a / (a + c) if (a + c) > 0 else np.nan       # Probability of Detection
    far  = b / (a + b) if (a + b) > 0 else np.nan       # False Alarm Ratio
    pofd = b / (b + d) if (b + d) > 0 else np.nan       # Prob. of False Detection
    sr   = a / (a + b) if (a + b) > 0 else np.nan       # Success Ratio (1 - FAR)
    csi  = a / (a + b + c) if (a + b + c) > 0 else np.nan   # Critical Success Index
    fbias = (a + b) / (a + c) if (a + c) > 0 else np.nan    # Frequency Bias

    # Equitable Threat Score (Gilbert Skill Score)
    denom = a + b + c - a_ref
    ets = (a - a_ref) / denom if denom != 0 else np.nan

    # Hanssen-Kuipers Discriminant (True Skill Statistic)
    hk = pod - pofd if not (np.isnan(pod) or np.isnan(pofd)) else np.nan

    # Heidke Skill Score
    hss_denom = (a + c) * (c + d) + (a + b) * (b + d)
    hss = 2 * (a * d - b * c) / hss_denom if hss_denom != 0 else np.nan

    # Accuracy (fraction correct)
    acc = (a + d) / n if n > 0 else np.nan

    return {
        'N': n, 'Hits': a, 'FA': b, 'Misses': c, 'CN': d,
        'POD': pod, 'FAR': far, 'POFD': pofd, 'SR': sr,
        'CSI': csi, 'ETS': ets, 'FBIAS': fbias,
        'HK': hk, 'HSS': hss, 'ACC': acc
    }
```

### Multi-Threshold Evaluation

```python
thresholds = [0.1, 1, 5, 10, 25, 50]  # mm for precipitation
results = []
for t in thresholds:
    a, b, c, d = contingency_table(fcst_precip, obs_precip, t)
    metrics = categorical_metrics(a, b, c, d)
    metrics['threshold_mm'] = t
    results.append(metrics)

df = pd.DataFrame(results)
print(df[['threshold_mm', 'POD', 'FAR', 'CSI', 'ETS', 'FBIAS']].to_string(index=False))

#  threshold_mm   POD    FAR    CSI    ETS  FBIAS
#           0.1  0.95   0.12   0.84   0.72   1.08
#           1.0  0.88   0.15   0.77   0.62   1.04
#           5.0  0.74   0.22   0.62   0.45   0.95
#          10.0  0.61   0.30   0.49   0.32   0.87
#          25.0  0.42   0.38   0.32   0.18   0.68
#          50.0  0.21   0.50   0.15   0.06   0.42
```

### Multi-Category Contingency Table

```python
from sklearn.metrics import confusion_matrix, classification_report

# Define categories (e.g., precipitation intensity)
def categorize_precip(values):
    """Categorize precipitation: None / Light / Moderate / Heavy."""
    cats = np.full(len(values), 'None', dtype=object)
    cats[values >= 0.1] = 'Light'
    cats[values >= 5.0] = 'Moderate'
    cats[values >= 25.0] = 'Heavy'
    return cats

fcst_cat = categorize_precip(fcst_precip)
obs_cat = categorize_precip(obs_precip)

labels = ['None', 'Light', 'Moderate', 'Heavy']
cm = confusion_matrix(obs_cat, fcst_cat, labels=labels)
print("Confusion Matrix:")
print(pd.DataFrame(cm, index=[f'Obs:{l}' for l in labels],
                       columns=[f'Fcst:{l}' for l in labels]))
```

---

## 6. Probabilistic Verification Metrics

> For metric definitions, see [Section 6.4 of the verification theory notes](forecast-verification-notes.md#64-methods-for-probabilistic-forecasts).

### Brier Score and Brier Skill Score

```python
def brier_score(prob_fcst, obs_binary):
    """Brier Score: mean squared error of probability forecasts.

    Parameters:
        prob_fcst:  array of forecast probabilities [0, 1]
        obs_binary: array of binary observations (0 or 1)

    Returns:
        Brier Score (0 = perfect, 1 = worst)
    """
    return np.mean((np.asarray(prob_fcst) - np.asarray(obs_binary))**2)

def brier_skill_score(prob_fcst, obs_binary):
    """BSS relative to climatological probability."""
    bs = brier_score(prob_fcst, obs_binary)
    climo_prob = np.mean(obs_binary)
    bs_ref = brier_score(np.full_like(prob_fcst, climo_prob), obs_binary)
    return 1 - bs / bs_ref if bs_ref > 0 else np.nan
```

### Brier Score Decomposition

```python
def brier_decomposition(prob_fcst, obs_binary, n_bins=10):
    """Decompose BS into reliability, resolution, and uncertainty.

    BS = reliability - resolution + uncertainty
    """
    p = np.asarray(prob_fcst, dtype=float)
    o = np.asarray(obs_binary, dtype=float)
    n = len(p)
    climo = np.mean(o)

    bin_edges = np.linspace(0, 1, n_bins + 1)
    reliability = 0.0
    resolution = 0.0

    for i in range(n_bins):
        mask = (p >= bin_edges[i]) & (p < bin_edges[i + 1])
        if i == n_bins - 1:  # include right edge in last bin
            mask = mask | (p == bin_edges[i + 1])
        nk = mask.sum()
        if nk == 0:
            continue
        ok = o[mask].mean()
        pk = p[mask].mean()
        reliability += nk * (pk - ok)**2
        resolution += nk * (ok - climo)**2

    reliability /= n
    resolution /= n
    uncertainty = climo * (1 - climo)

    return {
        'BS': reliability - resolution + uncertainty,
        'reliability': reliability,
        'resolution': resolution,
        'uncertainty': uncertainty
    }
```

### Reliability Diagram Data

```python
def reliability_data(prob_fcst, obs_binary, n_bins=10):
    """Compute data for a reliability diagram.

    Returns:
        bin_centers: forecast probability bin centers
        obs_freq:    observed relative frequency per bin
        counts:      number of forecasts per bin
    """
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    obs_freq = np.full(n_bins, np.nan)
    counts = np.zeros(n_bins, dtype=int)

    for i in range(n_bins):
        if i < n_bins - 1:
            mask = (prob_fcst >= bin_edges[i]) & (prob_fcst < bin_edges[i + 1])
        else:
            mask = (prob_fcst >= bin_edges[i]) & (prob_fcst <= bin_edges[i + 1])
        counts[i] = mask.sum()
        if counts[i] > 0:
            obs_freq[i] = obs_binary[mask].mean()

    return bin_centers, obs_freq, counts
```

### ROC Curve and AUC

```python
from sklearn.metrics import roc_curve, roc_auc_score

def compute_roc(prob_fcst, obs_binary):
    """Compute ROC curve and AUC.

    Returns:
        fpr: false positive rates
        tpr: true positive rates (= POD)
        auc: area under curve
    """
    fpr, tpr, _ = roc_curve(obs_binary, prob_fcst)
    auc = roc_auc_score(obs_binary, prob_fcst)
    return fpr, tpr, auc
```

### Ranked Probability Score (Multi-Category)

```python
def ranked_probability_score(prob_categories, obs_category, n_categories):
    """RPS for a single multi-category probabilistic forecast.

    Parameters:
        prob_categories: array of probabilities for each category (sums to 1)
        obs_category:    integer index of observed category (0-based)
        n_categories:    total number of categories

    Returns:
        RPS value (0 = perfect)
    """
    cum_fcst = np.cumsum(prob_categories)
    cum_obs = np.zeros(n_categories)
    cum_obs[obs_category:] = 1.0
    return np.sum((cum_fcst - cum_obs)**2) / (n_categories - 1)

# Example: 3 categories (Below / Near / Above normal)
rps = ranked_probability_score(
    prob_categories=[0.2, 0.5, 0.3],  # forecast probabilities
    obs_category=1,                     # observed: "Near normal"
    n_categories=3
)
print(f"RPS = {rps:.3f}")
```

---

## 7. Spatial Verification

> For the conceptual framework, see [Section 7.1 of the verification theory notes](forecast-verification-notes.md#71-spatial-forecast-verification).

### Fractions Skill Score (FSS)

FSS addresses the **double penalty problem** by comparing fractional coverage within neighborhoods rather than grid-point values.

```python
from scipy.ndimage import uniform_filter

def fractions_skill_score(fcst_field, obs_field, threshold, window_size):
    """Compute Fractions Skill Score.

    Parameters:
        fcst_field:  2D array of forecast values
        obs_field:   2D array of observation values
        threshold:   value threshold for binary conversion
        window_size: neighborhood size (odd integer, in grid points)

    Returns:
        FSS value (0 = no skill, 1 = perfect)
    """
    fcst_binary = (np.asarray(fcst_field) >= threshold).astype(float)
    obs_binary = (np.asarray(obs_field) >= threshold).astype(float)

    # Compute fractions using uniform (box) filter
    fcst_frac = uniform_filter(fcst_binary, size=window_size, mode='constant')
    obs_frac = uniform_filter(obs_binary, size=window_size, mode='constant')

    mse = np.mean((fcst_frac - obs_frac)**2)
    mse_ref = np.mean(fcst_frac**2) + np.mean(obs_frac**2)

    return 1 - mse / mse_ref if mse_ref > 0 else 1.0
```

### FSS as a Function of Scale

```python
def fss_by_scale(fcst_field, obs_field, threshold, window_sizes):
    """Compute FSS at multiple spatial scales.

    Returns:
        list of (window_size, FSS) tuples
    """
    return [(w, fractions_skill_score(fcst_field, obs_field, threshold, w))
            for w in window_sizes]

# Example
window_sizes = [1, 3, 5, 11, 21, 41, 81, 161]
fss_results = fss_by_scale(fcst_precip_grid, obs_precip_grid,
                           threshold=1.0, window_sizes=window_sizes)

for w, fss in fss_results:
    print(f"  Window {w:3d}: FSS = {fss:.3f}")
```

### Minimum Useful Scale

The forecast becomes "useful" at the scale where FSS exceeds `0.5 + f_o/2`, where `f_o` is the observed base rate.

```python
def minimum_useful_scale(fss_results, obs_field, threshold):
    """Find the minimum scale at which the forecast has useful skill.

    Parameters:
        fss_results: list of (window_size, FSS)
        obs_field:   2D observation array
        threshold:   threshold used for FSS computation

    Returns:
        minimum useful window size (or None if never reached)
    """
    f_o = np.mean(obs_field >= threshold)
    useful_threshold = 0.5 + f_o / 2

    for w, fss in fss_results:
        if fss >= useful_threshold:
            return w
    return None

min_scale = minimum_useful_scale(fss_results, obs_precip_grid, threshold=1.0)
print(f"Minimum useful scale: {min_scale} grid points")
```

### SAL (Structure-Amplitude-Location) — Amplitude Component

Full SAL requires object identification (complex). The **amplitude** component is straightforward:

```python
def sal_amplitude(fcst_field, obs_field):
    """SAL Amplitude component: normalized difference in domain-average.

    A = (D_fcst - D_obs) / (0.5 * (D_fcst + D_obs))
    Range: -2 to +2, perfect = 0
    Positive = forecast overestimates, negative = underestimates
    """
    d_f = np.mean(fcst_field)
    d_o = np.mean(obs_field)
    denom = 0.5 * (d_f + d_o)
    return (d_f - d_o) / denom if denom > 0 else 0.0
```

---

## 8. Ensemble Verification

> For rank histogram interpretation and ensemble diagnostics, see [Section 7.2 of the verification theory notes](forecast-verification-notes.md#72-ensembleprobabilistic-diagnostics).

### Rank Histogram

A flat rank histogram indicates a well-calibrated ensemble. U-shaped = underdispersive, dome-shaped = overdispersive.

```python
def rank_histogram(ensemble, obs):
    """Compute rank histogram from ensemble and observations.

    Parameters:
        ensemble: 2D array (n_times, n_members)
        obs:      1D array (n_times,)

    Returns:
        ranks: 1D array of counts for each rank (length = n_members + 1)
    """
    n_times, n_members = ensemble.shape
    ranks = np.zeros(n_members + 1, dtype=int)

    for t in range(n_times):
        combined = np.append(ensemble[t, :], obs[t])
        # Rank of observation among ensemble + obs
        rank = np.sum(ensemble[t, :] < obs[t])
        # Randomize ties
        n_ties = np.sum(ensemble[t, :] == obs[t])
        if n_ties > 0:
            rank += np.random.randint(0, n_ties + 1)
        ranks[rank] += 1

    return ranks
```

### Spread-Skill Relationship

```python
def spread_skill(ensemble, obs, n_bins=10):
    """Compute binned spread-skill relationship.

    Parameters:
        ensemble: 2D array (n_times, n_members)
        obs:      1D array (n_times,)
        n_bins:   number of spread bins

    Returns:
        DataFrame with columns: spread_bin, mean_spread, rmse, count
    """
    ens_mean = np.mean(ensemble, axis=1)
    ens_spread = np.std(ensemble, axis=1, ddof=1)
    errors = ens_mean - obs

    # Overall stats
    total_spread = np.mean(ens_spread)
    total_rmse = np.sqrt(np.mean(errors**2))

    # Binned relationship
    bin_edges = np.percentile(ens_spread, np.linspace(0, 100, n_bins + 1))
    results = []
    for i in range(n_bins):
        mask = (ens_spread >= bin_edges[i]) & (ens_spread < bin_edges[i + 1])
        if i == n_bins - 1:
            mask = mask | (ens_spread == bin_edges[i + 1])
        if mask.sum() > 0:
            results.append({
                'spread_bin': i + 1,
                'mean_spread': ens_spread[mask].mean(),
                'rmse': np.sqrt(np.mean(errors[mask]**2)),
                'count': mask.sum()
            })

    return pd.DataFrame(results), total_spread, total_rmse
```

### CRPS (Continuous Ranked Probability Score)

```python
# Using properscoring (recommended)
import properscoring as ps

def crps_ensemble(obs, ensemble):
    """Compute CRPS for ensemble forecasts.

    Parameters:
        obs:      1D array (n_times,) of observations
        ensemble: 2D array (n_times, n_members) of ensemble values

    Returns:
        mean_crps, individual_crps
    """
    crps_values = ps.crps_ensemble(obs, ensemble)
    return np.mean(crps_values), crps_values

# CRPS Skill Score relative to climatology
def crpss(obs, ensemble, climo_mean, climo_std):
    """CRPS Skill Score relative to Gaussian climatology."""
    crps_ens = np.mean(ps.crps_ensemble(obs, ensemble))
    crps_climo = np.mean(ps.crps_gaussian(obs, mu=climo_mean, sig=climo_std))
    return 1 - crps_ens / crps_climo if crps_climo > 0 else np.nan
```

### Manual CRPS (Without properscoring)

```python
def crps_ensemble_manual(obs, ensemble):
    """CRPS from ensemble members using the exact formulation.

    CRPS = E|X - y| - 0.5 * E|X - X'|
    where X, X' are independent draws from the ensemble distribution.
    """
    n_times, n_members = ensemble.shape
    crps = np.zeros(n_times)

    for t in range(n_times):
        ens = ensemble[t, :]
        # First term: mean absolute error
        mae = np.mean(np.abs(ens - obs[t]))
        # Second term: mean absolute difference between members
        diff = np.abs(ens[:, None] - ens[None, :])
        mad = np.mean(diff)
        crps[t] = mae - 0.5 * mad

    return np.mean(crps), crps
```

### PIT Histogram (Probability Integral Transform)

```python
def pit_histogram(ensemble, obs, n_bins=10):
    """Compute PIT histogram for ensemble calibration assessment.

    PIT = fraction of ensemble members below the observation.
    A calibrated ensemble produces a uniform PIT histogram.
    """
    n_times, n_members = ensemble.shape
    pit_values = np.zeros(n_times)

    for t in range(n_times):
        pit_values[t] = np.mean(ensemble[t, :] < obs[t])
        # Add half of ties for continuity
        pit_values[t] += 0.5 * np.mean(ensemble[t, :] == obs[t])

    counts, bin_edges = np.histogram(pit_values, bins=n_bins, range=(0, 1))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, counts
```

---

## 9. Visualization

### 9.1 Taylor Diagrams

The Taylor diagram summarizes correlation, standard deviation ratio, and centered RMSD on a single polar plot — see [theory notes](forecast-verification-notes.md#74-other-diagnostic-tools).

```python
def taylor_diagram(models, ref_std, ax=None, colors=None):
    """Plot a Taylor diagram.

    Parameters:
        models:  list of dicts with keys 'name', 'std', 'corrcoef'
        ref_std: standard deviation of the observations
        ax:      optional matplotlib polar axes
        colors:  optional list of colors
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, polar=True)

    # Limit to first quadrant (correlation 0 to 1)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)

    # Correlation as angle
    angles = np.arccos([m['corrcoef'] for m in models])

    # Normalized standard deviation as radius
    radii = [m['std'] / ref_std for m in models]

    # Plot reference point (observations)
    ax.plot(0, 1, 'ko', markersize=10, label='Observation')

    # Plot models
    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(models)))
    for i, m in enumerate(models):
        ax.plot(angles[i], radii[i], 'o', color=colors[i],
                markersize=8, label=m['name'])

    # RMSD contours (centered)
    rs = np.linspace(0, 2, 100)
    for rmsd in [0.25, 0.5, 0.75, 1.0, 1.5]:
        thetas = np.arccos(np.clip((1 + rs**2 - rmsd**2) / (2 * rs), -1, 1))
        valid = ~np.isnan(thetas)
        ax.plot(thetas[valid], rs[valid], '--', color='gray', alpha=0.4, linewidth=0.8)

    ax.set_rlabel_position(0)
    ax.set_xlabel('Standard Deviation Ratio')
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
    return ax

# Usage
models = [
    {'name': 'WRF-3km',  'std': 2.8, 'corrcoef': 0.95},
    {'name': 'WRF-9km',  'std': 2.5, 'corrcoef': 0.91},
    {'name': 'GFS-25km', 'std': 2.1, 'corrcoef': 0.85},
]
taylor_diagram(models, ref_std=3.0)
plt.title('Temperature Verification')
plt.savefig('taylor_diagram.png', dpi=150, bbox_inches='tight')
```

### 9.2 Performance Diagrams

Performance diagram (Roebber 2009): POD vs Success Ratio (1-FAR) with CSI contours and frequency bias lines.

```python
def performance_diagram(results_list, ax=None):
    """Plot a performance diagram.

    Parameters:
        results_list: list of dicts, each with keys:
            'name', 'thresholds', 'POD', 'SR' (Success Ratio = 1-FAR)
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # CSI contours
    sr_grid = np.linspace(0.01, 1, 100)
    pod_grid = np.linspace(0.01, 1, 100)
    SR, POD = np.meshgrid(sr_grid, pod_grid)
    CSI = 1 / (1/SR + 1/POD - 1)
    cs = ax.contour(SR, POD, CSI, levels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                    colors='gray', alpha=0.4, linewidths=0.8)
    ax.clabel(cs, fmt='%.1f', fontsize=8)

    # Frequency bias lines
    for bias in [0.25, 0.5, 1.0, 2.0, 4.0]:
        pod_line = bias * sr_grid
        valid = pod_line <= 1
        ax.plot(sr_grid[valid], pod_line[valid], '--', color='gray',
                alpha=0.3, linewidth=0.8)
        # Label
        idx = np.argmin(np.abs(pod_line - 0.95))
        if pod_line[idx] <= 1:
            ax.text(sr_grid[idx], pod_line[idx], f'{bias}',
                    fontsize=7, color='gray', ha='center')

    # Plot models
    colors = plt.cm.tab10(np.linspace(0, 1, len(results_list)))
    for i, r in enumerate(results_list):
        ax.plot(r['SR'], r['POD'], 'o-', color=colors[i], label=r['name'],
                markersize=6, linewidth=1.5)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Success Ratio (1 - FAR)')
    ax.set_ylabel('Probability of Detection (POD)')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax
```

### 9.3 Reliability Diagrams

```python
def plot_reliability(prob_fcst, obs_binary, n_bins=10, ax=None):
    """Plot a reliability diagram with sharpness histogram."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))

    centers, obs_freq, counts = reliability_data(prob_fcst, obs_binary, n_bins)
    climo = np.mean(obs_binary)

    # Perfect reliability line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Perfect reliability')

    # No-skill region (shaded)
    ax.fill_between([0, 1], [climo, climo], [0, 1],
                    alpha=0.08, color='gray', label='No skill')

    # Climatology line
    ax.axhline(climo, color='gray', alpha=0.3, linewidth=0.8)
    ax.axvline(climo, color='gray', alpha=0.3, linewidth=0.8)

    # Reliability curve
    valid = counts > 0
    ax.plot(centers[valid], obs_freq[valid], 'o-', color='C0',
            markersize=6, linewidth=1.5, label='Model')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Forecast Probability')
    ax.set_ylabel('Observed Relative Frequency')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Inset: sharpness histogram
    ax_inset = ax.inset_axes([0.55, 0.05, 0.4, 0.25])
    ax_inset.bar(centers, counts, width=1/n_bins * 0.8, color='C0', alpha=0.5)
    ax_inset.set_xlabel('Forecast Prob.', fontsize=8)
    ax_inset.set_ylabel('Count', fontsize=8)
    ax_inset.tick_params(labelsize=7)

    return ax
```

### 9.4 ROC Curves

```python
def plot_roc(prob_fcst, obs_binary, label='Model', ax=None):
    """Plot ROC curve with AUC annotation."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))

    fpr, tpr, auc = compute_roc(prob_fcst, obs_binary)

    ax.plot(fpr, tpr, '-', linewidth=2, label=f'{label} (AUC={auc:.3f})')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.4, label='No skill')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('False Positive Rate (POFD)')
    ax.set_ylabel('True Positive Rate (POD)')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax
```

### 9.5 Rank Histograms

```python
def plot_rank_histogram(ensemble, obs, ax=None):
    """Plot rank histogram with uniform reference line."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    ranks = rank_histogram(ensemble, obs)
    n_ranks = len(ranks)
    expected = len(obs) / n_ranks

    ax.bar(range(n_ranks), ranks, color='C0', alpha=0.7, edgecolor='black')
    ax.axhline(expected, color='red', linestyle='--', linewidth=1.5,
               label=f'Expected ({expected:.0f})')

    ax.set_xlabel('Rank')
    ax.set_ylabel('Count')
    ax.set_xticks(range(0, n_ranks, max(1, n_ranks // 10)))
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    return ax
```

### 9.6 Spatial Error Maps

```python
def plot_spatial_comparison(fcst_field, obs_field, lats, lons,
                            field_name='Temperature (K)', figsize=(18, 5)):
    """Plot forecast, observation, and error maps side by side."""
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    vmin = min(np.nanmin(fcst_field), np.nanmin(obs_field))
    vmax = max(np.nanmax(fcst_field), np.nanmax(obs_field))

    fig, axes = plt.subplots(1, 3, figsize=figsize,
                             subplot_kw={'projection': ccrs.PlateCarree()})

    # Forecast
    c1 = axes[0].pcolormesh(lons, lats, fcst_field, vmin=vmin, vmax=vmax,
                             cmap='RdYlBu_r', transform=ccrs.PlateCarree())
    axes[0].coastlines()
    axes[0].set_title('Forecast')
    plt.colorbar(c1, ax=axes[0], shrink=0.7, label=field_name)

    # Observation
    c2 = axes[1].pcolormesh(lons, lats, obs_field, vmin=vmin, vmax=vmax,
                             cmap='RdYlBu_r', transform=ccrs.PlateCarree())
    axes[1].coastlines()
    axes[1].set_title('Observation')
    plt.colorbar(c2, ax=axes[1], shrink=0.7, label=field_name)

    # Error (forecast - obs)
    error = fcst_field - obs_field
    emax = np.nanmax(np.abs(error))
    c3 = axes[2].pcolormesh(lons, lats, error, vmin=-emax, vmax=emax,
                             cmap='RdBu_r', transform=ccrs.PlateCarree())
    axes[2].coastlines()
    axes[2].set_title('Error (Fcst - Obs)')
    plt.colorbar(c3, ax=axes[2], shrink=0.7, label=field_name)

    plt.tight_layout()
    return fig, axes
```

---

## 10. End-to-End Example: WRF vs. Station Observations

A complete workflow verifying WRF 2-m temperature and 24-h precipitation against surface stations.

```
wrfout_d01 NetCDF  ──→  xr.open_dataset()
                              │
Station CSV obs   ──→  pd.read_csv()
                              │
                    ┌─────────┴──────────┐
                    │  KDTree matching    │
                    │  (fcst → stations)  │
                    └─────────┬──────────┘
                              │
              ┌───────────────┼───────────────┐
              v               v               v
     continuous_metrics  categorical_metrics  brier_score
              │               │               │
              v               v               v
     Taylor diagram    Performance diag   Reliability diag
```

### Step 1: Load Data

```python
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

# Load WRF output
wrf = xr.open_dataset('wrfout_d01_2024-01-15_12:00:00')
t2m_wrf = wrf['T2'].isel(Time=0)                         # K
rain_wrf = (wrf['RAINNC'] + wrf['RAINC']).isel(Time=0)   # mm (accumulated)
lats = wrf['XLAT'].isel(Time=0).values
lons = wrf['XLONG'].isel(Time=0).values

# Load station observations
obs = pd.read_csv('surface_obs_20240115.csv', parse_dates=['datetime'])
obs_12z = obs[obs['datetime'] == '2024-01-15 12:00:00'].copy()
print(f"Loaded {len(obs_12z)} stations")
```

### Step 2: Match Forecast to Stations

```python
# Build KD-tree from WRF grid
grid_points = np.column_stack([lats.ravel(), lons.ravel()])
tree = cKDTree(grid_points)

# Find nearest grid point for each station
station_points = np.column_stack([obs_12z['lat'].values, obs_12z['lon'].values])
dist, idx = tree.query(station_points)

# Extract WRF values at station locations
obs_12z['fcst_t2m'] = t2m_wrf.values.ravel()[idx]
obs_12z['fcst_rain'] = rain_wrf.values.ravel()[idx]

# Drop stations too far from any grid point (e.g., > 0.5 degrees)
obs_12z = obs_12z[dist < 0.5].copy()
print(f"Matched {len(obs_12z)} stations within grid")
```

### Step 3: Continuous Verification (Temperature)

```python
fcst_t = obs_12z['fcst_t2m'].values
obs_t = obs_12z['temperature_K'].values

# Compute metrics
t_metrics = continuous_metrics(fcst_t, obs_t)
t_decomp = mse_decomposition(fcst_t, obs_t)

print("=== 2-m Temperature Verification ===")
print(f"  N     = {t_metrics['N']}")
print(f"  ME    = {t_metrics['ME']:.2f} K")
print(f"  MAE   = {t_metrics['MAE']:.2f} K")
print(f"  RMSE  = {t_metrics['RMSE']:.2f} K")
print(f"  r     = {t_metrics['r']:.3f}")
print(f"  Bias² = {t_decomp['bias_squared']:.3f} K²")
print(f"  Amp.  = {t_decomp['amplitude_error']:.3f} K²")
print(f"  Phase = {t_decomp['phase_error']:.3f} K²")

# Bootstrap 95% CI for RMSE
rmse_func = lambda f, o: np.sqrt(np.mean((f - o)**2))
est, lo, hi = bootstrap_metric(fcst_t, obs_t, rmse_func)
print(f"  RMSE 95% CI: [{lo:.2f}, {hi:.2f}]")
```

### Step 4: Categorical Verification (Precipitation)

```python
fcst_p = obs_12z['fcst_rain'].values
obs_p = obs_12z['precip_mm'].values

# Multi-threshold evaluation
thresholds = [0.1, 1, 5, 10, 25]
precip_results = []
for t in thresholds:
    a, b, c, d = contingency_table(fcst_p, obs_p, t)
    m = categorical_metrics(a, b, c, d)
    m['threshold_mm'] = t
    precip_results.append(m)

df_precip = pd.DataFrame(precip_results)
print("\n=== Precipitation Verification ===")
print(df_precip[['threshold_mm', 'POD', 'FAR', 'CSI', 'ETS', 'FBIAS']].to_string(index=False))
```

### Step 5: Generate Plots

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Taylor diagram for temperature (comparing to a second model if available)
ax_taylor = fig.add_subplot(131, polar=True)
obs_std = np.std(obs_t)
models = [
    {'name': 'WRF', 'std': np.std(fcst_t), 'corrcoef': np.corrcoef(fcst_t, obs_t)[0, 1]}
]
taylor_diagram(models, ref_std=obs_std, ax=ax_taylor)
ax_taylor.set_title('Temperature', pad=20)

# Performance diagram for precipitation
ax_perf = fig.add_subplot(132)
perf_data = {
    'name': 'WRF',
    'SR': [1 - r['FAR'] for r in precip_results if not np.isnan(r['FAR'])],
    'POD': [r['POD'] for r in precip_results if not np.isnan(r['POD'])]
}
performance_diagram([perf_data], ax=ax_perf)
ax_perf.set_title('Precipitation')

# Scatter plot: forecast vs observed temperature
ax_scatter = fig.add_subplot(133)
ax_scatter.scatter(obs_t, fcst_t, alpha=0.5, s=15, edgecolors='none')
lims = [min(obs_t.min(), fcst_t.min()) - 1, max(obs_t.max(), fcst_t.max()) + 1]
ax_scatter.plot(lims, lims, 'k--', alpha=0.5)
ax_scatter.set_xlim(lims)
ax_scatter.set_ylim(lims)
ax_scatter.set_xlabel('Observed T2m (K)')
ax_scatter.set_ylabel('Forecast T2m (K)')
ax_scatter.set_aspect('equal')
ax_scatter.set_title(f'T2m: ME={t_metrics["ME"]:.2f}, RMSE={t_metrics["RMSE"]:.2f}')
ax_scatter.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('wrf_verification_summary.png', dpi=150, bbox_inches='tight')
print("\nSaved: wrf_verification_summary.png")
```

### Step 6: Verification by Lead Time (Multiple Init Times)

```python
# If you have multiple forecast lead times, loop and aggregate
lead_hours = [6, 12, 18, 24, 36, 48, 72]
lead_results = []

for lead in lead_hours:
    # Load forecast for this lead time (pseudo-code — adapt to your file naming)
    # fcst_data = load_wrf_forecast(init_time, lead)
    # obs_data = load_obs(valid_time)
    # matched = match_fcst_to_obs(fcst_data, obs_data)

    # Compute partial sums (for correct aggregation)
    # ps = partial_sums(matched['fcst'], matched['obs'])
    # ps['lead_h'] = lead
    # lead_results.append(ps)
    pass

# Aggregate by lead time
# df_leads = pd.DataFrame(lead_results)
# for lead in lead_hours:
#     subset = df_leads[df_leads['lead_h'] == lead]
#     agg = aggregate_partial_sums(subset.to_dict('records'))
#     print(f"Lead {lead:3d}h: ME={agg['ME']:.2f}, RMSE={agg['RMSE']:.2f}")
```

---

## 11. Tips & Best Practices

### Verification Principles → Python Implementation

| Principle | Python Implementation |
|---|---|
| **Use multiple metrics** | Compute `continuous_metrics()` + `categorical_metrics()` together |
| **Stratify results** | `df.groupby(['season', 'region', 'lead_h']).apply(...)` |
| **Use skill scores** | Pass climatology to `brier_skill_score()`, `crpss()` |
| **Bootstrap CIs** | `bootstrap_metric()` with `n_resamples >= 1000` |
| **Avoid double penalty** | Use `fractions_skill_score()` for precipitation |
| **Aggregate correctly** | Use `partial_sums()` + `aggregate_partial_sums()`, never average RMSE |
| **Homogeneous comparison** | `pd.merge(model1, model2, on='station_id')` to ensure common sample |
| **Verify at appropriate scales** | Use `fss_by_scale()` to find minimum useful scale |

### When to Graduate to MET/METplus

- You need **MODE/MTD** object-based verification
- You need **TC track** verification (TC-Pairs, TC-Stat)
- You're running **operational, production** verification
- You need **METviewer** for interactive database-backed exploration
- Datasets are **very large** and MET's C++ performance helps

### Memory Management for Large Datasets

```python
# Use xarray with dask for lazy loading
ds = xr.open_dataset('large_file.nc', chunks={'time': 10})

# Process in chunks
for t in range(0, len(ds.time), 10):
    chunk = ds.isel(time=slice(t, t + 10)).load()
    # ... compute metrics on chunk ...

# Or use dask directly
import dask.array as da
result = da.from_delayed(...)
```

### Aggregation Pitfall: Never Average RMSE

```python
# WRONG: averaging RMSE across days
daily_rmse = [1.5, 2.0, 1.8, 2.3, 1.6]
wrong_rmse = np.mean(daily_rmse)  # 1.84 — INCORRECT

# CORRECT: aggregate partial sums, then derive RMSE
daily_sums = [partial_sums(f, o) for f, o in daily_pairs]
correct = aggregate_partial_sums(daily_sums)
correct_rmse = correct['RMSE']  # May differ significantly from 1.84
```

---

## Appendix A: Metric Quick Reference

| Metric | Python Function | Theory Notes | MET Equivalent | Range | Perfect |
|---|---|---|---|---|---|
| ME | `continuous_metrics()` | Sec 6.3 | CNT:ME | -inf to +inf | 0 |
| MAE | `continuous_metrics()` | Sec 6.3 | CNT:MAE | 0 to +inf | 0 |
| RMSE | `continuous_metrics()` | Sec 6.3 | CNT:RMSE | 0 to +inf | 0 |
| r | `continuous_metrics()` | Sec 6.3 | CNT:PR_CORR | -1 to 1 | 1 |
| Skill Score | `continuous_metrics()` | Sec 6.3 | CNT:MSE_SS | -inf to 1 | 1 |
| POD | `categorical_metrics()` | Sec 6.1 | CTS:POD | 0 to 1 | 1 |
| FAR | `categorical_metrics()` | Sec 6.1 | CTS:FAR | 0 to 1 | 0 |
| CSI | `categorical_metrics()` | Sec 6.1 | CTS:CSI | 0 to 1 | 1 |
| ETS | `categorical_metrics()` | Sec 6.1 | CTS:GSS | -1/3 to 1 | 1 |
| FBIAS | `categorical_metrics()` | Sec 6.1 | CTS:FBIAS | 0 to +inf | 1 |
| HK | `categorical_metrics()` | Sec 6.1 | CTS:HK | -1 to 1 | 1 |
| HSS | `categorical_metrics()` | Sec 6.1 | CTS:HSS | -1 to 1 | 1 |
| BS | `brier_score()` | Sec 6.4 | PSTD:BRIER | 0 to 1 | 0 |
| BSS | `brier_skill_score()` | Sec 6.4 | PSTD:BSS | -inf to 1 | 1 |
| RPS | `ranked_probability_score()` | Sec 6.4 | RPS:RPS | 0 to 1 | 0 |
| FSS | `fractions_skill_score()` | Sec 7.1 | NBRCTS:FSS | 0 to 1 | 1 |
| CRPS | `crps_ensemble()` | Sec 7.2 | ECNT:CRPS | 0 to +inf | 0 |
| AUC | `compute_roc()` | Sec 6.4 | — | 0 to 1 | 1 |

---

## Appendix B: Useful Resources

| Resource | URL |
|---|---|
| Verification Theory (companion) | [forecast-verification-notes.md](forecast-verification-notes.md) |
| MET Practical Guide (companion) | [met-verification-guide.md](met-verification-guide.md) |
| xarray Documentation | https://docs.xarray.dev/ |
| cfgrib (GRIB2 for xarray) | https://github.com/ecmwf/cfgrib |
| properscoring (CRPS) | https://github.com/properscoring/properscoring |
| scikit-learn Metrics | https://scikit-learn.org/stable/modules/model_evaluation.html |
| cartopy Documentation | https://scitools.org.uk/cartopy/docs/latest/ |
| CAWCR Verification Methods | https://www.cawcr.gov.au/projects/verification/ |
| Wilks (2011) Statistical Methods in the Atmospheric Sciences | (textbook) |
| Jolliffe & Stephenson (2012) Forecast Verification | (textbook) |
