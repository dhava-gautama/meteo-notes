"""Synthetic-but-realistic model output for the post-processing notebooks.

Real WRF/ROMS/CROCO/SWAN output is multi-gigabyte and can't live in a notes
repo, but the *post-processing* (open file -> select variable -> plot) is
identical whether the data is real or synthetic. These generators build small
NetCDF files that look like the real thing — same variable names, dimensions,
and attributes — so the notebook code transfers verbatim to real output.

Everything is deterministic (fixed seed): re-running gives identical files.
The domains are placed over Indonesian waters to match the guides.
"""

from __future__ import annotations

import pathlib

import numpy as np
import xarray as xr

_DATA = pathlib.Path(__file__).resolve().parent / "data"
_RNG = np.random.default_rng(20240115)


def wrf_idealized_path():
    """Path to the committed *real* WRF output from the em_quarter_ss run.

    A compact extract (column-max updraft, accumulated precip, a W
    cross-section; 6 frames) from an actual WRFV4.7.1 ``em_quarter_ss``
    idealized supercell run — see the WRF guide's idealized-run section and
    ``models/notebooks/wrf_idealized_supercell.ipynb``.
    """
    return _DATA / "wrf_em_quarter_ss.nc"


def roms_upwelling_path():
    """Path to the committed *real* ROMS output from the UPWELLING run.

    A compact extract (SST, a cross-channel temperature section, bathymetry,
    21 frames) from an actual ROMS 4.3 ``UPWELLING`` idealized run — see the
    ROMS guide's idealized-run section and
    ``models/notebooks/roms_upwelling.ipynb``.
    """
    return _DATA / "roms_upwelling.nc"


def croco_basin_path():
    """Path to the committed *real* CROCO output from the BASIN run.

    A compact extract (free-surface SSH for 37 frames, surface temperature and
    depth-averaged currents) from an actual CROCO BASIN idealized run (the
    wind-driven double gyre) — see the CROCO guide's idealized-run section and
    ``models/notebooks/croco_basin.ipynb``.
    """
    return _DATA / "croco_basin.nc"


def ww3_propagation_path():
    """Path to the committed *real* WAVEWATCH III propagation-test output.

    A compact extract (Hs for 5 frames, mean period and direction) from an
    actual WW3 ``ww3_tp2.2`` 2-D propagation regression test — see the WW3
    guide's idealized-run section and ``models/notebooks/ww3_propagation.ipynb``.
    """
    return _DATA / "ww3_propagation.nc"

# np.trapz was renamed to np.trapezoid in NumPy 2.0; support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz


def _smooth(field, passes=4):
    """Cheap 3x3 box smoothing to make random fields look geophysical."""
    f = field.copy()
    for _ in range(passes):
        f = (f
             + np.roll(f, 1, 0) + np.roll(f, -1, 0)
             + np.roll(f, 1, 1) + np.roll(f, -1, 1)) / 5.0
    return f


def make_wrfout(path, nx=120, ny=100, ntime=4):
    """Write a WRF-like wrfout NetCDF over the Sunda Strait / West Java.

    Variables: T2, RAINNC, RAINC, U10, V10, XLAT, XLONG, Times — matching
    real wrfout naming so ``wrf-python``/``xarray`` code is identical.
    """
    lon = np.linspace(104.0, 109.0, nx)
    lat = np.linspace(-8.5, -4.5, ny)
    lon2d, lat2d = np.meshgrid(lon, lat)

    times = np.array([f"2024-01-15_{h:02d}:00:00" for h in (0, 6, 12, 18)][:ntime],
                     dtype="S19")

    # Diurnal 2-m temperature with a land/sea-ish gradient and noise.
    base = 300.0 - 0.6 * (lat2d + 6.5) ** 2 / 4.0
    t2 = np.empty((ntime, ny, nx), np.float32)
    rain_nc = np.zeros((ntime, ny, nx), np.float32)
    rain_c = np.zeros((ntime, ny, nx), np.float32)
    u10 = np.empty((ntime, ny, nx), np.float32)
    v10 = np.empty((ntime, ny, nx), np.float32)

    accum = np.zeros((ny, nx))
    for t in range(ntime):
        hour = (0, 6, 12, 18)[t]
        diurnal = 4.0 * np.cos((hour - 14) / 24 * 2 * np.pi)
        t2[t] = (base + diurnal + _smooth(_RNG.normal(0, 1.5, (ny, nx)))).astype(np.float32)
        # A convective blob that grows and drifts (typical MC afternoon storm).
        cx, cy = 106.0 + 0.2 * t, -6.5 + 0.1 * t
        blob = 30.0 * np.exp(-(((lon2d - cx) / 0.5) ** 2 + ((lat2d - cy) / 0.4) ** 2))
        accum = accum + blob * (hour >= 6)
        rain_nc[t] = _smooth(accum * 0.8).astype(np.float32)
        rain_c[t] = _smooth(accum * 0.2).astype(np.float32)
        # Monsoonal westerlies + noise.
        u10[t] = (5.0 + _smooth(_RNG.normal(0, 1.5, (ny, nx)))).astype(np.float32)
        v10[t] = (-2.0 + _smooth(_RNG.normal(0, 1.5, (ny, nx)))).astype(np.float32)

    ds = xr.Dataset(
        data_vars=dict(
            Times=(("Time",), times),
            T2=(("Time", "south_north", "west_east"), t2,
                {"description": "TEMP at 2 M", "units": "K"}),
            RAINNC=(("Time", "south_north", "west_east"), rain_nc,
                    {"description": "ACCUMULATED TOTAL GRID SCALE PRECIPITATION", "units": "mm"}),
            RAINC=(("Time", "south_north", "west_east"), rain_c,
                   {"description": "ACCUMULATED TOTAL CUMULUS PRECIPITATION", "units": "mm"}),
            U10=(("Time", "south_north", "west_east"), u10,
                 {"description": "U at 10 M", "units": "m s-1"}),
            V10=(("Time", "south_north", "west_east"), v10,
                 {"description": "V at 10 M", "units": "m s-1"}),
            XLAT=(("south_north", "west_east"), lat2d.astype(np.float32),
                  {"units": "degree_north"}),
            XLONG=(("south_north", "west_east"), lon2d.astype(np.float32),
                   {"units": "degree_east"}),
        ),
        attrs={"TITLE": "SYNTHETIC wrfout for meteo-notes examples",
               "GRID_ID": 1, "DX": 5000.0, "DY": 5000.0},
    )
    ds.to_netcdf(path)
    return path


def make_roms_history(path, xi=100, eta=90, s=20, ntime=3):
    """Write a ROMS/CROCO-like history file over the Indonesian Throughflow.

    Variables: temp (s_rho), zeta, u/v surface, lat_rho, lon_rho, s_rho,
    ocean_time — matching ROMS CF-ish naming for transferable post-processing.
    """
    lon = np.linspace(114.0, 120.0, xi)   # Makassar / Lombok Strait region
    lat = np.linspace(-9.0, -3.0, eta)
    lon2d, lat2d = np.meshgrid(lon, lat)
    s_rho = np.linspace(-0.975, -0.025, s)  # stretched vertical coordinate

    # SST: warm tropical surface with a cool throughflow tongue down a strait.
    sst = 29.5 - 1.5 * np.exp(-((lon2d - 117.0) / 0.6) ** 2) + _smooth(
        _RNG.normal(0, 0.2, (eta, xi)))

    temp = np.empty((ntime, s, eta, xi), np.float32)
    zeta = np.empty((ntime, eta, xi), np.float32)
    for t in range(ntime):
        drift = 0.1 * t
        surf = (sst + drift).astype(np.float32)
        # Thermocline: temperature decreasing with depth (s_rho from -1..0).
        for k, z in enumerate(s_rho):
            temp[t, k] = (surf - (1 - (z + 1)) * 12.0).astype(np.float32)
        zeta[t] = _smooth(_RNG.normal(0, 0.05, (eta, xi))).astype(np.float32)

    ds = xr.Dataset(
        data_vars=dict(
            temp=(("ocean_time", "s_rho", "eta_rho", "xi_rho"), temp,
                  {"long_name": "potential temperature", "units": "Celsius"}),
            zeta=(("ocean_time", "eta_rho", "xi_rho"), zeta,
                  {"long_name": "free-surface", "units": "meter"}),
            lon_rho=(("eta_rho", "xi_rho"), lon2d.astype(np.float32),
                     {"units": "degree_east"}),
            lat_rho=(("eta_rho", "xi_rho"), lat2d.astype(np.float32),
                     {"units": "degree_north"}),
            s_rho=(("s_rho",), s_rho.astype(np.float32),
                   {"long_name": "S-coordinate at RHO-points"}),
        ),
        coords={"ocean_time": np.arange(ntime) * 86400.0},
        attrs={"title": "SYNTHETIC ROMS history for meteo-notes examples",
               "type": "ROMS/CROCO history file"},
    )
    ds.to_netcdf(path)
    return path


def make_wave_spectrum(nfreq=30, ndir=36):
    """Return (freq, theta, E2d, Hs) for a JONSWAP-like directional spectrum.

    E2d has shape (nfreq, ndir) in m^2/Hz/rad; Hs is the integrated significant
    wave height. Mirrors what you'd read from SWAN/WW3 spectral output.
    """
    freq = np.linspace(0.04, 0.4, nfreq)          # Hz
    theta = np.linspace(0, 2 * np.pi, ndir, endpoint=False)

    # JONSWAP 1-D spectrum.
    fp, hs_target, gamma = 0.09, 2.5, 3.3
    alpha = 0.0081
    g = 9.81
    S = alpha * g ** 2 / (2 * np.pi) ** 4 / freq ** 5 * np.exp(-1.25 * (fp / freq) ** 4)
    sigma = np.where(freq <= fp, 0.07, 0.09)
    S *= gamma ** np.exp(-((freq - fp) ** 2) / (2 * sigma ** 2 * fp ** 2))

    # cos^2 directional spreading about a mean wave-from direction.
    mean_dir = np.deg2rad(225)   # swell from the SW
    D = np.maximum(np.cos((theta - mean_dir) / 2) ** 4, 0)
    D /= _trapz(D, theta)

    E2d = S[:, None] * D[None, :]
    # Scale so integrated Hs matches the target.
    m0 = _trapz(_trapz(E2d, theta, axis=1), freq)
    E2d *= (hs_target ** 2 / (16 * m0))
    m0 = _trapz(_trapz(E2d, theta, axis=1), freq)
    hs = 4 * np.sqrt(m0)
    return freq, theta, E2d, float(hs)
