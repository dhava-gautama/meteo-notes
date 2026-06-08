"""Data loaders for the example notebooks: live download + cached fallback.

Every loader follows the same contract so notebooks run **online or offline**:

1. If ``prefer_cache`` (or live download fails) and a cached copy exists in
   ``meteo_examples/data/``, return the cached copy.
2. Otherwise download live from a free, no-auth public source, write the result
   to the cache, and return it.

The small cached CSVs are committed to the repo so the notebooks are fully
reproducible without network access (e.g. in CI or on a plane). Refresh them by
deleting the file in ``meteo_examples/data/`` and re-running with network.

Sources (all free, no API key):
  * IEM ASOS         — surface station observations (Iowa Environmental Mesonet)
  * Open-Meteo       — archived NWP forecasts & ERA5 reanalysis (JSON)
  * NDBC             — moored-buoy marine observations
"""

from __future__ import annotations

import io
import pathlib

import pandas as pd
import requests

DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"
DATA_DIR.mkdir(exist_ok=True)

_TIMEOUT = 60


def _get(url, params=None):
    r = requests.get(url, params=params, timeout=_TIMEOUT,
                     headers={"User-Agent": "meteo-notes-examples/0.1"})
    r.raise_for_status()
    return r


def _cached_csv(cache_name, build_fn, prefer_cache=False, parse_dates=None):
    """Generic cache wrapper for DataFrame-producing loaders.

    build_fn() -> DataFrame is only called when a live fetch is needed.
    """
    path = DATA_DIR / cache_name
    if prefer_cache and path.exists():
        return pd.read_csv(path, parse_dates=parse_dates)
    try:
        df = build_fn()
        df.to_csv(path, index=False)
        return df
    except Exception as exc:  # network down, source moved, etc.
        if path.exists():
            print(f"  [dataio] live fetch failed ({exc!s:.80}); using cached {cache_name}")
            return pd.read_csv(path, parse_dates=parse_dates)
        raise RuntimeError(
            f"Live fetch failed and no cached fallback at {path}. "
            f"Connect to the network once to populate the cache."
        ) from exc


# --------------------------------------------------------------------------- #
# Surface station observations — IEM ASOS
# --------------------------------------------------------------------------- #


def iem_asos(station, start, end, network="ID__ASOS", prefer_cache=False):
    """Hourly ASOS observations for one station from the Iowa Env. Mesonet.

    Parameters
    ----------
    station : str   ICAO id, e.g. "WIII" (Jakarta Soekarno-Hatta).
    start, end : str/date  inclusive UTC date range "YYYY-MM-DD".
    network : str   IEM network code (e.g. "ID__ASOS" for Indonesia,
                    "AWOS"/"<state>_ASOS" for US states).

    Returns columns: valid (datetime, UTC), tmpc, dwpc, relh, sknt, drct, p01m, mslp.
    """
    start, end = pd.Timestamp(start), pd.Timestamp(end)
    cache = f"asos_{station}_{start:%Y%m%d}_{end:%Y%m%d}.csv"

    def build():
        url = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"
        params = {
            "station": station, "network": network,
            "data": ["tmpc", "dwpc", "relh", "sknt", "drct", "p01m", "mslp"],
            "year1": start.year, "month1": start.month, "day1": start.day,
            "year2": end.year, "month2": end.month, "day2": end.day,
            "tz": "UTC", "format": "onlycomma", "missing": "empty",
            "latlon": "no", "report_type": [3, 4],
        }
        text = _get(url, params).text
        df = pd.read_csv(io.StringIO(text))
        df["valid"] = pd.to_datetime(df["valid"])
        for c in ("tmpc", "dwpc", "relh", "sknt", "drct", "p01m", "mslp"):
            if c in df:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        return df

    return _cached_csv(cache, build, prefer_cache=prefer_cache, parse_dates=["valid"])


# --------------------------------------------------------------------------- #
# Archived NWP forecast / reanalysis — Open-Meteo
# --------------------------------------------------------------------------- #


def open_meteo_forecast(lat, lon, start, end, hourly="temperature_2m",
                        prefer_cache=False):
    """Archived deterministic NWP forecast (Open-Meteo Historical-Forecast API).

    Returns the *forecast that was actually issued* for past dates — ideal as
    the "model" side of a forecast-vs-obs comparison. Default model is the
    best-available global model at the requested location.

    Returns columns: time (datetime, UTC) + one column per requested variable.
    """
    cache = f"openmeteo_fc_{lat:.2f}_{lon:.2f}_{pd.Timestamp(start):%Y%m%d}_" \
            f"{pd.Timestamp(end):%Y%m%d}_{hourly.replace(',', '-')}.csv"

    def build():
        url = "https://historical-forecast-api.open-meteo.com/v1/forecast"
        params = {
            "latitude": lat, "longitude": lon,
            "start_date": str(pd.Timestamp(start).date()),
            "end_date": str(pd.Timestamp(end).date()),
            "hourly": hourly, "timezone": "UTC",
        }
        j = _get(url, params).json()["hourly"]
        df = pd.DataFrame(j)
        df["time"] = pd.to_datetime(df["time"])
        return df

    return _cached_csv(cache, build, prefer_cache=prefer_cache, parse_dates=["time"])


def open_meteo_ensemble(lat, lon, start, end, variable="temperature_2m",
                        model="gfs_seamless", prefer_cache=False):
    """Ensemble forecast (Open-Meteo Ensemble API) -> tidy long DataFrame.

    Returns columns: time, member, value. Members come from the requested
    ensemble system (default GEFS via "gfs_seamless").
    """
    cache = f"openmeteo_ens_{lat:.2f}_{lon:.2f}_{pd.Timestamp(start):%Y%m%d}_" \
            f"{pd.Timestamp(end):%Y%m%d}_{variable}_{model}.csv"

    def build():
        url = "https://ensemble-api.open-meteo.com/v1/ensemble"
        params = {
            "latitude": lat, "longitude": lon,
            "start_date": str(pd.Timestamp(start).date()),
            "end_date": str(pd.Timestamp(end).date()),
            "hourly": variable, "models": model, "timezone": "UTC",
        }
        hourly = _get(url, params).json()["hourly"]
        time = pd.to_datetime(hourly["time"])
        rows = []
        for key, vals in hourly.items():
            if not key.startswith(variable):
                continue
            # key is "temperature_2m" or "temperature_2m_memberNN"
            member = 0 if key == variable else int(key.rsplit("member", 1)[-1])
            rows.append(pd.DataFrame({"time": time, "member": member, "value": vals}))
        return pd.concat(rows, ignore_index=True)

    return _cached_csv(cache, build, prefer_cache=prefer_cache, parse_dates=["time"])


# --------------------------------------------------------------------------- #
# Marine observations — NDBC moored buoys
# --------------------------------------------------------------------------- #


def ndbc_realtime(buoy, prefer_cache=False):
    """Recent standard-meteorological observations for an NDBC buoy.

    Returns columns include: time (datetime, UTC), WSPD, WVHT (sig. wave
    height, m), DPD, WTMP (sea temp), PRES, ATMP. Realtime feed = last ~45 days.
    """
    cache = f"ndbc_{buoy}.csv"

    def build():
        url = f"https://www.ndbc.noaa.gov/data/realtime2/{buoy}.txt"
        text = _get(url).text
        lines = text.splitlines()
        # Header line is like "#YY  MM DD hh mm WDIR ..."; drop the leading '#'.
        names = lines[0].lstrip("#").split()
        # line 1 is units (starts with '#'); data starts at line 2
        df = pd.read_csv(io.StringIO("\n".join(lines[2:])),
                         sep=r"\s+", names=names, na_values="MM")
        df["time"] = pd.to_datetime(
            df[["YY", "MM", "DD", "hh", "mm"]].rename(
                columns={"YY": "year", "MM": "month", "DD": "day",
                         "hh": "hour", "mm": "minute"}), utc=True).dt.tz_localize(None)
        return df.sort_values("time").reset_index(drop=True)

    return _cached_csv(cache, build, prefer_cache=prefer_cache, parse_dates=["time"])
