"""Forecast-verification metrics used by the example notebooks.

These functions mirror the reference implementations in
``verification/python-verification-guide.md``. The guide explains the *why*;
this module is the importable *what* so notebooks read cleanly:

    from meteo_examples import verif
    stats = verif.continuous_metrics(fcst, obs)

If you change a formula here, update the guide to match (and vice-versa).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Continuous metrics
# --------------------------------------------------------------------------- #


def continuous_metrics(fcst, obs):
    """Standard continuous metrics: ME, MAE, RMSE, r, bias, skill score, N."""
    fcst, obs = np.asarray(fcst, dtype=float), np.asarray(obs, dtype=float)
    mask = ~(np.isnan(fcst) | np.isnan(obs))
    f, o = fcst[mask], obs[mask]

    err = f - o
    n = len(f)
    me = np.mean(err)
    mae = np.mean(np.abs(err))
    mse = np.mean(err ** 2)
    rmse = np.sqrt(mse)
    r = np.corrcoef(f, o)[0, 1] if n > 2 else np.nan
    mult_bias = np.mean(f) / np.mean(o) if np.mean(o) != 0 else np.nan

    mse_ref = np.var(o)
    ss = 1 - mse / mse_ref if mse_ref > 0 else np.nan

    return {
        "N": n, "ME": me, "MAE": mae, "RMSE": rmse,
        "r": r, "Mult_Bias": mult_bias, "Skill_Score": ss,
        "FBAR": np.mean(f), "OBAR": np.mean(o),
    }


def partial_sums(fcst, obs):
    """SL1L2 partial sums for correct aggregation (matches MET SL1L2)."""
    f, o = np.asarray(fcst, dtype=float), np.asarray(obs, dtype=float)
    mask = ~(np.isnan(f) | np.isnan(o))
    f, o = f[mask], o[mask]
    return {
        "N": len(f),
        "FBAR": np.mean(f), "OBAR": np.mean(o),
        "FFBAR": np.mean(f ** 2), "FOBAR": np.mean(f * o),
        "OOBAR": np.mean(o ** 2),
    }


def aggregate_partial_sums(sums_list):
    """Aggregate partial-sum dicts, then derive metrics (never average RMSE)."""
    total_n = sum(s["N"] for s in sums_list)
    fbar = sum(s["FBAR"] * s["N"] for s in sums_list) / total_n
    obar = sum(s["OBAR"] * s["N"] for s in sums_list) / total_n
    ffbar = sum(s["FFBAR"] * s["N"] for s in sums_list) / total_n
    fobar = sum(s["FOBAR"] * s["N"] for s in sums_list) / total_n
    oobar = sum(s["OOBAR"] * s["N"] for s in sums_list) / total_n

    me = fbar - obar
    mse = ffbar - 2 * fobar + oobar
    rmse = np.sqrt(mse) if mse >= 0 else np.nan
    return {"N": total_n, "ME": me, "RMSE": rmse, "FBAR": fbar, "OBAR": obar}


# --------------------------------------------------------------------------- #
# Wind
# --------------------------------------------------------------------------- #


def wind_speed(u, v):
    return np.sqrt(np.asarray(u) ** 2 + np.asarray(v) ** 2)


def wind_direction(u, v):
    """Meteorological wind direction (degrees, where wind comes FROM)."""
    return (270 - np.degrees(np.arctan2(v, u))) % 360


def wind_direction_error(fcst_dir, obs_dir):
    """Signed circular difference, wrapped to [-180, 180]."""
    diff = np.asarray(fcst_dir) - np.asarray(obs_dir)
    return (diff + 180) % 360 - 180


def wind_direction_mae(fcst_dir, obs_dir):
    return np.mean(np.abs(wind_direction_error(fcst_dir, obs_dir)))


# --------------------------------------------------------------------------- #
# Categorical
# --------------------------------------------------------------------------- #


def contingency_table(fcst, obs, threshold):
    """2x2 contingency table -> (hits, false_alarms, misses, correct_neg)."""
    f = np.asarray(fcst) >= threshold
    o = np.asarray(obs) >= threshold
    a = int(np.sum(f & o))
    b = int(np.sum(f & ~o))
    c = int(np.sum(~f & o))
    d = int(np.sum(~f & ~o))
    return a, b, c, d


def categorical_metrics(a, b, c, d):
    """POD/FAR/CSI/ETS/FBIAS/HK/HSS/ACC from contingency counts."""
    n = a + b + c + d
    a_ref = (a + b) * (a + c) / n if n > 0 else 0

    pod = a / (a + c) if (a + c) > 0 else np.nan
    far = b / (a + b) if (a + b) > 0 else np.nan
    pofd = b / (b + d) if (b + d) > 0 else np.nan
    sr = a / (a + b) if (a + b) > 0 else np.nan
    csi = a / (a + b + c) if (a + b + c) > 0 else np.nan
    fbias = (a + b) / (a + c) if (a + c) > 0 else np.nan

    denom = a + b + c - a_ref
    ets = (a - a_ref) / denom if denom != 0 else np.nan
    hk = pod - pofd if not (np.isnan(pod) or np.isnan(pofd)) else np.nan
    hss_denom = (a + c) * (c + d) + (a + b) * (b + d)
    hss = 2 * (a * d - b * c) / hss_denom if hss_denom != 0 else np.nan
    acc = (a + d) / n if n > 0 else np.nan

    return {
        "N": n, "Hits": a, "FA": b, "Misses": c, "CN": d,
        "POD": pod, "FAR": far, "POFD": pofd, "SR": sr,
        "CSI": csi, "ETS": ets, "FBIAS": fbias, "HK": hk, "HSS": hss, "ACC": acc,
    }


def multi_threshold_table(fcst, obs, thresholds):
    """Return a DataFrame of categorical metrics across thresholds."""
    rows = []
    for t in thresholds:
        a, b, c, d = contingency_table(fcst, obs, t)
        m = categorical_metrics(a, b, c, d)
        m["threshold"] = t
        rows.append(m)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Probabilistic
# --------------------------------------------------------------------------- #


def brier_score(prob_fcst, obs_binary):
    return np.mean((np.asarray(prob_fcst) - np.asarray(obs_binary)) ** 2)


def brier_skill_score(prob_fcst, obs_binary):
    bs = brier_score(prob_fcst, obs_binary)
    climo = np.mean(obs_binary)
    bs_ref = brier_score(np.full_like(np.asarray(prob_fcst, float), climo), obs_binary)
    return 1 - bs / bs_ref if bs_ref > 0 else np.nan


def reliability_data(prob_fcst, obs_binary, n_bins=10):
    """Bin centers, observed frequency, and counts for a reliability diagram."""
    prob_fcst = np.asarray(prob_fcst, float)
    obs_binary = np.asarray(obs_binary, float)
    edges = np.linspace(0, 1, n_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    obs_freq = np.full(n_bins, np.nan)
    counts = np.zeros(n_bins, dtype=int)
    for i in range(n_bins):
        if i < n_bins - 1:
            mask = (prob_fcst >= edges[i]) & (prob_fcst < edges[i + 1])
        else:
            mask = (prob_fcst >= edges[i]) & (prob_fcst <= edges[i + 1])
        counts[i] = mask.sum()
        if counts[i] > 0:
            obs_freq[i] = obs_binary[mask].mean()
    return centers, obs_freq, counts


# --------------------------------------------------------------------------- #
# Spatial
# --------------------------------------------------------------------------- #


def fractions_skill_score(fcst_field, obs_field, threshold, window_size):
    """Fractions Skill Score (0 = no skill, 1 = perfect)."""
    from scipy.ndimage import uniform_filter

    fb = (np.asarray(fcst_field) >= threshold).astype(float)
    ob = (np.asarray(obs_field) >= threshold).astype(float)
    ff = uniform_filter(fb, size=window_size, mode="constant")
    of = uniform_filter(ob, size=window_size, mode="constant")
    mse = np.mean((ff - of) ** 2)
    mse_ref = np.mean(ff ** 2) + np.mean(of ** 2)
    return 1 - mse / mse_ref if mse_ref > 0 else 1.0


# --------------------------------------------------------------------------- #
# Ensemble
# --------------------------------------------------------------------------- #


def rank_histogram(ensemble, obs, rng=None):
    """Rank-histogram counts (length n_members + 1)."""
    ensemble = np.asarray(ensemble, float)
    obs = np.asarray(obs, float)
    n_times, n_members = ensemble.shape
    rng = np.random.default_rng(0) if rng is None else rng
    ranks = np.zeros(n_members + 1, dtype=int)
    for t in range(n_times):
        rank = int(np.sum(ensemble[t, :] < obs[t]))
        n_ties = int(np.sum(ensemble[t, :] == obs[t]))
        if n_ties > 0:
            rank += int(rng.integers(0, n_ties + 1))
        ranks[rank] += 1
    return ranks


def crps_ensemble(obs, ensemble):
    """Mean CRPS and per-time CRPS via properscoring."""
    import properscoring as ps

    crps = ps.crps_ensemble(np.asarray(obs, float), np.asarray(ensemble, float))
    return float(np.mean(crps)), crps


def spread_skill(ensemble, obs):
    """Return (mean_spread, rmse_of_ensemble_mean)."""
    ensemble = np.asarray(ensemble, float)
    obs = np.asarray(obs, float)
    ens_mean = ensemble.mean(axis=1)
    ens_spread = ensemble.std(axis=1, ddof=1)
    rmse = np.sqrt(np.mean((ens_mean - obs) ** 2))
    return float(ens_spread.mean()), float(rmse)
