"""Regression tests for the meteo_examples helper package.

These are fast, offline, and deterministic — no network. They guard the metric
formulas (against hand-computed values), the native-format readers, and the
shapes/contents of the committed sample data the notebooks depend on.

Run:  pytest -q
"""

import numpy as np
import pytest

from meteo_examples import verif, sampledata, readers


# --------------------------------------------------------------------------- #
# Continuous metrics — checked against values computed by hand
# --------------------------------------------------------------------------- #


def test_continuous_metrics_known_values():
    fcst = np.array([1.0, 2.0, 3.0, 4.0])
    obs = np.array([1.0, 2.0, 3.0, 5.0])          # only last differs by -1
    m = verif.continuous_metrics(fcst, obs)
    assert m["N"] == 4
    assert m["ME"] == pytest.approx(-0.25)        # mean(0,0,0,-1)
    assert m["MAE"] == pytest.approx(0.25)
    assert m["RMSE"] == pytest.approx(0.5)        # sqrt(1/4)


def test_continuous_metrics_ignores_nan():
    f = np.array([1.0, np.nan, 3.0])
    o = np.array([1.0, 2.0, np.nan])
    m = verif.continuous_metrics(f, o)
    assert m["N"] == 1                            # only the first pair survives


def test_partial_sums_aggregation_matches_direct_rmse():
    rng = np.random.default_rng(0)
    f, o = rng.normal(size=200), rng.normal(size=200)
    direct = verif.continuous_metrics(f, o)["RMSE"]
    # split into two halves, aggregate partial sums, derive RMSE
    ps = [verif.partial_sums(f[:100], o[:100]), verif.partial_sums(f[100:], o[100:])]
    agg = verif.aggregate_partial_sums(ps)["RMSE"]
    assert agg == pytest.approx(direct, rel=1e-9)


# --------------------------------------------------------------------------- #
# Categorical
# --------------------------------------------------------------------------- #


def test_contingency_and_categorical_metrics():
    fcst = np.array([1.0, 1.0, 0.0, 0.0, 1.0])
    obs = np.array([1.0, 0.0, 0.0, 1.0, 1.0])
    a, b, c, d = verif.contingency_table(fcst, obs, threshold=0.5)
    assert (a, b, c, d) == (2, 1, 1, 1)           # hits, FA, miss, CN
    m = verif.categorical_metrics(a, b, c, d)
    assert m["POD"] == pytest.approx(2 / 3)
    assert m["FAR"] == pytest.approx(1 / 3)
    assert m["CSI"] == pytest.approx(2 / 4)


# --------------------------------------------------------------------------- #
# Wind (circular) and waves
# --------------------------------------------------------------------------- #


def test_wind_direction_mae_wraps():
    # 350 vs 10 deg should be 20 deg apart, not 340
    assert verif.wind_direction_mae([350.0], [10.0]) == pytest.approx(20.0)


def test_scatter_index_perfect_is_zero():
    o = np.array([1.0, 2.0, 3.0])
    assert verif.scatter_index(o, o) == pytest.approx(0.0)


def test_scatter_index_removes_bias():
    o = np.array([2.0, 2.0, 2.0, 2.0])
    f = o + 0.5                                   # pure bias, no scatter
    assert verif.scatter_index(f, o) == pytest.approx(0.0, abs=1e-9)


# --------------------------------------------------------------------------- #
# Probabilistic & ensemble
# --------------------------------------------------------------------------- #


def test_brier_score_and_skill():
    p = np.array([1.0, 0.0, 1.0, 0.0])
    o = np.array([1.0, 0.0, 1.0, 0.0])
    assert verif.brier_score(p, o) == pytest.approx(0.0)
    assert verif.brier_skill_score(p, o) == pytest.approx(1.0)


def test_rank_histogram_sums_to_ntimes():
    rng = np.random.default_rng(1)
    ens = rng.normal(size=(50, 10))
    obs = rng.normal(size=50)
    ranks = verif.rank_histogram(ens, obs)
    assert len(ranks) == 11                       # n_members + 1
    assert ranks.sum() == 50


def test_fss_perfect_field_is_one():
    field = np.zeros((20, 20))
    field[5:15, 5:15] = 10.0
    assert verif.fractions_skill_score(field, field, threshold=1.0, window_size=3) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# Native-format readers and committed sample data
# --------------------------------------------------------------------------- #


def test_swan_spec_reader_matches_table_hs():
    freq, theta, E2d, meta = readers.read_swan_spec(readers.swan_ideal_spec_path())
    assert E2d.shape == (len(freq), len(theta))
    # integrated Hs should match the SWAN-reported ~1.76 m at that point
    assert meta["hs"] == pytest.approx(1.76, abs=0.1)
    assert 4.0 < meta["tp"] < 7.0


@pytest.mark.parametrize("loader,must_have", [
    (sampledata.wrf_idealized_path, ("max_updraft", "rain")),
    (sampledata.roms_upwelling_path, ("sst", "temp_section")),
    (sampledata.croco_basin_path, ("zeta", "sst")),
    (sampledata.ww3_propagation_path, ("hs",)),
    (sampledata.schism_seiche_path, ("elevation", "node_x", "tri")),
])
def test_committed_samples_open_and_have_vars(loader, must_have):
    import xarray as xr
    ds = xr.open_dataset(loader())
    for v in must_have:
        assert v in ds.variables, f"{v} missing from {loader.__name__}"
        assert np.isfinite(np.asarray(ds[v])).any()


def test_synthetic_generators_run(tmp_path):
    wrf = sampledata.make_wrfout(tmp_path / "w.nc")
    roms = sampledata.make_roms_history(tmp_path / "r.nc")
    freq, theta, E2d, hs = sampledata.make_wave_spectrum()
    assert wrf.exists() and roms.exists()
    assert hs == pytest.approx(2.5, abs=0.3)     # target Hs of the JONSWAP spectrum
