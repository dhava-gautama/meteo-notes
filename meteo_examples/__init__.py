"""meteo_examples — shared helpers for the meteo-notes runnable examples.

This small package backs the Jupyter notebooks under each section's
``notebooks/`` directory. It deliberately mirrors the copy-paste functions in
the written guides (e.g. ``verification/python-verification-guide.md``) so the
notebooks stay short and the guides remain the single source of truth for the
explanations.

Submodules
----------
verif       Forecast-verification metrics (continuous, categorical,
            probabilistic, ensemble, spatial) — mirrors python-verification-guide.md
dataio      Live data download with cached-fallback (IEM ASOS, Open-Meteo,
            NDBC buoys). Notebooks call these so they run online *or* offline.
sampledata  Synthetic-but-realistic NetCDF generators for the model
            post-processing notebooks (WRF / ROMS-CROCO / wave-spectra), which
            can't ship multi-GB real model output.
"""

__version__ = "0.1.0"

from . import verif, dataio, sampledata  # noqa: E402,F401

__all__ = ["verif", "dataio", "sampledata"]
