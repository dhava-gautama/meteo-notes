"""Readers for native model-output text formats used by the example notebooks.

So far: SWAN ASCII spectral files (``SPECOUT ... SPEC2D``). These are the formats
that don't already have an xarray backend, so a small parser keeps the notebooks
clean and the code reusable on your own model output.
"""

from __future__ import annotations

import pathlib

import numpy as np

_DATA = pathlib.Path(__file__).resolve().parent / "data"


def swan_ideal_spec_path():
    """Path to the committed *real* SWAN 2-D spectrum from the idealized run."""
    return _DATA / "swan_ideal_P1.spc"


def swan_ideal_fetch_path():
    """Path to the committed centerline Hs/Tp-vs-fetch table from the SWAN run."""
    return _DATA / "swan_ideal_fetch.csv"


def read_swan_spec(path):
    """Parse a SWAN 2-D spectral file (``SPEC2D``) for a single location.

    Returns ``(freq, theta_rad, E2d, meta)`` where:
      * ``freq``      : 1-D frequencies (Hz)
      * ``theta_rad`` : 1-D nautical directions (radians), as written in the file
      * ``E2d``       : energy density (nfreq, ndir) in m²/Hz/deg
      * ``meta``      : dict with ``hs`` (integrated significant wave height, m),
                        ``tp`` (peak period, s), ``mean_dir_deg``.

    SWAN writes equal-width direction bins, so integrate over direction with a
    constant-bin sum (a periodic quantity — trapezoid over the sorted range would
    drop the wrap segment).
    """
    lines = pathlib.Path(path).read_text().splitlines()

    def _block(tag):
        i = next(k for k, l in enumerate(lines) if l.strip().startswith(tag))
        n = int(lines[i + 1].split()[0])
        return np.array([float(lines[i + 2 + j].split()[0]) for j in range(n)])

    freq = _block("AFREQ")
    dirs = _block("NDIR")                      # degrees, nautical
    fi = next(k for k, l in enumerate(lines) if l.strip() == "FACTOR")
    factor = float(lines[fi + 1])
    rows = [[float(p) for p in l.split()] for l in lines[fi + 2:]
            if len(l.split()) == len(dirs)]
    E2d = np.array(rows) * factor              # (nfreq, ndir), m²/Hz/deg
    E2d[E2d < 0] = 0.0                          # SWAN exception value (-99)

    dtheta = abs(dirs[1] - dirs[0])
    Ef = E2d.sum(axis=1) * dtheta              # 1-D spectrum, m²/Hz
    m0 = np.trapezoid(Ef, freq) if hasattr(np, "trapezoid") else np.trapz(Ef, freq)
    fp = freq[Ef.argmax()]
    a1 = (E2d * np.cos(np.deg2rad(dirs))[None, :]).sum()
    b1 = (E2d * np.sin(np.deg2rad(dirs))[None, :]).sum()
    meta = {
        "hs": float(4 * np.sqrt(m0)),
        "tp": float(1.0 / fp),
        "mean_dir_deg": float(np.rad2deg(np.arctan2(b1, a1)) % 360),
    }
    return freq, np.deg2rad(dirs), E2d, meta
