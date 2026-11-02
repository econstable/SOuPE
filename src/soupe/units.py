"""Convert energy values using exact SI h, c and Boltzmann constants."""

from __future__ import annotations

import numpy as np

# Joules per unit: cm^-1 is a spectroscopic wavenumber; K means E/k_B.
_JOULES = {
    "meV": 1.602176634e-22,
    "cm^-1": 6.62607015e-34 * 299792458 * 100,
    "K": 1.380649e-23,
}


def convert_energy(value, *, from_unit: str, to_unit: str):
    """Convert a finite real scalar or array between meV, cm^-1 and K.

    Parameters
    ----------
    value : real scalar or array_like
        Energy values or coefficients; negative values are permitted.
    from_unit, to_unit : {'meV', 'cm^-1', 'K'}
        Explicit units. K denotes energy divided by Boltzmann's constant,
        rather than a simulation temperature.

    Returns
    -------
    float or numpy.ndarray
        Converted values with unchanged shape; the input is not modified.

    Raises
    ------
    ValueError
        Units are unknown, or input/output contains nonfinite values.
    TypeError
        Input is not real numeric data.
    """
    if from_unit not in _JOULES or to_unit not in _JOULES:
        raise ValueError("Energy units must be 'meV', 'cm^-1' or 'K'.")
    data = np.asarray(value)
    if data.dtype.kind not in "iuf":
        raise TypeError("Energy values must be real numbers.")
    if not np.all(np.isfinite(data)):
        raise ValueError("Energy values must be finite.")
    with np.errstate(over="ignore", invalid="ignore"):
        result = data.astype(float) * (_JOULES[from_unit] / _JOULES[to_unit])
    if not np.all(np.isfinite(result)):
        raise ValueError("Converted energy values must be finite.")
    return float(result) if result.ndim == 0 else result
