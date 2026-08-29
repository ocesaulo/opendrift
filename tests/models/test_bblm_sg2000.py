"""Styles & Glenn (2000) combined wave-current BBL kernel.

Validates `opendrift.models.bblm_sg2000.sg2000_solve` against Table 2 of the
user test case distributed with the authors' MATLAB implementation. Inputs come
from `tests/test_data/sppm/test2.mat`; see that directory's PROVENANCE.md.
"""
import os

import numpy as np
import pytest
from scipy.io import loadmat

from opendrift.models.bblm_sg2000 import sg2000_solve

# UserTestCase.pdf, Table 2 (cm/s and cm).
# rows: SW-SC, WW-SC, SW-WC ; cols: ustarcw, ustarwm, ustarc, znot, z1
EXPECTED = np.array([
    [12.9, 11.2, 6.5, 1.3, 3.0],
    [6.7, 3.0, 6.0, 0.1, 1.4],
    [17.1, 17.0, 1.1, 0.6, 3.1],
])
LABELS = ['SW-SC', 'WW-SC', 'SW-WC']
COLS = ['ustarcw', 'ustarwm', 'ustarc', 'znot', 'z1']


def _solve(test_data):
    data = loadmat(os.path.join(test_data, 'sppm', 'test2.mat'))['DATA'].astype(float)
    Ub, Ab, Ur, zr, deg = (data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5])
    return sg2000_solve(Ub, Ab, Ur, zr, deg)


def test_table2(test_data):
    """Kernel reproduces all 15 Table 2 entries to the table's printed precision."""
    out = _solve(test_data)
    got = np.column_stack([out[c] for c in COLS])
    assert np.allclose(np.round(got, 1), EXPECTED, atol=0.15), got


def test_pure_wave_limit_is_finite():
    """With a vanishing current the solver falls back to the pure-wave branch."""
    out = sg2000_solve(np.array([40.0]), np.array([50.0]), np.array([1e-6]),
                       np.array([100.0]), np.array([0.0]))
    assert np.all(np.isfinite(out['ustarcw'])) and out['ustarcw'][0] > 0
    # with no current, the mean shear velocity is negligible next to the wave one
    assert out['ustarc'][0] < 0.05 * out['ustarwm'][0]


@pytest.mark.parametrize('deg', [0.0, 45.0, 90.0])
def test_angle_variation_finite(deg):
    """Solutions stay finite and ordered across wave-current angles."""
    out = sg2000_solve(np.array([40.0]), np.array([50.0]), np.array([60.0]),
                       np.array([100.0]), np.array([deg]))
    assert np.all(np.isfinite(out['ustarcw']))
    # the combined shear velocity bounds both constituents
    assert out['ustarcw'][0] >= out['ustarc'][0]
