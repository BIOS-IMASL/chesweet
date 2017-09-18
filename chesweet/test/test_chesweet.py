#from numpy.testing import assert_almost_equal
import numpy as np
from ..chesweet import _round_down_up, _nearest_chi

torsionals = [-10, 10., 180., -180., 60., 75., -185, 185]

def test_round_down_up():
    ref = [(-10, 0), (10, 20), (170, 180), (-180, -170), (60, 70), (70, 80),
           (170.0, 180.0), (-180.0, -170.0)]

    for idx, t in enumerate(torsionals):
        assert _round_down_up(t, 10) == ref[idx]
