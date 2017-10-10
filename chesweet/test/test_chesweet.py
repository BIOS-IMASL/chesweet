#from numpy.testing import assert_almost_equal
import numpy as np
from ..chesweet import *
from ..chesweet import _round_down_up, _nearest_chi

torsionals = [-10, 10., 180., -180., 60., 75., -185, 185]

def test_round_down_up():
    ref = [(-10, 0), (10, 20), (170, 180), (-180, -170), (60, 70), (70, 80),
           (170.0, 180.0), (-180.0, -170.0)]

    for idx, t in enumerate(torsionals):
        assert _round_down_up(t, 10) == ref[idx]

def test_compute_cs_all():
    che = CheSweet(path='../lut')
    che.compute_cs('a-D-Glcp-1-4-a-D-Glcp', 50, 60)

def test_compute_cs_few():
    che = CheSweet(path='../lut', disaccharides=['a-D-Glcp-1-1-a-D-Glcp','a-D-Glcp-1-4-a-D-Glcp','b-D-Galp-1-6-b-D-Galp'])
    che.compute_cs('a-D-Glcp-1-1-a-D-Glcp', 50, 60)
    che.compute_cs('a-D-Glcp-1-4-a-D-Glcp', 50, 60)
    che.compute_cs('b-D-Galp-1-6-b-D-Galp', 50, 60, 60)

def test_compute_cs_one():
    che = CheSweet(path='../lut', disaccharides=['b-D-GlcpNAc-1-4-b-D-Glcp1Me2NAc'])
    che.compute_cs('b-D-GlcpNAc-1-4-b-D-Glcp1Me2NAc', 50, 60)
