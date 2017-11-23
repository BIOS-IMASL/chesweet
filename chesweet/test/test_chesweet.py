#from numpy.testing import assert_almost_equal
import numpy as np
from ..chesweet import *
from ..chesweet import _round_down_up, _nearest_chi

torsionals = [-10, 10., 180., -180., 60., 75., -185, 185]

ef_corr = 183.4
disaccharides_list = ['a-D-Glcp-1-1-a-D-Glcp', 'a-D-Galp-1-3-b-D-Galp', 'b-D-Galp-1-6-b-D-Galp']

disaccharides_red = CheSweet(path='chesweet/lut', disaccharides = disaccharides_list)

disaccharides_full = CheSweet(path='chesweet/lut', disaccharides = disaccharides_list, full=True)

### Data to test compute_cs function

full_tors_1_3 = [(120, -120, 180, -60, 60), (-80, -130, 50, 50, 170), (0, 0, 60, 60, 60),
                 (105.7, 144.3, 65.3, 160.1, -45.6), (-55.5, -105.6, 78.9, -65.8, 46.79), 
                 (171.5, 65.2, 36.4, 58.9, 125.8)]
full_chem_1_3 = [(72.6629, 94.3636),(84.2208, 102.9435), (np.inf, np.inf), 
                 (84.18086, 101.497015), (84.6804, 112.0608), 
                 (np.inf, np.inf)]

red_tors_1_3 = [(120, -120), (-80, -130), (0, 0), (105.7, 144.3), 
                (-55.5, -105.6), (171.5, 65.2)]
red_chem_1_3 = [(71.9521, 91.3869),(82.612, 103.2848), (np.inf, np.inf), (84.386021, 103.482665),
                (84.4813, 111.5404), (np.inf, np.inf)]

full_tors_1_1 = [(-70, 100, 180, -60), (-105, 55, 60, 60), (-60, -60, 180, 64), 
                 (88.5, 154.3, 123.4, 70.5), (146.5, 5.8, 70.6, 60.3), (175.8, -15.4, 46.7, 58.9)]
full_chem_1_1 = [(82.9463, 87.7168),(82.4441, 90.0294), (np.inf, np.inf),
                 (78.73016, 75.417343), (90.4849, 84.0846), (np.inf, np.inf)]

red_tors_1_1 = [(-70, 120), (-105, 55), (-60, -60), 
                (88.5, 154.3), (146.5, 5.8), (175.8, -15.4)]
red_chem_1_1 = [(82.9808, 84.6705),(81.7204, 89.8081), (np.inf, np.inf), 
                (78.521635, 73.952066), (89.822, 83.7343), (np.inf, np.inf)]

full_tors_1_6 = [(-70, 120, 180, -60, -60), (135, 65, 60, 180, 60), (0, -10, 60, 60, 60), 
                 (-126.4, 46.7, 46.8, 120.6, 146.7), (-155.7, -76.4, 160.8, 37.8, -25.4),
                 (-175.89, -26.7, 78.6, -98.1, 64.5)]
full_chem_1_6 = [(72.9099, 113.8796), (76.7441, 112.9245), (np.inf, np.inf), 
                 (77.953199, 109.076416), (72.6276, 112.8706), 
                 (np.inf, np.inf)]

red_tors_1_6 = [(-70, 120, 180), (135, 65, 60), (0, -10, 60), 
                 (-126.4, 46.7, 46.8), (-155.7, -76.4, 160.8), (-175.89, -26.7, 78.6)]
red_chem_1_6 = [(72.4989, 113.5847),(77.6637, 113.626), (np.inf, np.inf), 
                (79.505894, 109.841949), (71.4614, 112.203), (np.inf, np.inf)]

all_tors_full = [full_tors_1_1, full_tors_1_3, full_tors_1_6]
all_chem_full = [full_chem_1_1, full_chem_1_3, full_chem_1_6]

all_tors_red = [red_tors_1_1, red_tors_1_3, red_tors_1_6]
all_chem_red = [red_chem_1_1, red_chem_1_3, red_chem_1_6]

### Data to test compute_tors function

chem_test_full_1_3 = [(72.6629, 94.3636), (20, 150), (72.6629, 150),
                      (20, 94.3636)]
tors_test_full_1_3 = [np.array([[ 100., -120.,  -60.,  -60.,   60.],
                                [ 100., -100.,  -60.,  -60.,   60.],
                                [ 100., -100.,  -60.,  -60.,  180.],
                                [ 110., -130.,  -60.,  -60.,   60.],
                                [ 110., -120.,  180.,  -60.,   60.],
                                [ 110., -120.,  180.,  -60.,  180.],
                                [ 110., -110.,  180.,  -60.,   60.],
                                [ 110., -110.,  180.,  -60.,  180.],
                                [ 110., -100.,  -60.,  -60.,   60.],
                                [ 110., -100.,  -60.,  -60.,  180.],
                                [ 120., -130.,  180.,  -60.,   60.],
                                [ 120., -120.,  180.,  -60.,   60.],
                                [ 120., -100.,  -60.,  -60.,   60.],
                                [ 140., -120.,  180.,  180.,  -60.]]),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5))
                     ]# np.zeros only used to completed the list

chem_test_full_1_1 = [(82.4441, 90.0294), (20, 150), (82.4441, 150),
                      (20, 90.0294)]
tors_test_full_1_1 = [np.array([[-120.,  100.,   60.,  180.],
                                [-110.,   90.,   60.,  180.],
                                [-110.,  130.,  180.,   60.],
                                [-100.,   60.,   60.,   60.],
                                [-100.,   60.,   60.,  180.],
                                [-100.,  110.,  180.,   60.],
                                [ -90.,   80.,   60.,  180.],
                                [ -90.,   90.,  180.,  180.],
                                [ -80.,   90.,   60.,   60.],
                                [ -80.,   90.,   60.,  -60.],
                                [ -70.,   80.,  180.,  -60.],
                                [ -70.,   90.,  180.,   60.]]),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5))
                     ]

chem_test_full_1_6 = [(68.7441, 112.9245), (20, 150), (68.7441, 150),
                      (20, 112.9245)]
tors_test_full_1_6 = [np.array([[  40.,  -20.,   60.,  180.,  -60.],
                                [-150.,  110.,  180.,  180.,  -60.],
                                [-140.,  -80.,  180.,  -60.,  -60.],
                                [ -80., -160.,  180.,  180.,  -60.],
                                [ -80., -120.,  180.,  180.,   60.],
                                [ -80., -120.,  180.,  180.,  180.],
                                [ -80., -120.,  180.,  180.,  -60.]]),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5))
                     ]

chem_test_red_1_3 = [(71.9521, 91.3869), (20, 150), (71.9521, 150),
                     (20, 91.3869)]
tors_test_red_1_3 = [np.array([[ 110., -120.],
                               [ 110., -110.],
                               [ 120., -130.],
                               [ 120., -120.]]),
                     np.zeros(shape=(0, 5)),
                     np.zeros(shape=(0, 5)),
                     np.zeros(shape=(0, 5))
                     ]

chem_test_red_1_1 = [(82.4441, 90.0294), (20, 150), (82.4441, 150),
                     (20, 90.0294)]
tors_test_red_1_1 = [np.array([[-110.,   90.],
                               [ -90.,  100.]]),
                     np.zeros(shape=(0, 5)),
                     np.zeros(shape=(0, 5)),
                     np.zeros(shape=(0, 5))
                    ]

chem_test_red_1_6 = [(68.7441, 112.9245), (20, 150), (68.7441, 150), 
                     (20, 112.9245)]
tors_test_red_1_6 = [np.array([[  20.,  -50.,  180.],
                               [  30.,  -50.,  180.],
                               [-170.,  150.,  180.],
                               [-170.,  -90.,  180.],
                               [-170.,  -80.,  180.],
                               [-140.,  100.,  180.],
                               [-140., -150.,  180.],
                               [-120.,  -70.,  180.],
                               [ -90.,  180.,  180.],
                               [ -90., -180.,  180.]]),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5)),
                      np.zeros(shape=(0, 5))
                    ]

all_chem_test_full = [chem_test_full_1_3, chem_test_full_1_1, chem_test_full_1_6]
all_tors_test_full = [tors_test_full_1_3, tors_test_full_1_1, tors_test_full_1_6]

all_chem_test_red = [chem_test_red_1_3, chem_test_red_1_1, chem_test_red_1_6]
all_tors_test_red = [tors_test_red_1_3, tors_test_red_1_1, tors_test_red_1_6]

# Tests

def test_round_down_up():
    ref = [(-10, 0), (10, 20), (170, 180), (-180, -170), (60, 70), (70, 80),
           (170.0, 180.0), (-180.0, -170.0)]
    for idx, t in enumerate(torsionals):
        assert _round_down_up(t, 10) == ref[idx]

def test_compute_cs_all():# load all and check if this run
    disaccharide_test = CheSweet(path='chesweet/lut')
    disaccharide_test.compute_cs(disaccharides_list[0], 50, 60)

def test_compute_cs_few():# load determinated disaccharides and check if this run
    disaccharides_test = CheSweet(path='chesweet/lut', disaccharides=disaccharides_list)
    disaccharides_test.compute_cs(disaccharides_list[0], 50, 60)
    disaccharides_test.compute_cs(disaccharides_list[1], 50, 60)
    disaccharides_test.compute_cs(disaccharides_list[2], 50, 60, 60)

def test_compute_cs_one():# load one disaccharide and check if this run
    che = CheSweet(path='chesweet/lut', disaccharides=[disaccharides_list[1]])
    che.compute_cs(disaccharides_list[1], 50, 60)


def test_compute_cs_full():
    for idx, disaccharide in enumerate(disaccharides_list):
        tors_list = all_tors_full[idx]
        chem_list = all_chem_full[idx]
        for i in range(0, len(tors_list)):
            if disaccharide != 'a-D-Glcp-1-1-a-D-Glcp':# 1-1 a Chi less
                chemC1_full, chemCx_full = disaccharides_full.compute_cs(disaccharide,
                                                                         tors_list[i][0], tors_list[i][1],
                                                                         tors_list[i][2], tors_list[i][3],
                                                                         tors_list[i][4])
            else:
                chemC1_full, chemCx_full = disaccharides_full.compute_cs(disaccharide,
                                                                         tors_list[i][0], tors_list[i][1],
                                                                         tors_list[i][2], tors_list[i][3])
            if chemC1_full != np.inf:
                np.testing.assert_almost_equal(ef_corr-chemC1_full, chem_list[i][0],
                                               decimal=4)
                np.testing.assert_almost_equal(ef_corr-chemCx_full, chem_list[i][1],
                                               decimal=4)
            else:# when the interpolation result is out of the caculated zone
                # not subtract ef_corr because convert inf into -inf
                np.testing.assert_almost_equal(chemC1_full, chem_list[i][0],
                                               decimal=4)
                np.testing.assert_almost_equal(chemCx_full, chem_list[i][1],
                                               decimal=4)


def test_compute_cs_red():
    for idx, disaccharide in enumerate(disaccharides_list):
        tors_list = all_tors_red[idx]
        chem_list = all_chem_red[idx]
        for i in range(0, len(tors_list)):
            if disaccharide != 'b-D-Galp-1-6-b-D-Galp':# 1-6 need omega
                chemC1_red, chemCx_red = disaccharides_red.compute_cs(disaccharide,
                                                                      tors_list[i][0], tors_list[i][1])
            else:
                chemC1_red, chemCx_red = disaccharides_red.compute_cs(disaccharide,
                                                                      tors_list[i][0], tors_list[i][1],
                                                                      tors_list[i][2])
            if chemC1_red != np.inf:
                np.testing.assert_almost_equal(ef_corr-chemC1_red, chem_list[i][0],
                                               decimal=4)
                np.testing.assert_almost_equal(ef_corr-chemCx_red, chem_list[i][1],
                                               decimal=4)
            else:# when the interpolation result is out of the caculated zone
                # not subtract ef_corr because convert inf into -inf
                np.testing.assert_almost_equal(chemC1_red, chem_list[i][0],
                                               decimal=4)
                np.testing.assert_almost_equal(chemCx_red, chem_list[i][1],
                                               decimal=4)


def test_compute_tors_full():
    for idx, disaccharide in enumerate(disaccharides_list):
        chem_list = all_chem_test_full[idx]
        tors_list = all_tors_test_full[idx]
        for i in range(0, len(chem_list)):
            cs1 = ef_corr-chem_list[i][0]
            cs2 = ef_corr-chem_list[i][1]
            x = tors_list[i]
            y = disaccharides_full.compute_tors('{}'.format(disaccharide),
                                                cs1, cs2, ef_corr=183.4, eps=0.5)
            if y.size != 0:# values into the distribution
                np.testing.assert_array_almost_equal(x, y, decimal=1, verbose=True)
            else:# values out of the distribution
                assert y.size == 0


def test_compute_tors_red():
    for idx, disaccharide in enumerate(disaccharides_list):
        chem_list = all_chem_test_red[idx]
        tors_list = all_tors_test_red[idx]
        for i in range(0, len(chem_list)):
            cs1 = ef_corr-chem_list[i][0]
            cs2 = ef_corr-chem_list[i][1]
            x = tors_list[i]
            y = disaccharides_red.compute_tors('{}'.format(disaccharide),
                                               cs1, cs2, ef_corr=183.4, eps=0.5)
            if y.size != 0:# values into the distribution
                np.testing.assert_array_almost_equal(x, y, decimal=1, verbose=True)
            else:# values out of the distribution
                assert y.size == 0



