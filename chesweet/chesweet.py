import math
import numpy as np
from scipy.interpolate import griddata


def _round_down_up(tor, step=10.):
    """
    Round tor down and up to multiples of step

    Parameters
    ----------
    tor : float
        angle to round in steps. tor can be outside the interval [-180, 180]
    step : float
        interval between the returned rounded values

    Returns
    ----------
    down_up : tuple
        rounded value of tor (down, up)
    """
    if tor > 180. or tor < -180.:
        tor_rad = tor * np.pi / 180.
        tor = math.atan2(math.sin(tor_rad), math.cos(tor_rad)) * 180. / np.pi
    down = math.floor(tor / step) * step
    up = down + step
    if up > 180:
        up = down
        down = 180 - step
    return down, up


def _nearest_chi(chi):
    """
    Compute nearest value in the look-up table for chi.

    Parameters
    ----------
    chi : float
        chi torsional angle

    Returns
    ----------
    nearest_chi : float
        The nearest value in the look-up table for chi.
    """
    chi_rotamers = np.array([-180., -60.,  60., 180.])
    indices = np.argmin(np.abs(chi_rotamers - chi))
    return chi_rotamers[indices]


def load(path, disaccharides, full=True):
    """
    Load CheSweet's look-up table as a dictionary of arrays.
    The keys are the name of the glycans
    """
    lut = {}
    if full:
        for disaccharide in disaccharides:
            try:
                lut[disaccharide] = np.fromfile(
                    '{}/{}'.format(path, disaccharide), sep=' ').reshape(-1, 8)
            except FileNotFoundError:
                print("The dissacharide {} is not calculated".format(disaccharide))
    else:
        for disaccharide in disaccharides:
            try:
                # Disaccharides with 1-6 glycosidics bond
                if '-1-6-' in disaccharide:  
                    lut[disaccharide] = np.fromfile(
                        '{}/{}_red'.format(path, disaccharide),
                        sep=' ').reshape(-1, 5)
                # Disaccharides with glycosidics bond diferent from 1-6
                else:
                    lut[disaccharide] = np.fromfile(
                        '{}/{}_red'.format(path, disaccharide),
                        sep=' ').reshape(-1, 4)
                    print('{}/{}_red'.format(path, disaccharide))
            except FileNotFoundError:
                print("The dissacharide {} is not calculated".format(disaccharide))
    return lut


def compute_cs(disaccharide, lt, phi, psi, chi1=None, chi2=None, chi3=None,
               full=True, ef_corr=183.4):
    """
    Compute the value of a chemical shift given a set of torsional angles,
    the name of a disaccharide and a look-up table.
    disaccharide : string
        blablablabla
    lt : Dictionary
        CheSweet's look-up table
    phi : float
        Phi torsional angle in degrees
    psi : float
        Psi torsional angle in degrees
    chi1 : float
        Chi1 torsional angle in degrees (optional)
    chi2 : float
        Chi2 torsional angle in degrees (optional)
    chi3 : float
        Chi2 torsional angle in degrees (optional)
    full: Boolean
        Whether to include chi's torsional angles (True) in the computation of
        chemical shifts or not (False)
    ef_corr : float
        Correction values used to turn shielding into chemical shifts. Default
        value is 183.4
    """

    # phi and psi angles in lt are compute using a 10 degree grid.
    phi_range = _round_down_up(phi, 10)
    psi_range = _round_down_up(psi, 10)
    # chi and omega angles in lt are compute in 3 position (60, -60, 180)
    if full:  # Check how this works with bonds 1-1 (a Chi less)
        chi1_n = _nearest_chi(chi1)
        chi2_n = _nearest_chi(chi2)
        chi3_n = _nearest_chi(chi3)
        s_db = lt[disaccharide][(lt[disaccharide][:, 0] >= phi_range[0]) &
                                (lt[disaccharide][:, 0] <= phi_range[1])]

        ss_db = s_db[(s_db[:, 1] >= psi_range[0]) &
                     (s_db[:, 1] <= psi_range[1])]
        data = ss_db[(ss_db[:, 2] == chi1_n) & (
            ss_db[:, 3] == chi2_n) & (ss_db[:, 4] == chi1_n)]
    else:
        s_lt = lt[disaccharide][(lt[disaccharide][:, 0] >= phi_range[0]) &
                                (lt[disaccharide][:, 0] <= phi_range[1])]
        ss_lt = s_lt[(s_lt[:, 1] >= psi_range[0]) &
                     (s_lt[:, 1] <= psi_range[1])]

        # phi, psi, chemical shift C1, chemical shift Cx
        if lt[disaccharide].shape[1] == 4:
            data = ss_lt
        # phi, psi, omega, chemical shift C1, chemical shift Cx
        elif lt[disaccharide].shape[1] == 5:
            omega_n = _nearest_chi(chi1)
            data = ss_lt[(ss_lt[:, 2] == omega_n)]

    d_shape = data.shape[0]
    # we are outside the zone of computed values
    if d_shape == 0:
        cs = np.array([np.inf, np.inf])
    # we hit a border of the computed values!
    elif d_shape < 4:
        cs = ef_corr - \
            griddata(data[:, :2],
                     data[:, -2:],
                     (phi, psi),
                     method='nearest')
    # life is sweet!
    elif d_shape == 4:
        cs = ef_corr - griddata(data[:, :2],
                                data[:, -2:],
                                (phi, psi),
                                method='linear')
    return cs
