import glob
from os.path import basename as bn
import pkg_resources as pkg
from os.path import join
import sys
import math
import numpy as np
from scipy.interpolate import griddata


class CheSweet():
    """
    Class to compute chemical shift or torsional angles of glycosidics bonds.
    """

    def __init__(self, path=None, full=False):
        """
        Parameters
        ----------
        path : string
            folder of lookup table, if None (default) CheSweet's lookup table will be used
        full : Boolean
            whether to include chi's torsional angles (True) in the computation
            of chemical shifts or not (False)
        """
        self.full = full
        
        if path is None:
            files = pkg.resource_filename(__name__, '/'.join(['lut']))
        else:
            files = '/'.join([path, 'lut'])
        
        if full:
            self.disaccharides = glob.glob('{}/*[!_red]'.format(files))
        else:
            self.disaccharides = glob.glob('{}/*_red'.format(files))

        self.lt = _load(self)


    def compute_cs(self, disaccharide, phi, psi, chi1=None, chi2=None,
                   chi3=None, ef_corr=183.4):
        """
        Compute the value of a chemical shift given a set of torsional angles,
        the name of a disaccharide and a look-up table.

        Parameters
        ----------
        disaccharide : string
            disaccharides names used as keys in lt dictionary
        lt : dictionary
            CheSweet's look-up table
        phi : float
            phi torsional angle in degrees
        psi : float
            psi torsional angle in degrees
        chi1 : float
            chi1 torsional angle in degrees (optional)
        chi2 : float
            chi2 torsional angle in degrees (optional)
        chi3 : float
            chi3 torsional angle in degrees (optional)
        full: Boolean
            whether to include chi's torsional angles (True) in the
            computation of chemical shifts or not (False)
        ef_corr : float
            correction values used to turn shielding into chemical shifts.
            Default value is 183.4

        Returns
        ----------
        cs : array
            interpolated chemical shifts of the first and second carbon in
            the glycosidic bond
        """

        lt = self.lt
        if not self.full:
            disaccharide = disaccharide + '_red'
        # phi and psi angles in lt are compute using a 10 degree grid.
        phi_range = _round_down_up(phi, 10)
        psi_range = _round_down_up(psi, 10)
        # chi and omega angles in lt are compute in 3 position (60, -60, 180)
        if self.full:
            chi1_n = _nearest_chi(chi1)
            chi2_n = _nearest_chi(chi2)
            if chi3 == None:# glycosidic bonds 1-1 (a Chi less)
                chi3_n = 0
            else:
                chi3_n = _nearest_chi(chi3)
            s_db = lt[disaccharide][(lt[disaccharide][:, 0] >= phi_range[0]) &
                                    (lt[disaccharide][:, 0] <= phi_range[1])]

            ss_db = s_db[(s_db[:, 1] >= psi_range[0]) &
                         (s_db[:, 1] <= psi_range[1])]
            
            data = ss_db[(ss_db[:, 2] == chi1_n) & (
                ss_db[:, 3] == chi2_n) & (ss_db[:, 4] == chi3_n)]
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
            cs = ef_corr - griddata(data[:, :2],
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

    def compute_tors(self, disaccharide, cs0, cs1, ef_corr=183.4, eps=0.5):
        """
        Compute the torsional angles given the chemical shift of the
        carbons in the glycosidic bond
        
        Parameters
        ----------
        dissacharide : string
            name of the dissacharide involved
        cs0, cs1 : float
            chemical shift of the first and second carbon on the glycosidic
            bond, respectively
        ef_corr : float
            correction values used to turn shielding into chemical shifts.
            Default value is 183.4
        eps: float
            chemical shift tolerance, increasing this value will, in general,
            increase the number of retrieved torsionals

        Returns
        ----------
        theoric_tors: array
            torsionals from the lookup table in the range of `eps`,
            if `full=True` returns torsionals phi, psi and omega/chi's
            otherwise return only phi and psi/omega.
            If the chemical shifts are outside the zone of computed values
            the function returns an empty array
        """
        # transform cs into shieldings
        cs0 = ef_corr - cs0
        cs1 = ef_corr - cs1

        if not self.full:
            disaccharide = disaccharide + '_red'

        x = self.lt[disaccharide]
        cond0 = (x[:,-2] < cs0 + eps) & (x[:,-2]  > cs0 - eps)
        cond1 = (x[:,-1] < cs1 + eps) & (x[:,-1]  > cs1 - eps)
        
        if self.full:
            if int(disaccharide.split('-')[4]) == 1:# bonds 1-1
                theoric_tors = x[:,:4][cond0 & cond1]
            else:# bonds different from 1-6
                theoric_tors = x[:,:5][cond0 & cond1]
        else:# reduced
            if x.shape[1] == 4:# bonds different from 1-6
                theoric_tors = x[:,:2][cond0 & cond1]
            else:# 1-6 bonds
                theoric_tors = x[:,:3][cond0 & cond1]
        return theoric_tors


def _load(self):
    """
    Load CheSweet's look-up table as a dictionary of arrays.

    Returns
    ----------
    lut : dictionary
        keys are the names of the dissacharides and values are arrays
        with the last two columns being the chemical shift pre-calculated 
        and the rest being torsinal angles. The number of columns in
        the arrays depends on wheter `full` is True or False
    """
    lut = {}
    if self.full:
        for disaccharide in self.disaccharides:
            lut[bn(disaccharide)] = np.fromfile(disaccharide, sep=' ').reshape(-1, 8)
    else:
        for disaccharide in self.disaccharides:
            # Disaccharides with 1-6 glycosidic bond
            if '-1-6-' in disaccharide:
                lut[bn(disaccharide)] = np.fromfile(disaccharide, sep=' ').reshape(-1, 5)
            # Disaccharides with glycosidic bond different from 1-6
            else:
                lut[bn(disaccharide)] = np.fromfile(disaccharide, sep=' ').reshape(-1, 4)

    if lut:
        return lut
    else:
        raise ValueError('I could not build a look-up table')


def _round_down_up(tor, step=10.):
    """
    Round tor down and up to multiples of step.

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
        the nearest value in the look-up table for chi
    """
    chi_rotamers = np.array([-180., -60.,  60., 180.])
    indices = np.argmin(np.abs(chi_rotamers - chi))
    return chi_rotamers[indices]
