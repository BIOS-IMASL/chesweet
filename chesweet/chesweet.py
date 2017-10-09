import glob
from os.path import basename as bn
import math
import numpy as np
from scipy.interpolate import griddata


class CheSweet():
    """
    Class to compute chemical shift or torsional angles of glycosidics bonds.
    """

    def __init__(self, path='lut', disaccharides=None, full=False):
        """
        Parameters
        ----------
        disaccharides : list of strings
            names of the disaccharides used as keys in the dictionary
        path : string
            folder of CheSweet's lookup table, by default is 'lut'
        full : Boolean
            whether to include chi's torsional angles (True) in the computation
            of chemical shifts or not (False)
        """
        self.path = path
        self.full = full
        if disaccharides is None:
            files = '{}/*[!_red]'.format(self.path)
            self.disaccharides = [bn(a) for a in glob.glob(files)]
        else:
            self.disaccharides = disaccharides

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
            chi2 torsional angle in degrees (optional)
        full: Boolean
            whether to include chi's torsional angles (True) in the computation of
            chemical shifts or not (False)
        ef_corr : float
            correction values used to turn shielding into chemical shifts. Default
            value is 183.4

        Returns
        ----------
        cs : array
            interpolated chemical shifts of the first and second carbon in
            the glycosidic bond
        """

        lt = self.lt
        # phi and psi angles in lt are compute using a 10 degree grid.
        phi_range = _round_down_up(phi, 10)
        psi_range = _round_down_up(psi, 10)
        # chi and omega angles in lt are compute in 3 position (60, -60, 180)
        if self.full:  # Check how this works with bonds 1-1 (a Chi less)
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
            the first columns is the phi angle and the second is psi.
            If the chemical shifts are outside the zone of computed values
            the function returns an empty array
        """
        # transform cs into shieldings
        cs0 = ef_corr - cs0
        cs1 = ef_corr - cs1
        x = self.lt[disaccharide]
        cond0 = (x[:,2] < cs0 + eps) & (x[:,2]  > cs0 - eps)
        cond1 = (x[:,3] < cs1 + eps) & (x[:,3]  > cs1 - eps)
        theoric_tors = x[:,:2][cond0 & cond1]
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
            try:
                lut[disaccharide] = np.fromfile('{}/{}'.format(self.path,
                                                disaccharide),
                                                sep=' ').reshape(-1, 8)
            except FileNotFoundError:
                print("The dissacharide {} is not calculated, "
                      "or the path {} is wrong".format(disaccharide, self.path))
    else:
        for disaccharide in self.disaccharides:
            try:
                # Disaccharides with 1-6 glycosidic bond
                if '-1-6-' in disaccharide:  
                    lut[disaccharide] = np.fromfile('{}/{}_red'.format(self.path,
                                                                       disaccharide),
                                                                       sep=' ').reshape(-1, 5)
                # Disaccharides with glycosidic bond different from 1-6
                else:
                    lut[disaccharide] = np.fromfile('{}/{}_red'.format(self.path,
                                                                       disaccharide),
                                                                       sep=' ').reshape(-1, 4)
            except FileNotFoundError:
                print("The dissacharide {} is not calculated, "
                      "or the path {} is wrong".format(disaccharide, self.path))

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
