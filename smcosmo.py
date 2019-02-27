"""
smcosmo: simulating SkyMapper (sm) supernova cosmology using SNCosmo.
u6253153@anu.edu.au,  2016 - 2019

Aims to:
- Generate a random, realistic set of type 1a supernovae.
- Simulate SkyMapper and Kepler Space Telescope observations of these supernovae, using SNCosmo.
- Fit SALT2 light-curve models to various combinations of these observations, to recover supernovae parameters (x1, x0, t0, c).
- Analyse how each 'combination' of observations performed (separate script)

Packaged with:
- converting_data.py (converts KST photometry to a readable format for this program)
- distance_wrapper.py (wrapper for luminosity_distance_fitter)
- filters.py (registers SkyMapper and KST transmission filters, requires filter data)
- luminosity_distance_fitter (calculates cosmological distances from fitted light-curve parameters)
- real_data_wrapper (wrapper for fitting real KST data, must be converted first - mainly an example script)
- statistical_analysis.py (produces residual plots etc. to evaluate quality of fits)
- wrapper.py (wrapper for this main script, smcosmo.py)

Requires installation of (not yet exhaustive):
- SNCosmo:  https://sncosmo.readthedocs.io/en/v1.4.x/index.html
- Emcee:
- sfdmap:
"""

import os
import pickle
import math
import numpy as np
import random
import datetime
import matplotlib.pyplot as plt
import sfdmap
import sncosmo
import copy
# filters.py should be in working directory
import filters
from shutil import copyfile
from sncosmo import mcmc_lc as mcmc_fit
from sncosmo import fit_lc as chi_fit
from astropy.time import Time
from numpy import random
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from astropy.cosmology import FlatLambdaCDM
from astropy.extern.six.moves import range
from astropy.table import Table, vstack
from scipy.stats import exponweib
from matplotlib.backends.backend_pdf import PdfPages

""" Important edit notes:
    - SM observable candidate range: RA(0-360), Dec(-90,10)
"""


# --- The parameters below are safe to change - make sure to check they are correct before running a big simulation! ---


# SkyMapper observing cadence (days)
# NOTE: for SM cadence to vary randomly between 1 and 2 days, cad_sm = 'well-sampled'
#       for SM cadence to vary randomly betwen 3 and 5 days, cad_sm = 'poorly-sampled'
cad_sm = 'well-sampled'

# Kepler observing cadence (days)
cad_k = 1./4. # 6 hours

# Number of observations in SkyMapper V filter (good seeing, centered around lightcurve peak)
# NOTE: this is currently disabled to mimic Pan-STARRS filter sets (g, r, i) - set v_obs > 0 to reactivate.
v_obs = 0

# Skymapper min and max observable redshift, representing redshift bounds for simulations.
zmin = 0.001
zmax = 0.15

# Parameter distribution to use for c and x1 (from Scolnic and Kessler 2016, Table 1)
# NOTE: options are 'lowz_g10', 'lowz_c11', 'ps1_g10', 'ps1_c11' (reflecting sample and intrinsic scatter model)
#       c bounds [-0.3, 0.3]; x bounds [-3, 3]
#       PS1 distributions need to be updated to values from Pantheon paper (Scolnic 2017)
lcdist = 'lowz_g10'

# Change path to location of dustmaps (must be downloaded first)
dustmap = sfdmap.SFDMap("/Users/macloan4/GTaylor/lightcurves/sfddata-master")
#dustmap = sfdmap.SFDMap("/home/georgie/sfddata-master")#thinkpad
#dustmap = sfdmap.SFDMap("/home/gtaylor/sfddata-master")#motley


# ---------------- You can change the below parameters if you must, but it'll probably break something -----------------


# Fitzpatrick 99 reddening law with Rv=3.1, as prescribed in Schlafly and Finkbeiner 2011.
dust = sncosmo.F99Dust()

# Salt2 model template (assuming most recent version i.e. JLA training)
model = sncosmo.Model(source='salt2',
                      effects=[dust, dust],
                      effect_names=['host', 'mw'],
                      effect_frames=['rest', 'obs'])


# -------------------------------------------------- Utility functions -------------------------------------------------

# Sets random number seed to system time
random.seed()

def save_obj(obj, name):
    """ Utility function for saving dictionary 'obj' as 'name'. """

    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    """ Utility function for opening dictionary. """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def ensure_dir(f):
    """ Checks if specified path exists, and creates it if not. """

    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def get_coords(nSNe, ra_range=0, dec_range=0):
    """ Generates random j2000 degree coordinates for SN
        Default values use SkyMapper observable range
    """

    if ra_range == 0:
        # Right ascension
        ra = np.random.uniform(0, 360, nSNe).tolist()
        # Declination
        dec = np.random.uniform(-90, 10, nSNe).tolist()
    else:
        # Right ascension (based on input values, probably from a K2 field)
        ra = np.random.uniform(ra_range[0], ra_range[1], nSNe).tolist()
        # Declination
        dec = np.random.uniform(dec_range[0], dec_range[1], nSNe).tolist()
    coords = [ra, dec]
    return coords

def get_skynoise(n):
    """ Calculates static skynoise in counts (background contribution to flux measurement error), assuming
        seeing-dominated gaussian psf with perfect psf photometry.
        seeing = fwhm of psf (pixels)
        sigma_psf = std. deviation of psf (pixels)
        sigma_pixel = background noise in single pixel (in counts) = bkgsig in data from Anais Moller
    """
    skynoise = []
    # Distributions recovered from SkyMapper data, don't change them without reason.
    for i in range(n):
        # Seeing = FWHM of psf, in arc sec.
        seeing = exponweib.rvs(461.87921971531955, 0.6860198725258686, loc=-0.869494090011727,
                               scale=0.43722954355185195, size=1)
        # Find probability of distribution giving generated seeing value (or less)
        percent = exponweib.cdf(seeing, 461.87921971531955, 0.6860198725258686, loc=-0.869494090011727,
                                scale=0.43722954355185195)
        # Randomly vary probability by 15%, staying within [0.01,1]
        # The 0.01 condition is to prevent negative values (see shape of distribution).
        if percent <= 0.15:
            percent_new = percent + np.random.uniform(-percent, percent)
        elif percent >= 0.85:
            percent_new = percent + np.random.uniform(-(1 - percent), (1 - percent))
        else:
            percent_new = percent + np.random.uniform(-0.15, 0.15)
        if percent_new < 0.01:
            percent_new = 0.01
        # sigma_pixel = background noise in a single pixel, in counts
        sigma_pixel = exponweib.ppf(percent_new, 26.47111428278329, 1.014405795857852, loc=-16.985179115012222,
                                    scale=9.81027950453267)
        # converting seeing from arc sec to pixels, with a scale of 0.5" per pixel - then converting to standard deviation.
        sigma_psf = 2 * seeing / 2.355
        x = 4 * math.pi * sigma_psf * sigma_pixel
        skynoise.append(x)
    # Flatten the nested list
    flat_sn = [item for sublist in skynoise for item in sublist]
    return flat_sn

def write_params(folder, sn):
    """ Writes observational parameters to a file for reference, not that important unless you forget what you did."""

    ensure_dir(folder)

    p_file = folder + 'observing_parameters.txt'
    pf = open(p_file, 'w')
    pf.write('Kepler cadence: %s days. \n\
    Number of SN: %s. \n\
    Fitting method: chi-squared initial guess passed to MCMC. \n\
    Distribution: %s.'% (cad_k, sn, lcdist))

    pf.close()

def modified_zdist(zmin, zmax, nsim, ratefunc=lambda z: 1.e-4,
                   cosmo=FlatLambdaCDM(H0=70.0, Om0=0.3)):
    """ *** Modified from original SNCosmo method to allow a user-defined number of simulations ***
         *** Majority of method is part of SNCosmo source code: https://sncosmo.readthedocs.io ***

    Generate a distribution of redshifts.

    Generates the correct redshift distribution and number of SNe, given
    the input volumetric SN rate, the cosmology, and the observed area and
    time.

    Parameters
    ----------
    zmin, zmax : float
        Minimum and maximum redshift.
    nsim : int
        Number of redshifts to generate
    ratefunc : callable
        A callable that accepts a single float (redshift) and returns the
        comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.
        The default is a function that returns ``1.e-4``.
    cosmo : `~astropy.cosmology.Cosmology`, optional
        Cosmology used to determine volume. The default is a FlatLambdaCDM
        cosmology with ``Om0=0.3``, ``H0=70.0``.
    """

    # Get comoving volume in each redshift shell.
    z_bins = 100  # Good enough for now.
    z_binedges = np.linspace(zmin, zmax, z_bins + 1)
    z_binctrs = 0.5 * (z_binedges[1:] + z_binedges[:-1])
    sphere_vols = cosmo.comoving_volume(z_binedges).value
    shell_vols = sphere_vols[1:] - sphere_vols[:-1]

    # SN / (observer year) in shell
    shell_snrate = np.array([shell_vols[i] *
                             ratefunc(z_binctrs[i]) / (1. + z_binctrs[i])
                             for i in range(z_bins)])

    # SN / (observer year) within z_binedges
    vol_snrate = np.zeros_like(z_binedges)
    vol_snrate[1:] = np.add.accumulate(shell_snrate)

    # Create a ppf (inverse cdf). We'll use this later to get
    # a random SN redshift from the distribution.
    snrate_cdf = vol_snrate / vol_snrate[-1]
    snrate_ppf = Spline1d(snrate_cdf, z_binedges, k=1)

    # The bit we changed #
    for i in range(nsim):
        yield float(snrate_ppf(random.random()))

def uniform_sample(xmin, xmax, ymin, ymax):
    """ Returns random value in 2d box, that helps with probabilities later """
    return np.random.random(2,)*np.array([xmax-xmin,ymax-ymin])+np.array([xmin, ymin])

def one_skew_rv(max_prob, sigma_m, sigma_p):
    """ Generate a distribution for c or x1, as defined in Scolnic and Kessler 2016.

    Returns a random float in the defined skew-normal distribution

    Parameters
    ----------

    max_prob: float
        value with maximum probability
    sigma_m: float
        standard deviation below max_prob
    sigma_p: float
        standard deviation above max_prob

    """
    while True:
        sample = uniform_sample(-3, 3, 0, 1)
        # x = random.uniform(-1,1)
        x = sample[0]
        numer = -((x-max_prob)**2)
        if x <= max_prob:
            prob =  math.exp(numer/(2*sigma_m**2))
        else:
            prob =  math.exp(numer/(2*sigma_p**2))
        if prob > sample[1]:
            break
    rv = sample[0]
    return rv

def one_skew_prob_c(param):
    """ Generate a probability for obtaining a certain param c, as defined in Scolnic and Kessler 2016.

    Returns the probability of getting 'param' in the defined distribution

    Parameters
    ----------

    max_prob: float
        value with maximum probability
    sigma_m: float
        standard deviation below max_prob
    sigma_p: float
        standard deviation above max_prob
    param: float
        value of param, for which you are trying to find the probability.
    normfactor: calculated as the area under the graph, using the trapz function.
        c:  lowz_g10: 0.21660889801728506;
            ps1_g10: 0.18780753712230902;
            lowz_c11: 0.1891463376131865;
            ps1_c11: 0.1654646632528336;
        x1: lowz_g10: 3.7561727085403804;
            ps1_g10: 1.74227465242822;
            lowz_c11: 3.739187213082589;
            ps1_c11: 1.7610474962637968;
    """

    # Handling the different normalisations for each lc distribution type (be careful, referencing global params)
    if lcdist == 'lowz_g10':
        normfactor = 0.21660889801728506
        max_prob = -0.055
        sigma_m = 0.023
        sigma_p = 0.150
    elif lcdist == 'ps1_g10':
        normfactor = 0.18780753712230902
        max_prob = -0.077
        sigma_m = 0.029
        sigma_p = 0.121
    elif lcdist == 'lowz_c11':
        normfactor = 0.1891463376131865
        max_prob = -0.069
        sigma_m = 0.003
        sigma_p = 0.148
    elif lcdist == 'ps1_c11':
        normfactor = 0.1654646632528336
        max_prob = -0.103
        sigma_m = 0.003
        sigma_p = 0.129
    else:
        raise ValueError('Invalid lcdist value - please use \'lowz_g10\', '
                         '\'ps1_g10\', \'lowz_c11\', or \'ps1_c11\'')
    x = param
    numer = -((x - max_prob) ** 2)
    # Low-z x1 condition (prevents math error for log of 0 probability case)
    if x <= -0.3 or x >= 0.3:
        prob = 0.000000000001
    elif x <= max_prob:
        prob =  (math.exp(numer/(2*sigma_m**2)))/normfactor
    else:
        prob =  (math.exp(numer/(2*sigma_p**2)))/normfactor
    return prob

def one_skew_prob_x1(param):
    """ Generate a probability for obtaining a certain param x1, as defined in Scolnic and Kessler 2016.

    Returns the probability of getting 'param' in the defined distribution

    Parameters
    ----------

    max_prob: float
        value with maximum probability
    sigma_m: float
        standard deviation below max_prob
    sigma_p: float
        standard deviation above max_prob
    param: float
        value of param, for which you are trying to find the probability.
    normfactor: calculated as the area under the graph, using the trapz function.
        c:  lowz_g10: 0.21660889801728506;
            ps1_g10: 0.18780753712230902;
            lowz_c11: 0.1891463376131865;
            ps1_c11: 0.1654646632528336;
        x1: lowz_g10: 3.7561727085403804;
            ps1_g10: 1.74227465242822;
            lowz_c11: 3.739187213082589;
            ps1_c11: 1.7610474962637968;
    """

    # Handling the different normalisations for each lc distribution type (be careful, referencing global params)
    if lcdist == 'lowz_g10':
        normfactor = 3.7561727085403804
        max_prob = 0.436
        sigma_m = 3.118
        sigma_p = 0.724
    elif lcdist == 'ps1_g10':
        normfactor = 1.74227465242822
        max_prob = 0.604
        sigma_m = 1.029
        sigma_p = 0.363
    elif lcdist == 'lowz_c11':
        normfactor = 3.739187213082589
        max_prob = 0.419
        sigma_m = 3.024
        sigma_p = 0.742
    elif lcdist == 'ps1_c11':
        normfactor = 1.7610474962637968
        max_prob = 0.589
        sigma_m = 1.026
        sigma_p = 0.381
    else:
        raise ValueError('Invalid lcdist value - please use \'lowz_g10\', '
                         '\'ps1_g10\', \'lowz_c11\', or \'ps1_c11\'')
    x = param
    numer = -((x - max_prob) ** 2)
    # Low-z x1 condition (prevents math error for log of 0 probability case)
    if x <= -3. or x >= 3.:
        prob = 0.000000000001
    elif x <= max_prob:
        prob =  (math.exp(numer/(2*sigma_m**2)))/normfactor
    else:
        prob =  (math.exp(numer/(2*sigma_p**2)))/normfactor
    return prob

def k2_info(campaign):
    """
    !!! Not really used anymore, but in theory it works.
    Determines field size (deg^2), coordinates, and observed dates, based on a given of K2 campaign
    Note: fields 9, 16, 17, 19, and 20 are forward-facing campaigns, for which simultaneous observations from the ground
     are possible throughout the duration of the campaign.
    :param campaign: K2 campaign to include in sample, of format 'C1'
    :return: observing period (list of bounds), right ascension (list of bounds), declination (list of bounds)
    """
    if campaign=='C0':
        # 2014 Mar 08 - 2014 May 27, MJD
        obs_period = [56724., 56804.]
        ra = [90, 106]
        dec = [14, 30]
        print "Warning: C0 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C1':
        # 2014 May 30 - 2014 Aug 21, MJD
        obs_period = [56807., 56890.]
        ra = [166, 182]
        dec = [-6, 10]
    elif campaign=='C2':
        # 2014 Aug 23 - 2014 Nov 13, MJD
        obs_period = [56892., 56974.]
        ra = [238, 256]
        dec = [-30, -14]
    elif campaign=='C3':
        # 2014 Nov 14 - 2015 Feb 03, MJD
        obs_period = [56975., 57056.]
        ra = [328, 346]
        dec = [-20, -2]
    elif campaign=='C4':
        # 2015 Feb 07 - 2015 Apr 23, MJD
        obs_period = [57060., 57135.]
        ra = [52, 68]
        dec = [12, 28]
        print "Warning: C4 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C5':
        # 2015 Apr 27 - 2015 Jul 10, MJD
        obs_period = [57139., 57213.]
        ra = [124, 138]
        dec = [10, 24]
        print "Warning: C5 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C6':
        # 2015 Jul 14 - 2015 Sep 30, MJD
        obs_period = [57217., 57295.]
        ra = [198, 212]
        dec = [-18, -4]
    elif campaign=='C7':
        # 2015 Oct 04 - 2015 Dec 26, MJD
        obs_period = [57299., 57382.]
        ra = [280, 296]
        dec = [-30, -16]
    elif campaign=='C8':
        # 2016 Jan 03 - 2016 Mar 23, MJD
        obs_period = [57390., 57470.]
        ra = [8, 24]
        dec = [-2, 14]
        print "Warning: C8 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C9':
        # 2016 Apr 21 - 2016 Jul 01, MJD
        obs_period = [57499., 57570.]
        ra = [262, 278]
        dec = [-14, -28]
    elif campaign=='C10':
        # 2016 Jul 06 - 2016 Sep 20, MJD
        obs_period = [57575., 57651.]
        ra = [180, 194]
        dec = [-10, 4]
    elif campaign=='C11':
        # 2016 Sep 24 - 2016 Dec 08, MJD
        obs_period = [57655., 57730.]
        ra = [252, 268]
        dec = [-32, -12]
    elif campaign=='C12':
        # 2016 Dec 15 - 2017 Mar 04, MJD
        obs_period = [57737., 57816.]
        ra = [344, 360]
        dec = [-12, 4]
    elif campaign=='C13':
        # 2017 Mar 08 - 2017 May 27, MJD
        obs_period = [57820., 57900.]
        ra = [64, 82]
        dec = [14, 28]
        print "Warning: C13 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C14':
        # 2017 May 31 - 2017 Aug 19, MJD
        obs_period = [57904., 57984.]
        ra = [154, 168]
        dec = [0, 14]
        print "Warning: C14 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C15':
        # 2017 Aug 23 - 2017 Nov 20, MJD
        obs_period = [57988., 58077]
        ra = [226, 242]
        dec = [-26, -12]
    elif campaign=='C16':
        # 2017 Dec 07 - 2018 Feb 25, MJD
        obs_period = [58094., 58174.]
        ra = [126, 142]
        dec = [12, 26]
        print "Warning: C16 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C17':
        # 2018 Mar 01 - 2018 May 08, MJD
        obs_period = [58178., 58246.]
        ra = [196, 212]
        dec = [-16, 0]
    elif campaign=='C18':
        # 2018 May 12 - 2018 Aug 02, MJD
        obs_period = [58250., 58332.]
        ra = [122, 138]
        dec = [10, 26]
        print "Warning: C18 field is outside field observable by SkyMapper.  Interpret results with caution."
    elif campaign=='C19':
        # 2018 Aug 06 - 2018 Oct 10, MJD
        obs_period = [58336., 58401.]
        ra = [340, 356]
        dec = [-12, 4]
    elif campaign=='C20':
        # 2018 Oct 15 - 2019 Jan 05, MJD
        obs_period = [58406., 58488.]
        ra = [60, 76]
        dec = [16, 30]
        print "Warning: C20 field is outside field observable by SkyMapper.  Interpret results with caution."
    else:
        print "Invalid K2 campaign entered.  Exiting."
        exit
    return obs_period, ra, dec


# ------------------------------------------ Generating Sets of Supernovae ---------------------------------------------

def simulate_sn_set(folder, nSNe=0, campaign=0):

    """ Generates full set of supernovae and observation parameters for all telescopes and filters, to be used by
    simulate_lc
    folder = folder to store data
    nSNe = number of SN to simulate - leave at 0 to trigger K2 simulations (requires campaign != 0)
    campaign = K2 campaign to simulate if nSNe = 0, e.g. 'C1'
    """

    ensure_dir(folder)
    params = []
    observations = []

    # ---------------------------------------------- Supernova parameters ----------------------------------------------
    # z = redshift
    if nSNe == 0:
        # Generate reasonable number of SN from sncosmo redshift distribution, based on some K2 campaign with an area of
        #  100 square degrees.
        obs_period, ra, dec = k2_info(campaign)
        tmin = obs_period[0]
        tmax = obs_period[1]
        z = list(sncosmo.zdist(zmin, zmax, time=(tmax - tmin), area=100))
        nSNe = len(z)
        # Equatorial coordinates
        coords = get_coords(nSNe, ra, dec)
        print coords
        print 'Obtained %s' % nSNe + ' redshifts from sncosmo distribution. \n'
    else:
        # Set number of SN and generate redshift (for statistical sample) from SNCosmo distribution using modified
        # method
        # Possible observing period: 365 days from today's date in mjd.
        now = Time([str(datetime.date.today())], format='iso', scale='utc')
        now_mjd = now.mjd
        tmin = now_mjd[0]
        tmax = tmin + 365
        z = list(modified_zdist(zmin, zmax, nSNe))
        # Equatorial coordinates
        coords = get_coords(nSNe)
        print 'Obtained %s' % nSNe + ' redshifts from SNCosmo distribution. \n'

    # t0 = observer-frame time corresponding to source's phase=0 (mjd).
    t0 = np.random.uniform(tmin, tmax, nSNe)

    # c = colour
    # distribution from Scolnic and Kessler 2016, Table 1
    c = []
    for i in range(nSNe):
        # Draws colour from distribution based on tag specified at start of script.
        if lcdist == 'lowz_g10':
            c.append(one_skew_rv(-0.055, 0.023, 0.150))
        elif lcdist == 'ps1_g10':
            c.append(one_skew_rv(-0.077, 0.029, 0.121))
        elif lcdist == 'lowz_c11':
            c.append(one_skew_rv(-0.069, 0.003, 0.148))
        elif lcdist == 'ps1_c11':
            c.append(one_skew_rv(-0.103, 0.003, 0.129))
        else:
            raise ValueError('Invalid lcdist value - please use \'lowz_g10\', '
                             '\'ps1_g10\', \'lowz_c11\', or \'ps1_c11\'')

    # x1 = stretch
    # distribution from Scolnic and Kessler 2016
    x1 = []
    for i in range(nSNe):
        if lcdist == 'lowz_g10':
            x1.append(one_skew_rv(0.436, 3.118, 0.724))
        elif lcdist == 'ps1_g10':
            x1.append(one_skew_rv(0.604, 1.029, 0.363))
        elif lcdist == 'lowz_c11':
            x1.append(one_skew_rv(0.419, 3.024, 0.742))
        elif lcdist == 'ps1_c11':
            x1.append(one_skew_rv(0.589, 1.026, 0.381))

    # x0 = scaling factor, set as prescribed by SNCosmo with some scatter
    for i in range(nSNe):
        # Absoute B-band magnitude of SN, with some scatter introduced.
        mabs = np.random.normal(-19.3, 0.3)
        model.set(z=z[i])
        model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        # x0 obtained from SNCosmo method
        x0 = model.get('x0')
        p = {'z': z[i], 't0': t0[i], 'x0': x0, 'x1': x1[i], 'c': c[i]}
        params.append(p)

    # ---------------------------------------------- Observing parameters ----------------------------------------------

    # times of each observation
    for t in range(nSNe):

        # -------------------------- SKYMAPPER --------------------------
        # Time of initial detection (mjd) (2 - 15 days BEFORE peak)
        td_sm = np.random.randint(-15, -2)
        t_det_sm = copy.deepcopy(td_sm)

        # Total observing period (days)
        to_sm = np.random.randint(25, 65)
        t_obs_sm = copy.deepcopy(to_sm)

        # Observation times
            # Well sampled
        time_hold_ws = params[t].get('t0') + td_sm
        t_sm_ws = []
        while time_hold_ws <= (params[t].get('t0') + td_sm + to_sm):
            t_sm_ws.append(time_hold_ws)
            time_hold_ws = time_hold_ws + random.randint(1, 2)
            # Poorly sampled
        time_hold_ps = params[t].get('t0') + td_sm
        t_sm_ps = []
        while time_hold_ps <= (params[t].get('t0') + td_sm + to_sm):
            t_sm_ps.append(time_hold_ps)
            time_hold_ps = time_hold_ps + random.randint(3, 5)

        n_obs_sm = [len(t_sm_ws), len(t_sm_ps)]
        time_sm = [t_sm_ws, t_sm_ps]

        # Zero points
        # For standard filter set (g, r, i) - used in 'bad seeing'
        zp_g_bad_ws = (np.random.normal(26.82, 0.79, len(t_sm_ws)))
        zp_i_bad_ws = (np.random.normal(25.21, 0.36, len(t_sm_ws)))
        zp_r_bad_ws = (np.random.normal(26.71, 0.76, len(t_sm_ws)))
        zp_g_bad_ps = (np.random.normal(26.82, 0.79, len(t_sm_ps)))
        zp_i_bad_ps = (np.random.normal(25.21, 0.36, len(t_sm_ps)))
        zp_r_bad_ps = (np.random.normal(26.71, 0.76, len(t_sm_ps)))
        zp_gri = [[zp_g_bad_ws, zp_i_bad_ws, zp_r_bad_ws], [zp_g_bad_ps, zp_i_bad_ps, zp_r_bad_ps]]

        # For extended filter set (g, r, i, v) - used in 'good seeing'
        zp_g_good_ws = (np.random.normal(26.87, 0.68, len(t_sm_ws)))
        zp_i_good_ws = (np.random.normal(25.85, 0.81, len(t_sm_ws)))
        zp_r_good_ws = (np.random.normal(26.63, 0.67, len(t_sm_ws)))
        zp_g_good_ps = (np.random.normal(26.87, 0.68, len(t_sm_ps)))
        zp_i_good_ps = (np.random.normal(25.85, 0.81, len(t_sm_ps)))
        zp_r_good_ps = (np.random.normal(26.63, 0.67, len(t_sm_ps)))
        zp_v_good = (np.random.normal(24.91, 0.70, v_obs))
        zp_griv = [[zp_g_good_ws, zp_i_good_ws, zp_r_good_ws, zp_v_good],[zp_g_good_ps, zp_i_good_ps, zp_r_good_ps, zp_v_good]]

        # Skynoise
        # For standard filter set (g, r, i) - used in 'bad seeing'
        sn_all_ws = [get_skynoise(len(t_sm_ws))]
        sn_all_ps = [get_skynoise(len(t_sm_ps))]
        sn_gri = [sn_all_ws*3, sn_all_ps*3]

        # For extended filter set (g, r, i, v) - used in 'good seeing'

        sn_griv = [sn_all_ws*4, sn_all_ps*4]


        # KEPLER--------------------
        if campaign == 0:
            # Time of init detection (before t0)
            td_k = np.random.randint(-17, -15)
            # Total observing period (days)
            to_k = np.random.randint(80, 85)
        else:
            # Time of init detection (before t0)
            td_k = params[t].get('t0')-tmin
            # Total observing period (days)
            to_k = tmax-tmin
        t_det_k = td_k
        t_obs_k = to_k

        # Observation times
        t_k = (params[t].get('t0')
               + np.arange(td_k, td_k + to_k, cad_k)).tolist()
        n_obs_k = len(t_k)
        time_k = t_k

        # Zero points
        zp_k = [25.47] * len(t_k)

        # Skynoise
        sn_k = [0.1] * len(t_k)

        # OUTPUTS -------------------
        # Possibly inefficient way of going about this...
        # List of observing properties for each SN, to return
        obs_util = [t_det_sm, t_det_k,
                    t_obs_sm, t_obs_k,
                    time_sm, time_k,
                    n_obs_sm, n_obs_k,
                    zp_gri, zp_griv, zp_k,
                    sn_gri, sn_griv, sn_k]
        observations.append(obs_util)

    # Save as dictionary ----------------------------------------------
    dict_out = {'Parameters': params, 'Coordinates': coords,
                'Observations': observations, 'nSNe': nSNe}
    save_obj(dict_out, folder + 'sn_dict')
    print 'Properties of simulation saved in %ssn_dict.pkl' % folder
    return nSNe

def make_sn_dict(folder, redshift, ra, dec):

    """ Parallel method to simulate_sn_set in the main build, to produce a dictionary of SN properties for real data
        redshift = list of redshifts
        ra = list of right ascensions
        dec = list of declinations
        make sure the order is correct :)
    """

    ensure_dir(folder)

    params = []
    observations = []

    # Supernova parameters ----------------------------------------------

    nSNe = len(redshift)
    # z = redshift
    coords = [ra,dec]

    # This is redundant, just in for now for consistency with real code
    for i in range(nSNe):
        p = {'z': redshift[i]}
        params.append(p)

    # Save as dictionary ----------------------------------------------
    dict_out = {'Parameters': params, 'Coordinates': coords,
                'Observations': observations, 'nSNe': nSNe}
    save_obj(dict_out, folder + 'sn_dict')
    print 'Properties of SN saved in %ssn_dict.pkl' % folder
    return nSNe


# B: SIMULATING OBSERVATIONS --------------------------------------------------

def simulate_lc(parent_folder, child_folder='TestFiles/', scope='sm', sm_cad='well-sampled'):
    """ Simulate 'observed' light-curves.
        Defaults:  scope = 'sm' - which telescope to use in simulations.
                            Accepted values are 'sm' or 'kst'
                   folder = 'Testfiles/' (path to store outputs)
                    observing properties
    """

    # Maintenance.
    child_folder = parent_folder + child_folder
    filters.register_filters()
    ensure_dir(child_folder)

    bands = []
    lcs = []
    time_sm = []
    time_k = []
    n_obs_sm = []
    n_obs_k = []
    zp_gri = []
    zp_griv = []
    zp_k = []
    sn_gri = []
    sn_griv = []
    sn_k = []

    # Setting bandpasses
    if scope == 'sm':
        if v_obs > 0:
            bands = ['smg', 'smr', 'smi', 'smv']
        else:
            bands = ['smg', 'smr', 'smi']
    elif scope == 'kst':
        bands = ['kst']
    else:
        print 'Hey, that\'s not a real telescope!  Please enter ' \
              'scope = \'sm\' for SkyMapper, or \'kst\' for Kepler. '

    # Load properties of simulation set
    sn_dict = load_obj(parent_folder+'sn_dict')
    params = sn_dict['Parameters']
    obs_in = sn_dict['Observations']
    nSNe = sn_dict['nSNe']

    if sm_cad == 'well-sampled':
        for n in range(nSNe):
            time_sm.append(obs_in[n][4][0])
            time_k.append(obs_in[n][5])
            n_obs_sm.append(obs_in[n][6][0])
            n_obs_k.append(obs_in[n][7])
            zp_gri.append(obs_in[n][8][0])
            zp_griv.append(obs_in[n][9][0])
            zp_k.append(obs_in[n][10])
            sn_gri.append(obs_in[n][11][0])
            sn_griv.append(obs_in[n][12][0])
            sn_k.append(obs_in[n][13])
    elif sm_cad == 'poorly-sampled':
        for n in range(nSNe):
            time_sm.append(obs_in[n][4][1])
            time_k.append(obs_in[n][5])
            n_obs_sm.append(obs_in[n][6][1])
            n_obs_k.append(obs_in[n][7])
            zp_gri.append(obs_in[n][8][1])
            zp_griv.append(obs_in[n][9][1])
            zp_k.append(obs_in[n][10])
            sn_gri.append(obs_in[n][11][1])
            sn_griv.append(obs_in[n][12][1])
            sn_k.append(obs_in[n][13])
    else:
        raise ValueError('Invalid SkyMapper cadence.  Please use \'well-sampled\' (1-2 days) or \'poorly-sampled\' '
                         '(3-5 days).')

    # RUNNING SIMULATION --------------------------------------

    true_file = child_folder + 'true_parameters.txt'
    tf = open(true_file, 'w')

    for t in range(nSNe):

        # time of observations for all filters (mjd)
        o_t = []
        observing_bands = []
        n_points = 0

        if scope == 'sm':

            # Adds 12s offset time between observations in different sm filters
            # (for slewing)
            for x in range(3):
                j = np.array(time_sm[t]) + 0.00013888888 * x
                o_t.extend(j.tolist())
                a = n_obs_sm[t] * [bands[x]]
                observing_bands.extend(a)
                n_points += n_obs_sm[t]

            # Handle v filter observations
            if v_obs > 0:
                peak = params[t].get('t0')

                # Find observation time closest to peak
                peak_hold = min(time_sm[t], key=lambda x: abs( - peak))
                peak_index = time_sm[t].index(peak_hold)

                # Select v_obs number of obervations after peak and add slew
                v_times = time_sm[t][peak_index:peak_index+v_obs]
                j = np.array(v_times) + 0.00013888888 * 3
                o_t.extend(j.tolist())
                a = v_obs * [bands[3]]
                observing_bands.extend(a)
                n_points += v_obs

                    # CHECK THESE BAD BOIS
            # Sets zp
            if v_obs > 0:
                zp = zp_griv[t]
            else:
                zp = zp_gri[t]

            # Sets skynoise
            if v_obs > 0:
                skynoise = sn_griv[t]
            else:
                skynoise = sn_gri[t]

        else:
            o_t.extend(time_k[t])
            k_hold = n_obs_k[t] * ['kst']
            observing_bands.extend(k_hold)
            n_points = n_obs_k[t]

            # Sets zp
            zp = [zp_k[t]]

            # Sets skynoise
            skynoise = [sn_k[t]]

        # Flattening lists of lists
        observing_time = o_t
        skynoise = [item for sublist in skynoise for item in sublist]
        zp = [item for sublist in zp for item in sublist]

        observing_dictionary = {'band': observing_bands,
                                'time': observing_time,
                                'zp': zp,
                                'zpsys': n_points*['ab'],
                                'gain': n_points*[1.0],
                                'skynoise': skynoise
                                }
        obs = Table(observing_dictionary)

        lc = sncosmo.realize_lcs(obs, model, [params[t]])
        lcs.append(lc[0])

        name = child_folder + 'observed_lc_%s.txt' % (t + 1)

        # Write observations.
        sncosmo.write_lc(lc[0], name)

        print 'Simulated observations for supernova %s saved in %s'\
                   %((t + 1), name)

        # Write true parameters in text file.
        tf.write('SN%s: '%(t+1))

        true_params = lc[0].meta

        for key in sorted(true_params):
            tf.write('%s:%s ' % (key, true_params[key]))

        tf.write('\n')

        print 'True parameters for supernova %s saved in %s \n'\
               %((t + 1), true_file)
    tf.close

    return lcs

def get_lc(filelist):

    """ Reads lightcurves from files
        File names should be in a list
        Returns list of light curve tables
        Note - for this method you'll need to generate SN coords using
        simulate_sn_set() before using fit_snlc()"""

    filters.register_filters()

    lcs = []

    for i in range(len(filelist)):
        print 'Reading light curve of supernova %s \n' %(i + 1)
        lc = sncosmo.read_lc(filelist[i])
        lcs.append(lc)

    return lcs

def combine_scopes(parent_folder, f_1, f_2, f_3, nSNe):

    """ Method to merge multiple observations into one light curve file.
        - Opens two 'source' folders (f_1 and f_2) and iteratively combines
        the light curve files 'observed_lc_i' into a new light curve file,
        which is saved in f3.
    """
    ensure_dir(parent_folder+f_3)
    lc_list = []

    # Copy true_parameters file to new folder.
    copyfile(parent_folder + f_1 + 'true_parameters.txt', parent_folder + f_3 + 'true_parameters.txt')

    for t in range(nSNe):

        lc_1 = sncosmo.read_lc(parent_folder + f_1 + 'observed_lc_%s.txt'%(
            t+1))
        lc_2 = sncosmo.read_lc(parent_folder + f_2 + 'observed_lc_%s.txt'%(
            t+1))

        lc_3 = vstack([lc_1, lc_2])

        ensure_dir(parent_folder + f_3)

        sncosmo.write_lc(lc_3, parent_folder + f_3 + 'observed_lc_%s.txt'%(
            t+1))

        lc_list.append(lc_3)

    return lc_list

# C: FITTING ------------------------------------------------------------------


def fit_snlc(lightcurve, parent_folder, child_folder='TestFiles/', t0_in=0):

    """ Utility function for fitting lightcurves to
    observations of multiple SN """

    folder = parent_folder + child_folder
    ensure_dir(folder)

    # Create text file for storing 'fitted' parameters for each SN.
    fitted_file = folder + 'fitted_parameters.txt'
    error_file = folder + 'error_sn.txt'
    snr_error_file = folder + 'snr_error_sn.txt'
    accuracy_file = folder + 'fitted_errors.txt'
    valueError_file = folder + 'valueError_error_params.txt'
    ff = open(fitted_file, 'w')
    ef = open(error_file, 'w')
    snrf = open(snr_error_file, 'w')
    af = open(accuracy_file, 'w')
    vf = open(valueError_file, 'w')
    nSNe = len(lightcurve)

    # Get coordinates
    sn_dict = load_obj(parent_folder+'sn_dict')
    coords_in = sn_dict['Coordinates']

    # Get redshift
    params = sn_dict['Parameters']

    # Hold t0 outputs (only used for Kepler)
    if t0_in == 0:
        t0_in = [0]*nSNe

    explosion_time = []

    for i in range(nSNe):
        try:
            coords_out = [el[i] for el in coords_in]
            z = params[i]['z']
            p, fitted_t0, err = fit_util_lc(lightcurve[i], i + 1, folder, coords_out, z, t0_in[i])

            explosion_time.append(fitted_t0)

            # Write fitted parameters in text file.
            ff.write('SN%s: ' %(i+1))
            af.write('SN%s: ' %(i+1))

            for key in sorted(p):
                ff.write('%s:%s ' % (key, p[key]))

            for key in sorted(err):
                af.write('%s:%s ' % (key, err[key]))

            ff.write('\n')
            af.write('\n')

            print 'Fitted parameters for supernova %s saved in %s \n'\
		        %((i + 1), fitted_file)

        except RuntimeError, e:  # working around NaN thrown with garbage emcee fitter.
            print 'Error:',e

            # List skipped SN in error file
            ef.write('SN%s: \n' %(i+1))

            # Add fitted parameters as 0 for all (to fix indexing issue when using fitted t0 values)
            ff.write('SN%s: c:0 t0:0 x0:0 x1:0 z:0 \n' %(i+1))
            af.write('SN%s: c:0 t0:0 x0:0 x1:0\n' %(i+1))
            explosion_time.append(0)

            pass
        except ValueError, e:  # working around NaN thrown with garbage emcee fitter.
            print 'Error:',e

            coords_out = [el[i] for el in coords_in]
            z = params[i]['z']

            # List skipped SN in error file
            ef.write('SN%s: \n' %(i+1))

            # Add fitted parameters as 0 for all (to fix indexing issue when using fitted t0 values)
            ff.write('SN%s: c:0 t0:0 x0:0 x1:0 z:0 \n' %(i+1))
            af.write('SN%s: c:0 t0:0 x0:0 x1:0\n' %(i+1))
            explosion_time.append(0)

            # Write params to try track source of value error
            # Not 100% yet - probably needs to be inside util_lc instead.
            vf.write('SN%s: c:0 t0:0 x0:0 x1:0 z:%s \n' %(i+1, z))

            pass
        except sncosmo.fitting.DataQualityError, e:
            print 'Error:', e

            # List skipped SN in error file
            snrf.write('SN%s: \n' % (i + 1))

            # Add fitted parameters as 0 for all (to fix indexing issue when using fitted t0 values)
            ff.write('SN%s: c:0 t0:0 x0:0 x1:0 z:0 \n' % (i + 1))
            af.write('SN%s: c:0 t0:0 x0:0 x1:0\n' % (i + 1))
            explosion_time.append(0)


    ff.close()
    af.close()
    ef.close()
    snrf.close()
    vf.close()

    print 'Process complete.'

    return explosion_time


def fit_util_lc(data, index, folder, coords_in, z, t0_in):

    """ Fits SALT2 light curve to a supernova observation
        Plots model, data and residuals """

    print 'Fitting light curve for supernova %s' % index

    plotname = folder + 'fitted_lc_%s' % index

    ebv = dustmap.ebv(coords_in[0], coords_in[1], frame = 'fk5j2000', unit='degree')

    model.set(mwebv=ebv)

    # Hold Z
    model.set(z=z)


    if t0_in == 0:
        # Fitting SALT2 model using chisquared and MCMC
        # and fitting for t0.

        result_1, fitted_model_1 = chi_fit(data, model,
                                           ['t0', 'x0', 'x1', 'c'],
                                           minsnr=3.0,
                                           # Bounds are rough at this fit step, full prior introduced in mcmc fit
                                           bounds={'x0': (0.0000000001, 10), 'x1': (-3, 3), 'c': (-0.3, 0.3)}
                                           )

        result, fitted_model = mcmc_fit(data, fitted_model_1,
                                            # Parameters of model to vary.
                                            ['t0', 'x0', 'x1', 'c'],
                                            minsnr=3.0,
                                            priors={'x1': one_skew_prob_x1, 'c': one_skew_prob_c},
                                            nwalkers=10, nburn=200, nsamples=1000,
                                            guess_t0=False,
                                            guess_amplitude=False,
                                            )


    else:
        # Fitting SALT2 model using chisquared and MCMC
        # and setting starting guess for t0 manually (should be fitted t0 from Kepler)
        # print "t0_in:%s"%t0_in
        model.set(t0=t0_in)

        result_1, fitted_model_1 = chi_fit(data, model,
                                           ['x0', 'x1', 'c'],
                                           minsnr=3.0,
                                           guess_t0=False,
                                           # Bounds are rough at this fit step, full prior introduced in mcmc fit
                                           bounds={'x0': (0.0000000001, 10), 'x1': (-3, 3), 'c': (-0.3, 0.3)},
                                           )

        result, fitted_model = mcmc_fit(data, fitted_model_1,
                                            # Parameters of model to vary.
                                            ['t0', 'x0', 'x1', 'c'],
                                            minsnr=3.0,
                                            guess_t0=False,
                                            guess_amplitude=False,
                                            nwalkers=10, nburn=200, nsamples=1000,
                                            priors={'x1': one_skew_prob_x1, 'c': one_skew_prob_c}
                                            )

    fitted_t0 = result.parameters[1]
    # print "fitted t0:%s"%fitted_t0

    fitted_params = dict([(result.param_names[0], result.parameters[0]),
                          (result.param_names[1], result.parameters[1]),
                          (result.param_names[2], result.parameters[2]),
                          (result.param_names[3], result.parameters[3]),
                          (result.param_names[4], result.parameters[4])
                          ])

    fitted_errors = result.errors

    # Use PdfPages if saving as a pdf (ensure to change plotname)
    #pp = PdfPages(plotname)

    fig = sncosmo.plot_lc(data, format='png')

    plt.savefig(plotname + '_data.png')
    plt.close(fig)

    fig = sncosmo.plot_lc(data, model=fitted_model_1,
                          errors=result_1.errors, format='png',
                          )

    plt.savefig(plotname + '_chi.png')
    plt.close(fig)

    fig = sncosmo.plot_lc(data, model=fitted_model,
                    errors=result.errors, format='png'
                    )

    plt.savefig(plotname +'_mcmc.png')
    plt.close(fig)

    # pp.close()

    print 'Fitted light curve plotted in ' + plotname

    return fitted_params, fitted_t0, fitted_errors
