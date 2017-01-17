"""

GT 05/12/16

Aims to:
a) generate SALT2 observation light curves for SN with randomly generated
   parameters and observations.
b) take an observation light curve and fit a SALT2 model to it.

Requires installation of:
- SNCosmo:  https://sncosmo.readthedocs.io/en/v1.4.x/index.html
- Emcee:
- sfdmap:

"""


import os
import csv
import pickle
import math
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages

import sfdmap
import sncosmo
from sncosmo import mcmc_lc as fit
#from tabulate import tabulate

# filters.py should be in working directory
import filters


# CONSTANTS -------------------------------------------------------------------

# SkyMapper observing cadence (days)
cad_sm = 4.


# Kepler observing cadence (30 minutes, in days)
cad_k = 1. / 48.


# SkyMapper field of view (square degrees)
AREA = 5.7


# Skymapper min and max observable redshift
zmin = 0.001
zmax = 0.1


# Possible observing period: all of 2017 (mjd)
tmin = 57754
tmax = 58118


Q0 = 0.2


# Speed of light (km/s)
LIGHT = 2.997*10**5


# (km/s/Mpc)
H0 = 70.00


# DUSTMAPS --------------------------------------------------------------------
dust = sncosmo.CCM89Dust()


# Change path to location of dustmaps
dustmap = sfdmap.SFDMap("/home/georgie/sfddata-master")


# SALT2 MODEL TEMPLATE --------------------------------------------------------
model = sncosmo.Model(source='salt2',
                      effects=[dust, dust],
                      effect_names=['host', 'mw'],
                      effect_frames=['rest', 'obs'])


# UTILITY FNS -----------------------------------------------------------------

def save_obj(obj, name):
    """ Utility function for saving dictionary. """

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


def get_coords(nSNe):
    """ Generates random galactic coordinates for SN """

    # Right ascension
    ra = np.random.uniform(0, 360, nSNe).tolist()

    # Declination (within SkyMapper bounds)
    dec = np.random.uniform(-90, 10, nSNe).tolist()

    coords = [ra, dec]

    return coords


def mu(z):
    """ Distance modulus formula used to obtain x0. """

    d_L = LIGHT * (z + 0.5 * (1 - Q0) * z ** 2) / H0

    return (5*np.log10(d_L) + 25)


def get_skynoise(filter, seeing):
    """ Calculates static skynoise (background contribution to flux measurement
        error), assuming seeing-dominated gaussian psf with perfect psf
        photometry.
        sigma_psf = std. deviation of psf (pixels)
        sigma_pixel = background noise in single pixel (in counts)
    """

    if filter == 'smg':
        if seeing == 'good':
            sigma_psf = 1.57
            sigma_pixel = 24.08
        else:
            sigma_psf = 2.03
            sigma_pixel = 22.79

    if filter == 'smi':
        if seeing == 'good':
            sigma_psf = 1.57
            sigma_pixel = 17.32
        else:
            sigma_psf = 2.28
            sigma_pixel = 13.54

    if filter == 'smr':
        if seeing == 'good':
            sigma_psf = 2.17
            sigma_pixel = 8.45
        else:
            sigma_psf = 1.95
            sigma_pixel = 1.25

    if filter == 'smv':
        sigma_psf = 1.96
        sigma_pixel = 6.23

    return 4*math.pi*sigma_psf*sigma_pixel


# GENERATING SIMULATION SET ---------------------------------------------------


def simulate_sn_set(folder, nSNe=0):

    """ Generates full set of supernovae and observation parameters for all
    telescopes and filters, to be used by simulate_lc"""

    ensure_dir(folder)

    params = []
    observations = []

    # Supernova parameters ----------------------------------------------

    # z = redshift
    if nSNe == 0:
        z = list(sncosmo.zdist(zmin, zmax, time=(tmax - tmin), area=AREA))
        nSNe = len(z)
        print 'Obtained %s' % nSNe + ' redshifts from sncosmo distribution. \n'
    else:
        z = np.random.uniform(zmin, zmax, size=nSNe)
        print 'Obtained %s' % nSNe + ' redshifts from uniform distribution. \n'

    # t0 = time at lightcurve peak (mjd).
    t0 = np.random.uniform(tmin, tmax, nSNe)

    # c = colour
    c = 0.3 * np.random.randn(nSNe)

    # x1 = stretch
    x1 = 3 * np.random.randn(nSNe)

    # x0 = amplitude (I think...)
    for i in range(nSNe):
        x0 = 10 ** ((29.69 - mu(z[i])) / 2.5)
        p = {'z': z[i], 't0': t0[i], 'x0': x0, 'x1': x1[i], 'c': c[i]}
        params.append(p)

    coords = get_coords(nSNe)

    # Observing parameters ----------------------------------------------

    t_det_sm = []
    t_det_k = []

    t_obs_sm = []
    t_obs_k = []

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

    for t in range(nSNe):

        # SKYMAPPER --------------------------
        # Time of init detection (mjd)
        td_sm = np.random.randint(-15, -2)
        t_det_sm.append(td_sm)

        # Total observing period (days)
        to_sm = np.random.randint(25, 65)
        t_obs_sm.append(to_sm)

        # Observation times
        t_sm = (params[t].get('t0')
                + np.arange(td_sm, td_sm + to_sm, cad_sm)).tolist()
        n_obs_sm.append(len(t_sm))
        time_sm.append(t_sm)

        # Zero points
        # For standard filter set (g, r, i) - used in 'bad seeing'
        zp_g_bad = (np.random.normal(26.82, 0.79, len(t_sm)))
        zp_i_bad = (np.random.normal(25.21, 0.36, len(t_sm)))
        zp_r_bad = (np.random.normal(26.71, 0.76, len(t_sm)))
        zp_gri = [zp_g_bad, zp_i_bad, zp_r_bad]

        # For extended filter set (g, r, i, v) - used in 'good seeing'
        zp_g_good = (np.random.normal(26.87, 0.68, len(t_sm)))
        zp_i_good = (np.random.normal(25.85, 0.81, len(t_sm)))
        zp_r_good = (np.random.normal(26.63, 0.67, len(t_sm)))
        zp_v_good = (np.random.normal(24.91, 0.70, len(t_sm)))
        zp_griv = [zp_g_good, zp_i_good, zp_r_good, zp_v_good]

        # Skynoise
        # For standard filter set (g, r, i) - used in 'bad seeing'
        sn_g_bad = [get_skynoise('smg', 'bad')] * len(t_sm)
        sn_i_bad = [get_skynoise('smi', 'bad')] * len(t_sm)
        sn_r_bad = [get_skynoise('smr', 'bad')] * len(t_sm)
        sn_gri = [sn_g_bad, sn_i_bad, sn_r_bad]

        # For extended filter set (g, r, i, v) - used in 'good seeing'
        sn_g_good = [get_skynoise('smg', 'good')] * len(t_sm)
        sn_i_good = [get_skynoise('smi', 'good')] * len(t_sm)
        sn_r_good = [get_skynoise('smr', 'good')] * len(t_sm)
        sn_v_good = [get_skynoise('smv', 'good')] * len(t_sm)
        sn_griv = [sn_g_good, sn_i_good, sn_r_good, sn_v_good]


        # KEPLER----------------------------------
        # Time of init detection (mjd)
        td_k = np.random.randint(-15, -2)
        t_det_k.append(td_k)

        # Total observing period (days)
        to_k = np.random.randint(80, 85)
        t_obs_k.append(to_k)

        # Observation times
        t_k = (params[t].get('t0')
               + np.arange(td_k, td_k + to_k, cad_k)).tolist()
        n_obs_k.append(len(t_k))
        time_k.append(t_k)

        # Zero points
        zp_k = [25.47] * len(t_k)

        # Skynoise
        sn_k = [0.] * len(t_k)

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
    return


# OBTAINING LC ----------------------------------------------------------------


def simulate_lc(parent_folder, child_folder='TestFiles/',
                scope='sm', follow_up=False
                ):
    """ Generate SN parameters and simulate 'observed' light curve
        Defaults:  scope = 'sm' - which telescope to use in simulations.
                            Accepted values are 'sm', 'kst' or 'both'.
                   folder = 'Testfiles/' (path to store outputs)
                   follow_up = False (if true, use 'good seeing'
                    observing properties
    """

    # Maintenance.
    child_folder = parent_folder + child_folder
    filters.register_filters()
    ensure_dir(child_folder)
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
        if follow_up:
            bands = ['smg', 'smr', 'smi', 'smv']
        else:
            bands = ['smg', 'smr', 'smi']
    else:
        if scope == 'kst':
            bands = ['kst']
        elif scope == 'both':
            if follow_up:
                bands = ['smg', 'smr', 'smi', 'smv', 'kst']
            else:
                bands = ['smg', 'smr', 'smi', 'kst']
        else:
            print 'Hey buddy, that\'s not a real telescope.  Please enter ' \
                  'scope = \'sm\' for SkyMapper, \'kst\' for Kepler, ' \
                  'or \'both\' for both.  Thanks.'

    # Load properties of simulation set
    sn_dict = load_obj(parent_folder+'sn_dict')
    params = sn_dict['Parameters']
    obs_in = sn_dict['Observations']
    nSNe = sn_dict['nSNe']

    for n in range(nSNe):
        time_sm.append(obs_in[n][4])
        time_k.append(obs_in[n][5])
        n_obs_sm.append(obs_in[n][6])
        n_obs_k.append(obs_in[n][7])
        zp_gri.append(obs_in[n][8])
        zp_griv.append(obs_in[n][9])
        zp_k.append(obs_in[n][10])
        sn_gri.append(obs_in[n][11])
        sn_griv.append(obs_in[n][12])
        sn_k.append(obs_in[n][13])

    # RUNNING SIMULATION --------------------------------------

    true_file = child_folder + 'true_parameters.txt'
    tf = open(true_file, 'w')

    for t in range(nSNe):

        # time of observations for all filters (mjd)
        o_t = []
        observing_bands = []

        if scope == 'sm':

            # Adds 12s offset time between observations in different sm filters
            # (for slewing)
            for x in range(len(bands)):
                j = np.array(time_sm[t]) + 0.00013888888*x
                o_t.append(j.tolist())
                a = n_obs_sm[t][0]*[bands[x]]
                observing_bands.extend(a)
            n_points = n_obs_sm[t][0] * len(bands)

            # Sets zp
            if follow_up:
                zp = zp_griv[t]
            else:
                zp = zp_gri[t]

            # Sets skynoise
            if follow_up:
                skynoise = sn_griv[t]
            else:
                skynoise = zp_gri[t]

        elif scope == 'kst':

            o_t.append(time_k[t])
            k_hold = n_obs_k[t][0] * ['kst']
            observing_bands.extend(k_hold)
            n_points = n_obs_k[t][0]

            # Sets zp
            zp = [zp_k[t]]

            # Sets skynoise
            skynoise = [sn_k[t]]

        else:

            # Adds 12s offset time between observations in different sm filters
            # (for slewing)
            for x in range(len(bands)-1):
                j = np.array(time_sm[t]) + 0.00013888888*x
                o_t.append(j.tolist())
                a = n_obs_sm[t][0]*[bands[x]]
                observing_bands.extend(a)

            # Adds kepler observation times to set
            o_t.append(time_k[t])
            k_hold = n_obs_k[t][0]*['kst']
            observing_bands.extend(k_hold)
            n_points = (n_obs_sm[t][0] * (len(bands) - 1)) + n_obs_k[t][0]

            # Sets zp
            if follow_up:
                zp = zp_griv[t]
            else:
                zp = zp_gri[t]
            zp.append(zp_k[t])

            # Sets skynoise
            if follow_up:
                skynoise = sn_griv[t]
            else:
                skynoise = sn_gri[t]
            skynoise.append(sn_k[t])


        # Flattening lists of lists
        observing_time = [item for sublist in o_t for item in sublist]
        observing_time = [item for sublist in observing_time for item in sublist]
        skynoise = [item for sublist in skynoise for item in sublist]
        zp = [item for sublist in zp for item in sublist]

        print len(observing_bands)
        print len(observing_time)
        print observing_time
        print zp
        print skynoise
        print len(zp)
        print len(skynoise)
        print n_points

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
        Returns list of light curve tables """

    filters.register_filters()

    lcs = []

    for i in range(len(filelist)):
        print 'Reading light curve of supernova %s \n' %(i + 1)
        lc = sncosmo.read_lc(filelist[i])
        lcs.append(lc)

    return lcs


# FITTING ---------------------------------------------------------------------


def fit_snlc(lightcurve, folder='TestFiles/', properties=''):

    """ Utility function for fitting lightcurves to
    observations of multiple SN """

    ensure_dir(folder)

    # Create text file for storing 'fitted' parameters for each SN.
    fitted_file = folder + 'fitted_parameters.txt'
    ff = open(fitted_file, 'w')
    nSNe = len(lightcurve)
    # Get coordinates
    if properties:
        sn_dict = load_obj(properties)
        coords_in = sn_dict['Coordinates']
    else:
        coords_in = get_coords(nSNe)

    for i in range(nSNe):
        coords_out = [el[i] for el in coords_in]
        p = fit_util_lc(lightcurve[i], i + 1, folder, coords_out)

        # Write fitted parameters in text file.
        ff.write('SN%s: ' %(i+1))

        for key in sorted(p):
            ff.write('%s:%s ' % (key, p[key]))

        ff.write('\n')

        print 'Fitted parameters for supernova %s saved in %s \n'\
            %((i + 1), fitted_file)

    ff.close()

    print 'Process complete.'

    return


def fit_util_lc(data, index, folder, coords_in):

    """ Fits SALT2 light curve to a supernova observation
        Plots model, data and residuals """

    print 'Fitting light curve for supernova %s' %index

    plotname = folder + 'fitted_lc_%s.pdf' %index

    ebv = dustmap.ebv(coords_in[0], coords_in[1])

    model.set(mwebv=ebv)

    # Fitting SALT2 model using MCMC.
    """result, fitted_model = fit(data, model,
                               # Parameters of model to vary.
                               ['z', 't0', 'x0', 'x1', 'c'],
                               bounds={'z': (0.001, 0.1)}, minsnr=3.0
                               )"""
    result, fitted_model = sncosmo.fit_lc(data, model,
            ['z', 't0', 'x0', 'x1', 'c'],  # parameters of model to vary
            bounds={'z': (0.001, 0.1)})

    fitted_params = dict([(result.param_names[0], result.parameters[0]),
                          (result.param_names[1], result.parameters[1]),
                          (result.param_names[2], result.parameters[2]),
                          (result.param_names[3], result.parameters[3]),
                          (result.param_names[4], result.parameters[4])
                          ])

    pp = PdfPages(plotname)

    sncosmo.plot_lc(data, model=fitted_model,
                    errors=result.errors, fname=pp, format='pdf'
                    )

    pp.close()

    print 'Fitted light curve plotted in ' + plotname

    return fitted_params