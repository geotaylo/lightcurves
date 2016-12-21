"""
GT 05/12/16
Aims to:
a) generate SALT2 observation light curves for SN with randomly 
   generated parameters and observations.
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

import filters


# CONSTANTS _----------------------------------------------------------


AREA = 5.7  # SkyMapper field of view.
Q0 = 0.2
LIGHT = 3*10**5  # km/s
H0 = 70.  # km/s/Mpc

smbands = ['smv', 'smg', 'smr', 'smi']  # SkyMapper bandpasses to use.
allbands = ['smv', 'smg', 'smr', 'smi', 'kst']  # SkyMapper + KST bands

dust = sncosmo.CCM89Dust()

# Change path to location of dustmaps
dustmap = sfdmap.SFDMap("/home/georgie/sfddata-master")

model = sncosmo.Model(source='salt2',
                      effects=[dust, dust],
                      effect_names=['host', 'mw'],
                      effect_frames=['rest', 'obs'])


# UTILITY FNS -----------------------------------------------------------------

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def ensure_dir(f):

    """ Checks if specified path exists, and creates it if not. """

    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def get_coords(nSNe):

    """ Generates random galactic coordinates for SN """

    # Galactic longitude.
    l = np.random.uniform(0, 360, nSNe).tolist()
    # Galactic latitude.
    b = np.random.uniform(-90, -30, nSNe).tolist()

    coords = [l, b]

    return coords


def mu(z):

    """ Distance modulus formula used to obtain x0. """

    d_L = LIGHT * (z + 0.5 * (1 - Q0) * z ** 2) / H0

    return 5*np.log10(d_L) + 25


def get_skynoise(sigma_psf, sigma_pixel):

    """ Calculates skynoise (background contribution to flux measurement
        error), assuming seeing-dominated gaussian psf with perfect psf
        photometry.
        sigma_psf = std. deviation of psf (pixels)
        sigma_pixel = background noise in single pixel (in counts)
    """

    return 4*math.pi*sigma_psf*sigma_pixel


# OBTAINING LC ----------------------------------------------------------------


def simulate_lc(nSNe=0, cadence=4, kpass=False, folder='TestFiles/',
                tmin=57754, tmax=58118, zmin=0.001, zmax=0.1,
                properties='', follow_up=False
               ):

    """ Generate SN parameters and simulate 'observed' light curve
        Defaults:  nSNe = 0 (number of SN)
                   cadence = 4 (days between observations)
                   kpass = False (if True, include kepler bandpass)
                   folder = 'Testfiles/' (path to store outputs)
                   tmin = 57754 (01/01/2017 in mjd)
                   tmax = 58118 (31/12/2017 in mjd)
                   zmin = 0.001 (min redshift detectable by SkyMapper)
                   zmax = 0.100 (max redshift detectable by SkyMapper)
                   properties = '' (path and name of dictionary of SN
                     parameters to use.  If blank, generates new params)
                   follow_up = False (if true, use 'good seeing'
                    observing properties
    """

    # Maintenance
    filters.register_filters()
    ensure_dir(folder)

    # Select filters.
    if kpass:
        bands = allbands
    else:
        bands = smbands

    lcs = []
    obs_out = []
    time = []
    n_obs = []
    zp_all = []
    gain = []
    skynoise = []

    # Load dictionary of SN properties
    if properties:
        sn_dict = load_obj(properties)
        params = sn_dict['Parameters']
        coords = sn_dict['Coordinates']
        obs_in = sn_dict['Observations']

        for n in range(nSNe):
            zp_sn = []
            time.append(obs_in[n][0])
            n_obs.append(obs_in[n][1])
            zp_sn.append(obs_in[n][2])
            gain = (obs_in[n][3])
            skynoise = (obs_in[n][4])
            s = set(obs_in[n][5])

            if 'kst' not in s:
                # Add kst zps to end
                zp_sn.append([25.47]*n_obs[n])
                zp_sn = [item for sublist in zp_sn for item in sublist]
            #flattening
            zp_all.append(zp_sn)

    # Generate SN properties
    else:
        # Set coordinates
        coords = get_coords(nSNe)

        # SN parameters
        params = []
        # Obtain redshifts
        if nSNe == 0:
            z = list(sncosmo.zdist(zmin, zmax,
                                  time=(tmax - tmin), area=AREA)
                                  )
            nSNe = len(z)
            print 'Obtained %s' %nSNe \
                  + ' redshifts from sncosmo distribution. \n'
        else:
            z = np.random.uniform(zmin, zmax, size=nSNe)
            print 'Obtained %s' %nSNe \
                  + ' redshifts from uniform distribution. \n'

        # Time at lightcurve peak (mjd).
        t0 = np.random.uniform(tmin, tmax, nSNe)
        c = 0.3*np.random.randn(nSNe)   # Colour
        x1 = 3*np.random.randn(nSNe)    # Stretch

        for i in range(nSNe):
            x0 = 10**((29.69 - mu(z[i]))/2.5)   # Amplitude (?)
            p = {'z': z[i], 't0': t0[i], 'x0': x0, 'x1': x1[i], 'c': c[i]}
            params.append(p)

        # Observations
        # ZP - distribuitons taken from SkYMapper obs paramteres
        for t in range(nSNe):
            tdet = np.random.randint(-15, -2)  # Time of init detection (mjd).
            tobs = np.random.randint(25, 65)  # Total observing period (days).
            tt = (params[t].get('t0')
                  + np.arange(tdet, tdet + tobs, cadence)).tolist()
            n_obs.append(len(tt))
            time.append(tt)
            gain = [1.]*len(bands)
            sn = [100.]*(len(bands)*len(tt))  # Filler value for skynoise

            if follow_up:
                zp_g = (np.random.normal(26.87, 0.68, len(tt)))
                zp_i = (np.random.normal(25.85, 0.81, len(tt)))
                zp_r = (np.random.normal(26.63, 0.67, len(tt)))
                zp_v = (np.random.normal(24.91, 0.70, len(tt)))
            else:
                zp_g = (np.random.normal(26.82, 0.79, len(tt)))
                zp_i = (np.random.normal(25.21, 0.36, len(tt)))
                zp_r = (np.random.normal(26.71, 0.76, len(tt)))
                # No v in bad seeing
                zp_v = (np.random.normal(24.91, 0.70, len(tt)))
            zp_k = [25.47]*len(tt)
            zp_sn = [zp_v, zp_g, zp_r, zp_i]
            # order is important
            if kpass:
                zp_sn.append(zp_k)
                sn[-len(tt):] = [0.]*len(tt)
            zp_sn = [item for sublist in zp_sn for item in sublist]
            zp_all.append(zp_sn)
            skynoise.append(sn)

    true_file = folder + 'true_parameters.txt'
    tf = open(true_file, 'w')

    for t in range(nSNe):
        o_t = []
        observing_bands = []

        # Adds 12s offset time between observations in different filters
        for x in range(len(bands)):
            j = np.array(time[t]) + 0.00013888888*x
            o_t.append(j.tolist())
            a = n_obs[t]*[bands[x]]
            observing_bands.extend(a)
        observing_time = [item for sublist in o_t for item in sublist]

        n_points = n_obs[t]*len(bands)

        # List of observing properties for each SN, to return
        obs_util = [time[t], n_obs[t], zp_all[t], gain, skynoise, bands]
        observing_dictionary = {'band': observing_bands,
                                'time': observing_time,
                                'zp': zp_all[t],
                                'zpsys': n_points*['ab'],
                                'gain': n_points*[gain[t]],
                                'skynoise': skynoise[t]
                                }

        obs_out.append(obs_util)
        obs = Table(observing_dictionary)

        lc = sncosmo.realize_lcs(obs, model, [params[t]])
        lcs.append(lc[0])

        name = folder + 'observed_lc_%s.txt'%(t + 1)

        # Write observations.
        sncosmo.write_lc(lc[0], name)

        if kpass:
            print 'Simulated observations for supernova %s saved in %s;\
                   Kepler filter included'%((t + 1), name)
        else:
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
    dict_out = {'Parameters': params, 'Coordinates': coords,
                 'Observations': obs_out}
    save_obj(dict_out, folder+'sn_dict')
    print 'Properties of simulation saved in %ssn_dict.pkl'%folder

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

    ebv = dustmap.ebv(coords_in[0], coords_in[1], frame='galactic',
                      unit='degree')

    model.set(mwebv=ebv)

    # Fitting SALT2 model using MCMC.
    result, fitted_model = fit(data, model,
                               # Parameters of model to vary.
                               ['z', 't0', 'x0', 'x1', 'c'],
                               bounds={'z': (0.001, 0.1)}, minsnr=2.0
                               )

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


def get_diff(folder):

    """ Reads text files of fitted and true parameters, and calculates
        diffences for each SN parameter """

    fp = folder + 'fitted_parameters.txt'
    tp = folder + 'true_parameters.txt'

    fitted_c = []
    fitted_t0 = []
    fitted_x0 = []
    fitted_x1 = []
    fitted_z = []

    with open(fp, 'rb') as f:

        reader = csv.reader(f, delimiter=' ')

        for row in reader:

            fitted_c.append(float(row[1].strip('c:')))
            fitted_t0.append(float(row[2].strip('t0:')))
            fitted_x0.append(float(row[3].strip('x0:')))
            fitted_x1.append(float(row[4].strip('x1:')))
            fitted_z.append(float(row[5].strip('z:')))

    #fitted_params = [fitted_c, fitted_t0, fitted_x0, fitted_x1, fitted_z]

    true_c = []
    true_t0 = []
    true_x0 = []
    true_x1 = []
    true_z = []

    with open(tp, 'rb') as file:

        reader = csv.reader(file, delimiter=' ')

        for row in reader:

            true_c.append(float(row[1].strip('c:')))
            true_t0.append(float(row[2].strip('t0:')))
            true_x0.append(float(row[3].strip('x0:')))
            true_x1.append(float(row[4].strip('x1:')))
            true_z.append(float(row[5].strip('z:')))

    #true_params = [true_c, true_t0, true_x0, true_x1, true_z]

    c_diff = 0
    t0_diff = 0
    x0_diff = 0
    x1_diff = 0
    z_diff = 0

    plotname = folder + 'differences.png'
    f, axarr = plt.subplots(5, sharex=True)
    axarr[0].plot(true_c, 'o')
    axarr[0].plot(fitted_c, 'v')
    axarr[0].set_ylabel('c')
    axarr[0].set_title('Fitted VS True Parameters')
    axarr[1].plot(true_t0, 'o')
    axarr[1].plot(fitted_t0, 'v')
    axarr[1].set_ylabel('t0')
    axarr[2].plot(true_x0, 'o')
    axarr[2].plot(fitted_x0, 'v')
    axarr[2].set_ylabel('x0')
    axarr[3].plot(true_x1, 'o')
    axarr[3].plot(fitted_x1, 'v')
    axarr[3].set_ylabel('x1')
    axarr[4].plot(true_z, 'o')
    axarr[4].plot(fitted_z, 'v')
    axarr[4].set_ylabel('z')
    f.savefig(plotname)

    for i in range(len(true_c)):
        c_diff = c_diff + abs(true_c[i] - fitted_c[i])
        t0_diff = t0_diff + abs(true_t0[i] - fitted_t0[i])
        x0_diff = x0_diff + abs(true_x0[i] - fitted_x0[i])
        x1_diff = x1_diff + abs(true_x1[i] - fitted_x1[i])
        z_diff = z_diff + abs(true_z[i] - fitted_z[i])

    return
