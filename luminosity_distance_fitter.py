""" A program to fit luminosity distances from light curve parameters."""

import csv
import pylab as py
import numpy as np
from scipy.optimize import curve_fit
import copy

# best fit values from Betoule paper
alpha = 0.141 # +- 0.006 (stat+sys)
beta = 3.101 # +- 0.075 (stat+sys)
M_1_B = -19.05 # +- 0.02 (stat+sys)
Delta_M = -0.07 # +- 0.023 (stat+sys)

# Hubble constant
H0 = 70.0

# Stellar mass - deal with this later
M_stellar = 1

# Cosmology - not clear on where these came from just yet, or why we need 'em.
my_Om0 = 0.3
Ode0 = 0.7
my_w0 = -1.0
wa = 0.0

def get_params(folder, set, fails=[]):
    """
    Calculates residuals (true - fitted value) of SN parameters for one filter set.

    :param folder: path to specific filter set data (e.g. '~\kst\')
    :param set: name of filter set, for labelling purposes (e.g. 'kst')
    :param fails: associated failed SN from another filter set.  For example, if kst has failed to fit a lightcurve to
                    SN10, the combined sm+kst fit will not be able to use kst's fitted t0 value as an initial guess,
                    and the fit will be affected.  SN10 should be counted as a 'fail' for sm+kst, as the fit has
                    underperformed.
    :param wipe_fails: if True, the 'failed' SN will not be included in the output of residuals, showing the best
                        possible performance of the code.  If False, the entire SN set is included, regardless of fit
                        success.
    :return: error set (list of 'failed' SN) and diffs, a table of absolute residuals.
    """

    # Check error SN and flag:
    error_file = folder + 'error_sn.txt'
    error_set = []
    f = open(error_file, "r")
    lines = f.readlines()

    for x in lines:
        error_set.append(int(x.strip('SN : \n')))
    f.close()

    # Import and strip fitted params
    fp = folder + 'fitted_parameters.txt'

    fitted_c = []
    fitted_t0 = []
    fitted_x0 = []
    fitted_x1 = []
    fitted_z = []

    with open(fp, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            fitted_c.append(float(row[1].replace('c:','')))
            fitted_t0.append(float(row[2].replace('t0:','')))
            fitted_x0.append(float(row[3].replace('x0:','')))
            fitted_x1.append(float(row[4].replace('x1:','')))
            fitted_z.append(float(row[5].replace('z:','')))

    # Double check for missed errors! - need to figure out why these are getting skipped
    for i in range(len(fitted_c)):
        if fitted_c[i] == 0 and fitted_t0[i] == 0 and fitted_x0[i] == 0 and fitted_x1[i] ==0 and fitted_z[i] ==0:
            if i+1 not in error_set:
                error_set.append(i+1)
    error_set.sort()

    # Create sn_num array
    sn_num = range(1, len(fitted_c)+1)
    for i in error_set:
        sn_num[i-1] = 'fit_error'+ str(i)

    # Flag kepler errors (only applicable for combined seeing where kst t0 couldn't be passed)
    for i in fails:
        sn_num[i-1] = 'kst_error' + str(i)

    total_fails = filter(lambda x:x in error_set, fails)

    # Handles failures in both kst and current fit
    for i in total_fails:
        sn_num[i-1] = 'fit_and_kst_error' + str(i)


    # remove fails from data
    for i in sorted(error_set+fails, reverse=True):
        del sn_num[i-1]
        del fitted_c[i-1]
        del fitted_t0[i-1]
        del fitted_x0[i-1]
        del fitted_x1[i-1]
        del fitted_z[i-1]

    params = [fitted_c, fitted_t0, fitted_x0, fitted_x1, fitted_z, sn_num]

    return error_set, params


def get_true(folder):
    """
    Returns true SN values for comparison with fitted params.

    :param folder: path to specific filter set data (e.g. '~\kst\') - SN will always be the same so only need to run on
    one - just use kst for now.

    :return: parameters
    """

    # Import and strip true params
    tp = folder + 'true_parameters.txt'

    true_c = []
    true_t0 = []
    true_x0 = []
    true_x1 = []
    true_z = []

    with open(tp, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            true_c.append(float(row[1].replace('c:','')))
            true_t0.append(float(row[2].replace('t0:','')))
            true_x0.append(float(row[3].replace('x0:','')))
            true_x1.append(float(row[4].replace('x1:','')))
            true_z.append(float(row[5].replace('z:','')))

    # Create sn_num array
    sn_num = range(1, len(true_c)+1)

    params = [true_c, true_t0, true_x0, true_x1, true_z, sn_num]

    return params



# def absolute_mag(M_stellar):
#     """ Calculates a 'nuisance parameter', absolute magnitude M_B, depending on stellar mass M_stellar
#         Reference: C11 (Conley et al 2011)"""
#     if M_stellar <= sun_factor:
#         M_b = M_1_B
#     else:
#         M_b = M_1_B + Delta_M
#     return M_b


def apparent_mag(x0):
    """ Calculates apparent peak b band magnitude (rest frame) from SALT2 parameter x0"""

    # Diemer et al 2013 - source for equation!
    m_b_star = -2.5*np.log10(x0)-10.095

    return m_b_star


def distance_modulus(x0, x1, c):
    """
    :param M_stellar: Mass of star (for calculating absolute magnitude
    :param x0: scaling parameter (units of flux?)
    :param x1: time stretching of LC
    :param c: SN colour at peak brightness, from LC fit
    :return: distance modulus, mu
    """

    # Absolute magnitude of SN
    # M_b = absolute_mag(M_stellar)
    M_b = M_1_B

    # Apparent magnitude of SN
    m_b_star = apparent_mag(x0)

    # distance modulus
    mu = m_b_star - M_b + (alpha*x1)-(beta*c)

    return mu


def distance(mu):
    """ Returns luminosity distance in Mpc """
    y = (mu -25)/5
    d_l = 10**(y)
    return d_l
