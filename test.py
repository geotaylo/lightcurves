import matplotlib.pyplot as plt

from scipy.stats.distributions import norm
import numpy as np
from numpy import random
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from astropy.cosmology import FlatLambdaCDM
from astropy.extern.six.moves import range
import vg
import math
import sncosmo
from numpy import trapz

lcdist = 'ps1_g10'
Q0 = -0.5
# Speed of light (km/s)
LIGHT = 2.997*10**5
# Hubble Constant (km/s/Mpc)
H0 = 70.00

model = sncosmo.Model(source='salt2')

# For z distributions - fully tested and implemented
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

    for i in range(nsim):
        yield float(snrate_ppf(random.random()))

# z = list(modified_zdist(0.001, 0.15, 1000))
# print z
# plt.hist(z)
# plt.show()


# Tested and implemented
def uniform_sample(xmin, xmax, ymin, ymax):
    return np.random.random(2,)*np.array([xmax-xmin,ymax-ymin])+np.array([xmin, ymin])

# Tested and implemented
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
    """ Generate a probability for obtaining a certain param c or x1, as defined in Scolnic and Kessler 2016.

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
    if x <= max_prob:
        prob =  (math.exp(numer/(2*sigma_m**2)))/normfactor
    else:
        prob =  (math.exp(numer/(2*sigma_p**2)))/normfactor
    return prob

def one_skew_prob_x1(param):
    """ Generate a probability for obtaining a certain param c or x1, as defined in Scolnic and Kessler 2016.

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
    # Low-z x1 condition:
    if x <= -3. or x >= 3.:
        prob = 0.000000000001
    elif x <= max_prob:
        prob =  (math.exp(numer/(2*sigma_m**2)))/normfactor
    else:
        prob =  (math.exp(numer/(2*sigma_p**2)))/normfactor
    return prob


# NormSamples = []
# for i in range(1000):
#     NormSamples.append(one_skew_rv(0.419, 3.024, 0.742))
#
fig, ax = plt.subplots()
# #norm_1, bins_1, patches_1 = plt.hist(NormSamples, bins=50, alpha = 0.4)
# norm, bins, patches = ax.hist(NormSamples, bins=50, alpha = 0.4, density=True)
# plt.xlim(-5, 5)
# # plt.show()
# x_sample = np.linspace(-20,20,1000)
# x_probs = []
# for i in range(len(x_sample)):
#     p = one_skew_prob_x1(x_sample[i])
#     x_probs.append(p)
#     print x_sample[i]
#     print p
#     logtest = math.log(p)
#
# # ax.scatter(x_sample, np.array(x_probs), color='k')
# #
# # plt.show()

zs = list(modified_zdist(0.001, 0.15, 100))

def mu(z):
    """ Distance modulus formula used to obtain x0. """

    d_L = LIGHT * (z + 0.5 * (1 - Q0) * z ** 2) / H0

    return 5*np.log10(d_L) + 25

xs = []
for i in range(100):
    x0 = 10 ** ((29.69 - mu(zs[i])) / 2.5)
    xs.append(x0)

ax.scatter(zs, xs, color='k')

mod_xs = []
for z in zs:
    # Absoute B-band magnitude of SN, with some scatter introduced.
    mabs = np.random.normal(-19.3, 0.3)
    model.set(z=z)
    model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    x0 = model.get('x0')
    mod_xs.append(x0)

ax.scatter(zs, mod_xs, color='m')
plt.show()