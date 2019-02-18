import matplotlib.pyplot as plt

import scipy.stats as st
import numpy as np
from numpy import random
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from astropy.cosmology import FlatLambdaCDM
from astropy.extern.six.moves import range
import math


# For z distributions
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

    for i in range(random.poisson(nsim)):
        yield float(snrate_ppf(random.random()))

# z = list(modified_zdist(0.001, 0.15, 1000))
# print z
# plt.hist(z)
# plt.show()


def uniform_sample(xmin, xmax, ymin, ymax):
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


# class my_pdf(st.rv_continuous):
#     def _pdf(self,x):
#         max_prob = -0.055
#         sigma_m = 0.023
#         sigma_p = 0.150
#         numer = -((x - max_prob) ** 2)
#         if x <= max_prob:
#             return math.exp(numer/(2*sigma_m**2))
#         else:
#             return math.exp(numer/(2*sigma_p**2))

# my_cv = my_pdf(a=-1., b=1., name='my_pdf')
#
# cs = my_cv.rvs(size=1000)

# probs = np.empty(1000)
# vals = np.empty(1000)
# for i in range(1000):
#     probs[i], vals[i] = triangular(-0.055, 0.023, 0.150)

NormSamples = []
for i in range(1000):
    NormSamples.append(one_skew_rv(0.436, 3.118, 0.724))
# print probs, vals

# plt.scatter(vals, probs)
# # plt.xlim(-0.2,0.2)
# plt.show()

plt.hist(NormSamples, bins=50)
plt.xlim(-2, 2)
plt.show()