import pandas as pd
import matplotlib.pyplot as plt
import os
import math
from scipy.optimize import curve_fit
import pylab as py
import scipy.stats as stats
import numpy as np
from numpy.polynomial import Polynomial
from itertools import chain
import pickle

# Import and plot histograms
df = pd.read_csv('info_Georgie_161209.list', delimiter=' ')

# mjd, filter, exptime, sp.field.id, sp.field.ra, sp.field.dec, sp.airmass, sp.seeing, si_bkgsig, si_maglim, nimg
df2 = pd.read_csv('survey_hist_bad_seeing_1216.txt', header=None, delim_whitespace=True, usecols=[1, 7, 8], names=['filters', 'seeing', 'bkgsig'])

df3 = pd.read_csv('survey_hist_good_seeing_1216.txt', header=None, delim_whitespace=True, usecols=[1, 7, 8], names=['filters', 'seeing', 'bkgsig'])

filters1 = df.filters.tolist()
seeing1 = df.seeing.tolist()
bkgsig1 = df.bkgsig.tolist()

filters2 = df2.filters.tolist()
seeing2 = df2.seeing.tolist()
bkgsig2 = df2.bkgsig.tolist()

filters3 = df3.filters.tolist()
seeing3 = df3.seeing.tolist()
bkgsig3 = df3.bkgsig.tolist()

filters = filters1 + filters2 + filters3
seeing = seeing1 + seeing2 + seeing3
bkgsig = bkgsig1 + bkgsig2 + bkgsig3

full_seeing = [i for i in seeing if not math.isnan(i)]
full_bkgsig = [i for i in bkgsig if not math.isnan(i)]
all_seeing = [i for i in full_seeing if i <= 25]
all_bkgsig = [i for i in full_bkgsig if i <= 100]

folder = 'skynoise/onefilter/'
if not os.path.isdir(folder):
    os.makedirs(folder)

weibull_params = []
figa, ax1 = plt.subplots(1)
data1 = ax1.hist(all_seeing, bins=100, normed=1)
x_fit = py.linspace(0, max(all_seeing), 200)
params = stats.exponweib.fit(all_seeing, loc=0)
weibull_params.append('all_seeing: %s, %s, %s, %s \n'%(params[0], params[1], params[2], params[3]))
ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')

# # Plot the results
ax1.set_title('Distribution of seeing in the SkyMapper filters')
plt.xlabel('Seeing')
plt.ylabel('Frequency')
plt.legend()
figa.savefig(folder+'all_seeing.png')
plt.close(figa)


figa, ax1 = plt.subplots(1)
data2 = ax1.hist(all_bkgsig, bins=100, normed=1)
x_fit = py.linspace(0, 100, 200)
params = stats.exponweib.fit(all_bkgsig, loc=0)
weibull_params.append('all_bkgsig: %s, %s, %s, %s \n'%(params[0], params[1], params[2], params[3]))

ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')

ax1.set_title('Distribution of bkgsig in the SkyMapper filters')
plt.xlabel('bkgsig')
plt.ylabel('Frequency')
figa.savefig(folder+'all_bkgsig.png')
plt.close(figa)


with open(folder + 'weibull_params.txt', 'w') as file:
    for item in weibull_params:
        file.write("%s" % item)
