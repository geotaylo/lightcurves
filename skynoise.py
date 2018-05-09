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

g_seeing = []
i_seeing = []
v_seeing = []
r_seeing = []

g_bkgsig = []
i_bkgsig = []
v_bkgsig = []
r_bkgsig = []

# good_g_seeing = []
# good_i_seeing = []
# good_v_seeing = []
# good_r_seeing = []
#
# good_g_bkgsig = []
# good_i_bkgsig = []
# good_v_bkgsig = []
# good_r_bkgsig = []

def f(x, a, b, c):
    return a * py.exp(-(x - b) ** 2.0 / (2 * c ** 2))

def f_parab(x, a, b, c):
    return a * x**2 + b * x + c

# # Sort
# for i in range(len(seeing)):
#     # Bad seeing
#     if not math.isnan(bkgsig[i]):
#         if seeing[i] >= 2.35:
#             if filters[i] == 'g':
#                 bad_g_seeing.append(seeing[i])
#                 bad_g_bkgsig.append(bkgsig[i])
#             elif filters[i] == 'i':
#                 bad_i_seeing.append(seeing[i])
#                 bad_i_bkgsig.append(bkgsig[i])
#             elif filters[i] == 'v':
#                 bad_v_seeing.append(seeing[i])
#                 bad_v_bkgsig.append(bkgsig[i])
#             elif filters[i] == 'r':
#                 bad_r_seeing.append(seeing[i])
#                 bad_r_bkgsig.append(bkgsig[i])
#         elif seeing[i] <= 2.35:
#             if filters[i] == 'g':
#                 good_g_seeing.append(seeing[i])
#                 good_g_bkgsig.append(bkgsig[i])
#             elif filters[i] == 'i':
#                 good_i_seeing.append(seeing[i])
#                 good_i_bkgsig.append(bkgsig[i])
#             elif filters[i] == 'v':
#                 good_v_seeing.append(seeing[i])
#                 good_v_bkgsig.append(bkgsig[i])
#             elif filters[i] == 'r':
#                 good_r_seeing.append(seeing[i])
#                 good_r_bkgsig.append(bkgsig[i])


for i in range(len(seeing)):
    # Bad seeing
    if not math.isnan(bkgsig[i]):
        if not math.isnan(seeing[i]):
            if filters[i] == 'g':
                g_seeing.append(seeing[i])
                g_bkgsig.append(bkgsig[i])
            elif filters[i] == 'i':
                i_seeing.append(seeing[i])
                i_bkgsig.append(bkgsig[i])
            elif filters[i] == 'v':
                v_seeing.append(seeing[i])
                v_bkgsig.append(bkgsig[i])
            elif filters[i] == 'r':
                r_seeing.append(seeing[i])
                r_bkgsig.append(bkgsig[i])


folder = 'skynoise/'
if not os.path.isdir(folder):
    os.makedirs(folder)

# figa, ax1 = plt.subplots(1)
# data1 = ax1.hist(g_seeing, bins=100, range=[0,25], normed=1)
# x_fit = py.linspace(0, 25, 200)
# m, s = stats.norm.fit(g_seeing) # get mean and standard deviation
# pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
# ax1.plot(x_fit, pdf_g, label="Norm") # plot it
#
# ag,bg,cg = stats.gamma.fit(g_seeing)
# pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
# ax1.plot(x_fit, pdf_gamma, label="Gamma")
#
# params = stats.exponweib.fit(g_seeing, loc=0)
# ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
#
# # # Plot the results
# ax1.set_title('Distribution of seeing in the SkyMapper g filter')
# plt.xlabel('Seeing')
# plt.ylabel('Frequency')
# plt.legend()
# figa.savefig(folder+'g_seeing.png')
# plt.close(figa)
#
#
#
# g_bkgsig2=[a for a in g_bkgsig if (a >= 1 and a <= 100)]
#
# figa, ax1 = plt.subplots(1)
# data2 = ax1.hist(g_bkgsig2, bins=100, range=[0,100], normed=1)
# x_fit = py.linspace(0, 100, 200)
# m, s = stats.norm.fit(g_bkgsig2) # get mean and standard deviation
# pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
# ax1.plot(x_fit, pdf_g, label="Norm") # plot it
#
# ag,bg,cg = stats.gamma.fit(g_bkgsig2)
# pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
# ax1.plot(x_fit, pdf_gamma, label="Gamma")
#
# params = stats.exponweib.fit(g_bkgsig2, loc=0)
# ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
#
# ax1.set_title('Distribution of bkgsig in the SkyMapper g filter')
# plt.xlabel('bkgsig')
# plt.ylabel('Frequency')
# figa.savefig(folder+'g_bkgsig.png')
# plt.close(figa)






# figa, ax1 = plt.subplots(1)
# data1 = ax1.hist(i_seeing, bins=100, range=[0,20], normed=1)
# x_fit = py.linspace(0, 20, 200)
# m, s = stats.norm.fit(i_seeing) # get mean and standard deviation
# pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
# ax1.plot(x_fit, pdf_g, label="Norm") # plot it
#
# ag,bg,cg = stats.gamma.fit(i_seeing)
# pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
# ax1.plot(x_fit, pdf_gamma, label="Gamma")
#
# params = stats.exponweib.fit(i_seeing, loc=0)
# ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
# ax1.set_title('Distribution of seeing in the SkyMapper i filter')
# plt.xlabel('Seeing')
# plt.ylabel('Frequency')
# figa.savefig(folder+'i_seeing.png')
# plt.close(figa)




i_bkgsig2=[a for a in i_bkgsig if (a >= 1 and a <= 100)]


figa, ax1 = plt.subplots(1)
data1 = ax1.hist(i_bkgsig2, bins=100, range=[0,100], normed=1)
x_fit = py.linspace(0, 100, 200)
m, s = stats.norm.fit(i_bkgsig2) # get mean and standard deviation
pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
ax1.plot(x_fit, pdf_g, label="Norm") # plot it

ag,bg,cg = stats.gamma.fit(i_bkgsig2)
pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
ax1.plot(x_fit, pdf_gamma, label="Gamma")

params = stats.exponweib.fit(i_bkgsig2, loc=0)
ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
ax1.set_title('Distribution of bkgsig in the SkyMapper i filter')
plt.xlabel('bkgsig')
plt.ylabel('Frequency')
figa.savefig(folder+'i_bkgsig.png')
plt.close(figa)




#
#
# figa, ax1 = plt.subplots(1)
# data1 = ax1.hist(r_seeing, bins=100, range=[0,25], normed=1)
# x_fit = py.linspace(0, 25, 200)
# m, s = stats.norm.fit(r_seeing) # get mean and standard deviation
# pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
# ax1.plot(x_fit, pdf_g, label="Norm") # plot it
#
# ag,bg,cg = stats.gamma.fit(r_seeing)
# pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
# ax1.plot(x_fit, pdf_gamma, label="Gamma")
#
# params = stats.exponweib.fit(r_seeing, loc=0)
# ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
# ax1.set_title('Distribution of seeing in the SkyMapper r filter')
# plt.xlabel('Seeing')
# plt.ylabel('Frequency')
# figa.savefig(folder+'r_seeing.png')
# plt.close(figa)
#



r_bkgsig2=[a for a in r_bkgsig if (a >= 1 and a <= 100)]


figa, ax1 = plt.subplots(1)
data1 = ax1.hist(r_bkgsig2, bins=100, range=[0,100], normed=1)
x_fit = py.linspace(0, 100, 200)
m, s = stats.norm.fit(r_bkgsig2) # get mean and standard deviation
pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
ax1.plot(x_fit, pdf_g, label="Norm") # plot it

ag,bg,cg = stats.gamma.fit(r_bkgsig2)
pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
ax1.plot(x_fit, pdf_gamma, label="Gamma")

params = stats.exponweib.fit(r_bkgsig2, loc=0)
ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
ax1.set_title('Distribution of bkgsig in the SkyMapper r filter')
plt.xlabel('bkgsig')
plt.ylabel('Frequency')
figa.savefig(folder+'r_bkgsig.png')
plt.close(figa)





#
# figa, ax1 = plt.subplots(1)
# data1 = ax1.hist(v_seeing, bins=100, range=[0,25], normed=1)
# x_fit = py.linspace(0, 20, 200)
# m, s = stats.norm.fit(v_seeing) # get mean and standard deviation
# pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
# ax1.plot(x_fit, pdf_g, label="Norm") # plot it
#
# ag,bg,cg = stats.gamma.fit(v_seeing)
# pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
# ax1.plot(x_fit, pdf_gamma, label="Gamma")
#
# params = stats.exponweib.fit(v_seeing, loc=0)
# ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
# ax1.set_title('Distribution of seeing in the SkyMapper v filter')
# plt.xlabel('Seeing')
# plt.ylabel('Frequency')
# figa.savefig(folder+'v_seeing.png')
# plt.close(figa)
#



v_bkgsig2=[a for a in v_bkgsig if (a >= 1 and a <= 50)]


figa, ax1 = plt.subplots(1)
data1 = ax1.hist(v_bkgsig2, bins=100, range=[0,50], normed=1)
x_fit = py.linspace(0, 50, 200)
m, s = stats.norm.fit(v_bkgsig2) # get mean and standard deviation
pdf_g = stats.norm.pdf(x_fit, m, s) # now get theoretical values in our interval
ax1.plot(x_fit, pdf_g, label="Norm") # plot it

ag,bg,cg = stats.gamma.fit(v_bkgsig2)
pdf_gamma = stats.gamma.pdf(x_fit, ag, bg,cg)
ax1.plot(x_fit, pdf_gamma, label="Gamma")

params = stats.exponweib.fit(v_bkgsig, loc=0)
ax1.plot(x_fit,stats.exponweib.pdf(x_fit,*params),label='Weibull')
ax1.set_title('Distribution of bkgsig in the SkyMapper v filter')
plt.xlabel('bkgsig')
plt.ylabel('Frequency')
figa.savefig(folder+'v_bkgsig.png')
plt.close(figa)