#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
GT 28/11/16
WIP version of simulation_new.py translation into sncosmo.

Aims to:
    a) generate SALT2 observation light curves for SN with randomly generated parameters and observations.
    b) take an observation light curve and fit a SALT2 model to it.
    c) provide cosmoloy

Requires SNCosmo to be installed: https://sncosmo.readthedocs.io/en/v1.4.x/index.html
"""
import sm_filters_v1_1    # Registers bandpasses
import sncosmo
import numpy as np
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages

#  "CONSTANTS"

q0 = 0.2
c_light = 3*10**5     #  Speed of light, km/s
H0 = 70.              #  Hubble constant, km/s/Mpc
area = 5.7            #  Skymapper field of view
filter_col = ['indigo', 'forestgreen', 'r', 'gold']     #  Colours for plotting bandpasses
model = sncosmo.Model(source='salt2')   #  Default SALT2 model


### Make this an input parameter
smbands = ['smv', 'smg', 'smr', 'smi']  #  SkyMapper bandpasses


#  "METHODS"

#  Distance modulus formula for obtaining X0
def mu(z):
        d_L = c_light*(z + 0.5*(1-q0)*z**2)/H0
        return 5*np.log10(d_L) + 25

#  Run 3 step program.

#  Generate light curve from random SNe and observations
def gen_lc(nSNe=0, tmin=57754, tmax=58118, zmin=0.001, zmax=0.07):
    
    #  Generate random parameters for SN
    #  Redshift:
    if nSNe == 0:
        z = list(sncosmo.zdist(zmin, zmax, time=(tmax-tmin), area=area))
        nSNe = len(z)
        print 'Obtained %s'%nSNe + ' redshifts from sncosmo distribution'
    else:
        z = np.random.uniform(zmin, zmax, size=nSNe)
        print 'Obtained %s'%nSNe + ' redshifts from uniform distribution'
       
    t0 = np.random.uniform(tmin, tmax)   #  Time at lightcurve peak
    c = 0.3*np.random.randn(nSNe)        #  SALT2 PARAMETER
    x1 = 3*np.random.randn(nSNe)         #  SALT2 PARAMETER
    
    #  Assign params to model for each SN
    params = []
    for i in range(nSNe):
        x0 = 10**((29.69-mu(z[i]))/2.5)             #  SALT2 PARAMETER
        p = {'z':z[i], 't0':t0, 'x0':x0, 'x1': x1[i], 'c': c[i]}
        params.append(p)
        
        
    #  Set observations
    nObs = 10        #  Number of observations for each filter, for each SN
    nPoints = nObs*len(smbands)     #  Total number of observations for each SN
    zp = 25.        #  Zero point for scaling flux, currently default
    tdet = np.random.randint(-15, -2)    #  Time of initial detection, relative to t0
    tobs = np.random.randint(25, 65)     #  Total observing period (days) 
    timearr = np.linspace(t0+tdet, t0+tdet+tobs,nPoints)
    time = timearr.tolist()
        
    obs = Table({'band':np.array(nObs * smbands),
             'time': time,
             'zp':nPoints*[zp],
             'zpsys':nPoints*['ab'],
             'gain':nPoints*[1.],
             'skynoise':nPoints*[160.40]})
        
        
    #  Generate light curve for each SN observation
    lcs = sncosmo.realize_lcs(obs, model, params)
    for k in range(nSNe):
        k_true=k+1
        name = 'TestFiles/observed_lc_%s.txt'%k_true
        sncosmo.write_lc(lcs[k], name)
        print 'Simulated observations for supernova %s saved in %s'%(k_true,name) 
    return lcs
    
#  Fit light curve to observations and plot
def fit_util_lc(lightcurve, index):
    print 'Fitting light curve to supernova %s'%index
    observed_lc = lightcurve
    plotname = 'TestFiles/fitted_lc_%s.pdf'%index
    result, fitted_model = sncosmo.mcmc_lc(observed_lc, model,
                                          ['z', 't0', 'x0', 'x1', 'c'],           # parameters of model to vary
                                          bounds={'z':(0.001, 0.1)}, minsnr=1.0)  # bounds on parameters (if any)
    #plot
    ### Colours not done
    pp = PdfPages(plotname)
    sncosmo.plot_lc(observed_lc, model=fitted_model, errors=result.errors, fname = pp, format ='pdf')
    pp.close()
    print 'Fitted light curve plotted in '+plotname
    return
    
#  Function for fitting multiple sn
def fit_sn_lc(lightcurve):
    for i in range(len(lightcurve)):
        fit_util_lc(lightcurve[i], i+1)
    return

    
#  Function for reading lightcurves from files
#  Note filenames should be in a list
def get_lc(filename):
    lc = []
    for i in range(len(filename)):
        print 'Reading light curve of supernova %s'%(i+1)
        x = sncosmo.read_lc(filename[i])
        lc.append(x)
    return lc
