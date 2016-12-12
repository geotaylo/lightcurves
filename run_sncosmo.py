"""
GT 05/12/16
Aims to:
a) generate SALT2 observation light curves for SN with randomly 
   generated parameters and observations.
b) take an observation light curve and fit a SALT2 model to it.
c) provide cosmoloy

Requires installation of:
- SNCosmo:  https://sncosmo.readthedocs.io/en/v1.4.x/index.html
- Emcee:
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages

import sfdmap
import sncosmo
from sncosmo import mcmc_lc as fit

import filters


# CONSTANTS _-------------------------------------------------------


AREA = 5.7  # SkyMapper field of view.

Q0 = 0.2

LIGHT = 3*10**5  # km/s

H0 = 70.  # km/s/Mpc

smbands = ['smv', 'smg', 'smr', 'smi']  # SkyMapper bandpasses to use.

allbands = ['smv', 'smg', 'smr', 'smi','kst']  # SkyMapper + KST bands

dust = sncosmo.CCM89Dust()

# Change path to location of dustmaps
dustmap = sfdmap.SFDMap("/home/georgie/sfddata-master")

model = sncosmo.Model(source='salt2',
                      effects=[dust, dust],
                      effect_names=['host', 'mw'],
                      effect_frames=['rest', 'obs'])

# USEFUL FNS ----------------------------------------------------


def simulate_lc(nSNe=0, cadence=4, kpass=False, folder='TestFiles/',
                tmin=57754, tmax=58118, zmin=0.001, zmax=0.1, 
                params=[], obs_in=[]
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
                   params = [] (SN parameters, generates if blank)
    """
        
    # Maintenance stuff
    filters.register_filters()
    ensure_dir(folder)
    
    # Select filters.           
    if kpass:
        bands = allbands
    else:
        bands = smbands 
        
    lcs=[]
    observing_list = []
    # GENERATING SUPERNOVAE PARAMETERS (IF NOT PROVIDED)
    if (params == []):
        # Obtain redshifts                    
        if nSNe == 0:
            z = list(sncosmo.zdist(zmin, zmax, 
                                  time=(tmax - tmin), area=AREA)
                                  )
            nSNe = len(z)
            print 'Obtained %s'%nSNe \
                  + ' redshifts from sncosmo distribution. \n'     
        else:
            z = np.random.uniform(zmin, zmax, size=nSNe)
            print 'Obtained %s'%nSNe \
                  + ' redshifts from uniform distribution. \n'
                      
        t0 = np.random.uniform(tmin, tmax, nSNe)  # Time at lightcurve peak (mjd).
        c = 0.3*np.random.randn(nSNe)   # Colour
        x1 = 3*np.random.randn(nSNe)    # Stretch
    
        for i in range(nSNe):
            x0 = 10**((29.69 - mu(z[i]))/2.5)   # Amplitude (?)
            p = {'z': z[i], 't0': t0[i], 'x0': x0, 'x1': x1[i], 'c': c[i]}
            params.append(p)        
    
    true_file = folder + 'true_parameters.txt'
    tf = open(true_file, 'w')
    
    for t in range(nSNe):
        # CREATING TABLE OF OBSERVATIONS            
        if (obs_in == []):
            tdet = np.random.randint(-15, -2)  # Time of init detection (mjd).
            tobs = np.random.randint(25, 65)  # Total observing period (days).
            time = (params[t].get('t0')
                    + np.arange(tdet, tdet + tobs, cadence)).tolist()
            nObs = len(time)
            ZP = 25. # Zero point for scaling (filler val)
            gain = 1. # Filler value
            skynoise = 100. # Filler value
        else:
            time = obs_in[t][0]
            nObs = obs_in[t][1]
            ZP = obs_in[t][2]
            gain = obs_in[t][3]
            skynoise = obs_in[t][4]
            
        o_t = []
        # Adds 20s offset time between observations in different filters
        for x in range(len(bands)):
            j = np.array(time) + 0.00023148148*x
            o_t.append(j.tolist())
        observing_time = [item for sublist in o_t for item in sublist]
            
        nPoints = nObs*len(bands)
        obs_util = [time, nObs, ZP, gain, skynoise] 
        observing_dictionary = {'band':np.array(nObs*bands),
                     'time': sorted(observing_time),
                     'zp':nPoints*[ZP],
                     'zpsys':nPoints*['ab'],
                     'gain':nPoints*[gain],
                     'skynoise':nPoints*[skynoise]
                     }  
        observing_list.append(obs_util) 
        obs = Table(observing_dictionary)
        lc = sncosmo.realize_lcs(obs, model, [params[t]])
        lcs.append(lc[0])
        
        name = folder + 'observed_lc_%s.txt'%(t + 1)
        
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
    return lcs, params, observing_list
                   
                   
def fit_snlc(lightcurve, folder='TestFiles/'):
                       
    """ Utility function for fitting lightcurves to 
    observations of multiple SN """
                 
    # Create text file for storing 'fitted' parameters for each SN.
    fitted_file = folder + 'fitted_parameters.txt'
    ff = open(fitted_file, 'w')
      
    for i in range(len(lightcurve)):
        p = fit_util_lc(lightcurve[i], i + 1, folder)
        
        # Write fitted parameters in text file.
        ff.write('SN%s: '%(i+1))
        
        for key in sorted(p):
            ff.write('%s:%s ' % (key, p[key]))
            
        ff.write('\n')
        
        print 'Fitted parameters for supernova %s saved in %s \n'\
               %((i + 1), fitted_file)

    ff.close()
    
    print 'Process complete.'
                       
    return
                           
                           
def get_lc(filelist):
                               
    """ Reads lightcurves from files
        File names should be in a list
        Returns list of light curve tables """
                               
    lc = []

    for i in range(len(filelist)):
        print 'Reading light curve of supernova %s \n'%(i + 1)
        x = sncosmo.read_lc(filelist[i])
        lc.append(x)
        
    return lc
    
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

    with open(fp, 'rb') as file:        
        
        reader = csv.reader(file,delimiter=' ')
        
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
        
        reader = csv.reader(file,delimiter=' ')
        
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
    

# UTILITY FNS --------------------------------------------------------


def ensure_dir(f):
    
    """ Checks if specified path exists, and creates it if not. """
    
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        

def mu(z):
    
    """ Distance modulus formula used to obtain x0. """
            
    d_L = LIGHT*(z + 0.5*(1 - Q0)*z **2)/H0
    return 5*np.log10(d_L) + 25

def fit_util_lc(data, index, folder):
                   
    """ Fits SALT2 light curve to a supernova observation
        Plots model, data and residuals """
                   
    print 'Fitting light curve for supernova %s'%index
                   
    plotname = folder + 'fitted_lc_%s.pdf'%index
    
    # Generating random coords for dustmap.
    
    # Galactic longitude.
    l = np.random.uniform(0, 360)
    
    # Galactic latitude.
    b = np.random.uniform(-90, -30)
    
    ebv = dustmap.ebv(l, b, frame='galactic', unit='degree')
    
    model.set(mwebv=ebv)
    
                   
    # Fitting SALT2 model using MCMC. 
    result, fitted_model = fit(data, model,
                               # Parameters of model to vary.
                               ['z', 't0', 'x0', 'x1', 'c'], 
                               bounds={'z':(0.001, 0.1)}, minsnr=1.0
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
                