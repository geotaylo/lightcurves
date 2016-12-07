"""
GT 05/12/16
Aims to:
a) generate SALT2 observation light curves for SN with randomly 
   generated parameters and observations.
b) take an observation light curve and fit a SALT2 model to it.
c) provide cosmoloy

Requires SNCosmo to be installed: 
https://sncosmo.readthedocs.io/en/v1.4.x/index.html
"""

import os
import numpy as np
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages

import sfdmap
import sncosmo
from sncosmo import mcmc_lc as fit

import filters  # Registers bandpasses.


AREA = 5.7  # SkyMapper field of view.

smbands = ['smv', 'smg', 'smr', 'smi']  # Bandpasses to use.

allbands = ['smv', 'smg', 'smr', 'smi','kst']

dust = sncosmo.CCM89Dust()

dustmap = sfdmap.SFDMap("/home/georgie/sfddata-master")

model = sncosmo.Model(source='salt2',
                      effects=[dust, dust],
                      effect_names=['host', 'mw'],
                      effect_frames=['rest', 'obs'])


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        

def mu(z):
    
    """ Distance modulus formula used to obtain x0 """
    
    Q0 = 0.2
    LIGHT = 3*10**5  # km/s
    H0 = 70.  # km/s/Mpc
            
    d_L = LIGHT*(z + 0.5*(1 - Q0)*z **2)/H0

    return 5*np.log10(d_L) + 25


def simulate_lc(nSNe=0, cadence=4, kpass=False, folder='TestFiles/',
                tmin=57754, tmax=58118, zmin=0.001, zmax=0.1
                ):
    
    """ Generate SN parameters and simulate 'observed' light curve 
        Defaults:  nSNe = 0 (number of SN)
                   cadence = 4 (days between observations)
                   kpass = false (include kepler bandpass)
                   folder = 'Testfiles/' (path to store outputs)
                   tmin = 57754 (01/01/2017 in mjd)
                   tmax = 58118 (31/12/2017 in mjd)
                   zmin = 0.001 (min redshift detectable by SkyMapper)
                   zmax = 0.100 (max redshift detectable by SkyMapper)
    """
    
    ensure_dir(folder)
               
    if kpass:
        bands = allbands
    else:
        bands = smbands 
                           
    if nSNe == 0:
        
        # Obtain redshifts.
        z = list(sncosmo.zdist(zmin, zmax, 
                              time=(tmax - tmin), area=AREA)
                              )
        nSNe = len(z)
        
        print 'Obtained %s'%nSNe \
              + ' redshifts from sncosmo distribution'
              
    else:
                  
        # Obtain redshifts.
        z = np.random.uniform(zmin, zmax, size=nSNe)
                  
        print 'Obtained %s'%nSNe \
              + ' redshifts from uniform distribution'
                  
    t0 = np.random.uniform(tmin, tmax)  # Time at lightcurve peak (mjd).
    
    c = 0.3*np.random.randn(nSNe)
    
    x1 = 3*np.random.randn(nSNe)

    params = []

    for i in range(nSNe):
        x0 = 10**((29.69 - mu(z[i]))/2.5) 
        p = {'z': z[i], 't0': t0, 'x0': x0, 'x1': x1[i], 'c': c[i]}
        params.append(p)        
        
    # Time of initial detection (mjd).
    tdet = np.random.randint(-15, -2) 
    
    # Total observing period (days).
    tobs = np.random.randint(25, 65)  
    
    timearr = t0 + np.arange(tdet, tdet + tobs, cadence)
    
    time = timearr.tolist()
    
    nObs = len(time)
    
    nPoints = nObs*len(bands)
    
    # Zero point for scaling.
    ZP = 25.  
    
        
    obs = Table({'band':np.array(nObs*bands),
                 'time': sorted(time*len(bands)),
                 'zp':nPoints*[ZP],
                 'zpsys':nPoints*['ab'],
                 'gain':nPoints*[1.],
                 'skynoise':nPoints*[160.40]
                 })
        
    # Generate light curve for each SN observation.
    lcs = sncosmo.realize_lcs(obs, model, params)
    
    # Create text file for storing 'true' parameters for each SN.
    true_file = folder + 'true_parameters.txt'
    tf = open(true_file, 'w')
    
    for k in range(nSNe):
        name = folder + 'observed_lc_%s.txt'%(k + 1)
        
        # Write observations.
        sncosmo.write_lc(lcs[k], name)
        
        print 'Simulated observations for supernova %s saved in %s'\
               %((k + 1), name)
               
        # Write true parameters in text file.
        tf.write('SN%s: '%(k+1))
        
        true_params = lcs[k].meta

        for key in sorted(true_params):
            tf.write('%s:%s ' % (key, true_params[key]))
            
        tf.write('\n')
        
        print 'True parameters for supernova %s saved in %s'\
               %((k + 1), true_file)
              
    tf.close 
    return lcs
               
               
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
        
        print 'Fitted parameters for supernova %s saved in %s'\
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
        print 'Reading light curve of supernova %s'%(i + 1)
        x = sncosmo.read_lc(filelist[i])
        lc.append(x)
        
    return lc
                