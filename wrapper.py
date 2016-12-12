"""
GT 05/12/2016

Wrapper script to run sncosmo light curve simulations from 
run_sncosmo.py.  Can be edited to suit user.

> Generates observed light curve from SN parameters (random or 
  specified);
> OR imports observed light curve;
> THEN fits SALT2 model to observations.
"""

import run_sncosmo as run

# Edit these as needed

# Paths to store simulations and results (make sure / is at end)
dir1 = 'TestFiles/'  # Without kepler filter (1st rn)
dir2 = 'randomDir2/'  # With kepler filter (2nd run)

# Number of SN to simulate (if using sncosmo distribution, use 0).
nSNe = 4    


# Generates random SN and simulates observed data.

# Without kepler filter (1st run).
"""run1 = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir1)

# With kepler filter (2nd run).  To use the same set of SN as 1st run,
# make sure to use dir1 as properties path.
run2 = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir2, kpass=True, 
                       properties=dir1+'sn_dict')"""


# Import list of lightcurves from files (if not simulating above).

run3 = run.get_lc(['TestFiles/observed_lc_1.txt', 
                     'TestFiles/observed_lc_2.txt'
                     ])

                        
# Fit SALT2 models to list of observed light curves.

"""run.fit_snlc(run1, folder=dir1, properties=dir1+'sn_dict')

run.fit_snlc(run2, folder=dir2, properties=dir2+'sn_dict')"""

run.fit_snlc(run3, folder='Testfiles2/')
