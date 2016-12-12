"""
GT 05/12/2016

Wrapper script to run sncosmo light curve simulations from sim_new.pp.
This is an example; simulations can also be run directly from command
line.

Generate observed light curve from SN parameters (random or specified);
OR import observed light curve;
THEN fit SALT2 model to observations.
"""

import run_sncosmo as run

dir1 = 'randomDir/'
dir2 = 'randomDir2/'
nSNe = 1

# Generates random SN and simulates observed data.
no_k_lc = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir1)
k_lc = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir2, kpass=True, 
                       properties=dir1+'sn_dict')



# Import list of lightcurves from files (if not simulating above).
#test_1 = run.get_lc(['TestFiles/observed_lc_1.txt', 
#                     'TestFiles/observed_lc_2.txt'
#                     ])
# Get random coordinates for imported lc SN (if not provided)
#coords3 = run.get_coords(nSNe=2)

# Fit SALT2 models to list of observed light curves.
run.fit_snlc(no_k_lc, folder=dir1, properties=dir1+'sn_dict')
run.fit_snlc(k_lc, folder=dir2, properties=dir2+'sn_dict')

# Calculate differences between true and fitted paramaters
#run.get_diff(dir)