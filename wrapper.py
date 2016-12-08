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

dir = '100sim/noK/cad4/'
dir2 = '100sim/K/cad4/'

# Generates random SN and simulates observed data.
no_k_lc, params, obs = run.simulate_lc(nSNe=100, cadence=4, folder=dir)
k_lc, params2, obs2 = run.simulate_lc(nSNe=100, cadence=4, folder=dir2, 
                                kpass=True, params=params, obs_in=obs)


# Import list of lightcurves from files (if not simulating above).
#test_1 = run.get_lc(['TestFiles/observed_lc_1.txt', 
#                     'TestFiles/observed_lc_2.txt'
#                     ])


# Fit SALT2 models to list of observed light curves.
res1 = run.fit_snlc(no_k_lc, folder=dir)
res2 = run.fit_snlc(k_lc, folder=dir2)

# Calculate differences between true and fitted paramaters
#run.get_diff(dir)