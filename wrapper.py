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


# Generates random SN and simulates observed data.
test = run.simulate_lc(nSNe=1, cadence=3, kpass=True)


# Import list of lightcurves from files (if not simulating above).
#test_1 = run.get_lc(['TestFiles/observed_lc_1.txt', 
#                     'TestFiles/observed_lc_2.txt'
#                     ])


# Fit SALT2 models to list of observed light curves.
res = run.fit_snlc(test)

