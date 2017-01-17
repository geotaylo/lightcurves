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


# Number of SN to simulate (if using sncosmo distribution, use 0).
nSNe = 1


# Path to store info about whole observing set, used by each run.
# ENSURE / is at end of path!
parent_folder = 'Simulations/testing/flipping_methods/1/'


# Paths to store individual runs
# ENSURE / is at end of path!
child_folder_1 = 'sm_bad_seeing/'
child_folder_2 = 'sm_good_seeing/'
child_folder_3 = 'k_bad_seeing/'
child_folder_4 = 'k_good_seeing/'
child_folder_5 = 'both_bad_seeing/'
child_folder_6 = 'both_good_seeing'


# Generates all info about SN and observing parameters, to be used by each
# light curve simulation and fit.
# MUST BE RUN BEFORE simulate_lc!
run.simulate_sn_set(parent_folder, nSNe)


# Generates random SN and simulates observed data.
print'attempting run 1:'
run1 = run.simulate_lc(parent_folder, child_folder=child_folder_1,
                       scope='sm', follow_up=False)

print'attempting run 2:'
run2 = run.simulate_lc(parent_folder, child_folder=child_folder_2,
                       scope='sm', follow_up=True)

print'attempting run 3:'
run3 = run.simulate_lc(parent_folder, child_folder=child_folder_3,
                       scope='kst', follow_up=False)

print'attempting run 4:'
run4 = run.simulate_lc(parent_folder, child_folder=child_folder_4,
                       scope='kst', follow_up=True)

print'attempting run 5:'
run5 = run.simulate_lc(parent_folder, child_folder=child_folder_5,
                       scope='both', follow_up=False)

print'attempting run 6:'
run6 = run.simulate_lc(parent_folder, child_folder=child_folder_6,
                       scope='both', follow_up=True)

# Import list of lightcurves from files (if not simulating above).

#run3 = run.get_lc(['best_zps/noK/observed_lc_1.txt',
#                     'best_zps/noK/observed_lc_2.txt'
#                  ])

                        
# Fit SALT2 models to list of observed light curves.

print'attempting fit 1:'
#run.fit_snlc(run1, folder=dir1, properties=dir1+'sn_dict')

