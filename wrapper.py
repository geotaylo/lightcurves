"""
GT 05/12/2016

Wrapper script to run sncosmo light curve simulations from run_sncosmo.py.  
Can be edited to suit user.
Cadences etc. must be changed inside run_sncosmo.py script, but number of 
supernovae and which telescope to simulate data from can happen in here.

(1) Generates set of supernovae parameters for a range of observations;
(2) Simulates observed light curves of supernovae for the SkyMapper or Kepler 
    telescopes (or both), with parameters sampled from (1)...
    ... OR imports observed light curve;
(3) Fits SALT2 model to observations.

Contact Georgie: g.l.taylor@outlook.com with any questions!
"""

import run_sncosmo as run


# Number of supernovae to simulate (if using sncosmo distribution, use 0).
nSNe = 30


# Path to store info about whole observing set, used by each run.
# -- ENSURE / is at end of path! --
parent_folder = 'Simulations/cad_sm_4d/cad_k_6h/30sn/'


# Paths to store individual runs
# -- ENSURE / is at end of path! --
child_folder_1 = 'sm_bad_seeing/'
child_folder_2 = 'sm_good_seeing/'
child_folder_3 = 'k_bad_seeing/'
child_folder_4 = 'both_bad_seeing/'
child_folder_5 = 'both_good_seeing/'


# ---------- Generates all info about SN and observing parameters, ----------
#             to be used by each light curve simulation and fit.

# -- MUST BE RUN BEFORE simulate_lc! --
run.simulate_sn_set(parent_folder, nSNe)


# ---------- Simulates observations of supernovae from above set. ----------

print 'Attempting to simulate observation 1: SkyMapper in bad seeing.'
run1 = run.simulate_lc(parent_folder, child_folder=child_folder_1, scope='sm', 
                       follow_up=False)


print 'Attempting to simulate observation 2: SkyMapper in good seeing.'
run2 = run.simulate_lc(parent_folder, child_folder=child_folder_2, scope='sm', 
                       follow_up=True)


print 'Attempting to simulate observation 3: Kepler.'
run3 = run.simulate_lc(parent_folder, child_folder=child_folder_3, scope='kst',
                       follow_up=False)


print 'Attempting to simulate observation 4: SkyMapper and Kepler in bad \
       seeing.'
run4 = run.simulate_lc(parent_folder, child_folder=child_folder_4, 
                       scope='both', follow_up=False)


print 'Attempting to simulate observation 5: SkyMapper and Kepler in good \
       seeing.'
run5 = run.simulate_lc(parent_folder, child_folder=child_folder_5,
                       scope='both', follow_up=True)


# ---------- Imports photometric data (if not simulated above) ----------

# Uncomment if needed!
# Note - parent folder and simulate_sn_set() is still needed for fit_snlc()

#run6 = run.get_lc(['best_zps/noK/observed_lc_1.txt', 
#                   'best_zps/noK/observed_lc_2.txt'
#                  ])


# ---------- Fit SALT2 models to light curves ----------

print 'Attempting fit 1: SkyMapper in bad seeing.'
run.fit_snlc(run1, parent_folder, child_folder=child_folder_1)


print 'Attempting fit 2: SkyMapper in good seeing.'
run.fit_snlc(run2, parent_folder, child_folder=child_folder_2)


print 'Attempting fit 3: Kepler.'
run.fit_snlc(run3, parent_folder, child_folder=child_folder_3)


print 'Attempting fit 4: SkyMapper and Kepler in bad seeing.'
run.fit_snlc(run4, parent_folder, child_folder=child_folder_4)


print 'Attempting fit 5: SkyMapper and Kepler in good seeing.'
run.fit_snlc(run5, parent_folder, child_folder=child_folder_5)
