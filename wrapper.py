"""
GT 03/04/2016
Test wrapper script to run sncosmo light curve simulations from
run_sncosmo.py.  Can be edited to suit user.

> Generates observed light curve from SN parameters (random or
  specified);
> OR imports observed light curve;
> THEN fits SALT2 model to observations.
"""

import run_sncosmo as run

# Number of SN to simulate (if using sncosmo distribution, use 0).
nSNe = 100

# Path to store info about whole observing set, used by each run.
# ENSURE / is at end of path!
parent_folder = 'Honours_data_sets/4_080218/1_stat_sample \
                '/Kepler_6hours/SM_5day/vObs_2/100_SN/'


# Paths to store individual runs
# ENSURE / is at end of path!
child_folder_1 = 'sm_bad_seeing/'
child_folder_2 = 'sm_good_seeing/'
child_folder_3 = 'kst/'
child_folder_4 = 'both_bad_seeing/'
child_folder_5 = 'both_good_seeing/'


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

print 'attempting bad seeing combination'
run4 = run.combine_scopes(parent_folder, child_folder_1, child_folder_3,
                         child_folder_4, nSNe)

print 'attempting good seeing combination'
run5 = run.combine_scopes(parent_folder, child_folder_2, child_folder_3,
                         child_folder_5, nSNe)


# Import list of lightcurves from files (if not simulating above).
# Note - parent folder is still needed for fit_snlc()

# run1 = run.get_lc([parent_folder+child_folder_1+'observed_lc_1.txt',
#                    parent_folder+child_folder_1+'observed_lc_2.txt',
#                    parent_folder+child_folder_1+'observed_lc_3.txt',
#                    parent_folder+child_folder_1+'observed_lc_4.txt',
#                    parent_folder+child_folder_1+'observed_lc_5.txt',
#                  ])

# run2 = run.get_lc([parent_folder+child_folder_2+'observed_lc_1.txt',
#                     parent_folder+child_folder_2+'observed_lc_2.txt',
#                     parent_folder+child_folder_2+'observed_lc_3.txt',
#                     parent_folder+child_folder_2+'observed_lc_4.txt',
#                     parent_folder+child_folder_2+'observed_lc_5.txt',
#                   ])

# run3 = run.get_lc([parent_folder+child_folder_3+'observed_lc_1.txt',
#                    parent_folder+child_folder_3+'observed_lc_2.txt',
#                    parent_folder+child_folder_3+'observed_lc_3.txt',
#                    parent_folder+child_folder_3+'observed_lc_4.txt',
#                    parent_folder+child_folder_3+'observed_lc_5.txt',
#                  ])

# run4 = run.get_lc([parent_folder+child_folder_4+'observed_lc_1.txt',
#                    parent_folder+child_folder_4+'observed_lc_2.txt',
#                    parent_folder+child_folder_4+'observed_lc_3.txt',
#                    parent_folder+child_folder_4+'observed_lc_4.txt',
#                    parent_folder+child_folder_4+'observed_lc_5.txt',
#                  ])

# run5 = run.get_lc([parent_folder+child_folder_5+'observed_lc_1.txt',
#                     parent_folder+child_folder_5+'observed_lc_2.txt',
#                     parent_folder+child_folder_5+'observed_lc_3.txt',
#                     parent_folder+child_folder_5+'observed_lc_4.txt',
#                     parent_folder+child_folder_5+'observed_lc_5.txt',
#                   ])

# Fit SALT2 models to list of observed light curves.
print'attempting fit 1:'
run.fit_snlc(run1, parent_folder, child_folder=child_folder_1)

print'attempting fit 2:'
run.fit_snlc(run2, parent_folder, child_folder=child_folder_2)

print'attempting fit 3:'
kepler_time = run.fit_snlc(run3, parent_folder, child_folder=child_folder_3)

####### 4 and 5 (combined scopes) need kepler time!

print'attempting fit 4:'
run.fit_snlc(run4, parent_folder, child_folder=child_folder_4, t0=kepler_time)

print'attempting fit 5:'
run.fit_snlc(run5, parent_folder, child_folder=child_folder_5, t0=kepler_time)
