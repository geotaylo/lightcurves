"""
GT 03/04/2016
Test wrapper script to run sncosmo light curve simulations from
run_sncosmo.py.  Can be edited to suit user.

> Generates observed light curve from SN parameters (random or
  specified);
> OR imports observed light curve;
> THEN fits SALT2 model to observations.
"""

import heaven_and_earth as run

# Number of SN to simulate (if using sncosmo distribution, use 0).
nSNe = 1

# Kepler campaign (if using SNCosmo dist.)
campaign = 0

# Path to store info about whole observing set, used by each run.
# ENSURE / is at end of path!
parent_folder = 'Honours_data_sets/test/'

# Paths to store individual runs
# ENSURE / is at end of path!
child_folder_1 = 'sm_ws/'
child_folder_2 = 'sm_ps/'
child_folder_3 = 'kst/'
child_folder_4 = 'combined_ws/'
child_folder_5 = 'combined_ps/'


# Generates all info about SN and observing parameters, to be used by each
# light curve simulation and fit.
# MUST BE RUN BEFORE simulate_lc!
nSNe = run.simulate_sn_set(parent_folder, nSNe, campaign)

# Generates random SN and simulates observed data.
print 'attempting well-sampled SM simulations:'

run1 = run.simulate_lc(parent_folder, child_folder=child_folder_1,
                      scope='sm', sm_cad='well-sampled', follow_up=False)

print 'attempting poorly-sampled SM simulations'
run2 = run.simulate_lc(parent_folder, child_folder=child_folder_2,
                      scope='sm', sm_cad='poorly-sampled', follow_up=False)

print'attempting KST simulations:'
run3 = run.simulate_lc(parent_folder, child_folder=child_folder_3,
                      scope='kst', follow_up=False)

print 'attempting well-sampled SM + KST combined simulations'
run4 = run.combine_scopes(parent_folder, child_folder_1, child_folder_3,
                         child_folder_4, nSNe)

print 'attempting poorly-sampled SM + KST combined simulations'
run5 = run.combine_scopes(parent_folder, child_folder_2, child_folder_3,
                         child_folder_5, nSNe)

# Write observational parameters
run.write_params(parent_folder, nSNe)



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
#
# run3 = run.get_lc([parent_folder+child_folder_3+'observed_lc_1.txt',
#                    parent_folder+child_folder_3+'observed_lc_2.txt',
#                    parent_folder+child_folder_3+'observed_lc_3.txt',
#                    parent_folder+child_folder_3+'observed_lc_4.txt',
#                    parent_folder+child_folder_3+'observed_lc_5.txt',
#                  ])


# Fit SALT2 models to list of observed light curves.
print'attempting well-sampled SM fits:'
run.fit_snlc(run1, parent_folder, child_folder=child_folder_1)

print'attempting poorly-sampled SM fits:'
run.fit_snlc(run2, parent_folder, child_folder=child_folder_2)

print'attempting KST fits:'
kepler_time = run.fit_snlc(run3, parent_folder, child_folder=child_folder_3)

print'attempting well-sampled combined fits:'
run.fit_snlc(run4, parent_folder, child_folder=child_folder_4, t0_in=kepler_time)

print'attempting poorly-sampled combined fits:'
run.fit_snlc(run5, parent_folder, child_folder=child_folder_5, t0_in=kepler_time)