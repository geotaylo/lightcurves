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
dir1 = 'Simulations/run2/cad4/30SN/no_k_bad_seeing/'
dir2 = 'Simulations/run2/cad4/30SN/no_k_good_seeing/'
dir3 = 'Simulations/run2/cad4/30SN/k_bad_seeing/'
dir4 = 'Simulations/run2/cad4/30SN/k_good_seeing/'

dir1_5 = 'Simulations/run2/cad4/50SN/no_k_bad_seeing/'
dir2_5 = 'Simulations/run2/cad4/50SN/no_k_good_seeing/'
dir3_5 = 'Simulations/run2/cad4/50SN/k_bad_seeing/'
dir4_5 = 'Simulations/run2/cad4/50SN/k_good_seeing/'

dir1_1 = 'Simulations/run2/cad4/100SN/no_k_bad_seeing/'
dir2_1 = 'Simulations/run2/cad4/100SN/no_k_good_seeing/'
dir3_1 = 'Simulations/run2/cad4/100SN/k_bad_seeing/'
dir4_1 = 'Simulations/run2/cad4/100SN/k_good_seeing/'

folders = [dir1, dir2, dir3, dir4]
folders_5 = [dir1_5, dir2_5, dir3_5, dir4_5]
folders_1 = [dir1_1, dir2_1, dir3_1, dir4_1]


# Number of SN to simulate (if using sncosmo distribution, use 0).
nSNe = 100


# Generates random SN and simulates observed data.

run1 = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir1, kpass=False, follow_up=False)
run2 = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir2, kpass=False, follow_up=True)
run3 = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir3, kpass=True, follow_up=False,
                       properties=dir1+'sn_dict')
run4 = run.simulate_lc(nSNe=nSNe, cadence=4, folder=dir4, kpass=True, follow_up=True,
                       properties=dir2+'sn_dict')


# Import list of lightcurves from files (if not simulating above).

#run3 = run.get_lc(['best_zps/noK/observed_lc_1.txt',
#                     'best_zps/noK/observed_lc_2.txt'
#                  ])

                        
# Fit SALT2 models to list of observed light curves.

run.fit_snlc(run1, folder=dir1, properties=dir1+'sn_dict')
run.fit_snlc(run2, folder=dir2, properties=dir2+'sn_dict')
run.fit_snlc(run3, folder=dir3, properties=dir3+'sn_dict')
run.fit_snlc(run4, folder=dir4, properties=dir4+'sn_dict')
"""

# Find residuals (or at least something)

x = run.get_diff(folders, 'Simulations/run1/cad4/30SN/')
y = run.get_diff(folders_5, 'Simulations/run1/cad4/50SN/')
z = run.get_diff(folders_1, 'Simulations/run1/cad4/100SN/')
"""
