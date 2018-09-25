"""
GT 21/09/2018
Fitting real data from SNCosmo rebuild.
"""

import run_sncosmo_k2fields as run

# Path to store info about whole observing set, used by each run.
# ENSURE / is at end of path!
parent_folder = 'Real_data/KEGS/2017k/'

# Names of SN
name = '2017k'

# List of redshifts
z = [0.13182]

# List of RA (in J2000 frame and degrees)
ra = [run.dms_to_deg(1.,58.,6.)]

# List of dec (in J2000 frame and degrees)
dec = [run.dms_to_deg(13.,46.,44.2)]

# Make SN dictionary
run.make_sn_dict(parent_folder, z, ra, dec)

# Import lightcurve from files (it must have already been converted using converting_data.py.
# Can eventually import a list of lcs if you want to - not right now though, it's not ready.
run1 = run.get_lc([parent_folder+name+'_lc.txt'])


# Fit SALT2 models to list of observed light curves.
print'attempting fits:'
run.fit_snlc(run1, parent_folder, child_folder='')
