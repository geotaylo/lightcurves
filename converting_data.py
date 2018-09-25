import sncosmo
import run_sncosmo_k2fields as run
import pandas as pd
from astropy.table import Table
import filters
import matplotlib.pyplot as plt
import numpy as np
import sfdmap


# INPUTS ---------------------------------------------------------------------------------------------------------------
# Only change these - unless something goes horribly wrong, in which case you're on your own.
# Important - change NaN to 0 entries first, because Python is a stupid language.
# Time is in BJD!

# Which dust law?
dust = sncosmo.CCM89Dust()

# change location to map
dustmap = sfdmap.SFDMap("C:\Users\gltay\Documents\ANU 2018\Honours Project\lightcurves\sfddata-master\sfddata-master")

# Folder to read/write to.
folder = 'Real_data/KEGS/'

# Name of SN
name = '2017k'

# redshift (from K2 spreadsheet)
z = 0.13182

# Count offset - to subtract for galaxy.  Estimated from calibrated ground-based photometry on K2 drive.
galaxy_counts = 1288

# skiplines - check data file for each SN (inconsistent).
sl = [0,2]

# Coordinate frame
frame = 'fk5j2000'

# ra, in degrees, minutes, seconds
ra_d = 01.
ra_m = 58.
ra_s = 06.

ra = run.dms_to_deg(ra_d, ra_m, ra_s)

# dec, in degrees, minutes, seconds
dec_d = 13.
dec_m = 46.
dec_s = 44.2

dec = run.dms_to_deg(dec_d, dec_m, dec_s)



# CODE -----------------------------------------------------------------------------------------------------------------
# Need to register filters to get zps
filters.register_filters()


# Read in csv and convert to ascii -------------------------------------------------------------------------------------

# Name of file to convert - eventually loop over all files in directory
target = folder+name+'.csv'

# Load data frame
# Check skip rows are correct!
df = pd.read_csv(target, skiprows=sl)
df = df[df.Raw != 0]

# select relevant data
df_picky = df[['K2_Time', 'Raw', 'Error']]
df_picky.apply(pd.to_numeric)

# Remove nan rows
# # df_picky.replace(["NaN", 'NaT'], np.nan, inplace = True)
# df_numbers = df_picky.dropna()

#df_numbers = df_picky[~df_picky.isin(['NaN', 'NaT']).any(axis=1)]
t = Table.from_pandas(df_picky)

# add zpsys column and rename dflux in line with sncosmo reqs.  - correct filter entries
t['K2_Time'].name = 'time'
t['Raw'].name = 'flux'
t['zp'] = 25.47
t['zpsys'] = 'ab'
t['Error'].name = 'flux_err'
t['filter'] = 'kst'


# Convert kepler time to MJD
# Assume K2_time = BJD - 2454833 (source = Brad, https://archive.stsci.edu/prepds/kegs/#dataformat) - we're going to
# leave it in BJD and convert everything else to BJD as well
t['time'] = t['time'] - 2454833


# Convert counts, error, and galaxy offset to flux
ab = sncosmo.get_magsystem('ab')
# This is old conversion stuff, don't worry about it.

# First convert to KST_mag = -2.5 log_10 (counts) + 25.47 (from Brad)
# t['flux'] = -2.5*np.log10(t['flux'])+25.47
# t['flux_err'] = -2.5*np.log10(t['flux_err'])+25.47
# galaxy_mag = -2.5*np.log10(galaxy_counts)+25.47
# Then to flux = zpband*10^(-0.4*mag) (from SNCOSMO)
# t['flux'] = ab.band_mag_to_flux(t['flux'], 'kst')
# t['flux_err'] = ab.band_mag_to_flux(t['flux_err'], 'kst')
# galaxy_flux = ab.band_mag_to_flux(galaxy_mag, 'kst')
# Finally, subtract galaxy flux from counts
sky_percent = (100/t['flux'])*galaxy_counts
t['flux'] = t['flux'] - galaxy_counts
t['flux_err'] = np.abs(t['flux_err']*(100-sky_percent)/100)


# Maybe save in different folder?
sncosmo.write_lc(t, folder+name+'_lc.txt')

# Fitting --------------------------------------------------------------------------------------------------------------

# run1 = run.get_lc([folder+name+'_lc.txt')
# print'attempting well-sampled SM fits:'
# run.fit_snlc(run1, parent_folder, child_folder=child_folder_1)
dat = sncosmo.read_lc(folder+name+'_lc.txt')

# create a model
model = sncosmo.Model(source='salt2',
                      effects=[dust, dust],
                      effect_names=['host', 'mw'],
                      effect_frames=['rest', 'obs'])


# Fit for dust
ebv = dustmap.ebv(ra, dec, frame = frame, unit='degree')

model.set(mwebv=ebv)


# Don't fit redshift - it's known
model.set(z=z, x0=0.1, x1=0.1, c=0.1)

result, fitted_model = sncosmo.fit_lc(
    dat, model,
    ['t0', 'x0', 'x1', 'c'],  # parameters of model to vary,
    # bounds={'z':(0.1,0.15)},
    guess_amplitude=False,
    minsnr=3)  # bounds on parameters (if any)

fig = sncosmo.plot_lc(dat, model=fitted_model,
                    errors=result.errors, format='png'
                    )

plt.savefig(folder+name+'_fit.png')
plt.close(fig)

fig2 = sncosmo.plot_lc(dat, format='png')

plt.savefig(folder+name+'_data.png')
plt.close(fig2)


