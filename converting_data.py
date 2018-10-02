import sncosmo
import heaven_and_earth as run
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
name = '2018j'

# redshift (from K2 spreadsheet)
z = 0.01

# Count offset - to subtract for galaxy.  Estimated from calibrated ground-based photometry on K2 drive.
galaxy_counts = 9330

# skiplines - check data file for each SN (inconsistent).
sl = [1]

# Coordinate frame
frame = 'fk5j2000'


# CODE -----------------------------------------------------------------------------------------------------------------
# Need to register filters to get zps
filters.register_filters()


# Read in csv and convert to ascii -------------------------------------------------------------------------------------

# Name of file to convert - eventually loop over all files in directory
target = folder+name+'/'+name+'.csv'

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
sncosmo.write_lc(t, folder+name+'/'+name+'_lc.txt')
