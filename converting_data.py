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
name = '2017i'

# Count offset - to subtract for galaxy.  Estimated from calibrated ground-based photometry on K2 drive.
galaxy_counts = 4425

# skiplines - check data file for each SN (inconsistent).
sl = [0]

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
df = df[df.Counts != 0]

# select relevant data
df_picky = df[['K2_Time', 'Counts', 'Error']]
df_picky.apply(pd.to_numeric)

t = Table.from_pandas(df_picky)

# add zpsys column and rename dflux in line with sncosmo reqs.  - correct filter entries
t['K2_Time'].name = 'time'
t['Counts'].name = 'flux'
t['zp'] = 25.47
t['zpsys'] = 'ab'
t['Error'].name = 'flux_err'
t['filter'] = 'kst'

# Convert kepler time to MJD
# Assume K2_time = BJD - 2454833 (source = Brad, https://archive.stsci.edu/prepds/kegs/#dataformat) - we're going to
# leave it in BJD and convert everything else to BJD as well
t['time'] = t['time'] + 2454833


# Subtract galaxy background
sky_percent = (100/t['flux'])*galaxy_counts
t['flux'] = t['flux'] - galaxy_counts
t['flux_err'] = np.abs(t['flux_err']*(100-sky_percent)/100)

# Maybe save in different folder?
sncosmo.write_lc(t, folder+name+'/'+name+'_lc.txt')

# Test conversion by reading in and plotting data
data = sncosmo.read_lc(folder+name+'/'+name+'_lc.txt')
fig = sncosmo.plot_lc(data,format='png')
plt.savefig(folder+name+'/'+name+'_data.png')
