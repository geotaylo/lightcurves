# lightcurves
Codes for simulating and fitting light curves using SNCosmo, for ANU Summership 2016-2017


Code purpose:
* This project is designed to simulate observations of supernovae of known parameters, and then fit a light curve model to the   observations to approximately obtain these parameters.  This allows insight into the recovery of SN parameters under a range of observing conditions.
* Simulating observations:  The code randomly generate parameters of supernovae that we could reasonably expect to observe as part of the   SkyMapper Transient Survey.  Observing properties are generated based on user-inputted filters and cadence.  Observations are then combined with the SN parameters and a 'blank' SALT2 model using SNCosmo's 'realize_lc()' method to simulate an observed light curve, which is stored as a observed_lc.txt file (as are the 'true' parameters of each SN, in true_parameters.txt).  The features of the simulation (SN parameters, coordinates, and relevant observing properties) are saved in sn_dict.pkl so that the same set of SN can be 'observed' again with different properties (e.g. cadence, filter set) if desired.
* Fitting a light curve:  Using observed light curves (either pre-existing or returned from the simulation above), the code uses SNCosmo's   'mcmc_lc()' method to fit a SALT2 model to the observations.  Prior to this, the extinction of the model is set using SNCosmo's dustmap   funtionality.  The plot of model vs. data is saved as fitted_lc.pdf (as are the 'fitted' parameters of each SN, in fitted_parameters.txt).


Installation requirements:
* astropy: http://www.astropy.org/
* emcee: http://dan.iel.fm/emcee/current/
* matplotlib: http://matplotlib.org/
* sfdmap: https://github.com/kbarbary/sfdmap
* sncosmo: http://sncosmo.readthedocs.io/en/v1.4.x/install.html


Basic instructions:

1. Complete required installations.

2. Edit the wrapper.py code for your desired simulations.  Available methods include:

   *  simulate_lc(nSNe=0, cadence=4, kpass=False, folder='TestFiles/', tmin=57754, tmax=58118, zmin=0.001, zmax=0.1, properties=''): Returns list of 'observed' lightcurves.  
   SN parameters are saved in true_parameters.txt; simulated observations for SN 'n' are saved in observed_lc_n.txt; simulation features are saved in sn_dict.pkl.
   Inputs:
   - nSNe: number of SN to simulate (if 0, takes expected SN rate from SNCosmo based on SkyMapper field of view).
   - cadence: number of days between simulated observations.
   - kpass: whether the kepler filter is included in simulation.
   - folder: directory to store outputs.
   - tmin & tmax: possible interval covered by simulation, in mjd (default is all of 2017).
   - zmin & zmax: possible redshift range of SN.
   - properties: file containing set of features to 're-observe'.  If blank, random features are generated.
   
   * get_lc(filelist): If you're fitting models to observations from an existing ascii file (rather than simulating observations above), this will return the a list of lightcurves to pass to the fitting method.
   Inputs:
   - filelist: Python list of observation filename strings.
   
   * fit_snlc(lightcurve, folder='TestFiles/', properties=''): Fits SALT2 models to list of observed light curves.
   Fitted model parameters are saved in fitted_parameters.txt; plots of simulated observations vs fitted model are saved in fitted_lc_n.pdf.
   Inputs:
   - lightcurve: list of observed lightcurves
   - folder: directory to store outputs.
   - properties: file containing set of features of simulation, used to pass coordinates to dust maps.  If blank (i.e. if using get_lc() instead of a simulation), random SN cordinates are generated.
   
3. Save wrapper.py file and run from terminal.  In current directory: python wrapper.py

4. You will need an internet connection (to load dustmaps) and time - each SN model takes ~3 minutes to fit.


About files:

* filters.py reads wavelength/transmission files for SkyMapper and Kepler filters, and create and registers bandpasses for use in SNCosmo   code.  It's currently designed to work only for the following filters:
  - 'Kepler2cut', where cut refers to the upper wavelength limit being reduced from 9500Å to 9200Å in order to fit in SALT2 model's wavelength range.  NOTE that this cut currently needs to be performed manually in the ascii file, if using updated filter sets.
  - 'SkyMapperg'
  - 'SkyMapperi'
  - 'SkyMapperr'
  - 'SkyMapperv'
  All of these filters then become available for the simulation code to use as desired.  To include additional filters, register them in filters.py and then add them to filter lists at top of run_sncosmo.py
  
* run_sncosmo.py does all of the heavy lifting, and should only be changed with caution.  A brief overview of the code follows:
  - Cosmology and telescope parameters are set (field of view, filter sets, h0, ...)
  - 

* wrapper.py is a wrapper script, used to run a complete simulation and fitting of a set of SN.  The example script currently runs simulated observations of the same 5 SN at a cadence of 4 days, both with and without kepler filters, and then fits SALT2 models to the observations.

