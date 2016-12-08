# lightcurves
Codes for simulating and fitting light curves using SNCosmo, for ANU Summership 2016-2017

Code purpose:
* This project is designed to simulate observations of supernovae of known parameters, and then fit a light curve model to the               observations to approximately obtain these parameters.  This allows insight into the recovery of SN parameters under a range of           observing conditions.
* Simulating observations:  The code randomly generate parameters of supernovae that we could reasonably expect to observe as part of the   SkyMapper Transient Survey.  Observing properties are generated based on user-inputted filters and cadence.  Observations are then         combined with the SN parameters and a 'blank' SALT2 model using SNCosmo's 'realize_lc()' method to simulate an observed light curve,       which is stored as a observed_lc.txt file (as are the 'true' parameters of each SN, in true_parameters.txt).
* Fitting a light curve:  Using observed light curves (either pre-existing or returned from the simulation above), the code uses SNCosmo's   'mcmc_lc()' method to fit a SALT2 model to the observations.  Prior to this, the extinction of the model is set using SNCosmo's dustmap   funtionality.  The plot of model vs. data is saved as fitted_lc.pdf (as are the 'fitted' parameters of each SN, in                         fitted_parameters.txt).

Installation requirements:

Basic instructions:
* After completing the required installations listed above, you can run the simulation code directly from your command line, or use the     wrapper.py script which performs an example simulation.

About files:
* filters.py reads wavelength/transmission files for SkyMapper and Kepler filters, and create and registers bandpasses for use in SNCosmo   code.  It's currently designed to work only for the following filters:
  - 'Kepler2cut', where cut refers to the upper wavelength limit being reduced from 9500Å to 9200Å in order to fit in SALT2 model's           wavelength range.
  - 'SkyMapperg'
  - 'SkyMapperi'
  - 'SkyMapperr'
  - 'SkyMapperv'
  All of these filters then become available for the simulation code to use as desired.  To include additional filters, register them in     filters.py and then add them to run_sncosmo.py
* run_sncosmo.py
* wrapper.py
