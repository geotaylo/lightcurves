#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
GT 30/11/2016

Test run of a wrapper script to take inputs and run section A of sncosmo situation
.
"""

### Keep updating as needed - might be able to change to cmd eventually
import sim_new as sim

#  Simulate observed lightcurve
#  Defaults: nSNe=0, tmin=57754 (01Jan2017), tmax=58118 (31Dec2017), zmin=0.001, zmax=0.07
#  Note: nSNe=0 obtains number of SN and z values from sncosmo.zdist based on SkyMapper field of view
#  and a flat lambdacdm cosmology.
#  Outputs list of light curves
test1 = sim.gen_lc(nSNe=2)

#  Get list of lightcurves from files - if not simulating above
test = sim.get_lc(['TestFiles/observed_lc_1.txt', 'TestFiles/observed_lc_2.txt'])


#  Fit SALT2 model to (multiple) light curves
#  'Test' parameter = list of light curves
res = sim.fit_sn_lc(test)
