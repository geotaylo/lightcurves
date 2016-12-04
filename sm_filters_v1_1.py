#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
GT 30/11/16

Program for creating and registering bandpasses of the SkyMapper Telescope
and Kepler Space Telescope in SNcosmo.
For use with sim_new_v1_1.py
Files must be in wd
"""

from astropy.io import ascii
import sncosmo

#  KST Filter
kst_data=ascii.read('Kepler2', data_start=1)
kst_wavelength = kst_data['col1']
kst_transmission = kst_data['col2']
kstband = sncosmo.Bandpass(kst_wavelength, kst_transmission, name='kst')   #  Creating bandpass
sncosmo.registry.register(kstband, 'kst', force=True)      #  Registering bandpass

#  SKYMAPPER v Filter
smv_data=ascii.read('SkyMapperv', data_start=1)
smv_wavelength = smv_data['col1']
smv_transmission = smv_data['col2']
smvband = sncosmo.Bandpass(smv_wavelength, smv_transmission, name='smv')   #  Creating bandpass
sncosmo.registry.register(smvband, 'smv', force=True)      #  Registering bandpass

#  SKYMAPPER r Filter
smr_data=ascii.read('SkyMapperr', data_start=1)
smr_wavelength = smr_data['col1']
smr_transmission = smr_data['col2']
smrband = sncosmo.Bandpass(smr_wavelength, smr_transmission, name='smr')   #  Creating bandpass
sncosmo.registry.register(smrband, 'smr', force=True)      #  Registering bandpass

#  SKYMAPPER g Filter
smg_data=ascii.read('SkyMapperg', data_start=1)
smg_wavelength = smg_data['col1']
smg_transmission = smg_data['col2']
smgband = sncosmo.Bandpass(smg_wavelength, smg_transmission, name='smg')   #  Creating bandpass
sncosmo.registry.register(smgband, 'smg', force=True)      #  Registering bandpass

#  SKYMAPPER i Filter
smi_data=ascii.read('SkyMapperi', data_start=1)
smi_wavelength = smi_data['col1']
smi_transmission = smi_data['col2']
smiband = sncosmo.Bandpass(smi_wavelength, smi_transmission, name='smi')   #  Creating bandpass
sncosmo.registry.register(smiband, 'smi', force=True)      #  Registering bandpass