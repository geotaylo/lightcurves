"""
GT 09/12/16
Program for creating and registering bandpasses of the SkyMapper 
Telescope and Kepler Space Telescope in SNcosmo.
For use with run_sncosmo.py
Files must be in wd (or update path as needed)
"""

from astropy.io import ascii

import sncosmo


def register_filters():
    
    """  Registers SkyMapper and Kepler filters to SNCosmo, for use in 
         simulating supernova observations.  MUST be executed at least 
         once before filters are available to use, but there's no harm
         in running it before every simulation.
         
         NOTE:  Current build requires specified ascii filters files to 
         be in working directory:
             - 'Kepler2cut'  << CUT to 9200A to fit in SALT2 range.
             - 'SkyMapperg'
             - 'SkyMapperi'
             - 'SkyMapperr'
             - 'SkyMapperv'
    """
    
    print 'Registering filters.\n'
    
    # KST filter.

    kst_data = ascii.read('Kepler2cut', data_start=1)

    kst_wavelength = kst_data['col1']

    kst_transmission = kst_data['col2']

    kstband = sncosmo.Bandpass(kst_wavelength, kst_transmission, 
                               name='kst')

    sncosmo.registry.register(kstband, 'kst', force=True)
        
    # SKYMAPPER v filter.

    smv_data = ascii.read('SkyMapperv', data_start=1)

    smv_wavelength = smv_data['col1']

    smv_transmission = smv_data['col2']

    smvband = sncosmo.Bandpass(smv_wavelength, smv_transmission,
                               name='smv') 

    sncosmo.registry.register(smvband, 'smv', force=True)
        
    # SKYMAPPER r filter.

    smr_data = ascii.read('SkyMapperr', data_start=1)

    smr_wavelength = smr_data['col1']

    smr_transmission = smr_data['col2']

    smrband = sncosmo.Bandpass(smr_wavelength, smr_transmission,
                               name='smr') 

    sncosmo.registry.register(smrband, 'smr', force=True)    
    
    # SKYMAPPER g filter.

    smg_data = ascii.read('SkyMapperg', data_start=1)

    smg_wavelength = smg_data['col1']

    smg_transmission = smg_data['col2']

    smgband = sncosmo.Bandpass(smg_wavelength, smg_transmission,
                               name='smg')
    sncosmo.registry.register(smgband, 'smg', force=True)
        
    # SKYMAPPER i filter.

    smi_data = ascii.read('SkyMapperi', data_start=1)

    smi_wavelength = smi_data['col1']

    smi_transmission = smi_data['col2']
    
    smiband = sncosmo.Bandpass(smi_wavelength, smi_transmission,
                               name='smi')
    
    sncosmo.registry.register(smiband, 'smi', force=True)
    
    return
