"""
GT 09/12/16
Program for creating and registering bandpasses of the SkyMapper 
Telescope and Kepler Space Telescope in SNcosmo.
For use with run_sncosmo.py
Files must be in wd (or update path as needed)
"""

from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.lines import Line2D

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
                               name='Kepler')

    sncosmo.registry.register(kstband, 'Kepler', force=True)

        
    # SKYMAPPER v filter.

    smv_data = ascii.read('SkyMapperv', data_start=1)

    smv_wavelength = smv_data['col1']

    smv_transmission = smv_data['col2']

    smvband = sncosmo.Bandpass(smv_wavelength, smv_transmission,
                               name='SMv')

    sncosmo.registry.register(smvband, 'SMv', force=True)
        
    # SKYMAPPER r filter.

    smr_data = ascii.read('SkyMapperr', data_start=1)

    smr_wavelength = smr_data['col1']

    smr_transmission = smr_data['col2']

    smrband = sncosmo.Bandpass(smr_wavelength, smr_transmission,
                               name='SMr')

    sncosmo.registry.register(smrband, 'SMr', force=True)
    
    # SKYMAPPER g filter.

    smg_data = ascii.read('SkyMapperg', data_start=1)

    smg_wavelength = smg_data['col1']

    smg_transmission = smg_data['col2']

    smgband = sncosmo.Bandpass(smg_wavelength, smg_transmission,
                               name='SMg')
    sncosmo.registry.register(smgband, 'SMg', force=True)
        
    # SKYMAPPER i filter.

    smi_data = ascii.read('SkyMapperi', data_start=1)

    smi_wavelength = smi_data['col1']

    smi_transmission = smi_data['col2']
    
    smiband = sncosmo.Bandpass(smi_wavelength, smi_transmission,
                               name='SMi')
    
    sncosmo.registry.register(smiband, 'SMi', force=True)

    # Plot the filter curves
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig_width_pt = 240.0  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0 / 72.27  # Convert pt to inches
    golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean  # height in inches
    fig_size = [fig_width, fig_height]


    fig = plt.figure(figsize=(1.5 * fig_width, 1.5 * fig_height))

    plt.plot(kst_wavelength, kst_transmission, label=r'$Kepler$', linewidth=1)
    plt.plot(smg_wavelength, smg_transmission, label=r'$SkyMapper_g$', linewidth=1)
    plt.plot(smi_wavelength, smi_transmission, label=r'$SkyMapper_i$', linewidth=1)
    plt.plot(smr_wavelength, smr_transmission, label=r'$SkyMapper_r$', linewidth=1)

    plt.xlabel(r'Wavelength $\AA$', fontsize=13)
    plt.ylabel('Transmission', fontsize=13)


    legend_elements = [Line2D([0], [0], marker='o', color='w', label=r'$Kepler$',
                              markerfacecolor='g', markersize=5),
                       Line2D([0], [0], marker='o', color='w', label=r'$SkyMapper_g$',
                              markerfacecolor='b', markersize=5),
                       Line2D([0], [0], marker='o', color='w', label=r'$SkyMapper_r$',
                              markerfacecolor='y', markersize=5),
                       Line2D([0], [0], marker='o', color='w', label=r'$SkyMapper_i$',
                              markerfacecolor='m', markersize=5)
                       ]

    plt.legend(handles=legend_elements, loc='upper right', ncol=1, fancybox=True, framealpha=0.5)

    plt.savefig('filters.pdf', dpi=200, format='pdf', bbox_inches='tight')

    return

if __name__ == "__main__":
    register_filters()
