"""
Module for testing the dt_fusion module.
"""

import os
import netCDF4
from python_scripts import my_gfile_reader, dt_fusion

def test_calc_fusion_power():
    """
    Function to test the calc_fusion_power function.
    """

    spr_string = 'SPR-045-16'
    current_dir = os.path.dirname(__file__)
    input_data_dir = os.path.join(current_dir, '..', 'input_data')
    cdf_filename = os.path.join(input_data_dir, f'profiles_{spr_string}.CDF')
    #pylint: disable=no-member
    with netCDF4.Dataset(cdf_filename, 'r') as profile_cdf:
        #pylint: enable=no-member
        psin = profile_cdf.variables['XPSI'][-1, :]
        ti = profile_cdf.variables['TI'][-1, :]
        nd = profile_cdf.variables['NID'][-1, :]
        nt = profile_cdf.variables['NIT'][-1, :]
    gfile_filename = f'{spr_string}.eqdsk'
    gfile_path = os.path.join(input_data_dir, gfile_filename)
    gfile = my_gfile_reader.getGfile(gfile_path)
    fusion_power, _ = dt_fusion.calc_fusion_power(psin, nd, nt, ti, gfile)
    assert abs(fusion_power - 1.66e9) < 0.02e9
