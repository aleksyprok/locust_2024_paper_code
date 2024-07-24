"""
This module contains unit tests for the prepare_profiles module.

The prepare_profiles module provides functions for preparing profiles for a given SPR string.
"""

import os
import netCDF4
from python_scripts import prepare_profiles

def test_get_cdf_filename():
    """
    Test case for the get_cdf_filename function.

    This test case verifies that the get_cdf_filename function returns the expected filename
    for a given SPR string.
    """
    spr_string = 'SPR-045-16'
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data")
    expected_filename = os.path.join(
        repo_path,
        'input_data',
        f'profiles_{spr_string}.CDF'
    )
    acutal_filename = prepare_profiles.get_cdf_filename(spr_string, input_dir)
    assert os.path.abspath(acutal_filename) == expected_filename

def test_read_cdf_file():
    """
    Test case for the read_cdf_file function.

    This test case verifies that the read_cdf_file function reads the CDF file correctly
    and returns the expected data.
    """
    spr_string = 'SPR-045-16'
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data")
    cdf_filename = os.path.join(
        input_dir,
        f'profiles_{spr_string}.CDF'
    )
    psin, ti, te, ne, nd, nt = prepare_profiles.read_cdf_file(cdf_filename)

    # Check that psin[0] is approximately 0 i.e. < 1e-4
    assert abs(psin[0]) < 1e-4
    # Check that ti, te and ne are greater in the core than the edge
    assert ti[0] > ti[-1]
    assert te[0] > te[-1]
    assert ne[0] > ne[-1]
    assert nd[0] > nd[-1]
    assert nt[0] > nt[-1]
    # Check ne > nd > nt
    assert ne[0] > nd[0] > nt[0]

def test_get_gfile():
    """
    Test case for the get_gfile function.

    This test case verifies that the get_gfile function returns the expected gfile path
    for a given gfile filename.
    """
    gfile_filename = 'SPR-045-16.eqdsk'
    current_dir = os.path.dirname(__file__)
    input_data_dir = os.path.join(current_dir, '..', 'input_data')
    gfile_path = os.path.join(input_data_dir, gfile_filename)
    gfile = prepare_profiles.get_gfile(gfile_path)
    assert abs(gfile.rmaxis - 4.381458618) < 1e-6
    assert abs(gfile.zmaxis - 0.0) < 1e-6

def test_calculate_number_of_impurities():
    """
    Test case for the calculate_number_of_impurities function.

    This test case verifies that the calculate_number_of_impurities function returns the expected
    number of impurities for a given CDF file.
    """
    spr_string = 'SPR-045-16'
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data")
    cdf_filename = os.path.join(
        input_dir,
        f'profiles_{spr_string}.CDF'
    )
    num_of_impurities = prepare_profiles.calculate_number_of_impurities(cdf_filename)
    assert num_of_impurities == 2

def test_get_impurity_data():
    """
    Test case for the get_impurity_data function.
    """
    spr_string = 'SPR-045-16'
    current_dir = os.path.dirname(__file__)
    input_data_dir = os.path.join(current_dir, '..', 'input_data')
    cdf_filename = os.path.join(input_data_dir, f'profiles_{spr_string}.CDF')
    impurity_names, zia, nim = prepare_profiles.get_impurity_data(cdf_filename, 2)
    assert impurity_names == ['AXe', 'AHe4']
    assert zia == [54.0, 2.0]
    assert nim[0][0] < nim[1][0]

def test_calculate_ion_fractions():
    """
    Test case for the calculate_ion_fractions function.
    """
    spr_string = 'SPR-045-16'
    current_dir = os.path.dirname(__file__)
    input_data_dir = os.path.join(current_dir, '..', 'input_data')
    cdf_filename = os.path.join(input_data_dir, f'profiles_{spr_string}.CDF')
    #pylint: disable=no-member
    with netCDF4.Dataset(cdf_filename, 'r') as profile_cdf:
        #pylint: enable=no-member
        ne = profile_cdf.variables['NE'][-1, :]
        nd = profile_cdf.variables['NID'][-1, :]
        nt = profile_cdf.variables['NIT'][-1, :]
        nim1 = profile_cdf.variables['NIM1'][-1, :]
        nim2 = profile_cdf.variables['NIM2'][-1, :]
    fd, ft, fi = prepare_profiles.calculate_ion_fractions(ne, nd, nt, [nim1, nim2])
    # check charge neutrality
    assert abs(fd + ft + fi[0] * 54 + fi[1] * 2 - 1) < 1e-6

def test_main():
    """
    Test case for the main function.
    """
    num_markers = 10
    spr_string = 'SPR-045-16'
    current_dir = os.path.dirname(__file__)
    input_data_dir = os.path.join(current_dir, '..', 'input_data')
    output_dir = os.path.join(current_dir, 'output_data')
    prepare_profiles.main(spr_string, input_data_dir, output_dir, num_markers)
    assert os.path.exists(os.path.join(output_dir, f'profile_{spr_string}_ne.dat'))
    assert os.path.exists(os.path.join(output_dir, f'profile_{spr_string}_Te.dat'))
    assert os.path.exists(os.path.join(output_dir, f'profile_{spr_string}_Ti.dat'))
    assert os.path.exists(os.path.join(output_dir, f'ion_info_{spr_string}.dat'))
    assert os.path.exists(os.path.join(output_dir, f'{spr_string}_markers_10.dat'))
