"""
Test module for the log.py module.
"""

import os
from python_scripts import log

def test_log_init():
    """
    Test the Log class initialization.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_path = os.path.join(repo_path, "input_data")

    log_path_axisymmetric = os.path.join(input_path, 'LOCUST_SPR-045-14_OutputFiles/axisymmetric/'
                                                      'gpu-q-41/LOG_13-12-2023_16-51-52.811.out')
    log_path_ripple = os.path.join(input_path, 'LOCUST_SPR-045-14_OutputFiles/vary_ripple/'
                                                'gpu-q-18/LOG_15-08-2023_20-34-29.636.out')
    log_path_rmp = os.path.join(input_path, 'LOCUST_SPR-045-14_OutputFiles/vary_dave_bplasma/'
                                             'gpu-q-9/LOG_11-08-2023_06-06-19.733.out')
    log_path_rwm = os.path.join(input_path, 'LOCUST_SPR-045-14_OutputFiles/rwm_tesla_scan/lrdn1185/'
                                            'LOG_20-11-2023_16-04-07.624.out')
    log_axisymmetric = log.Log(log_path_axisymmetric)
    log_ripple = log.Log(log_path_ripple)
    log_rmp = log.Log(log_path_rmp)
    log_rwm = log.Log(log_path_rwm)

    assert log_axisymmetric.log_path == log_path_axisymmetric
    assert log_ripple.log_path == log_path_ripple
    assert log_rmp.log_path == log_path_rmp
    assert log_rwm.log_path == log_path_rwm

    assert log_axisymmetric.analytic_ripple is False
    assert log_ripple.analytic_ripple is True
    assert log_rmp.analytic_ripple is False
    assert log_rwm.analytic_ripple is False

    assert log_axisymmetric.bplasma_file is None
    assert log_ripple.bplasma_file is None
    assert log_rmp.bplasma_file == 'BPLASMA_efcc_response=1_current=100_100x200_phase=225.0_n2'
    assert log_rwm.bplasma_file == 'BPLASMA_cylindrical_tesla_G=0.0_bscale=1_n1'

    assert log_axisymmetric.bplasma is False
    assert log_ripple.bplasma is False
    assert log_rmp.bplasma is True
    assert log_rwm.bplasma is True

    assert log_axisymmetric.axisymmetric is True
    assert log_ripple.axisymmetric is False
    assert log_rmp.axisymmetric is False
    assert log_rwm.axisymmetric is False

    assert log_axisymmetric.rcoil is None
    assert abs(log_ripple.rcoil - 8.50) < 1e-10
    assert log_rmp.rcoil is None
    assert log_rwm.rcoil is None

    assert log_axisymmetric.ncoil is None
    assert log_ripple.ncoil is 12
    assert log_rmp.ncoil is None
    assert log_rwm.ncoil is None

    assert log_axisymmetric.coil_set is None
    assert log_ripple.coil_set is None
    assert log_rmp.coil_set == 'exterior_rmp'
    assert log_rwm.coil_set is 'rwm_acc'

    assert log_axisymmetric.bplasma_n is None
    assert log_ripple.bplasma_n is None
    assert log_rmp.bplasma_n == 2
    assert log_rwm.bplasma_n == 1

    assert log_axisymmetric.rmp_current is None
    assert log_ripple.rmp_current is None
    assert log_rmp.rmp_current == 100
    assert log_rwm.rmp_current is None

    assert log_axisymmetric.rmp_phase is None
    assert log_ripple.rmp_phase is None
    assert log_rmp.rmp_phase == 225
    assert log_rwm.rmp_phase is None#

    assert log_axisymmetric.rmp_response is None
    assert log_ripple.rmp_response is None
    assert log_rmp.rmp_response is True
    assert log_rwm.rmp_response is None

    assert log_axisymmetric.rwm_gain is None
    assert log_ripple.rwm_gain is None
    assert log_rmp.rwm_gain is None
    assert log_rwm.rwm_gain == 0

    assert log_axisymmetric.rwm_bscale is None
    assert log_ripple.rwm_bscale is None
    assert log_rmp.rwm_bscale is None
    assert log_rwm.rwm_bscale == 1

    assert log_axisymmetric.total_stopped_power == 0.26120333
    assert log_ripple.total_stopped_power == 0.27798763
    assert log_rmp.total_stopped_power == 7.7370629
    assert log_rwm.total_stopped_power == 0.26159607

    assert log_axisymmetric.total_stopped_power_error == 0.0025217316
    assert log_ripple.total_stopped_power_error == 0.0028677862
    assert log_rmp.total_stopped_power_error == 0.043398268
    assert log_rwm.total_stopped_power_error is None

    # assert log_axisymmetric.max_energy_flux == 11.74194
    # assert log_ripple.max_energy_flux == 4.61870
    # assert log_rmp.max_energy_flux == 1.95909
    # assert log_rwm.max_energy_flux == 11.74194

    assert log_axisymmetric.pinj == 338
    assert log_ripple.pinj == 338
    assert log_rmp.pinj == 338
    assert log_rwm.pinj == 338

    assert log_axisymmetric.simulation_time == 29869.8
    assert log_ripple.simulation_time == 0.410031E+05
    assert log_rmp.simulation_time == 0.369179E+05
    assert log_rwm.simulation_time == 0.412703E+05
