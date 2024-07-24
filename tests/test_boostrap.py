"""
Test the bootstrap module.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from python_scripts import bootstrap, flux, run

def test_bootstrap_resample_stopped():
    """
    Test the bootstrap_resample function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    test_run.init_markers()
    s_phi, s_theta, energy, weight = bootstrap.bootstrap_resample_stopped(test_run)
    assert len(s_phi) == 41927
    assert len(s_theta) == 41927
    assert len(energy) == 41927
    assert len(weight) == 41927

def test_bootstrap_resample_all():
    """
    Test the bootstrap_resample function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    test_run.init_markers()
    energy, energy0, stopped_weight, all_weight = bootstrap.bootstrap_resample_all(test_run)
    assert len(energy) < 524288
    assert len(energy0) == 524288
    assert len(stopped_weight) < 524288
    assert len(all_weight) == 524288
def test_calc_be_1d_array():
    """
    Test the calc_be_kde_1d_array function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join("tests", "output_plots")
    output_path = os.path.join(output_dir, 'be_array_1d.png')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points_1d=1000,
                       num_bootstraps=1024)
    test_run.flux.h_theta_1d = 0.1
    test_run.flux.energy_1d = flux.calc_energy_flux_1d(test_run)
    be_array_1d = bootstrap.calc_be_1d_array(test_run)

    # Make histogram of the BE array
    fig, ax = plt.subplots()
    ax.hist(be_array_1d, bins=100, density=True)
    ax.set_xlabel('BE [MW/m^2]')
    ax.set_ylabel('Probability')
    ax.set_title(f'Max energy flux 2d = {np.max(test_run.flux.energy_1d):.1e} MW/m^2')
    fig.savefig(output_path,
                bbox_inches='tight', dpi=300)
    plt.close(fig)

def test_calc_be_2d_array():
    """
    Test the calc_be_kde_2d_array function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join("tests", "output_plots")
    output_path = os.path.join(output_dir, 'be_array_2d.png')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points_2d=1000,
                       num_bootstraps=16)
    test_run.flux.h_phi = 1
    test_run.flux.h_theta_2d = 0.1
    test_run.flux.energy_2d = flux.calc_energy_flux_2d(test_run)
    be_array_2d = bootstrap.calc_be_2d_array(test_run)

    # Make histogram of the BE array
    fig, ax = plt.subplots()
    ax.hist(be_array_2d, bins=100, density=True)
    ax.set_xlabel('BE [MW/m^2]')
    ax.set_ylabel('Probability')
    ax.set_title(f'Max energy flux 2d = {np.max(test_run.flux.energy_2d):.1e} MW/m^2')
    fig.savefig(output_path,
                bbox_inches='tight', dpi=300)
    plt.close(fig)

def test_calc_be_total_array():
    """
    Test the calc_be_total_array function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join("tests", "output_plots")
    output_path = os.path.join(output_dir, 'be_array_total.png')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    test_run.init_markers()
    test_run.init_flux(num_bootstraps=128)
    be_array_total = bootstrap.calc_be_total_array(test_run)

    # Make histogram of the BE array
    fig, ax = plt.subplots()
    ax.hist(be_array_total, bins=100, density=True)
    ax.set_xlabel('BE [MW]')
    ax.set_ylabel('Probability')
    ax.set_title(f'Total energy flux = {test_run.flux.total_energy:.1e} MW')
    fig.savefig(output_path,
                bbox_inches='tight', dpi=300)
    plt.close(fig)
