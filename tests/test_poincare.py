"""
Test the poincare module.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from python_scripts import poincare, run

def test_read_poincare_file():
    """
    Test the read_poincare_file function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    poincare_path = os.path.join(repo_path, 'output_data', 'poincare_test',
                                 'Poincare_map_12-07-2024_16-54-18.698_formatted.dat')
    plot_dir = os.path.join(repo_path, 'tests', 'output_plots')
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, 'poincare_plot_rz.png')
    r_array, z_array = poincare.read_poincare_file(poincare_path)
    fig, ax = plt.subplots()
    fig_size = fig.get_size_inches()
    fig_size *= 50
    fig.set_size_inches(fig_size)
    ax.scatter(r_array, z_array,
               marker='.')
    ax.set_aspect('equal')
    fig.savefig(plot_path,
                bbox_inches='tight',
                dpi=100)

def test_map_grid_to_rz_array():
    """
    Test the map grid to rz array function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    poincare_path = os.path.join(repo_path, 'output_data', 'poincare_test',
                                 'Poincare_map_12-07-2024_16-54-18.698_formatted.dat')
    poincare.read_poincare_file(poincare_path)
    values = []
    with open(poincare_path, 'r', encoding='utf-8') as poincare_file:
        poincare_file.readline()
        line = poincare_file.readline()
        nrp, nzp = [int(x) for x in line.split()[:2]]
        for _ in range(6):
            poincare_file.readline()
        line = poincare_file.readline()
        r_min, r_max, z_min, z_max = [float(x) for x in line.split()]
        for _ in range(nrp * nzp // 6):
            line = poincare_file.readline()
            values.extend(map(float, line.split()))
    poincare_map = np.array(values).reshape(nzp, nrp)

    r_array, z_array = poincare.map_grid_to_rz_array(poincare_map, r_min, r_max, z_min, z_max)

    assert len(r_array) == 194352
    assert len(z_array) == 194352

def test_calc_theta_psin_values():
    """
    Test the calc psi theta values function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    run_dir = os.path.join(repo_path, 'output_data', 'poincare_test')
    run_tag = '12-07-2024_16-54-18.698'
    run0 = run.Run(run_dir, run_tag)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    run0.init_gfile(gfile_path)
    run0.init_poincare()

    assert len(run0.poincare.psin_array) == 194352
    assert len(run0.poincare.theta_array) == 194352

    fig, ax = plt.subplots()
    fig_size = fig.get_size_inches()
    fig_size *= 5
    fig.set_size_inches(fig_size)
    ax.scatter(run0.poincare.theta_array, run0.poincare.psin_array,
               marker='.')
    fig.savefig(os.path.join(repo_path, 'tests', 'output_plots', 'poincare_psi_theta.png'),
                bbox_inches='tight', dpi=100)
