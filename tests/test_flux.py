"""
Test the flux module.
"""
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from python_scripts import flux, run

def test_calc_total_energy():
    """
    Test the calc_total_energy function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux()

    assert np.isclose(test_run.flux.total_energy,
                      test_run.log.total_stopped_power,
                      atol=1e-2)

def test_calc_energy_flux_1d():
    """
    Test the calc_energy_flux_1d function.
    """

    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join(repo_path, 'tests', 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points_1d=10**4)
    test_run.flux.energy_1d = flux.calc_energy_flux_1d(test_run, h_theta_1d=0.1)
    fig, ax = plt.subplots()
    ax.plot(test_run.flux.s_theta_1d, test_run.flux.energy_1d)
    ax.set_xlabel('s_theta [m]')
    ax.set_ylabel('Energy flux [MW/m^2]')
    fig.savefig(os.path.join(output_dir, 'axisymmetric_energy_flux_1d.png'))
    plt.close(fig)

def test_calc_energy_flux_2d():
    """
    Test the calc_energy_flux_2d function.
    """

    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join(repo_path, 'tests', 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points_2d=10**3)
    test_run.flux.energy_2d = flux.calc_energy_flux_2d(test_run, h_phi=1, h_theta_2d=0.1)

    # Make an imshow of the energy flux
    fig, ax = plt.subplots()
    ax.imshow(test_run.flux.energy_2d.T, origin='lower',
               extent=[test_run.wall.s_phi_min, test_run.wall.s_phi_max,
                       test_run.wall.s_theta_min, test_run.wall.s_theta_max])
    ax.set_xlabel('s_phi [m]')
    ax.set_ylabel('s_theta [m]')
    ax.set_title('Energy flux [MW/m^2]\n'
                 f'axisymmetric, h_theta_2d = {test_run.flux.h_theta_2d}'
                 f', h_phi = {test_run.flux.h_phi}')
    fig.savefig(os.path.join(output_dir, 'axisymmetric_energy_flux_2d.png'),
                bbox_inches='tight', dpi=300)
    plt.close(fig)

def test_calc_optimum_bandwidth_1d():
    """
    Test the calc_optimum_bandwidth_1d function.
    """

    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join(repo_path, 'tests', 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points_1d=1000)
    test_run.flux.h_theta_1d_array = np.logspace(-2, 0, 64)
    flux.calc_optimum_bandwidth_1d(test_run)

    fig, ax = plt.subplots()
    ax.plot(test_run.flux.h_theta_1d_array, test_run.flux.amise_1d)
    ax.plot(test_run.flux.h_theta_1d, test_run.flux.amise_1d.min(),
            'r.', markersize=10)
    ax.set_xscale('log')
    ax.set_xlabel('h_theta_1d [m]')
    ax.set_ylabel('AMISE')
    fig.savefig(os.path.join(output_dir, 'axisymmetric_amise_1d.png'))
    plt.close(fig)
    assert test_run.flux.h_theta_1d == 0.0372759372031494

def test_calc_optimum_bandwidth_2d():
    """
    Test the calc_optimum_bandwidth_2d function.
    """

    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join(repo_path, 'tests', 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points_2d=10**3)
    test_run.flux.h_phi_array = np.logspace(-2, 0, 8)
    test_run.flux.h_theta_2d_array = np.logspace(-2, 0, 8)
    flux.calc_optimum_bandwidth_2d(test_run)

    hx_mesh, hy_mesh = np.meshgrid(test_run.flux.h_phi_array,
                                   test_run.flux.h_theta_2d_array)
    fig, ax = plt.subplots()
    im = ax.pcolormesh(hx_mesh, hy_mesh, test_run.flux.amise_2d,
                       norm=LogNorm())
    ax.plot(test_run.flux.h_phi, test_run.flux.h_theta_2d,
            'r.', markersize=10)
    ax.set_title('Asymptotic mean integrated square error')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.colorbar(im, ax=ax)
    ax.set_xlabel(r'$h_\phi$ [m]')
    ax.set_ylabel(r'$h_\theta$ [m]')
    fig.savefig(output_dir + '/axisymmetric_amise_2d_array.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)

    h_theta_index = np.where(test_run.flux.h_theta_2d_array == test_run.flux.h_theta_2d)[0][0]
    h_phi_index = np.where(test_run.flux.h_phi_array == test_run.flux.h_phi)[0][0]
    axs[0].plot(test_run.flux.h_phi_array, test_run.flux.amise_2d[h_theta_index, :])
    axs[0].plot(test_run.flux.h_phi, test_run.flux.amise_2d[h_theta_index, h_phi_index],
                'r.', markersize=10)
    axs[0].set_xlabel(r'$h_\phi$ [m]')
    axs[0].set_ylabel(r'$AMISE$')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_title(r'$h_\theta$ = ' f'{test_run.flux.h_theta_2d} m')

    axs[1].plot(test_run.flux.h_theta_2d_array, test_run.flux.amise_2d[:, h_phi_index])
    axs[1].plot(test_run.flux.h_theta_2d, test_run.flux.amise_2d[h_theta_index, h_phi_index],
                'r.', markersize=10)
    axs[1].set_xlabel(r'$h_\theta$ [m]')
    axs[1].set_ylabel(r'$AMISE$')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_title(r'$h_\phi$ = ' f'{test_run.flux.h_phi} m')
    fig.savefig(output_dir + '/axisymmetric_amise_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)
