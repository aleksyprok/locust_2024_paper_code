"""
Module contains routines for plotting histograms of the final pitch of the
stopped markers, final energy and stopped time as well as
the thermalization time.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from KDEpy import FFTKDE
from python_scripts import run

def plot_energy_histogram(run0, plot_dir):
    """
    Plots histograms of the final energy.
    """
    fig, ax = plt.subplots()
    ax.hist(run0.markers.stopped.energy / 1e6 * 6.6464731e-27 / 1.602e-19,
            weights=run0.markers.stopped.weight,
            histtype='step',
            bins=100)
    ax.set_xlabel('Energy [MeV]')
    ax.set_ylabel('Weighted Number of Markers')
    ax.set_title('Energy distribution of stopped markers')
    fig.savefig(os.path.join(plot_dir, 'stopped_energy_distribution.png'),
                bbox_inches='tight',
                dpi=300)
    
def plot_thermalization_time_histogram(run0, plot_dir):
    """
    Plots histograms of the thermalisation time.
    """
    fig, ax = plt.subplots()
    ax.hist(run0.markers.thermal.t,
            weights=run0.markers.thermal.weight,
            histtype='step',
            bins=100)
    ax.set_xlabel('Thermalization Time [s]')
    ax.set_ylabel('Weighted Number of Markers')
    ax.set_xscale('log')
    fig.savefig(os.path.join(plot_dir, 'thermalisation_time.png'),
                bbox_inches='tight',
                dpi=300)
    
def plot_kde_of_stopped_marker_initial_positions(run0, plot_dir):
    """
    Plots the KDE of the initial position of stopped markers.
    """

    bw = 0.1
    num_grid_points = 1000
    fft_kde = FFTKDE(kernel='gaussian', bw=bw)
    coords = np.vstack([run0.markers.stopped.r0, run0.markers.stopped.z0]).T
    poloidal_grid, pdf = \
        fft_kde.fit(coords, weights=run0.markers.stopped.weight *
                                    run0.markers.stopped.energy).evaluate(num_grid_points)
    pdf = pdf.reshape(num_grid_points, num_grid_points).T
    r_unique, z_unique = np.unique(poloidal_grid[:, 0]), np.unique(poloidal_grid[:, 1])
    grid_r, grid_z = np.meshgrid(r_unique, z_unique)

    fig, ax = plt.subplots(1, 1)
    cs = ax.pcolormesh(grid_r, grid_z, pdf)
    plt.colorbar(cs)
    ax.plot(run0.gfile.R_bnd, run0.gfile.Z_bnd, 'r')
    ax.plot(run0.wall.r, run0.wall.z, 'k')
    ax.set_aspect('equal')
    ax.set_xlim([min(run0.gfile.R_bnd), max(run0.gfile.R_bnd)])
    ax.set_ylim([min(run0.gfile.Z_bnd), max(run0.gfile.Z_bnd)])
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    fig.savefig(os.path.join(plot_dir, 'stopped_marker_initial_positions.png'),
                bbox_inches='tight',
                dpi=300)
    
def check_toroidal_symmetry(run0, plot_dir):
    """
    Checks if the markers are toroidally symmetric.
    """
    fig, ax = plt.subplots(1, 3)
    fig_size = fig.get_size_inches()
    fig.set_size_inches(fig_size[0] * 3, fig_size[1] * 1)
    ax[0].hist(run0.markers.stopped.phi0 * 180 / np.pi,
               weights=run0.markers.stopped.weight * run0.markers.stopped.energy,
               bins=100)
    ax[0].set_xlabel('Initial Toroidal Angle [deg]')
    ax[0].set_ylabel('Weighted Number of Markers')
    ax[0].set_title('Initial Toroidal Angle Distribution')
    ax[1].hist(run0.markers.stopped.phi * 180 / np.pi,
               weights=run0.markers.stopped.weight * run0.markers.stopped.energy,
               bins=100)
    ax[1].set_xlabel('Final Toroidal Angle [deg]')
    ax[1].set_ylabel('Weighted Number of Markers')
    ax[1].set_title('Final Toroidal Angle Distribution')
    ax[2].hist(run0.markers.stopped.s_phi,
               weights=run0.markers.stopped.weight * run0.markers.stopped.energy,
               bins=100)
    ax[2].set_xlabel('Final toroidal distance [m]')
    ax[2].set_ylabel('Weighted Number of Markers')
    ax[2].set_title('Final toroidal distance distribution')

    fig.savefig(os.path.join(plot_dir, 'toroidal_symmetry.png'),
                bbox_inches='tight',
                dpi=300)
    
if __name__ == '__main__':
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    runs_path = os.path.join(repo_path, "output_data",
                             "spr_068_spr_045_ripple_scan")
    all_runs = run.create_runs_list(runs_path)
    for run_i in all_runs:
        plots_dir_i = os.path.join(run_i.dir_path, "plots")
        if "SPR-045-14" in run_i.log.eqdsk_fname:
            spr_string = "SPR-045-14"
            wall_path = os.path.join(repo_path, "input_data", "SPP-001_wall.dat")
            gfile_path = os.path.join(repo_path, "input_data", "SPR-045-14.eqdsk")
        elif "SPR-045-16" in run_i.log.eqdsk_fname:
            spr_string = "SPR-045-16"
            wall_path = os.path.join(repo_path, "input_data", "SPP-001_wall.dat")
            gfile_path = os.path.join(repo_path, "input_data", "SPR-045-16.eqdsk")
        elif "SPR-068-7" in run_i.log.eqdsk_fname:
            spr_string = "SPR-068-7"
            wall_path = os.path.join(repo_path, "input_data", "SPR-068_wall.dat")
            gfile_path = os.path.join(repo_path, "input_data", "SPR-068-7.eqdsk")
        else:
            raise ValueError("Could not determine run type")
        if run_i.log.axisymmetric:
            output_dir = os.path.join(repo_path, "plots", "spr_068_spr_045_axisymmetric_scan")
            output_dir_i = os.path.join(output_dir, spr_string)
            run_i.init_gfile(gfile_path)
            run_i.init_wall(wall_path)
            run_i.init_markers()
            plot_energy_histogram(run_i, output_dir_i)
            plot_thermalization_time_histogram(run_i, output_dir_i)
            plot_kde_of_stopped_marker_initial_positions(run_i, output_dir_i)
            check_toroidal_symmetry(run_i, output_dir_i)
        else:
            output_dir = os.path.join(repo_path, "plots", "spr_068_spr_045_ripple_scan")
            output_dir_i = os.path.join(output_dir, spr_string,
                            f"rcoil_{run_i.log.rcoil}_ncoil_{run_i.log.ncoil}")
            run_i.init_gfile(gfile_path)
            run_i.init_wall(wall_path)
            run_i.init_markers()
            plot_energy_histogram(run_i, output_dir_i)
            plot_thermalization_time_histogram(run_i, output_dir_i)
            plot_kde_of_stopped_marker_initial_positions(run_i, output_dir_i)
            check_toroidal_symmetry(run_i, output_dir_i)
