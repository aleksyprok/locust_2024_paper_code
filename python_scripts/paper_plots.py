"""
This file contains functions used to make the plots that appears in the IAEA FEC 2024
paper.

Note that the code in this file assumes that the LOCUST runs have completed and
the data files are in the output_data/FEC_2024 directory.
"""
import os
import pickle
import time
from matplotlib.colors import Normalize, LinearSegmentedColormap
# pylint: disable=no-name-in-module
from matplotlib.cm import ScalarMappable, viridis
# pylint: enable=no-name-in-module
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from python_scripts import bootstrap, dt_fusion, flux, my_gfile_reader, paper_plots_extra, \
                           paper_plots_3d, prepare_profiles, ripple_check, run, wall

plt.rcParams.update({'font.size': 14})
REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_DIRECTORY = os.path.join(REPOSITORY_PATH, "output_data",
                                "FEC_2024")
GFILE_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPR-045-16.eqdsk")
WALL_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
NUM_GRID_POINTS_1D = 10**5
NUM_GRID_POINTS_2D = 10**3
NUM_BOOTSTRAPS = 128
CONFIDENCE_LEVEL = 0.95
H_PHI_ARRAY = np.logspace(-2, 0, 8)
H_THETA_1D_ARRAY = np.logspace(-3, 0, 64)
H_THETA_2D_ARRAY = np.logspace(-2, 0, 8)
SPECIAL_NODES = (16, 55, 138, 178)

def calc_energy_flux(run_i, output_dir_i,
                     remap_phi_n=None,
                     previous_run=None,
                     free_space=True,
                     wall_path=WALL_PATH,
                     special_nodes=SPECIAL_NODES):
    """
    This function calculates the energy flux for a given run.
    And plots some of the data along the way.
    """
    os.makedirs(output_dir_i, exist_ok=True)
    run_i.init_wall(wall_path,
                    special_nodes=special_nodes)
    run_i.init_markers(remap_phi_n=remap_phi_n)
    run_i.init_flux(num_grid_points_1d=NUM_GRID_POINTS_1D,
                    num_grid_points_2d=NUM_GRID_POINTS_2D,
                    num_bootstraps=NUM_BOOTSTRAPS,
                    h_phi_array=H_PHI_ARRAY,
                    h_theta_1d_array=H_THETA_1D_ARRAY,
                    h_theta_2d_array=H_THETA_2D_ARRAY)
    # run_i.flux.h_phi = 1
    # run_i.flux.h_theta_2d = 0.01
    flux.calc_optimum_bandwidth_1d(run_i)
    flux.calc_optimum_bandwidth_2d(run_i)
    if previous_run is not None:
        # If the optimum bandwidth finder fails then use the value from the previous run.
        if np.isclose(run_i.flux.h_phi, H_PHI_ARRAY[-1], atol=1e-3) and \
            np.isclose(run_i.flux.h_theta_2d, H_THETA_2D_ARRAY[-1], atol=1e-3):
            run_i.flux.h_phi = previous_run.flux.h_phi
            run_i.flux.h_theta_2d = previous_run.flux.h_theta_2d
    paper_plots_extra.plot_amise_1d(run_i, output_dir_i)
    paper_plots_extra.plot_amise_2d(run_i, output_dir_i)
    run_i.flux.energy_1d = flux.calc_energy_flux_1d(run_i)
    run_i.flux.energy_2d = flux.calc_energy_flux_2d(run_i)
    paper_plots_extra.plot_energy_flux_1d(run_i, output_dir_i)
    paper_plots_extra.plot_energy_flux_2d(run_i, output_dir_i)
    be_array_1d = bootstrap.calc_be_1d_array(run_i)
    be_array_2d = bootstrap.calc_be_2d_array(run_i)
    run_i.flux.conf_band_1d = np.quantile(be_array_1d, CONFIDENCE_LEVEL)
    run_i.flux.conf_band_2d = np.quantile(be_array_2d, CONFIDENCE_LEVEL)
    be_array_total = bootstrap.calc_be_total_array(run_i)
    run_i.flux.conf_band_total = np.quantile(be_array_total, CONFIDENCE_LEVEL)
    # run_i.flux.conf_band_2d = 0.1
    # run_i.flux.conf_band_total = 0.1
    run_i.flux.max_energy_1d = np.max(run_i.flux.energy_1d)
    run_i.flux.max_energy_2d = np.max(run_i.flux.energy_2d)
    if free_space:
        run_i.free_space()
    paper_plots_extra.save_attributes_to_file(run_i, output_dir_i,
                                                indent=4)

def plot_ripple_runs(all_runs):
    """
    This function plots the ripple runs to produce Figure 2 of the paper.
    """

    def create_csv():
        # First we need to filter the ripple runs
        runs = []
        for run_i in all_runs:
            if run_i.log.analytic_ripple:
                runs.append(run_i)
        # Next sort the runs by rcoil and ncoil
        runs.sort(key=lambda x: (x.log.ncoil, x.log.rcoil))

        runs_metadata = []
        for i, run_i in enumerate(runs):
            output_dir_i = os.path.join(output_dir,
                                        f"rcoil_{run_i.log.rcoil}_ncoil_{run_i.log.ncoil}")
            if i==0:
                calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=run_i.log.ncoil)
            else:
                calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=run_i.log.ncoil,
                                 previous_run=runs[i-1])
            runs_metadata.append([run_i.log.ncoil,
                                  run_i.log.rcoil,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_2d])
        columns = ['ncoil', 'rcoil', 'max_energy_flux', 'total_energy_flux', 'conf_band_2d',
                   'conf_band_total', 'h_phi', 'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'ripple_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "ripple_runs")
    make_csv = False
    save_axisymmetric = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'ripple_runs.csv'))
    for run_i in all_runs:
        if run_i.log.axisymmetric:
            run_axisymmetric = run_i
    output_dir_i = os.path.join(output_dir, "axisymmetric")
    if save_axisymmetric:
        calc_energy_flux(run_axisymmetric, output_dir_i)
        pickle.dump(run_axisymmetric, open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'wb'))
    else:
        run_axisymmetric = pickle.load(open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'rb'))
    linestyles = ['-', '--', ':']
    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    for i, ncoil in enumerate([12, 16, 18]):
        rcoils = df[df.ncoil == ncoil].sort_values(by='rcoil').rcoil.values
        max_energy_flux = df[df.ncoil == ncoil].sort_values(by='rcoil').max_energy_flux.values
        total_energy_flux = \
            df[df.ncoil == ncoil].sort_values(by='rcoil').total_energy_flux.values / \
            run_axisymmetric.log.pinj * 100
        conf_band_2d = df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_2d.values
        conf_band_total = \
            df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_total.values / \
            run_axisymmetric.log.pinj * 100
        axs[0].errorbar(rcoils, max_energy_flux, yerr=conf_band_2d,
                        label=r"$N_{coil}$ = "f"{ncoil}",
                        linestyle=linestyles[i])
        axs[0].axhline(y=run_axisymmetric.flux.max_energy_2d,
                       color='k',
                       linestyle=':')
        axs[0].legend()
        axs[1].errorbar(rcoils, total_energy_flux, yerr=conf_band_total,
                        label=r"$N_{coil}$ = "f"{ncoil}",
                        linestyle=linestyles[i])
        axs[1].axhline(y=run_axisymmetric.flux.total_energy / run_axisymmetric.log.pinj * 100,
                       color='k',
                       linestyle=':')
        axs[1].legend()
    axs[0].set_ylabel(r'Max Alpha Particle Energy Flux [MW m$^{-2}$]')
    axs[1].set_ylabel(r'Alpha Power Lost [%]')
    fig.suptitle('TF Ripple Field Results')
    for i in range(2):
        axs[i].set_xlabel(r'Major radius of TF coil outer limb ($R_{outer}$) [m]')
        axs[i].set_yscale('log')
    output_path = os.path.join(output_dir, 'max_and_total_flux_vs_rcoil')
    fig.savefig(output_path + ".pdf", bbox_inches='tight')
    fig.savefig(output_path + ".png", bbox_inches='tight',
                dpi=300)
    plt.close(fig)

def plot_rmp_runs(all_runs):
    """
    This function plots the RMP runs to produce Figures 3 and 4 of the paper.
    """

    def create_csv():
        # First we need to filter the rmp runs
        runs = []
        for run_i in all_runs:
            if run_i.log.bplasma:
                if run_i.log.bplasma_n > 1:
                    runs.append(run_i)
        # Next sort the runs by rcoil and ncoil
        runs.sort(key=lambda x: x.log.rmp_phase, reverse=False)
        runs.sort(key=lambda x: x.log.rmp_response, reverse=False)
        runs.sort(key=lambda x: x.log.rmp_current, reverse=False)
        runs.sort(key=lambda x: x.log.bplasma_n, reverse=False)
        runs.sort(key=lambda x: x.log.coil_set, reverse=True)
        for run_i in runs:
            print(run_i.log.rmp_phase, run_i.log.rmp_response, run_i.log.rmp_current,
                  run_i.log.bplasma_n, run_i.log.coil_set)

        runs_metadata = []
        for run_i in runs:
            output_dir_i = os.path.join(output_dir,
                                        f"{run_i.log.coil_set}",
                                        f"{run_i.log.bplasma_n}",
                                        f"{run_i.log.rmp_current}_{run_i.log.rmp_response}",
                                        f"{run_i.log.rmp_phase}")
            # calc_energy_flux(run_i, output_dir_i,
            #                  remap_phi_n=run_i.log.bplasma_n)
            calc_energy_flux(run_i, output_dir_i,
                             remap_phi_n=None)
            runs_metadata.append([run_i.log.coil_set,
                                  run_i.log.bplasma_n,
                                  run_i.log.rmp_current,
                                  run_i.log.rmp_response,
                                  run_i.log.rmp_phase,
                                  run_i.flux.max_energy_1d,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_1d,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_1d,
                                  run_i.flux.h_theta_2d])
        columns = ['coil_set', 'bplasma_n', 'rmp_current', 'rmp_response', 'rmp_phase',
                   'max_energy_flux_1d', 'max_energy_flux_2d', 'total_energy_flux',
                   'conf_band_1d', 'conf_band_2d', 'conf_band_total', 'h_phi', 'h_theta_1d',
                   'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'rmp_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "rmp_runs")
    make_csv = False
    save_axisymmetric = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'rmp_runs.csv'))
    for run_i in all_runs:
        if run_i.log.axisymmetric:
            run_axisymmetric = run_i
    output_dir_i = os.path.join(output_dir, "axisymmetric")
    if save_axisymmetric:
        calc_energy_flux(run_axisymmetric, output_dir_i)
        pickle.dump(run_axisymmetric, open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'wb'))
    else:
        run_axisymmetric = pickle.load(open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'rb'))
    linestyles = ['--', '-', '--', '-']
    clrs = ['tab:blue', 'tab:blue', 'tab:orange', 'tab:orange']
    for coil_set in ["interior_rmp", "exterior_rmp"]:
        fig, axs = plt.subplots(3, 2)
        fig_size = fig.get_size_inches()
        fig_size[0] *= 2
        fig_size[1] *= 3
        fig.set_size_inches(fig_size)
        for i, bplasma_n in enumerate([2, 3, 4]):
            if coil_set == "interior_rmp":
                if bplasma_n == 2:
                    current0 = 30
                    optimum_phase = 265
                elif bplasma_n == 3:
                    current0 = 50
                    optimum_phase = 173
                elif bplasma_n == 4:
                    current0 = 80
                    optimum_phase = 67
            elif coil_set == "exterior_rmp":
                if bplasma_n == 2:
                    current0 = 50
                    optimum_phase = 61
                elif bplasma_n == 3:
                    current0 = 90
                    optimum_phase = 20
                elif bplasma_n == 4:
                    current0 = 150
                    optimum_phase = 321
            parameters = []
            for rmp_current in [current0, 2 * current0]:
                for rmp_response in [False, True]:
                    parameters.append((rmp_current, rmp_response))
            for j, (rmp_current, rmp_response) in enumerate(parameters):
                df_i = df.loc[(df['coil_set'] == coil_set) & (df['bplasma_n'] == bplasma_n) &
                              (df['rmp_current'] == rmp_current) &
                              (df['rmp_response'] == rmp_response)]
                rmp_phases = df_i['rmp_phase'].values
                max_energy_flux = df_i.max_energy_flux_1d.values
                total_energy_flux = df_i.total_energy_flux.values / run_axisymmetric.log.pinj * 100
                conf_band_1d = df_i.conf_band_1d.values
                conf_band_total = df_i.conf_band_total.values / run_axisymmetric.log.pinj * 100
                axs[i, 0].errorbar(rmp_phases, max_energy_flux, yerr=conf_band_1d,
                                   linestyle=linestyles[j], color=clrs[j])
                axs[i, 0].axhline(y=run_axisymmetric.flux.max_energy_1d,
                                  color='k',
                                  linestyle=':')
                axs[i, 0].axvline(x=optimum_phase, color='k', linestyle=':')
                axs[i, 1].errorbar(rmp_phases, total_energy_flux, yerr=conf_band_total,
                                   linestyle=linestyles[j], color=clrs[j])
                axisymmetric_total_energy_flux = run_axisymmetric.flux.total_energy \
                                               / run_axisymmetric.log.pinj * 100
                axs[i, 1].axhline(y=axisymmetric_total_energy_flux,
                                  color='k',
                                  linestyle=':')
                axs[i, 1].axvline(x=optimum_phase, color='k', linestyle=':')
                axs[i, 0].get_shared_y_axes().joined(axs[i, 0], axs[0, 0])
                axs[i, 1].get_shared_y_axes().joined(axs[i, 1], axs[0, 1])
            blue_line = mlines.Line2D([], [],
                                      color='tab:blue',
                                      label=f'Current: {current0} kAt')
            orange_line = mlines.Line2D([], [],
                                        color='tab:orange',
                                        label=f'Current: {2 * current0} kAt')
            dashed_line = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        label='Plasma response')
            solid_line = mlines.Line2D([], [],
                                       color='black',
                                       linestyle='--',
                                       label='Vacuum')
            leg1 = axs[i, 0].legend(handles=[orange_line, blue_line])
            leg2 = axs[i, 1].legend(handles=[dashed_line, solid_line])
            axs[i, 0].add_artist(leg1)
            axs[i, 1].add_artist(leg2)
            for j in range(2):
                axs[i, j].set_title(f'n = {bplasma_n}')
                axs[i, j].set_yscale('log')
            axs[i, 0].set_ylabel(r'Max Alpha Particle Energy Flux [MW m$^{-2}$]')
            axs[i, 1].set_ylabel(r'Alpha Power Lost [%]')
        for j in range(2):
            axs[2, j].set_xlabel(r'Phase shift ($\Delta \phi$) [degrees]')
        fig.tight_layout()
        output_path = os.path.join(output_dir,
                                    f"max_and_total_flux_vs_phase_{coil_set}")
        if coil_set == "interior_rmp":
            fig.suptitle("In-Vessel ELM Suppression Results",
                            y=1)
        elif coil_set == "exterior_rmp":
            fig.suptitle("Out-of-Vessel ELM Suppression Results",
                            y=1)
        fig.savefig(output_path + '.pdf', bbox_inches='tight')
        fig.savefig(output_path + '.png', bbox_inches='tight',
                    dpi=300)
        plt.close(fig)

def plot_rmp_distribution(all_runs):
    """
    This function plots the the energy distribution runs to produce Figure 5 of the paper.
    """
    for run_i in all_runs:
        if run_i.log.bplasma and \
           run_i.log.coil_set == "exterior_rmp" and \
           run_i.log.bplasma_n == 3 and \
           run_i.log.rmp_current == 90 and \
           run_i.log.rmp_response and \
           run_i.log.rmp_phase == 20:
            print("Run found!")
            run0 = run_i
    print(run0.log.bplasma_n, run0.log.rmp_current, run0.log.rmp_response,
            run0.log.rmp_phase)
    run0.init_wall(WALL_PATH,
                    special_nodes=SPECIAL_NODES)
    run0.init_markers(remap_phi_n=None)
    run0.init_flux(num_grid_points_1d=NUM_GRID_POINTS_1D,
                   num_grid_points_2d=NUM_GRID_POINTS_2D,
                   num_bootstraps=NUM_BOOTSTRAPS,
                   h_phi_array=H_PHI_ARRAY,
                   h_theta_1d_array=H_THETA_1D_ARRAY,
                   h_theta_2d_array=H_THETA_2D_ARRAY)
    flux.calc_optimum_bandwidth_1d(run0)
    run0.flux.energy_1d = flux.calc_energy_flux_1d(run0)
    # Make copies of run0.flux.energy_1d and run0.flux.s_theta_1d
    # with shorter length using interpolation
    new_length = 1000
    num_edges = new_length - 1
    new_s_theta = np.interp(np.linspace(0, 1, new_length),
                            np.linspace(0, 1, len(run0.flux.s_theta_1d)),
                            run0.flux.s_theta_1d)
    new_energy_flux = np.interp(np.linspace(0, 1, new_length),
                                np.linspace(0, 1, len(run0.flux.energy_1d)),
                                run0.flux.energy_1d)
    x_points, y_points = wall.get_rz_from_s_theta(new_s_theta,
                                                  run0.wall.r,
                                                  run0.wall.z)
    new_energy_flux_mid = 0.5 * (new_energy_flux[1:] + new_energy_flux[:-1])

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)

    norm = Normalize(vmin=0, vmax=np.max(run0.flux.energy_1d))
    scalar_map = ScalarMappable(norm=norm, cmap=viridis,)
    for i in range(num_edges):
        color = scalar_map.to_rgba(new_energy_flux_mid[i])
        axs[0].plot(x_points[i:i+2], y_points[i:i+2],
                    color=color)
    plt.colorbar(scalar_map, ax=axs[0], orientation='vertical')
    axs[0].set_aspect('equal')
    axs[0].set_xlim(0, 7.5)
    axs[0].set_ylim(-10, 10)
    axs[0].set_xlabel('R [m]')
    axs[0].set_ylabel('Z [m]')
    clrs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    symbols = ['+', 'x', 'o', 's']
    axs[1].plot(run0.flux.s_theta_1d, run0.flux.energy_1d, 'k')
    y_mid = 0.5 * (axs[1].get_ylim()[1] - axs[1].get_ylim()[0])
    for i, s_nod_i in enumerate(run0.wall.special_nodes):
        axs[0].plot(run0.wall.r[s_nod_i], run0.wall.z[s_nod_i],
                    color=clrs[i],
                    marker=symbols[i])
        axs[1].axvline(x=run0.wall.s_nodes[s_nod_i],
                       color=clrs[i])
        axs[1].plot(run0.wall.s_nodes[s_nod_i], y_mid,
                    color=clrs[i],
                    marker=symbols[i])

    axs[1].set_xlabel(r'Poloidal Distance along the Inner Wall ($s_\theta$) [m]')
    fig.suptitle(r'Max Alpha Particle Energy Flux [MW m$^{-2}$]',
                 x=0.6)
    output_dir = os.path.join(REPOSITORY_PATH, "plots")
    output_path = os.path.join(output_dir,
                               "energy_flux_distribution_rmp")
    fig.savefig(output_path + '.pdf', bbox_inches='tight')
    fig.savefig(output_path + '.png', bbox_inches='tight',
                dpi=300)
    plt.close(fig)

def plot_rwm_runs(all_runs):
    """
    This function plots the RWM runs to produce Figure 6 of the paper.
    """

    def create_csv():
        # First we need to filter the rwm runs
        runs = []
        for run_i in all_runs:
            if run_i.log.bplasma:
                if run_i.log.bplasma_n == 1:
                    runs.append(run_i)
        # Next sort the runs by bscale
        runs.sort(key=lambda x: x.log.rwm_bscale, reverse=False)

        runs_metadata = []
        for run_i in runs:
            output_dir_i = os.path.join(output_dir,
                                        f"bscale_{run_i.log.rwm_bscale}")
            calc_energy_flux(run_i, output_dir_i)
            runs_metadata.append([run_i.log.rwm_bscale,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_2d])
        columns = ['bscale', 'max_energy_flux', 'total_energy_flux', 'conf_band_2d',
                   'conf_band_total', 'h_phi', 'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'rwm_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "rwm_runs")
    make_csv = False
    save_axisymmetric = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'rwm_runs.csv'))
    output_dir_i = os.path.join(output_dir, "axisymmetric")
    if save_axisymmetric:
        for run_i in all_runs:
            if run_i.log.axisymmetric:
                run_axisymmetric = run_i
        calc_energy_flux(run_axisymmetric, output_dir_i)
        pickle.dump(run_axisymmetric,
                    open(os.path.join(output_dir_i, 'axisymmetric_run.pickle'), 'wb'))
    else:
        run_axisymmetric = pickle.load(open(os.path.join(output_dir_i,
                                                         'axisymmetric_run.pickle'), 'rb'))

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    bscales = df.sort_values(by='bscale').bscale.values
    max_energy_flux = df.sort_values(by='bscale').max_energy_flux.values
    total_energy_flux = df.sort_values(by='bscale').total_energy_flux.values \
                      / run_axisymmetric.log.pinj * 100
    conf_band_2d = df.sort_values(by='bscale').conf_band_2d.values
    conf_band_total = df.sort_values(by='bscale').conf_band_total.values \
                    / run_axisymmetric.log.pinj * 100
    axs[0].errorbar(bscales / 1e4, max_energy_flux, yerr=conf_band_2d)
    axs[1].errorbar(bscales / 1e4, total_energy_flux, yerr=conf_band_total)
    axs[0].set_ylabel(r'Max Alpha Particle Energy Flux [MW m$^{-2}$]')
    axs[1].set_ylabel('Alpha Power Lost [%]')
    for i in range(2):
        axs[i].set_xlabel('Magnetic field strength of RWM at sensors [T]')
        axs[i].set_xscale('log')
        axs[i].set_yscale('log')
    fig.suptitle('RWM field results')

    output_path = os.path.join(output_dir, 'max_and_total_flux_vs_bscale')
    fig.savefig(output_path + '.pdf', bbox_inches='tight')
    fig.savefig(output_path + '.png', bbox_inches='tight',
                dpi=300)
    plt.close(fig)

def spr_045_14_vs_spr_045_16():
    """
    This function plots the energy flux distribution for an SPR-045-14 run vs an
    SPR-045-16 run.
    """

    runs_path = os.path.join(REPOSITORY_PATH, "output_data",
                             "full_3d_field_spr_045_14_and_16_processed")
    runs = run.create_runs_list(runs_path)
    output_dir = os.path.join(REPOSITORY_PATH, "plots", "spr_045_14_vs_spr_045_16")
    energy_flux_dict = {}
    for run_i in runs:
        if "SPR-045-14" in run_i.log.eqdsk_fname:
            spr_string = "SPR-045-14"
        elif "SPR-045-16" in run_i.log.eqdsk_fname:
            spr_string = "SPR-045-16"
        else:
            raise ValueError("Could not determine run type")
        output_dir_i = os.path.join(output_dir, spr_string)
        os.makedirs(output_dir_i, exist_ok=True)
        calc_energy_flux(run_i, output_dir_i,
                         free_space=False)
        energy_flux_dict[spr_string] = run_i.flux.energy_1d

    fig, ax = plt.subplots()
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    for spr_string in ['SPR-045-14', 'SPR-045-16']:
        energy_flux = energy_flux_dict[spr_string]
        if spr_string == 'SPR-045-14':
            linestyle = '-'
            color = 'tab:blue'
            alpha = 1
            label = 'Dirtier plasma scenario'
        elif spr_string == 'SPR-045-16':
            linestyle = '--'
            color = 'tab:orange'
            alpha = 1
            label = 'Baseline plasma scenario'
        ax.plot(runs[0].flux.s_theta_1d, energy_flux,
                linestyle=linestyle,
                color=color,
                alpha=alpha,
                label=label)
    ax.legend()
    clrs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    symbols = ['+', 'x', 'o', 's']
    y_mid = 0.5 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    for i, s_nod_i in enumerate(runs[0].wall.special_nodes):
        ax.axvline(x=runs[0].wall.s_nodes[s_nod_i],
                   color=clrs[i])
        ax.plot(runs[0].wall.s_nodes[s_nod_i], y_mid,
                color=clrs[i],
                marker=symbols[i])
    ax.annotate(r'Collisional losses',
                xy=(14, 0.06),
                xytext=(19, 0.06),
                arrowprops=dict(facecolor='black', arrowstyle='->',
                lw=1.5))
    ax.annotate(r'Prompt losses',
                xy=(40, 0.025),
                xytext=(20, 0.025),
                arrowprops=dict(facecolor='black', arrowstyle='->',
                lw=1.5))
    ax.set_xlabel(r'$s_\theta$ [m]')
    ax.set_ylabel(r'Max Alpha Particle Energy Flux [MW m$^{-2}$]')
    fig.savefig(output_dir + '/energy_flux_spr_045_14_vs_spr_045_16.png',
                bbox_inches='tight', dpi=300)
    fig.savefig(output_dir + '/energy_flux_spr_045_14_vs_spr_045_16.pdf',
                bbox_inches='tight')
    plt.close('all')

    # Now make SPR-045-16 plots

    spr_string = "SPR-045-16"
    energy_flux = energy_flux_dict[spr_string]
    run0 = runs[0]

    new_length = 1000
    num_edges = new_length - 1
    new_s_theta = np.interp(np.linspace(0, 1, new_length),
                            np.linspace(0, 1, len(run0.flux.s_theta_1d)),
                            run0.flux.s_theta_1d)
    new_energy_flux = np.interp(np.linspace(0, 1, new_length),
                                np.linspace(0, 1, len(run0.flux.energy_1d)),
                                run0.flux.energy_1d)
    x_points, y_points = wall.get_rz_from_s_theta(new_s_theta,
                                                  run0.wall.r,
                                                  run0.wall.z)
    new_energy_flux_mid = 0.5 * (new_energy_flux[1:] + new_energy_flux[:-1])

    fig = plt.figure()
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig_size[1] *= 3
    fig.set_size_inches(fig_size)
    gs = gridspec.GridSpec(3, 2, width_ratios=[0.05, 1])
    ax1 = fig.add_subplot(gs[0:2, 1])
    ax2 = fig.add_subplot(gs[2, :])
    cmap = LinearSegmentedColormap.from_list(
        'custom_red',
        [(1, 1, 1), (1, 0.5, 1)],
        N=256
    )
    norm = Normalize(vmin=0, vmax=0.03)
    scalar_map = ScalarMappable(norm=norm,
                                cmap=cmap)
    ax1.plot(run0.wall.r, run0.wall.z, 'k',
            linewidth=1)
    ax1.annotate(r'$s_\theta=0$',
                xy=(run0.wall.r[0], run0.wall.z[0]),
                xytext=(8, run0.wall.z[0]),
                arrowprops=dict(facecolor='black', arrowstyle='->',
                                linewidth=1.5))
    ax1.annotate(r'Prompt losses',
                xy=(2.6, 6.3),
                xytext=(8, 6.3),
                arrowprops=dict(facecolor='black', arrowstyle='->',
                                lw=1.5))
    ax1.annotate(r'Collisional losses',
                xy=(6, -8.5),
                xytext=(8, -8.5),
                arrowprops=dict(facecolor='black', arrowstyle='->',
                                lw=1.5))
    for i in range(num_edges):
        color = scalar_map.to_rgba(new_energy_flux_mid[i])
        alpha = new_energy_flux_mid[i]/np.max(new_energy_flux_mid)
        alpha = np.max([0, alpha])
        ax1.plot(x_points[i:i+2], y_points[i:i+2],
                alpha=alpha,
                linewidth=6,
                color=color)
    plt.colorbar(scalar_map, ax=ax1, orientation='vertical',
                 location='left',
                 pad=0.3,
                 label=r'Max Alpha Particle Energy Flux [MW m$^{-2}$]')

    for i, s_nod_i in enumerate(run0.wall.special_nodes):
        ax1.plot(run0.wall.r[s_nod_i], run0.wall.z[s_nod_i],
                color=clrs[i],
                marker=symbols[i])
    ax1.set_aspect('equal')
    ax1.set_xlim(0, 7.5)
    ax1.set_ylim(-10, 10)
    ax1.set_xlabel('R [m]')
    ax1.set_ylabel('Z [m]')

    # Regular line plot of energy flux
    ax2.plot(run0.flux.s_theta_1d, energy_flux, 'k')
    ax2.set_xlabel(r'$s_\theta$ [m]')
    ax2.set_ylabel(r'Max Alpha Particle Energy Flux [MW m$^{-2}$]')
    y_mid = 0.5 * (ax2.get_ylim()[1] - ax2.get_ylim()[0])
    for i, s_nod_i in enumerate(runs[0].wall.special_nodes):
        ax2.axvline(x=runs[0].wall.s_nodes[s_nod_i],
                   color=clrs[i])
        ax2.plot(runs[0].wall.s_nodes[s_nod_i], y_mid,
                color=clrs[i],
                marker=symbols[i])
    ax2.annotate(r'Collisional losses',
                 xy=(14, 0.025),
                 xytext=(19, 0.015),
                 arrowprops=dict(facecolor='black', arrowstyle='->',
                 lw=1.5))
    ax2.annotate(r'Prompt losses',
                 xy=(41, 0.025),
                 xytext=(25, 0.025),
                 arrowprops=dict(facecolor='black', arrowstyle='->',
                 lw=1.5))
    output_path = os.path.join(output_dir, "energy_flux_full_3d")
    fig.suptitle('TF ripple + RWM + Out-of-Vessel ELM Suppression Results',
                 y=0.93)
    fig.savefig(output_path + '.pdf', bbox_inches='tight')
    fig.savefig(output_path + '.png', bbox_inches='tight', dpi=300)
    plt.close('all')

def spr_068_spr_045_ripple_scan():
    """
    This function plots the ripple runs for SPR-068-7 and SPR-045-14 and SPR-045-16.
    """
    def create_csv_ripple():
        # First we need to filter the ripple runs
        runs = []
        for run_i in all_runs:
            print(run_i.dir_path)
            if run_i.log.analytic_ripple:
                runs.append(run_i)
        # Next sort the runs by rcoil and ncoil
        runs.sort(key=lambda x: (x.log.eqdsk_fname,x.log.ncoil, x.log.rcoil))
        # runs.sort(key=lambda x: (x.log.ncoil, x.log.rcoil))
        # runs.sort(key=lambda x: x.log.eqdsk_fname, reverse=True)

        runs_metadata = []
        for i, run_i in enumerate(runs):
            if "SPR-045-14" in run_i.log.eqdsk_fname:
                spr_string = "SPR-045-14"
                wall_path = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
                special_nodes = (16, 55, 138, 178)
            elif "SPR-045-16" in run_i.log.eqdsk_fname:
                spr_string = "SPR-045-16"
                special_nodes = (16, 55, 138, 178)
                wall_path = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
            elif "SPR-068-7" in run_i.log.eqdsk_fname:
                spr_string = "SPR-068-7"
                special_nodes = (4, 27, 40, 63)
                wall_path = os.path.join(REPOSITORY_PATH, "input_data", "SPR-068_wall.dat")
            else:
                raise ValueError("Could not determine run type")
            output_dir_i = os.path.join(output_dir, spr_string,
                                        f"rcoil_{run_i.log.rcoil}_ncoil_{run_i.log.ncoil}")
            if i==0:
                calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=run_i.log.ncoil,
                                 wall_path=wall_path,
                                 special_nodes=special_nodes)
            else:
                calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=run_i.log.ncoil,
                                 previous_run=runs[i-1],
                                 wall_path=wall_path,
                                 special_nodes=special_nodes)
            runs_metadata.append([run_i.log.ncoil,
                                  run_i.log.rcoil,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_2d,
                                  spr_string,
                                  run_i.log.pinj])
        columns = ['ncoil', 'rcoil', 'max_energy_flux', 'total_energy_flux', 'conf_band_2d',
                   'conf_band_total', 'h_phi', 'h_theta_2d', 'spr_string', 'pinj']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'spr_068_spr_045_ripple_scan.csv'))

    def create_csv_axisymmetric():
        # First we need to filter the axisymmetric runs
        runs = []
        for run_i in all_runs:
            if run_i.log.axisymmetric:
                runs.append(run_i)
        runs.sort(key=lambda x: (x.log.eqdsk_fname), reverse=True)

        runs_metadata = []
        for i, run_i in enumerate(runs):
            if "SPR-045-14" in run_i.log.eqdsk_fname:
                spr_string = "SPR-045-14"
                wall_path = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
                special_nodes = (16, 55, 138, 178)
            elif "SPR-045-16" in run_i.log.eqdsk_fname:
                spr_string = "SPR-045-16"
                wall_path = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
                special_nodes = (16, 55, 138, 178)
            elif "SPR-068-7" in run_i.log.eqdsk_fname:
                spr_string = "SPR-068-7"
                special_nodes = (4, 27, 40, 63)
                wall_path = os.path.join(REPOSITORY_PATH, "input_data", "SPR-068_wall.dat")
            else:
                raise ValueError("Could not determine run type")
            output_dir_i = os.path.join(output_dir,
                                        spr_string)
            if i==0:
                calc_energy_flux(run_i, output_dir_i,
                                 wall_path=wall_path,
                                 special_nodes=special_nodes)
            else:
                calc_energy_flux(run_i, output_dir_i,
                                 previous_run=runs[i-1],
                                 wall_path=wall_path,
                                 special_nodes=special_nodes)
            runs_metadata.append([run_i.log.ncoil,
                                  run_i.log.rcoil,
                                  run_i.flux.max_energy_1d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_1d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_1d,
                                  spr_string,
                                  run_i.log.pinj])
        columns = ['ncoil', 'rcoil', 'max_energy_flux', 'total_energy_flux', 'conf_band_2d',
                   'conf_band_total', 'h_phi', 'h_theta_1d', 'spr_string', 'pinj']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'spr_068_spr_045_axisymmetric_scan.csv'))

    runs_path = os.path.join(REPOSITORY_PATH, "output_data",
                             "spr_068_spr_045_ripple_scan")
    all_runs = run.create_runs_list(runs_path)
    make_csv = False
    save_axisymmetric = False
    output_dir = os.path.join(REPOSITORY_PATH, "plots", "spr_068_spr_045_ripple_scan")
    if make_csv:
        create_csv_ripple()
    df_r = pd.read_csv(os.path.join(output_dir, 'spr_068_spr_045_ripple_scan.csv'))
    output_dir = os.path.join(REPOSITORY_PATH, "plots", "spr_068_spr_045_axisymmetric_scan")
    if save_axisymmetric:
        create_csv_axisymmetric()
    df_a = pd.read_csv(os.path.join(output_dir, 'spr_068_spr_045_axisymmetric_scan.csv'))
    linestyles = ['-', '--', ':']
    colours = ['tab:blue', 'tab:orange', 'tab:green']
    spr_strings = ['SPR-045-14', 'SPR-045-16', 'SPR-068-7']
    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    for i, spr_string in enumerate(spr_strings):
        rcoils = df_r[df_r.spr_string == spr_string].sort_values(by='rcoil').rcoil.values
        max_energy_flux = \
            df_r[df_r.spr_string == spr_string].sort_values(by='rcoil').max_energy_flux.values
        total_energy_flux = \
            df_r[df_r.spr_string == spr_string].sort_values(by='rcoil').total_energy_flux.values /\
            df_r[df_r.spr_string == spr_string].sort_values(by='rcoil').pinj.values[0] * 100
        conf_band_2d = \
            df_r[df_r.spr_string == spr_string].sort_values(by='rcoil').conf_band_2d.values
        conf_band_total = \
            df_r[df_r.spr_string == spr_string].sort_values(by='rcoil').conf_band_total.values / \
            df_r[df_r.spr_string == spr_string].sort_values(by='rcoil').pinj.values[0] * 100
        axs[0].errorbar(rcoils, max_energy_flux, yerr=conf_band_2d,
                        label=spr_string,
                        linestyle=linestyles[i],
                        color=colours[i])
        axs[0].axhline(y=df_a[df_a.spr_string == spr_string].max_energy_flux.values[0],
                       color='k',
                       linestyle=linestyles[i])
        axs[0].legend()
        axs[1].errorbar(rcoils, total_energy_flux, yerr=conf_band_total,
                        label=spr_string,
                        linestyle=linestyles[i],
                        color=colours[i])
        axs[1].axhline(y=df_a[df_a.spr_string == spr_string].total_energy_flux.values[0] /
                         df_a[df_a.spr_string == spr_string].pinj.values[0] * 100,
                       color='k',
                       linestyle=linestyles[i])
        axs[1].legend()
    axs[0].set_ylabel(r'Max Alpha Particle Energy Flux [MW m$^{-2}$]')
    axs[1].set_ylabel(r'Alpha Power Lost [%]')
    fig.suptitle('TF Ripple Field Results')
    for i in range(2):
        axs[i].set_xlabel(r'Major radius of TF coil outer limb ($R_{outer}$) [m]')
        axs[i].set_yscale('log')
        axs[i].grid()
        axs[i].minorticks_on()
    output_path = os.path.join(output_dir, 'max_and_total_flux_vs_rcoil')
    fig.savefig(output_path + ".pdf", bbox_inches='tight')
    fig.savefig(output_path + ".png", bbox_inches='tight',
                dpi=300)
    plt.close(fig)

def plot_background_plasma_curves():
    """
    This function plots the background plasma curves to produce the 
    density temperature, safety factor and reaction rate plots in the paper.
    """
    spr_string = 'SPR-068-RV'
    cdf_filename = os.path.join(REPOSITORY_PATH, "input_data", f"profiles_{spr_string}.CDF")
    psin, ti, te, ne, nd, nt = prepare_profiles.read_cdf_file(cdf_filename)
    gfile_path = os.path.join(REPOSITORY_PATH, "input_data", f'{spr_string}.eqdsk')
    gfile = prepare_profiles.get_gfile(gfile_path)
    num_impurities = prepare_profiles.calculate_number_of_impurities(cdf_filename)
    _, _, nim = prepare_profiles.get_impurity_data(cdf_filename, num_impurities)
    fd, ft, fim = prepare_profiles.calculate_ion_fractions(ne, nd, nt, nim)
    reactivity = dt_fusion.reactivity(ti*1e-3) # convert from eV to keV
    reactivity *= 1e-6  # convert from cm^3/s to m^3/s
    reaction_rate = reactivity * ne * fd * ne * ft

    fig, axs = plt.subplots(2, 2)
    fig_size = fig.get_size_inches()
    fig_size *= 2
    fig.set_size_inches(fig_size)
    axs[0, 0].plot(psin, ti * 1e-3, label=r'$T_i$')
    axs[0, 0].plot(psin, te * 1e-3, label=r'$T_e$')
    axs[0, 0].legend()
    axs[0, 0].set_ylabel('Temperature [keV]')
    axs[0, 1].plot(gfile.psin, gfile.q)
    axs[0, 1].set_ylabel('Safety Factor')
    axs[1, 0].plot(psin, ne, label=r'$n_e$')
    axs[1, 0].plot(psin, fd * ne, label=r'$n_D$')
    axs[1, 0].plot(psin, ft * ne, label=r'$n_T$')
    axs[1, 0].plot(psin, fim[0] * ne, label=r'$n_{Xe}$')
    axs[1, 0].plot(psin, fim[1] * ne, label=r'$n_{He}$')
    axs[1, 0].set_ylabel('Density of particle species [m$^{-3}$]')
    axs[1, 0].set_yscale('log')
    axs[1, 0].legend()
    axs[1, 1].plot(psin, reaction_rate)
    axs[1, 1].set_ylabel(r'DT Fusion Reaction Rate [m$^{-3}$ s$^{-1}$]')
    for i in range(2):
        for j in range(2):
            axs[i, j].set_xlabel(r'$\psi_N$')
    output_dir = os.path.join(REPOSITORY_PATH, "plots")
    output_path = os.path.join(output_dir, 'background_plasma_curves_' + spr_string)
    fig.savefig(output_path + '.pdf',
                bbox_inches='tight')
    fig.savefig(output_path + '.png',
                bbox_inches='tight',
                dpi=300)
    plt.close('all')

    # Plot comparing SPR-045-14 and SPR-045-16
    cdf_filename_14 = os.path.join(REPOSITORY_PATH, "input_data", "profiles_SPR-045-14.CDF")
    psin_14, ti_14, te_14, ne_14, nd_14, nt_14 = prepare_profiles.read_cdf_file(cdf_filename_14)
    gfile_path_14 = os.path.join(REPOSITORY_PATH, "input_data", 'SPR-045-14.eqdsk')
    gfile_14 = prepare_profiles.get_gfile(gfile_path_14)
    num_impurities_14 = prepare_profiles.calculate_number_of_impurities(cdf_filename_14)
    _, _, nim_14 = prepare_profiles.get_impurity_data(cdf_filename_14, num_impurities_14)
    fd_14, ft_14, fim_14 = prepare_profiles.calculate_ion_fractions(ne_14, nd_14, nt_14, nim_14)
    reactivity_14 = dt_fusion.reactivity(ti_14*1e-3) # convert from eV to keV
    reactivity_14 *= 1e-6  # convert from cm^3/s to m^3/s
    reaction_rate_14 = reactivity_14 * ne_14 * fd_14 * ne_14 * ft_14
    fig, axs = plt.subplots(2, 2)
    fig_size = fig.get_size_inches()
    fig_size *= 2
    fig.set_size_inches(fig_size)
    axs[0, 0].plot(psin_14, ti_14 * 1e-3, label=r'$T_i$')
    axs[0, 0].plot(psin_14, te_14 * 1e-3, label=r'$T_e$')
    axs[0, 0].legend()
    axs[0, 0].set_ylabel('Temperature [keV]')
    axs[0, 1].plot(gfile_14.psin, gfile_14.q, label='Dirtier plasma scenario')
    axs[0, 1].set_ylabel('Safety Factor')
    axs[1, 0].plot(psin_14, ne_14, label=r'$n_e$')
    axs[1, 0].plot(psin_14, fd_14 * ne_14, label=r'$n_D$')
    axs[1, 0].plot(psin_14, ft_14 * ne_14, label=r'$n_T$')
    axs[1, 0].plot(psin_14, fim_14[0] * ne_14, label=r'$n_{Xe}$')
    axs[1, 0].plot(psin_14, fim_14[1] * ne_14, label=r'$n_{He}$')
    axs[1, 0].plot(psin_14, fim_14[2] * ne_14, label=r'$n_{Ar}$')
    axs[1, 0].set_ylabel('Density of particle species [m$^{-3}$]')
    axs[1, 0].set_yscale('log')
    axs[1, 0].legend()
    axs[1, 1].plot(psin_14, reaction_rate_14)
    axs[1, 1].set_ylabel(r'DT Fusion Reaction Rate [m$^{-3}$ s$^{-1}$]')
    for i in range(2):
        for j in range(2):
            axs[i, j].set_prop_cycle(None)
    axs[0, 0].plot(psin, ti * 1e-3, alpha=0.5, linestyle='--')
    axs[0, 0].plot(psin, te * 1e-3, alpha=0.5, linestyle='--')
    axs[0, 1].plot(gfile.psin, gfile.q, alpha=0.5, linestyle='--', label='Baseline plasma scenario')
    axs[0, 1].legend()
    axs[1, 0].plot(psin, ne, label=r'$n_e$', alpha=0.5, linestyle='--')
    axs[1, 0].plot(psin, fd * ne, label=r'$n_D$', alpha=0.5, linestyle='--')
    axs[1, 0].plot(psin, ft * ne, label=r'$n_T$', alpha=0.5, linestyle='--')
    axs[1, 0].plot(psin, fim[0] * ne, label=r'$n_{Xe}$', alpha=0.5, linestyle='--')
    axs[1, 0].plot(psin, fim[1] * ne, label=r'$n_{He}$', alpha=0.5, linestyle='--')
    axs[1, 1].plot(psin, reaction_rate, alpha=0.5, linestyle='--')
    for i in range(2):
        for j in range(2):
            axs[i, j].set_xlabel(r'$\psi_N$')
    output_dir = os.path.join(REPOSITORY_PATH, "plots")
    output_path = os.path.join(output_dir, 'background_plasma_curves_14_vs_16')
    fig.savefig(output_path + '.pdf',
                bbox_inches='tight')
    fig.savefig(output_path + '.png',
                bbox_inches='tight',
                dpi=300)
    plt.close('all')

def plot_hotpost_distributions():

    class Circle:
        """Class to represent a circle with a center and radius."""
        def __init__(self, centre_r, centre_z, radius):
            self.centre_r = centre_r
            self.centre_z = centre_z
            self.radius = radius

        def generate_circle(self, num_points=100):
            """Generates the r and z coordinates of the circle for plotting."""
            circle_s = np.linspace(0, 2 * np.pi, num_points)
            circle_r = self.centre_r + self.radius * np.cos(circle_s)
            circle_z = self.centre_z + self.radius * np.sin(circle_s)
            return circle_r, circle_z

    # Define the hotspots using the Circle class
    hotspots = {
        "top": Circle(centre_r=2.5, centre_z=6.5, radius=1),
        "bottom": Circle(centre_r=6.0, centre_z=-8.0, radius=1),
        "left": Circle(centre_r=1.8, centre_z=2.5, radius=1)
    }

    # Paths and run configuration
    runs_path = os.path.join(REPOSITORY_PATH, "output_data", "full_3d_field_spr_045_14_and_16_processed")
    plot_dir = os.path.join(REPOSITORY_PATH, "plots", "hotspot_distributions")
    os.makedirs(plot_dir, exist_ok=True)
    full_3d_run_dir = os.path.join(runs_path, "gpu-q-52")
    full_3d_run_tag = "25-07-2024_12-17-50.555"
    full_3d_run = run.Run(full_3d_run_dir, full_3d_run_tag)
    full_3d_run.init_gfile(GFILE_PATH)
    full_3d_run.init_markers()

    # Plot markers and the hotspot circles
    # fig, ax = plt.subplots()
    # ax.scatter(full_3d_run.markers.stopped.r, full_3d_run.markers.stopped.z, s=0.1)
    # # Plot each hotspot circle using a loop
    # for label, hotspot in hotspots.items():
    #     circle_r, circle_z = hotspot.generate_circle()
    #     ax.plot(circle_r, circle_z, 'k--')
    # ax.set_aspect('equal')
    # ax.set_xlabel('R [m]')
    # ax.set_ylabel('Z [m]')

    # Find the indices of the markers inside the hotspot circles
    hotspot_indices = {}
    for label, hotspot in hotspots.items():
        hotspot_indices[label] = np.where(
            (full_3d_run.markers.stopped.r < hotspot.centre_r + hotspot.radius) &
            (full_3d_run.markers.stopped.r > hotspot.centre_r - hotspot.radius) &
            (full_3d_run.markers.stopped.z < hotspot.centre_z + hotspot.radius) &
            (full_3d_run.markers.stopped.z > hotspot.centre_z - hotspot.radius)
        )[0]
    
    # Plot markers inside each hotspot
    # fig, ax = plt.subplots()
    # for label, indices in hotspot_indices.items():
    #     ax.scatter(full_3d_run.markers.stopped.r[indices], full_3d_run.markers.stopped.z[indices], s=0.1)
    # ax.set_aspect('equal')
    # ax.set_xlabel('R [m]')
    # ax.set_ylabel('Z [m]')

    # Make a histogram of the energy distribution for each hotspot
    # Convert energy from J to MeV
    fig, ax = plt.subplots()
    for label, indices in hotspot_indices.items():
        # energy = energy * mass of helium * charge of electron
        ax.hist(full_3d_run.markers.stopped.energy[indices] / 1e6 * 6.6464731e-27 / 1.602e-19,
                weights=full_3d_run.markers.stopped.weight[indices],
                histtype='step',
                label = label + ' hotspot',
                bins=100)
    ax.legend()
    ax.set_xlabel('Energy [MeV]')
    ax.set_ylabel('Weighted Number of Markers')
    ax.set_title('Energy distribution of SPR-045-16 Hotspots')
    fig.savefig(os.path.join(plot_dir, 'hotspot_energy.png'),
                bbox_inches='tight',
                dpi=300)

    # Make a histogram of the initial v_|| / v for each hotspot
    fig, ax = plt.subplots()
    for label, indices in hotspot_indices.items():
        ax.hist(full_3d_run.markers.stopped.v_parallel0[indices] /
                np.sqrt(2 * full_3d_run.markers.stopped.energy0[indices]),
                weights=full_3d_run.markers.stopped.weight[indices],
                histtype='step',
                label = label + ' hotspot',
                bins=100)
    ax.legend()
    ax.set_xlabel(r'$v_{||}$ / |v|')
    ax.set_title('Initial pitch angle distribution of SPR-045-16 Hotspots')
    ax.set_ylabel('Weighted Number of Markers')
    fig.savefig(os.path.join(plot_dir, 'hotspot_v_parallel.png'),
                bbox_inches='tight',
                dpi=300)

    # Make a histogram of initial arctan(Z, R - gfile.rmaxis) for each hotspot
    # fig, ax = plt.subplots()
    # for label, indices in hotspot_indices.items():
    #     ax.hist(np.arctan2(full_3d_run.markers.stopped.z0[indices],
    #                        full_3d_run.markers.stopped.r0[indices] - full_3d_run.gfile.rmaxis),
    #             weights=full_3d_run.markers.stopped.weight[indices],
    #             range = (-np.pi, np.pi),
    #             histtype='step',
    #             bins=100)
    # ax.set_xlabel('arctan(Z, R - gfile.rmaxis)')
    # ax.set_ylabel('Weighted Number of Markers')
    # plt.show()

def plot_magnetic_flux_surfaces():
    """
    Plot the magnetic flux surfaces using the eqdsk files
    for SPR-045-14, SPR-045-16 and SPR-068-7.
    """

    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gfile_paths = []
    gfile_paths.append(os.path.join(repo_path, "input_data", "SPR-045-14.eqdsk"))
    gfile_paths.append(os.path.join(repo_path, "input_data", "SPR-045-16.eqdsk"))
    gfile_paths.append(os.path.join(repo_path, "input_data", "SPR-068-7.eqdsk"))
    gfile_labels = ["SPR-045-14 LCFS", "SPR-045-16 LCFS", "SPR-068-7 LCFS"]
    gfiles = []
    for gfile_path in gfile_paths:
        gfile = my_gfile_reader.getGfile(gfile_path)
        gfiles.append(gfile)
    
    wall_paths = []
    wall_paths.append(os.path.join(repo_path, "input_data", "SPP-001_wall.dat"))
    wall_paths.append(os.path.join(repo_path, "input_data", "SPR-068_wall.dat"))
    wall_labels = ["SPR-045 Wall", "SPR-068 Wall"]

    fig, ax = plt.subplots()
    fig_size = fig.get_size_inches()
    fig_size *= 2
    fig.set_size_inches(fig_size)
    for i, gfile_path in enumerate(gfile_paths):
        gfile = my_gfile_reader.getGfile(gfile_path)
        ax.plot(gfile.R_bnd, gfile.Z_bnd,
                linestyle='--',
                label=gfile_labels[i])
    for i, wall_path in enumerate(wall_paths):
        rz_coords = np.loadtxt(wall_path)
        r_coords = rz_coords[:, 0]
        z_coords = rz_coords[:, 1]
        print(np.shape(rz_coords))
        ax.plot(r_coords, z_coords,
                label=wall_labels[i])
    xlim = ax.get_xlim()
    xlim = [0, xlim[1]]
    ax.set_xlim(xlim)
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    ax.grid(True)
    ax.legend()
    ax.set_aspect('equal')
    output_dir = os.path.join(REPOSITORY_PATH, "plots")
    output_path = os.path.join(output_dir, 'lcfs_and_walls.png')
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    start_time = time.time()
    RUNS = run.create_runs_list(RUNS_DIRECTORY)
    # paper_plots_3d.coil_plot_3d(gfile_path=GFILE_PATH)
    # plot_ripple_runs(RUNS)
    # plot_rmp_runs(RUNS)
    # plot_rmp_distribution(RUNS)
    # plot_rwm_runs(RUNS)
    # plot_background_plasma_curves()
    # ripple_check.plot_ripple_field()
    # paper_plots_extra.tf_coil_inner_limb_scan()
    # spr_045_14_vs_spr_045_16()
    # plot_hotpost_distributions()
    # plot_magnetic_flux_surfaces()
    spr_068_spr_045_ripple_scan()

    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.2e} seconds")
