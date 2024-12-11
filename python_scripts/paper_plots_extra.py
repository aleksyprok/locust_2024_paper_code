"""
This file contains plots which are not shown in the IAEA FEC 2024 paper,
but are generated at the same time as the paper plots.
"""
import os
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import pandas as pd
import python_scripts.run as python_scripts_run
from python_scripts import my_gfile_reader, paper_plots

def plot_amise_1d(run, output_dir):
    """
    Plot the Asymptotic Mean Integrated Square error array.
    """
    fig, ax = plt.subplots()
    ax.plot(run.flux.h_theta_1d_array, run.flux.amise_1d)
    h_theta_index = np.where(run.flux.h_theta_1d_array == run.flux.h_theta_1d)[0][0]
    ax.plot(run.flux.h_theta_1d, run.flux.amise_1d[h_theta_index],
            'r.', markersize=10)
    ax.set_title('Asymptotic mean integrated square error')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'$h_\theta$ [m]')
    fig.savefig(output_dir + '/amise_1d_array.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

def plot_amise_2d(run, output_dir):
    """
    Plot the Asymptotic Mean Integrated Square error array.
    """
    hx_mesh, hy_mesh = np.meshgrid(run.flux.h_phi_array,
                                   run.flux.h_theta_2d_array)
    fig, ax = plt.subplots()
    im = ax.pcolormesh(hx_mesh, hy_mesh, run.flux.amise_2d,
                       norm=LogNorm())
    ax.plot(run.flux.h_phi, run.flux.h_theta_2d,
            'r.', markersize=10)
    ax.set_title('Asymptotic mean integrated square error')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.colorbar(im, ax=ax)
    ax.set_xlabel(r'$h_\phi$ [m]')
    ax.set_ylabel(r'$h_\theta$ [m]')
    fig.savefig(output_dir + '/amise_2d_array.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)

    h_theta_index = np.where(run.flux.h_theta_2d_array == run.flux.h_theta_2d)[0][0]
    h_phi_index = np.where(run.flux.h_phi_array == run.flux.h_phi)[0][0]
    axs[0].plot(run.flux.h_phi_array, run.flux.amise_2d[h_theta_index, :])
    axs[0].plot(run.flux.h_phi, run.flux.amise_2d[h_theta_index, h_phi_index],
                'r.', markersize=10)
    axs[0].set_xlabel(r'$h_\phi$ [m]')
    axs[0].set_ylabel(r'$AMISE$')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_title(r'$h_\theta$ = ' f'{run.flux.h_theta_2d} m')

    axs[1].plot(run.flux.h_theta_2d_array, run.flux.amise_2d[:, h_phi_index])
    axs[1].plot(run.flux.h_theta_2d, run.flux.amise_2d[h_theta_index, h_phi_index],
                'r.', markersize=10)
    axs[1].set_xlabel(r'$h_\theta$ [m]')
    axs[1].set_ylabel(r'$AMISE$')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_title(r'$h_\phi$ = ' f'{run.flux.h_phi} m')
    fig.savefig(output_dir + '/amise_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

def plot_energy_flux_1d(run, output_dir):
    """
    Plot the energy flux array.
    """
    fig, ax = plt.subplots()
    ax.plot(run.flux.s_theta_1d, run.flux.energy_1d)
    ax.set_title('Energy flux [MW/m^2]\n'
                 f'h_phi = {run.flux.h_phi:.2e}, h_theta_2d = {run.flux.h_theta_2d:.2e}')
    for s_nod_i in run.wall.special_nodes:
        ax.axvline(x=run.wall.s_nodes[s_nod_i], color='k')
    ax.set_xlabel(r's_theta [m]')
    fig.savefig(output_dir + '/energy_flux_1d_array.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

def plot_energy_flux_2d(run, output_dir):
    """
    Plot the energy flux array.
    """
    fig, ax = plt.subplots()
    ax.imshow(run.flux.energy_2d.T, origin='lower',
               extent=[run.wall.s_phi_min, run.wall.s_phi_max,
                       run.wall.s_theta_min, run.wall.s_theta_max])
    ax.set_xlabel('s_phi [m]')
    ax.set_ylabel('s_theta [m]')
    ax.set_title('Energy flux [MW/m^2]\n'
                 f'h_phi = {run.flux.h_phi:.2e}, h_theta_2d = {run.flux.h_theta_2d:.2e}')
    fig.savefig(output_dir + '/energy_flux_2d.png',
                bbox_inches='tight', dpi=300)

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    axs[0].plot(run.flux.s_theta_2d, np.max(run.flux.energy_2d, axis=0))
    axs[1].plot(run.flux.s_phi, np.max(run.flux.energy_2d, axis=1))
    for s_nod_i in run.wall.special_nodes:
        axs[0].axvline(x=run.wall.s_nodes[s_nod_i], color='k')
    for i in range(2):
        axs[i].set_title('Max Energy flux [MW/m^2]\n'
                            f'h_phi = {run.flux.h_phi:.2e}, h_theta_2d = '
                            f'{run.flux.h_theta_2d:.2e}')
    axs[0].set_xlabel('s_theta [m]')
    axs[1].set_xlabel('s_phi [m]')
    fig.savefig(output_dir + '/max_energy_flux_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

def save_attributes_to_file(obj, output_dir, indent=0):
    """
    Recursively saves attributes and their values for a given class instance to a file.
    If an attribute is a class instance, it does the same for that class.

    Parameters:
    obj (object): The class instance.
    output_dir (str): Directory to the output file.
    indent (int): The indentation level for pretty printing.
    """
    file_path = os.path.join(output_dir, 'run_i_attributes.txt')
    with open(file_path, "w", encoding="utf-8") as file:
        def write_attributes(obj, indent=0):
            indent_str = '    ' * indent

            if not hasattr(obj, "__dict__"):
                file.write(f"{indent_str}{obj} is not a class instance.\n")
                return

            for attr, value in vars(obj).items():
                if hasattr(value, "__dict__"):
                    file.write(f"{indent_str}{attr}:\n")
                    write_attributes(value, indent + 1)
                else:
                    file.write(f"{indent_str}{attr}: {value}\n")

        write_attributes(obj, indent=indent)

def save_percent_power_lost_to_file(output_dir):
    """
    Saves the percent power lost to a file.
    """

    def get_percent_power_lost(conditions_list):
        for run_i in runs:
            conditions_met = True
            for condition_i in conditions_list:
                key = condition_i[0]
                value = condition_i[1]
                if not run_i.log.__dict__[key] == value:
                    conditions_met = False
            if conditions_met:
                return run_i.log.total_stopped_power * 100 / run_i.log.pinj
    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    runs_directory = os.path.join(repository_path, "output_data", "FEC_2024")
    runs = python_scripts_run.create_runs_list(runs_directory)

    master_conditions_list = [
        [('axisymmetric', True)],
        [('coil_set', 'exterior_rmp'), ('rmp_response', True),
         ('rmp_current', 90), ('rmp_phase', 20)],
        [('coil_set', 'exterior_rmp'), ('rmp_response', True),
         ('rmp_current', 180), ('rmp_phase', 20)],
        [('analytic_ripple', True), ('rcoil', 7), ('ncoil', 12)],
        [('analytic_ripple', True), ('rcoil', 8), ('ncoil', 12)],
        [('analytic_ripple', True), ('rcoil', 8.5), ('ncoil', 12)],
        [('analytic_ripple', True), ('rcoil', 7), ('ncoil', 16)],
        [('analytic_ripple', True), ('rcoil', 7), ('ncoil', 18)],
    ]

    file_path = os.path.join(output_dir, 'percent_power_lost.txt')
    with open(file_path, "w", encoding="utf-8") as file:
        for conditions_list in master_conditions_list:
            percent_power_lost = get_percent_power_lost(conditions_list)
            for condition_i in conditions_list:
                key = condition_i[0]
                value = condition_i[1]
                file.write(f"{key} = {value}, ")
            file.write(f"percent_power_lost = {percent_power_lost}\n")

def tf_coil_inner_limb_scan():
    """
    Routine to read in the scan_tf_rcoil_inner_processed data and plot
    the total power vs. the inner radius of the plasma.
    """

    def create_csv():
        runs = python_scripts_run.create_runs_list(runs_directory)
        runs.sort(key=lambda run_i: (run_i.log.rcoil_inner, run_i.log.rcoil))

        runs_metadata = []
        for run_i in runs:
            output_dir_i = os.path.join(output_dir,
                                        f"rcoil_inner_{run_i.log.rcoil_inner}"
                                        f"_rcoil_{run_i.log.rcoil}")
            paper_plots.calc_energy_flux(run_i, output_dir_i)
            runs_metadata.append([run_i.log.rcoil,
                                  run_i.log.rcoil_inner,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_2d,
                                  run_i.log.pinj])
        columns = ['rcoil', 'rcoil_inner', 'max_energy_flux', 'total_energy_flux',
                   'conf_band_2d', 'conf_band_total', 'h_phi', 'h_theta_2d', 'pinj']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'rcoil_inner_runs.csv'))

    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    runs_directory = os.path.join(repository_path, "output_data", "scan_tf_rcoil_inner_processed")
    output_dir = os.path.join(repository_path, "plots", "scan_tf_rcoil_inner")
    make_csv = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'rcoil_inner_runs.csv'))
    rcoil_array = np.array([7.50, 100.0])

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)

    for rcoil in rcoil_array:
        df_i = df[df['rcoil'] == rcoil]
        axs[0].errorbar(df_i['rcoil_inner'], df_i['max_energy_flux'],
                        yerr=df_i['conf_band_2d'],
                        label=r"$R_\text{outer}$" f" = {rcoil} m")
        total_energy_flux = df_i['total_energy_flux'] * 100 / df_i['pinj']
        conf_band_total = df_i['conf_band_total'] * 100 / df_i['pinj']
        axs[1].errorbar(df_i['rcoil_inner'], total_energy_flux,
                        yerr=conf_band_total,
                        label=r"$R_\text{outer}$" f" = {rcoil} m")
    for ax in axs:
        ax.set_xlabel(r"$R_\text{inner}$ [m]")
        ax.legend()
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"Max Alpha Particle Energy Flux [MW m$^{-2}$]")
    axs[1].set_ylabel(r"Alpha Power lost [%]")

    fig.suptitle('TF Ripple Field Results')

    output_path = os.path.join(output_dir, 'max_and_total_energy_flux_vs_rcoil_inner')
    fig.savefig(output_path + '.png',
                bbox_inches='tight', dpi=300)
    fig.savefig(output_path + '.pdf',
                bbox_inches='tight')
    plt.close('all')

def plot_poincare_theta_psi(run_i, output_dir):
    """
    Plot a poincare plot for a poincare run.

    Note that run is assumed to already be initialized with the
    log and poincare attributes.
    """
    fig, ax = plt.subplots()
    fig_size = fig.get_size_inches()
    fig_size *= 5
    fig.set_size_inches(fig_size)
    ax.scatter(run_i.poincare.theta_array * 180 / np.pi,
               run_i.poincare.psin_array,
               marker='.')
    for item in ([ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        font_size = item.get_fontsize()
        item.set_fontsize(font_size * 3)
    ax.set_xlabel(r'$\theta$ [deg]')
    ax.set_ylabel(r'$\psi_N$')
    titlefont_size = ax.title.get_fontsize() * 3
    title = (f'n={run_i.log.bplasma_n}_phase={run_i.log.rmp_phase}_'
             f'current={run_i.log.rmp_current}kAt'
             f'_response={run_i.log.rmp_response}_coil_set={run_i.log.coil_set}'
             f'_i3dr={run_i.log.i3dr}')
    ax.set_title(title, fontsize=titlefont_size)
    fig.savefig(os.path.join(output_dir, f'poincare_psi_theta_{title}.png'),
                bbox_inches='tight', dpi=100)
    plt.close('all')

def poincare_rmp_toggle_i3dr_scan():
    """
    Routine to read in the poincare_rmp_toggle_i3dr_processed data and plot
    poincare plots for each simulation.
    """
    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    runs_directory = os.path.join(repository_path, "output_data",
                                  "scan_poincare_toggle_i3dr_processed")
    gfile_path = os.path.join(repository_path, "input_data", "SPR-045-16.eqdsk")
    output_dir = os.path.join(repository_path, "plots", "poincare", "rmp_toggle_i3dr")
    os.makedirs(output_dir, exist_ok=True)
    runs = python_scripts_run.create_runs_list(runs_directory, poincare_run=True)
    for run_i in runs:
        run_i.init_gfile(gfile_path)
        run_i.init_poincare()
        plot_poincare_theta_psi(run_i, output_dir)
        run_i.free_space()

def csv_of_key_values():
    """
    Saves the following key values to a csv file:
      - Absolute value of the magnetic field at the magnetic axis
      - Major radius at the magnetic axis
      - Total plasma current.
      - Total fusion power.
      - Total alpha power.
      - Zeff
    """

    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gfile_path = os.path.join(repository_path, "input_data", "SPR-045-16.eqdsk")
    gfile = my_gfile_reader.getGfile(gfile_path)
    cdf_path = os.path.join(repository_path, "input_data", "profiles_SPR-045-16.CDF")
    #pylint: disable=no-member
    with netCDF4.Dataset(cdf_path, 'r') as cdf:
        #pylint: enable=no-member
        zeff = cdf.variables['ZEFF'][0][0]
    current = 20.10
    total_fusion_power = 1.66
    total_alpha_power = 0.33

    # Write csv file
    output_dir = os.path.join(repository_path, "plots", "list_of_key_values.csv")
    with open(output_dir, "w", encoding="utf-8") as file:
        file.write(f"R_axis = {gfile.rmaxis:.2f} m\n")
        file.write(f"B_axis = {gfile.FSpline(0)/gfile.rmaxis:.2f} T\n")
        file.write(f"Plasma current = {current:.2f} MA\n")
        file.write(f"Zeff = {zeff:.2f}\n")
        file.write(f"Total fusion power = {total_fusion_power:.2f} GW\n")
        file.write(f"Total alpha power = {total_alpha_power:.2f} GW\n")

    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gfile_path = os.path.join(repository_path, "input_data", "SPR-045-14.eqdsk")
    gfile = my_gfile_reader.getGfile(gfile_path)
    cdf_path = os.path.join(repository_path, "input_data", "profiles_SPR-045-14.CDF")
    #pylint: disable=no-member
    with netCDF4.Dataset(cdf_path, 'r') as cdf:
        #pylint: enable=no-member
        zeff = cdf.variables['ZEFF'][0][0]
    current = 20.46
    total_fusion_power = 1.42
    total_alpha_power = 0.28

    # Write csv file
    output_dir = os.path.join(repository_path, "plots", "list_of_key_values_spr_045_14.csv")
    with open(output_dir, "w", encoding="utf-8") as file:
        file.write(f"R_axis = {gfile.rmaxis:.2f} m\n")
        file.write(f"B_axis = {gfile.FSpline(0)/gfile.rmaxis:.2f} T\n")
        file.write(f"Plasma current = {current:.2f} MA\n")
        file.write(f"Zeff = {zeff:.2f}\n")
        file.write(f"Total fusion power = {total_fusion_power:.2f} GW\n")
        file.write(f"Total alpha power = {total_alpha_power:.2f} GW\n")

def plot_divb_pcolourmesh():
    """
    Plot the real and imaginary components of the divergence of the magnetic field (divB)
    and visualize field components as pcolor meshes.
    """
    # Determine the repository and data paths
    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rwm_control_path = os.path.join(repository_path, "input_data", "BPLASMA", "RWM_control")

    # Get full paths of all files that start with "BPLASMA" in the rwm_control_path directory
    control_files = [os.path.join(rwm_control_path, file) \
                     for file in os.listdir(rwm_control_path) if file.startswith("BPLASMA")]

    # Grid resolution for radial and vertical directions (assuming fixed)
    num_r_points = 100
    num_z_points = 200
    toroidal_mode_number = -1

    # Initialize arrays to store data
    radial_grid = np.zeros((num_z_points, num_r_points))
    vertical_grid = np.zeros((num_z_points, num_r_points))

    # Iterate through each file to read and process data
    for filename in control_files:
        # Load the numerical data, skipping the first 3 header lines
        data = np.loadtxt(filename, skiprows=3)

        # Extract columns corresponding to R, Z, and B-field components
        radius = data[:, 0]
        height = data[:, 1]
        br_real = data[:, 2]
        br_imag = data[:, 3]
        bz_real = data[:, 4]
        bz_imag = data[:, 5]
        bt_real = data[:, 6]
        bt_imag = data[:, 7]

        # Reshape R, Z, and the magnetic fields onto a grid
        radial_grid = radius.reshape(num_r_points, num_z_points).T
        vertical_grid = height.reshape(num_r_points, num_z_points).T
        br_real_grid = br_real.reshape(num_r_points, num_z_points).T
        br_imag_grid = br_imag.reshape(num_r_points, num_z_points).T
        bz_real_grid = bz_real.reshape(num_r_points, num_z_points).T
        bz_imag_grid = bz_imag.reshape(num_r_points, num_z_points).T
        bt_real_grid = bt_real.reshape(num_r_points, num_z_points).T
        bt_imag_grid = bt_imag.reshape(num_r_points, num_z_points).T

        # Create a figure with subplots to visualize all the field components
        fig, axes = plt.subplots(1, 3)
        fig_size = fig.get_size_inches()
        fig_size[0] *= 2
        fig.set_size_inches(fig_size)

        # Plot divB Real
        divb_real = (1 / radial_grid) * np.gradient(radial_grid * br_real_grid, axis=1) \
                    / np.gradient(radial_grid, axis=1) \
                    - toroidal_mode_number * bt_imag_grid / radial_grid \
                    + np.gradient(bz_real_grid, axis=0) \
                    / np.gradient(vertical_grid, axis=0)
        mesh1 = axes[0].pcolormesh(radial_grid, vertical_grid, divb_real, shading='auto', cmap='viridis')
        fig.colorbar(mesh1, ax=axes[0])
        axes[0].set_title('divB Real')

        # Plot divB Imaginary (calculated from provided formulas)
        divb_imag = (1 / radial_grid) * np.gradient(radial_grid * br_imag_grid, axis=1) \
                    / np.gradient(radial_grid, axis=1) \
                    + toroidal_mode_number * bt_real_grid  / radial_grid \
                    + np.gradient(bz_imag_grid, axis=0) / \
                    np.gradient(vertical_grid, axis=0)
        mesh2 = axes[1].pcolormesh(radial_grid, vertical_grid, divb_imag, shading='auto', cmap='viridis')
        fig.colorbar(mesh2, ax=axes[1])
        axes[1].set_title('divB Imaginary')
        for i in range(3):
            axes[i].set_xlabel('Radius (m)')
            axes[i].set_ylabel('Height (m)')
            axes[i].set_aspect('equal')

        # Plot absolute value of B
        B_abs_avg = 0.5 * (
                    np.sqrt(br_real_grid**2 + bz_real_grid**2 + bt_real_grid**2) +
                    np.sqrt(br_imag_grid**2 + bz_imag_grid**2 + bt_imag_grid**2))
        mesh3 = axes[2].pcolormesh(radial_grid, vertical_grid, B_abs_avg, shading='auto', cmap='viridis')
        fig.colorbar(mesh3, ax=axes[2])
        axes[2].set_title('<|B|>')                

        # Create figure suptitle with last part of the filename
        fig.suptitle(filename.split('/')[-1])     

        # Save figure to output directory
        # Make a divb directory if it doesn't exist
        divb_directory = os.path.join(repository_path, 'plots', 'divb')
        os.makedirs(divb_directory, exist_ok=True)
        fig.savefig(divb_directory + '/' + filename.split('/')[-1] + '.png',
                    bbox_inches='tight', dpi=300)

        # Make new plot of just the magnetic fields
        fig, axes = plt.subplots(1, 6)
        fig_size = fig.get_size_inches()
        fig_size[0] *= 4
        fig.set_size_inches(fig_size)

        mesh = axes[0].pcolormesh(radial_grid, vertical_grid, bt_real_grid, shading='auto', cmap='viridis')
        fig.colorbar(mesh, ax=axes[0])
        axes[0].set_title('Bt Real')
        mesh = axes[1].pcolormesh(radial_grid, vertical_grid, bt_imag_grid, shading='auto', cmap='viridis')
        fig.colorbar(mesh, ax=axes[1])
        axes[1].set_title('Bt Imag')
        mesh = axes[2].pcolormesh(radial_grid, vertical_grid, br_real_grid, shading='auto', cmap='viridis')
        fig.colorbar(mesh, ax=axes[2])
        axes[2].set_title('Br Real')
        mesh = axes[3].pcolormesh(radial_grid, vertical_grid, br_imag_grid, shading='auto', cmap='viridis')
        fig.colorbar(mesh, ax=axes[3])
        axes[3].set_title('Br Imag')
        mesh = axes[4].pcolormesh(radial_grid, vertical_grid, bz_real_grid, shading='auto', cmap='viridis')
        fig.colorbar(mesh, ax=axes[4])
        axes[4].set_title('Bz Real')
        mesh = axes[5].pcolormesh(radial_grid, vertical_grid, bz_imag_grid, shading='auto', cmap='viridis')
        fig.colorbar(mesh, ax=axes[5])
        axes[5].set_title('Bz Imag')
        for i in range(6):
            axes[i].set_xlabel('Radius (m)')
            axes[i].set_ylabel('Height (m)')
            axes[i].set_aspect('equal')
        fig.suptitle(filename.split('/')[-1])
        fig.savefig(divb_directory + '/' + filename.split('/')[-1] + '_fields.png',
                    bbox_inches='tight', dpi=300)
        plt.close('all')

if __name__ == "__main__":

    # tf_coil_inner_limb_scan()
    # poincare_rmp_toggle_i3dr_scan()
    csv_of_key_values()
    # plot_divb_pcolourmesh()
