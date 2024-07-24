"""
This repo reads in a CDF file and outputs data about the ions and electrons for LOCUST.
"""

import os
import time
import netCDF4
import numpy as np
from python_scripts import dt_fusion, generate_alphas, my_gfile_reader

def get_cdf_filename(spr_string, input_data_dir):
    """
    Return the CDF filename for the given SPR string.
    """
    return os.path.join(input_data_dir, f'profiles_{spr_string}.CDF')

def read_cdf_file(cdf_filename):
    """
    Read the CDF file and return the profiles.
    """
    #pylint: disable=no-member
    with netCDF4.Dataset(cdf_filename, 'r') as profile_cdf:
        #pylint: enable=no-member
        psin = profile_cdf.variables['XPSI'][-1, :]
        ti = profile_cdf.variables['TI'][-1, :]
        te = profile_cdf.variables['TE'][-1, :]
        ne = profile_cdf.variables['NE'][-1, :]
        nd = profile_cdf.variables['NID'][-1, :]
        nt = profile_cdf.variables['NIT'][-1, :]
    return psin, ti, te, ne, nd, nt

def save_profiles(psin, profiles, output_dir, profile_names, spr_string):
    """
    Save the profiles to a file.
    """
    for profile, profile_name in zip(profiles, profile_names):
        filepath = os.path.join(output_dir, f"profile_{spr_string}_{profile_name}.dat")
        np.savetxt(filepath, np.transpose([psin, profile]), fmt='%13.5e',
                   header=str(len(psin)), comments='')

def get_impurity_data(cdf_filename, num_of_impurities):
    """
    Return the impurity names, atomic numbers and number densities for the given CDF file.
    """
    impurity_names = []
    zia = []
    nim = []

    #pylint: disable=no-member
    with netCDF4.Dataset(cdf_filename, 'r') as profile_cdf:
        #pylint: enable=no-member
        for i in range(num_of_impurities):
            impurity_data = profile_cdf.variables[f'NIM{i+1}'][-1, :]
            ziai = float(profile_cdf.variables[f'ZIA{i+1}'][-1, 0].data)
            zia.append(ziai)
            nim.append(impurity_data)

            impurity_name = get_impurity_name(ziai)
            if impurity_name is None:
                raise ValueError(f'Unknown impurity with atomic number: {ziai}')
            impurity_names.append(impurity_name)

    return impurity_names, zia, nim

def calculate_number_of_impurities(cdf_filename):
    """
    Return the number of impurities for the given CDF file.
    """
    num_of_impurities = 0
    impurity_data_exists = True
    impurity_variable = None  # Assign the impurity variable to None
    #pylint: disable=no-member
    with netCDF4.Dataset(cdf_filename, 'r') as profile_cdf:
        #pylint: enable=no-member
        while impurity_data_exists:
            try:
                impurity_variable = profile_cdf.variables['NIM' + str(num_of_impurities + 1)][-1, :]
                if len(impurity_variable) == 0:
                    raise ValueError('Impurity data is empty')
                num_of_impurities += 1
            except KeyError:
                impurity_data_exists = False
        return num_of_impurities

def get_impurity_name(ziai):
    """
    Return the impurity name for the given atomic number.
    """
    if int(np.round(ziai)) == 2:
        return 'AHe4'
    elif int(np.round(ziai)) == 18:
        return 'AAr'
    elif 50 <= ziai <= 54:
        return 'AXe'
    else:
        return None

def calculate_ion_fractions(ne, nd, nt, nim):
    """
    Calculate the ion fractions. 
    """
    fid = nd[0] / ne[0]
    fit = nt[0] / ne[0]
    fim = [nimi[0] / ne[0] for nimi in nim]
    return fid, fit, fim

def write_ion_info_file(filename, num_of_impurities, alpha_power, fid, fit, fim, zia,
                        impurity_names):
    """
    Write the ion info to a file.
    """
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f'{num_of_impurities + 2}\n')
        f.write(f'{alpha_power:.1e}\n')
        f.write(f'{1:11.8f} {fid:11.8f} AD\n')
        f.write(f'{1:11.8f} {fit:11.8f} AT\n')
        for ziai, fimi, impurity_name in zip(zia, fim, impurity_names):
            f.write(f'{ziai:11.8f} {fimi:11.8f} {impurity_name}\n')

def create_dirs(input_dir, output_dir):
    """
    Create the input and output directories if they do not exist.
    """
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def get_gfile(gfile_path):
    """
    Return the gfile object for the given gfile path.
    """
    return my_gfile_reader.getGfile(gfile_path)

def main(spr_string, input_dir, output_dir, num_markers):
    """
    Main function.

    Note that LOCUST assumes a constant Zeff so for consistency we will
    also use a constant Zeff. Namely, we use the Zeff at the magnetic
    axis.
    """
    start_time = time.time()

    create_dirs(input_dir, output_dir)
    cdf_filename = get_cdf_filename(spr_string, input_dir)
    psin, ti, te, ne, nd, nt = read_cdf_file(cdf_filename)
    gfile_path = os.path.join(input_dir, f'{spr_string}.eqdsk')
    gfile = get_gfile(gfile_path)
    num_impurities = calculate_number_of_impurities(cdf_filename)
    impurity_names, zia, nim = get_impurity_data(cdf_filename, num_impurities)
    fd, ft, fim = calculate_ion_fractions(ne, nd, nt, nim)
    fusion_power, alpha_power = dt_fusion.calc_fusion_power(psin, fd * ne, ft * ne,
                                                            ti, gfile)
    print(f'Fusion power: {fusion_power:.3e} W')
    print(f'Alpha power: {alpha_power:.3e} W')
    ion_info_path = os.path.join(output_dir, f'ion_info_{spr_string}.dat')
    write_ion_info_file(ion_info_path, num_impurities, alpha_power, fd, ft, fim, zia,
                        impurity_names)
    profiles = [te, ti, ne]
    profile_names = ["Te", "Ti", "ne"]
    save_profiles(psin, profiles, output_dir, profile_names, spr_string)

    # Generate markers uniformly across the plasma volume
    # with weights according to the DT reaction rate
    marker_coords = generate_alphas.gen_marker_coords(num_markers, gfile)
    marker_weights = generate_alphas.gen_marker_weights(marker_coords, gfile, psin, ne * fd,
                                                        ne * ft, ti)
    e_alpha_std_coords = dt_fusion.calc_e_alpha_std(marker_coords, gfile, ti, psin)
    marker_velocities = generate_alphas.gen_marker_velocities(e_alpha_std_coords)
    marker_filename = f"{spr_string}_markers_{num_markers:d}.dat"
    marker_filepath = os.path.join(output_dir, marker_filename)
    generate_alphas.save_markers(marker_filepath, marker_coords, marker_velocities,
                                 marker_weights)

    end_time = time.time()
    print(f'Time taken for run main: {end_time - start_time:.2e} seconds')

if __name__ == "__main__":
    SPR_STRING = 'SPR-045-16'
    CURRENT_DIR = os.path.dirname(__file__)
    INPUT_DATA_DIR = os.path.join(CURRENT_DIR, '..', 'input_data')
    OUTPUT_DIR = INPUT_DATA_DIR
    NUM_MARKERS = int(1e6)
    main(SPR_STRING, INPUT_DATA_DIR, OUTPUT_DIR, NUM_MARKERS)
