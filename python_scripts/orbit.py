"""
This file contains routines related to reading ORBIT files and
producing plots of the particle orbits. Note that the
DPLOT flag is needed to produce the ORBIT files in LOCUST.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from  python_scripts import run

class Orbit:
    def __init__(self,
                 orbit_path: str):
        self.orbit_path = orbit_path
        self.num_timesteps = 0
        self.num_markers = 0
        self.r = np.zeros(0)
        self.phi = np.zeros(0)
        self.z = np.zeros(0)

def read_orbit(orbit_path):

    # Read the file and process the data in chunks
    with open(orbit_path, 'r', encoding='utf-8') as file:
        num_markers = int(file.readline().strip())

    # Read last line
    with open(orbit_path, "rb") as file:
        try:
            file.seek(-2, os.SEEK_END)
            while file.read(1) != b'\n':
                file.seek(-2, os.SEEK_CUR)
        except OSError:
            file.seek(0)
        num_timesteps = int(file.readline().decode())

    data = np.loadtxt(orbit_path, skiprows=1, max_rows=num_timesteps*num_markers)
    data = data.reshape((num_timesteps, 3, num_markers))

    return data

def generate_markers():
    """
    Function to generate a ptcles.dat which contains the initial coordinates
    of markers that hit hotspots.

    Take the 50 highest weighted stopped markers with final z greater than 6m
    and 50 highest weighted stopped markers with final z less than -6m.
    """
    def get_top_markers_indices(z_values, weights, z_condition, top_n=50):
        """
        Get the indices of the top N highest weighted markers based on a z condition.

        Parameters:
        z_values (numpy array): The z coordinates of the markers.
        weights (numpy array): The weights of the markers.
        z_condition (callable): A function that takes z_values and returns a boolean array.
        top_n (int): Number of top markers to return. Default is 50.

        Returns:
        numpy array: Indices of the top N highest weighted markers meeting the z condition.
        """
        # Find the indices of markers meeting the z condition
        condition_indices = np.where(z_condition(z_values))[0]
        # Get the weights of these markers
        condition_weights = weights[condition_indices]
        # Sort the weights in descending order and get the top N indices
        top_indices_sorted = np.argsort(condition_weights)[::-1][:top_n]
        # Get the original indices of these top N markers
        top_markers_indices = condition_indices[top_indices_sorted]
        return top_markers_indices
    
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    run_dir = os.path.join(repo_path, "output_data",
                           "full_3d_field_spr_045_14_and_16_processed",
                           "gpu-q-52")
    run_tag = "25-07-2024_12-17-50.555"
    run0 = run.Run(run_dir, run_tag)
    run0.init_log()
    run0.init_markers()

    probabilities = run0.markers.stopped.weight
    probabilities *= (run0.markers.stopped.z > 6)
    probabilities /= np.sum(probabilities)
    upper_indices = np.random.choice(len(values), size=3, replace=False, p=probabilities)
    

if __name__ == "__main__":
    # REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # ORBIT_PATH = os.path.join(REPOSITORY_PATH, "output_data", "fec_2024_orbit_runs",
    #                           "gpu-q-57", "ORBIT_25-07-2024_14-33-04.121.dat")
    # ORBIT_PATH = os.path.join(REPOSITORY_PATH, "output_data", "fec_2024_orbit_runs",
    #                           "gpu-q-33", "ORBIT_25-07-2024_17-51-57.537.dat")
    # orbit_data = read_orbit(ORBIT_PATH)
    # marker_index = 3
    # fig, ax = plt.subplots()
    # # Scatter plot in R Z with colour given by timestep
    # timesteps = np.arange(orbit_data.shape[0])
    # scatter = ax.scatter(orbit_data[::1, 0, marker_index], orbit_data[::1, 2, marker_index],
    #                      c=timesteps, cmap='viridis', s=0.1)
    # fig.colorbar(scatter)
    # ax.set_xlabel("r")
    # ax.set_ylabel("z")
    # plt.show()
    # ax.set_aspect('equal')
    generate_markers()