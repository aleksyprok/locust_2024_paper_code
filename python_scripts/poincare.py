"""
This module contains routines for updating a Run object from the
runs.py module with information from a LOCUST Poincare_map_*_formatted file.
"""

import numpy as np

class Poincare:
    """
    This class contains routines for updating a Poincare object which stores information
    about a LOCUST run with the -DPOINCARE=3 flag enabled.
    """

    def __init__(self,
                 poincare_path: str,
                 gfile=None):
        """
        Args:
            poincare_path: path to the Poincare map file
        """

        self.poincare_path = poincare_path
        self.r_array, self.z_array = read_poincare_file(poincare_path)
        self.psin_array = None
        self.theta_array = None
        if gfile is not None:
            self.calc_theta_psin_values(gfile)

    def calc_theta_psin_values(self, gfile):
        """
        Calculate the psi and theta values from the r_array and z_array values
        and using the gfile object.

        Args:
            gfile: gfile object

        Returns:
            psi_array, theta_array: arrays of psi and theta values
        """
        self.theta_array = np.arctan2(self.z_array, self.r_array - gfile.rmaxis)
        psi_array = gfile.PsiSpline(self.r_array, self.z_array, grid=False)
        self.psin_array = (psi_array - gfile.simag) / (gfile.sibry - gfile.simag)

def map_grid_to_rz_array(poincare_map,
                         r_min, r_max, z_min, z_max):
    """
    Map the grid to a list of r and z values. We only record the values where 
    the grid value is greater than 1

    Input

    poincare_map: 2D array of shape (nzp, nrp)

    r_min, r_max, z_min, z_max: float

    Output:

    r_array, z_array: arrays of r and z values which correspond
    to the cell centres of the cells of poincare_map where the entry is
    greater than 0.
    """
    nzp, nrp = poincare_map.shape
    r_vertices = np.linspace(r_min, r_max, nrp+1)
    z_vertices = np.linspace(z_min, z_max, nzp+1)
    r_centres = 0.5 * (r_vertices[:-1] + r_vertices[1:])
    z_centres = 0.5 * (z_vertices[:-1] + z_vertices[1:])
    r_list = []
    z_list = []
    for i in range(nzp):
        for j in range(nrp):
            if poincare_map[i, j] > 0:
                r_list.append(r_centres[j])
                z_list.append(z_centres[i])
    r_array = np.array(r_list)
    z_array = np.array(z_list)
    return r_array, z_array

def read_poincare_file(poincare_path):
    """
    Read the Poincare map file and return corresponding lists of R and Z.
    """
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

    return map_grid_to_rz_array(poincare_map, r_min, r_max, z_min, z_max)
