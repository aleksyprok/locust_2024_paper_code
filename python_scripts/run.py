"""
This module contains routines for updating a Run object which stores information
about a LOCUST run.

The run object has four classes as attributes:
- log: contains information obtained by reading the LOG*.out file
- wall: contains routines from converting from R, Phi, Z coordinates on the wall to
        wall coordinates, s_phi, s_theta.
- markers: contains routines related to the particle groups from reading the FINAL_STATE*.dat file.
- flux: contains routines related to the energy flux on the PFCs.
Note that the log wall and markers classes are standalone classes, but many of the methods in the
flux class require the log, wall and markers classes.
Note that to save memory we will often set the markers attribute to None.
"""
import glob
import os
from typing import Optional
from python_scripts import flux, log, markers, my_gfile_reader, poincare, wall

class Run:
    """
    This class contains routines for updating a Run object which stores information
    about a LOCUST run.
    """
    def __init__(self, dir_path, tag):
        self.log: Optional[log.Log] = None
        self.wall: Optional[wall.Wall] = None
        self.markers: Optional[markers.Markers] = None
        self.gfile: Optional[my_gfile_reader.getGfile] = None
        self.flux: Optional[flux.Flux] = None
        self.poincare: Optional[poincare.Poincare] = None
        self.dir_path = dir_path
        self.tag = tag
        self.log_path = self.dir_path + f'/LOG_{self.tag}.out'
        self.fstate_path = self.dir_path + f'/FINAL_STATE_{self.tag}.dat'

    def init_log(self):
        """
        Initialize the Run object with information from the LOCUST .log file.
        """
        self.log = log.Log(self.log_path)

    def init_wall(self, wall_path,
                  special_nodes=(0, 1, 2, 3)):
        """
        Initialize the Run object with information from the LOCUST .log file.
        """
        self.wall = wall.Wall(wall_path,
                              special_nodes=special_nodes)

    def init_markers(self,
                     remap_phi_n=None):
        """
        Initialize the Run object with information from the LOCUST FINAL_STATE*.dat file.
        """
        self.markers = markers.Markers(self.fstate_path)
        if self.wall is not None:
            markers.get_s_phi_s_theta_from_r_z_phi(self,
                                                   remap_phi_n=remap_phi_n)

    def init_gfile(self, gfile_path):
        """
        Initialize the Run object with information from the eqdsk file.
        """
        self.gfile = my_gfile_reader.getGfile(gfile_path)

    def init_flux(self,
                  num_grid_points_1d=10**4,
                  num_grid_points_2d=10**3,
                  num_bootstraps=128,
                  h_theta_1d_array=None,
                  h_theta_2d_array=None,
                  h_phi_array=None):
        """
        Initialize the Run object with the flux class based on
        desired num_grid_points, num_h_theta_1d, num_h_theta_2d,
        num_h_phi.
        """
        self.flux = flux.Flux(num_grid_points_1d=num_grid_points_1d,
                              num_grid_points_2d=num_grid_points_2d,
                              num_bootstraps=num_bootstraps)
        if self.wall is not None:
            flux.calc_s_theta_s_phi(self)
        if self.markers is not None:
            self.flux.total_energy = flux.calc_total_energy(self)
        if h_theta_1d_array is not None:
            self.flux.h_theta_1d_array = h_theta_1d_array
        if h_theta_2d_array is not None:
            self.flux.h_theta_2d_array = h_theta_2d_array
        if h_phi_array is not None:
            self.flux.h_phi_array = h_phi_array

    def init_poincare(self):
        """
        Initialize the Run object with information from the Poincare map file.
        """
        poincare_path = self.dir_path + f'/Poincare_map_{self.tag}_formatted.dat'
        self.poincare = poincare.Poincare(poincare_path,
                                          gfile=self.gfile)

    def free_space(self):
        """
        Free the memory used by the Run object.
        """
        self.markers = None
        self.gfile = None
        self.wall = None
        if self.flux is not None:
            self.flux.s_phi = None
            self.flux.s_theta_1d = None
            self.flux.s_theta_2d = None
            self.flux.energy_1d = None
            self.flux.energy_2d = None
        self.poincare = None

def create_runs_list(dir_path,
                     poincare_run=False):
    """
    Create a list of Run objects.

    Args:
        dir_path: str
            The path to the directory containing the runs.
    
    Returns:
        list
            A list of Run objects.
    """
    runs = []
    if not poincare_run:
        for file in glob.iglob(f'{dir_path}/**/FINAL_STATE*.dat', recursive=True):
            tag = os.path.splitext(os.path.basename(file))[0].split('FINAL_STATE_')[-1]
            run_i = Run(os.path.dirname(file), tag)
            run_i.init_log()
            runs.append(run_i)
    else:
        for file in glob.iglob(f'{dir_path}/**/Poincare_map*.dat', recursive=True):
            tag = os.path.splitext(os.path.basename(file))[0].split('Poincare_map_')[-1]
            tag = tag.split('_formatted')[0]
            run_i = Run(os.path.dirname(file), tag)
            run_i.init_log()
            runs.append(run_i)
    return runs
