"""
The module contains routines associated with the energy flux on the PFCs.
It also contains information about the distribution of the energy flux on the PFCs.
"""

import time
from typing import Optional
import numpy as np
from scipy.ndimage import minimum_filter
from python_scripts import pkde

class Flux:
    """
    This class contains routines for calculating the energy flux on the PFCs
    """
    def __init__(self,
                 num_grid_points_1d=10**4,
                 num_grid_points_2d=10**3,
                 num_bootstraps=128):
        self.num_grid_points_1d = num_grid_points_1d
        self.num_grid_points_2d = num_grid_points_2d
        self.num_bootstraps = num_bootstraps
        self.s_theta_1d: Optional[np.ndarray] = None
        self.s_theta_2d: Optional[np.ndarray] = None
        self.s_phi: Optional[np.ndarray] = None
        self.total_energy: Optional[float] = None # Store the total energy flux in MW on the PFCs
        self.h_phi: Optional[float] = None
        self.h_theta_1d: Optional[float] = None
        self.h_theta_2d: Optional[float] = None
        self.energy_1d: Optional[np.ndarray] = None
        self.energy_2d: Optional[np.ndarray] = None
        self.h_phi_array: Optional[np.ndarray] = None
        self.h_theta_1d_array: Optional[np.ndarray] = None
        self.h_theta_2d_array: Optional[np.ndarray] = None
        self.amise_1d: Optional[np.ndarray] = None
        self.amise_2d: Optional[np.ndarray] = np.array([[None]])
        self.conf_band_1d: Optional[float] = None
        self.conf_band_2d: Optional[float] = None
        self.conf_band_total: Optional[float] = None
        self.max_energy_1d: Optional[float] = None
        self.max_energy_2d: Optional[float] = None

def calc_s_theta_s_phi(run):
    """
    Calculate the s_theta and s_phi attributes assuming that
    the wall attribute of the run object is populated.
    """
    run.flux.s_phi = np.linspace(run.wall.s_phi_min,
                                 run.wall.s_phi_max,
                                 run.flux.num_grid_points_2d)
    run.flux.s_theta_1d = np.linspace(run.wall.s_theta_min,
                                      run.wall.s_theta_max,
                                      run.flux.num_grid_points_1d)
    run.flux.s_theta_2d = np.linspace(run.wall.s_theta_min,
                                      run.wall.s_theta_max,
                                      run.flux.num_grid_points_2d)

def calc_total_energy(run,
                      stopped_energy=None,
                      all_energy0=None,
                      stopped_weight=None,
                      all_weight=None):
    """
    Calculate the total energy flux on the Plasma-Facing Components (PFCs) in megawatts (MW).

    This method relies on specific attributes and states of the provided `run` object. 
    Before calling this method, ensure that the `run` object meets the following conditions:
    
    1. It has been initialized with the Flux class attribute.
    2. The `markers.stopped` attribute has been populated with relevant data.
    3. The `markers.all` attribute has been populated with relevant data.
    4. The `log` attribute has been populated and includes the `pinj` attribute, 
    which specifies the total power injected in MW.

    Args:
        run: A Run object
            An instance of the Run class that meets the above requirements.

    Returns:
        float: The total energy flux on the PFCs in MW. If the requirements are not 
            met, the behavior of this method is undefined and may result in errors.
    """
    if stopped_energy is None:
        stopped_energy = run.markers.stopped.energy
    if stopped_weight is None:
        stopped_weight = run.markers.stopped.weight
    if all_energy0 is None:
        all_energy0 = run.markers.all.energy0
    if all_weight is None:
        all_weight = run.markers.all.weight
    return run.log.pinj * np.sum(stopped_energy * stopped_weight) \
                        / np.sum(all_energy0 * all_weight)

def calc_energy_flux_1d(run,
                        stopped_s_theta=None,
                        stopped_energy=None,
                        stopped_weight=None,
                        h_theta_1d=None):
    """
    Calculate the energy flux on the PFCs using the distribution of the
    particle flux on the PFCs.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_theta attribute, num_grid_points and the
            h_theta_1d attribute (if h_theta_1d is None).

    Returns:
        Flux object
            The Flux object with new attributes, namely the energy_fun_1d
            and energy_arr_1d which give the energy flux along the wall
            in MW/m^2.

    """
    if stopped_s_theta is None:
        stopped_s_theta = run.markers.stopped.s_theta
    if stopped_energy is None:
        stopped_energy = run.markers.stopped.energy
    if stopped_weight is None:
        stopped_weight = run.markers.stopped.weight
    if h_theta_1d is None:
        h_theta_1d = run.flux.h_theta_1d
    kde_fun = pkde.periodic_kde_1d(stopped_s_theta,
                                   run.wall.s_theta_min,
                                   run.wall.s_theta_max,
                                   stopped_weight * stopped_energy,
                                   h_theta_1d,
                                   kernel='gaussian',
                                   num_grid_points=run.flux.num_grid_points_1d)
    kde_array = kde_fun(run.flux.s_theta_1d)

    # Calculate bias corrected KDE
    asymptotic_bias = \
        pkde.calc_asymptotic_bias_1d(kde_fun, run.flux.s_theta_1d, h_theta_1d)
    kde_array = kde_array - asymptotic_bias

    # Calculate unstretched KDE
    kde_array *= run.wall.scale_factor_1d(run.flux.s_theta_1d)

    # Calculate energy flux
    energy_flux = kde_array * run.flux.total_energy

    return energy_flux

def calc_energy_flux_2d(run,
                        stopped_s_phi=None,
                        stopped_s_theta=None,
                        stopped_energy=None,
                        stopped_weight=None,
                        h_phi=None,
                        h_theta_2d=None):
    """
    Calculate the energy flux on the PFCs using the distribution of the
    particle flux on the PFCs.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_phi and s_theta attributes, num_grid_points
            and the h_phi and h_theta_2d attributes (if h_theta_2d is None).

    Returns:
        Flux object
            The Flux object with new attributes, namely the energy_fun_2d
            and energy_arr_2d which give the energy flux along the wall
            in MW/m^2.

    """

    if stopped_s_phi is None:
        stopped_s_phi = run.markers.stopped.s_phi
    if stopped_s_theta is None:
        stopped_s_theta = run.markers.stopped.s_theta
    if stopped_energy is None:
        stopped_energy = run.markers.stopped.energy
    if stopped_weight is None:
        stopped_weight = run.markers.stopped.weight
    if h_theta_2d is None:
        h_theta_2d = run.flux.h_theta_2d
    if h_phi is None:
        h_phi = run.flux.h_phi

    kde_fun = pkde.periodic_kde_2d(stopped_s_phi,
                                   stopped_s_theta,
                                   run.wall.s_phi_min,
                                   run.wall.s_phi_max,
                                   run.wall.s_theta_min,
                                   run.wall.s_theta_max,
                                   stopped_weight * stopped_energy,
                                   h_phi,
                                   h_theta_2d,
                                   num_grid_points=run.flux.num_grid_points_2d)
    kde_array = kde_fun(run.flux.s_phi, run.flux.s_theta_2d)

    # Calculate bias corrected KDE
    asymptotic_bias = \
        pkde.calc_asymptotic_bias_2d(kde_fun, run.flux.s_phi, run.flux.s_theta_2d,
                                     h_phi, h_theta_2d)
    kde_array = kde_array - asymptotic_bias

    # Calculate unstretched KDE
    kde_array *= run.wall.scale_factor_2d(run.flux.s_phi, run.flux.s_theta_2d)

    # Calculate energy flux
    energy_flux = kde_array * run.flux.total_energy
    return energy_flux

def calc_optimum_bandwidth_1d(run, h_theta_1d_array=None):
    """
    Calculate the optimal bandwidth for the 1D KDE. Note that we assume
    that a unique optimum bandwidth is possible.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_theta attribute and num_grid_points

        h_theta_1d_array: array
            The set of bandwidths to be tested for.

    Returns:
        h_theta_1d: float
            The optimal bandwidth
        amise_array: array
            The AMISE for each of the bandwidths.
    """
    if h_theta_1d_array is None:
        h_theta_1d_array = run.flux.h_theta_1d_array
    amise_array = pkde.calc_amise_1d_array(run.markers.stopped.s_theta,
                                           run.flux.s_theta_1d,
                                           run.markers.stopped.weight * run.markers.stopped.energy,
                                           h_theta_1d_array,
                                           num_grid_points=run.flux.num_grid_points_1d)
    h_theta_1d = h_theta_1d_array[np.argmin(amise_array)]
    run.flux.amise_1d = amise_array
    run.flux.h_theta_1d = h_theta_1d

def calc_optimum_bandwidth_2d(run, h_phi_array=None, h_theta_2d_array=None):
    """
    Calculate the optimal bandwidth for the 2D KDE. Note that we assume
    that a unique optimum bandwidth combination is possible. Note we take
    the optimum bandwidth to be the minima of the AMISE with the largest
    h_phi^2 + h_theta_2d^2.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_phi and s_theta attributes, num_grid_points
            and the h_phi and h_theta_2d attributes (if h_theta_2d is None).

        h_phi_array: array
            The set of bandwidths to be tested for.
        h_theta_2d_array: array
            The set of bandwidths to be tested for.

    Returns:
        h_phi: float
            The optimal bandwidth
        h_theta_2d: float
            The optimal bandwidth
        amise_array: array
            The AMISE for each of the bandwidths.
    """
    print('Calculating optimum 2d bandwidths')
    start_time = time.time()
    if h_phi_array is None:
        h_phi_array = run.flux.h_phi_array
    if h_theta_2d_array is None:
        h_theta_2d_array = run.flux.h_theta_2d_array
    amise_array = pkde.calc_amise_2d_array(run.markers.stopped.s_phi,
                                           run.markers.stopped.s_theta,
                                           run.flux.s_phi,
                                           run.flux.s_theta_2d,
                                           run.markers.stopped.weight * run.markers.stopped.energy,
                                           h_phi_array,
                                           h_theta_2d_array,
                                           num_grid_points=run.flux.num_grid_points_2d)
    # min_indices = np.unravel_index(np.argmin(amise_array), amise_array.shape)
    minima = minimum_filter(amise_array, size=3, mode = 'nearest') == amise_array
    min_h_coords = np.argwhere(minima)
    min_h_coord = min_h_coords[np.argmin(min_h_coords[:, 0]**2 + min_h_coords[:, 1]**2)]
    h_phi = h_phi_array[min_h_coord[1]]
    h_theta_2d = h_theta_2d_array[min_h_coord[0]]
    run.flux.amise_2d = amise_array
    run.flux.h_phi = h_phi
    run.flux.h_theta_2d = h_theta_2d
    print(f'Optimum 2d bandwidths calculated in {time.time() - start_time} seconds')
