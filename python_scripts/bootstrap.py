"""
This module contains the routines needed to calculate the bootstrap error
for the maximum energy flux and total energy flux on the PFCs.
"""
import time
# pylint: disable=no-name-in-module
from multiprocess import Process, Manager, cpu_count
# pylint: enable=no-name-in-module
import numpy as np
from sklearn.utils import resample
from python_scripts import flux

def bootstrap_resample_stopped(run):
    """
    Resample the run object using the resample function from sklearn.
    This function is intended to be used in calc_be_1d_array() and
    calc_be_2d_array().

    Args:
        run: Run object which needs to have been initialized with
            the markers.stopped and wall attribute to ensure
            markers.stopped.s_phi and markers.stopped.s_theta are not None.
    Returns:
        s_phi, s_theta, energy, weight
    """
    indexes = resample(np.arange(len(run.markers.stopped.r),
                                 dtype = int))
    s_phi = run.markers.stopped.s_phi[indexes]
    s_theta = run.markers.stopped.s_theta[indexes]
    weight = run.markers.stopped.weight[indexes]
    energy = run.markers.stopped.energy[indexes]
    return s_phi, s_theta, energy, weight

def bootstrap_resample_all(run):
    """
    Resample the run object using the resample function from sklearn.
    This function is intended to be used in calc_be_total_array()

    Args:
        run: Run object which needs to have been initialized with
            the markers.all

    Returns:
        energy, energy0, stopped_weight, all_weight
    """
    indexes = resample(np.arange(len(run.markers.all.r),
                                 dtype = int))
    statuses = run.markers.all.s[indexes]
    stopped_indices = np.where(statuses == -5)[0]
    energy = run.markers.all.energy[indexes]
    energy = energy[stopped_indices]
    energy0 = run.markers.all.energy0[indexes]
    stopped_weight = run.markers.all.weight[indexes]
    stopped_weight = stopped_weight[stopped_indices]
    all_weight = run.markers.all.weight[indexes]
    return energy, energy0, stopped_weight, all_weight

def calc_be_1d_array(run,
                     num_bootstraps=None):
    """
    Calculate an array of energy flux 1D values using bootstrap resampling so we can
    construct confidence intervals.

    Note that BE stands for bootstrap error.

    Note that the run object needs to have been initialized with
    the flux class attribute, also the wall and markers.stopped attributes need to
    have been populated. The flux object needs to have been initialized with
    the s_theta and num_grid_points attributes and the energy_flux_1d attribute as well
    as the h_theta_1d attribute.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_theta attribute and num_grid_points

        num_bootstraps: int
            The number of bootstrap samples to be used to calculate BE.

    Returns:
        be_array: array
            The array of BE values
    """

    start_time = time.time()

    if num_bootstraps is None:
        num_bootstraps = run.flux.num_bootstraps

    def calc_be(n_array, be_list):

        np.random.seed()

        for _ in n_array:

            _, s_theta, energy, weight = bootstrap_resample_stopped(run)
            energy_flux_1d = flux.calc_energy_flux_1d(run,
                                                      stopped_s_theta=s_theta,
                                                      stopped_energy=energy,
                                                      stopped_weight=weight)
            be_list.append(np.max(np.abs(energy_flux_1d - run.flux.energy_1d)))


    be_list = []
    n_array_list = np.array_split(np.arange(num_bootstraps), cpu_count())

    with Manager() as manager:
        be_list = manager.list()
        processes = []

        for n_array in n_array_list:
            # pylint: disable=not-callable
            p = Process(target = calc_be, \
                args = (n_array, be_list))
            # pylint: enable=not-callable
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

        be_list = list(be_list)

    end_time = time.time()
    print(f"Time taken to calculate BE 1D array: {end_time - start_time:.2f} seconds")

    return np.array(be_list)

def calc_be_2d_array(run,
                     num_bootstraps=None):
    """
    Calculate an array of energy flux 2D values using bootstrap resampling so we can
    construct confidence intervals.

    Note that BE stands for bootstrap error.

    Note that the run object needs to have been initialized with
    the flux class attribute, also the wall and markers.stopped attributes need to
    have been populated. The flux object needs to have been initialized with
    the s_theta and num_grid_points attributes and the energy_flux_1d attribute as well
    as the h_theta_1d attribute.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_theta attribute and num_grid_points

        num_bootstraps: int
            The number of bootstrap samples to be used to calculate BE.

    Returns:
        be_array: array
            The array of BE values
    """

    start_time = time.time()

    if num_bootstraps is None:
        num_bootstraps = run.flux.num_bootstraps

    def calc_be(n_array, be_list):

        np.random.seed()

        for _ in n_array:

            s_phi, s_theta, energy, weight = bootstrap_resample_stopped(run)
            energy_flux_2d = flux.calc_energy_flux_2d(run,
                                                      stopped_s_phi=s_phi,
                                                      stopped_s_theta=s_theta,
                                                      stopped_energy=energy,
                                                      stopped_weight=weight)
            be_list.append(np.max(np.abs(energy_flux_2d - run.flux.energy_2d)))


    be_list = []
    n_array_list = np.array_split(np.arange(num_bootstraps), cpu_count())

    with Manager() as manager:
        be_list = manager.list()
        processes = []

        for n_array in n_array_list:
            # pylint: disable=not-callable
            p = Process(target = calc_be, \
                args = (n_array, be_list))
            # pylint: enable=not-callable
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

        be_list = list(be_list)

    end_time = time.time()
    print(f"Time taken to calculate BE 2D array: {end_time - start_time:.2f} seconds")

    return np.array(be_list)

def calc_be_total_array(run,
                        num_bootstraps=None):
    """
    Calculate an array of total energy flux values using bootstrap resampling so we can
    construct confidence intervals.

    Note that BE stands for bootstrap error.

    Note that the run object needs to have been initialized with
    the flux class attribute, also the wall and markers.stopped attributes need to
    have been populated. The flux object needs to have been initialized with
    the s_theta and num_grid_points attributes and the energy_flux_1d attribute as well
    as the h_theta_1d attribute.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_theta attribute and num_grid_points

        num_bootstraps: int
            The number of bootstrap samples to be used to calculate BE.

    Returns:
        be_array: array
            The array of BE values
    """

    start_time = time.time()

    if num_bootstraps is None:
        num_bootstraps = run.flux.num_bootstraps

    def calc_be(n_array, be_list):

        np.random.seed()

        for _ in n_array:

            energy, energy0, stopped_weight, all_weight = bootstrap_resample_all(run)
            total_energy_flux = flux.calc_total_energy(run,
                                                       stopped_energy=energy,
                                                       all_energy0=energy0,
                                                       stopped_weight=stopped_weight,
                                                       all_weight=all_weight)
            be_list.append(np.max(np.abs(total_energy_flux - run.flux.total_energy)))


    be_list = []
    n_array_list = np.array_split(np.arange(num_bootstraps), cpu_count())

    with Manager() as manager:
        be_list = manager.list()
        processes = []

        for n_array in n_array_list:
            # pylint: disable=not-callable
            p = Process(target = calc_be, \
                args = (n_array, be_list))
            # pylint: enable=not-callable
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

        be_list = list(be_list)

    end_time = time.time()
    print(f"Time taken to calculate BE total array: {end_time - start_time:.2f} seconds")

    return np.array(be_list)
