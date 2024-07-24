"""
This module contains routines for calculating kernel desnity estimates
over 1D and 2D domains. With periodic boundary conditions.
"""
# pylint: disable=no-name-in-module
from multiprocess import Process, Manager, cpu_count, shared_memory
# pylint: enable=no-name-in-module
import numpy as np
from scipy import interpolate
from KDEpy import FFTKDE

def extend_coords(x, y, x_min, x_max, y_min, y_max, weights):
    """
    Repeats the coordinates in each dimension by assuming periodicity
    to ensure KDE is well resolved at the boundaries.

    Parameters
    ----------
    x : array
        1D array of x coordinates
    y : array
        1D array of y coordinates
    x_min : float
        Minimum x coordinate
    x_max : float
        Maximum x coordinate
    y_min : float
        Minimum y coordinate
    y_max : float
        Maximum y coordinate
    weights : array
        1D array of weights

    Returns
    -------
    x_extended : array
        1D array of extended x coordinates
    y_extended : array 
        1D array of extended y coordinates
    weights_extended : array
        1D array of extended weights
    """

    period_x = x_max - x_min
    period_y = y_max - y_min

    # Repeat array in each dimension
    x_block =  np.concatenate((x, x, x))
    x_extended = np.concatenate((x_block - period_x, x_block, x_block + period_x))

    y_block = np.concatenate((y - period_y, y, y + period_y))
    y_extended = np.concatenate((y_block, y_block, y_block))

    weights_block = np.concatenate((weights, weights, weights))
    weights_extended = np.concatenate((weights_block, weights_block, weights_block))

    delta_x = period_x / 2
    delta_y = period_y / 2

    reduced_indices = np.where(
        (x_min - delta_x <= x_extended) *
                           (x_extended <= x_max + delta_x) *
        (y_min - delta_y <= y_extended) *
                           (y_extended <= y_max + delta_y))
    x_extended = x_extended[reduced_indices]
    y_extended = y_extended[reduced_indices]
    weights_extended = weights_extended[reduced_indices]

    return x_extended, y_extended, weights_extended

def periodic_kde_1d(x, x_min, x_max, weights, hx,
                    kernel='gaussian',
                    num_grid_points=10**3):
    """
    Calculates a 1D kernel density estimate using FFTs with periodic
    boundary conditions. We simulate periodic boundary conditions by repeating
    the data in each dimension.

    Parameters
    ----------
    x : array
        1D array of x coordinates (corresponding to the final position of markers
        on the wall of the tokamak.)
    x_min : float
        Minimum x coordinate
    x_max : float
        Maximum x coordinate
    weights : array
        1D array of weights
    hx : float
        Bandwidth in metres
    kernel : str
        KDE kernel
    num_grid_points : int
        Number of grid points

    Returns
    -------
    pdf_fun : function
        KDE estimate.
    """

    # Extending coordiantes
    period_x = x_max - x_min
    x_norm = np.concatenate((x - period_x, x, x + period_x))
    weights_norm = np.concatenate((weights, weights, weights))
    delta_x = period_x / 2
    reduced_indices = np.where(
        (x_min - delta_x <= x_norm) *
                           (x_norm <= x_max + delta_x))
    x_norm = x_norm[reduced_indices]
    weights_norm = weights_norm[reduced_indices]
    # Double number of grid points to account for extention
    num_grid_points *= 2

    # Stretching coordiantes
    x_norm = x_norm / hx

    # Calculating extended and stretched KDE
    fft_kde = FFTKDE(kernel = kernel, bw = 1)
    grid_norm_x, pdf = \
        fft_kde.fit(x_norm, weights = weights_norm).evaluate(num_grid_points)
    pdf *= np.sum(weights_norm) / np.sum(weights) / hx # Renormalise

    # Evaluating extended and stretched KDE on an unstretched grid
    grid_x = hx * grid_norm_x
    pdf_fun = interpolate.UnivariateSpline(grid_x, pdf,
                                           s=0, ext='raise')

    return pdf_fun

def periodic_kde_2d(x, y, x_min, x_max, y_min, y_max, weights,
                    hx, hy,
                    kernel = 'gaussian',
                    num_grid_points = 10**3):
    """
    Calculates a 2D kernel density estimate using FFTs with periodic
    boundary conditions. We simulate periodic boundary conditions by repeating
    the data in each dimension.

    Parameters
    ----------
    x : array
        1D array of x coordinates
    y : array
        1D array of y coordinates
    x_min : float
        Minimum x coordinate
    x_max : float
        Maximum x coordinate
    y_min : float
        Minimum y coordinate
    y_max : float
        Maximum y coordinate
    weights : array
        1D array of weights
    hx : float
        Bandwidth in metres
    hy : float
        Bandwidth in metres
    kernel : str
        KDE kernel
    num_grid_points : int
        Number of grid points

    Returns
    -------
    pdf_fun : array
        2D function  which interpolates the FFTKDE
    """

    # Extending coordiantes
    x_norm, y_norm, weights_norm = extend_coords(x, y, x_min, x_max, y_min,
                                                 y_max, weights)
    # Double number of grid points to account for extention
    num_grid_points *= 2

    # Stretching coordiantes
    x_norm = x_norm / hx
    y_norm = y_norm / hy

    # Calculating extended and stretched KDE
    fft_kde = FFTKDE(kernel=kernel, bw=1)
    coords_norm = np.vstack([x_norm, y_norm]).T
    grid_norm, pdf = \
        fft_kde.fit(coords_norm, weights = weights_norm).evaluate(num_grid_points)
    grid_norm_x, grid_norm_y = np.unique(grid_norm[:, 0]), np.unique(grid_norm[:, 1])
    pdf = pdf.reshape(num_grid_points, num_grid_points) \
        * np.sum(weights_norm) / np.sum(weights) / (hx * hy) # Renormalise

    # Evaluating extended and stretched KDE on an unstretched grid
    grid_x = grid_norm_x * hx
    grid_y = grid_norm_y * hy
    pdf_fun = interpolate.RectBivariateSpline(grid_x, grid_y, pdf,
                                              kx=3, ky=3, s=0.0)
    return pdf_fun

def calc_asymptotic_bias_1d(pdf_fun, grid_x, hx):
    """
    Calculates the asymptotic bias of the KDE estimate using a formula which was derived
    following the following reference:
    https://users.ssc.wisc.edu/~bhansen/718/NonParametrics1.pdf

    Parameters
    ----------
    pdf_fun : function
        1D function which interpolates the KDE
    grid_x : array
        1D array of grid points
    hx : float
        Bandwidth in metres

    Returns
    -------
    bias : array
        1D array of asymptotic bias
    """

    pdf_xx_array = pdf_fun(grid_x, nu=2)

    return 0.5 * pdf_xx_array * hx**2

def calc_asymptotic_variance_1d(pdf_fun, grid_x, hx, weights):
    """
    Calculates the asymptotic variance of the KDE estimate using a formula which was derived
    following the following reference:
    https://users.ssc.wisc.edu/~bhansen/718/NonParametrics1.pdf

    Parameters
    ----------
    pdf_fun : function
        1D function which interpolates the KDE
    grid_x : array
        1D array of grid points
    hx : float
        Bandwidth in metres
    weights : array
        1D array of weights

    Returns
    -------
    variance : array
        1D array of asymptotic variance
    """

    pdf_array = pdf_fun(grid_x)
    neff = np.sum(weights)**2 / np.sum(weights**2)

    return pdf_array / (np.sqrt(4 * np.pi) * neff * hx)
def calc_amse_1d(pdf_fun, grid_x, hx, weights):
    """
    Calculates the asymptotic mean square error of the KDE estimate using the formula
    amse = bias^2 + variance

    Parameters
    ----------
    pdf_fun : function
        1D function which interpolates the KDE
    grid_x : array
        1D array of grid points
    hx : float
        Bandwidth in metres
    weights : array
        1D array of weights

    Returns
    -------
    amse : array
        1D array of asymptotic mean square error
    """

    asymptotic_bias = calc_asymptotic_bias_1d(pdf_fun, grid_x, hx)
    asymptotic_variance = calc_asymptotic_variance_1d(pdf_fun, grid_x, hx, weights)

    return asymptotic_bias**2 + asymptotic_variance

def calc_amise_1d_array(x, grid_x, weights, hx_array,
                        num_grid_points = 10**3):
    """
    Calculate the asymptotic mean integrated square error for a set of bandwidths.
    Multiprocessing is used to speed up the calculation.

    Parameters
    ----------
    x : array
        1D array of data
    grid_x : array
        1D array of grid points
    weights : array
        1D array of weights
    hx_array : array
        1D array of bandwidths
    num_grid_points : int
        Number of grid points

    Returns
    -------
    amise_array : array
        1D array of asymptotic mean integrated square error
    """

    def calc_amise(n_array):

        for n in n_array:

            hx = hx_array[n]
            pdf_fun = periodic_kde_1d(x, x_min, x_max, weights, hx,
                                      num_grid_points=num_grid_points)
            amse_array = calc_amse_1d(pdf_fun, grid_x, hx, weights)
            amise = np.trapz(amse_array, x=grid_x)
            amise_array[n] = amise

    x_min = np.min(grid_x)
    x_max = np.max(grid_x)
    shm1 = shared_memory.SharedMemory(create = True,
                                      size = np.zeros(len(hx_array)).nbytes)
    amise_array = np.ndarray(len(hx_array), dtype = float, buffer = shm1.buf)
    n_array_list = np.array_split(np.arange(len(hx_array)), cpu_count())

    with Manager():
        processes = []

        for n_array in n_array_list:
            # pylint: disable=not-callable
            p = Process(target = calc_amise,
                args = (n_array,))
            # pylint: enable=not-callable
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

    amise_array = np.array(amise_array)
    shm1.close()
    shm1.unlink()

    return amise_array

def calc_asymptotic_bias_2d(pdf_fun, grid_x, grid_y, hx, hy):
    """
    Calculates the asymptotic bias of the KDE estimate using a formula which was derived
    following the following reference:
    https://users.ssc.wisc.edu/~bhansen/718/NonParametrics1.pdf

    Parameters
    ----------
    pdf_fun : function
        2D function which interpolates the KDE
    grid_x : array
        1D array of grid points
    grid_y : array
        1D array of grid points
    hx : float
        Bandwidth
    hy : float
        Bandwidth

    Returns
    -------
    bias : array
        2D array of asymptotic bias
    """
    pdf_xx_array = pdf_fun(grid_x, grid_y, dx=2)
    pdf_yy_array = pdf_fun(grid_x, grid_y, dy=2)
    # pdf_array = pdf_fun(grid_x, grid_y)
    # left = np.roll(pdf_array, 1, axis=0)
    # right = np.roll(pdf_array, -1, axis=0)
    # delta_x = grid_x[1] - grid_x[0]
    # pdf_xx_array = (right - 2 * pdf_array + left) / delta_x**2
    # up = np.roll(pdf_array, 1, axis=1)
    # down = np.roll(pdf_array, -1, axis=1)
    # delta_y = grid_y[1] - grid_y[0]
    # pdf_yy_array = (down - 2 * pdf_array + up) / delta_y**2

    return 0.5 * (pdf_xx_array * hx**2 + pdf_yy_array * hy**2)

def calc_asymptotic_variance_2d(pdf_fun, grid_x, grid_y, hx, hy, weights):
    """
    Calculates the asymptotic variance of the KDE estimate using a formula which was derived
    following the following reference:
    https://users.ssc.wisc.edu/~bhansen/718/NonParametrics1.pdf

    Parameters
    ----------
    pdf_fun : function
        2D function which interpolates the KDE
    grid_x : array
        1D array of grid points
    grid_y : array
        1D array of grid points
    hx : float
        Bandwidth
    hy : float
        Bandwidth
    weights : array
        1D array of weights

    Returns
    -------
    variance : array
    """

    pdf_array = pdf_fun(grid_x, grid_y)
    neff = np.sum(weights)**2 / np.sum(weights**2)

    return pdf_array / (4 * np.pi * neff * hx * hy)

def calc_amse_2d(pdf_fun, grid_x, grid_y, hx, hy, weights):
    """
    Calculates the asymptotic mean square error of the KDE estimate using the formula
    amse = bias^2 + variance

    Parameters
    ----------
    pdf_fun : function
        2D function which interpolates the KDE
    grid_x : array
        1D array of grid points
    grid_y : array
        1D array of grid points
    hx : float
        Bandwidth
    hy : float
        Bandwidth
    weights : array
        1D array of weights

    Returns
    -------
    amse : array
        2D array of asymptotic mean square error
    """

    asymptotic_bias = calc_asymptotic_bias_2d(pdf_fun, grid_x, grid_y, hx, hy)
    asymptotic_variance = calc_asymptotic_variance_2d(pdf_fun, grid_x, grid_y,
                                                      hx, hy, weights)

    return asymptotic_bias**2 + asymptotic_variance

def calc_amise_2d_array(x, y, grid_x, grid_y, weights,
                        hx_array, hy_array,
                        num_grid_points = 10**3):
    """
    Calculates the asymptotic mean integral square error for a set of bandwidths.
    Multiprocessing is used to speed up the calculation.

    Parameters
    ----------
    pdf_fun : function
        2D function which interpolates the KDE
    grid_x : array
        1D array of grid points
    grid_y : array
        1D array of grid points
    hx_array : array
        1D array of bandwidths
    hy_array : array
        1D array of bandwidths
    weights : array
        1D array of weights
    num_grid_points : int
        Number of grid points

    Returns
    -------
    amise : array
        2D array of the asymptotic mean integrated square error
        for each bandwidth combination
    """

    def calc_amise(n_array):

        for n in n_array:

            i = i_array[n]
            j = j_array[n]

            hx = hx_array[i]
            hy = hy_array[j]

            pdf_fun = periodic_kde_2d(x, y, x_min, x_max, y_min,
                                      y_max, weights, hx, hy,
                                      num_grid_points=num_grid_points)

            amse_array = calc_amse_2d(pdf_fun, grid_x, grid_y, hx, hy, weights)

            amise = np.trapz(np.trapz(amse_array, x=grid_x, axis=0), x=grid_y)
            amise_array[j, i] = amise

    x_min = np.min(grid_x)
    x_max = np.max(grid_x)
    y_min = np.min(grid_y)
    y_max = np.max(grid_y)

    shm1 = shared_memory.SharedMemory(create = True,
                                      size = np.zeros((len(hy_array), len(hx_array))).nbytes)
    amise_array = np.ndarray((len(hy_array), len(hx_array)), dtype = float, buffer = shm1.buf)
    i_array = np.zeros(len(hx_array) * len(hy_array), dtype = int)
    j_array = np.zeros(len(hx_array) * len(hy_array), dtype = int)
    count = 0
    for i in range(len(hx_array)):
        for j in range(len(hy_array)):
            i_array[count] = i
            j_array[count] = j
            count += 1
    n_array_list = np.array_split(np.arange(len(hx_array) * len(hy_array), dtype = int),
                                  cpu_count())

    with Manager():
        processes = []

        for n_array in n_array_list:
            # pylint: disable=not-callable
            p = Process(target = calc_amise,
                        args = (n_array,))
            # pylint: enable=not-callable
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

    amise_array = np.array(amise_array)
    shm1.close()
    shm1.unlink()

    return amise_array
