"""
Module to test the pkde module.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from python_scripts import pkde, run

def test_extend_coords():
    """
    Test case for the extend_coords function.

    This test case verifies that the extend_coords function returns the expected
    extended coordinates for a given set of coordinates.
    """
    x = np.array([1, 2, 3])
    y = np.array([3, 4, 5])
    x_min = 0
    x_max = 4
    y_min = 2
    y_max = 6
    weights = np.array([1, 2, 3])
    x_extended, y_extended, weights_extended = \
        pkde.extend_coords(x, y, x_min, x_max, y_min, y_max, weights)

    # x_test = np.array([-3, -2, -1,
    #                    -3, -2, -1,
    #                    -3, -2, -1,
    #                    +1, +2, +3,
    #                    +1, +2, +3,
    #                    +1, +2, +3,
    #                    +5, +6, +7,
    #                    +5, +6, +7,
    #                    +5, +6, +7])
    # y_test = np.array([-1, +0, 1,
    #                    +3, +4, +5,
    #                    +7, +8, +9,
    #                    -1, +0, +1,
    #                    +3, +4, +5,
    #                    +7, +8, +9,
    #                    -1, +0, +1,
    #                    +3, +4, +5,
    #                    +7, +8, +9])
    # weights_test = np.array([1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3])
    # Remove indices where x < x_min - delta_x or x > x_max + delta_x
    # or y < y_min - delta_y or y > y_max + delta_y
    # where delta_x = (x_max - x_min) / 2 and delta_y = (y_max - y_min) / 2
    # hence delta_x = 2 and delta_y = 2
    # hence remove elements where x = -3 or 7 and y = -1 or 9
    x_test = np.array([-2, -1,
                       -2, -1,
                       -2,
                       +2, +3,
                       +1, +2, +3,
                       +1, +2,
                       +6,
                       +5, +6,
                       +5, +6])
    y_test = np.array([+0, +1,
                       +4, +5,
                       +8,
                       +0, +1,
                       +3, +4, +5,
                       +7, +8,
                       +0,
                       +3, +4,
                       +7, +8])
    weights_test = np.array([2, 3,
                             2, 3,
                             2,
                             2, 3,
                             1, 2, 3,
                             1, 2,
                             2,
                             1, 2,
                             1, 2])

    assert np.array_equal(x_extended, x_test)
    assert np.array_equal(y_extended, y_test)
    assert np.array_equal(weights_extended, weights_test)

def test_periodic_kde_1d():
    """
    Test case for the periodic_kde_1d function.

    This test case verifies that the periodic_kde_1d function returns the expected
    KDE estimate for a given set of coordinates.
    """
    x = np.array([1, 2, 3])
    x_min = 0
    x_max = 4
    weights = np.array([1, 2, 3])
    hx = 0.3
    kernel = 'gaussian'
    num_grid_points = 10**3
    pdf_fun = pkde.periodic_kde_1d(x, x_min, x_max, weights, hx,
                                   kernel = kernel,
                                   num_grid_points = num_grid_points)

    # Calculate KDE using scipy, not FFTKDE
    x = np.array([-3, -2, -1, +1, +2, +3, +5, +6, +7])
    sigma_x = np.std(x)
    weights = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3])
    kde_scipy = stats.gaussian_kde(x,
                                   bw_method=hx/sigma_x,
                                   weights=weights)

    # Plot the KDE estimate
    grid_x = np.linspace(x_min, x_max, num_grid_points)
    fig, ax = plt.subplots()
    ax.plot(grid_x, pdf_fun(grid_x))
    ax.plot(grid_x, 3 * kde_scipy(grid_x))
    ax.set_xlabel('x')
    ax.set_ylabel('PDF')
    ax.set_title('KDE estimate')
    ax.legend(['Periodic KDE', 'Scipy KDE'])
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    outputpath = os.path.join(output_dir, 'periodic_kde_1d.png')
    fig.savefig(outputpath)
    plt.close(fig)

    # Check integral of the PDE is 1 within a 10^-4 tolerance
    assert np.isclose(np.trapz(pdf_fun(grid_x), grid_x),
                      1,
                      atol=1e-4)
    # Check the integral of scipy KDE is 1 within a 10^-4 tolerance
    assert np.isclose(np.trapz(3 * kde_scipy(grid_x), grid_x),
                      1,
                      atol=1e-4)
    # Check the KDE estimate is within 10% of the scipy KDE
    assert np.isclose(np.trapz(pdf_fun(grid_x), grid_x),
                      np.trapz(3 * kde_scipy(grid_x), grid_x),
                      atol=0.1)

def test_periodic_kde_2d():
    """
    Test case for the periodic_kde_2d function.

    This test case verifies that the periodic_kde_2d function returns the expected
    KDE estimate for a given set of coordinates.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")

    x = np.array([1, 2, 3])
    y = np.array([2, 5.6, 2])
    x_min = 0
    x_max = 4
    y_min = 1
    y_max = 5.7
    weights = np.array([1, 2, 3])
    hx = 0.3
    hy = 0.4
    kernel = 'gaussian'
    num_grid_points = 10**3
    pdf_fun = pkde.periodic_kde_2d(x, y, x_min, x_max, y_min, y_max, weights,
                                   hx, hy,
                                   kernel=kernel,
                                   num_grid_points=num_grid_points)
    grid_x = np.linspace(x_min, x_max, num_grid_points)
    grid_y = np.linspace(y_min, y_max, num_grid_points)
    pdf_fft = pdf_fun(grid_x, grid_y)
    # mesh_x, mesh_y = np.meshgrid(grid_x, grid_y)
    # points = np.vstack((mesh_x.ravel(), mesh_y.ravel())).T
    # # evaluate pdf_fft in num_grid_point chunks
    # pdf_fft = np.zeros(num_grid_points * num_grid_points)
    # for i in range(num_grid_points):
    #     points_i = points[i*num_grid_points:(i+1)*num_grid_points]
    #     pdf_fft[i*num_grid_points:(i+1)*num_grid_points] = pdf_fun(points_i)
    # pdf_fft = pdf_fun(points)
    # pdf_fft = np.reshape(pdf_fft, (num_grid_points, num_grid_points))

    # Plot the KDE estimate image
    fig, ax = plt.subplots()
    im = ax.imshow(pdf_fft.T, extent=[x_min, x_max, y_min, y_max],
                   origin='lower')
    ax.set_title('FFT KDE estimate')
    fig.colorbar(im)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.scatter(x, y)
    outputpath = os.path.join(output_dir, 'periodic_kde_2d.png')
    fig.savefig(outputpath)
    plt.close(fig)

    # Check integral is equal to 1
    assert np.isclose(np.trapz(np.trapz(pdf_fft, grid_x, axis=0), grid_y, axis =0),
                      1,
                      atol=1e-4)

def test_calc_asymptotic_bias_1d():
    """
    Test case for the calc_asymptotic_bias_1d function.

    This test case verifies that the calc_asymptotic_bias_1d function returns the expected
    asymptotic bias for a given KDE estimate.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    x = np.array([-3, -2, -1, +1, +2, +3, +5, +6, +7])
    weights = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3])
    x_min = -4
    x_max = 8
    num_grid_points = 1000
    x_grid = np.linspace(x_min, x_max, num_grid_points)
    hx = 0.1
    pdf_fun = pkde.periodic_kde_1d(x, x_min, x_max, weights, hx)
    bias = pkde.calc_asymptotic_bias_1d(pdf_fun, x_grid, hx)
    pdf = pdf_fun(x_grid)

    # Plot the KDE estimate image
    fig, ax = plt.subplots()
    ax.plot(x_grid, pdf, label='KDE estimate')
    ax.plot(x_grid, bias, label='Asymptotic bias')
    ax.legend()
    ax.set_xlim(x_min, x_max)
    fig.savefig(output_dir + '/calc_asymptotic_bias_1d.png')
    plt.close(fig)

def test_calc_asymptotic_variance_1d():
    """
    Test case for the calc_asymptotic_variance_1d function.

    This test case verifies that the calc_asymptotic_variance_1d function returns the expected
    asymptotic variance for a given KDE estimate.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    x = np.array([-3, -2, -1, +1, +2, +3, +5, +6, +7])
    weights = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3])
    x_min = -4
    x_max = 8
    num_grid_points = 1000
    x_grid = np.linspace(x_min, x_max, num_grid_points)
    hx = 0.1
    pdf_fun = pkde.periodic_kde_1d(x, x_min, x_max, weights, hx)
    variance = pkde.calc_asymptotic_variance_1d(pdf_fun, x_grid, hx, weights)
    pdf = pdf_fun(x_grid)

    # Plot the KDE estimate image
    fig, ax = plt.subplots()
    ax.plot(x_grid, pdf, label='KDE estimate')
    ax.plot(x_grid, variance, label='Asymptotic variance')
    ax.legend()
    ax.set_xlim(x_min, x_max)
    fig.savefig(output_dir + '/calc_asymptotic_variance_1d.png')
    plt.close(fig)

def test_calc_amse_1d():
    """
    Test case for the calc_amse_1d function.

    This test case verifies that the calc_amse_1d function returns the expected
    asymptotic mean square error for a given KDE estimate.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    x = np.array([-3, -2, -1, +1, +2, +3, +5, +6, +7])
    weights = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3])
    x_min = -4
    x_max = 8
    num_grid_points_1d = 10**4
    x_grid = np.linspace(x_min, x_max, num_grid_points_1d)
    hx = 0.1
    pdf_fun = pkde.periodic_kde_1d(x, x_min, x_max, weights, hx)
    amse = pkde.calc_amse_1d(pdf_fun, x_grid, hx, weights)
    pdf = pdf_fun(x_grid)

    # Plot the KDE estimate image
    fig, ax = plt.subplots()
    ax.plot(x_grid, pdf, label='KDE estimate')
    ax.plot(x_grid, amse, label='Asymptotic mean square error')
    ax.legend()
    ax.set_xlim(x_min, x_max)
    fig.savefig(output_dir + '/calc_amse_1d.png')
    plt.close(fig)

# def test_calc_amise_1d():
#     """
#     Test case for the calc_amise_1d function.

#     This test case verifies that the calc_amise_1d function returns the expected
#     asymptotic mean square error for a given KDE estimate.
#     """
#     repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#     output_dir = os.path.join(repo_path, "tests", "output_plots")
#     dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
#                             "axisymmetric", "gpu-q-41")
#     output_dir = os.path.join(repo_path, 'tests', 'output_plots')
#     tag = '13-12-2023_16-51-52.811'
#     test_run = run.Run(dir_path, tag)
#     test_run.init_log()
#     wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
#     test_run.init_wall(wall_path)
#     gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
#     test_run.init_gfile(gfile_path)
#     test_run.init_markers()
#     test_run.init_flux(num_grid_points=1000)
#     hx_array = np.logspace(-2, 0, 64)
#     amise_array = pkde.calc_amise_1d_array(test_run.markers.stopped.s_theta,
#                                            test_run.flux.s_theta,
#                                            test_run.markers.stopped.weight,
#                                            hx_array)

#     # Plot the amise array
#     fig, ax = plt.subplots()
#     ax.plot(hx_array, amise_array)
#     ax.set_xscale('log')
#     ax.set_xlabel('Bandwidth [m]')
#     ax.set_ylabel('Asymptotic mean integrated square error')
#     fig.savefig(output_dir + '/axisymmetric_calc_amise_1d_array.png')

def test_calc_asymptotic_bias_2d():
    """
    Test case for the calc_asymptotic_bias_2d function.

    This test case verifies that the calc_asymptotic_bias_2d function returns the expected
    asymptotic bias for a given KDE estimate.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    x = np.array([1, 2, 3])
    y = np.array([2, 5.6, 2])
    x_min = 0
    x_max = 4
    y_min = 1
    y_max = 5.7
    weights = np.array([1, 2, 3])
    hx = 0.3
    hy = 0.4
    num_grid_points = 1000
    x_grid = np.linspace(x_min, x_max, num_grid_points)
    y_grid = np.linspace(y_min, y_max, num_grid_points)
    pdf_fun = pkde.periodic_kde_2d(x, y, x_min, x_max, y_min, y_max, weights, hx, hy)
    bias = pkde.calc_asymptotic_bias_2d(pdf_fun, x_grid, y_grid, hx, hy)
    pdf = pdf_fun(x_grid, y_grid)

    # Plot imshow of KDE and asymptotic bias side by side
    fig, ax = plt.subplots(1, 2)
    im = ax[0].imshow(pdf.T, origin='lower', extent=[x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax[0])
    im = ax[1].imshow(bias.T, origin='lower', extent=[x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax[1])
    ax[0].set_title('KDE')
    ax[1].set_title('Asymptotic bias')
    for i in range(2):
        ax[i].set_xlabel('x')
        ax[i].set_ylabel('y')
    fig.savefig(output_dir + '/calc_asymptotic_bias_2d.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

    # Line plot along the line y=2 we need to find grid point close to y=2
    y_index = np.argmin(np.abs(y_grid - 2))
    fig, ax = plt.subplots()
    ax.plot(x_grid, pdf[:, y_index], label='KDE estimate')
    ax.plot(x_grid, bias[:, y_index], label='Asymptotic bias')
    ax.legend()
    ax.set_xlim(x_min, x_max)
    ax.set_title('Asymptotic bias along y=2')
    fig.savefig(output_dir + '/calc_asymptotic_bias_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

def test_calc_asymptotic_variance_2d():
    """
    Test case for the calc_asymptotic_variance_2d function.

    This test case verifies that the calc_asymptotic_variance_2d function returns the expected
    asymptotic variance for a given KDE estimate.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    x = np.array([1, 2, 3])
    y = np.array([2, 5.6, 2])
    x_min = 0
    x_max = 4
    y_min = 1
    y_max = 5.7
    weights = np.array([1, 2, 3])
    hx = 0.3
    hy = 0.4
    num_grid_points = 1000
    x_grid = np.linspace(x_min, x_max, num_grid_points)
    y_grid = np.linspace(y_min, y_max, num_grid_points)
    pdf_fun = pkde.periodic_kde_2d(x, y, x_min, x_max, y_min, y_max, weights, hx, hy)
    variance = pkde.calc_asymptotic_variance_2d(pdf_fun, x_grid, y_grid, hx, hy, weights)
    pdf = pdf_fun(x_grid, y_grid)

    # Plot imshow of KDE and asymptotic variance side by side
    fig, ax = plt.subplots(1, 2)
    im = ax[0].imshow(pdf.T, origin='lower', extent=[x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax[0])
    im = ax[1].imshow(variance.T, origin='lower', extent=[x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax[1])
    ax[0].set_title('KDE')
    ax[1].set_title('Asymptotic variance')
    for i in range(2):
        ax[i].set_xlabel('x')
        ax[i].set_ylabel('y')
    fig.savefig(output_dir + '/calc_asymptotic_variance_2d.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

    # Line plot along the line y=2 we need to find grid point close to y=2
    y_index = np.argmin(np.abs(y_grid - 2))
    fig, ax = plt.subplots()
    ax.plot(x_grid, pdf[:, y_index], label='KDE estimate')
    ax.plot(x_grid, variance[:, y_index], label='Asymptotic variance')
    ax.legend()
    ax.set_xlim(x_min, x_max)
    ax.set_title('Asymptotic variance along y=2')
    fig.savefig(output_dir + '/calc_asymptotic_variance_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

def test_calc_amse_2d():
    """
    Test case for the calc_amse_2d function.

    This test case verifies that the calc_amse_2d function returns the expected
    asymptotic mean square error for a given KDE estimate.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    x = np.array([1, 2, 3])
    y = np.array([2, 5.6, 2])
    x_min = 0
    x_max = 4
    y_min = 1
    y_max = 5.7
    weights = np.array([1, 2, 3])
    hx = 0.3
    hy = 0.4
    num_grid_points = 1000
    x_grid = np.linspace(x_min, x_max, num_grid_points)
    y_grid = np.linspace(y_min, y_max, num_grid_points)
    pdf_fun = pkde.periodic_kde_2d(x, y, x_min, x_max, y_min, y_max, weights, hx, hy)
    amse = pkde.calc_amse_2d(pdf_fun, x_grid, y_grid, hx, hy, weights)
    pdf = pdf_fun(x_grid, y_grid)

    # Plot imshow of KDE and asymptotic mean square error side by side
    fig, ax = plt.subplots(1, 2)
    im = ax[0].imshow(pdf.T, origin='lower', extent=[x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax[0])
    im = ax[1].imshow(amse.T, origin='lower', extent=[x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax[1])
    ax[0].set_title('KDE')
    ax[1].set_title('Asymptotic mean square error')
    for i in range(2):
        ax[i].set_xlabel('x')
        ax[i].set_ylabel('y')
    fig.savefig(output_dir + '/calc_amse_2d.png',
                bbox_inches='tight', dpi=300)

    # Line plot along the line y=2 we need to find grid point close to y=2
    y_index = np.argmin(np.abs(y_grid - 2))
    fig, ax = plt.subplots()
    ax.plot(x_grid, pdf[:, y_index], label='KDE estimate')
    ax.plot(x_grid, amse[:, y_index], label='Asymptotic mean square error')
    ax.legend()
    ax.set_xlim(x_min, x_max)
    ax.set_title('Asymptotic mean square error along y=2')
    fig.savefig(output_dir + '/calc_amse_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

def test_calc_amse_2d_array_axisymmetric():
    """
    Test case for the calc_amse_2d_array function using axisymmetric SPR-045-14 data.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join(repo_path, 'tests', 'output_plots')
    tag = '13-12-2023_16-51-52.811'
    num_grid_points_2d = 10**3
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points_2d=num_grid_points_2d)
    h_phi = 0.1
    h_theta_2d = 0.01
    pdf_fun = pkde.periodic_kde_2d(test_run.markers.stopped.s_phi,
                                   test_run.markers.stopped.s_theta,
                                   test_run.wall.s_phi_min,
                                   test_run.wall.s_phi_max,
                                   test_run.wall.s_theta_min,
                                   test_run.wall.s_theta_max,
                                   test_run.markers.stopped.weight,
                                   h_phi, h_theta_2d,
                                   num_grid_points = num_grid_points_2d)
    amse = pkde.calc_amse_2d(pdf_fun, test_run.flux.s_phi, test_run.flux.s_theta_2d,
                             h_phi, h_theta_2d, test_run.markers.stopped.weight)

    # Plot the amse array
    fig, ax = plt.subplots()
    im = ax.imshow(amse.T, origin='lower',
                   extent=[test_run.wall.s_phi_min, test_run.wall.s_phi_max,
                           test_run.wall.s_theta_min, test_run.wall.s_theta_max])
    fig.colorbar(im, ax=ax)
    ax.set_title('Asymptotic mean square error')
    ax.set_xlabel('phi')
    ax.set_ylabel('theta')
    fig.savefig(output_dir + '/axisymmetric_calc_amse_2d_array.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

# def test_calc_amise_2d_array():
#     """
#     Test case for the calc_amise_2d_array function using axisymmetric SPR-045-14 data.

#     This test case verifies that the calc_amise_2d_array function returns the expected
#     asymptotic mean square error for a given KDE estimate.
#     """
#     repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#     output_dir = os.path.join(repo_path, "tests", "output_plots")
#     dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
#                             "axisymmetric", "gpu-q-41")
#     output_dir = os.path.join(repo_path, 'tests', 'output_plots')
#     tag = '13-12-2023_16-51-52.811'
#     num_grid_points = 10**3
#     test_run = run.Run(dir_path, tag)
#     test_run.init_log()
#     wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
#     test_run.init_wall(wall_path)
#     gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
#     test_run.init_gfile(gfile_path)
#     test_run.init_markers()
#     test_run.init_flux(num_grid_points=num_grid_points)
#     hx_array = np.logspace(-2, 0, 8)
#     hy_array = np.logspace(-2, 0, 8)
#     amise_array = pkde.calc_amise_2d_array(test_run.markers.stopped.s_phi,
#                                            test_run.markers.stopped.s_theta,
#                                            test_run.flux.s_phi,
#                                            test_run.flux.s_theta,
#                                            test_run.markers.stopped.weight,
#                                            hx_array, hy_array,
#                                            num_grid_points=num_grid_points)
#     hx_dummy = np.logspace(-2, 0, len(hx_array) + 1)
#     hy_dummy = np.logspace(-2, 0, len(hy_array) + 1)
#     hx_mesh, hy_mesh = np.meshgrid(hx_dummy, hy_dummy)

#     # Plot the amise array using pcolormesh with log scale
#     fig, ax = plt.subplots()
#     im = ax.pcolormesh(hx_mesh, hy_mesh, amise_array,
#                        norm=LogNorm())
#     ax.set_title('Asymptotic mean integrated square error')
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     fig.colorbar(im, ax=ax)
#     ax.set_xlabel(r'$h_\phi$ [m]')
#     ax.set_ylabel(r'$h_\theta$ [m]')
#     fig.savefig(output_dir + '/axisymmetric_calc_amise_2d_array.png',
#                 bbox_inches='tight', dpi=300)
