"""
Module for testing the wall module.
"""
import numpy as np
from python_scripts import wall

def test_calc_wall_centroid():
    """
    Test the calc_wall_centroid function.
    """
    r_wall = np.array([0, 1, 1, 0, 0])
    z_wall = np.array([0, 0, 1, 1, 0])
    assert wall.calc_wall_centroid(r_wall, z_wall) == (0.5, 0.5)

def test_get_s_theta_from_rz():
    """
    Test the get_s_theta_from_rz function.
    """
    r_coord = np.array([0.300, 1.001, 0.700, 0.001])
    z_coord = np.array([0.001, 0.200, 0.999, 0.200])
    r_wall = np.array([0, 1, 1, 0, 0])
    z_wall = np.array([0, 0, 1, 1, 0])
    s_theta = np.array([0, 0, 0, 0, 0])
    s_theta = wall.get_s_theta_from_rz(r_coord, z_coord, r_wall, z_wall)
    s_theta_test = np.array([0.3, 1.2, 2.3, 3.8])
    assert np.allclose(s_theta, s_theta_test, atol=1e-8)

def test_get_s_nodes():
    """
    Test the get_s_nodes function.
    """
    r_wall = np.array([0, 1, 1, 0, 0])
    z_wall = np.array([0, 0, 1, 1, 0])
    s_nodes = np.array([0, 1, 2, 3, 4])
    assert np.allclose(wall.get_s_nodes(r_wall, z_wall), s_nodes, atol=1e-8)

def test_get_rz_from_s_theta():
    """
    Test the get_rz_from_s_theta function.
    """
    s_theta = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5])
    r_wall = np.array([0, 1, 1, 0, 0])
    z_wall = np.array([0, 0, 1, 1, 0])
    r_array, z_array = wall.get_rz_from_s_theta(s_theta, r_wall, z_wall)
    r_array_test = np.array([0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0])
    z_array_test = np.array([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.5])
    assert np.allclose(r_array, r_array_test, atol=1e-8)
    assert np.allclose(z_array, z_array_test, atol=1e-8)

def test_get_rz_from_s_theta_funs():
    """
    Test the get_rz_from_s_theta_funs function.
    """
    r_wall = np.array([0, 1, 1, 0, 0])
    z_wall = np.array([0, 0, 1, 1, 0])
    r_fun, z_fun = wall.get_rz_from_s_theta_funs(r_wall, z_wall)
    s_theta_array = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5])
    r_array_test = np.array([0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0])
    z_array_test = np.array([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.5])
    r_array = r_fun(s_theta_array)
    z_array = z_fun(s_theta_array)
    assert np.allclose(r_array, r_array_test, atol=1e-8)
    assert np.allclose(z_array, z_array_test, atol=1e-8)

def test_scalefactor1d():
    """
    Test the ScaleFactor1D class.
    """
    r_wall = np.array([1, 2, 2, 1, 1])
    z_wall = np.array([0, 0, 1, 1, 0])
    s_theta = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5])
    r_array = np.array([1, 1.5, 2, 2, 2, 1.5, 1, 1])
    sf = wall.ScaleFactor1D(r_wall, z_wall)
    sf_array = sf(s_theta)
    sf_array_test = 1 / (2 * np.pi * r_array)
    assert np.allclose(sf_array, sf_array_test, atol=1e-8)

def test_get_s_phi_from_phi():
    """
    Test the get_s_phi_from_phi function.
    """
    phi = np.array([-2, 2, 4, 6, 2 * np.pi + 2])
    r_wall = np.array([1, 2, 2, 1, 1])
    z_wall = np.array([0, 0, 1, 1, 0])
    cr = 1.5
    s_phi = wall.get_s_phi_from_phi(phi, r_wall, z_wall)
    s_phi_test = cr * (phi % (2 * np.pi))
    assert np.allclose(s_phi, s_phi_test, atol=1e-8)
