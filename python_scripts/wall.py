"""
Module contains the routines to parameterize the wall.
"""

import time
# pylint: disable=no-name-in-module
from multiprocess import Process, Manager, cpu_count, shared_memory
# pylint: enable=no-name-in-module
import numpy as np
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import nearest_points
from scipy.spatial.distance import cdist
from scipy import interpolate

class Wall:
    """
    Class for storing information about the wall. Specifically, converting from
    cylindrical to coordinates at the wall to s_phi, s_theta wall coordiantes
    and vice versa.
    """
    def __init__(self,
                 wall_path: str,
                 special_nodes=(0, 1, 2, 3)):
        """
        Args:
            r: array of the r coordinates of the wall
            z_wall: array of the z coordinates of the wall
            special_nodes: list of the special node indices, we use this to
                           mark where boundaries of the outer wall
                           upper divertor, inner wall and lower divertor.
        """
        wall_rz = np.loadtxt(wall_path)
        self.r = wall_rz[:, 0]
        self.z = wall_rz[:, 1]
        self.special_nodes = special_nodes


        self.cr, self.cz = calc_wall_centroid(self.r, self.z)
        self.s_nodes = get_s_nodes(self.r, self.z)

        self.s_phi_min = 0
        self.s_phi_max = 2 * np.pi * self.cr

        self.s_theta_min = 0
        self.s_theta_max = self.s_nodes[-1]

        self.scale_factor_1d = ScaleFactor1D(self.r, self.z)
        self.scale_factor_2d = ScaleFactor2D(self.r, self.z)

def calc_wall_centroid(r_wall, z_wall):
    """
    calculate the centroid of the wall using formula given here:
    https://en.wikipedia.org/wiki/centroid#of_a_polygon
    """

    a = 0.5 * np.sum(r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1])

    cr = np.sum((r_wall[:-1] + r_wall[1:]) * (r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1]))
    cr /= 6 * a
    cz = np.sum((z_wall[:-1] + z_wall[1:]) * (r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1]))
    cz /= 6 * a

    return cr, cz

def get_s_theta_from_rz(r_coord, z_coord, r_wall, z_wall):
    """
    Calculate the s_theta coordinate from a set of 
    r and z coords which lie on the wall

    Args:
        r_coord: np.ndarray
            The r coordinates of the markers
        z_coord: np.ndarray
            The z coordinates of the markers
        r_wall: np.ndarray
            The r coordinates of the wall
        z_wall: np.ndarray
            The z coordinates of the wall
    """
    def calculate_s_theta(n_array):

        for n in n_array:

            # Calculate r_new, z_new which give the closest points on the wall
            point = Point(r_coord[n], z_coord[n])
            p1, _ = nearest_points(poly.boundary, point)
            r_new = p1.x
            z_new = p1.y

            # Determine the wall vertices either side of r_new, z_new
            nearest_node = np.argmin(cdist([[r_new, z_new]], wall_coords))
            if nearest_node == 0:
                coords0 = [r_wall[-2], z_wall[-2]]
            else:
                coords0 = [r_wall[nearest_node - 1], z_wall[nearest_node - 1]]
            coords1 = [r_wall[nearest_node], z_wall[nearest_node]]
            coords2 = [r_wall[nearest_node + 1], z_wall[nearest_node + 1]]

            # Check if the point is on line connecting coords0 with coords1
            if LineString((coords0, coords1)).distance(Point(r_new, z_new)) < 1e-8:
                if nearest_node == 0:
                    s_theta[n] = s_theta_nodes[-1] \
                            + cdist([coords0], [[r_new, z_new]])[0][0]
                else:
                    s_theta[n] = s_theta_nodes[nearest_node - 1] \
                            + cdist([coords0], [[r_new, z_new]])[0][0]

            # Check if the point is on line connecting coords1 with coords2
            elif LineString((coords1, coords2)).distance(Point(r_new, z_new)) < 1e-8:
                s_theta[n] = s_theta_nodes[nearest_node] \
                        + cdist([coords1], [[r_new, z_new]])[0][0]

            # In the unlikely event that the particle is not on the wall connected to the
            # nearest vertices then we check all walls. Note that this is more
            # computationally expensive.
            else:
                index = -1
                for coords1, coords2 in zip(wall_coords[:-1], wall_coords[1:]):
                    index += 1
                    if LineString((coords1, coords2)).distance(Point(r_new, z_new)) < 1e-8:
                        s_theta[n] = s_theta_nodes[index] \
                                + cdist([coords1], [[r_new, z_new]])[0][0]
                        # print('This point did not lie on walls connected to nearest vertex:', \
                        #       r_coord[n], z_coord[n])
                        break
                    if index == len(r_wall) - 2:
                        print('Error: Unable to calculate s_theta.')

    print('Calculating s_theta from r_coords and z_coords.')
    start_time = time.time()

    wall_coords = np.vstack([r_wall, z_wall]).T
    poly = Polygon(wall_coords)

    ds = np.sqrt(np.diff(r_wall)**2 + np.diff(z_wall)**2)
    s_theta_nodes = np.zeros(len(r_wall) - 1)
    s_theta_nodes[1:] = np.cumsum(ds[:-1])

    shm = shared_memory.SharedMemory(create = True, \
                                     size = np.zeros(len(r_coord)).nbytes)
    s_theta = np.ndarray(len(r_coord), dtype = float, buffer = shm.buf)
    n_array_list = np.array_split(np.arange(len(r_coord)), cpu_count())
    with Manager():
        processes = []

        for n_array in n_array_list:
            # pylint: disable=not-callable
            p = Process(target = calculate_s_theta, \
                args = (n_array,))
            # pylint: enable=not-callable
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

    s_theta = np.array(s_theta)
    shm.close()
    shm.unlink()

    print(f'Time taken: {time.time() - start_time:.2f} seconds')

    return s_theta

def get_s_phi_from_phi(phi, r_wall, z_wall):
    """
    Calculate the s_phi coordinate from the phi coordinate
    of the wall.

    Args:
        phi: np.ndarray
            The phi coordinate markers that hit the wall
        r_wall: np.ndarray
            The r coordinates of the wall
        z_wall: np.ndarray
            The z coordinates of the wall

    Returns:
        s_phi: np.ndarray
            The s_phi coordinates
    """
    cr, _ = calc_wall_centroid(r_wall, z_wall)
    return cr * (phi % (2 * np.pi))

def get_s_nodes(r_wall, z_wall):
    """
    Calculate the s coordinate from the r and z coordinates
    of the wall.

    Args:
        r_wall: np.ndarray
            The r coordinates of the wall
        z_wall: np.ndarray
            The z coordinates of the wall

    Returns:
        s_nodes: np.ndarray
            The s_theta coordinates of the wall
            nodes
    """

    dr = np.diff(r_wall)
    dz = np.diff(z_wall)
    ds = np.sqrt(dr**2 + dz**2)
    s_nodes = np.zeros(len(r_wall))
    s_nodes[1:] = np.cumsum(ds)

    return s_nodes

def get_rz_from_s_theta(s_theta, r_wall, z_wall):
    """
    Calculate the r and z coordinates from s_theta
    coordinates on the wall.

    Args:
        s_theta: np.ndarray
            The s_theta coordinates of the wall
        r_wall: np.ndarray
            The r coordinates of the wall
        z_wall: np.ndarray
            The z coordinates of the wall

    Returns:
        r_array: np.ndarray
            The r coordinates of the wall
        z_array: np.ndarray
            The z coordinates of the wall
    """

    dr = np.diff(r_wall)
    dz = np.diff(z_wall)
    ds = np.sqrt(dr**2 + dz**2)
    s_nodes = np.zeros(len(r_wall))
    s_nodes[1:] = np.cumsum(ds)

    # if s < 0 or s > s_max then impose periodic condition to get a value which
    # lies in the range [0, s_max].
    s_theta_norms = s_theta % s_nodes[-1]

    r_array = np.zeros(len(s_theta))
    z_array = np.zeros(len(s_theta))

    for i, s_theta_norm in enumerate(s_theta_norms):
        nearest_node_index_min = np.max(np.where(s_nodes <= s_theta_norm))
        index_float_part = (s_theta_norm - s_nodes[nearest_node_index_min]) \
                         / ds[nearest_node_index_min]
        r_array[i] = r_wall[nearest_node_index_min] + index_float_part * dr[nearest_node_index_min]
        z_array[i] = z_wall[nearest_node_index_min] + index_float_part * dz[nearest_node_index_min]

    return r_array, z_array

def get_rz_from_s_theta_funs(r_wall, z_wall):
    """
    Return the interpolating functions for the r and z coordinates
    """

    s_nodes = get_s_nodes(r_wall, z_wall)

    return interpolate.interp1d(s_nodes, r_wall, kind = 'linear'), \
           interpolate.interp1d(s_nodes, z_wall, kind = 'linear')

class ScaleFactor1D:
    """
    Class to calculate the scale factor as a function of s_theta.
    This can be thought of as a Jacobian determinant to give a measure
    of the size of an area element as function of s_theta.
    """

    def __init__(self, r_wall, z_wall):

        self.s_nodes = get_s_nodes(r_wall, z_wall)
        self.r_from_s_theta_fun, self.z_from_s_theta_fun = \
            get_rz_from_s_theta_funs(r_wall, z_wall)

    def __call__(self, s_theta):

        r_array = self.r_from_s_theta_fun(s_theta % self.s_nodes[-1])
        return 1 / (2 * np.pi * r_array)

class ScaleFactor2D:
    """
    Class to calculate the scale factor as a function of s_theta.
    This can be thought of as a Jacobian determinant to give a measure
    of the size of an area element as function of s_theta.
    """

    def __init__(self, r_wall, z_wall):

        self.s_nodes = get_s_nodes(r_wall, z_wall)
        self.r_from_s_theta_fun, self.z_from_s_theta_fun = \
            get_rz_from_s_theta_funs(r_wall, z_wall)
        self.cr, self.cz = calc_wall_centroid(r_wall, z_wall)

    def __call__(self, s_phi, s_theta):

        r_array = self.r_from_s_theta_fun(s_theta % self.s_nodes[-1])
        r_array = np.repeat(r_array[:, np.newaxis], len(s_phi), axis=1)
        return self.cr / r_array
