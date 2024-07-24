"""
This file contains the routines used to generate the initial position and velocity of the alpha
particles.

At the moment, this file contains only routines to generate markers uniformly across the plasma
volume and we have put the routines to generate markers from a KDE in a separate file, but the
plan is to combine the two in the future.
"""
import numpy as np
from matplotlib.path import Path
from python_scripts import dt_fusion

def gen_marker_coords(num_markers, gfile):
    """
    Generate the initial marker coordinates uniformly across the plasma volume.

    Args:
        num_markers: int
            The number of markers to generate.

    Returns:
        np.ndarray
            The marker positions with shape (3, num_markers).
            The 
    """

    r_coords = np.array([])
    z_coords = np.array([])
    r_min = min(gfile.R_bnd)
    r_max = max(gfile.R_bnd)
    z_min = min(gfile.Z_bnd)
    z_max = max(gfile.Z_bnd)
    while len(r_coords) < num_markers:
        r_coords = np.sqrt(np.random.uniform(low=r_min**2, high=r_max**2,
                                             size=num_markers * 2))
        z_coords = np.random.uniform(low=z_min, high=z_max,
                                     size=num_markers * 2)
        # remove markers outside LCFS
        polygon_path = Path(np.column_stack((gfile.R_bnd, gfile.Z_bnd)))
        coords = np.column_stack((r_coords, z_coords))
        inside_polygon = polygon_path.contains_points(coords)
        r_coords = r_coords[inside_polygon]
        z_coords = z_coords[inside_polygon]
        # psin_values = (gfile.PsiSpline(r_coords, z_coords, grid=False) - gfile.simag) \
        #             / (gfile.sibry - gfile.simag)
        # indices = np.where((0 <= psin_values) & (psin_values <= 1) &
        #                    (r_min <= r_coords) & (r_coords <= r_max) &
        #                    (z_min <= z_coords) & (z_coords <= z_max))[0]
        # r_coords = r_coords[indices]
        # z_coords = z_coords[indices]
    r_coords = r_coords[:num_markers]
    z_coords = z_coords[:num_markers]
    phi_sample = np.random.uniform(low=0, high=2 * np.pi, size=num_markers)
    coords = np.vstack([r_coords, phi_sample, z_coords])
    return coords

def gen_marker_weights(coords, gfile, psin, nd, nt, ti):
    """
    Generate the initial marker weights. We weight them accoding to the
    DT reaction rate in W/cm^3.

    Args:
        num_markers: int
            The number of markers to generate.

    Returns:
        np.ndarray
            The marker weights.
    """
    interp_nd = dt_fusion.ProfileCylindrical(gfile, nd, psin)
    interp_nt = dt_fusion.ProfileCylindrical(gfile, nt, psin)
    interp_ti = dt_fusion.ProfileCylindrical(gfile, ti, psin)

    ti_coords = interp_ti(coords[0], coords[2])
    ti_coords *= 1e-3  # convert from eV to keV
    reactivity_coords = dt_fusion.reactivity(ti_coords)
    reactivity_coords *= 1e-6  # convert from cm^3/s to m^3/s
    energy_per_fusion = 17.6e6 * 1.602e-19  # J
    reaction_rate_coords = reactivity_coords * interp_nd(coords[0], coords[2])
    reaction_rate_coords *= interp_nt(coords[0], coords[2]) * energy_per_fusion
    # convert from W/m^3 to MW/m^3
    reaction_rate_coords *= 1e-6

    return reaction_rate_coords

def generate_random_3d_vector(v_alpha):
    """
    Generate a random 3D vector with a magnitude of v_alpha.

    Args:
        v_alpha: 
            np.ndarray with shape (num_samples,)
    
    Returns:
        np.ndarray
            The random 3D vectors with shape (3, num_samples).
    """
    num_samples = v_alpha.size

    # Random azimuthal angle
    theta = np.random.uniform(0, 2 * np.pi, num_samples)

    # Random polar angle with correct distribution
    cos_phi = np.random.uniform(-1, 1, num_samples)
    phi = np.arccos(cos_phi)

    # Spherical to Cartesian conversion
    vx = v_alpha * np.sin(phi) * np.cos(theta)
    vy = v_alpha * np.sin(phi) * np.sin(theta)
    vz = v_alpha * np.cos(phi)

    return np.vstack([vx, vy, vz])

def gen_marker_velocities(e_std):
    """
    Generate the initial marker velocities.

    Args:
        e_std: np.ndarray
            The standard deviation of the alpha energy distribution
            in Joules.

    Returns:
        np.ndarray
            The marker velocities with shape (3, num_markers)
            in m/s.
    """
    num_markers = e_std.size
    # Alpha particle energy has mean 3.5 MeV
    mean_e_alpha = 3.5e6 * 1.602e-19  # J
    e_alpha = np.random.normal(loc=mean_e_alpha, scale=e_std, size=num_markers)
    # mass of alpha particle = 6.644657230(82) × 10-27 Kg
    m_alpha = 6.644657230e-27
    v_alpha = np.sqrt(2 * e_alpha / m_alpha)
    return generate_random_3d_vector(v_alpha)

def save_markers(output_path, coords, velocities, weights):
    """
    Save the markers to a file using np.savetxt.

    Args:
        output_path: str
            The path to the output file.
        coords: np.ndarray
            The marker positions with shape (3, num_markers).
        velocities: np.ndarray
            The marker velocities with shape (3, num_markers)
            in m/s.
        weights: np.ndarray
            The marker weights.
    """
    num_markers = coords.shape[1]
    data = np.zeros((num_markers, 7))
    data[:, :3] = coords.T
    data[:, 3:6] = velocities.T
    data[:, 6] = weights

    fmt = ['%13.6E'] * data.shape[1]

    header = '1.0\n1.0'

    np.savetxt(output_path, data, fmt=fmt, header=header, comments='')
