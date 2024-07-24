"""
Calculate the total fusion power in the plasma using a .CDF file

We use a monte carlo integration technique.

We sample  markers uniformly from an annulus volume
with inner radius R_MIN, outer radius R_MAX, and height Z_MAX - Z_MIN.

We then calculate the ion temperature at each marker using a
cylindrical profile interpolation.

We then calculate the DT reactivity at each marker using the ion
temperature.

We then calculate the fusion power at each marker using the DT
reactivity.

We then calculate the total fusion power by summing the fusion power
at each marker and dividing by the number of markers.
"""

import time
import numpy as np
from scipy import interpolate

class ProfileCylindrical:
    """
    Profile class to provide a way to calculate the ion temperature
    or density from R, Z coordinates.
    """
    def __init__(self, input_gfile, input_profile, input_psin):
        self.gfile = input_gfile
        self.interp_profile = interpolate.interp1d(input_psin, input_profile,
                                                   bounds_error=False,
                                                   fill_value=1e-10)

    def __call__(self, r_cyl, z_cyl):
        """
        Calculate the ion temperature or density at a given R, Z coordinate
        """
        psi_at_rz = self.gfile.PsiSpline(r_cyl, z_cyl, grid=False)
        psin_at_rz = (psi_at_rz - self.gfile.simag) / (self.gfile.sibry - self.gfile.simag)
        return self.interp_profile(psin_at_rz)

def reactivity(temp_ions):
    """
    This function computes the DT reactivity. Using this paper:
    H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611

    Input:
      - temp_ions: temperature of the ions in keV

    Output:
      - reactivity: DT reactivity in cm^3/s
    """

    gamov_constant = 34.382  # in sqrt(kev)
    mrc2 = 1124656  # in keV
    const_1 = 1.17302e-9
    const_2 = 1.51361e-2
    const_3 = 7.51886e-2
    const_4 = 4.60643e-3
    const_5 = 1.35000e-2
    const_6 = -1.06750e-4
    const_7 = 1.36600e-5

    # Calculate theta
    denominator = -temp_ions * (const_2 +
                                temp_ions * (const_4 +
                                             temp_ions * const_6))
    denominator /= 1 + temp_ions * (const_3 +
                                    temp_ions * (const_5 +
                                                 temp_ions * const_7))
    denominator += 1
    theta = temp_ions / denominator
    xi_greek = (0.25 * gamov_constant**2 / theta)**(1/3)
    sigma_v = const_1 * theta * np.sqrt(xi_greek / (mrc2 * temp_ions**3))
    sigma_v *= np.exp(-3 * xi_greek)
    return sigma_v

def calc_fusion_power(psin, nd, nt, ti, gfile):
    """
    Calculate the fusion power using a monte carlo integration technique.
    """

    start_time = time.time()

    r_min = np.min(gfile.R_bnd)
    r_max = np.max(gfile.R_bnd)
    z_min = np.min(gfile.Z_bnd)
    z_max = np.max(gfile.Z_bnd)
    num_markers = int(1e6)
    r_coords = np.sqrt(np.random.uniform(low=r_min**2, high=r_max**2,
                                        size=num_markers))
    z_coords = np.random.uniform(low=z_min, high=z_max,
                                size=num_markers)

    interp_nd = ProfileCylindrical(gfile, nd, psin)
    interp_nt = ProfileCylindrical(gfile, nt, psin)
    interp_ti = ProfileCylindrical(gfile, ti, psin)

    annulus_area = np.pi * (r_max**2 - r_min**2)
    annulus_area_times_height = annulus_area * (z_max - z_min)

    ti_coords = interp_ti(r_coords, z_coords)
    ti_coords *= 1e-3  # convert from eV to keV
    reactivity_coords = reactivity(ti_coords)
    reactivity_coords *= 1e-6  # convert from cm^3/s to m^3/s
    energy_per_fusion = 17.6e6 * 1.602e-19  # J
    reaction_rate = reactivity_coords * interp_nd(r_coords, z_coords)
    reaction_rate *= interp_nt(r_coords, z_coords) * energy_per_fusion
    fusion_power = np.sum(reaction_rate) * annulus_area_times_height / num_markers
    alpha_power = fusion_power * 3.5 / 17.6

    end_time = time.time()
    print(f'Time taken to calculate fusion power: {end_time - start_time:.2f} s')

    return fusion_power, alpha_power

def calc_e_alpha_std(coords, gfile, ti, psin):
    """
    Calculate the standard deviation of the alpha energy distribution
    using Eq 36 of H. Brysk, Plasma Physics 15, 611 (1973)

    Args:
        coords: np.ndarray
            The marker positions with shape (3, num_markers).
        gfile: gfile object
            The gfile object.
        ti: np.ndarray
            The ion temperature profile in eV
        psin: np.ndarray
            The normalised poloidal flux profile.
    
    Returns:
        np.ndarray
            The standard deviation of the alpha energy distribution
            with shape (num_markers,) in Joules.
    """

    interp_ti = ProfileCylindrical(gfile, ti, psin)
    ti_coords = interp_ti(coords[0], coords[2])
    # Convert from eV to J
    ti_coords *= 1.602e-19
    # e_alpha = 3.5 MeV
    # Convert from MeV to J
    e_alpha = 3.5 * 1e6 * 1.602e-19
    # mass of alpha particle = 6.644657230(82) × 10-27 Kg
    m_alpha = 6.644657230e-27
    # mass of a neutron = 6.644657230(82) × 10-27 Kg
    m_n = 1.67492749804e-27
    variance_coords = 2 * m_alpha * ti_coords * e_alpha \
                    / (m_alpha + m_n)
    return np.sqrt(variance_coords)
