"""
This file check the ripple_check.py module.
"""
import numpy as np
import matplotlib.pyplot as plt
from python_scripts import ripple_check

def test_field_vwire():
    """
    Test case for the field_vwire function.

    We plot the magnetic field to check it's as expected
    """

    z_min = -1 # m
    z_max = 1 # m
    current = 1e6 # Amps
    y_vals = [-1, 0, 1]
    z_vals = [1, 0, -1]
    labels = [r"$B_x$", r"$B_y$", r"$B_z$"]
    fig, axs = plt.subplots(3, 3)
    fig_size = fig.get_size_inches()
    fig_size *= 3
    fig.set_size_inches(fig_size)
    for i, z_val in enumerate(z_vals):
        for j, y_val in enumerate(y_vals):
            x = np.concatenate((np.linspace(-5, -1, 100), np.linspace(1, 5, 100)))
            y = np.zeros_like(x) + y_val
            z = np.zeros_like(x) + z_val
            coords = [x, y, z]
            b_field = ripple_check.field_vwire(coords,
                                            z_min=z_min,
                                            z_max=z_max,
                                            current=current)
            for k, label in enumerate(labels):
                axs[i, j].plot(x, b_field[k],
                               label=label)
            axs[i, j].set_xlabel("x [m]")
            axs[i, j].set_ylabel("Magnetic field [T]")
            axs[i, j].set_title(f"y = {y_val}, z = {z_val} m")
            axs[i, j].set_ylim(-0.15, 0.15)
            axs[i, j].legend()

    fig.savefig('tests/output_plots/field_vwire.png')

def test_rotate_vector():
    """
    Test case for the rotate_vector function.

    We plot the magnetic field to check it's as expected
    """

    vector = [0, 0, 1]
    phi = np.pi/6
    new_vector = ripple_check.rotate_vector(vector, phi)
    new_vector_check = [np.sqrt(3)/2, 0.5, 0]
    assert np.allclose(new_vector, new_vector_check)

    vector = [1/np.sqrt(2), 1/np.sqrt(2), 0]
    phi = np.pi/4
    new_vector = ripple_check.rotate_vector(vector, phi,
                                            inverse=True)
    new_vector_check = [0, 0, 1]
    assert np.allclose(new_vector, new_vector_check)

def test_rotate_vector_plot():
    """
    Test case for the rotate_vector function.

    We plot the magnetic field to check it's as expected
    """

    z_min = 0 # m
    z_max = 2 # m
    current = 1e6 # Amps
    phi = np.pi/6
    wire_z_coord = 5
    new_y_plot_vals = [-1, 0, 1]
    new_z_plot_vals = [2, 1, 0]
    labels = [r"$B_x$", r"$B_y$", r"$B_z$"]
    fig, axs = plt.subplots(3, 3)
    fig_size = fig.get_size_inches()
    fig_size *= 3
    fig.set_size_inches(fig_size)
    for i, new_z_plot_val in enumerate(new_z_plot_vals):
        for j,  new_y_plot_val in enumerate(new_y_plot_vals):
            z = np.concatenate((np.linspace(wire_z_coord - 5, wire_z_coord - 1, 100),
                                np.linspace(wire_z_coord + 1, wire_z_coord + 5, 100)))
            x = np.zeros_like(z)
            y = np.zeros_like(z)
            coords = [x, y, z]
            new_coords = ripple_check.rotate_vector(coords, phi,
                                                    inverse=True)
            new_coords[0] += wire_z_coord
            new_coords[1] += new_y_plot_val
            new_coords[2] += new_z_plot_val
            new_plot_vals = [0, new_y_plot_val, new_z_plot_val]
            plot_vals = ripple_check.rotate_vector(new_plot_vals, phi,
                                                   inverse=False)
            b_field =  ripple_check.field_vwire(new_coords,
                                                z_min=z_min,
                                                z_max=z_max,
                                                current=current)
            b_field = ripple_check.rotate_vector(b_field, phi)
            for k, label in enumerate(labels):
                axs[i, j].plot(z, b_field[k],
                               label=label)
            if i == 2:
                axs[i, j].set_xlabel("z [m]")
            axs[i, j].set_ylabel("Magnetic field [T]")
            axs[i, j].set_title(f"new_y = {new_y_plot_val}, new_z = {new_z_plot_val}, "
                                "\n"
                                f"x = {plot_vals[0]:.2f}, "
                                f"y = {plot_vals[1]:.2f}, "
                                f"z = {plot_vals[2]:.2f}")
            axs[i, j].set_ylim(-0.15, 0.15)
            axs[i, j].legend()

    fig.suptitle(f"phi = {phi}")
    fig.savefig('tests/output_plots/rotate_vector.png')

def test_field_coil():
    """
    Test case for the field_coil function.

    We check that Ampere's law is satisfied
    """

    height = 2
    current = 1e6
    phi = 0
    r_inner = 4
    r_outer = 8
    mu0 = 4*np.pi*1e-7
    n_coords = int(1e4)

    r0 = 7
    phi_coords = np.linspace(0, 2*np.pi, n_coords)
    r_coords = np.zeros_like(phi_coords) + r0
    z_coords = np.zeros_like(phi_coords)
    coords = [r_coords * np.cos(phi_coords),
              r_coords * np.sin(phi_coords),
              z_coords]

    b_field = ripple_check.field_coil(coords, phi,
                                      r_inner=r_inner,
                                      r_outer=r_outer,
                                      height=height,
                                      current=current)

    # Check that Ampere's law is satisfied
    # By integrating the magnetic field in phi
    b_phi = np.cos(phi_coords) * b_field[1] \
          - np.sin(phi_coords) * b_field[0]
    b_phi_int = r0 * np.trapz(b_phi, phi_coords)

    print("b_phi_int = ", b_phi_int)
    print("mu0 * current = ", mu0 * current)
    print("% error = ", 100*(b_phi_int - mu0*current)/mu0/current)
    assert np.isclose(b_phi_int, mu0*current)

def test_total_field():
    """
    Test case for the total_field function.

    We check that Ampere's law is satisfied
    """

    height = 2
    current = 1e6
    r_inner = 4
    r_outer = 8
    mu0 = 4*np.pi*1e-7
    n_coords = int(1e4)
    num_coils = 16

    r0 = 7
    phi_coords = np.linspace(0, 2*np.pi, n_coords)
    r_coords = np.zeros_like(phi_coords) + r0
    z_coords = np.zeros_like(phi_coords)
    coords = [r_coords * np.cos(phi_coords),
              r_coords * np.sin(phi_coords),
              z_coords]

    b_field = ripple_check.total_field(coords, num_coils,
                                       r_inner=r_inner,
                                       r_outer=r_outer,
                                       height=height,
                                       current=current)

    # Check that Ampere's law is satisfied
    # By integrating the magnetic field in phi
    b_phi = np.cos(phi_coords) * b_field[1] \
          - np.sin(phi_coords) * b_field[0]
    b_phi_int = r0 * np.trapz(b_phi, phi_coords)

    assert np.isclose(b_phi_int, mu0*current*num_coils)

if __name__ == "__main__":
    test_field_coil()
