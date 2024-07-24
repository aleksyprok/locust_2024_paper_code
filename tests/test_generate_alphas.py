"""
Tests for the generate_alphas module.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from python_scripts import dt_fusion, generate_alphas, my_gfile_reader, prepare_profiles

def test_gen_marker_coords():
    """
    Test the gen_marker_coords function.

    This test case verifies that the gen_marker_coords function returns an array of marker
    coordinates with the correct shape.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gfile_path = os.path.join(repo_path, "input_data", 'SPR-045-16.eqdsk')
    num_markers = 1000
    gfile = my_gfile_reader.getGfile(gfile_path)
    marker_coords = generate_alphas.gen_marker_coords(num_markers, gfile)
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots()
    ax.scatter(marker_coords[0], marker_coords[2])
    ax.plot(gfile.R_bnd, gfile.Z_bnd, 'k')
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    ax.set_title('Initial marker coordinates')
    ax.set_aspect('equal')
    fig.savefig(output_dir + '/marker_coords_rz.png',
                bbox_inches='tight', dpi=300)

    fig, ax = plt.subplots()
    # make a scatter of the X, Y coordinates
    x_coords = marker_coords[0] * np.cos(marker_coords[1])
    y_coords = marker_coords[0] * np.sin(marker_coords[1])
    ax.scatter(x_coords, y_coords)
    # add circles with radius min(R_bnd) and max(R_bnd) both with centre at (0, 0)
    ax.add_artist(plt.Circle((0, 0), min(gfile.R_bnd), fill=False, color='k'))
    ax.add_artist(plt.Circle((0, 0), max(gfile.R_bnd), fill=False, color='k'))
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title('Initial marker coordinates')
    ax.set_aspect('equal')
    fig.savefig(output_dir + '/marker_coords_xy.png',
                bbox_inches='tight', dpi=300)

    assert marker_coords.shape == (3, num_markers)
    assert np.all(min(gfile.R_bnd) <= marker_coords[0])
    assert np.all(marker_coords[0] <= max(gfile.R_bnd))
    assert np.all(min(gfile.Z_bnd) <= marker_coords[2])
    assert np.all(marker_coords[2] <= max(gfile.Z_bnd))

def test_gen_marker_weights():
    """
    Test the gen_marker_weights function.

    This test case verifies that the gen_marker_weights function returns an array of marker
    weights with the correct shape.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gfile_path = os.path.join(repo_path, "input_data", 'SPR-045-16.eqdsk')
    input_dir = os.path.join(repo_path, "input_data")
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    os.makedirs(output_dir, exist_ok=True)
    spr_string = 'SPR-045-16'
    num_markers = 10
    gfile = my_gfile_reader.getGfile(gfile_path)
    cdf_filename = prepare_profiles.get_cdf_filename(spr_string, input_dir)
    psin, ti, _, _, nd, nt = prepare_profiles.read_cdf_file(cdf_filename)
    marker_coords = generate_alphas.gen_marker_coords(num_markers, gfile)
    weights = generate_alphas.gen_marker_weights(marker_coords, gfile, psin, nd, nt, ti)

    fig, ax = plt.subplots()
    ax.scatter(marker_coords[0], marker_coords[2])
    ax.plot(gfile.R_bnd, gfile.Z_bnd, 'k')
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    ax.set_title('Initial marker coordinates')
    ax.set_aspect('equal')
    for i, txt in enumerate(weights):
        ax.annotate(f'{txt:.2f}', (marker_coords[0, i], marker_coords[2, i]))
    fig.savefig(output_dir + '/marker_coords_rz_weights.png',
                bbox_inches='tight', dpi=300)

    assert weights.shape == (num_markers,)
    assert np.all(weights >= 0)

def test_generate_random_3d_vector():
    """
    Test the generate_random_3d_vector function.

    This test case verifies that the generate_random_3d_vector function returns an array of
    random 3D vectors with the correct shape and distribution.
    """

    num_samples = 1000
    v_alpha = np.ones(num_samples)
    vectors = generate_alphas.generate_random_3d_vector(v_alpha)

    # Make a histogram of the vectors
    fig, axs = plt.subplots(1, 3)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 3
    fig.set_size_inches(fig_size)
    for i, ax in enumerate(axs):
        ax.hist(vectors[i], bins=100, density=True)
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Probability')
    axs[0].set_title('X component')
    axs[1].set_title('Y component')
    axs[2].set_title('Z component')
    fig.savefig('tests/output_plots/random_3d_vectors.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

    assert vectors.shape == (3, num_samples)
    assert np.all(vectors >= -v_alpha)
    assert np.all(vectors <= v_alpha)
    # Check magnitude of vectors is correct
    magnitudes = np.sqrt(np.sum(vectors**2, axis=0))
    assert np.all(magnitudes >=  0)
    # They should all be equal to v_alpha within machine precision
    assert np.allclose(magnitudes, v_alpha, atol=1e-15)

def test_gen_marker_velocities():
    """
    Test the gen_marker_velocities function.

    This test case verifies that the gen_marker_velocities function returns an array of marker
    velocities with the correct shape and distribution.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    gfile_path = os.path.join(repo_path, "input_data", 'SPR-045-16.eqdsk')
    input_dir = os.path.join(repo_path, "input_data")
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    os.makedirs(output_dir, exist_ok=True)
    spr_string = 'SPR-045-16'
    num_markers = 1000
    gfile = my_gfile_reader.getGfile(gfile_path)
    cdf_filename = prepare_profiles.get_cdf_filename(spr_string, input_dir)
    psin, ti, _, _, _, _ = prepare_profiles.read_cdf_file(cdf_filename)
    marker_coords = generate_alphas.gen_marker_coords(num_markers, gfile)
    e_alpha_std_coords = dt_fusion.calc_e_alpha_std(marker_coords, gfile, ti, psin)
    marker_velocities = generate_alphas.gen_marker_velocities(e_alpha_std_coords)
    # mass of alpha particle = 6.644657230(82) × 10-27 Kg
    m_alpha = 6.644657230e-27
    marker_energies = 0.5 * m_alpha * np.sum(marker_velocities**2, axis=0)
    # Convert to MeV
    marker_energies /= 1.602e-13
    mean_e_alpha = 3.5 # MeV
    std_e_alpha = np.mean(e_alpha_std_coords) / 1.602e-13 # MeV

    fig, axs = plt.subplots(1, 3)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 3
    fig.set_size_inches(fig_size)
    for i, ax in enumerate(axs):
        ax.hist(marker_velocities[i], bins=100, density=True)
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Probability')
    axs[0].set_title('X component')
    axs[1].set_title('Y component')
    axs[2].set_title('Z component')
    fig.savefig(output_dir + '/marker_velocities.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(marker_energies, bins=100, density=True)
    xlim = ax.get_xlim()
    x = np.linspace(xlim[0], xlim[1], 100)
    y = (1 / (std_e_alpha * np.sqrt(2 * np.pi)) *
         np.exp(-0.5 * ((x - mean_e_alpha) / std_e_alpha)**2))
    ax.set_xlim(xlim)
    ax.plot(x, y, 'r')
    ax.set_xlabel('Energy [MeV]')
    ax.set_ylabel('Probability')
    ax.set_title('Alpha particle energy distribution')
    fig.savefig(output_dir + '/alpha_energy_distribution.png',
                bbox_inches='tight', dpi=300)
    plt.close(fig)

    assert marker_velocities.shape == (3, num_markers)
