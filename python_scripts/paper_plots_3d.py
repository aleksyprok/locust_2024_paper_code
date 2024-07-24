"""
This file contains functions used to make the 3D plots that appears in the IAEA FEC 2024.
"""

import os
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from python_scripts import my_gfile_reader

def coil_plot_3d(gfile_path=None):
    """
    Plot the 3D coil plot in the paper.
    """

    gfile = my_gfile_reader.getGfile(gfile_path)

    # Create a mesh of points for revolution
    phi_plasma = np.linspace(0, 2 * np.pi, 100)  # angles for revolution
    r_mg, phi_mg = np.meshgrid(gfile.R_bnd, phi_plasma)
    x_mg = r_mg * np.cos(phi_mg)
    y_mg = r_mg * np.sin(phi_mg)
    r_mg = np.outer(np.ones_like(phi_plasma), gfile.Z_bnd)

    # Create the 3D plot
    fig = plt.figure()
    fig_size = fig.get_size_inches()
    fig.set_size_inches((fig_size[0] * 2, fig_size[1] * 2))
    ax = fig.add_subplot(111, projection='3d')

    color_array = np.full_like(r_mg, fill_value='tab:orange', dtype=np.dtype('U10'))
    ax.plot_surface(x_mg, y_mg, r_mg, facecolors=color_array, edgecolor='none', alpha=0.2)

    # EFCC

    num_elm_coils = 16

    elm_coil_r00 = 8.0
    elm_coil_r01 = 8.0
    elm_coil_z00 = 2.5
    elm_coil_z01 = 5.5
    r_elm_0 = [elm_coil_r00, elm_coil_r01]
    z_elm_0 = [elm_coil_z00, elm_coil_z01]

    elm_coil_r10 = 8.0
    elm_coil_r11 = 8.0
    elm_coil_z10 = -2.5
    elm_coil_z11 = -5.5
    r_elm_1 = [elm_coil_r10, elm_coil_r11]
    z_elm_1 = [elm_coil_z10, elm_coil_z11]

    arc_length = 17.5  # in degrees
    gap_length = 5  # in degrees
    start_angle = 2.5  # in degrees

    # Compute arc points
    def compute_arc_points(r, z, start_phi, arc_length, resolution=100):
        phi_values = np.linspace(np.radians(start_phi),
                                 np.radians(start_phi + arc_length),
                                 resolution)
        x_values = r * np.cos(phi_values)
        y_values = r * np.sin(phi_values)
        return x_values, y_values, np.array([z for _ in phi_values])


    # Loop to plot all the arcs
    for i in range(num_elm_coils):

        phi0 = start_angle + i * (arc_length + gap_length)
        phi1 = phi0 + arc_length

        for r_elm, z_elm in [(r_elm_0, z_elm_0), (r_elm_1, z_elm_1)]:
            for r, z in zip(r_elm, z_elm):
                x, y, z = compute_arc_points(r, z, phi0, arc_length)
                ax.plot(x, y, z, c='tab:red')
            x00 = r_elm[0] * np.cos(np.radians(phi0))
            x01 = r_elm[1] * np.cos(np.radians(phi0))
            y00 = r_elm[0] * np.sin(np.radians(phi0))
            y01 = r_elm[1] * np.sin(np.radians(phi0))
            ax.plot([x00, x01], [y00, y01], z_elm, c='tab:red')
            x10 = r_elm[0] * np.cos(np.radians(phi1))
            x11 = r_elm[1] * np.cos(np.radians(phi1))
            y10 = r_elm[0] * np.sin(np.radians(phi1))
            y11 = r_elm[1] * np.sin(np.radians(phi1))
            ax.plot([x10, x11], [y10, y11], z_elm, c='tab:red')

    # RWM

    num_elm_coils = 16

    elm_coil_r00 = 6.22
    elm_coil_r01 = 5.13
    elm_coil_z00 = 3.5
    elm_coil_z01 = 5.92
    r_elm_0 = [elm_coil_r00, elm_coil_r01]
    z_elm_0 = [elm_coil_z00, elm_coil_z01]

    elm_coil_r10 = 6.22
    elm_coil_r11 = 5.13
    elm_coil_z10 = -3.5
    elm_coil_z11 = -5.92
    r_elm_1 = [elm_coil_r10, elm_coil_r11]
    z_elm_1 = [elm_coil_z10, elm_coil_z11]

    arc_length = 17.5  # in degrees
    gap_length = 5  # in degrees
    start_angle = 2.5  # in degrees

    # Loop to plot all the arcs
    for i in range(num_elm_coils):

        phi0 = start_angle + i * (arc_length + gap_length)
        phi1 = phi0 + arc_length

        for r_elm, z_elm in [(r_elm_0, z_elm_0), (r_elm_1, z_elm_1)]:
            for r, z in zip(r_elm, z_elm):
                x, y, z = compute_arc_points(r, z, phi0, arc_length)
                ax.plot(x, y, z, c='tab:blue')
            x00 = r_elm[0] * np.cos(np.radians(phi0))
            x01 = r_elm[1] * np.cos(np.radians(phi0))
            y00 = r_elm[0] * np.sin(np.radians(phi0))
            y01 = r_elm[1] * np.sin(np.radians(phi0))
            ax.plot([x00, x01], [y00, y01], z_elm, c='tab:blue')
            x10 = r_elm[0] * np.cos(np.radians(phi1))
            x11 = r_elm[1] * np.cos(np.radians(phi1))
            y10 = r_elm[0] * np.sin(np.radians(phi1))
            y11 = r_elm[1] * np.sin(np.radians(phi1))
            ax.plot([x10, x11], [y10, y11], z_elm, c='tab:blue')

    # RWM ACC

    num_elm_coils = 8

    elm_coil_r00 = 7.46
    elm_coil_r01 = 7.46
    elm_coil_z00 = 1.81
    elm_coil_z01 = -1.81
    r_elm_0 = [elm_coil_r00, elm_coil_r01]
    z_elm_0 = [elm_coil_z00, elm_coil_z01]

    elm_coil_r10 = 7.48
    elm_coil_r11 = 7.46
    elm_coil_z10 = 1.81
    elm_coil_z11 = -1.81
    r_elm_1 = [elm_coil_r10, elm_coil_r11]
    z_elm_1 = [elm_coil_z10, elm_coil_z11]

    arc_length = 40  # in degrees
    gap_length = 5  # in degrees
    start_angle = 2.5  # in degrees

    # Loop to plot all the arcs
    for i in range(num_elm_coils):

        phi0 = start_angle + i * (arc_length + gap_length)
        phi1 = phi0 + arc_length

        for r_elm, z_elm in [(r_elm_0, z_elm_0), (r_elm_1, z_elm_1)]:
            for r, z in zip(r_elm, z_elm):
                x, y, z = compute_arc_points(r, z, phi0, arc_length)
                ax.plot(x, y, z, c='tab:green')
            x00 = r_elm[0] * np.cos(np.radians(phi0))
            x01 = r_elm[1] * np.cos(np.radians(phi0))
            y00 = r_elm[0] * np.sin(np.radians(phi0))
            y01 = r_elm[1] * np.sin(np.radians(phi0))
            ax.plot([x00, x01], [y00, y01], z_elm, c='tab:green')
            x10 = r_elm[0] * np.cos(np.radians(phi1))
            x11 = r_elm[1] * np.cos(np.radians(phi1))
            y10 = r_elm[0] * np.sin(np.radians(phi1))
            y11 = r_elm[1] * np.sin(np.radians(phi1))
            ax.plot([x10, x11], [y10, y11], z_elm, c='tab:green')

    # define parameters
    num_tf_coils = 16
    tf_coil_r_min = 0
    tf_coil_r_max = 9
    tf_coil_z_min = -11
    tf_coil_z_max = 11

    def get_coil_vertices(r_min, r_max, z_min, z_max, phi):
        """
        Function to get the X, Y, Z coordinates of the coil for a given phi
        """
        # Points in cylindrical coordinates
        r_points = [r_min, r_max, r_max, r_min, r_min]
        z_points = [z_min, z_min, z_max, z_max, z_min]

        # Convert cylindrical to cartesian coordinates
        x_coords = [r * np.cos(phi) for r in r_points]
        y_coords = [r * np.sin(phi) for r in r_points]
        z_coords = z_points

        return x_coords, y_coords, z_coords

    for i in range(num_tf_coils):
        phi_tf = 2 * np.pi * i / num_tf_coils
        x_tf, y_tf, z_tf = get_coil_vertices(tf_coil_r_min, tf_coil_r_max,
                                            tf_coil_z_min, tf_coil_z_max, phi_tf)
        ax.plot(x_tf, y_tf, z_tf, color='k')

    # equal aspect ratio
    axis_limits = [tf_coil_z_min, tf_coil_z_max]
    ax.set_xlim(axis_limits)
    ax.set_ylim(axis_limits)
    ax.set_zlim(axis_limits)
    ax.grid(False)
    ax.axis('off')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

    orange_patch = mpatches.Patch(color='tab:orange',
                                  label='Last Closed Flux Surface')
    legend_line_blue = Line2D([0], [0], color='tab:blue',
                            label='In-vessel ELM' +
                            ' suppression coils/\nPassive RWM suppression coils')
    legend_line_green = Line2D([0], [0], color='tab:green',
                            label='RWM Active' +
                            ' control coils')
    legend_line_red = Line2D([0], [0], color='tab:red',
                            label='Ex-vessel ELM' +
                            ' suppression coils/\nError field correction coils')
    legend_line_black = Line2D([0], [0], color='k', label='TF coils')
    leg = ax.legend(handles=[orange_patch, legend_line_blue, legend_line_green,
                            legend_line_red, legend_line_black],
                    loc='center left', bbox_to_anchor=(0.75, 0.5))

    ax.view_init(elev=7, azim=0)
    ax.set_box_aspect([1,1,1])

    fig.savefig('plots/coil_plot_3d.png',
                bbox_extra_artists=(leg,),
                bbox_inches=mtransforms.Bbox.from_bounds(4.5, 2.25, 7.5, 5),
                dpi=300)
    fig.savefig('plots/coil_plot_3d.pdf',
                bbox_extra_artists=(leg,),
                bbox_inches=mtransforms.Bbox.from_bounds(4.5, 2.25, 7.5, 5))
    plt.close('all')

if __name__ == '__main__':

    repo_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    gfile_path_name = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    coil_plot_3d(gfile_path=gfile_path_name)
