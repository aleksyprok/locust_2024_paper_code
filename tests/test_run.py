"""
Test module for the run.py module.
"""

import os
import matplotlib.pyplot as plt
from python_scripts import run

def test_run_init():
    """
    Test the Run class initialization.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    assert test_run.dir_path == dir_path
    assert test_run.tag == tag
    assert test_run.log is None
    assert test_run.markers is None
    assert test_run.log_path == dir_path + f'/LOG_{tag}.out'
    assert test_run.fstate_path == dir_path + f'/FINAL_STATE_{tag}.dat'

def test_create_runs_list():
    """
    Test the create_runs_list function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles")
    runs = run.create_runs_list(dir_path)
    assert len(runs) == 7

def test_stopped_scatter():
    """
    Check we can read a run and make a scatter plot of the stopped particles.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()

    fig, ax = plt.subplots()
    ax.scatter(test_run.markers.stopped.r, test_run.markers.stopped.z)
    ax.plot(test_run.wall.r, test_run.wall.z, 'k')
    ax.plot(test_run.gfile.R_bnd, test_run.gfile.Z_bnd, 'r')
    ax.set_aspect('equal')
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    output_dir = os.path.join(repo_path, 'tests', 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'stopped_particles.png')
    fig.savefig(output_path,
                bbox_inches='tight', dpi=300)
    plt.close(fig)
