"""
Test the markers module.
"""
import os
import numpy as np
from python_scripts import markers, run

def test_particle_group_init():
    """
    Test the ParticleGroup class initialization.
    """
    test_particle_group = markers.ParticleGroup()
    assert test_particle_group.r.size < 1
    assert test_particle_group.phi.size < 1
    assert test_particle_group.z.size < 1
    assert test_particle_group.vr.size < 1
    assert test_particle_group.vphi.size < 1
    assert test_particle_group.vz.size < 1
    assert test_particle_group.t.size < 1
    assert test_particle_group.s.size < 1
    assert test_particle_group.particle_id.size < 1
    assert test_particle_group.r0.size < 1
    assert test_particle_group.phi0.size < 1
    assert test_particle_group.z0.size < 1
    assert test_particle_group.vr0.size < 1
    assert test_particle_group.vphi0.size < 1
    assert test_particle_group.vz0.size < 1
    assert test_particle_group.weight.size < 1

def test_particle_group_add_particles():
    """
    Test the ParticleGroup class add_particles method.
    """
    test_particle_group = markers.ParticleGroup()
    data = np.array([
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0],
        [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
    ])
    test_particle_group.add_particles(data)
    assert np.all(test_particle_group.r == np.array([1.0, 2.0]))
    assert np.all(test_particle_group.phi == np.array([2.0, 3.0]))
    assert np.all(test_particle_group.z == np.array([3.0, 4.0]))
    assert np.all(test_particle_group.vr == np.array([4.0, 5.0]))
    assert np.all(test_particle_group.vphi == np.array([5.0, 6.0]))
    assert np.all(test_particle_group.vz == np.array([6.0, 7.0]))
    assert np.all(test_particle_group.t == np.array([7.0, 8.0]))
    assert np.all(test_particle_group.s == np.array([8.0, 9.0]))
    assert np.all(test_particle_group.particle_id == np.array([9.0, 10.0]))
    assert np.all(test_particle_group.r0 == np.array([10.0, 11.0]))
    assert np.all(test_particle_group.phi0 == np.array([11.0, 12.0]))
    assert np.all(test_particle_group.z0 == np.array([12.0, 13.0]))
    assert np.all(test_particle_group.vr0 == np.array([13.0, 14.0]))
    assert np.all(test_particle_group.vphi0 == np.array([14.0, 15.0]))
    assert np.all(test_particle_group.vz0 == np.array([15.0, 16.0]))
    assert np.all(test_particle_group.weight == np.array([16.0, 17.0]))

def test_markers_init():
    """
    Test the Markers class initialization.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    fstate_path = os.path.join(input_dir, "FINAL_STATE_13-12-2023_16-51-52.811.dat")
    axisymmetric_markers = markers.Markers(fstate_path)
    assert axisymmetric_markers.fstate_path == fstate_path
    assert len(axisymmetric_markers.moving.r) == 18
    assert len(axisymmetric_markers.thermal.r) == 482343
    assert len(axisymmetric_markers.stopped.r) == 41927
    assert len(axisymmetric_markers.all.r) == 524288
    assert len(axisymmetric_markers.unresolved.r) == 0
    assert abs(axisymmetric_markers.stopped.r[0] - 5.77736) < 1e-5
    assert abs(axisymmetric_markers.stopped.phi[0] + 1.56278) < 1e-5
    assert abs(axisymmetric_markers.stopped.z[0] - 8.53939) < 1e-5
    assert abs(axisymmetric_markers.stopped.vr[0] + 114536.0) < 1e-5
    assert abs(axisymmetric_markers.stopped.vphi[0] + 577659.0) < 1e-5
    assert abs(axisymmetric_markers.stopped.vz[0] - 347993.0) < 1e-5
    assert abs(axisymmetric_markers.stopped.t[0] - 0.833556) < 1e-5
    assert abs(axisymmetric_markers.stopped.s[0] + 5.0) < 1e-5
    assert int(axisymmetric_markers.stopped.particle_id[0]) == 64954
    assert abs(axisymmetric_markers.stopped.r0[0] - 5.12387) < 1e-5
    assert abs(axisymmetric_markers.stopped.phi0[0] - 2.37106) < 1e-5
    assert abs(axisymmetric_markers.stopped.z0[0] - 2.5327) < 1e-5
    assert abs(axisymmetric_markers.stopped.vr0[0] - 8201490.0) < 1e-5
    assert abs(axisymmetric_markers.stopped.vphi0[0] + 8432840.0) < 1e-5
    assert abs(axisymmetric_markers.stopped.vz0[0] + 5514400.0) < 1e-5
    assert abs(axisymmetric_markers.stopped.weight[0] - 23987000000.0) < 1e-5

def test_get_s_phi_s_theta_from_r_z_phi():
    """
    Test the get_s_phi_s_theta_from_r_z_phi function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    test_run.init_markers()
    markers.get_s_phi_s_theta_from_r_z_phi(test_run)
    assert len(test_run.markers.stopped.s_phi) == 41927
    assert len(test_run.markers.stopped.s_theta) == 41927
