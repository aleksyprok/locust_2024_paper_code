"""
This script contains routines to reduce the velocity of markers
in a markers file.
"""

import os
import numpy as np

# if __name__ == "__main__":
#     repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
#     partciles_path = os.path.join(repo_dir, "input_data",
#                                   "SPR-068-RV_markers_1000000.dat")
#     data = np.loadtxt(partciles_path,
#                       skiprows=2)
#     data[:, 3:6] *= 0.8
#     np.savetxt(partciles_path, data)
