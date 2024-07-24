# nuclear_fusion_2024_code

The codebase associated with our research, as detailed in the paper available for review and use at the following GitHub repository: <https://github.com/aleksyprok/Prokopyszyn_et_al_Nuclear_Fusion_IAEA_2023>.

### External Data Storage

To keep this repository manageable, we have stored large datasets on Google Drive and OneDrive. You can access the necessary data through the following links:

- **Input Data and Plots**: Available on [Google Drive](https://drive.google.com/drive/folders/19O-Q3VZ2t9bzik9bUeZ-YEsVlpDDwWD7?usp=share_link).
- **Simulation Output Data**: Available on [OneDrive](https://stfc365-my.sharepoint.com/:f:/g/personal/alexander_prokopyszyn_stfc_ac_uk/Ev5nc3ZpORVPqvzkBBQJYCUBEnaOCg2FhIZq7Fq3rOzMSA?email=alex.prokopyszyn%40ukaea.uk&e=ZBfYRV). Note that access to this OneDrive folder requires permission, as it is stored on a work account and exceeds 30 GB.

## Repository Structure

This repository is organized as follows:

`/input_data` - Input datasets required by the scripts are stored here.

  - `BPLASMA/` - Contains the numerically calculated 3D magnetic fields. These were calculated using MARS-F and were provided by Dave Ryan and Guoliang Xia.
  - `LOCUST_SPR-045-14_OutputFiles/` - Contains output data from LOCUST simulations. We use this in `tests/` to run our unit tests for the modules in `python_scripts/`.
  - `base.f90` - This is a copy of `prec_mod.f90` from the original LOCUST code, it contains placeholder values, e.g. `file_eqm = 'eqdsk_file.eqdsk' ! apkp` which we later change using `sed` commands in `shell_scripts/compile.sh`.
  - `ion_info_SPR-045-16.dat` - This is produced by `python_scripts/prepare_profiles.py` and contains the information for the plasma ions.
  - `makefile_template` - The is a copy of `makefile` from the original LOCUST code, it contains placeholder values e.g. `CUDALIB = ccxx,cudaxx.x` which we later change using `sed` commands in `shell_scripts/compile.sh`.
  - `profile_SPR-045-16_ne.dat` - This is produced by `python_scripts/prepare_profiles.py` and contains the profile data and give the number density of the electrons as a function of $\psi_n$.
  - `profile_SPR-045-16_Te.dat` - Similar to `profile_SPR-045-16_ne.dat` but for the electron temperatures.
  - `profile_SPR-045-16_Ti.dat` - Similar to `profile_SPR-045-16_ne.dat` but for the ion temperatures.
  - `profiles_SPR-045-16.CDF` - Outputted data from a JETTO run in CDF format, specifically SPR-045-16 see (https://simdb.step.ukaea.uk/alias/smars/jetto/step/88888/dec1023/seq-1). `python_scripts/prepare_profiles.py` reads this file in to produce the `profile_SPR-045-16_ne.dat` and `profile_SPR-045-16_Te.dat` files.
  - `SPP-001_wall.dat` - Contains a 2D trace of the inner wall of the tokamak in the poloidal plane. We use this to make plots later.
  - `SPP-001-1.cdb.locust` - The full tetrahedral volumetric wall mesh used by LOCUST.
  - `SPR-045-16_markers_1000000.dat` - The initial position and velocity of the alpha particle markers. This file is produced by running `python_scripts/prepare_profiles.py` on the `profiles_SPR-045-16.CDF` file.
  - `SPR-045-16.eqdsk` - The magnetic field which we took from the JETTO run above.

`/output_data/FEC-2024/` - Contains the output data from the LOCUST simulations.

`/plots` - Contains plots produced by running `python_scripts/paper_plots.py`. Some of these plots are also available in the paper.

`/python_scripts` - Core Python scripts for processing and analysis:

  - `__init__.py` - Indicates that this directory is treated as a Python package.
  - `bootstrap.py` - Module for calculating the bootstrap errors, i.e. for estimating the errors due to Monte Carlo noise.
  - `dt_fusion.py` - dt here stand for Deuterium-Tritium. This module contains routines for calculating the fusion power and reaction rate.
  - `flux.py` - This module contains routines for calculating the alpha particle energy flux on the inner wall of the tokamak.
  - `generate_alphas.py` - This module contains routines for calculating the initial position and velocity of the alpha particle markers and is used in `prepare_profiles.py`.
  - `log.py` - This module contains routines for reading in the `LOG` files from the LOCUST simulations and storing useful information as attributes in the `log` object for use in the other modules.
  - `markers.py` - This module reads in the `FINAL_STATE` files from the LOCUST simulations and stores useful information as attributes in the `markers` object for use in the other modules.
  - `my_gfile_reader.py` - Module for reading in eqdsk files (not written by me).
  - `paper_plots_3d.py` - Contain the code for producing the 3D plot which appears in the paper. This module is used in `paper_plots.py`.
  - `paper_plots_extra.py` - Contains code for producing the extra plots which don't appear in the paper but are still useful for checking the results. This module is used in `paper_plots.py`.
  - `paper_plots.py` - This module contains the code for producing the plots that appear in the paper.
  - `pkde.py` - Here PKDE stands for Periodic Kernel Density Estimation. This routine is a self-contained designed to calculate periodic Kernel Density Estimates in 1D and 2D for a set of points with weights. Note that `flux.py` uses this to calculate the alpha particle energy flux on the inner wall of the tokamak.
  - `prepare_profiles.py` - The main script that processes profile data for simulations. It produces data that LOCUST needs to run.
  - `run.py` - This module contains code for the `run` object which is a class which stores information about an individual LOCUST simulation. It has attributes called `log`, `wall`, `markers`, `gfile`, `flux` which store information about the `LOG` files, wall, `FINAL_STATE files`, `.eqdsk` files, and alpha particle energy flux on the wall. The `run.py` module calls on the `flux`, `log`, `wall`, `markers` and  `my_gfile_reader` modules.
  - `wall.py` - This module contains routines for reading in the `SPP-001_wall.dat` file and storing useful information as attributes in the `wall` object for use in the other modules.

`shell_scripts/` Contains shell scripts intended to be run on a GPU cluster to compile, run and post process the LOCUST simulations.

  - `compile.sh` - This complies LOCUST and produces a binary for each set of parameters we need to make the plots in the paper.
  - `extract_log_and_fstate.sh` - This extracts the data we need from the `LOG` files and `FINAL_STATE` files.
  - `run_interactive.sh` - This is only needed if you intend to run a LOCUST simulation in interactive mode.
  - `run.sh` - This is intended to be run after compiling and executes the binaries produced by `compile.sh`.


`/tests` - Contains unit tests for the modules in `python_scripts/`.

`/.gitattributes` - Specifies attributes for path names and handling of line endings across different environments. This is also where we specify the large files for the Git Large File Storage (LFS).

`/.gitignore` - Lists files and directories that are to be ignored by version control.

`LICENSE` - Contains the LICENSE file. It says that I am happy for anyone to use this code as they wish, however, I'm not sure what UKAEA's policy is!

`/README.md` - Provides an overview and documentation for the project.

 `/setup.py` - Used for packaging and distributing the project, defining the package name, version, and included packages.

 ## Installation instructions:

  - To run the python scripts you need to have python3.9+ installed.
  - Setup a virtual environment with e.g. "python3 -m venv venv".
  - Activate the environment with e.g. "source venv/bin/activate"
  - Install python packages with "pip install -r requirements.txt"
  - Install the repo itself with "pip install -e .".

## How to produce data for the paper

### Step 1: Produce density/temperature profiles as well as alpha particle initial position and velocity.

The first thing we did was produce the profile data. This was done using the `prepare_profiles.py` script. This script takes the raw data from the input directory and processes it into a format that can be used for simulations. Note that I have already done this to produce:
  - `input_data/ion_info_SPR-045-16.dat`
  - `input_data/profile_SPR-045-16_ne.dat`
  - `input_data/profile_SPR-045-16_Te.dat`
  - `input_data/profile_SPR-045-16_Ti.dat`
  - `input_data/SPR-045-16_markers_1000000.dat`

### Step 2: ssh into GPU cluster

For this work I used the CSD3 cluster:
https://docs.hpc.cam.ac.uk/hpc/index.html

### Step 3: Clone locust code

For the shell scripts to work the `locust` code must be cloned into the home space and you need to run
``git checkout prok``
to get the same branch of the code. The LOCUST repo is stored on GitLab here:
[https://git.ccfe.ac.uk/fast-particles/locust](https://git.ccfe.ac.uk/fast-particles/locust)

### Step 4: Clone this repository

Clone this repo to a directory of your choosing on CSD3. Make sure you have `git-lfs`, I load these modules in my `.bashrc` to ensure that I have `git-lfs`. Note that this also sets up `slurm` which you will need to run jobs:
```
module purge
module load gcc/5
module load rhel8/default-amp
module load slurm-current-gcc-5.4.0-6idu76o
module load git-lfs-2.3.0-gcc-5.4.0-oktvmkw
```

### Step 5: Create output directories

Create the file structure for LOCUST to output its files to. Your file structure needs to be:
- `~/locust.STEP/OutputFiles`
- `~/locust.STEP/InputFiles`
- `~/locust.STEP/CacheFiles`

I recommend using symbolic links to ensure you files are saved in e.g. `/rds/project/iris_vol2/rds-ukaea-ap001/` as you only have a limited amount of space in your home drive.

### Step 6: Compile the code.

Navigate to `nuclear_fusion_2024_code/shell_scripts/compile.sh`. Edit the the file and ensure `device="csd3"` (assuming you are running the code on CSD3). Then create a session on a GPU node with e.g.:
```
sintr -t 36:0:0 -N1 -A UKAEA-AP001-GPU -p ampere --qos=INTR --gres=gpu:1
```
Then run `./compile.sh` to compile the code.

### Step 7: Run the code.

After the code has finished compiling. Edit `nuclear_fusion_2024_code/shell_scripts/run.sh` and ensure the settings are ok. You may need to edit the lines that read
```
#SBATCH --ntasks=248
#SBATCH --array=0-247
```
to 
```
#SBATCH --ntasks=124
#SBATCH --array=0-123
```
as sometime CSD3 won't let you execute so many jobs as a single array job. So you need to submit two separate array jobs.

### Step 8: Filter out essential data.

After the simulations have finished, you should be able to see the outputted data in `~/locust.STEP/OutputFiles`.
We now need to run `shell_scripts/extract_log_and_fstate.sh` to get the data we need. I recommend putting all the output data in a directory e.g.
`~/locust.STEP/OutputFiles/output_data_raw`, then setting
```
SRC_DIR=$HOME"/locust.STEP/OutputFiles/output_data_raw"
DEST_DIR=$HOME"/locust.STEP/OutputFiles/output_data_processed"
```
in the shell script before running it. This script extracts only the data we need, by copying the `FINAL_STATE` and `LOG` files as well as removing redundant lines
of text from the `LOG` file. This saves storage space and significantly reduces the time it take to the read the `LOG` files later.

### Step 9: Download data to your local machine.

Now download the data from CSD3 to your local machine. I put the data in `output_data/FEC_2024`.

### Step 10: Produce plots.

To produce plots simply execute `python_scripts/paper_plots.py`. Note that you need to ensure that you have followed the installation instructions above first. The code should take about one minute to run with the `make_csv` and `save_axisymmetric` values set to `False`. With these to `True` the code will take about 6 hours to run. With `make_csv` and `save_axisymmetric` set to `True` the code will perform a full run of the code and calculate the optimum bandwidth to use and calculate the bootstrap errors. However, I have already calculated these values and saved them in `csv` files in:
  - `plots/ripple_runs.csv`.
  - `plots/rmp_runs.csv`.
  - `plots/rwm_runs.csv`.

