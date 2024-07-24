#!/bin/bash

device="csd3"
# device="leonardo"
# device="sdcc"
tokamak="STEP"

run_name="FEC_2024"

if [[ $device == "csd3" ]]; then
    account="ukaea-ap001-GPU"
    partition="ampere"
    time="36:00:00"
    ngpu=4
elif [[ $device == "leonardo" ]]; then
    account="FUAL7_UKAEA_ML"
    partition="boost_fua_prod"
    time="24:00:00"
    ngpu=4
elif [[ $device == "sdcc" ]]; then
    account="default"
    partition="gpu_p100_titan"
    time="99-00:00:00"
    ngpu=1
else
    echo "Invalid device."
    exit 1
fi

# Wrapper for sbatch script
cat > job.sbatch << EOF
#!/bin/bash
#SBATCH -J gpujob
#SBATCH -A $account
#SBATCH --nodes=1
#SBATCH --gres=gpu:$ngpu
#SBATCH --time=$time
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alex.prokopyszyn@ukaea.uk
#SBATCH --output=test_dmtcp_%A_%a.out
#SBATCH --error=test_dmtcp_%A_%a.err
#SBATCH --ntasks=248
#SBATCH --array=0-247
#SBATCH --exclusive
#SBATCH -p $partition

device=$device
if [[ \$device == "csd3" ]]; then
    . /etc/profile.d/modules.sh
    module purge
    module load rhel8/default-amp
    module load nvhpc/22.3/gcc-9.4.0-ywtqynx
    module load hdf5/1.10.7/openmpi-4.1.1/nvhpc-22.3-strpuv5
    export HDF5_DIR="/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/nvhpc-22.3/hdf5-1.10.7-strpuv55e7ggr5ilkjrvs2zt3jdztwpv"
elif [[ \$device == "leonardo" ]]; then
    module purge
    module load nvhpc/23.1
elif [[ \$device == "sdcc" ]]; then
    module purge
    module load IMAS
fi

export OMP_NUM_THREADS=$ngpu

workdir="\$SLURM_SUBMIT_DIR"
cd \$workdir

export HDF5_DIR="/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/nvhpc-22.3/hdf5-1.10.7-strpuv55e7ggr5ilkjrvs2zt3jdztwpv"
export LIBRARY_PATH=$LIBRARY_PATH:"$HDF5_DIR/lib"
export CFLAGS="-I$HDF5_DIR/include"
export FFLAGS="-I$HDF5_DIR/include"
export LDFLAGS="-L$HDF5_DIR/lib"

ulimit -s 2000000
export OMP_STACKSIZE=102400
export CUDA_CACHE_DISABLE=1

echo "OMP_NUM_THREADS="\$OMP_NUM_THREADS

# Clear CacheFiles
echo $HOSTNAME
rm -vf $HOME"/locust."$tokamak"/CacheFiles/"\$HOSTNAME"/"*
# Execute Locust
$HOME"/locust/locust_"$run_name"_"\$SLURM_ARRAY_TASK_ID
EOF

sbatch job.sbatch