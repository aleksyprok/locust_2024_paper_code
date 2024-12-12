#!/bin/bash
# The is a shell script to run the locust code in interactive mode.

#device="csd3"
device="leonardo"
tokamak="STEP"

run_name="FEC_2024"

if [[ $device == "aws_v100" ]]; then
    export PATH=$PATH:/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/compilers/bin
    ngpu=1
elif [[ $device == "csd3" ]]; then
    . /etc/profile.d/modules.sh
    module purge
    module load rhel8/default-amp
    module load nvhpc/22.3/gcc-9.4.0-ywtqynx
    module load hdf5/1.10.7/openmpi-4.1.1/nvhpc-22.3-strpuv5
    export HDF5_DIR="/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/nvhpc-22.3/hdf5-1.10.7-strpuv55e7ggr5ilkjrvs2zt3jdztwpv"
    ngpu=1
elif [[ $device == "leonardo" ]]; then
    account="FUAL8_UKAEA_ML"
    partition="boost_fua_prod"
    time="24:00:00"
    ngpu=4
else
    echo "Invalid device."
    exit 1
fi

export OMP_NUM_THREADS=$ngpu
ulimit -s 2000000
export OMP_STACKSIZE=102400
export CUDA_CACHE_DISABLE=1

echo "OMP_NUM_THREADS="$OMP_NUM_THREADS

# Clear CacheFiles
echo $HOSTNAME
rm -vf $HOME"/locust."$tokamak"/CacheFiles/"$HOSTNAME"/"*
$HOME"/locust/locust_"$run_name"_"4
