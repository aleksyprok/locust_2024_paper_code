#!/bin/bash

# This shell script will compile and produce a range of binary files for LOCUST
# to execute. We compile every simulation used for the IAEA FEC 2024 paper.
# This includes simulations where we have an axisymmetric simulation,
# where we vary the ripple field, where we vary the RMP field and where
# we vary the amplitude of a RWM.

absolute_value() {
    if [ $1 -lt 0 ]; then
      echo $((-$1))
    else
      echo $1
    fi
}

device="leonardo"
tokamak="STEP"
input_dir="$(pwd)/../input_data"
home_dir=$HOME
prec_file=$home_dir"/locust/prec_mod.f90"
mesh_file="SPP-001-1.cdb.locust"
run_name="FEC_2024"
spr_strings_unique=("SPR-045-14" "SPR-045-16" "SPR-068-7")

niter=1
threadsPerBlock=256
blocksPerGrid=256

# Axisymmetric simulations
rcoil=0
rcoil_inner=0
ncoil=0
bripple=0
toroidal_mode=0
coil_set="none"
current=0
response=0
phase=0
gain_value=0
bscale=0
bplasma=0
rwm_control=0
rmp=0
spr_string="none"
dplot=0
timax="10.0_gpu"
unbor=1000
dt0="1.0e-06_gpu"
run_category=1
for spr_string in "${spr_strings_unique[@]}"; do
    rcoils+=("$rcoil")
    rcoils_inner+=("$rcoil_inner")
    ncoils+=("$ncoil")
    bripples+=("$bripple")
    toroidal_modes+=("$toroidal_mode")
    coil_sets+=("$coil_set")
    currents+=("$current")
    responses+=("$response")
    phases+=("$phase")
    gain_values+=("$gain_value")
    bscales+=("$bscale")
    bplasmas+=("$bplasma")
    rwm_controls+=("$rwm_control")
    rmps+=("$rmp")
    spr_strings+=("$spr_string")
    dplots+=("$dplot")
    timaxs+=("$timax")
    unbors+=("$unbor")
    dt0s+=("$dt0")
    run_categories+=("$run_category")
done

# Ripple simulations
rcoil=0
rcoil_inner=0.75
ncoil=0
bripple=1
toroidal_mode=0
coil_set="none"
current=0
response=0
phase=0
gain_value=0
bscale=0
bplasma=0
rwm_control=0
rmp=0
spr_string="SPR-045-16"
dplot=0
timax="10.0_gpu"
unbor=1000
dt0="1.0e-06_gpu"
run_category=2
# ncoils_unique=(12 16 18)
ncoils_unique=(12)
# rcoils_unique=("7.0" "7.25" "7.5" "7.75" "8.0" "8.25" "8.5" "8.75" "9.0")
rcoils_unique=("7.5" "7.75" "8.0" "8.25" "8.5" "8.75" "9.0" "9.25" "9.50")
for ncoil in "${ncoils_unique[@]}"; do
    for rcoil in "${rcoils_unique[@]}"; do
        for spr_string in "${spr_strings_unique[@]}"; do
            ncoils+=("$ncoil")
            rcoils+=("$rcoil")
            rcoils_inner+=("$rcoil_inner")
            bripples+=("$bripple")
            toroidal_modes+=("$toroidal_mode")
            coil_sets+=("$coil_set")
            currents+=("$current")
            responses+=("$response")
            phases+=("$phase")
            gain_values+=("$gain_value")
            bscales+=("$bscale")
            bplasmas+=("$bplasma")
            rwm_controls+=("$rwm_control")
            rmps+=("$rmp")
            spr_strings+=("$spr_string")
            dplots+=("$dplot")
            timaxs+=("$timax")
            unbors+=("$unbor")
            dt0s+=("$dt0")
            run_categories+=("$run_category")
        done
    done
done

# RMP simulations
rcoil=0
rcoil_inner=0
ncoil=0
bripple=0
toroidal_mode=0
coil_set="none"
current=0
response=0
phase=0
gain_value=0
bscale=0
bplasma=1
rwm_control=0
rmp=1
spr_string="SPR-045-16"
dplot=0
timax="10.0_gpu"
unbor=1000
dt0="1.0e-06_gpu"
run_category=3
phases_unique=("000.00" "045.00" "090.00" "135.00"
               "180.00" "225.00" "270.00" "315.00")
toroidal_modes_unique=(-2 -3 -4)
coil_sets_unique=("efcc" "rwm")
responses_unique=(0 1)
current_multipliers_unique=(1 2)
for coil_set in "${coil_sets_unique[@]}"; do
    for toroidal_mode in "${toroidal_modes_unique[@]}"; do
        for current_multiplier in "${current_multipliers_unique[@]}"; do
            if [[ $coil_set == "efcc" ]]; then
                if [[ $(($toroidal_mode * $toroidal_mode)) -eq 4 ]]; then
                    current=$(($current_multiplier * 50))
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 9 ]]; then
                    current=$(($current_multiplier * 90))
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 16 ]]; then
                    current=$(($current_multiplier * 150))
                else
                    echo "Invalid toroidal number"
                    exit 1
                fi
            elif [[ $coil_set == "rwm" ]]; then
                if [[ $(($toroidal_mode * $toroidal_mode)) -eq 4 ]]; then
                    current=$(($current_multiplier * 30))
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 9 ]]; then
                    current=$(($current_multiplier * 50))
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 16 ]]; then
                    current=$(($current_multiplier * 80))
                else
                    echo "Invalid toroidal number"
                    exit 1
                fi
            else
                echo "Invalid coil set"
                exit 1
            fi
            for response in "${responses_unique[@]}"; do
                for phase in "${phases_unique[@]}"; do
                    ncoils+=("$ncoil")
                    rcoils+=("$rcoil")
                    rcoils_inner+=("$rcoil_inner")
                    bripples+=("$bripple")
                    toroidal_modes+=("$toroidal_mode")
                    coil_sets+=("$coil_set")
                    currents+=("$current")
                    responses+=("$response")
                    phases+=("$phase")
                    gain_values+=("$gain_value")
                    bscales+=("$bscale")
                    bplasmas+=("$bplasma")
                    rwm_controls+=("$rwm_control")
                    rmps+=("$rmp")
                    spr_strings+=("$spr_string")
                    dplots+=("$dplot")
                    timaxs+=("$timax")
                    unbors+=("$unbor")
                    dt0s+=("$dt0")
                    run_categories+=("$run_category")
                done
            done
        done
    done
done

# RMP special phases simulations
rcoil=0
rcoil_inner=0
ncoil=0
bripple=0
toroidal_mode=0
coil_set="none"
current=0
response=0
phase=0
gain_value=0
bscale=0
bplasma=1
rwm_control=0
rmp=1
spr_string="SPR-045-16"
dplot=0
timax="10.0_gpu"
unbor=1000
dt0="1.0e-06_gpu"
run_category=4
toroidal_modes_unique=(-2 -3 -4)
coil_sets_unique=("efcc" "rwm")
responses_unique=(0 1)
current_multipliers_unique=(1 2)
for coil_set in "${coil_sets_unique[@]}"; do
    for toroidal_mode in "${toroidal_modes_unique[@]}"; do
        for current_multiplier in "${current_multipliers_unique[@]}"; do
            if [[ $coil_set == "efcc" ]]; then
                if [[ $(($toroidal_mode * $toroidal_mode)) -eq 4 ]]; then
                    current=$(($current_multiplier * 50))
                    phase="061.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 9 ]]; then
                    current=$(($current_multiplier * 90))
                    phase="020.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 16 ]]; then
                    current=$(($current_multiplier * 150))
                    phase="321.0"
                else
                    echo "Invalid toroidal number"
                    exit 1
                fi
            elif [[ $coil_set == "rwm" ]]; then
                if [[ $(($toroidal_mode * $toroidal_mode)) -eq 4 ]]; then
                    current=$(($current_multiplier * 30))
                    phase="265.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 9 ]]; then
                    current=$(($current_multiplier * 50))
                    phase="173.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 16 ]]; then
                    current=$(($current_multiplier * 80))
                    phase="067.0"
                else
                    echo "Invalid toroidal number"
                    exit 1
                fi
            else
                echo "Invalid coil set"
                exit 1
            fi
            for response in "${responses_unique[@]}"; do
                ncoils+=("$ncoil")
                rcoils+=("$rcoil")
                rcoils_inner+=("$rcoil_inner")
                bripples+=("$bripple")
                toroidal_modes+=("$toroidal_mode")
                coil_sets+=("$coil_set")
                currents+=("$current")
                responses+=("$response")
                phases+=("$phase")
                gain_values+=("$gain_value")
                bscales+=("$bscale")
                bplasmas+=("$bplasma")
                rwm_controls+=("$rwm_control")
                rmps+=("$rmp")
                spr_strings+=("$spr_string")
                dplots+=("$dplot")
                timaxs+=("$timax")
                unbors+=("$unbor")
                dt0s+=("$dt0")
                run_categories+=("$run_category")
            done
        done
    done
done

# RWM control simulations
rcoil=0
rcoil_inner=0
ncoil=0
bripple=0
toroidal_mode=-1
coil_set="none"
current=90
response=1
phase="020.0"
gain_value=0
bscale=0
bplasma=1
rwm_control=0
rmp=1
spr_string="SPR-045-16"
dplot=0
timax="10.0_gpu"
unbor=1000
dt0="1.0e-06_gpu"
run_category=5
bscales_unique=(1 10 100 1000)
for bscale in "${bscales_unique[@]}"; do
    ncoils+=("$ncoil")
    rcoils+=("$rcoil")
    rcoils_inner+=("$rcoil_inner")
    bripples+=("$bripple")
    toroidal_modes+=("$toroidal_mode")
    coil_sets+=("$coil_set")
    currents+=("$current")
    responses+=("$response")
    phases+=("$phase")
    gain_values+=("$gain_value")
    bscales+=("$bscale")
    bplasmas+=("$bplasma")
    rwm_controls+=("$rwm_control")
    rmps+=("$rmp")
    spr_strings+=("$spr_string")
    dplots+=("$dplot")
    timaxs+=("$timax")
    unbors+=("$unbor")
    dt0s+=("$dt0")
    run_categories+=("$run_category")
done

# Ripple simulations - scan in TF_Rcoil_inner
rcoil=0
rcoil_inner=0
ncoil=16
bripple=1
toroidal_mode=0
coil_set="none"
current=0
response=0
phase=0
gain_value=0
bscale=0
bplasma=0
rwm_control=0
rmp=0
spr_string="SPR-045-16"
dplot=0
timax="10.0_gpu"
unbor=1000
dt0="1.0e-06_gpu"
run_category=6
rcoils_unique=("7.5" "100.0")
rcoils_inner_unique=("0.0" "0.25" "0.5" "0.75" "1.0" "1.25" "1.50")
for rcoil in "${rcoils_unique[@]}"; do
    for rcoil_inner in "${rcoils_inner_unique[@]}"; do
        ncoils+=("$ncoil")
        rcoils+=("$rcoil")
        rcoils_inner+=("$rcoil_inner")
        bripples+=("$bripple")
        toroidal_modes+=("$toroidal_mode")
        coil_sets+=("$coil_set")
        currents+=("$current")
        responses+=("$response")
        phases+=("$phase")
        gain_values+=("$gain_value")
        bscales+=("$bscale")
        bplasmas+=("$bplasma")
        rwm_controls+=("$rwm_control")
        rmps+=("$rmp")
        spr_strings+=("$spr_string")
        dplots+=("$dplot")
        timaxs+=("$timax")
        unbors+=("$unbor")
        dt0s+=("$dt0")
        run_categories+=("$run_category")
    done
done

# TF Ripple + RMP + RWM
rcoil=9.0
rcoil_inner=0.75
ncoil=16
bripple=1
toroidal_mode=-3
coil_set="efcc"
current=90
response=0
phase="020.0"
gain_value=0
bscale=0
bplasma=1
rwm_control=0
rmp=1
spr_string="none"
dplot=0
timax="10.0_gpu"
unbor=1000
dt0="1.0e-06_gpu"
run_category=7
for spr_string in "${spr_strings_unique[@]}"; do
    ncoils+=("$ncoil")
    rcoils+=("$rcoil")
    rcoils_inner+=("$rcoil_inner")
    bripples+=("$bripple")
    toroidal_modes+=("$toroidal_mode")
    coil_sets+=("$coil_set")
    currents+=("$current")
    responses+=("$response")
    phases+=("$phase")
    gain_values+=("$gain_value")
    bscales+=("$bscale")
    bplasmas+=("$bplasma")
    rwm_controls+=("$rwm_control")
    rmps+=("$rmp")
    spr_strings+=("$spr_string")
    dplots+=("$dplot")
    timaxs+=("$timax")
    unbors+=("$unbor")
    dt0s+=("$dt0")
    run_categories+=("$run_category")
done

# Make DPLOTS of the full 3D simulations
rcoil=9.0
rcoil_inner=0.75
ncoil=16
bripple=1
toroidal_mode=-3
coil_set="efcc"
current=90
response=0
phase="020.0"
gain_value=0
bscale=0
bplasma=1
rwm_control=0
rmp=1
spr_string="SPR-045-16"
dplot=1
dt0="none"
unbor=0
timax="none"
run_category=8
dt0s_unique=("1.0e-09_gpu" "1.0e-08_gpu" "1.0e-07_gpu" "1.0e-06_gpu" "1.0e-05_gpu")
unbors_unique=(1 10 100 1000 10000)
timaxs_unique=("1.0e-04_gpu" "1.0e-03_gpu" "1.0e-02_gpu" "1.0e-01_gpu" "1.0e+00_gpu")
num_dt0s=${#dt0s_unique[@]}
for ((i = 0; i < num_dt0s; i++)); do
    dt0="${dt0s_unique[$i]}"
    unbor="${unbors_unique[$i]}"
    timax="${timaxs_unique[$i]}"
    ncoils+=("$ncoil")
    rcoils+=("$rcoil")
    rcoils_inner+=("$rcoil_inner")
    bripples+=("$bripple")
    toroidal_modes+=("$toroidal_mode")
    coil_sets+=("$coil_set")
    currents+=("$current")
    responses+=("$response")
    phases+=("$phase")
    gain_values+=("$gain_value")
    bscales+=("$bscale")
    bplasmas+=("$bplasma")
    rwm_controls+=("$rwm_control")
    rmps+=("$rmp")
    spr_strings+=("$spr_string")
    dplots+=("$dplot")
    timaxs+=("$timax")
    unbors+=("$unbor")
    dt0s+=("$dt0")
    run_categories+=("$run_category")
done

if [[ $device == "csd3" ]]; then
    module purge
    module load rhel8/default-amp
    module load nvhpc/22.3/gcc-9.4.0-ywtqynx
    module load hdf5/1.10.7/openmpi-4.1.1/nvhpc-22.3-strpuv5
    export HDF5_DIR="/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/nvhpc-22.3/hdf5-1.10.7-strpuv55e7ggr5ilkjrvs2zt3jdztwpv"
    export LIBRARY_PATH=$LIBRARY_PATH:"$HDF5_DIR/lib"
    export CFLAGS="-I$HDF5_DIR/include"
    export FFLAGS="-I$HDF5_DIR/include"
    export LDFLAGS="-L$HDF5_DIR/lib"
    root_dir="/home"
    cc="80"
    cuda="11.6"
    nohdf5=0
elif [[ $device == "leonardo" ]]; then
    module purge
    module load nvhpc/23.1
    root_dir="/leonardo/home/userexternal"
    cc="80"
    cuda="11.8"
    nohdf5=1
elif [[ $device == "aws_v100" ]]; then
    export PATH=$PATH:/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/compilers/bin
    root_dir="/home"
    cc="70"
    cuda="11.0"
    nohdf5=1
else
    echo "Invalid device."
    exit 1
fi

cp -vf \
$input_dir/"makefile_template" \
$home_dir"/locust/makefile"
SRC="ccxx,cudaxx\.x"
DST="cc"$cc",cuda"$cuda
echo "nohdf5 = "$nohdf5
sed -i "s/$SRC/$DST/g" "$home_dir/locust/makefile"
if [[ $nohdf5 == 1 ]]; then
    SRC="-lstdc++ -lhdf5_fortran -lhdf5"
    DST="-lstdc++"
    sed -i "s/$SRC/$DST/g" "$home_dir/locust/makefile"
fi
diff $input_dir/"makefile_template" "$home_dir/locust/makefile"

rsync -avh \
$input_dir"/"$mesh_file \
$home_dir"/locust."$tokamak"/InputFiles/."

cd $home_dir"/locust"

num_runs=${#rcoils[@]}
for ((n=0; n<num_runs; n++)); do

    echo "n="$n
    # If run category is not x skip
    if [[ ${run_categories[$n]} -ge 3 ]]; then
        continue
    fi

    rcoil=${rcoils[$n]}
    rcoil_inner=${rcoils_inner[$n]}
    ncoil=${ncoils[$n]}
    bripple=${bripples[$n]}
    toroidal_mode=${toroidal_modes[$n]}
    coil_set=${coil_sets[$n]}
    current=${currents[$n]}
    response=${responses[$n]}
    phase=${phases[$n]}
    gain_value=${gain_values[$n]}
    bscale=${bscales[$n]}
    bplasma=${bplasmas[$n]}
    rwm_control=${rwm_controls[$n]}
    rmp=${rmps[$n]}
    spr_string=${spr_strings[$n]}
    unbor=${unbors[$n]}
    dt0=${dt0s[$n]}
    timax=${timaxs[$n]}
    dplot=${dplots[$n]}
    run_category=${run_categories[$n]}
    echo "ncoil="$ncoil
    echo "rcoil="$rcoil
    echo "rcoil_inner"=$rcoil_inner
    echo "bripple="$bripple
    echo "toroidal_mode="$toroidal_mode
    echo "coil_set="$coil_set
    echo "current="$current
    echo "response="$response
    echo "phase="$phase
    echo "gain_value="$gain_value
    echo "bscale="$bscale
    echo "bplasma="$bplasma
    echo "rwm_control="$rwm_control
    echo "rmp="$rmp
    echo "spr_string="$spr_string
    echo "unbor="$unbor
    echo "dt0="$dt0
    echo "timax="$timax
    echo "dplot="$dplot
    echo "run_category="$run_category

    toroidal_mode_abs=$(absolute_value $toroidal_mode)
    echo "toroidal_mode_abs="$toroidal_mode_abs


    FLAGS_BASE="-DCONLY -DPFCMOD -DTOKHEAD -DFSTATE -DLEIID=6 -DSTDOUT \
                -DSMALLEQ -DOPENTRACK -DOPENTERM -DPSIT=0.7 \
                -DNOTUNE -DUNBOR="$unbor" \
                -DTETALL -DSOLCOL \
                -DRFORCE -DBP -DTIMAX="$timax" -DWREAL -DWLIST"
    if [[ $bplasma == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DB3D -DB3D_EX"
    fi
    if [[ $nohdf5 == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DNOHDF5"
    fi
    if [[ $tokamak == "ITER" ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DTOKAMAK=1"
    elif [[ $tokamak == "STEP" ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DTOKAMAK=10"
    else
        echo "Invalid tokamak."
        exit 1
    fi
    if [[ $bripple == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DBRIPPLE"
    fi
    if [[ $dplot == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DPLOT=100"
    fi
    echo $FLAGS_BASE

    cp -vf \
    $input_dir"/base.f90" \
    $prec_file 
    cp -f \
    $input_dir"/profile_"$spr_string"_ne.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$spr_string"_ne.dat"
    cp -f \
    $input_dir"/profile_"$spr_string"_Ti.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$spr_string"_Ti.dat"
    cp -f \
    $input_dir"/profile_"$spr_string"_Te.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$spr_string"_Te.dat"
    eqdsk_file=$spr_string".eqdsk"
    cp -f \
    $input_dir"/"$eqdsk_file \
    $home_dir"/locust."$tokamak"/InputFiles/."
    ptcle_file=$spr_string"_markers_1000000.dat"
    cp -f \
    $input_dir"/"$ptcle_file \
    $home_dir"/locust."$tokamak"/InputFiles/."
    if [[ $dplot == 1 ]]; then
        ptcle_file_dplot=$spr_string"_markers_dplot_100.dat"
        cp -f \
        $input_dir"/"$ptcle_file \
        $home_dir"/locust."$tokamak"/InputFiles/."
    fi

    SRC="TF_Ncoil = 00"
    DST="TF_Ncoil = "$ncoil
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="TF_Rcoil = 0.00_gpu"
    DST="TF_Rcoil = "$rcoil"_gpu"
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="TF_Rcoil_inner = 0.00_gpu"
    DST="TF_Rcoil_inner = "$rcoil_inner"_gpu"
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="'eqdsk_file.eqdsk' ! apkp"
    DST="'"$eqdsk_file"'"
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="'mesh_file.cdb.locust' ! apkp"
    DST="'"$mesh_file"'"
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="niter  = 0 ! apkp - Needs changing"
    DST="niter  = "$niter
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="'ptcle_file.dat' ! apkp - Needs changing"
    DST="'"$ptcle_file"'"
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="threadsPerBlock = 256"
    DST="threadsPerBlock = "$threadsPerBlock
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="blocksPerGrid   = 512"
    DST="blocksPerGrid   = "$blocksPerGrid
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="dt0    = 0.0e-00_gpu ! apkp - Needs changing"
    DST="dt0    = "$dt0
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="root = '\/home' ! apkp - Needs changing on Marconi"
    DST="root = '"$root_dir"'"
    sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_ne = 'profile_ne_file' ! apkp"
    DST="file_profile_ne = 'profile_"$spr_string"_ne.dat'"
    sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_Te = 'profile_Te_file' ! apkp"
    DST="file_profile_Te = 'profile_"$spr_string"_Te.dat'"
    sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_Ti = 'profile_Ti_file' ! apkp"
    DST="file_profile_Ti = 'profile_"$spr_string"_Ti.dat'"
    sed -i "s|$SRC|$DST|g" $prec_file
    if [[ $spr_string == "SPR-068-7" ]] then
        SRC="R0F    = 1.2_gpu !0.85_gpu"
        DST="R0F    = 1.8_gpu"
        sed -i "s/$SRC/$DST/g" $prec_file
        SRC="R1F    = 6.0_gpu !4.15_gpu"
        DST="R1F    = 7.0_gpu"
        sed -i "s/$SRC/$DST/g" $prec_file
        SRC="Z0F    = -7.0_gpu !-4.30_gpu"
        DST="Z0F    = -7.25_gpu"
        sed -i "s/$SRC/$DST/g" $prec_file
        SRC="Z1F    = 7.0_gpu !4.30_gpu"
        DST="Z1F    = 7.25_gpu"
        sed -i "s/$SRC/$DST/g" $prec_file
    fi

    # BPLASMA code

    if [[ $bplasma == 1 ]]; then

        SRC="nnum   = \[16, -3\] ! apkp - Needs changing"
        DST="nnum   = ["$toroidal_mode"]"
        sed -i "s/$SRC/$DST/g" $prec_file

        SRC="nmde   = 2 ! apkp - Needs changing"
        DST="nmde   = 1"
        sed -i "s/$SRC/$DST/g" $prec_file

        SRC="phase  = \[0.0e0_gpu, 0.0e0_gpu\] ! apkp - Needs changing"
        DST="phase  = [0.0e0_gpu]"
        sed -i "s/$SRC/$DST/g" $prec_file

        # rwm_control and rmp cannot both equal 1
        # Check this:
        if [[ $rwm_control == 1 ]] && [[ $rmp == 1 ]]; then
            echo "rwm_control and rmp cannot both equal 1."
            exit 1
        fi
        if [[ $rwm_control == 1 ]]; then
            BPLASMA_parts=(
                "BPLASMA_cylindrical_tesla_G="$gain_value
                "_bscale="$bscale
                "_n")
        elif [[ $rmp == 1 ]]; then
            BPLASMA_parts=(
                "BPLASMA_"$coil_set
                "_response="$response
                "_current="$(printf "%03d" $current)
                "_100x200"
                "_phase="$(printf "%05.1f" $phase)
                "_n")
        fi
        bplasma_file=$(printf "%s" "${BPLASMA_parts[@]}")
        SRC="'bplasma_file' ! apkp"
        DST="'"$bplasma_file"'"
        sed -i "s|$SRC|$DST|g" $prec_file

        if [[ $rwm_control == 1 ]]; then
            BPLASMA_directory="RWM_control"
        elif [[ $rmp == 1 ]]; then
            BPLASMA_directory="DAVE_v2_100x200/"$coil_set"_n"$toroidal_mode
        fi
        cp -vf \
        $input_dir"/BPLASMA/"$BPLASMA_directory"/"$bplasma_file$toroidal_mode_abs \
        $home_dir"/locust."$tokamak"/InputFiles/."

    fi

    # Read and process pdep_fi.dat file
    pdep_fi_file="$input_dir/ion_info_"$spr_string".dat"
    nion=$(sed -n '1p' "$pdep_fi_file")
    pdep=$(sed -n '2p' "$pdep_fi_file")
    echo $pdep

    # Prepare arrays for fi, Ai, and Zi
    declare -a fi_array=() Ai_array=() Zi_array=()
    while IFS=' ' read -r zi fi ion_name; do
        fi_array+=("$fi")
        Zi_array+=("$zi")
        Ai_array+=("$ion_name") # Directly use the ion name string
    done < <(tail -n +3 "$pdep_fi_file")

    # Convert arrays to comma-separated strings for Fortran
    fi_str=$(IFS=,; echo "${fi_array[*]}")
    Ai_str=$(IFS=,; echo "${Ai_array[*]}")
    Zi_str=$(IFS=,; echo "${Zi_array[*]}")

    # Update prec_file with the new values
    sed -i "s/real( gpu )                                   :: Pdep   = x/real( gpu )                                   :: Pdep   = $pdep/" $prec_file
    sed -i "s/integer,                            parameter :: nion   = x/integer,                            parameter :: nion   = $nion/" $prec_file
    sed -i "s/real( gpu ),  dimension(nion),      parameter :: fi     = \[x\]/real( gpu ),  dimension(nion),      parameter :: fi     = \[$fi_str\]/" $prec_file
    sed -i "s/real( gpu ),  dimension(nion),      parameter :: Ai     = \[x\]/real( gpu ),  dimension(nion),      parameter :: Ai     = \[$Ai_str\]/" $prec_file
    sed -i "s/real( gpu ),  dimension(nion),      parameter :: Zi     = \[x\]/real( gpu ),  dimension(nion),      parameter :: Zi     = \[$Zi_str\]/" $prec_file

    make clean
    make FLAGS="$FLAGS_BASE" -j
    mv -f "locust" "locust_"$run_name"_"$n

done

# Print length of arrays
echo "Length of rcoils: ${#rcoils[@]}"
echo "Length of rcoils_inner: ${#rcoils_inner[@]}"
echo "Length of ncoils: ${#ncoils[@]}"
echo "Length of bripples: ${#bripples[@]}"
echo "Length of toroidal_modes: ${#toroidal_modes[@]}"
echo "Length of coil_sets: ${#coil_sets[@]}"
echo "Length of currents: ${#currents[@]}"
echo "Length of responses: ${#responses[@]}"
echo "Length of phases: ${#phases[@]}"
echo "Length of gain_values: ${#gain_values[@]}"
echo "Length of bscales: ${#bscales[@]}"
echo "Length of bplasmas: ${#bplasmas[@]}"
echo "Length of rwm_controls: ${#rwm_controls[@]}"
echo "Length of rmps: ${#rmps[@]}"
echo "Length of spr_strings: ${#spr_strings[@]}"
echo "Length of unbors: ${#unbors[@]}"
echo "Length of dt0s: ${#dt0s[@]}"
echo "Length of timaxs: ${#timaxs[@]}"
echo "Length of dplots: ${#dplots[@]}"
echo "Length of run_categories: ${#run_categories[@]}"

