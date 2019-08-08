#!/bin/bash
#SBATCH --job-name=pafrica
#SBATCH --partition=allq
#SBATCH --time=6:00:00
###########################################################
# USER PARAMETERS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
###########################################################
#SBATCH -o output_%j.txt
#SBATCH -e stderr_%j.txt

module purge
module load userspace/all
module load intel-compiler/64/2018.0.128
module load intel-mpi/64/2018.0.128

### Your script ####

HOME=/quanta1/home/rhys

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib

BASE=$HOME/software/GeneralVoronoiCartesian/

LBIN=$BASE/tides/post_likelihood_mpi
BIN=$BASE/tides/post_mean_mpi
TBIN=$BASE/tides/post_mean_tide_mpi
RBIN=$BASE/tides/post_resolution_mpi

#
# Input observations
#
SCRATCH=/scratch/seismology/rhys
DATABASE=$SCRATCH/data/TidesData/
INPUT=$DATABASE/africa_combined.txt

#
# Results dir
#
RESULTSBASE=$SCRATCH/Inversions/tides

RESULTS=$RESULTSBASE/africa_pt_d/

OUTPUT=$RESULTS/post
mkdir -p $OUTPUT

#
# Prior
#
cd $OUTPUT

BASEOPTIONS="-i ../ch.dat \
        -A 1 -A 1 \
        -x -30.0 -X 60.0 \
        -y -40.0 -Y 40.0 \
        -T 500 \
        -t 10"

OPTIONS="$BASEOPTIONS \
        -s 500000 \
        -W 200 -H 200"

SEAOPTIONS="$OPTIONS \
        -D stddev_sea.txt \
        -m median_sea.txt \
        -g histogram_sea.txt \
        -z -5.0 -Z 5.0 -b 400"
        
LANDOPTIONS="$OPTIONS \
        -D stddev_land.txt \
        -m median_land.txt \
        -g histogram_land.txt \
        -z -10.0 -Z 10.0 -b 400"

TIDEOPTIONS="$OPTIONS \
        -D stddev_tide.txt \
        -m median_tide.txt \
        -g histogram_tide.txt \
        -z -10.0 -Z 10.0 -b 400"

# To run an mpi code with Omnipath network
mpirun -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $LBIN $BASEOPTIONS -o like.txt -H hierarchical.txt -K khistory.txt
mpirun -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $SEAOPTIONS -I 0 -o mean_sea.txt
mpirun -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $LANDOPTIONS -I 1 -o mean_land.txt
mpirun -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $TBIN $TIDEOPTIONS -o mean_tide.txt
mpirun -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $RBIN $OPTIONS -I 0 -o res_mean_sea.txt -D res_std_sea.txt
mpirun -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $RBIN $OPTIONS -I 1 -o res_mean_land.txt -D res_std_land.txt


exit
