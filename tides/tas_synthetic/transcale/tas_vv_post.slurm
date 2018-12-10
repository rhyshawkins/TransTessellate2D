#!/bin/bash
#SBATCH --job-name=posttasvv
#SBATCH --partition=transcale
#SBATCH --time=6:00:00
###########################################################
# USER PARAMETERS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
###########################################################
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH -o output_%j.txt
#SBATCH -e stderr_%j.txt
#
# Intel 17.0.4
source /softs/intel/l_ics/2017_update4/compilers_and_libraries_2017.4.196/linux/bin/compilervars.sh intel64
source /softs/intel/l_ics/2017_update4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpivars.sh intel64
export INTEL_LICENSE_FILE=/softs/intel/l_ics/license:$INTEL_LICENSE_FILE

### Your script ####

BASE=/home_nfs/lgltpe/rhys.hawkins/software/GeneralVoronoiCartesian/
LBIN=$BASE/tides/post_likelihood_mpi
BIN=$BASE/tides/post_mean_mpi
TBIN=$BASE/tides/post_mean_tide_mpi

#
# Input observations
#
DATABASE=$BASE/tides/tas_synthetic
INPUT=$DATABASE/syntheticobs_50_50_50_tas.txt

#
# Results dir
#
RESULTSBASE=/erc_transcale/rhys.hawkins/tides/tas_synthetic

RESULTS=$RESULTSBASE/results_vv/

OUTPUT=$RESULTS/post
mkdir -p $OUTPUT

#
# Prior
#
cd $OUTPUT

BASEOPTIONS="-i ../ch.dat \
        -x 144.0 -X 149.0 \
        -y -44.0 -Y -40.0 \
        -T 200 \
        -t 100"

OPTIONS="$BASEOPTIONS \
	-s 500000 \
        -W 200 -H 160"

SEAOPTIONS="$OPTIONS \
        -D stddev_sea.txt \
        -m median_sea.txt \
        -g histogram_sea.txt \
        -z -5.0 -Z 5.0 -b 200"
        
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
mpirun -PSM2 -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $LBIN $BASEOPTIONS -o like.txt -H hierarchical.txt -K khistory.txt
mpirun -PSM2 -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $SEAOPTIONS -I 0 -o mean_sea.txt
mpirun -PSM2 -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $LANDOPTIONS -I 1 -o mean_land.txt
mpirun -PSM2 -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $TBIN $TIDEOPTIONS -o mean_tide.txt


exit
