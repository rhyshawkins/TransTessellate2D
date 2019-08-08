#!/bin/bash
#SBATCH --job-name=southamericahmc
#SBATCH --partition=allq
#SBATCH --time=12:00:00
###########################################################
# USER PARAMETERS
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --threads-per-core=1
###########################################################
#SBATCH -o output_%j.txt
#SBATCH -e stderr_%j.txt
#

module purge
module load load userspace/all
module load intel-compiler/64/2018.0.128
module load intel-mpi/64/2018.0.128


### Your script ####

HOME=/quanta1/home/rhys

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib

BASE=$HOME/software/GeneralVoronoiCartesian/
BIN=$BASE/tides/tideshmcpt

#
# Input observations
#
SCRATCH=/scratch/seismology/rhys
DATABASE=$SCRATCH/data/TidesData/
INPUT=$DATABASE/southamerica_sealand.txt

#
# Output dir
#
OUTPUTBASE=$SCRATCH/Inversions/tides
mkdir -p $OUTPUTBASE

OUTPUT=$OUTPUTBASE/southamerica_pt_d_sealand/
mkdir -p $OUTPUT

cd $OUTPUT

#
# Priors
#
echo "prior Uniform" > prior_sea.txt
echo "-10.0 10.0" >> prior_sea.txt
echo "proposal Gaussian" >> prior_sea.txt
echo "0.2" >> prior_sea.txt

echo "prior Uniform" > prior_land.txt
echo "-10.0 10.0" >> prior_land.txt
echo "proposal Gaussian" >> prior_land.txt
echo "0.3" >> prior_land.txt

#
# Position priors
#
echo "positionprior UniformPosition" > positionprior_sea.txt
echo "positionproposal GaussianPosition" >> positionprior_sea.txt
echo "0.5" >> positionprior_sea.txt

echo "positionprior UniformPosition" > positionprior_land.txt
echo "positionproposal GaussianPosition" >> positionprior_land.txt
echo "0.5" >> positionprior_land.txt

#
# Hierarchical priors
#
echo "prior Jeffreys" > hierarchical_prior_tide.txt
echo "0.01 10.0" >> hierarchical_prior_tide.txt
echo "proposal LogGaussian" >> hierarchical_prior_tide.txt
echo "0.1" >> hierarchical_prior_tide.txt

echo "prior Jeffreys" > hierarchical_prior_land.txt
echo "0.01 10.0" >> hierarchical_prior_land.txt
echo "proposal LogGaussian" >> hierarchical_prior_land.txt
echo "0.025" >> hierarchical_prior_land.txt

echo "prior Jeffreys" > hierarchical_prior_sea.txt
echo "0.01 10.0" >> hierarchical_prior_sea.txt
echo "proposal LogGaussian" >> hierarchical_prior_sea.txt
echo "0.1" >> hierarchical_prior_sea.txt

OPTIONS="-i $INPUT \
        -P prior_sea.txt -P prior_land.txt \
        -M positionprior_sea.txt -M positionprior_land.txt \
        -H hierarchical_prior_tide.txt -H hierarchical_prior_land.txt -H hierarchical_prior_sea.txt \
        -x -95.0 -X -25.0 \
        -y -60.0 -Y 15.0 \
        -C 100 \
        -T 500 \
        -t 1000000 \
        -v 10000 \
        -z 10000 \
        -Z 5.0 \
        -c 24 \
        -A 1 -A 1 \
        -K 4 -m 5.0 -e 10 \
        -S 1.0e-2 "


# To run an mpi code with Omnipath network
mpirun -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $OPTIONS

exit
