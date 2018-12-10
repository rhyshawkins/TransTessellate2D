#!/bin/bash
#SBATCH --job-name=tascv
#SBATCH --partition=transcale
#SBATCH --time=6:00:00
###########################################################
# USER PARAMETERS
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-socket=14
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
BIN=$BASE/tides/tideshmcpt

#
# Input observations
#
DATABASE=$BASE/tides/tas_synthetic
INPUT=$DATABASE/syntheticobs_50_50_50_tas.txt

#
# Output dir
#
OUTPUTBASE=/erc_transcale/rhys.hawkins/tides/tas_synthetic
mkdir -p $OUTPUTBASE

OUTPUT=$OUTPUTBASE/results_cv/

mkdir -p $OUTPUT

cd $OUTPUT

#
# Priors
#
echo "prior Uniform" > prior_sea.txt
echo "-5.0 5.0" >> prior_sea.txt
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
echo "0.3" >> positionprior_sea.txt

echo "positionprior UniformPosition" > positionprior_land.txt
echo "positionproposal GaussianPosition" >> positionprior_land.txt
echo "0.2" >> positionprior_land.txt

#
# Hierarchical priors
#
echo "prior Uniform" > hierarchical_prior_tide.txt
echo "0.1 10.0" >> hierarchical_prior_tide.txt
echo "proposal Gaussian" >> hierarchical_prior_tide.txt
echo "0.1" >> hierarchical_prior_tide.txt

echo "prior Uniform" > hierarchical_prior_land.txt
echo "0.1 10.0" >> hierarchical_prior_land.txt
echo "proposal Gaussian" >> hierarchical_prior_land.txt
echo "0.04" >> hierarchical_prior_land.txt

echo "prior Uniform" > hierarchical_prior_sea.txt
echo "0.1 10.0" >> hierarchical_prior_sea.txt
echo "proposal Gaussian" >> hierarchical_prior_sea.txt
echo "0.07" >> hierarchical_prior_sea.txt

OPTIONS="-i $INPUT \
        -P prior_sea.txt -P prior_land.txt \
        -M positionprior_sea.txt -M positionprior_land.txt \
        -H hierarchical_prior_tide.txt -H hierarchical_prior_land.txt -H hierarchical_prior_sea.txt \
        -x 144.0 -X 149.0 \
        -y -44.0 -Y -40.0 \
        -C 20 \
        -T 200 \
        -t 1000000 \
        -v 10000 \
        -z 10000 \
        -Z 5.0 \
        -c 28 \
        -A 2 -A 0 \
        -K 4 -m 5.0 -e 10"


# To run an mpi code with Omnipath network
mpirun -PSM2 -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $OPTIONS

exit
