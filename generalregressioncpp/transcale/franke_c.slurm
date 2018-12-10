#!/bin/bash
#SBATCH --job-name=frankec
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
BIN=$BASE/generalregressioncpp/hmcpt

#
# Input observations
#
DATABASE=$BASE/generalregressioncpp/example_franke_hmc_pt
INPUT=$DATABASE/syntheticobs_franke.txt

#
# Output dir
#
OUTPUTBASE=/erc_transcale/rhys.hawkins/regression/

mkdir -p $OUTPUTBASE

OUTPUT=$OUTPUTBASE/results_frankec/

mkdir -p $OUTPUT

cd $OUTPUT

#
# Priors
#
echo "prior Uniform" > prior.txt
echo "-0.5 2.0" >> prior.txt
echo "proposal Gaussian" >> prior.txt
echo "0.05" >> prior.txt

#
# Position priors
#
echo "positionprior UniformPosition" > positionprior.txt
echo "positionproposal GaussianPosition" >> positionprior.txt
echo "0.1" >> positionprior.txt

#
# Hierarchical priors
#
echo "prior Uniform" > hierarchical_prior.txt
echo "0.5 5.0" >> hierarchical_prior.txt
echo "proposal Gaussian" >> hierarchical_prior.txt
echo "0.07" >> hierarchical_prior.txt

OPTIONS="-i $INPUT \
        -P prior.txt \
        -M positionprior.txt \
        -H hierarchical_prior.txt \
        -x 0.0 -X 1.0 \
        -y 0.0 -Y 1.0 \
        -A 2 \
        -C 20 \
        -T 200 \
        -t 1000000 \
        -v 10000 \
        -z 10000 \
        -Z 5.0 \
        -c 28 \
        -K 4 -m 5.0 -e 10"


# To run an mpi code with Omnipath network
mpirun -PSM2 -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $OPTIONS

exit
