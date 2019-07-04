#!/bin/bash
# Job name:
#SBATCH --job-name=EMI-tune
#
# Project:
#SBATCH --account=nn9279k
#
# Wall clock limit:
#SBATCH --time=01:00:00
#
# Max memory usage per task:
#SBATCH --mem-per-cpu=8192M
#
# Number of tasks (cores):
#SBATCH --nodes=24

#SBATCH --ntasks-per-node=4

# Set up job environment:
# module load openmpi.intel
source /cluster/bin/jobsetup
source ~johannr/fenics-2018.2.0.dev0+zampini-2019.02.04.abel.intel.conf
# module purge   # clear any inherited modules
set -o errexit # exit on errors

# Set up input and output files:
cp -r . $SCRATCH

chkfile "./results"

cd $SCRATCH
mpirun python3 emi_tune_solver.py
