#!/bin/bash 
# caller_arrays.sbatch
#
#SBATCH -J caller # A single job name for the array
#SBATCH -n 1 # Number of cores requested 
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 450 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to 
#SBATCH --mem=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append 
#SBATCH -o ./outputs/caller_arrays_%A_%a.out # Standard out goes to this file
#SBATCH -e ./outputs/caller_arrays_%A_%a.err # Standard err goes to this filehostname hostname

module load Anaconda/5.0.1-fasrc02

python /n/holyscratch01/murray_lab/fbarber/20_growth_rate_project/200326_revision_sizer_no_penalty_model/simulation_caller.py