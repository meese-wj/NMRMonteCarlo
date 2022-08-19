#!/usr/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=10
#SBATCH --mem=2g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH -o %x-%j.out
#=
    srun julia $(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
    exit
# =#

using DrWatson; @quickactivate :NMRMonteCarlo
names(NMRMonteCarlo)