#!/usr/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=10
#SBATCH --mem=2g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH -o %x-%j.out
#=
    pwd
    srun julia clean_nmr_temp_sweep.jl
    exit
=#

using DrWatson; @quickactivate :NMRMonteCarlo
names(NMRMonteCarlo)