#!/usr/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --job-name="nmr-L-64_%a.out"
#SBATCH -o %x-%j.out
#=
    pwd
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    srun julia --threads=$SLURM_CPUS_PER_TASK clean_nmr_temp_sweep.jl
    exit
=#

using DrWatson; @quickactivate :NMRMonteCarlo
@show Threads.nthreads()
@show names(NMRMonteCarlo)
