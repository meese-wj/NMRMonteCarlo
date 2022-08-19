#!/usr/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --array=1-50
#SBATCH --job-name=nmr-L-64_%A_%a.out
#SBATCH -o %x-%j.out
#=
    pwd
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    echo
    srun julia --threads=$SLURM_CPUS_PER_TASK clean_nmr_temp_sweep.jl
    exit
=#

using DrWatson; @quickactivate :NMRMonteCarlo

const slurm_arr_length::Int = parse(Int, ENV["SLURM_ARRAY_TASK_COUNT"])

beta_vals = LinRange( 1/3, 1/2, slurm_arr_length )

const my_index = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
const my_beta::Float64 = beta_vals[my_index] 

@info "Setting up a small run first."
sim = CleanNMRATMSimulation(; Lx = 4, βvalue= my_beta )

println(sim)
@timev simulate!(sim)
@info "End of the small simulation."

@info "Starting real run."
sim = CleanNMRATMSimulation(; Lx = 64, βvalue= my_bet, Ntherm = 2^18)
@show sim

@timev simulate!(sim)
@info "End of real run."

datapath = savename("clean_temp_sweep", SimulationParameters(sim), "jld2")
@info "Find the data at: $( datadir(datapath) )"

using JLD2
save_object( datadir(datapath), sim )

DrWatson.safesave( datadir(datapath), sim )


