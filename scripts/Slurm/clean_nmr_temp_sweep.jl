#!/usr/bin/bash -l
#SBATCH --time=17:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --array=1-50
#SBATCH --job-name=nmr-L-32_%A_%a.out
#SBATCH -o %x-%j.out
#=
    pwd
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    echo
    srun julia --threads=$SLURM_CPUS_PER_TASK clean_nmr_temp_sweep.jl
    exit
=#
using Pkg
using DrWatson
@quickactivate :NMRMonteCarlo
# Pkg.instantiate()

const slurm_arr_length::Int = parse(Int, ENV["SLURM_ARRAY_TASK_COUNT"])

@show const Lvalue = 32
@show const Jex = 1.0
@show const Kex = 0.5
@show const Tc = critical_temperature(Jex, Kex)
@show const βc = 1 / Tc
const dTLow = 0.125 * Tc
const dTHigh = 0.375 * Tc
@show const betaHigh = 1 / (Tc - dTLow)
@show const betaLow = 1 / (Tc + dTHigh)

beta_vals = LinRange( betaLow, betaHigh, slurm_arr_length )

const my_index = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
const my_beta::Float64 = beta_vals[my_index] 
const nmr_type::Type = Out_of_Plane

@info "Setting up a small run first."
sim = CleanNMRATMSimulation( Lx = 4, βvalue = my_beta, nmr_spin_type = nmr_type)

println(sim)
@timev simulate!(sim)
@info "End of the small simulation."

@info "Starting real run."
save_sim = CleanNMRATMSimulation(Jex = Jex, Kex = Kex, 
                            Lx = Lvalue, βvalue= my_beta,
                            Ntherm = 2^20, Nmeas = 2^18, 
                            nmr_spin_type = nmr_type)
println(save_sim)                            
@show const nreplicas = Threads.nthreads()
all_local_susc = zeros(Float64, nreplicas, 2 * Lvalue^2)
Threads.@threads for rep_dx ∈ 1:nreplicas
    sim = CleanNMRATMSimulation(Jex = Jex, Kex = Kex, 
                                Lx = Lvalue, βvalue= my_beta,
                                Ntherm = 2^20, Nmeas = 2^18, 
                                nmr_spin_type = nmr_type)
    @show rep_dx, sim

    @timev simulate!(sim)
    rep_local_chi_vals = collect_hyperfine_susceptibilites(sim)
    for atom_idx ∈ eachindex(all_local_susc[rep_dx, :])
        all_local_susc[rep_dx, atom_idx] = rep_local_chi_vals[rep_dx].val
    end
end
@info "End of real run."

flush(stdout)

using Measurements, Statistics
local_chi_vals = similar(all_local_susc[begin, :], Measurement{Float64})
for atom_idx ∈ eachindex(all_local_susc[begin, :])
    mean_χ = mean( @view all_local_susc[:, atom_idx] )
    err_χ  = std( @view all_local_susc[:, atom_idx]; mean = mean_χ )
    local_chi_vals[atom_idx] = measurement(mean_χ, err_χ)
end
@show local_chi_vals
println()

chi_path = savename("threaded_hyperfine_susceptibilites_$(nmr_type)", SimulationParameters(save_sim), "jld2")
datapath = savename("clean_temp_sweep_$(nmr_type)", SimulationParameters(save_sim), "jld2")
@info "Find the data at: $( datadir(chi_path) )"
# @info "Find the data at: $( datadir(datapath) )"

using JLD2
save_object( datadir(chi_path), local_chi_vals )
# save_object( datadir(datapath), save_sim )

DrWatson.safesave( datadir(chi_path), local_chi_vals )
# DrWatson.safesave( datadir(datapath), save_sim )


