using Pkg
Pkg.activate(".")

using BenchmarkTools, Plots, LaTeXStrings

include("src/ParametersJSON.jl")
include("src/Lattices/CubicLattice2D.jl")
include("src/Ashkin_Teller/AT_Hamiltonian.jl")
include("src/Ashkin_Teller/AT_State.jl")
include("src/Monte_Carlo_Core/MonteCarloCore.jl")
include("src/SimulationParameters.jl")
include("src/File_Handling/FileHandling.jl")

const import_directory = "NMR_Simulation_Data/Ashkin-Teller/May-10-2022/"

sim_params, mc_states = import_simulation(import_directory, AT_2DCL_Metro_Params{Float64}, "job-1")
latt = (reciprocal_type(sim_params.latt_params))(sim_params.latt_params)
ham  = (reciprocal_type(sim_params.ham_params))(latt, sim_params.ham_params)

(num_dof, num_states) = size(mc_states)
num_AT_sites = num_dof / NUM_AT_COLORS

data_t = eltype(mc_states)
sweep_indices = sweeps_per_export(sim_params.mc_params) .* (1:1:num_states)

# Energies, sigma and tau mags
energies = zeros(data_t, num_states)
mags = zeros(NUM_AT_COLORS, num_states)
@inbounds for idx ∈ eachindex(1:num_states)
    energies[idx] = AT_total_energy(mc_states[:, idx], ham.params, latt) / num_AT_sites
    mags[1, idx] = sum(mc_states[begin:2:end, idx]) / num_AT_sites
    mags[2, idx] = sum(mc_states[2:2:end, idx]) / num_AT_sites
end

@show mags[1,:] == mags[2,:]

plot( sweep_indices, energies,
      color=:blue, xlabel = "MC Sweeps",
      ylabel = "Energy / Site", 
      label="\$\\beta = $(sim_params.mc_params.β)\\,J^{-1}\$")

sigma_plt = plot( sweep_indices, mags[1, :], 
                  color=:blue,
                  ylabel = L"$N^{-1}\sum_i \sigma_i$", 
                  label=nothing)

tau_plt   = plot( sweep_indices, mags[2, :], 
                  color=:orange, xlabel = "MC Sweeps",
                  ylabel = L"$N^{-1}\sum_i \tau_i$", 
                  label=nothing)


plot(sigma_plt, tau_plt, layout=(2,1), link = :xy,
     plot_title="\$\\beta = $(sim_params.mc_params.β)\\,J^{-1}\$")
