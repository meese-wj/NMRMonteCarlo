using Pkg
Pkg.activate(@__DIR__)

include("src/ParametersJSON.jl")
include("src/Lattices/CubicLattice2D.jl")
include("src/Ashkin_Teller/AT_Hamiltonian.jl")
include("src/Ashkin_Teller/AT_State.jl")
include("src/Monte_Carlo_Core/MonteCarloCore.jl")

function run_default_metropolis()
    latt = CubicLattice2D()
    ham = AT_Hamiltonian{Float64}(latt)
    mc_params = MetropolisParameters{Float64}()
    mc_model = Model(latt, ham)
    thermalize!(mc_model, mc_params, metropolis_sweep!)
    sweep_and_measure!(mc_model, mc_params, metropolis_sweep!, nothing)
end

using BenchmarkTools
latt = CubicLattice2D()
ham = AT_Hamiltonian{Float64}(latt)
mc_params = MetropolisParameters{Float64}()
mc_model = Model(latt, ham)

# @time run_default_metropolis()
mc_states = build_state_container(typeof(mc_model.ham), num_DoF(ham), mc_params.total_measurements)

@time thermalize!(mc_model, mc_params, metropolis_sweep!)
@time thermalize!(mc_model, mc_params, metropolis_sweep!)

@time sweep_and_measure!(mc_model, mc_params, metropolis_sweep!, mc_states)
@time sweep_and_measure!(mc_model, mc_params, metropolis_sweep!, mc_states)
export_states(mc_states, joinpath(@__DIR__, "NMR_Simulation_Data"))