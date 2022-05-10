using Pkg
Pkg.activate(@__DIR__)

include("src/ParametersJSON.jl")
include("src/File_Handling/FileHandling.jl")
include("src/Lattices/CubicLattice2D.jl")
include("src/Ashkin_Teller/AT_Hamiltonian.jl")
include("src/Ashkin_Teller/AT_State.jl")
include("src/Monte_Carlo_Core/MonteCarloCore.jl")

struct AT_2DCL_Metro_Params{T <: AbstractFloat} <: SimulationParameters
    latt_params::CubicLattice2DParams
    ham_params::AT_Parameters{T}
    mc_params::MetropolisParameters{T}
end
StructTypes.StructType(::Type{AT_2DCL_Metro_Params{T}}) where {T} = StructTypes.Struct()
num_exports(sim_params::AT_2DCL_Metro_Params) = mc_params.total_measurements

using BenchmarkTools
latt = CubicLattice2D()
ham = AT_Hamiltonian{Float64}(latt)
mc_params = MetropolisParameters{Float64}()
mc_model = Model(latt, ham)

simulation_parameters = AT_2DCL_Metro_Params( latt.params, ham.params, mc_params )

# @time run_default_metropolis()
mc_states = build_state_container(typeof(mc_model.ham), num_DoF(ham), mc_params.total_measurements)

@time thermalize!(mc_model, mc_params, metropolis_sweep!)
@time thermalize!(mc_model, mc_params, metropolis_sweep!)

@time sweep_and_measure!(mc_model, mc_params, metropolis_sweep!, mc_states)
@time sweep_and_measure!(mc_model, mc_params, metropolis_sweep!, mc_states)
# export_states(mc_states, joinpath(@__DIR__, "NMR_Simulation_Data"))

# Export simulation
data_dir = build_data_directory(; model_name = "Ashkin-Teller")
export_simulation(data_dir, mc_states, simulation_parameters, "job-1")