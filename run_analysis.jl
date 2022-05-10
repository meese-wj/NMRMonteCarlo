using Pkg
Pkg.activate(".")

using BenchmarkTools

include("src/ParametersJSON.jl")
include("src/Lattices/CubicLattice2D.jl")
include("src/Ashkin_Teller/AT_Hamiltonian.jl")
include("src/Ashkin_Teller/AT_State.jl")
include("src/Monte_Carlo_Core/MonteCarloCore.jl")
include("src/SimulationParameters.jl")
include("src/File_Handling/FileHandling.jl")

const import_directory = "NMR_Simulation_Data/Ashkin-Teller/May-10-2022/"

sim_params, mc_states = import_simulation(import_directory, AT_2DCL_Metro_Params{Float64}, "job-1")

