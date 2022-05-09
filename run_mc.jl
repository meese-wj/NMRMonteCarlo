using Pkg
Pkg.activate(@__DIR__)

include("src/ParametersJSON.jl")
include("src/Lattices/CubicLattice2D.jl")
include("src/Ashkin_Teller/AT_Hamiltonian.jl")
include("src/Monte_Carlo_Core/MonteCarloCore.jl")

using BenchmarkTools
ham_file = "defined_AT_Params.jl"
latt = CubicLattice2D()
ham = AT_Hamiltonian{Float64}(latt)
mc_params = MetropolisParameters{Float64}()
mc_model = Model(latt, ham)

@time thermalize!(mc_model, mc_params, metropolis_sweep!)
@time thermalize!(mc_model, mc_params, metropolis_sweep!)

@time sweep_and_measure!(mc_model, mc_params, metropolis_sweep!, nothing)
@time sweep_and_measure!(mc_model, mc_params, metropolis_sweep!, nothing)