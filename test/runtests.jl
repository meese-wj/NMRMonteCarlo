using DrWatson
@quickactivate :NMRMonteCarlo

using Test, BenchmarkTools

test_timer = @timed @testset "NMRMonteCarlo" begin
    include(joinpath("Lattices", "latticetests.jl"))
    include(joinpath("Hamiltonians", "test_AshkinTellerHamiltonian.jl"))
    include(joinpath("MonteCarloMethods", "test_MonteCarloMethods.jl"))
    include(joinpath("SimulatingNMR", "test_SimulatingNMR.jl"))
end
@info "Test timing: $(test_timer.time) seconds" 