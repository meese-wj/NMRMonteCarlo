using DrWatson
@quickactivate :NMRMonteCarlo

using Test

test_timer = @timed @testset "NMRMonteCarlo" begin
    include(joinpath("Lattices", "latticetests.jl"))
    include(joinpath("Hamiltonians", "test_AshkinTellerHamiltonian.jl"))
end
@info "Test timing: $(test_timer.time) seconds" 