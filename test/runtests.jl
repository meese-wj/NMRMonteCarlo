using DrWatson
@quickactivate :NMRMonteCarlo

using Test

test_timer = @timed @testset "NMRMonteCarlo" begin
    include(joinpath("Lattices", "latticetests.jl"))
end
@info "Test timing: $(test_timer.time) seconds" 