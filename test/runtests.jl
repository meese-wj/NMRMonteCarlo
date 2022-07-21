using DrWatson
@quickactivate "NMR_Monte_Carlo"

using Test

test_timer = @timed @testset "NMR_Monte_Carlo" begin
    include(joinpath("Lattices", "latticetests.jl"))
end
@info "Test timing: $(test_timer.time) seconds" 