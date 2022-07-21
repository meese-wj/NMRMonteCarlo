using DrWatson
@quickactivate "NMR_Monte_Carlo"
include(srcdir("Lattices", "Lattices.jl"))

using Test

@time @testset "NMR_Monte_Carlo" begin
    
    @time @testset "CubicLattice2D" begin
        latt = CubicLattice2D( 4, 4 )
        @test size(latt.neighbors) == (16, 4)
        @test nearest_neighbors(latt, 1) == [13, 4, 2, 5]
    end
    
end