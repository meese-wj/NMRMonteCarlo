# Assumes that @quickactivate "NMR_Monte_Carlo" has been 
# called in the parent scope!

include(srcdir("Lattices", "Lattices.jl"))
using .Lattices


@info "Testing CubicLattice2D..."
@time @testset "CubicLattice2D" begin
    latt = CubicLattice2D( 4, 4 )

    @time @testset "Structure" begin
        @test Lattices.num_sites(latt) == 16
        @test size(latt.neighbors) == (16, 4)
    end

    @time @testset "Neighborhood" begin
        @testset "y = 1" begin
            @test nearest_neighbors(latt, 1) == [13, 4, 2, 5]
            @test nearest_neighbors(latt, 2) == [14, 1, 3, 6]
            @test nearest_neighbors(latt, 3) == [15, 2, 4, 7]
            @test nearest_neighbors(latt, 4) == [16, 3, 1, 8]
        end
        
        @testset "y = 2" begin
            @test nearest_neighbors(latt, 5) == [1, 8, 6, 9]
            @test nearest_neighbors(latt, 6) == [2, 5, 7, 10]
            @test nearest_neighbors(latt, 7) == [3, 6, 8, 11]
            @test nearest_neighbors(latt, 8) == [4, 7, 5, 12]
        end
        
        @testset "y = 3" begin
            @test nearest_neighbors(latt, 9)  == [5, 12, 10, 13]
            @test nearest_neighbors(latt, 10) == [6, 9, 11, 14]
            @test nearest_neighbors(latt, 11) == [7, 10, 12, 15]
            @test nearest_neighbors(latt, 12) == [8, 11, 9, 16]
        end
        
        @testset "y = 4" begin
            @test nearest_neighbors(latt, 13) == [9, 16, 14, 1]
            @test nearest_neighbors(latt, 14) == [10, 13, 15, 2]
            @test nearest_neighbors(latt, 15) == [11, 14, 16, 3]
            @test nearest_neighbors(latt, 16) == [12, 15, 13, 4]
        end
    end
end
@info "End of CubicLattice2D tests."