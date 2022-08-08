function create_latt_ham(Lx=16, Ly=Lx, Jex=1.0, Kex=0.0)
    latt = CubicLattice2D(Lx, Ly)
    ham = AshkinTellerHamiltonian( latt, AshkinTellerParameters(Jex, Kex) )
    return latt, ham
end

@info "Testing AshkinTellerHamiltonian..."
@testset "AshkinTellerHamiltonian" begin
    
    @time @testset "Structure" begin
        latt, ham = create_latt_ham()
        @test eltype(ham) === Float64
        @test num_sites(ham) == num_sites(latt)
        @test num_DoF(ham) == num_colors(ham) * num_sites(ham)
    end

    @time @testset "Ground State Properties" begin
        size_list = [ (128, 64), (64, 128), (128, 128) ]
        K_list = LinRange(0, 1, 5) 
        for (Lx, Ly) ∈ size_list, Kex ∈ K_list            
            latt, ham = create_latt_ham(Lx, Ly, 1.0, Kex)

            spin_states = (one(eltype(ham)), -one(eltype(ham)))
            for (σval, τval) ∈ Iterators.product(spin_states, spin_states)
                # Set the spin state to one of the ferromagnetic ground states
                for (site, vals) ∈ enumerate(Hamiltonians.IterateBySite, ham)
                    ham[site, Hamiltonians.AT_sigma] = σval
                    ham[site, Hamiltonians.AT_tau] = τval
                end
                @test energy(ham, latt) == -0.5 * Lattices.NN_SQUARE_LATT * ( num_colors(ham) * ham.params.Jex + ham.params.Kex ) * num_sites(ham) 
                # @test round(energy(ham, latt)) == round( -0.5 * Lattices.NN_SQUARE_LATT * ( num_colors(ham) * ham.params.Jex + ham.params.Kex ) * num_sites(ham) )
                @test sum(Hamiltonians.sigma_values(ham)) == σval * num_sites(ham)
                @test sum(Hamiltonians.tau_values(ham)) == τval * num_sites(ham)
                @test sum( idx -> Hamiltonians.site_Baxter(ham, idx), (1:num_sites(ham)) ) == (σval * τval) * num_sites(ham)
            end
        end
    end

    @testset "Iterators" begin
        
        @testset "IterateByDefault" begin
            
        end
        
        @testset "IterateBySite" begin
            
        end
        
        @testset "IterateByDoFType" begin
            
        end

    end

end
@info "End of AshkinTellerHamiltonian tests."