function create_latt_ham(Lx=16, Ly=Lx, Jex=1.0, Kex=0.0)
    latt = CubicLattice2D(Lx, Ly)
    ham = BasicAshkinTellerHamiltonian( latt, AshkinTellerParameters(Jex, Kex) )
    return latt, ham
end

@info "Testing BasicAshkinTellerHamiltonian..."
@time @testset "BasicAshkinTellerHamiltonian" begin
    
    println("  Testing Structure")
    @time @testset "Structure" begin
        latt, ham = create_latt_ham()
        @test eltype(ham) === Float64
        @test num_sites(ham) == num_sites(latt)
        @test num_DoF(ham) == num_colors(ham) * num_sites(ham)
    end
    println()

    println("  Testing Ground State Properties")
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
    println()

    println("  Testing Iterators")
    @time @testset "Iterators" begin
        
        @inline _alloctest_Iterator(itertype, ham) = for tup ∈ enumerate( itertype, ham ) end

        @time @testset "IterateByDefault" begin
            latt, ham = create_latt_ham()
            # First, check that the iterator is non-allocating
            bm = @benchmark $_alloctest_Iterator( Hamiltonians.IterateByDefault, $ham)
            @test bm.allocs == zero(bm.allocs)
            @test bm.memory == zero(bm.memory)
            # Second, check that its indices and values are correct.
            indices = Int[]
            sites = Int[]
            values = eltype(ham)[]
            for (idx, site_val) ∈ enumerate( Hamiltonians.IterateByDefault, ham ) 
                push!(indices, idx)
                push!(sites, site_val[1])
                push!(values, site_val[2])
            end
            @test all( indices .== collect(1:num_DoF(ham)) )
            @test all( sites .== [ Hamiltonians.site_index(ham, dof_idx) for dof_idx ∈ 1:num_DoF(ham) ] )
            @test all( values .== spins(ham) )
        end
        
        @time @testset "IterateBySite" begin
            latt, ham = create_latt_ham()
            # First, check that the iterator is non-allocating
            bm = @benchmark $_alloctest_Iterator( Hamiltonians.IterateBySite, $ham)
            @test bm.allocs == zero(bm.allocs)
            @test bm.memory == zero(bm.memory)
            # Second, check that its iterations and values are correct.
            iterations = Int[]
            sites = Int[]
            σvals = eltype(ham)[]
            τvals = eltype(ham)[]
            baxvals = eltype(ham)[]
            for (idx, site_vals) ∈ enumerate( Hamiltonians.IterateBySite, ham ) 
                push!(iterations, idx)
                push!(sites, site_vals[1])
                push!(σvals, site_vals[2][1])
                push!(τvals, site_vals[2][2])
                push!(baxvals, site_vals[2][3])
            end
            @test all( iterations .== collect(1:num_sites(ham)) )
            @test all( sites .== collect(1:num_sites(ham)) )
            @test all( σvals .== Hamiltonians.sigma_values(ham) )
            @test all( τvals .== Hamiltonians.tau_values(ham) )
            @test all( baxvals .== [ Hamiltonians.site_Baxter(ham, site) for site ∈ sites ] )
        end
        
        @testset "IterateByDoFType" begin
            latt, ham = create_latt_ham()
            # First, check that the iterator is non-allocating
            bm = @benchmark $_alloctest_Iterator( Hamiltonians.IterateByDoFType, $ham)
            @test bm.allocs == zero(bm.allocs)
            @test bm.memory == zero(bm.memory)
            # Second, check that its iterations and values are correct.
            iterations = Int[]
            sites = Int[]
            σvals = eltype(ham)[]
            τvals = eltype(ham)[]
            for (idx, site_vals) ∈ enumerate( Hamiltonians.IterateByDoFType, ham ) 
                push!(iterations, idx)
                push!(sites, site_vals[1])
                if idx <= num_sites(ham)
                    push!(σvals, site_vals[2])
                else
                    push!(τvals, site_vals[2])
                end
            end
            @test all( iterations .== collect(1:num_DoF(ham)) )
            @test all( sites .== collect( Iterators.flatten((1:num_sites(ham), 1:num_sites(ham))) ) )
            @test all( σvals .== Hamiltonians.sigma_values(ham) )
            @test all( τvals .== Hamiltonians.tau_values(ham) )
        end

    end
    println()

    println("  Total testset timing:")
end
@info "End of BasicAshkinTellerHamiltonian tests."