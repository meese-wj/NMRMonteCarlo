using DrWatson
@quickactivate :NMRMonteCarlo
using Test, BenchmarkTools
using .Lattices: site_index

struct NMRTests end
# Using Lx = Ly = 3 means it is easiest to test the central unit cell
# (site number 5 = (2, 2)).
function create_test_suite(::Type{NMRTests}, Lx = 3, Ly=Lx, Jex=1.0, Kex= 0.)
    latt = CubicLattice2D(Lx, Ly)
    atparams = AshkinTellerParameters(Jex, Kex)
    ham = BasicAshkinTellerHamiltonian(latt, atparams)
    return latt, atparams, ham
end
TestSite(::Type{NMRTests}) = 5

reset_spins!(ham, value) = for (iteration, val) ∈ enumerate(Hamiltonians.IterateByDefault, ham) spins(ham)[iteration] = value end

@info "Beginning SimulatingNMR tests..."
@time @testset "SimulatingNMR" begin

    println("  Testing BaFe2As2 Constructions")
    @time @testset "BaFe2As2 Constructions tests" begin
        latt, atparams, ham = create_test_suite(NMRTests)
        # Set spin state to be uniform +1, +1
        reset_spins!(ham, one(eltype(ham)))

        for construct ∈ [ SimulatingNMR.Out_of_Plane, SimulatingNMR.Easy_Axis_In_Plane, SimulatingNMR.Spin_Orbit_Coupling ]
            for color ∈ (Hamiltonians.AT_sigma, Hamiltonians.AT_tau)
                eval(quote 
                    @testset "$($construct): $($color))" begin
                        bm = @benchmark SimulatingNMR.mag_vector($construct, $ham, TestSite(NMRTests), $color)
                        @test bm.allocs == zero(bm.allocs)
                        @test bm.memory == zero(bm.memory)
                    end
                end)
            end
        end
    end
    println()

    println("  Testing BaFe2As2 OOP & EAIP values")
    @time @testset "BaFe2As2 OOP & EAIP value tests" begin
        OOP_vals =  [ 2.9584, 1.4792, 0.0 ]
        OOP_degeneracies = [4, 8, 4]
        
        EAIP_vals = [ 2.9584, 2.4820, 0.0, 6.9696 ]
        EAIP_degeneracies = [2, 8, 4, 2]

        latt, atparams, ham = create_test_suite(NMRTests)

        function degeneracy_equality(a_arr, b_arr, deg_a, deg_b)
            if size(deg_a) != size(deg_b) return false end
            bool_arr = zeros(Bool, size(deg_a))
            for (idx, val) ∈ enumerate(a_arr)
                bool_arr[idx] = deg_a[idx] == deg_b[ findfirst(≈(val), b_arr) ]
            end
            return all(bool_arr)
        end
        
        @testset "As⁺ atom tests" begin
            reset_spins!(ham, one(eltype(ham)))
            
            Ωval_idx = 1
            # Check (σ0, τ0, σ(1,0), τ(0,1)) = (+1, +1, +1, +1)
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[1]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[1]
            
            # Check (σ0, τ0, σ(1,0), τ(0,1)) = (+1, +1, +1, -1)
            ham[site_index(latt, TestSite(NMRTests), (0, 1)), Hamiltonians.AT_tau] *= -1
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[2]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[2]
            reset_spins!(ham, one(eltype(ham)))
            
            # Check (σ0, τ0, σ(1,0), τ(0,1)) = (+1, -1, +1, -1)
            ham[site_index(latt, TestSite(NMRTests), (0, 0)), Hamiltonians.AT_tau] *= -1
            ham[site_index(latt, TestSite(NMRTests), (0, 1)), Hamiltonians.AT_tau] *= -1
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[1]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[3]
            reset_spins!(ham, one(eltype(ham)))

            # Check (σ0, τ0, σ(1,0), τ(0,1)) = (+1, +1, -1, -1)
            ham[site_index(latt, TestSite(NMRTests), (1, 0)), Hamiltonians.AT_sigma] *= -1
            ham[site_index(latt, TestSite(NMRTests), (0, 1)), Hamiltonians.AT_tau] *= -1
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[3]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[4]
            reset_spins!(ham, one(eltype(ham)))

            # Now test all others...
            test_spins = (one(eltype(ham)), -one(eltype(ham)))
            calc_Ωvals_OOP  = eltype(ham)[]
            calc_Ωvals_EAIP = eltype(ham)[]
            for (σ0, τ0, σ1p0, τ0p1) ∈ Iterators.product(test_spins, test_spins, test_spins, test_spins)
                ham[site_index(latt, TestSite(NMRTests), (0, 0)), Hamiltonians.AT_sigma]  = σ0
                ham[site_index(latt, TestSite(NMRTests), (0, 0)), Hamiltonians.AT_tau]    = τ0
                ham[site_index(latt, TestSite(NMRTests), (1, 0)), Hamiltonians.AT_sigma] = σ1p0
                ham[site_index(latt, TestSite(NMRTests), (0, 1)), Hamiltonians.AT_sigma] = τ0p1
                
                Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
                push!(calc_Ωvals_OOP, round(Ωvals[Ωval_idx], digits = 4))
                Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
                push!(calc_Ωvals_EAIP, round(Ωvals[Ωval_idx], digits = 4))
            end
            unique_OOP = unique(calc_Ωvals_OOP)
            unique_EAIP = unique(calc_Ωvals_EAIP)
            
            counts_OOP  = map( x -> count(==(x), calc_Ωvals_OOP), unique_OOP )
            counts_EAIP = map( x -> count(==(x), calc_Ωvals_EAIP), unique_EAIP )
            @test unique_OOP ⊆ OOP_vals && OOP_vals ⊆ unique_OOP
            @test unique_EAIP ⊆ EAIP_vals && EAIP_vals ⊆ unique_EAIP
            @test degeneracy_equality(OOP_vals, unique_OOP, OOP_degeneracies, counts_OOP)
            @test degeneracy_equality(EAIP_vals, unique_EAIP, EAIP_degeneracies, counts_EAIP)
        end

        @testset "As⁻ atom tests" begin
            reset_spins!(ham, one(eltype(ham)))

            Ωval_idx = 2
            # Check (σ0, τ0, σ(0,-1), τ(-1,0)) = (+1, +1, +1, +1)
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[1]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[1]
            
            # Check (σ0, τ0, σ(0,-1), τ(-1,0)) = (+1, +1, +1, -1)
            ham[site_index(latt, TestSite(NMRTests), (-1, 0)), Hamiltonians.AT_tau] *= -1
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[2]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[2]
            reset_spins!(ham, one(eltype(ham)))
            
            # Check (σ0, τ0, σ(0,-1), τ(-1,0)) = (+1, -1, +1, -1)
            ham[site_index(latt, TestSite(NMRTests), (0, 0)), Hamiltonians.AT_tau] *= -1
            ham[site_index(latt, TestSite(NMRTests), (-1, 0)), Hamiltonians.AT_tau] *= -1
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[1]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[3]
            reset_spins!(ham, one(eltype(ham)))

            # Check (σ0, τ0, σ(0,-1), τ(-1,0)) = (+1, +1, -1, -1)
            ham[site_index(latt, TestSite(NMRTests), (0, -1)), Hamiltonians.AT_sigma] *= -1
            ham[site_index(latt, TestSite(NMRTests), (-1, 0)), Hamiltonians.AT_tau] *= -1
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ OOP_vals[3]
            Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
            @test Ωvals[Ωval_idx] ≈ EAIP_vals[4]
            reset_spins!(ham, one(eltype(ham)))

            # Now test all others...
            test_spins = (one(eltype(ham)), -one(eltype(ham)))
            calc_Ωvals_OOP  = eltype(ham)[]
            calc_Ωvals_EAIP = eltype(ham)[]
            for (σ0, τ0, σ0m1, τ1m0) ∈ Iterators.product(test_spins, test_spins, test_spins, test_spins)
                ham[site_index(latt, TestSite(NMRTests), (0, 0)), Hamiltonians.AT_sigma]  = σ0
                ham[site_index(latt, TestSite(NMRTests), (0, 0)), Hamiltonians.AT_tau]    = τ0
                ham[site_index(latt, TestSite(NMRTests), (0, -1)), Hamiltonians.AT_sigma] = σ0m1
                ham[site_index(latt, TestSite(NMRTests), (-1, 0)), Hamiltonians.AT_sigma] = τ1m0
                
                Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
                push!(calc_Ωvals_OOP, round(Ωvals[Ωval_idx], digits = 4))
                Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Easy_Axis_In_Plane, ham, latt, TestSite(NMRTests))
                push!(calc_Ωvals_EAIP, round(Ωvals[Ωval_idx], digits = 4))
            end
            unique_OOP = unique(calc_Ωvals_OOP)
            unique_EAIP = unique(calc_Ωvals_EAIP)
            
            counts_OOP  = map( x -> count(==(x), calc_Ωvals_OOP), unique_OOP )
            counts_EAIP = map( x -> count(==(x), calc_Ωvals_EAIP), unique_EAIP )
            @test unique_OOP ⊆ OOP_vals && OOP_vals ⊆ unique_OOP
            @test unique_EAIP ⊆ EAIP_vals && EAIP_vals ⊆ unique_EAIP
            @test degeneracy_equality(OOP_vals, unique_OOP, OOP_degeneracies, counts_OOP)
            @test degeneracy_equality(EAIP_vals, unique_EAIP, EAIP_degeneracies, counts_EAIP)
        end

    end
    println()
   
    println("  Total testset timing:")
end
@info "End of SimulatingNMR tests."

