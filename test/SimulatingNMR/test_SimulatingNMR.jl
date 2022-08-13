using DrWatson
@quickactivate :NMRMonteCarlo
using Test

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

@info "Beginning SimulatingNMR tests..."
@time @testset "SimulatingNMR" begin

    println("  Testing BaFe2As2 Out_of_Plane values")
    @time @testset "Out_of_Plane tests" begin
        latt, atparams, ham = create_test_suite(NMRTests)
        # Set spin state to be uniform +1, +1
        for (iteration, val) ∈ enumerate(Hamiltonians.IterateByDefault, ham) spins(ham)[iteration] = one(eltype(ham)) end
        
        Ωvals = SimulatingNMR.inst_hyperfine_fluctuations(SimulatingNMR.Out_of_Plane, ham, latt, TestSite(NMRTests))
        @show Ωvals
    end
    println()
   
    println("  Total testset timing:")
end
@info "End of SimulatingNMR tests."

