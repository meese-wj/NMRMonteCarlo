using DrWatson
@quickactivate :NMRMonteCarlo
using Profile

struct ATModel <: MonteCarloMethods.AbstractModel
    lattice::CubicLattice2D
    hamiltonian::AshkinTellerHamiltonian{Float64}
end

function create_test_suite()
    latt = CubicLattice2D(64, 64)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)

    model = ATModel(latt, ham)
    return latt, atparams, ham, model
end


function memtest()
    latt, atparams, ham, model = create_test_suite()

    timer = @timed MonteCarloMethods.metropolis_update!(model, 0.5, 5)
    timer = @timed MonteCarloMethods.metropolis_sweep!(model, 0.5)
    println("$(timer.time) seconds")
    println("$(timer.bytes) bytes allocated")
    println("$(round(num_DoF(ham) / timer.time; sigdigits = 4)) updates per second")
end

# let 
#     latt, atparams, ham, model = create_test_suite()
#     @code_warntype Hamiltonians.DoF_energy_change( ham, latt, 5, Hamiltonians.AT_sigma )
# end

# Compile everything
memtest()
# Clear the compilation memory
Profile.clear_malloc_data()
# Test the memory
memtest()