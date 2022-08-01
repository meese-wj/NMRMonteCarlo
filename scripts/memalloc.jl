using DrWatson
@quickactivate :NMRMonteCarlo
using Profile

struct ATModel <: MonteCarloMethods.AbstractModel
    lattice::CubicLattice2D
    hamiltonian::AshkinTellerHamiltonian{Float64}
end

function memtest()
    latt = CubicLattice2D(64, 64)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)

    model = ATModel(latt, ham)

    MonteCarloMethods.metropolis_update!(model, 0.5, 5)
    timer = @timed MonteCarloMethods.metropolis_sweep!(model, 0.5)
    println("$(timer.time) seconds")
    println("$(timer.bytes) bytes allocated")
    println("$(round(num_DoF(ham) / timer.time; sigdigits = 4)) updates per second")
end

# Compile everything
memtest()
# Clear the compilation memory
Profile.clear_malloc_data()
# Test the memory
memtest()