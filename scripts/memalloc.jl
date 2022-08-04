using DrWatson
@quickactivate :NMRMonteCarlo
using Profile

function create_test_suite()
    latt = CubicLattice2D(2*8192, 2*8192)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)
    metroparams = MetropolisParameters{Float64}([0.3], 1024, 1024, 1024)

    model = CleanAshkinTellerModel( latt.params.Lx, latt.params.Ly, atparams.Jex, atparams.Kex, metroparams.total_measurements )
    return metroparams, model
end


function memtest()
    metroparams, model = create_test_suite()

    timer = @timed metropolis_sweep!(model, metroparams.βvalues[begin])
    println("$(timer.time) seconds")
    println("$(timer.bytes) bytes allocated")
    println("$(round(num_DoF(Hamiltonian(model)) / timer.time; sigdigits = 4)) updates per second")
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