using DrWatson
@quickactivate :NMRMonteCarlo
using Profile

function create_test_suite()
    latt = CubicLattice2D(64, 64)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)
    metroparams = MetropolisParameters{Float64}([0.3], 1024, 1024, 1024)

    model = CleanAshkinTellerModel( latt.params.Lx, latt.params.Ly, atparams.Jex, atparams.Kex, metroparams.total_measurements )
    return latt, atparams, ham, model
end


function memtest()
    latt, atparams, ham, model = create_test_suite()

    timer = @timed update_observables!(model)
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