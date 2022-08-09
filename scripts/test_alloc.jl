using DrWatson
@quickactivate :NMRMonteCarlo
using Profile

function create_test_suite(temperature = 2.5, Lx=32, Ly=Lx)
    latt = CubicLattice2D(Lx, Ly)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)
    metroparams = MetropolisParameters{Float64}([1/temperature], 2^10, 2^10, 2^10)

    model = CleanAshkinTellerModel( latt.params.Lx, latt.params.Ly, atparams.Jex, atparams.Kex, metroparams.total_measurements )
    return metroparams, model
end

function test()
    metroparams, model = create_test_suite()
    sweep_and_measure!(model, metroparams.βvalues[begin], metroparams, MonteCarloMethods.metropolis_sweep!)
end

test()
Profile.clear_malloc_data()
test()