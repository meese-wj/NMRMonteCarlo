# Simulate a single temperature
using DrWatson
@quickactivate :NMRMonteCarlo
using MonteCarloMeasurementUncertainty

function create_test_suite(temperature = 2.5, Lx=16, Ly=Lx)
    latt = CubicLattice2D(Lx, Ly)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)
    metroparams = MetropolisParameters{Float64}([1/temperature], 2^12, 2^12, 2^12)

    model = CleanAshkinTellerModel( latt.params.Lx, latt.params.Ly, atparams.Jex, atparams.Kex, metroparams.total_measurements )
    return metroparams, model
end

function simulate!( model, metroparams )
    @info "Start of simulation."

    beta = metroparams.βvalues[begin]
    thermalize!(model, beta, metroparams, metropolis_sweep!)
    sweep_and_measure!(model, beta, metroparams, metropolis_sweep!)
    
    for obs ∈ Observables(model)
        @show name(obs)
        @show binning_analysis(obs)
    end

    @info "End of simulation."
end

metroparams, model = create_test_suite()
timer = @timed simulate!(model, metroparams)
println()
@info "$(round( (metroparams.therm_sweeps + metroparams.measure_sweeps) * num_DoF(Hamiltonian(model)) / timer.time, sigdigits = 4) ) updates/second"
@info "$( timer.bytes / 1000 ) KiB allocated"