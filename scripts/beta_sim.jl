# Simulate a single temperature
using DrWatson
@quickactivate :NMRMonteCarlo
using MonteCarloMeasurementUncertainty
import OnlineLogBinning: BinningAnalysisResult

function create_test_suite(temperature = 2.5, Lx=32, Ly=Lx)
    latt = CubicLattice2D(Lx, Ly)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)
    metroparams = MetropolisParameters{Float64}([1/temperature], 2^16, 2^16, 2^16)

    model = CleanAshkinTellerModel( latt.params.Lx, latt.params.Ly, atparams.Jex, atparams.Kex, metroparams.total_measurements )
    return metroparams, model
end

function simulate!( model, metroparams )
    @info "Start of simulation."

    beta = metroparams.βvalues[begin]
    println("\nBeginning thermalization.")
    @time thermalize!(model, beta, metroparams, metropolis_sweep!)
    println("End of thermalization.")

    println("\nBeginning measurement sweeps.")
    @time sweep_and_measure!(model, beta, metroparams, metropolis_sweep!)
    println("End of measurement sweeps.")
    
    println("\nBeginning statistical analysis.")
    result_statistics = BinningAnalysisResult[]
    for obs ∈ Observables(model)
        @show name(obs)
        @show res = binning_analysis(obs)
        push!(result_statistics, res)
    end
    println("End of statistical analysis.\n")

    @info "End of simulation."
    return result_statistics
end

metroparams, model = create_test_suite()
timer = @timed results = simulate!(model, metroparams)
@info "Simulation time: $(round(timer.time; sigdigits=4)) seconds"
@info "$(round( (metroparams.therm_sweeps + metroparams.measure_sweeps) * num_DoF(Hamiltonian(model)) / timer.time; sigdigits = 4) ) updates/second"
@info "$( timer.bytes / 1000 ) KiB allocated"