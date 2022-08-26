# Simulate a single temperature
using DrWatson
@quickactivate :NMRMonteCarlo
using MonteCarloMeasurementUncertainty
import OnlineLogBinning: BinningAnalysisResult

function create_test_suite(temperature = 2.269, Lx=8, Ly=Lx)
    latt = CubicLattice2D(Lx, Ly)
    atparams = AshkinTellerParameters(1., 0.)
    ham = BasicAshkinTellerHamiltonian(latt, atparams)
    metroparams = MetropolisParameters{Float64}([1/temperature], 2^16, 2^20, 2^16)

    model = CleanNMRAshkinTellerModel( latt.params.Lx, latt.params.Ly, atparams.Jex, atparams.Kex, metroparams.total_measurements, SimulatingNMR.Out_of_Plane )
    return metroparams, model
end

function analyze_results(model)
    time_statistics = BinningAnalysisResult[]
    accum_statistics = BinningAnalysisResult[]
    for obstype ∈ (Models.TimeSeriesObservables, Models.AccumulatedSeriesObservables), obs ∈ obstype(model)
        @show name(obs)
        @show res = binning_analysis(obs)
        obstype === Models.TimeSeriesObservables ? push!(time_statistics, res) : push!(accum_statistics, res)
    end
    return time_statistics, accum_statistics
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
    result_statistics = analyze_results(model)
    println("End of statistical analysis.\n")

    @info "End of simulation."
    return result_statistics
end

metroparams, model = create_test_suite()
timer = @timed results = simulate!(model, metroparams)
@info "Simulation time: $(round(timer.time; sigdigits=4)) seconds"
@info "$(round( (metroparams.therm_sweeps + metroparams.measure_sweeps) * num_DoF(Hamiltonian(model)) / timer.time; sigdigits = 4) ) updates/second"
@info "$( timer.bytes / 1000 ) KiB allocated"