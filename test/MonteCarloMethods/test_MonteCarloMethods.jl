
function create_test_suite(temperature = 2.5, Lx=32, Ly=Lx)
    latt = CubicLattice2D(Lx, Ly)
    atparams = AshkinTellerParameters(1., 0.)
    ham = AshkinTellerHamiltonian(latt, atparams)
    metroparams = MetropolisParameters{Float64}([1/temperature], 2^16, 2^16, 2^16)

    model = CleanAshkinTellerModel( latt.params.Lx, latt.params.Ly, atparams.Jex, atparams.Kex, metroparams.total_measurements )
    return metroparams, model
end

@info "Testing MonteCarloMethods..."
@time @testset "MonteCarloMethods" begin

    println("  Testing memory allocations")
    @time @testset "Memory allocations" begin
        metroparams, model = create_test_suite()

        bm = @benchmark MonteCarloMethods.metropolis_sweep!($model, $(metroparams.βvalues[begin]))
        @test bm.allocs == zero(bm.allocs)
        @test bm.memory == zero(bm.memory)

        timer = @timed MonteCarloMethods.sweep_and_measure!(model, metroparams.βvalues[begin], metroparams, MonteCarloMethods.metropolis_sweep!)
        timer = @timed MonteCarloMethods.sweep_and_measure!(model, metroparams.βvalues[begin], metroparams, MonteCarloMethods.metropolis_sweep!)
        @test timer.bytes == zero(timer.bytes)
        @test timer.gctime == zero(timer.gctime)
    end
    
    println("  Total testset timing:")
end
@info "End of MonteCarloMethods tests."
