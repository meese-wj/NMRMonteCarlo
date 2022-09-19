
import ..MonteCarloMethods: thermalize!, sweep_and_measure!, Lattice, Hamiltonian, Observables
import ..Models: collect_hyperfine_susceptibilites
import Base: show
using MonteCarloMeasurementUncertainty
import OnlineLogBinning: BinningAnalysisResult

export 
#      Base overloads
       show,
#      NMRMonteCarlo Convenience overloads
       collect_hyperfine_susceptibilites,
#      AbstractSimulations exports
       SimulationParameters, SimulationModel, SimulationMethod, simulate!, Lattice, Hamiltonian, Observables, analyze

"""
    abstract type AbstractMCMCSimulation end

Supertype for all simulations that will be run.

# Required fields

There are only two required fields for a `<: AbstractMCMCSimulation`:
    
1. `parameters`: Encoding all information needed to run a particular simulation. Must be `<:` [`AbstractMonteCarloParameters`](@ref).
1. `model`: The actual `<:` [`AbstractModel`](@ref) to be [`simulate!`](@ref)d.

# Required Interface Methods

* `Base.iterate`  ← this method must be overloaded for each simulation's observables
* [`SimulationMethod`](@ref)

# Default Interface Methods

## Base overloads

* [`show`](@ref)

## [`NMRMonteCarlo`](@ref) definitions

* [`SimulationParameters`](@ref)
* [`SimulationModel`](@ref)
* [`SimulationMethod`](@ref)
* [`Lattice`](@ref) (forwarded from [`MonteCarloMethods`](@ref))
* [`Hamiltonian`](@ref) (forwarded from [`MonteCarloMethods`](@ref))
* [`Observables`](@ref) (forwarded from [`MonteCarloMethods`](@ref))
* [`simulate!`](@ref)
* [`analyze`](@ref)

"""
abstract type AbstractMCMCSimulation end

"""
    SimulationParameters(::AbstractMCMCSimulation) = sim.parameters

Extract the parameters for a given `<:` [`AbstractMCMCSimulation`](@ref).
"""
SimulationParameters(sim::AbstractMCMCSimulation) = sim.parameters
"""
    SimulationModel(::AbstractMCMCSimulation) = sim.model

Extract the model for a given `<:` [`AbstractMCMCSimulation`](@ref).
"""
SimulationModel(sim::AbstractMCMCSimulation) = sim.model
"""
    Lattice(sim::AbstractMCMCSimulation)

Convenience wrapper around `Lattice(SimulationModel(sim))`.
"""
Lattice(sim::AbstractMCMCSimulation) = Lattice(SimulationModel(sim))
"""
    Hamiltonian(sim::AbstractMCMCSimulation)

Convenience wrapper around `Hamiltonian(SimulationModel(sim))`.
"""
Hamiltonian(sim::AbstractMCMCSimulation) = Hamiltonian(SimulationModel(sim))
"""
    Observables(sim::AbstractMCMCSimulation)

Convenience wrapper around `Observables(SimulationModel(sim))`.
"""
Observables(sim::AbstractMCMCSimulation) = Observables(SimulationModel(sim))
"""
    SimulationMethod(::AbstractMCMCSimulation)

Returns a function method that implements a specific Monte Carlo sweep technique 
to be used in a given `<:` [`AbstractMCMCSimulation`](@ref). Returns `nothing` 
if no method for a simulation subtype has been explicitly implemented.
"""
function SimulationMethod(sim::AbstractMCMCSimulation) end

"""
    show([io::IO = stdout], ::AbstractMCMCSimulation)

`Base` overload to display [`AbstractMCMCSimulation`](@ref)s.
"""
function show(io::IO, sim::AbstractMCMCSimulation)
    tab = "  "
    params = SimulationParameters(sim)
    model  = SimulationModel(sim)
    println(io, "$(typeof(sim)):")
    println()
    println(io, tab, "Simulation Parameters::$(typeof(params))" )
    for field ∈ fieldnames(typeof(params))
        println(io, tab * tab, String(field) * "::$(getfield(params, field) |> typeof) = $(getfield(params, field))")
    end
    println()
    println(io, tab, "Simulation Model::$(typeof(model))")
    for field ∈ fieldnames(typeof(model))
        println(io, tab * tab, String(field) * "::$( getfield(model, field) |> typeof )")
    end
    println()
    println(io, tab, "Simulation Method: $(SimulationMethod(sim))")
end
show(sim::AbstractMCMCSimulation) = show(stdout, sim)

"""
    simulate!(::AbstractMCMCSimulation, [thermalization = true])

Default simulation that will may [`thermalize!`](@ref) and *will* perform
Monte Carlo sweeps according to the [`SimulationMethod`](@ref) rule.
"""
function simulate!(sim::AbstractMCMCSimulation, thermalization = true)
    @info "Start of simulation."

    params = SimulationParameters(sim)
    model = SimulationModel(sim)

    beta = params.βvalue
    if thermalization
        println("\nBeginning thermalization.")
        @time thermalize!(model, beta, params, SimulationMethod(sim))
        println("End of thermalization.")
    else
        println("\nNo thermalization to be performed.")
    end

    println("\nBeginning measurement sweeps.")
    @time sweep_and_measure!(model, beta, params, SimulationMethod(sim))
    println("End of measurement sweeps.")
    
    println("\nBeginning statistical analysis.")
    # result_statistics = analyze_results(model)
    println("End of statistical analysis.\n")

    @info "End of simulation."
    return sim
end
"""
    analyze(::AbstractMCMCSimulation)

Run a `binning_analysis` on each of the [`Observables`](@ref) for a given
`sim <:` [`AbstractMCMCSimulation`](@ref).

Return **all** resulting `BinningAnalysisResult`s as a **single** `Vector`.
It is up to the user to separate the indices appropriately for the different
types of observables.

!!! note
    `binning_analysis` is taken from `MonteCarloMeasurementUncertainty.jl` and
    `BinningAnalysisResult` comes from `OnlineLogBinning.jl`.
"""
function analyze(sim::AbstractMCMCSimulation)
    results = BinningAnalysisResult{eltype(Hamiltonian(sim))}[]
    for obs ∈ Observables(sim)
        push!(results, binning_analysis(obs))
    end
    return results
end

collect_hyperfine_susceptibilites(sim::AbstractMCMCSimulation) = collect_hyperfine_susceptibilites(SimulationModel(sim))