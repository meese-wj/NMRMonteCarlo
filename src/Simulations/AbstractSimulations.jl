
import ..MonteCarloMethods: thermalize!, sweep_and_measure!, Lattice, Hamiltonian, Observables
import Base: show

export 
#      Base overloads
       show,
#      AbstractSimulations exports
       SimulationParameters, SimulationModel, SimulationMethod, simulate!, Lattice, Hamiltonian, Observables

abstract type AbstractMCMCSimulation end

SimulationParameters(sim::AbstractMCMCSimulation) = sim.parameters
SimulationModel(sim::AbstractMCMCSimulation) = sim.model
Lattice(sim::AbstractMCMCSimulation) = Lattice(SimulationModel(sim))
Hamiltonian(sim::AbstractMCMCSimulation) = Hamiltonian(SimulationModel(sim))
Observables(sim::AbstractMCMCSimulation) = Observables(SimulationModel(sim))
function SimulationMethod(sim::AbstractMCMCSimulation) end

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

function simulate!(sim::AbstractMCMCSimulation)
    @info "Start of simulation."

    params = SimulationParameters(sim)
    model = SimulationModel(sim)

    beta = params.βvalue
    println("\nBeginning thermalization.")
    @time thermalize!(model, beta, params, SimulationMethod(sim))
    println("End of thermalization.")

    println("\nBeginning measurement sweeps.")
    @time sweep_and_measure!(model, beta, params, SimulationMethod(sim))
    println("End of measurement sweeps.")
    
    println("\nBeginning statistical analysis.")
    # result_statistics = analyze_results(model)
    println("End of statistical analysis.\n")

    @info "End of simulation."
    return sim
end