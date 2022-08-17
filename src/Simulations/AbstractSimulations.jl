
import ..MonteCarloMethods: thermalize!, sweep_and_measure!, Lattice, Hamiltonian, Observables

export SimulationParameters, SimulationModel, SimulationMethod, simulate!, Lattice, Hamiltonian, Observables

abstract type AbstractMCMCSimulation end

SimulationParameters(sim::AbstractMCMCSimulation) = sim.parameters
SimulationModel(sim::AbstractMCMCSimulation) = sim.model
Lattice(sim::AbstractMCMCSimulation) = Lattice(SimulationModel(sim))
Hamiltonian(sim::AbstractMCMCSimulation) = Hamiltonian(SimulationModel(sim))
Observables(sim::AbstractMCMCSimulation) = Observables(SimulationModel(sim))
function SimulationMethod(sim::AbstractMCMCSimulation) end

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