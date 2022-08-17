module Simulations

include("AbstractSimulations.jl")

# Automatically submodulizes this code
import ..MonteCarloMethods: AbstractModel, AbstractMonteCarloParameters, thermalization_sweeps, sampling_sweeps, total_measurements, metropolis_sweep!
import ..Models: CleanNMRAshkinTellerModel

export CleanNMRATMParameters, CleanNMRATMSimulation, thermalization_sweeps, sampling_sweeps, total_measurements


struct CleanNMRATMParameters{T <: AbstractFloat} <: AbstractMonteCarloParameters
    Lx::Int
    Ly::Int
    Jex::T
    Kex::T
    Ntherm::Int
    Lτ::Int
    Nmeas::Int
    βvalue::T
end
@inline thermalization_sweeps(params::CleanNMRATMParameters) = params.Ntherm
@inline sampling_sweeps(params::CleanNMRATMParameters) = params.Lτ
@inline total_measurements(params::CleanNMRATMParameters) = params.Nmeas

struct CleanNMRATMSimulation{T <: AbstractFloat} <: AbstractMCMCSimulation
    parameters::CleanNMRATMParameters{T}
    model::CleanNMRAshkinTellerModel{T}

    CleanNMRATMSimulation{T}(p::CleanNMRATMParameters{T}, m::CleanNMRAshkinTellerModel{T}) where {T} = new{T}(p, m)
    CleanNMRATMSimulation(args...) = CleanNMRATMSimulation{Float64}(args...)
end

function CleanNMRATMSimulation(; Lx = 8, Ly = Lx, Jex = 1.0, Kex = 0.0, 
                                 Ntherm = 2^16, Lτ = 2^18, Nmeas = 2^16, βvalue = 1/2.269 )
    Jex, Kex, βvalue = promote(Jex, Kex, βvalue)
    params = CleanNMRATMParameters{typeof(Jex)}( Lx, Ly, Jex, Kex, Ntherm, Lτ, Nmeas, βvalue )
    model = CleanNMRAshkinTellerModel{typeof(Jex)}( Lx, Ly, Jex, Kex, Nmeas )
    return CleanNMRATMSimulation{typeof(Jex)}( params, model )
end

SimulationMethod(sim::CleanNMRATMSimulation) = metropolis_sweep!

end # Simulations