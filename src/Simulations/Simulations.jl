module Simulations

include("AbstractSimulations.jl")

# Automatically submodulizes this code
import ..MonteCarloMethods: AbstractModel
import ..Models: CleanNMRAshkinTellerModel

export CleanNMRATMParameters, CleanNMRATMSimulation

struct CleanNMRATMParameters{T <: AbstractFloat}
    Lx::Int
    Ly::Int
    Jex::T
    Kex::T
    Ntherm::Int
    Lτ::Int
    Nmeas::Int
    βvalue::T
end

struct CleanNMRATMSimulation{T <: AbstractFloat} <: AbstractMCMCSimulation
    parameters::CleanNMRATMParameters{T}
    model::CleanNMRAshkinTellerModel{T}

    CleanNMRATMSimulation(args...) = CleanNMRATMSimulation{Float64}(args...)
end

function CleanNMRATMSimulation(; Lx = 8, Ly = Lx, Jex = 1.0, Kex = 0.0, 
                                 Ntherm = 2^16, Lτ = 2^18, Nmeas = 2^16, βvalue = 2.269 )
    Jex, Kex, βvalue = promote(Jex, Kex, βvalue)
    params = CleanNMRATMParameters{typeof(Jex)}( Lx, Ly, Jex, Kex, Ntherm, Lτ, Nmeas, βvalue )
    model = CleanNMRAshkinTellerModel{typeof(Jex)}( Lx, Ly, Jex, Kex, Nmeas )
    return CleanNMRATMSimulation{typeof(Jex)}( params, model )
end

end # Simulations