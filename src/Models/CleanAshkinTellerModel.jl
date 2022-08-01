
import ..Lattices: CubicLattice2D
import ..Hamiltonians: AshkinTellerHamiltonian, AshkinTellerParameters, 
                       energy, num_sites, spins, AT_sigma, AT_tau, num_colors,
                       site_Baxter
import ..MonteCarloMethods: AbstractModel, Lattice, Hamiltonian, Observables
import StaticArrays: @SVector
import Base: getindex, setindex!
using MonteCarloMeasurementUncertainty

export CleanAshkinTellerModel, CleanAshkinTellerModelParameters, update_observables!

abstract type CATMObservable end
const CATM_observables = ["Energy", "Energy2", "Sigma", "Tau", "Baxter"]
for (idx, obs) ∈ enumerate(CATM_observables)
    obs_symb = Symbol(obs)
    @eval struct $obs_symb <: CATMObservable end
    @eval ObservableIndex(::Type{$obs_symb}) = $idx
    @eval getindex(vec::AbstractArray, ::Type{$obs_symb}) = vec[ObservableIndex($obs_symb)]
    @eval getindex(vec::AbstractArray, ::$obs_symb) = vec[$obs_symb]
    @eval ObservableString(::Type{$obs_symb}) = String(Symbol($obs_symb))
    @eval ObservableString(::$obs_symb) = ObservableString($obs_symb)
    @eval ObservableType(::Type{CATMObservable}, ::Type{Val{$idx}}) = $obs_symb
end
ObservableString(s::Symbol) = ObservableString(eval(s))
ObservableType(type::Type{CATMObservable}, idx) = ObservableType(type, Val{idx})

struct CleanAshkinTellerModelParameters{T <: AbstractFloat}
    Lx::Int
    Ly::Int
    Jex::T
    Kex::T
    num_measurements::Int
    
    CleanAshkinTellerModelParameters{T}(Lx, Ly, Jex, Kex, n) where T = new{T}(Int(Lx), Int(Ly), T(Jex), T(Kex), Int(n))
    CleanAshkinTellerModelParameters(Lx, Ly, Jex, Kex, n) = CleanAshkinTellerModelParameters{Float64}(Lx, Ly, Jex, Kex, n)
end

struct CleanAshkinTellerModel{T <: AbstractFloat} <: AbstractModel
    lattice::CubicLattice2D
    hamiltonian::AshkinTellerHamiltonian{T}
    observables::Vector{TimeSeries}

    CleanAshkinTellerModel(args...) = CleanAshkinTellerModel{Float64}(args...)
    function CleanAshkinTellerModel{T}( Lx, Ly, Jex = 1, Kex = 0, num_meas = 0 ) where T
        latt = CubicLattice2D(Int(Lx), Int(Ly))
        atparams = AshkinTellerParameters( T(Jex), T(Kex) )
        ham = AshkinTellerHamiltonian(latt, atparams)
        obs_list = TimeSeries{T}[]
        for obs ∈ CATM_observables
            push!(obs_list, TimeSeries{T}( obs, num_meas ))
        end
        return new{T}( latt, ham, obs_list )
    end
    CleanAshkinTellerModel{T}(params::CleanAshkinTellerModelParameters{T}) where T = CleanAshkinTellerModel( params.Lx, params.Ly, params.Jex, params.Kex, params.num_measurements )
end

function update_observables!(model::CleanAshkinTellerModel)
    for idx ∈ eachindex(Observables(model))
        update_observable!(model, ObservableType(CATMObservable, idx))
    end
end

# Individual observable definitions
function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Energy})
    @show en = energy(Hamiltonian(model))
    push!( obs, energy(Hamiltonian(model)) / num_sites(Hamiltonian(model)) )
end

function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Energy2})
    push!( obs, (energy(Hamiltonian(model)) / num_sites(Hamiltonian(model)))^2 )
end

function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Sigma})
    ham = Hamiltonian(model)
    push!( obs, sum( @view spins(ham)[ColorIndex(AT_sigma):num_colors(ham):end] ) )
end

function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Tau})
    ham = Hamiltonian(model)
    push!( obs, mean( @view spins(ham)[ColorIndex(Tau):num_colors(ham):end] ) )
end

function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Baxter})
    ham = Hamiltonian(model)
    push!( obs, mean( site_idx -> site_Baxter(ham, site_idx), (1:num_sites(ham)) ) )
end

# Dispatch on observable types
for obs_str ∈ CATM_observables
    obs_symb = Symbol(obs_str)
    @eval update_observable!( model::CleanAshkinTellerModel, ::Type{$obs_symb} ) = update_observable!( Observables(model)[$obs_symb], model, $obs_symb )
end

