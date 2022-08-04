
import ..Lattices: CubicLattice2D
import ..Hamiltonians: AshkinTellerHamiltonian, AshkinTellerParameters, 
                       energy, num_sites, spins, AT_sigma, AT_tau, num_colors,
                       site_Baxter, sigma_values, tau_values
import ..MonteCarloMethods: AbstractModel, Lattice, Hamiltonian, Observables, update_observables!
import StaticArrays: @SVector
import Statistics: mean
import Base: getindex, setindex!
using MonteCarloMeasurementUncertainty

export CleanAshkinTellerModel, CleanAshkinTellerModelParameters, update_observables!

abstract type CATMObservable end
const CATM_observables = @SVector ["Energy", "Energy2", "Sigma", "Tau", "Baxter"]
CATM_temp_obs_types = []
for (idx, obs) ∈ enumerate(CATM_observables)
    obs_symb = Symbol(obs)
    @eval struct $obs_symb <: CATMObservable end
    @eval push!(CATM_temp_obs_types, $obs_symb)
    @eval ObservableIndex(::Type{$obs_symb}) = $idx
    @eval getindex(vec::AbstractArray, ::Type{$obs_symb}) = vec[ObservableIndex($obs_symb)]
    @eval getindex(vec::AbstractArray, ::$obs_symb) = vec[$obs_symb]
    @eval ObservableString(::Type{$obs_symb}) = String(Symbol($obs_symb))
    @eval ObservableString(::$obs_symb) = ObservableString($obs_symb)
    @eval ObservableType(::Type{CATMObservable}, ::Type{Val{$idx}}) = $obs_symb()
end
const CATM_observable_types = @SVector [ type for type ∈ CATM_temp_obs_types ]
ObservableString(s::Symbol) = ObservableString(eval(s))

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
    observables::Vector{TimeSeries{T}}

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

"""
    @generated update_observables!(::CleanAshkinTellerModel)

Function to update all of the (possibly-many) observable 
types for the CleanAshkinTellerModel.

!!! note
    The naive way of doing this is type-unstable as the iterable
    observable `Type`s are different, and so the iterator itself is 
    type-unstable. So a generated function is needed to to manually
    unroll the for-loop over the observables. The naive function 
    defintion looks something like the following:

    ```julia
    function update_observables!(model::CleanAshkinTellerModel)
        for type ∈ CATM_observable_types
            update_observable!(model, type)
        end
        return Observables(model)
    end
    ```
"""
@generated function update_observables!(model::CleanAshkinTellerModel)
    expr = :()
    for type ∈ CATM_observable_types
        expr = :( $expr; update_observable!(model, $type) )
    end
    expr = :( $expr; return model )
    return expr
end

# Individual observable definitions
@inline function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Energy})
    push!( obs, energy(Hamiltonian(model), Lattice(model)) / num_sites(Hamiltonian(model)) )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Energy2})
    push!( obs, (energy(Hamiltonian(model), Lattice(model)) / num_sites(Hamiltonian(model)))^2 )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Sigma})
    ham = Hamiltonian(model)
    push!( obs, mean( sigma_values(Hamiltonian(model)) ) )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Tau})
    ham = Hamiltonian(model)
    push!( obs, mean( tau_values(Hamiltonian(model)) ) )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::CleanAshkinTellerModel, ::Type{Baxter})
    ham = Hamiltonian(model)
    push!( obs, mean( site_idx -> site_Baxter(ham, site_idx), (1:num_sites(ham)) ) )
end

# Dispatch on observable types
for obs_str ∈ CATM_observables
    obs_symb = Symbol(obs_str)
    @eval @inline update_observable!( model::CleanAshkinTellerModel, ::Type{$obs_symb} ) = update_observable!( Observables(model)[$obs_symb], model, $obs_symb )
end

