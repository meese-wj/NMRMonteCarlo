
import ..MonteCarloMethods: AbstractModel, Lattice, Hamiltonian, Observables, update_observables!
import StaticArrays: @SVector
import Statistics: mean
import Base: getindex, setindex!
using MonteCarloMeasurementUncertainty

export update_observables!

"""
    abstract type AbstractAshkinTellerModel end

Super type for all [`AbstractModel`](@ref)s that implement any 
[`AbstractAshkinTellerHamiltonian`](@ref) as its [`Hamiltonian`](@ref)
type.
"""
abstract type AbstractAshkinTellerModel <: AbstractModel end
# TODO: Keep this in mind?
# abstract type AbstractTwoColorAshkinTellerModel <: AbstractAshkinTellerModel end
"""
    TimeSeriesObservables(::AbstractAshkinTellerModel)

Return the collection of `TimeSeries` `MonteCarloMeasurement`s.
"""
function TimeSeriesObservables(::AbstractAshkinTellerModel) end
"""
    AccumulatedSeriesObservables(::AbstractAshkinTellerModel)

Return the collection of `AccumulatedSeries` `MonteCarloMeasurement`s.
"""
function AccumulatedSeriesObservables(::AbstractAshkinTellerModel) end

abstract type BaseATMObservable end
const BaseATMTimeSeriesObs = @SVector ["Energy", "Energy2", "Sigma", "Tau", "Baxter"]
BaseATM_temp_obs_types = []
for (idx, obs) ∈ enumerate(BaseATMTimeSeriesObs)
    obs_symb = Symbol(obs)
    @eval struct $obs_symb <: BaseATMObservable end
    @eval push!(BaseATM_temp_obs_types, $obs_symb)
    @eval ObservableIndex(::Type{$obs_symb}) = $idx
    @eval getindex(vec::AbstractArray, ::Type{$obs_symb}) = vec[ObservableIndex($obs_symb)]
    @eval getindex(vec::AbstractArray, ::$obs_symb) = vec[$obs_symb]
    @eval ObservableString(::Type{$obs_symb}) = String(Symbol($obs_symb))
    @eval ObservableString(::$obs_symb) = ObservableString($obs_symb)
    @eval ObservableType(::Type{BaseATMObservable}, ::Type{Val{$idx}}) = $obs_symb()
end
const BaseATM_observable_types = @SVector [ type for type ∈ BaseATM_temp_obs_types ]
ObservableString(s::Symbol) = ObservableString(eval(s))

"""
    @generated update_TimeSeriesObservables!(::AbstractAshkinTellerModel)

Function to update all of the (possibly-many) observable 
types for the AbstractAshkinTellerModel.

!!! note
    The naive way of doing this is type-unstable as the iterable
    observable `Type`s are different, and so the iterator itself is 
    type-unstable. So a generated function is needed to to manually
    unroll the for-loop over the observables. The naive function 
    defintion looks something like the following:

    ```julia
    function update_TimeSeriesObservables!(model::AbstractAshkinTellerModel)
        for type ∈ BaseATM_observable_types
            update_observable!(model, type)
        end
        return Observables(model)
    end
    ```
"""
@generated function update_TimeSeriesObservables!(model::AbstractAshkinTellerModel)
    expr = :()
    for type ∈ BaseATM_observable_types
        expr = :( $expr; update_observable!(model, $type) )
    end
    expr = :( $expr; return model )
    return expr
end

# Individual observable definitions
@inline function update_observable!(obs::MonteCarloMeasurement, model::AbstractAshkinTellerModel, ::Type{Energy})
    push!( obs, energy(Hamiltonian(model), Lattice(model)) / num_sites(Hamiltonian(model)) )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::AbstractAshkinTellerModel, ::Type{Energy2})
    push!( obs, (energy(Hamiltonian(model), Lattice(model)) / num_sites(Hamiltonian(model)))^2 )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::AbstractAshkinTellerModel, ::Type{Sigma})
    ham = Hamiltonian(model)
    push!( obs, mean( sigma_values(Hamiltonian(model)) ) )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::AbstractAshkinTellerModel, ::Type{Tau})
    ham = Hamiltonian(model)
    push!( obs, mean( tau_values(Hamiltonian(model)) ) )
end

@inline function update_observable!(obs::MonteCarloMeasurement, model::AbstractAshkinTellerModel, ::Type{Baxter})
    ham = Hamiltonian(model)
    push!( obs, mean( site_idx -> site_Baxter(ham, site_idx), (1:num_sites(ham)) ) )
end

# Dispatch on observable types
for obs_str ∈ BaseATMTimeSeriesObs
    obs_symb = Symbol(obs_str)
    @eval @inline update_observable!( model::AbstractAshkinTellerModel, ::Type{$obs_symb} ) = update_observable!( TimeSeriesObservables(model)[$obs_symb], model, $obs_symb )
end
