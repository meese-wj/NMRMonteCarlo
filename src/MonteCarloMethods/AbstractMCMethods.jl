using Parameters2JSON

import ..Hamiltonians: Hamiltonian   # Automatically submodulizes this code
export sweeps_per_export, Lattice, Hamiltonian, Observables, update_observables!

abstract type AbstractMonteCarloParameters end
sweeps_per_export(params::AbstractMonteCarloParameters) = params.measure_sweeps <= params.total_measurements ? 1 : params.measure_sweeps ÷ params.total_measurements

# # TODO: Change this to an abstract type
# struct Model{L <: AbstractLattice, H <: AbstractHamiltonian}
#     latt::L
#     ham::H
# end
"""
    abstract type AbstractModel end

Supertype for all model `Type`s to be used. We take the approach that any instance
of an `<: AbstractModel` is uniquely defined by _three_ fields:

1. A [`lattice`](@ref) (or some other _geometry_)
1. A [`Hamiltonian`](@ref) (demarking the _system_ to be simulated)
1. A set of [`observables`](@ref) (defined by the _system_ and the _geometry_ it lives in)

# Required Interface Methods

These `methods` _must_ be defined for every `<: AbstractModel`. By default they
`throw` `MethodError`s if not implemented.

- [`update_observables!`](@ref)

# Default Interface Methods

These `methods` are implemented by default for any `<: AbstractModel`. One is free 
to overload any of them for any peculiar `subtype`s.

- [`Lattice`](@ref)
- [`Hamiltonian`](@ref)
- [`Observables`](@ref)
"""
abstract type AbstractModel end

# Required Interface Methods
function update_observables!(::AbstractModel) end

# Default Interface Methods
"""
    Lattice(::AbstractModel) -> MethodError

Returns the lattice or geometry of an [`AbstractModel`](@ref) `subtype`.
"""
# Lattice(model::AbstractModel) = throw(MethodError(lattice, model))
Lattice(model::AbstractModel) = model.lattice

"""
    Hamiltonian(::AbstractModel) -> MethodError

Returns the system in an [`AbstractModel`](@ref) `subtype`.
"""
# Hamiltonian(model::AbstractModel) = throw(MethodError(Hamiltonian, model))
Hamiltonian(model::AbstractModel) = model.hamiltonian

"""
    Observables(::AbstractModel) -> MethodError

Returns the set of defined observables for an [`AbstractModel`](@ref) `subtype`.
"""
Observables(model::AbstractModel) = model.observables