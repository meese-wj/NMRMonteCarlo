
import Base: enumerate
export energy, energy_change, DoF_energy, num_DoF, enumerate

"""
    abstract type AbstractHamiltonian end 

Supertype of all Hamiltonians. This `Type` defines the 
interface for each system defined by its Hamiltonian.

# Required Interface Methods

These `methods` _must_ be defined for every `<: AbstractHamiltonian`. By default they
`throw` `MethodError`s if not implemented.

- [`energy`](@ref)
- [`energy_change`](@ref)
- [`DoF_energy`](@ref)
- [`num_DoF`](@ref)
- [`enumerate`](@ref)

# Default Interface Methods

These `methods` are implemented by default for any `<: AbstractModel`. One is free 
to overload any of them for any peculiar `subtype`s.
"""
abstract type AbstractHamiltonian end

# Required Interface Methods
"""
    energy(::AbstractHamiltonian) -> MethodError

Return the total energy of a given [`AbstractHamiltonian`](@ref).
"""
energy(ham::AbstractHamiltonian) = throw(MethodError(energy, ham))
"""
    energy_change(::AbstractHamiltonian, args...) -> MethodError

Return the change in energy of a given [`AbstractHamiltonian`](@ref).
The `args` are variadic and can change for a given system or method.
"""
energy_change(ham::AbstractHamiltonian, args...) = throw(MethodError(energy_change, ham, args...))
"""
    DoF_energy(::AbstractHamiltonian, args...) -> MethodError

Return the energy for a given degree of freedom (DoF) of a given [`AbstractHamiltonian`](@ref).
The specific DoF should be specified by the `args`.
"""
DoF_energy(ham::AbstractHamiltonian, args...) = throw(MethodError(DoF_energy, ham, args...))
"""
    num_DoF(::AbstractHamiltonian) -> MethodError

Return the number of degrees of freedom of a given [`AbstractHamiltonian`](@ref).
"""
num_DoF(ham::AbstractHamiltonian) = throw(MethodError(num_DoF, ham))
"""
    enumerate(::AbstractHamiltonian, ::Type{<: HamiltonianIterationScheme}) -> MethodError

Define a set of enumeration iterators to traverse an [`AbstractHamiltonian`](@ref) with 
according to a specific [`HamiltonianIterationScheme`](@ref).

Should always return `(idx, value)`-`Tuple`s.
"""
enumerate(ham::AbstractHamiltonian, type::Type{<: HamiltonianIterationScheme}) = throw(MethodError(enumerate_DoFs, ham, type))
"""
    enumerate(::AbstractHamiltonian, ::Type{IterateByDoFType}) -> MethodError

Define a set of enumeration iterators to traverse an [`AbstractHamiltonian`](@ref) with 
according to [`IterateByDoFType`](@ref) scheme. This is _required_ for defining generic
sweeps through a system.

Should always return `(idx, value)`-`Tuple`s.
"""
enumerate(ham::AbstractHamiltonian, type::Type{IterateByDoFType}) = throw(MethodError(enumerate_DoFs, ham, type))

# Default Interface Methods

