
export energy, energy_change, DoF_energy, num_DoF

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

In addition to these required methods, one must define the `iterate` interface for a 
given `<:` [`AbstractHamiltonian`](@ref) using a [`HamiltonianIterator`](@ref).

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

# Default Interface Methods

