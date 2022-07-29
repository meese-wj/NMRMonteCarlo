
export energy, energy_change, DoF_energy, DoF_energy_change, num_DoF

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
function energy(::AbstractHamiltonian) end
"""
    energy_change(::AbstractHamiltonian, args...) -> MethodError

Return the change in energy of a given [`AbstractHamiltonian`](@ref).
The `args` are variadic and can change for a given system or method.
"""
function energy_change(::AbstractHamiltonian, args...) end
"""
    DoF_energy(::AbstractHamiltonian, args...) -> MethodError

Return the energy for a given degree of freedom (DoF) of a given [`AbstractHamiltonian`](@ref).
The specific DoF should be specified by the `args`.
"""
function DoF_energy(::AbstractHamiltonian, args...) end
"""
    num_DoF(::AbstractHamiltonian) -> MethodError

Return the number of degrees of freedom of a given [`AbstractHamiltonian`](@ref).
"""
function num_DoF(::AbstractHamiltonian) end

# Default Interface Methods

