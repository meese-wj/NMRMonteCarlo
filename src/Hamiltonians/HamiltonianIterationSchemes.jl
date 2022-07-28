
export traverse, enumerate

"""
abstract type AbstractEnumerationScheme end 

Interface for different [`AbstractHamiltonian`](@ref) iteration/enumeration `methods`.
This is necessary because, depending on the specific `method` context, the specific 
way one needs to traverse a system's degrees of freedom (DoF) may change.
"""
abstract type HamiltonianTraversalScheme end
"""
    ByDefault <: HamiltonianTraversalScheme

Easy iteration scheme. Simply traverse through the `<:` [`AbstractHamiltonian`](@ref)
as it is stored by default in whatever `iter`able container is used.
"""
struct ByDefault <: HamiltonianTraversalScheme end
"""
    BySite <: HamiltonianTraversalScheme

Traverse the degrees of freedom by outputting a `Tuple` of values at each define site
or within each unit cell.
"""
struct BySite <: HamiltonianTraversalScheme end
"""
    ByDoFType <: HamiltonianTraversalScheme

Traverse the degrees of freedom by their type. Useful for cases where sweeps are 
defined with one variable type while all others are _quenched_. 
"""
struct ByDoFType <: HamiltonianTraversalScheme end

"""
    TraversalScheme(::Type{T}) where T <: HamiltonianTraversalScheme = T()

Create an iterator trait for various `subtype`s of [`HamiltonianTraversalScheme`](@ref).
"""
TraversalScheme(::Type{T}) where {T <: HamiltonianTraversalScheme} = T()