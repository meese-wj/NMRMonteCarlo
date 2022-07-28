
export traverse, enumerate

"""
abstract type AbstractEnumerationScheme end 

Interface for different [`AbstractHamiltonian`](@ref) iteration/enumeration `methods`.
This is necessary because, depending on the specific `method` context, the specific 
way one needs to traverse a system's degrees of freedom (DoF) may change.
"""
abstract type HamiltonianTraversalScheme end
"""
TraverseByDefault <: HamiltonianTraversalScheme

Easy iteration scheme. Simply traverse through the `<:` [`AbstractHamiltonian`](@ref)
as it is stored by default in whatever `iter`able container is used.
"""
struct TraverseByMemory <: HamiltonianTraversalScheme end
"""
TraverseBySite <: HamiltonianTraversalScheme

Traverse the degrees of freedom by outputting a `Tuple` of values at each define site
or within each unit cell.
"""
struct TraverseBySite <: HamiltonianTraversalScheme end
"""
TraverseByDoFType <: HamiltonianTraversalScheme

Traverse the degrees of freedom by their type. Useful for cases where sweeps are 
defined with one variable type while all others are _quenched_. 
"""
struct TraverseByDoFType <: HamiltonianTraversalScheme end

"""
    TraversalScheme(::Type{T}) where T <: HamiltonianTraversalScheme = T()

Create an iterator trait for various `subtype`s of [`HamiltonianTraversalScheme`](@ref).
"""
TraversalScheme(::Type{T}) where {T <: HamiltonianTraversalScheme} = T()

traverse(iter::Tuple{<: HamiltonianTraversalScheme, <: Any}, args...) = traverse(iter..., args...)
traverse(a::Any, b::Type{<: HamiltonianTraversalScheme}, args...) = traverse((b, a), args...)
enumerate(a::Any, b::Type{<: HamiltonianTraversalScheme}) = enumerate((b, a))