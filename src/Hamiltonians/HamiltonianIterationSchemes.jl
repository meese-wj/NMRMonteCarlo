
import Base: iterate, enumerate
export iterate, enumerate

"""
abstract type AbstractEnumerationScheme end 

Interface for different [`AbstractHamiltonian`](@ref) iteration/enumeration `methods`.
This is necessary because, depending on the specific `method` context, the specific 
way one needs to traverse a system's degrees of freedom (DoF) may change.
"""
abstract type HamiltonianIterationScheme end
"""
IterateByDefault <: HamiltonianIterationScheme

Easy iteration scheme. Simply traverse through the `<:` [`AbstractHamiltonian`](@ref)
as it is stored by default in whatever `iter`able container is used.
"""
struct IterateByMemory <: HamiltonianIterationScheme end
"""
IterateBySite <: HamiltonianIterationScheme

Iterate the degrees of freedom by outputting a `Tuple` of values at each define site
or within each unit cell.
"""
struct IterateBySite <: HamiltonianIterationScheme end
"""
IterateByDoFType <: HamiltonianIterationScheme

Iterate the degrees of freedom by their type. Useful for cases where sweeps are 
defined with one variable type while all others are _quenched_. 
"""
struct IterateByDoFType <: HamiltonianIterationScheme end

"""
    IterationScheme(::Type{T}) where T <: HamiltonianIterationScheme = T()

Create an iterator trait for various `subtype`s of [`HamiltonianIterationScheme`](@ref).
"""
IterationScheme(::Type{T}) where {T <: HamiltonianIterationScheme} = T()

iterate(iter::Tuple{<: HamiltonianIterationScheme, <: Any}, args...) = iterate(iter..., args...)
iterate(a::Any, b::Type{<: HamiltonianIterationScheme}, args...) = iterate((b, a), args...)
enumerate(a::Any, b::Type{<: HamiltonianIterationScheme}) = enumerate((b, a))