
export iterate, enumerate

"""
    abstract type HamiltonianIterationScheme end 

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
struct IterateByDefault <: HamiltonianIterationScheme end
"""
    IterateBySite <: HamiltonianIterationScheme

Traverse the degrees of freedom by outputting a `Tuple` of values at each define site
or within each unit cell.
"""
struct IterateBySite <: HamiltonianIterationScheme end
"""
    IterateByDoFType <: HamiltonianIterationScheme

Traverse the degrees of freedom by their type. Useful for cases where sweeps are 
defined with one variable type while all others are _quenched_. 
"""
struct IterateByDoFType <: HamiltonianIterationScheme end

"""
    IterationScheme(::Type{T}) where T <: HamiltonianIterationScheme = T()

Create an iterator trait for various `subtype`s of [`HamiltonianIterationScheme`](@ref).
"""
IterationScheme(::Type{T}) where {T <: HamiltonianIterationScheme} = T()

"""
    HamiltonianIterator{H, S <: HamiltonianIterationScheme}

A wrapper iterator around a given Hamiltonian of type `H`. 

# Fields 

- `iter`: a pointer to the reference Hamiltonian

# Constructor

    HamiltonianIterator(ham::H, ::Type{S})

Creates a new `HamiltonianIterator` with respect to the Hamiltonian 
`ham` with a given `S <:`[`HamiltonianIterationScheme`](@ref).
"""
struct HamiltonianIterator{H, S <: HamiltonianIterationScheme} 
    iter::H
    HamiltonianIterator(ham::H, ::Type{S}) where {H, S} = new{H, S}(ham)
end
"""
    IterationScheme(::HamiltonianIterator{H, S}) → S 

Return the [`HamiltonianIterationScheme`](@ref) from a [`HamiltonianIterator`](@ref).
"""
IterationScheme(::HamiltonianIterator{H, S}) where {H, S} = S()
"""
    HamiltonianType(::HamiltonianIterator{H, S}) → H 

Return the Hamiltonian `Type` from a [`HamiltonianIterator`](@ref).
"""
HamiltonianType(::HamiltonianIterator{H, S}) where {H, S} = H
"""
    Hamiltonian(ham::HamiltonianIterator) → ham.iter

Return the specific Hamiltonian referenced by `ham.iter`.
"""
Hamiltonian(ham::HamiltonianIterator) = ham.iter
"""
    iterate(scheme::Type{<: HamiltonianIterationScheme}, ham, args...)

Wrapper for `iterate(HamiltonianIterator(ham, scheme), args...)`.
"""
iterate(scheme::Type{<: HamiltonianIterationScheme}, ham, args...) = iterate( HamiltonianIterator(ham, scheme), args...)
"""
    enumerate(scheme::Type{<: HamiltonianIterationScheme}, ham)

Wrapper for `enumerate(HamiltonianIterator(ham, scheme))`.
"""
enumerate(scheme::Type{<: HamiltonianIterationScheme}, ham) = enumerate(HamiltonianIterator(ham, scheme))
