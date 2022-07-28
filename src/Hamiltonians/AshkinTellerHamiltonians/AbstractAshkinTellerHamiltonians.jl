
import StaticArrays: @SVector, @MVector
import Base: getindex, setindex!
export 
# Base overloads
       getindex, setindex!, eltype, length,
# Ashkin-Teller functionality
       num_colors, spins, AT_sigma, AT_tau, Index

getTypename(::Type{T}) where T = isempty(T.parameters) ? T : T.name.wrapper

############################################################################
#             Abstract Ashkin-Teller Hamiltonian Interface                 #
############################################################################

abstract type AbstractAshkinTeller <: AbstractHamiltonian end
num_colors(type::Type{<: AbstractAshkinTeller}) = throw(MethodError(num_colors, type))
num_colors(::ATType) where ATType = num_colors(ATType)
spins(ham::AbstractAshkinTeller) = ham.spins
color_update(ham::AbstractAshkinTeller) = ham.color_update

Base.eltype(ham::AbstractAshkinTeller) = eltype(spins(ham))
Base.length( ham::AbstractAshkinTeller ) = length(spins(ham))
num_sites( ham::AbstractAshkinTeller ) = length(ham) ÷ num_colors(ham)
# Base.size( ham::AbstractAshkinTeller ) = ( num_sites(ham), num_sites(ham) )
num_DoF( ham::AbstractAshkinTeller ) = length(ham)

abstract type AshkinTellerColor end
const ATColorList = @SVector [ :AT_sigma, :AT_tau ]
for (idx, col) ∈ enumerate(ATColorList)
       @eval struct $col <: AshkinTellerColor end
       @eval @inline Index(::Type{$col}) = $idx
       @eval @inline Index(AT::Type{<: AbstractAshkinTeller}, site, ::Type{$col}) = num_colors(AT) * (site - one(Int)) + Index($col) 
       @eval @inline Index(::ATType, site, type::Type{$col}) where ATType = Index(ATType, site, type) 
       @eval getindex(ham::AbstractAshkinTeller, site, ::Type{$col}) = spins(ham)[Index(ham, site, $col)]
       @eval setindex!(ham::AbstractAshkinTeller, value, site, ::Type{$col}) = spins(ham)[Index(ham, site, $col)] = value
end

############################################################################
#        Abstract Two-Color Ashkin-Teller Hamiltonian Interface            #
############################################################################

abstract type AbstractTwoColorAshkinTellerHamiltonian <: AbstractAshkinTeller end
const TwoC_ATH = AbstractTwoColorAshkinTellerHamiltonian
@inline num_colors(::Type{<: TwoC_ATH}) = 2

"""
    iterate(::HamiltonianIterator{<: AbstractTwoColorAshkinTellerHamiltonian, IterateByDefault}, [state])

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) as it is laid out in memory.
"""
function iterate(iter::HamiltonianIterator{<: TwoC_ATH, IterateByDefault}, state = (one(Int),))
    ham = Hamiltonian(iter)
    spin_idx, = state
    return spin_idx <= length(ham) ? (spins(ham)[spin_idx], (spin_idx + one(Int),)) : nothing
end
"""
    iterate(::HamiltonianIterator{<: AbstractTwoColorAshkinTellerHamiltonian, IterateByDoFType}, [state])

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) by the spin types separately.
"""
function iterate(iter::HamiltonianIterator{<: TwoC_ATH, IterateByDoFType}, state = (one(Int), AT_sigma))
    ham = Hamiltonian(iter)
    site_idx, color = state
    ham.color_update = color
    next = site_idx + one(Int), color
    iteration_complete = false
    if site_idx == num_sites(ham)
        switch_color_update!(ham)
        next = one(Int), ifelse( color === AT_sigma, AT_tau, AT_sigma )
        iteration_complete = color === AT_tau
    end
    return iteration_complete ? nothing : ( ham[site_idx, color], next )
end
