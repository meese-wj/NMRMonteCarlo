
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
num_colors(::Type{T}) where {T <: AbstractAshkinTeller} = throw(MethodError(num_colors, T))
num_colors(::ATType) where ATType = num_colors(ATType)
spin_index(::Type{T}, site_idx, color_idx) where {T <: AbstractAshkinTeller} = num_colors(T) * (site_idx - one(Int)) + color_idx
spin_index(::ATType, site_idx, color_idx) where ATType = num_colors(ATType) * (site_idx - one(Int)) + color_idx
site_index(::Type{T}, spin_idx) where {T <: AbstractAshkinTeller} = one(Int) + spin_idx ÷ num_colors(T)
site_index(::ATType, spin_idx) where ATType = site_index(ATType, spin_idx)
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
    @eval @inline spin_index(::ATType, site, type::Type{$col}) where ATType = spin_index(ATType, site, Index(type)) 
    @eval getindex(ham::AbstractAshkinTeller, site, ::Type{$col}) = spins(ham)[spin_index(ham, site, $col)]
    @eval setindex!(ham::AbstractAshkinTeller, value, site, ::Type{$col}) = spins(ham)[spin_index(ham, site, $col)] = value
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
    dof_location_value = site_index(ham, spin_idx), spins(ham)[spin_idx]
    return spin_idx <= length(ham) ? (dof_location_value, (spin_idx + one(Int),)) : nothing
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
        next = one(Int), ifelse( color === AT_sigma, AT_tau, AT_sigma )
        iteration_complete = color === AT_tau
    end
    dof_location_value = site_idx, ham[site_idx, color]
    return iteration_complete ? nothing : ( dof_location_value, next )
end
