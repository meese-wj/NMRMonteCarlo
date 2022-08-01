
import StaticArrays: @SVector, @MVector
import Base: getindex, setindex!
import ..Lattices: num_sites  # Automatically submodulizes this code
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
spin_index(::ATType, site_idx, color_idx) where ATType = spin_index(ATType, site_idx, color_idx)

site_index(::Type{T}, spin_idx) where {T <: AbstractAshkinTeller} = one(Int) + (spin_idx - one(Int)) ÷ num_colors(T)
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
    @eval color_update(::Type{$col}) = $col()
    @eval color_update(::$col) = $col()
    @eval @inline Index(::Type{$col}) = $idx
    @eval @inline Index(::$col) = Index($col)
    @eval @inline Index(::Type{AshkinTellerColor}, ::Type{Val{$idx}}) = $col()
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
@inline site_Baxter( ham::AbstractTwoColorAshkinTellerHamiltonian, site ) = ham[site, AT_sigma] * ham[site, AT_tau]

"""
    iterate(::HamiltonianIterator{<: AbstractTwoColorAshkinTellerHamiltonian, IterateByDefault}, [state])

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) as it is laid out in memory.
"""
function iterate(iter::HamiltonianIterator{<: TwoC_ATH, IterateByDefault}, state = (one(Int),))
    ham = Hamiltonian(iter)
    spin_idx, = state
    next = (spin_idx + one(Int),)
    return spin_idx <= length(ham) ? ((site_index(ham, spin_idx), spins(ham)[spin_idx]), next) : nothing
end

"""
    iterate(::HamiltonianIterator{<: AbstractTwoColorAshkinTellerHamiltonian, IterateByDoFType}, [state])

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) by the spin types separately.
"""
function iterate(iter::HamiltonianIterator{<: TwoC_ATH, IterateByDoFType}, state = (one(Int), one(Int), one(Int)))
    ham = Hamiltonian(iter)
    site_idx = state[1]
    iterations = state[2]
    color = Index(AshkinTellerColor, Val{state[end]})
    ham.color_update = color
    next_site = site_idx + one(Int)
    next_iter = iterations + one(Int)
    next_color_int = Index(color)
    if site_idx == num_sites(ham)
        next_site = one(Int)
        next_color_int = ifelse( color === AT_sigma, Index(AT_tau), Index(AT_sigma) )
    end
    return iterations <= length(ham) ? ( (site_idx, ham[site_idx, color]), (next_site, next_iter, next_color_int) ) : nothing
end
