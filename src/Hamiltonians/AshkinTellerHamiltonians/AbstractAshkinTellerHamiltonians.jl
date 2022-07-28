
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
    traverse(ByMemory, ::AbstractTwoColorAshkinTellerHamiltonian)

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) as it is
laid out in memory.
"""
function traverse(::Type{ByDefault}, ham::TwoC_ATH)
    iter = MVector{Tuple{Int, eltype(ham)}, num_DoF(ham)}(undef)
end
"""
    traverse(TraverseByDoFType, ::AbstractTwoColorAshkinTellerHamiltonian)

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) by the 
spin types separately.
"""
traverse(::Type{ByDoFType}, ham::TwoC_ATH) = length(ham) >= num_colors(ham) ? ( ham.color_update = AT_sigma; ( one(Int), ham[one(Int), color_update(ham)] )) : nothing

function _next_state!(::Type{TraverseByMemory}, ham::TwoC_ATH, state)
    idx = state[begin]
    next = idx == length(ham) ? nothing : (idx + one(Int), spins[idx + one(Int)])
    return next
end

function _next_state!(::Type{TraverseByDoFType}, ham::TwoC_ATH, state)
    site_index = state[begin]
    next_site = site_index + one(Int)
    end_of_iteration = false
    if site_index == num_sites(ham)
        switch_color_update!(ham)
        next_site = one(Int)
        end_of_iteration = color_update(ham) === AT_tau
    end
    return end_of_iteration ? nothing : (next_site, ham[next_site, color_update(ham)])
end

traverse(::Type{T}, ham::TwoC_ATH, state) where T <: HamiltonianTraversalScheme = _next_state!(T, ham, state)