
import StaticArrays: @SVector
import Base: getindex, setindex!, iterate, enumerate
export 
# Base overloads
       getindex, setindex!, eltype, length, iterate, enumerate,
# Ashkin-Teller functionality
       num_colors, spins, AT_sigma, AT_tau, Index, 
       AshkinTellerHamiltonian, AshkinTellerParameters

getTypename(::Type{T}) where T = isempty(T.parameters) ? T : T.name.wrapper

############################################################################
#             Abstract Ashkin-Teller Hamiltonian Interface                 #
############################################################################

abstract type AbstractAshkinTeller <: AbstractHamiltonian end
num_colors(type::Type{<: AbstractAshkinTeller}) = throw(MethodError(num_colors, type))
num_colors(::ATType) where ATType = num_colors(ATType)
spins(ham::AbstractAshkinTeller) = ham.spins

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

"""
iterate(::Type{IterateByMemory}, ham::TwoC_ATH) = length(ham) > zero(Int) ? (one(Int), spins(ham)[begin]) : nothing
function iterate(::Type{IterateByMemory}, ham::TwoC_ATH, state)
    idx = state[begin]
    next = idx == length(ham) ? nothing : (idx + one(Int), spins[idx + one(Int)])
    return next
end