
import StaticArrays: @SVector
import Base: getindex, setindex!
export 
# Base overloads
       getindex, setindex!, eltype, length,
# Ashkin-Teller functionality
       num_colors, spins, AT_sigma, AT_tau, Index, 
       AshkinTellerHamiltonian, AshkinTellerParameters

getTypename(::Type{T}) where T = isempty(T.parameters) ? T : T.name.wrapper

abstract type AbstractAshkinTeller <: AbstractHamiltonian end
num_colors(type::Type{<: AbstractAshkinTeller}) = throw(MethodError(num_colors, type))
num_colors(::ATType) where ATType = num_colors(ATType)
spins(ham::AbstractAshkinTeller) = ham.spins

Base.eltype(ham::AbstractAshkinTeller) = eltype(spins(ham))
Base.length( ham::AbstractAshkinTeller ) = length(spins(ham))
num_sites( ham::AbstractAshkinTeller ) = length(ham) ÷ num_colors(ham)
# Base.size( ham::AbstractAshkinTeller ) = ( num_sites(ham), num_sites(ham) )
num_DoF( ham::AbstractAshkinTeller ) = length(ham)

abstract type AbstractTwoColorAshkinTellerHamiltonian <: AbstractAshkinTeller end
@inline num_colors(::Type{<: AbstractTwoColorAshkinTellerHamiltonian}) = 2

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