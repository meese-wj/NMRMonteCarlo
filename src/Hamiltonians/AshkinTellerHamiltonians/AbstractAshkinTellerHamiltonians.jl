
import StaticArrays: @SVector, @MVector
import Base: getindex, setindex!, to_index
import ..Lattices: num_sites  # Automatically submodulizes this code
export 
# Base overloads
       getindex, setindex!, eltype, length, to_index,
# Ashkin-Teller functionality
       num_colors, spins, AT_sigma, AT_tau, ColorIndex

getTypename(::Type{T}) where T = isempty(T.parameters) ? T : T.name.wrapper

############################################################################
#             Abstract Ashkin-Teller Hamiltonian Interface                 #
############################################################################

abstract type AbstractAshkinTeller <: AbstractHamiltonian end

num_colors(::Type{T}) where {T <: AbstractAshkinTeller} = throw(MethodError(num_colors, T))
num_colors(::ATType) where ATType = num_colors(ATType)

spin_index(::Type{T}, site_idx, color_idx) where {T <: AbstractAshkinTeller} = num_colors(T) * (site_idx - one(Int)) + to_index(color_idx)
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

@enum AshkinTellerColor AT_sigma=1 AT_tau
@inline to_index(color::AshkinTellerColor) = Int(color)
@inline to_AshkinTellerColor(idx) = instances(AshkinTellerColor)[idx]
@inline getindex(ham::T, site, color) where {T <: AbstractAshkinTeller} = spins(ham)[spin_index(T, site, color)]
@inline setindex!(ham::T, value, site, color) where {T <: AbstractAshkinTeller} = spins(ham)[spin_index(T, site, color)] = value
# const ATColorList = @SVector [ :AT_sigma, :AT_tau ]
# for (idx, col) ∈ enumerate(ATColorList)
#     # ******************************************************************************
#     # Conversions between indices and singletons
#     # ******************************************************************************
#     # Define the singleton type
#     @eval struct $col <: AshkinTellerColor end
#     # Convert the singleton type to an index
#     @eval @inline ColorIndex(::Type{$col}) = $idx
#     # Dispatch from the Type to values of that singleton type (for stability)
#     @eval @inline ColorIndex(::$col) = ColorIndex($col)
#     # Invert the ColorIndex function from integers to singletons
#     @eval @inline ColorIndex(::Type{AshkinTellerColor}, ::Type{Val{$idx}}) = $col()
#     # Create an equivalent inversion for the AbstractTwoColorAshkinTellerHamiltonian
#     @eval @inline ColorIndex(::Type{TwoC_ATH}, ::Type{Val{$idx}}) = ColorIndex(AshkinTellerColor, Val{$idx})
#     # ******************************************************************************

#     # ******************************************************************************
#     # Array access functionality
#     # ******************************************************************************
#     # Map a singleton type to a specific index and access the spin from it
#     @eval @inline spin_index(::ATType, site, type::Type{$col}) where ATType = spin_index(ATType, site, ColorIndex(type)) 
#     # Overload getindex for an AbstractAshkinTeller Hamiltonian with a singleton type index
#     @eval getindex(ham::AbstractAshkinTeller, site, ::Type{$col}) = spins(ham)[spin_index(ham, site, $col)]
#     # Overload getindex for an AbstractAshkinTeller Hamiltonian with a singleton index (for type stability)
#     @eval getindex(ham::AbstractAshkinTeller, site, ::$col) = ham[site, $col]
#     # Overload setindex! for an AbstractAshkinTeller Hamiltonian with a singleton type index
#     @eval setindex!(ham::AbstractAshkinTeller, value, site, ::Type{$col}) = spins(ham)[spin_index(ham, site, $col)] = value
#     # Overload setindex! for an AbstractAshkinTeller Hamiltonian with a singleton index (for type stability)
#     @eval setindex!(ham::AbstractAshkinTeller, value, site, ::$col) = ham[site, $col] = value
#     # ******************************************************************************
# end
# # Have a cutoff for the ATColorList for AbstractTwoColorAshkinTellerHamiltonian
# ColorIndex(::Type{AbstractTwoColorAshkinTellerHamiltonian}, ::Type{Val{num_colors(TwoC_ATH) + one(Int)}}) = ATDefaultColor()
# # Overload the ColorIndex inverse to take any idx index to a specific singleton
# ColorIndex(::Type{T}, idx::Int) where {T <: AshkinTellerColor} = ColorIndex(T, Val{idx})
# # Create a default value for the Ashkin Teller colors
# ATDefaultColor() = ColorIndex(AshkinTellerColor, one(Int))

############################################################################
#        Abstract Two-Color Ashkin-Teller Hamiltonian Interface            #
############################################################################

abstract type AbstractTwoColorAshkinTellerHamiltonian <: AbstractAshkinTeller end
const TwoC_ATH = AbstractTwoColorAshkinTellerHamiltonian
@inline num_colors(::Type{<: TwoC_ATH}) = 2
@inline site_Baxter( ham::AbstractTwoColorAshkinTellerHamiltonian, site ) = ham[site, AT_sigma] * ham[site, AT_tau]
@inline ATDefaultColor(::Type{<: TwoC_ATH}) = AT_sigma

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
    iterate(::HamiltonianIterator{<: AbstractTwoColorAshkinTellerHamiltonian, IterateBySite}, [state])

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) and return an SVector of sigma, tau, and
Baxter values.
"""
function iterate(iter::HamiltonianIterator{<: TwoC_ATH, IterateBySite}, state = (one(Int),))
    ham = Hamiltonian(iter)
    site, = state
    next = (site + one(Int), )
    return site <= num_sites(ham) ? ((site, (@SVector [ ham[site, AT_sigma], ham[site, AT_tau], site_Baxter(ham, site) ])), next) : nothing
end

"""
    iterate(::HamiltonianIterator{<: AbstractTwoColorAshkinTellerHamiltonian, IterateByDoFType}, [state])

Traverse a `<:`[`AbstractTwoColorAshkinTellerHamiltonian`](@ref) by the spin types separately.
"""
function iterate(iter::HamiltonianIterator{T, IterateByDoFType}, state = (one(Int), one(Int), to_index(ATDefaultColor(T)))) where {T <: TwoC_ATH}
    ham = Hamiltonian(iter)
    site_idx = state[1]
    iterations = state[2]
    # @show state
    # color = ifelse( state[end] <= num_colors(ham), ColorIndex(TwoC_ATH, Val{state[end]}), ColorIndex(TwoC_ATH, Val{state[end]}) )
    # color = state[end]
    ham.color_update = state[end] <= num_colors(ham) ? to_AshkinTellerColor(state[end]) : ATDefaultColor(T)
    next_iter = iterations + one(Int)
    next_site = one(Int) + iterations % num_sites(ham)
    next_color_int = one(Int) + iterations ÷ num_sites(ham)
    # if site_idx == num_sites(ham)
    #     next_site = one(Int)
    #     @show next_color_int = ifelse( color === AT_sigma, ColorIndex(AT_tau), ColorIndex(AT_sigma) )
    # end
    return iterations <= length(ham) ? 
               ( (site_idx, ham[site_idx, color_update(ham)]), (next_site, next_iter, next_color_int) ) : 
               nothing
end
