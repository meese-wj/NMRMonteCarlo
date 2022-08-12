
import StaticArrays: @SVector, @MVector
import Base: getindex, setindex!, to_index
import ..Lattices: num_sites  # Automatically submodulizes this code
export 
# Base overloads
       getindex, setindex!, eltype, length, to_index,
# Ashkin-Teller functionality
       num_colors, spins, AT_sigma, AT_tau, ColorIndex

############################################################################
#             Abstract Ashkin-Teller Hamiltonian Interface                 #
############################################################################

abstract type AbstractAshkinTeller <: AbstractHamiltonian end
abstract type AbstractAshkinTellerParameters <: AbstractHamiltonianParameters end

@inline num_colors(::Type{T}) where {T <: AbstractAshkinTeller} = throw(MethodError(num_colors, T))
@inline num_colors(::ATType) where ATType = num_colors(ATType)

spin_index(::Type{T}, site_idx, color_idx) where {T <: AbstractAshkinTeller} = num_colors(T) * (site_idx - one(Int)) + to_index(color_idx)
spin_index(::ATType, site_idx, color_idx) where ATType = spin_index(ATType, site_idx, color_idx)

@inline site_index(::Type{T}, spin_idx) where {T <: AbstractAshkinTeller} = one(Int) + (spin_idx - one(Int)) ÷ num_colors(T)
@inline site_index(::ATType, spin_idx) where ATType = site_index(ATType, spin_idx)

@inline spins(ham::AbstractAshkinTeller) = ham.spins
@inline color_update(ham::AbstractAshkinTeller) = ham.color_update

@inline Base.eltype(ham::AbstractAshkinTeller) = eltype(spins(ham))
@inline Base.length( ham::AbstractAshkinTeller ) = length(spins(ham))
@inline num_sites( ham::AbstractAshkinTeller ) = length(ham) ÷ num_colors(ham)
# Base.size( ham::AbstractAshkinTeller ) = ( num_sites(ham), num_sites(ham) )
@inline num_DoF( ham::AbstractAshkinTeller ) = length(ham)

@enum AshkinTellerColor AT_sigma=1 AT_tau
@inline to_index(color::AshkinTellerColor) = Int(color)
@inline to_AshkinTellerColor(idx) = instances(AshkinTellerColor)[idx]
@inline getindex(ham::T, site, color) where {T <: AbstractAshkinTeller} = spins(ham)[spin_index(T, site, color)]
@inline setindex!(ham::T, value, site, color) where {T <: AbstractAshkinTeller} = spins(ham)[spin_index(T, site, color)] = value

############################################################################
#        Abstract Two-Color Ashkin-Teller Hamiltonian Interface            #
############################################################################

abstract type AbstractTwoColorAshkinTellerHamiltonian <: AbstractAshkinTeller end
abstract type AbstractTwoColorAshkinTellerParameters <: AbstractAshkinTellerParameters end
const TwoC_ATH = AbstractTwoColorAshkinTellerHamiltonian
const TwoC_AT_Params = AbstractTwoColorAshkinTellerParameters

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
    ham.color_update = state[end] <= num_colors(ham) ? to_AshkinTellerColor(state[end]) : ATDefaultColor(T)
    next_iter = iterations + one(Int)
    next_site = one(Int) + iterations % num_sites(ham)
    next_color_int = one(Int) + iterations ÷ num_sites(ham)
    return iterations <= length(ham) ? 
               ( (site_idx, ham[site_idx, color_update(ham)]), (next_site, next_iter, next_color_int) ) : 
               nothing
end

"""
    sigma_values(::AbstractTwoColorAshkinTellerHamiltonian)

Return a `view` into the σ degrees of freedom.
"""
@inline sigma_values(ham::AbstractTwoColorAshkinTellerHamiltonian) = @view spins(ham)[to_index(AT_sigma):num_colors(ham):end]

"""
    tau_values(::AbstractTwoColorAshkinTellerHamiltonian)

Return a `view` into the τ degrees of freedom.
"""
@inline tau_values(ham::AbstractTwoColorAshkinTellerHamiltonian)   = @view spins(ham)[to_index(AT_tau):num_colors(ham):end]

"""
    neighbor_fields(::AbstractTwoColorAshkinTellerHamiltonian, ::AshkinTellerParameters, latt, site)

Return the fields at a given `site` generated by the neighbors in a `latt`ice.
"""
@inline function neighbor_fields(ham::TwoC_ATH, hamparams::TwoC_AT_Params, latt, site)
    near_neighbors = nearest_neighbors(latt, site)
    σ_field = zero(eltype(ham))
    τ_field = σ_field
    bax_field = σ_field
    @inbounds for nn ∈ near_neighbors
        σ_field += ham[nn, AT_sigma]
        τ_field += ham[nn, AT_tau]
        bax_field += site_Baxter(ham, nn)
    end
    return @SVector [ hamparams.Jex * σ_field, hamparams.Jex * τ_field, hamparams.Kex * bax_field ]
end

@doc raw"""
    _base_site_energy(::AbstractTwoColorAshkinTellerHamiltonian, latt, site, site_values)

Calculate the energy for a single site based on the `3` degrees of freedom (DoF)
present (two primary, σ and τ, and one composite, [`site_Baxter`](@ref)). This 
_base_ energy is for a clean model defined by 

```math
H = -\sum_{\langle ij \rangle} \left[ J\left( \sigma_i\sigma_j + \tau_i\tau_j \right) + K\sigma_i\tau_i\sigma_j\tau_j \right].
```
"""
@inline function _base_site_energy( ham::TwoC_ATH, latt, site, site_values )
    effective_fields = neighbor_fields(ham, ham.params, latt, site)    
    en = zero(eltype(ham))
    @inbounds for (idx, eff_field) ∈ enumerate(effective_fields)
        en += site_values[idx] * eff_field
    end
    return -en
end

@doc raw"""
    _base_energy(::TwoC_ATH, latt)

Calculate the base energy for a pure Ashkin-Teller Hamiltonian given
on a given `latt`ice by the following Hamiltonian:

```math
H = -\sum_{\langle ij \rangle} \left[ J\left( \sigma_i\sigma_j + \tau_i\tau_j \right) + K\sigma_i\tau_i\sigma_j\tau_j \right].
```
"""
@inline function _base_energy( ham::TwoC_ATH, latt )
    en = zero(eltype(ham))
    @inbounds for (iter, dof_location_val) ∈ enumerate(IterateBySite, ham)
        en += site_energy( ham, latt, dof_location_val...)
    end
    return 0.5 * en
end

"""
    _base_DoF_energy_change(::AbstractTwoColorAshkinTellerHamiltonian, latt, site, [color = color_update(ham)])

Calculate the change in energy for a single `color` at a single site for a base two-color Ashkin-Teller Hamiltonian.
This function works by calling [`_base_site_energy`](@ref) while setting temporarily one of the primary field `site_values`
opposite to itself. Then, the _difference_ between the temporary copy of the DoF values and the real ones are used as the
new `site_values`.
"""
@inline function _base_DoF_energy_change(ham::TwoC_ATH, latt, site, color = color_update(ham))
    old_σ, old_τ, old_bax = ham[site, AT_sigma], ham[site, AT_tau], site_Baxter(ham, site)
    cond = color === AT_sigma
    σ_value = ham[site, AT_sigma] * ( cond ? -one(eltype(ham)) : one(eltype(ham)) )  # color === AT_sigma, cond == true
    τ_value = ham[site, AT_tau]   * ( !cond ? -one(eltype(ham)) : one(eltype(ham)) ) # color === AT_tau,   cond == false
    bax_val = σ_value * τ_value
    return site_energy( ham, latt, site, @SVector [ σ_value - old_σ, τ_value - old_τ, bax_val - old_bax ] )
end

"""
    _base_accept_move!(::Bool, ::AbstractTwoColorAshkinTellerHamiltonian, site)

Defines a binary operation of how to accept a single-DoF spin flip move if the `Bool`ean is 
evaluated to `true`. This `Bool`ean may be, for example, the Metropolis acceptance condition.
"""
@inline function _base_accept_move!( condition::Bool, ham::TwoC_ATH, site )
    ham[site, color_update(ham)] = condition ? -ham[site, color_update(ham)] : ham[site, color_update(ham)]
    return nothing
end