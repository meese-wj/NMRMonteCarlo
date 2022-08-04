# Main code for concrete AT Hamiltonians

using Parameters2JSON
import StaticArrays: @SVector, SVector
import ..Lattices: nearest_neighbors
export AshkinTellerHamiltonian, AshkinTellerParameters, sigma_values, tau_values

@jsonable struct AshkinTellerParameters{T <: AbstractFloat}
    Jex::T  # Ising exchanges in the AT model. Jex > 0 is ferromagnetic
    Kex::T  # Baxter exchange measured in units of Jex
end
AshkinTellerParameters{T}(; J::T = 1., K::T = 0. ) where {T <: AbstractFloat} = AshkinTellerParameters{T}( J, K )
reciprocal_type(ty::Type{T}) where {T} = error("\nNo corresponding object defined for $(typeof(ty)) types.\n")
reciprocal_type(obj) = reciprocal_type(typeof(obj))
reciprocal_type(ty::Type{AshkinTellerParameters{T}}) where {T} = AshkinTellerHamiltonian{T}

mutable struct AshkinTellerHamiltonian{T <: AbstractFloat} <: AbstractTwoColorAshkinTellerHamiltonian
    color_update::AshkinTellerColor
    params::AshkinTellerParameters{T}
    spins::Vector{T}

    function AshkinTellerHamiltonian(latt, params::AshkinTellerParameters{T}) where T
        ndofs = num_colors(AshkinTellerHamiltonian) * num_sites(latt)
        return new{T}(ATDefaultColor(AshkinTellerHamiltonian), params, rand([one(T), -one(T)], ndofs) )
    end
end

@inline sigma_values(ham::AbstractTwoColorAshkinTellerHamiltonian) = @view spins(ham)[to_index(AT_sigma):num_colors(ham):end]
@inline tau_values(ham::AbstractTwoColorAshkinTellerHamiltonian)   = @view spins(ham)[to_index(AT_tau):num_colors(ham):end]

@inline function neighbor_fields(ham::TwoC_ATH, hamparams::AshkinTellerParameters{T}, latt, site) where {T}
    near_neighbors = nearest_neighbors(latt, site)
    σ_field = zero(T)
    τ_field = σ_field
    bax_field = σ_field
    @inbounds for nn ∈ near_neighbors
        σ_field += ham[nn, AT_sigma]
        τ_field += ham[nn, AT_tau]
        bax_field += site_Baxter(ham, nn)
    end
    return @SVector [ hamparams.Jex * σ_field, hamparams.Jex * τ_field, hamparams.Kex * bax_field ]
end

@inline function DoF_energy( ham::AshkinTellerHamiltonian{T}, latt, site, site_values::SVector{3} ) where {T}
    effective_fields::SVector{3, T} = neighbor_fields(ham, ham.params, latt, site)    
    en = zero(T)
    @inbounds for (idx, eff_field) ∈ enumerate(effective_fields)
        en += site_values[idx] * eff_field
    end
    return -en
end

@inline function energy( ham::AshkinTellerHamiltonian{T}, latt ) where {T}
    en = zero(T)
    @inbounds for (iter, dof_location_val) ∈ enumerate(IterateBySite, ham)
        en += DoF_energy( ham, latt, dof_location_val...)
    end
    return 0.5 * en
end

@inline function DoF_energy_change(ham::AshkinTellerHamiltonian{T}, latt, site, color = color_update(ham)) where T
    old_σ, old_τ, old_bax = ham[site, AT_sigma], ham[site, AT_tau], site_Baxter(ham, site)
    # σ_value = ifelse( color == AT_sigma, -ham[site, AT_sigma], ham[site, AT_sigma] )
    # τ_value = ifelse( color == AT_tau, -ham[site, AT_tau], ham[site, AT_tau] )
    cond = color === AT_sigma
    σ_value = ham[site, AT_sigma] * ( cond ? -one(T) : one(T) )  # color === AT_sigma, cond == true
    τ_value = ham[site, AT_tau]   * ( !cond ? -one(T) : one(T) ) # color === AT_tau,   cond == false
    # σ_value = -ham[site, AT_sigma]
    # τ_value = ham[site, AT_tau]
    # if color === AT_tau
    #     # σ_value *= -one(eltype(ham))
    #     # τ_value *= -one(eltype(ham))
    #     σ_value *= -one(T)
    #     τ_value *= -one(T)
    # end
    bax_val = σ_value * τ_value
    # return DoF_energy( ham, latt, site, @SVector [ σ_value - ham[site, AT_sigma], τ_value - ham[site, AT_tau], bax_val - site_Baxter(ham, site) ] )
    return DoF_energy( ham, latt, site, @SVector [ σ_value - old_σ, τ_value - old_τ, bax_val - old_bax ] )
end

@inline function site_flip!( condition::Bool, ham::AshkinTellerHamiltonian, site )
    ham[site, color_update(ham)] = condition ? -ham[site, color_update(ham)] : ham[site, color_update(ham)]
    return nothing
end

function export_state!(state_container::AbstractArray, ham::AshkinTellerHamiltonian, export_index)
    state_container[:, export_index] .= ham.colors
    return nothing 
end

load_state!(ham::AshkinTellerHamiltonian, state) = ham.colors .= state
