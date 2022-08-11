# Main code for concrete AT Hamiltonians

using Parameters2JSON
import StaticArrays: @SVector, SVector
import ..Lattices: nearest_neighbors
export AshkinTellerHamiltonian, AshkinTellerParameters, sigma_values, tau_values

@jsonable struct AshkinTellerParameters{T <: AbstractFloat} <: AbstractTwoColorAshkinTellerParameters
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

@inline site_energy(ham::AshkinTellerHamiltonian{T}, latt, site, site_values::SVector{3}) where {T} = _base_site_energy(ham, latt, site, site_values)

@inline function energy( ham::AshkinTellerHamiltonian{T}, latt ) where {T}
    en = zero(T)
    @inbounds for (iter, dof_location_val) ∈ enumerate(IterateBySite, ham)
        en += site_energy( ham, latt, dof_location_val...)
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
    # return site_energy( ham, latt, site, @SVector [ σ_value - ham[site, AT_sigma], τ_value - ham[site, AT_tau], bax_val - site_Baxter(ham, site) ] )
    return site_energy( ham, latt, site, @SVector [ σ_value - old_σ, τ_value - old_τ, bax_val - old_bax ] )
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
