# Main code for concrete AT Hamiltonians

using Parameters2JSON
import StaticArrays: @SVector, SVector
import ..Lattices: nearest_neighbors
export BasicAshkinTellerHamiltonian, AshkinTellerParameters, sigma_values, tau_values

@jsonable struct AshkinTellerParameters{T <: AbstractFloat} <: AbstractTwoColorAshkinTellerParameters
    Jex::T  # Ising exchanges in the AT model. Jex > 0 is ferromagnetic
    Kex::T  # Baxter exchange measured in units of Jex
end
AshkinTellerParameters{T}(; J::T = 1., K::T = 0. ) where {T <: AbstractFloat} = AshkinTellerParameters{T}( J, K )
reciprocal_type(ty::Type{T}) where {T} = error("\nNo corresponding object defined for $(typeof(ty)) types.\n")
reciprocal_type(obj) = reciprocal_type(typeof(obj))
reciprocal_type(ty::Type{AshkinTellerParameters{T}}) where {T} = BasicAshkinTellerHamiltonian{T}

mutable struct BasicAshkinTellerHamiltonian{T <: AbstractFloat} <: AbstractTwoColorAshkinTellerHamiltonian
    color_update::AshkinTellerColor
    params::AshkinTellerParameters{T}
    spins::Vector{T}

    function BasicAshkinTellerHamiltonian(latt, params::AshkinTellerParameters{T}) where T
        ndofs = num_colors(BasicAshkinTellerHamiltonian) * num_sites(latt)
        return new{T}(ATDefaultColor(BasicAshkinTellerHamiltonian), params, rand([one(T), -one(T)], ndofs) )
    end
end

@inline site_energy(ham::BasicAshkinTellerHamiltonian, latt, site, site_values) = _base_site_energy(ham, latt, site, site_values)

@inline energy(ham::BasicAshkinTellerHamiltonian, latt) = _base_energy(ham, latt)

@inline DoF_energy_change(ham::BasicAshkinTellerHamiltonian, latt, site, color = color_update(ham)) = _base_DoF_energy_change(ham, latt, site, color)

@inline accept_move!(condition::Bool, ham::BasicAshkinTellerHamiltonian, site) = _base_accept_move!(condition, ham, site)

function export_state!(state_container::AbstractArray, ham::BasicAshkinTellerHamiltonian, export_index)
    state_container[:, export_index] .= ham.colors
    return nothing 
end

load_state!(ham::BasicAshkinTellerHamiltonian, state) = ham.colors .= state
