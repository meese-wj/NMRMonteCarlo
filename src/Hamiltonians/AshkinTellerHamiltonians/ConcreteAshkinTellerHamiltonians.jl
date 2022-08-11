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

@inline site_energy(ham::AshkinTellerHamiltonian, latt, site, site_values) = _base_site_energy(ham, latt, site, site_values)

@inline energy(ham::AshkinTellerHamiltonian, latt) = _base_energy(ham, latt)

@inline DoF_energy_change(ham::AshkinTellerHamiltonian, latt, site, color = color_update(ham)) = _base_DoF_energy_change(ham, latt, site, color)

@inline function site_flip!( condition::Bool, ham::AshkinTellerHamiltonian, site )
    ham[site, color_update(ham)] = condition ? -ham[site, color_update(ham)] : ham[site, color_update(ham)]
    return nothing
end

function export_state!(state_container::AbstractArray, ham::AshkinTellerHamiltonian, export_index)
    state_container[:, export_index] .= ham.colors
    return nothing 
end

load_state!(ham::AshkinTellerHamiltonian, state) = ham.colors .= state
