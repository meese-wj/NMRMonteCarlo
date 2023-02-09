# Main code for concrete AT Hamiltonians

using Parameters2JSON
import StaticArrays: @SVector, SVector
import ..Lattices: nearest_neighbors
export BasicAshkinTellerHamiltonian, AshkinTellerParameters, sigma_values, tau_values

@jsonable struct AshkinTellerParameters{T <: AbstractFloat} <: AbstractTwoColorAshkinTellerParameters
    Jex::T  # Ising exchanges in the AT model. Jex > 0 is ferromagnetic
    Kex::T  # Baxter exchange measured in units of Jex
end

struct RandomBaxterFieldParameters{T <: AbstractFloat, D <: DisorderDistribution} <: AbstractRandomBaxterFieldParameters
    Jex::T  # Ising exchanges in the AT model. Jex > 0 is ferromagnetic
    Kex::T  # Baxter exchange measured in units of Jex
    Δε::T   # Width of the Baxter field disorder distribution
    disorder_distribution::D # Type of Baxter field distribution
end
AshkinTellerParameters{T}(; J::T = 1., K::T = 0. ) where {T <: AbstractFloat} = AshkinTellerParameters{T}( J, K )
RandomBaxterFieldParameters(; J = 1., K = 0., Δε = 0. ) = ( args = promote(J, K, Δε); AshkinTellerParameters{eltype(args)}( args... ) )

# begin DEPRECATION
reciprocal_type(ty::Type{T}) where {T} = error("\nNo corresponding object defined for $(typeof(ty)) types.\n")
reciprocal_type(obj) = reciprocal_type(typeof(obj))
reciprocal_type(ty::Type{AshkinTellerParameters{T}}) where {T} = BasicAshkinTellerHamiltonian{T}
# End

mutable struct BasicAshkinTellerHamiltonian{T <: AbstractFloat} <: AbstractTwoColorAshkinTellerHamiltonian
    color_update::AshkinTellerColor
    params::AshkinTellerParameters{T}
    spins::Vector{T}

    function BasicAshkinTellerHamiltonian(latt, params::AshkinTellerParameters{T}) where T
        ndofs = num_colors(BasicAshkinTellerHamiltonian) * num_sites(latt)
        return new{T}(ATDefaultColor(BasicAshkinTellerHamiltonian), params, rand([one(T), -one(T)], ndofs) )
    end
end

mutable struct RandomBaxterFieldHamiltonian{T <: AbstractFloat} <: AbstractRandomBaxterFieldHamiltonian
    color_update::AshkinTellerColor
    params::RandomBaxterFieldParameters{T}
    spins::Vector{T}
    random_baxter_fields::Vector{T}

    function RandomBaxterFieldHamiltonian(latt, params::RandomBaxterFieldParameters{T}, rng::Random.AbstractRNG = Random.GLOBAL_RNG) where T
        ndofs = num_colors(RandomBaxterFieldHamiltonian) * num_sites(latt)
        qgd = QuenchedDisorderGenerator(zero(params.Δε), params.Δε, rng, params.disorder_distribution)
        fields = generate_disorder(qgd, num_sites(latt))
        return new{T}(ATDefaultColor(RandomBaxterFieldHamiltonian), params, rand([one(T), -one(T)], ndofs), fields)
    end
end

@inline site_energy(ham::BasicAshkinTellerHamiltonian, latt, site, site_values) = _base_site_energy(ham, latt, site, site_values)
@inline site_energy(ham::RandomBaxterFieldHamiltonian, latt, site, site_values) = _base_site_energy(ham, latt, site, site_values) + _Baxter_site_energy(ham, site, site_values[end])

@inline energy(ham::BasicAshkinTellerHamiltonian, latt) = _base_energy(ham, latt)
function energy(ham::RandomBaxterFieldHamiltonian, latt)
    en = _base_energy(ham, latt)
    @inbounds for (site, dof_location_val) ∈ enumerate(IterateBySite, ham)
        en += _Baxter_site_energy(ham, dof_location_val[1], dof_location_val[2][end])
    end
    return en
end

@inline DoF_energy_change(ham::BasicAshkinTellerHamiltonian, latt, site, color = color_update(ham)) = _base_DoF_energy_change(ham, latt, site, color)
@inline DoF_energy_change(ham::RandomBaxterFieldHamiltonian, latt, site, color = color_update(ham)) = _base_DoF_energy_change(ham, latt, site, color) + _Baxter_DoF_energy_change(ham, site)

@inline accept_move!(condition::Bool, ham::BasicAshkinTellerHamiltonian, site) = _base_accept_move!(condition, ham, site)

function export_state!(state_container::AbstractArray, ham::BasicAshkinTellerHamiltonian, export_index)
    state_container[:, export_index] .= ham.colors
    return nothing 
end

load_state!(ham::BasicAshkinTellerHamiltonian, state) = ham.colors .= state
