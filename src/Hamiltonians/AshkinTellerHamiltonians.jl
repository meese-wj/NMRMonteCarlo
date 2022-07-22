
using Parameters2JSON
import StaticArrays: @SVector, SVector
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

@jsonable struct AshkinTellerParameters{T <: AbstractFloat}
    Jex::T  # Ising exchanges in the AT model. Jex > 0 is ferromagnetic
    Kex::T  # Baxter exchange measured in units of Jex
end
AshkinTellerParameters{T}(; J::T = 1., K::T = 0. ) where {T <: AbstractFloat} = AshkinTellerParameters{T}( J, K )
reciprocal_type(ty::Type{T}) where {T} = error("\nNo corresponding object defined for $(typeof(ty)) types.\n")
reciprocal_type(obj) = reciprocal_type(typeof(obj))
reciprocal_type(ty::Type{AshkinTellerParameters{T}}) where {T} = AshkinTellerHamiltonian{T}

mutable struct AshkinTellerHamiltonian{T <: AbstractFloat} <: AbstractTwoColorAshkinTellerHamiltonian
    color_update::Type{<: AshkinTellerColor}
    params::AshkinTellerParameters{T}
    spins::Vector{T}
end

AshkinTellerHamiltonian{T}( num_dof::Int; params = AshkinTellerParameters{T}() ) where {T <: AbstractFloat} = AshkinTellerHamiltonian( AT_sigma, params, rand([one(T), -one(T)], num_dof) )
# AshkinTellerHamiltonian{T}( latt::AbstractLattice; params = AshkinTellerParameters{T}() ) where {T <: AbstractFloat} = AshkinTellerHamiltonian{T}( num_colors(AshkinTellerHamiltonian) * num_sites(latt); params = params )
# AshkinTellerHamiltonian{T}( latt::AbstractLattice, params::AshkinTellerParameters{T} ) where {T <: AbstractFloat} = AshkinTellerHamiltonian{T}( num_colors(AshkinTellerHamiltonian) * num_sites(latt); params = params )
AshkinTellerHamiltonian{T}( latt, params::AshkinTellerParameters{T} ) where {T <: AbstractFloat} = AshkinTellerHamiltonian{T}( num_colors(AshkinTellerHamiltonian) * num_sites(latt); params = params )

const default_params = String(joinpath(@__DIR__, "default_AT_Params.jl"))
function AshkinTellerHamiltonian{T}( latt; # ::AbstractLattice
                            params_file::String = joinpath(@__DIR__, "default_AT_Params.jl"),
                            display_io::IO = stdout ) where {T}
    at_params = import_json_and_display(params_file, AshkinTellerParameters{T}, display_io )
    return AshkinTellerHamiltonian{T}(latt, at_params)
end
  
site_Baxter( ham::AbstractTwoColorAshkinTellerHamiltonian, site ) = ham[site, AT_sigma] * ham[site, AT_tau]
function switch_color_update!(ham::AbstractTwoColorAshkinTellerHamiltonian)
    # if ham.color_update == AT_sigma
    #     ham.color_update = AT_tau
    #     return nothing
    # end
    # ham.color_update = AT_sigma
    ifelse(ham.color_update == AT_sigma, AT_tau, AT_sigma)
    return nothing
end

# TODO: get rid of this maybe? Only use the ham::AbstractTwoColorAshkinTellerHamiltonian version?
# function neighbor_fields(colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt::AbstractLattice, site) where {T}
function neighbor_fields(colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt, site) where {T}
    near_neighbors = nearest_neighbors(latt, site)
    σ_field = zero(T)
    τ_field = σ_field
    bax_field = σ_field
    @inbounds for nn ∈ near_neighbors
        σ_field += colors[nn, AT_sigma]
        τ_field += colors[nn, AT_tau]
        bax_field += site_Baxter(colors, nn)
    end
    return @SVector [ hamparams.Jex * σ_field, hamparams.Jex * τ_field, hamparams.Kex * bax_field ]
end

# neighbor_fields(ham::AshkinTellerHamiltonian, latt::AbstractLattice, site) = neighbor_fields(ham.colors, ham.params, latt, site)
neighbor_fields(ham::AshkinTellerHamiltonian, latt, site) = neighbor_fields(ham.colors, ham.params, latt, site)

# TODO: get rid of this maybe? Only use the ham::AbstractTwoColorAshkinTellerHamiltonian version?
# function DoF_energy( colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt::AbstractLattice, site, site_values::SVector{3} ) where {T}
function DoF_energy( colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt, site, site_values::SVector{3} ) where {T}
    effective_fields::SVector = neighbor_fields(colors, hamparams, latt, site)    
    en = zero(T)
    @inbounds for idx ∈ 1:length(effective_fields)
        en += site_values[idx] * effective_fields[idx]
    end
    return -en
end

# DoF_energy(ham::AshkinTellerHamiltonian, latt::AbstractLattice, site, site_values::SVector{3}) = DoF_energy(ham.colors, ham.params, latt, site, site_values)
DoF_energy(ham::AbstractAshkinTeller, latt, site, site_values::SVector{3}) = DoF_energy(ham.colors, ham.params, latt, site, site_values)
# TODO: get rid of this maybe? Only use the ham::AbstractTwoColorAshkinTellerHamiltonian version?
# DoF_energy(colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt::AbstractLattice, site) where {T} = DoF_energy(colors, hamparams, latt, site, @SVector [ colors[site, AT_sigma], colors[site, AT_tau], site_Baxter(colors, site) ])
DoF_energy(colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt, site) where {T} = DoF_energy(colors, hamparams, latt, site, @SVector [ colors[site, AT_sigma], colors[site, AT_tau], site_Baxter(colors, site) ])

# TODO: get rid of this maybe? Only use the ham::AbstractTwoColorAshkinTellerHamiltonian version?
# function energy( colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt::AbstractLattice ) where {T}
function energy( colors::AbstractVector, hamparams::AshkinTellerParameters{T}, latt ) where {T}
    en = DoF_energy( colors, hamparams, latt, 1 )
    for site ∈ 2:num_sites(latt)
        en += DoF_energy( colors, hamparams, latt, site)
    end
    return 0.5 * en
end

# DoF_energy(ham::AbstractAshkinTeller, latt::AbstractLattice, site) = DoF_energy( ham, latt, site, @SVector [ ham[site, AT_sigma], ham[site, AT_tau], site_Baxter(ham, site) ] )
DoF_energy(ham::AbstractAshkinTeller, latt, site) = DoF_energy( ham, latt, site, @SVector [ ham[site, AT_sigma], ham[site, AT_tau], site_Baxter(ham, site) ] )
# energy(ham::AbstractAshkinTeller, latt::AbstractLattice) = energy(ham.colors, ham.params, latt)
energy(ham::AbstractAshkinTeller, latt) = energy(ham.colors, ham.params, latt)

# function DoF_energy_change(ham::AshkinTellerHamiltonian, latt::AbstractLattice, site, color = AT_sigma)
function DoF_energy_change(ham::AshkinTellerHamiltonian, latt, site, color = AT_sigma)
    σ_value = -ham[site, AT_sigma]
    τ_value = ham[site, AT_tau]
    if color == AT_tau
        σ_value *= -one(eltype(ham))
        τ_value *= -one(eltype(ham))
    end
    bax_val = σ_value * τ_value
    return DoF_energy( ham, latt, site, @SVector [ σ_value - ham[site, AT_sigma], τ_value - ham[site, AT_tau], bax_val - site_Baxter(ham, site) ] )
end

function site_flip!( ham::AshkinTellerHamiltonian, site )
    ham[site, ham.color_update] *= -one(ham[site, ham.color_update])
    return nothing
end

function export_state!(state_container::AbstractArray, ham::AshkinTellerHamiltonian, export_index)
    state_container[:, export_index] .= ham.colors
    return nothing 
end

load_state!(ham::AshkinTellerHamiltonian, state) = ham.colors .= state
