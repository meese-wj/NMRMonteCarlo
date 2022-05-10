using StaticArrays

abstract type AbstractHamiltonian end

@enum AT_Colors AT_σ=1 AT_τ=2 AT_enum_end
Base.getindex(collection, color::AT_Colors) = collection[Int(color)]
Base.getindex(collection::AbstractVector, color::AT_Colors) = collection[Int(color)]
Base.setindex!(collection, value, color::AT_Colors) = collection[Int(color)] = value
const NUM_AT_COLORS = Int(AT_enum_end) - Int(AT_σ)

struct AT_Parameters{T <: AbstractFloat}
    Jex::T  # Ising exchanges in the AT model. Jex > 0 is ferromagnetic
    Kex::T  # Baxter exchange measured in units of Jex
end
AT_Parameters{T}(; J::T = 1., K::T = 0. ) where {T <: AbstractFloat} = AT_Parameters{T}( J, K )
StructTypes.StructType(::Type{AT_Parameters{T}}) where {T} = StructTypes.Struct()
reciprocal_type(ty::Type{T}) where {T} = error("\nNo corresponding object defined for $(typeof(ty)) types.\n")
reciprocal_type(obj) = reciprocal_type(typeof(obj))
reciprocal_type(ty::Type{AT_Parameters{T}}) where {T} = AT_Hamiltonian{T}

mutable struct AT_Hamiltonian{T <: AbstractFloat} <: AbstractHamiltonian
    color_update::AT_Colors
    params::AT_Parameters{T}
    colors::Vector{T}
end

AT_Hamiltonian{T}( num_dof::Int; params = AT_Parameters{T}() ) where {T <: AbstractFloat} = AT_Hamiltonian( AT_σ, params, ones(T, num_dof) )
# AT_Hamiltonian{T}( latt::AbstractLattice; params = AT_Parameters{T}() ) where {T <: AbstractFloat} = AT_Hamiltonian{T}( NUM_AT_COLORS * num_sites(latt); params = params )
AT_Hamiltonian{T}( latt::AbstractLattice, params::AT_Parameters{T} ) where {T <: AbstractFloat} = AT_Hamiltonian{T}( NUM_AT_COLORS * num_sites(latt); params = params )

const default_params = String(joinpath(@__DIR__, "default_AT_Params.jl"))
function AT_Hamiltonian{T}( latt::AbstractLattice; 
                            params_file::String = joinpath(@__DIR__, "default_AT_Params.jl"),
                            display_io::IO = stdout ) where {T}
    at_params = import_json_and_display(params_file, AT_Parameters{T}, display_io )
    return AT_Hamiltonian{T}(latt, at_params)
end

Base.eltype(ham::AT_Hamiltonian) = eltype(ham.colors)
Base.length( ham::AT_Hamiltonian ) = length( ham.colors )
num_sites( ham::AT_Hamiltonian ) = length(ham) ÷ NUM_AT_COLORS
Base.size( ham::AT_Hamiltonian ) = ( num_sites(ham), num_sites(ham) )
num_DoF( ham::AT_Hamiltonian ) = length(ham)
color_index( site, col ) = convert(Int, NUM_AT_COLORS * (site - 1) + Int(col) )
Base.getindex( colors::AbstractVector, site, col::AT_Colors ) = getindex(colors, color_index(site, col))
Base.getindex( ham::AT_Hamiltonian, site, col::AT_Colors ) = ham.colors[ site, col ]
Base.setindex!( colors::AbstractVector, value, site, col::AT_Colors ) = setindex!(colors, value, color_index(site, col))
Base.setindex!(ham::AT_Hamiltonian, value, site, col::AT_Colors) = ham.colors[ site, col ] = value  
site_Baxter( colors::AbstractVector, site ) = colors[site, AT_σ] * colors[site, AT_τ]
site_Baxter( ham::AT_Hamiltonian, site ) = ham[site, AT_σ] * ham[site, AT_τ]
function switch_color_update!(ham::AT_Hamiltonian)
    if ham.color_update == AT_σ
        ham.color_update = AT_τ
        return nothing
    end
    ham.color_update = AT_σ
    return nothing
end

function AT_neighbor_fields(colors::AbstractVector, hamparams::AT_Parameters{T}, latt::AbstractLattice, site) where {T}
    near_neighbors = nearest_neighbors(latt, site)
    σ_field = zero(T)
    τ_field = σ_field
    bax_field = σ_field
    @inbounds for nn ∈ near_neighbors
        σ_field += colors[nn, AT_σ]
        τ_field += colors[nn, AT_σ]
        bax_field += site_Baxter(colors, nn)
    end
    return @SVector [ hamparams.Jex * σ_field, hamparams.Jex * τ_field, hamparams.Kex * bax_field ]
end

AT_neighbor_fields(ham::AT_Hamiltonian, latt::AbstractLattice, site) = AT_neighbor_fields(ham.colors, ham.params, latt, site)

function AT_site_energy( colors::AbstractVector, hamparams::AT_Parameters{T}, latt::AbstractLattice, site, site_values::SVector{3} ) where {T}
    effective_fields::SVector = AT_neighbor_fields(colors, hamparams, latt, site)    
    en = zero(T)
    @inbounds for idx ∈ 1:length(effective_fields)
        en += site_values[idx] * effective_fields[idx]
    end
    return -en
end

AT_site_energy(ham::AT_Hamiltonian, latt::AbstractLattice, site, site_values::SVector{3}) = AT_site_energy(ham.colors, ham.params, latt, site, site_values)

function AT_total_energy( colors::AbstractVector, hamparams::AT_Parameters{T}, latt::AbstractLattice ) where {T}
    en = AT_site_energy( colors, hamparams, latt, 1 )
    for site ∈ 2:num_sites(latt)
        en += AT_site_energy( colors, hamparams, latt, site)
    end
    return -0.5 * en
end

AT_total_energy(ham::AT_Hamiltonian, latt::AbstractLattice) = AT_total_energy(ham.colors, ham.params, latt)

AT_site_energy(ham::AT_Hamiltonian, latt::AbstractLattice, site) = AT_site_energy( ham, latt, site, @SVector [ ham[site, AT_σ], ham[site, AT_τ], site_Baxter(ham, site) ] )
function AT_site_energy_change(ham::AT_Hamiltonian, latt::AbstractLattice, site, color = AT_σ)
    σ_value = -ham[site, AT_σ]
    τ_value = ham[site, AT_τ]
    if color == AT_τ
        σ_value *= -one(eltype(ham))
        τ_value *= -one(eltype(ham))
    end
    bax_val = σ_value * τ_value
    return AT_site_energy( ham, latt, site, @SVector [ σ_value - ham[site, AT_σ], τ_value - ham[site, AT_τ], bax_val - site_Baxter(ham, site) ] )
end

function AT_site_flip!( ham::AT_Hamiltonian, site )
    ham[site, ham.color_update] *= -one(ham[site, ham.color_update])
    return nothing
end

function export_state!(state_container::AbstractArray, ham::AT_Hamiltonian, export_index)
    state_container[:, export_index] .= ham.colors
    return nothing 
end

load_state!(ham::AT_Hamiltonian, state) = ham.colors .= state
