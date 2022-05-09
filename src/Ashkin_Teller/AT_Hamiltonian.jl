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
Base.getindex( ham::AT_Hamiltonian, site, col ) = ham.colors[ color_index(site, col) ]
Base.setindex!(ham::AT_Hamiltonian, value, site, col) = ham.colors[ color_index(site, col) ] = value  
site_Baxter( ham::AT_Hamiltonian, site ) = ham[site, AT_σ] * ham[site, AT_τ]
function switch_color_update!(ham::AT_Hamiltonian)
    if ham.color_update == AT_σ
        ham.color_update = AT_τ
        return nothing
    end
    ham.color_update = AT_σ
    return nothing
end

using StaticArrays

function AT_neighbor_fields(ham::AT_Hamiltonian, latt::AbstractLattice, site)
    near_neighbors = nearest_neighbors(latt, site)
    # (σ_field::Float64, τ_field::Float64, bax_field::Float64) = ( zero(Float64), zero(Float64), zero(Float64) )
    σ_field = zero(eltype(ham))
    τ_field = σ_field
    bax_field = σ_field
    for nn ∈ near_neighbors
        σ_field += ham[nn, AT_σ]
        τ_field += ham[nn, AT_σ]
        bax_field += site_Baxter(ham, nn)
    end
    return @SVector [ ham.params.Jex * σ_field, ham.params.Jex * τ_field, ham.params.Kex * bax_field ]
    # return ( ham.params.Jex * σ_field, ham.params.Jex * τ_field, ham.params.Kex * bax_field )
end

function AT_site_energy( ham::AT_Hamiltonian, latt::AbstractLattice, site, site_values::SVector{3} )
    effective_fields::SVector = AT_neighbor_fields(ham, latt, site)    
    en = zero(eltype(effective_fields))
    @inbounds for idx ∈ 1:length(effective_fields)
        en += site_values[idx] * effective_fields[idx]
    end
    # en = site_values[AT_σ] * effective_fields[AT_σ] + site_values[AT_τ] * effective_fields[AT_τ] + site_values[end] * effective_fields[end]
    return -en
end

function AT_total_energy( ham, latt )
    en = AT_site_energy( ham, latt, 1 )
    for site ∈ 2:num_sites(latt)
        en += AT_site_energy( ham, latt, site)
    end
    return -0.5 * en
end

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

state_file_name(ham::AT_Hamiltonian, sweep_number) = "sweep-$sweep_number.bin"

function write_state(ham::AT_Hamiltonian, sweep_number, data_dir::String)
    data_path = joinpath(data_dir, state_file_name(ham, sweep_number))
    write(data_path, ham.colors)
    return nothing
end