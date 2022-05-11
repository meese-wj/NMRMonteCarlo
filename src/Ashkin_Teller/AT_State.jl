
build_state_container(ham_t::Type, num_dof::Int, num_state_writes::Int ) = error("\nNo container designed for Hamiltonians of type $(typeof(ham_t)).")

function build_state_container(ham_t::Type{AT_Hamiltonian{T}}, num_dof::Int, num_state_writes::Int ) where {T}
    return zeros(T, (num_dof, num_state_writes))
end

const int8spinminus = -one(Int8)
const int8spinplus = one(Int8)

function ising_to_int8(arraylike)
    int8vector = Vector{typeof(int8spinplus)}(undef, length(arraylike))
    @inbounds @simd for idx in eachindex(arraylike)
        int8vector[idx] = arraylike[idx] == one(arraylike[idx]) ? int8spinplus : int8spinminus
    end
    return int8vector
end

function int8_to_ising!(ising_arraylike, arraylike)
    length(ising_arraylike) == length(arraylike) ? nothing : error("Ising vector has length $(length(ising_arraylike)) but the Int8 vector has length $(length(arraylike)).")
    @inbounds @simd for idx in eachindex(ising_arraylike)
        ising_arraylike[idx] = arraylike[idx] == int8spinplus ? one(ising_arraylike[idx]) : -one(ising_arraylike[idx])
    end
    return nothing
end

function int8_to_ising(arraylike, type::Type = Float64)
    ising_vector = Vector{type}(undef, length(arraylike))
    int8_to_ising!(ising_vector, arraylike)
    return ising_vector
end

state_file_name(identifiers...) = "states$(concat_ids(identifiers...)).bin"

function export_states(state_container::AbstractArray, data_dir::String, identifiers...)
    data_path = joinpath(data_dir, state_file_name(identifiers...))
    open(data_path, "w") do io
        write(io, ising_to_int8(state_container))
    end
    return nothing
end

import_states(ty::Type, data_file, num_dofs, num_states) = error("\nNo method defined for importing states for types $(ty).")

function import_states(ty::Type{AT_Hamiltonian{T}}, data_file, num_dofs, num_states) where {T}
    states = zeros(Int8, num_dofs * num_states)
    read!(data_file, states)
    ising_states = int8_to_ising(states, T)
    return reshape(ising_states, num_dofs, num_states)
end

function separate_colors!(::Type{AT_Hamiltonian{T}}, colors::Matrix{T}, state::AbstractVector) where {T}
    length(colors[:,Int(AT_sigma)]) == length(state) ? nothing : error("\nDestination matrix and state have different lengths: 2 × $(length(colors[:,Int(AT_sigma)])) ≂̸ $(length(state)).")
    @inbounds for site ∈ 1:(length(state) ÷ NUM_AT_COLORS)
        colors[site, Int(AT_sigma)] = state[ color_index(site, AT_sigma) ]
        colors[site, Int(AT_tau)] = state[ color_index(site, AT_tau) ]
    end
    return nothing
end

function separate_colors_from_states(ty::Type{AT_Hamiltonian{T}}, states::AbstractArray) where {T}
    num_color, num_states = size(states)
    state_colors = Array{T, 3}(undef, num_color ÷ NUM_AT_COLORS, NUM_AT_COLORS, num_states)
    @inbounds for state_idx ∈ 1:num_states
        separate_colors!(ty, state_colors[:,:,state_idx], states[:, state_idx])
    end
    return state_colors
end