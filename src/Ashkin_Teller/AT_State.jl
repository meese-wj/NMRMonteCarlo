
build_state_container(ham_t::Type, num_dof::Int, num_state_writes::Int ) = error("\nNo container designed for Hamiltonians of type $(typeof(ham_t)).")

function build_state_container(ham_t::Type{AT_Hamiltonian{T}}, num_dof::Int, num_state_writes::Int ) where {T}
    return zeros(T, (num_dof, num_state_writes))
end