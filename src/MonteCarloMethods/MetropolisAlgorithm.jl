using Parameters2JSON

export MetropolisParameters, thermalize!, sweep_and_measure!

@jsonable struct MetropolisParameters{T <: AbstractFloat} <: AbstractMonteCarloParameters
    βvalues::Vector{T}
    therm_sweeps::Int
    measure_sweeps::Int
    total_measurements::Int
end

function MetropolisParameters{T}(; params_file = joinpath(@__DIR__, "default_Metropolis_Params.jl"),
                                   display_io::IO = stdout) where {T}
    return import_json_and_display(params_file, MetropolisParameters{T}, display_io)
end

metropolis_accepted(ΔE, β) = ( ΔE < zero(ΔE) || rand() < exp(-β * ΔE) )::Bool

# function metropolis_update!( model::Model{L, AT_Hamiltonian{T}}, site, beta ) where {L <: AbstractLattice, T <: AbstractFloat}
function metropolis_update!( model::AbstractModel, site, beta )
    # ΔE = AT_site_energy_change( model.ham, model.latt, site, model.ham.color_update )
    ham = hamiltonian(model)
    ΔE = AT_site_energy_change( ham, lattice(model), site, ham.color_update )
    if metropolis_accepted( ΔE, beta )
        # AT_site_flip!( model.ham, site )
        AT_site_flip!( ham, site )
    end
    return nothing
end

# function metropolis_sweep!( model::Model{L, H}, beta ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
function metropolis_sweep!( model::AbstractModel, beta )

    # # TODO: Generalize this loop to something like this:
    # for (dof_idx, dof_val) ∈ enumerate( hamiltonian(model), TraverseByDoF )
    #     metropolis_update!(model, beta, dof_idx, dof_val)
    # end

    for dof_color ∈ 1:NUM_AT_COLORS
        # for site ∈ 1:num_sites(model.latt)
        for site ∈ 1:num_sites(lattice(model))
            metropolis_update!( model, site, beta )
        end
        # switch_color_update!(model.ham)
        switch_color_update!(hamiltonian(model))
    end
    return nothing
end

# function thermalize!( model::Model{L, H}, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
function thermalize!( model::AbstractModel, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function )
    total_sweeps::Int64 = mc_params.therm_sweeps
    @inbounds for sweep ∈ 1:total_sweeps
        mc_sweep(model, beta)
    end
    return nothing
end

# function sweep_and_measure!( model::Model{L, H}, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function, state_container::AbstractArray = [] ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
function sweep_and_measure!( model::AbstractModel, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function, state_container::AbstractArray = [] )
    total_exports = mc_params.total_measurements
    spe = sweeps_per_export(mc_params)
    @inbounds for write ∈ (1:total_exports)
        @inbounds for sweep ∈ (1:spe)
            mc_sweep(model, beta)
        end

        # Only write out observables if the state container is defined
        if length(state_container) > 0
            # export_state!( state_container, model.ham, write )
            export_state!( state_container, hamiltonian(model), write )
        end
    end
    return nothing
end