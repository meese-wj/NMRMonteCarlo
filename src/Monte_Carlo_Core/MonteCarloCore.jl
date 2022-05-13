struct Model{L <: AbstractLattice, H <: AbstractHamiltonian}
    latt::L
    ham::H
end

abstract type MonteCarloParameters end
sweeps_per_export(params::MonteCarloParameters) = params.measure_sweeps <= params.total_measurements ? 1 : params.measure_sweeps ÷ params.total_measurements

struct MetropolisParameters{T <: AbstractFloat} <: MonteCarloParameters
    βvalues::Vector{T}
    therm_sweeps::Int
    measure_sweeps::Int
    total_measurements::Int
end
StructTypes.StructType(::Type{MetropolisParameters{T}}) where {T} = StructTypes.Struct()

function MetropolisParameters{T}(; params_file = joinpath(@__DIR__, "default_Metropolis_Params.jl"),
                                   display_io::IO = stdout) where {T}
    return import_json_and_display(params_file, MetropolisParameters{T}, display_io)
end

metropolis_accepted(ΔE, β) = ( ΔE < zero(ΔE) || rand() < exp(-β * ΔE) )::Bool

function metropolis_update!( model::Model{L, AT_Hamiltonian{T}}, site, beta ) where {L <: AbstractLattice, T <: AbstractFloat}
    ΔE = AT_site_energy_change( model.ham, model.latt, site, model.ham.color_update )
    if metropolis_accepted( ΔE, beta )
        AT_site_flip!( model.ham, site )
    end
    return nothing
end

function metropolis_sweep!( model::Model{L, H}, beta ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
    for dof_color ∈ 1:NUM_AT_COLORS
        for site ∈ 1:num_sites(model.latt)
            metropolis_update!( model, site, beta )
        end
        switch_color_update!(model.ham)
    end
    return nothing
end

function thermalize!( model::Model{L, H}, beta, mc_params::MonteCarloParameters, mc_sweep::Function ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
    total_sweeps::Int64 = mc_params.therm_sweeps
    @inbounds for sweep ∈ 1:total_sweeps
        mc_sweep(model, beta)
    end
    return nothing
end

function sweep_and_measure!( model::Model{L, H}, beta, mc_params::MonteCarloParameters, mc_sweep::Function, state_container::AbstractArray = [] ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
    total_exports = mc_params.total_measurements
    sue = sweeps_per_export(mc_params)
    @inbounds for write ∈ (1:total_exports)
        @inbounds for sweep ∈ (1:sue)
            mc_sweep(model, beta)
        end

        # Only write out observables if the state container is defined
        if length(state_container) > 0
            export_state!( state_container, model.ham, write )
        end
    end
    return nothing
end