struct Model{L <: AbstractLattice, H <: AbstractHamiltonian}
    latt::L
    ham::H
end

abstract type MonteCarloParameters end

struct MetropolisParameters{T <: AbstractFloat} <: MonteCarloParameters
    β::T
    therm_sweeps::Int
    measure_sweeps::Int
    sweeps_until_write::Int
end
StructTypes.StructType(::Type{MetropolisParameters{T}}) where {T} = StructTypes.Struct()

function MetropolisParameters{T}(; params_file = joinpath(@__DIR__, "default_Metropolis_Params.jl"),
                                   display_io::IO = stdout) where {T}
    return import_json_and_display(params_file, MetropolisParameters{T}, display_io)
end

metropolis_accepted(ΔE, β) = ( ΔE < zero(ΔE) || rand() < exp(-β * ΔE) )::Bool

function metropolis_update!( model::Model{L, AT_Hamiltonian{T}}, site, mc_params::MetropolisParameters ) where {L <: AbstractLattice, T <: AbstractFloat}
    ΔE = AT_site_energy_change( model.ham, model.latt, site, model.ham.color_update )
    if metropolis_accepted( ΔE, mc_params.β )
        AT_site_flip!( model.ham, site )
    end
    return nothing
end

function metropolis_sweep!( model::Model{L, H}, mc_params::MonteCarloParameters ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
    for dof_color ∈ 1:NUM_AT_COLORS
        for site ∈ 1:num_sites(model.latt)
            metropolis_update!( model, site, mc_params )
        end
        switch_color_update!(model.ham)
    end
    return nothing
end

function thermalize!( model::Model{L, H}, mc_params::MonteCarloParameters, mc_sweep::Function ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
    total_sweeps::Int = mc_params.therm_sweeps
    @inbounds for sweep ∈ 1:total_sweeps
        mc_sweep(model, mc_params)
    end
    return nothing
end

write_state( model, sweep_number, data_dir ) = nothing

function sweep_and_measure!( model::Model{L, H}, mc_params::MonteCarloParameters, mc_sweep::Function, data_dir ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
    total_writes = mc_params.measure_sweeps ÷ mc_params.sweeps_until_write
    @inbounds for write ∈ (1:total_writes)
        @inbounds for sweep ∈ (1:mc_params.sweeps_until_write)
            mc_sweep(model, mc_params)
        end
        write_state( model, write * mc_params.sweeps_until_write, data_dir )
    end
    return nothing
end