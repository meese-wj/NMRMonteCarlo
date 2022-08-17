using Parameters2JSON

# Automatically submodulizes this code
import ..Hamiltonians: IterateByDoFType, Hamiltonian, DoF_energy_change, accept_move!
export MetropolisParameters, metropolis_sweep!, thermalize!, sweep_and_measure!

@jsonable struct MetropolisParameters{T <: AbstractFloat} <: AbstractMonteCarloParameters
    βvalues::Vector{T}
    therm_sweeps::Int
    measure_sweeps::Int
    total_measurements::Int
end
thermalization_sweeps(params::MetropolisParameters) = params.therm_sweeps
sampling_sweeps(params::MetropolisParameters) = params.measure_sweeps
total_measurements(params::MetropolisParameters) = params.total_measurements

function MetropolisParameters{T}(; params_file = joinpath(@__DIR__, "default_Metropolis_Params.jl"),
                                   display_io::IO = stdout) where {T}
    return import_json_and_display(params_file, MetropolisParameters{T}, display_io)
end

@inline metropolis_accepted(ΔE, β) = ( ΔE < zero(ΔE) || rand() < @fastmath exp(-β * ΔE) )::Bool

@inline function metropolis_update!( model::AbstractModel, beta, site )
    ΔE = DoF_energy_change( Hamiltonian(model), Lattice(model), site )
    accept_move!( metropolis_accepted( ΔE, beta ), Hamiltonian(model), site )
    return nothing
end

# function metropolis_sweep!( model::Model{L, H}, beta ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
@inline function metropolis_sweep!( model::AbstractModel, beta )
    for (iteration, dof_site_val) ∈ enumerate( IterateByDoFType, Hamiltonian(model) )
        metropolis_update!(model, beta, dof_site_val[1])
    end
    return nothing
end

# function thermalize!( model::Model{L, H}, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
function thermalize!( model::AbstractModel, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function )
    total_sweeps = thermalization_sweeps(mc_params)
    @inbounds for sweep ∈ 1:total_sweeps
        mc_sweep(model, beta)
    end
    return nothing
end

# function sweep_and_measure!( model::Model{L, H}, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function, state_container::AbstractArray = [] ) where {L <: AbstractLattice, H <: AbstractHamiltonian}
function sweep_and_measure!( model::AbstractModel, beta, mc_params::AbstractMonteCarloParameters, mc_sweep::Function )
    total_exports = total_measurements(mc_params)
    spe = sweeps_per_export(mc_params)
    @inbounds for write ∈ (1:total_exports)
        @inbounds for sweep ∈ (1:spe)
            mc_sweep(model, beta)
        end
        update_observables!(model)
        # Only write out observables if the state container is defined
        # if length(state_container) > 0
        #     # export_state!( state_container, model.ham, write )
        #     export_state!( state_container, Hamiltonian(model), write )
        # end
    end
    return nothing
end