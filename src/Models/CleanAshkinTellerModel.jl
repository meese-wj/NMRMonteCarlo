
import Base: iterate, length
import ..Lattices: CubicLattice2D
import ..Hamiltonians: BasicAshkinTellerHamiltonian, AshkinTellerParameters, 
                       energy, num_sites, spins, AT_sigma, AT_tau, num_colors,
                       site_Baxter, sigma_values, tau_values, IterateBySite

using ..SimulatingNMR
using MonteCarloMeasurementUncertainty

export 
# Base overloads
       iterate, length,
# Concrete type exports
       CleanAshkinTellerModel, CleanNMRAshkinTellerModel, CleanAshkinTellerModelParameters

struct CleanAshkinTellerModelParameters{T <: AbstractFloat}
    Lx::Int
    Ly::Int
    Jex::T
    Kex::T
    num_measurements::Int
    
    CleanAshkinTellerModelParameters{T}(Lx, Ly, Jex, Kex, n) where T = new{T}(Int(Lx), Int(Ly), T(Jex), T(Kex), Int(n))
    CleanAshkinTellerModelParameters(Lx, Ly, Jex, Kex, n) = CleanAshkinTellerModelParameters{Float64}(Lx, Ly, Jex, Kex, n)
end

"""
    CleanAshkinTellerModel{T <: AbstractFloat} <: AbstractAshkinTellerModel

Basic _clean_ implementation of an Ashkin-Teller model. The [`Observables`](@ref)
are only `TimeSeries` observables defined by [`BaseATM_observable_types`](@ref).

!!! note 
    `observables` iteration is implemented by default since it is a `Vector`.
"""
struct CleanAshkinTellerModel{T <: AbstractFloat} <: AbstractAshkinTellerModel
    lattice::CubicLattice2D
    hamiltonian::BasicAshkinTellerHamiltonian{T}
    observables::Vector{TimeSeries{T}}

    CleanAshkinTellerModel(args...) = CleanAshkinTellerModel{Float64}(args...)
    function CleanAshkinTellerModel{T}( Lx, Ly, Jex = 1, Kex = 0, num_meas = 0 ) where T
        latt = CubicLattice2D(Int(Lx), Int(Ly))
        atparams = AshkinTellerParameters( T(Jex), T(Kex) )
        ham = BasicAshkinTellerHamiltonian(latt, atparams)
        obs_list = TimeSeries{T}[]
        for obs ∈ BaseATMTimeSeriesObs
            push!(obs_list, TimeSeries{T}( obs, num_meas ))
        end
        return new{T}( latt, ham, obs_list )
    end
    CleanAshkinTellerModel{T}(params::CleanAshkinTellerModelParameters{T}) where T = CleanAshkinTellerModel( params.Lx, params.Ly, params.Jex, params.Kex, params.num_measurements )
end

@inline TimeSeriesObservables(model::CleanAshkinTellerModel) = Observables(model)
@inline update_observables!(model::CleanAshkinTellerModel) = update_TimeSeriesObservables!(model)

"""
    NMRObservables{T <: Number}

Wrapper `struct` around two `Vector`s of `MonteCarloMeasurement`s:

* `time_series_obs::Vector{TimeSeries{T}}`
* `acc_series_obs::Vector{AccumulatedSeries{T}}`

!!! note 
    The `Base` `iterate` interface is explicitly defined.
"""
struct NMRObservables{T <: Number}
    time_series_obs::Vector{TimeSeries{T}}
    acc_series_obs::Vector{AccumulatedSeries{T}}
end
"""
    length(::NMRObservables)

Return the total number of [`NMRObservables`](@ref) contained.
"""
length(all_obs::NMRObservables) = length(all_obs.time_series_obs) + length(all_obs.acc_series_obs)
"""
    iterate(::NMRObservables, [state = 1])

`Base` iteration interface method to traverse the [`NMRObservables`](@ref).
"""
function iterate(all_obs::NMRObservables, state = 1)
    # TODO: Maybe there is a better way to do this with Iterators.flatten?
    index = one(state) + (state - one(state)) % length(all_obs.time_series_obs)
    return state <= length(all_obs) ? 
          ( state <= length(all_obs.time_series_obs) ? all_obs.time_series_obs[index] : all_obs.acc_series_obs[index], state + 1 ) : 
          nothing
end

"""
    CleanNMRAshkinTellerModel{T <: AbstractFloat} <: AbstractAshkinTellerModel

Basic _clean_ implementation of an Ashkin-Teller model used to simulate NMR spectra.
The [`Observables`](@ref) are comprised of a [`NMRObservables`](@ref) wrapper `struct`
of the [`BaseATM_observable_types`](@ref) as `TimeSeries` and the set of local NMR 
spectra in each unit cell (Ω⁺, Ω⁻) as `AccumulatedSeries`.
"""
struct CleanNMRAshkinTellerModel{T <: AbstractFloat} <: AbstractAshkinTellerModel
    nmr_spin_type::Type
    lattice::CubicLattice2D
    hamiltonian::BasicAshkinTellerHamiltonian{T}
    observables::NMRObservables{T}

    CleanNMRAshkinTellerModel(args...) = CleanNMRAshkinTellerModel{Float64}(args...)
    function CleanNMRAshkinTellerModel{T}( Lx, Ly, Jex = 1, Kex = 0, num_meas = 0, spin_type = SimulatingNMR.Easy_Axis_In_Plane ) where T
        latt = CubicLattice2D(Int(Lx), Int(Ly))
        atparams = AshkinTellerParameters( T(Jex), T(Kex) )
        ham = BasicAshkinTellerHamiltonian(latt, atparams)
        ts_list = TimeSeries{T}[]
        for obs ∈ BaseATMTimeSeriesObs
            push!(ts_list, TimeSeries{T}( obs, num_meas ))
        end
        acc_list = AccumulatedSeries{T}[]
        for (iter, val) ∈ enumerate(IterateBySite, ham)
            # push!(acc_list, AccumulatedSeries{T}( "Site $(iter): Ω⁺", num_meas ))
            # push!(acc_list, AccumulatedSeries{T}( "Site $(iter): Ω⁻", num_meas ))
            
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): |hx⁺|", num_meas ))
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): |hy⁺|", num_meas ))
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): (hx⁺)²", num_meas ))
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): (hy⁺)²", num_meas ))
            
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): |hx⁻|", num_meas ))
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): |hy⁻|", num_meas ))
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): (hx⁻)²", num_meas ))
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): (hy⁻)²", num_meas ))
        end
        return new{T}( spin_type, latt, ham, NMRObservables{T}( ts_list, acc_list ) )
    end
    CleanNMRAshkinTellerModel{T}(params::CleanAshkinTellerModelParameters{T}, spin_type = SimulatingNMR.Easy_Axis_In_Plane) where T = CleanNMRAshkinTellerModel( params.Lx, params.Ly, params.Jex, params.Kex, params.num_measurements, spin_type )
end

@inline TimeSeriesObservables(model::CleanNMRAshkinTellerModel) = Observables(model).time_series_obs
@inline AccumulatedSeriesObservables(model::CleanNMRAshkinTellerModel) = Observables(model).acc_series_obs

@inline nmr_observable_index(atom_idx, ob) = NMR_OBS_PER_AS * (atom_idx - one(NMR_OBS_PER_AS)) + ob
@inline nmr_observable_index(site, As_sign::As_atoms, ob) = NMR_OBS_PER_AS * (As_atom_index(site, As_sign) - one(NMR_OBS_PER_AS)) + ob

function update_NMR_values!(model::CleanNMRAshkinTellerModel, ty)
    for (site, val) ∈ enumerate(IterateBySite, Hamiltonian(model))
        hyp_obs_tup = inst_hyperfine_observables(ty, Hamiltonian(model), Lattice(model), site)
        # Ωvals = inst_hyperfine_fluctuations(ty, Hamiltonian(model), Lattice(model), site)
        # push!( AccumulatedSeriesObservables(model)[As_atom_index(site, As_plus)],  Ωvals[As_plus] )
        # push!( AccumulatedSeriesObservables(model)[As_atom_index(site, As_minus)], Ωvals[As_minus] )
        for (hyp_obs, As_sign) ∈ zip(hyp_obs_tup, (As_plus, As_minus)), ob ∈ (1:4)  
            push!( AccumulatedSeriesObservables(model)[nmr_observable_index(site, As_sign, ob)],  hyp_obs[ob] )
        end
    end
    return AccumulatedSeriesObservables(model)
end

function update_observables!(model::CleanNMRAshkinTellerModel, ty = model.nmr_spin_type)
    update_TimeSeriesObservables!(model)
    update_NMR_values!(model, ty)
    return Observables(model)
end

function collect_hyperfine_susceptibilites(model::AbstractAshkinTellerModel)
    hyp_chi_vals = []
    num_chi = length( AccumulatedSeriesObservables(model) ) ÷ NUM_OBS_PER_AS
    for atom_idx ∈ (1:num_chi)
        hyp_tuple = AccumulatedSeriesObservables(model)[ nmr_observable_index(atom_idx, 1) : nmr_observable_index(atom_idx, NUM_OBS_PER_AS) ]
        measurements = []
        for obs ∈ hyp_tuple
            push!(measurements, measurement(obs))
        end
        push!(hyp_chi_vals, hyperfine_field_susceptibility(measurements))
    end
    return hyp_chi_vals
end