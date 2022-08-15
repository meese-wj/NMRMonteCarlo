
import ..Lattices: CubicLattice2D
import ..Hamiltonians: BasicAshkinTellerHamiltonian, AshkinTellerParameters, 
                       energy, num_sites, spins, AT_sigma, AT_tau, num_colors,
                       site_Baxter, sigma_values, tau_values, IterateBySite

using ..SimulatingNMR

export CleanAshkinTellerModel, CleanNMRAshkinTellerModel, CleanAshkinTellerModelParameters

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
"""
struct NMRObservables{T <: Number}
    time_series_obs::Vector{TimeSeries{T}}
    acc_series_obs::Vector{AccumulatedSeries{T}}
end

"""
    CleanNMRAshkinTellerModel{T <: AbstractFloat} <: AbstractAshkinTellerModel

Basic _clean_ implementation of an Ashkin-Teller model used to simulate NMR spectra.
The [`Observables`](@ref) are comprised of a [`NMRObservables`](@ref) wrapper `struct`
of the [`BaseATM_observable_types`](@ref) as `TimeSeries` and the set of local NMR 
spectra in each unit cell (Ω⁺, Ω⁻) as `AccumulatedSeries`.
"""
struct CleanNMRAshkinTellerModel{T <: AbstractFloat} <: AbstractAshkinTellerModel
    lattice::CubicLattice2D
    hamiltonian::BasicAshkinTellerHamiltonian{T}
    observables::NMRObservables{T}

    CleanNMRAshkinTellerModel(args...) = CleanNMRAshkinTellerModel{Float64}(args...)
    function CleanNMRAshkinTellerModel{T}( Lx, Ly, Jex = 1, Kex = 0, num_meas = 0 ) where T
        latt = CubicLattice2D(Int(Lx), Int(Ly))
        atparams = AshkinTellerParameters( T(Jex), T(Kex) )
        ham = BasicAshkinTellerHamiltonian(latt, atparams)
        ts_list = TimeSeries{T}[]
        for obs ∈ BaseATMTimeSeriesObs
            push!(ts_list, TimeSeries{T}( obs, num_meas ))
        end
        acc_list = AccumulatedSeries{T}[]
        for (iter, val) ∈ enumerate(IterateBySite, ham)
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): Ω⁺", num_meas ))
            push!(acc_list, AccumulatedSeries{T}( "Site $(iter): Ω⁻", num_meas ))
        end
        return new{T}( latt, ham, NMRObservables{T}( ts_list, acc_list ) )
    end
    CleanNMRAshkinTellerModel{T}(params::CleanAshkinTellerModelParameters{T}) where T = CleanNMRAshkinTellerModel( params.Lx, params.Ly, params.Jex, params.Kex, params.num_measurements )
end

@inline TimeSeriesObservables(model::CleanNMRAshkinTellerModel) = Observables(model).time_series_obs
@inline AccumulatedSeriesObservables(model::CleanNMRAshkinTellerModel) = Observables(model).acc_series_obs

function update_NMR_values!(model::CleanNMRAshkinTellerModel, ty)
    for (site, val) ∈ enumerate(IterateBySite, Hamiltonian(model))
        Ωvals = inst_hyperfine_fluctuations(ty, Hamiltonian(model), Lattice(model), site)
        push!( AccumulatedSeriesObservables(model)[site], Ωvals[1] )
        push!( AccumulatedSeriesObservables(model)[site], Ωvals[2] )
    end
    return AccumulatedSeriesObservables(model)
end

function update_observables!(model::CleanNMRAshkinTellerModel, ty = SimulatingNMR.Easy_Axis_In_Plane)
    update_TimeSeriesObservables!(model)
    update_NMR_values!(model, ty)
    return Observables(model)
end