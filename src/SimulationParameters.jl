abstract type SimulationParameters end
num_exports(simparams::SimulationParameters) = error("\nNo method defined for $(typeof(sim_params)) types.")

struct AT_2DCL_Metro_Params{T <: AbstractFloat} <: SimulationParameters
    latt_params::CubicLattice2DParams
    ham_params::AT_Parameters{T}
    mc_params::MetropolisParameters{T}
end
StructTypes.StructType(::Type{AT_2DCL_Metro_Params{T}}) where {T} = StructTypes.Struct()
num_exports(sim_params::AT_2DCL_Metro_Params) = sim_params.mc_params.total_measurements