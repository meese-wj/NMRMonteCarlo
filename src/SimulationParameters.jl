using Parameters2JSON

abstract type SimulationParameters end
num_exports(simparams::SimulationParameters) = error("\nNo method defined for $(typeof(simparams)) types.")

@jsonable struct AT_2DCL_Metro_Params{T <: AbstractFloat} <: SimulationParameters
    latt_params::CubicLattice2DParams
    ham_params::AT_Parameters{T}
    mc_params::MetropolisParameters{T}
end
num_exports(sim_params::AT_2DCL_Metro_Params) = sim_params.mc_params.total_measurements