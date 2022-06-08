using Dates

const default_mc_root = abspath(joinpath(@__DIR__,"..","..","NMR_Simulation_Data"))
const mc_date_format = "u-dd-yyyy"

eligible_models() = ["Ashkin-Teller", "RBFM"]

function concat_ids(identifiers...)
    if length(identifiers) == 0
        error("\nIdentifiers must be provided to save the simulation files.\n")
    end

    ids = ""
    for id in identifiers
        ids *= "_$id"
    end
    return ids
end

simulation_date_format( date, form = mc_date_format ) = Dates.format(date, form)

simulation_date_directory( model_name, date, data_root = default_mc_root, form = mc_date_format ) = joinpath(data_root, model_name, simulation_date_format(date, form))

function build_data_directory(; model_name, data_root = default_mc_root)
    model_name ∈ eligible_models() ? nothing : error("\nModel name $(model_name) not an eligible type.\n")
    data_path = simulation_date_directory(model_name, Dates.now(), data_root, mc_date_format)
    return mkpath(data_path)
end

function export_parameters(data_path, simparams::SimulationParameters, identifiers...)
    ids = concat_ids(identifiers...)
    file_name = "$(typeof(simparams))$ids.params"
    export_json(simparams, joinpath(data_path, file_name))
    return nothing 
end

function export_simulation(data_path, states::AbstractArray, 
                           simparams::SimulationParameters, identifiers...)
    # Assumes data_path exists...
    export_states(states, data_path, identifiers...)
    # export_parameters(data_path, simparams, identifiers...)
    return nothing
end

function import_parameters(data_path, simparam_t::Type, identifiers...)
    param_file = joinpath(data_path, "$(simparam_t)$(concat_ids(identifiers...)).params")
    isfile(param_file) ? nothing : error("\nSimulation parameters file:\n\t$param_file\nnot found.")
    return import_json( param_file, simparam_t )::SimulationParameters
end

function import_simulation(data_path, sim_params, identifiers...)
    # Get the parameters
    # sim_params = import_parameters(data, simparam_t, identifiers...)

    # Get the states with a temporary Hamiltonian
    temp_latt = (reciprocal_type(sim_params.latt_params))(sim_params.latt_params)
    temp_ham = (reciprocal_type(sim_params.ham_params))(temp_latt, sim_params.ham_params)
    states_file = joinpath(data_path, "states$(concat_ids(identifiers...)).bin")
    states = import_states( reciprocal_type(sim_params.ham_params), states_file, num_DoF(temp_ham), num_exports(sim_params) )
    # return (sim_params, states)
    return states
end