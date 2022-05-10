using Dates

const default_mc_root = abspath(joinpath(@__DIR__,"..","..","NMR_Simulation_Data"))
const mc_date_format = "u-dd-yyyy"

eligible_models() = ["Ashkin-Teller", "RBFM"]

function build_data_directory(; model_name, data_root = default_mc_root)
    model_name ∈ eligible_models() ? nothing : error("\nModel name $(model_name) not an eligible type.\n")
    data_path = joinpath(data_root, model_name, Dates.format(Dates.now(), mc_date_format))
    return mkpath(data_path)
end

function export_parameters(data_path, params...)
    for param_struct ∈ params
        file_name = "$(Base.typename(param_struct).wrapper).json"
        println(file_name)
        # export_json(param_struct, joinpath(data_path, "$file_name.json"))
    end
    return nothing 
end