using Dates

const default_mc_root = abspath(joinpath(@__DIR__,"..","..","NMR_Simulation_Data"))
const mc_date_format = "u-dd-yyyy"

function build_data_directory(; model_name, data_root = default_mc_root)
    data_path = joinpath(data_root, model_name, Dates.format(Dates.now(), mc_date_format))
    return mkpath(data_path)
end