using JSON3, StructTypes

function import_json(fl, ::Type{T}) where {T}
    return JSON3.read( read(fl, String), T )
end

function export_json(mystruct, fl; openmode="w")
    open(fl, openmode) do io 
        JSON3.pretty(io, mystruct)
    end
    nothing
end

function pretty_display(io::IO, params)
    println(io, "\nValues for $(typeof(params)):")
    JSON3.pretty(io, params)
    println(io, "\n")
    return nothing 
end

pretty_display(params) = pretty_display(stdout, params)

function import_json_and_display(fl, ::Type{T}, io::IO = stdout; illusion=true) where {T}
    println("\nImporting values from $fl...")
    if illusion
        sleep(1.5)
    end
    params = import_json(fl, T)
    pretty_display(io, params)
    return params
end