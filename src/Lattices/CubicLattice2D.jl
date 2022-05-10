abstract type AbstractLattice end
construct_lattice!( latt::AbstractLattice ) = error("No implementation defined for $(typeof(latt)) types.")
num_sites( latt::AbstractLattice ) = error("No implementation defined for $(typeof(latt)) types.")
site_index( latt::AbstractLattice, indices ) = error("No implementation defined for $(typeof(latt)) and $(typeof(indices)) types.")
nearest_neighbors( latt::AbstractLattice, site ) = error("No implementation defined for $(typeof(latt)) and $(typeof(site)) types.")

const NN_SQUARE_LATT = 4
struct CubicLattice2DParams
    Lx::Int
    Ly::Int
end
StructTypes.StructType(::Type{CubicLattice2DParams}) = StructTypes.Struct()
corresponding_object(::Type{CubicLattice2DParams}) = CubicLattice2D

struct CubicLattice2D <: AbstractLattice
    params::CubicLattice2DParams
    neighbors::Matrix{Int32}
end
num_sites_CL2D(Lx, Ly) = Lx * Ly
num_sites_CL2D(lattparams) = (lattparams.Lx * lattparams.Ly)::Int
num_sites( latt::CubicLattice2D ) = num_sites_CL2D( latt.params.Lx, latt.params.Ly )
Base.getindex( latt::CubicLattice2D, site, neighbor ) = latt.neighbors[ site, neighbor ]
Base.setindex!( latt::CubicLattice2D, value, site, neighbor ) = latt.neighbors[ site, neighbor ] = value 
site_index( latt::CubicLattice2D, indices::Tuple{Int, Int} ) = latt.params.Lx * ( indices[2] - 1 ) + indices[1]
function pbc_add(x1, x2, Lsize) 
    output = x1 + x2
    while output > Lsize
        output -= Lsize
    end
    while output < 1
        output += Lsize
    end
    return output
end

function construct_lattice!( latt::CubicLattice2D )
    for ydx ∈ 1:latt.params.Ly, xdx ∈ 1:latt.params.Lx
        site = site_index(latt, (xdx, ydx))
        
        # Order of neighbors: (x, y-1), (x-1, y), (x+1, y), (x, y+1)
        latt[site, 1] = site_index(latt, ( xdx, pbc_add(ydx, -1, latt.params.Ly) ) )
        latt[site, 2] = site_index(latt, ( pbc_add(xdx, -1, latt.params.Lx), ydx ) )
        latt[site, 3] = site_index(latt, ( pbc_add(xdx, 1, latt.params.Lx), ydx ) )
        latt[site, 4] = site_index(latt, ( xdx, pbc_add(ydx, 1, latt.params.Ly) ) )
    end
    return nothing
end

function CubicLattice2D(; params_file::String = joinpath(@__DIR__,"default_CL2D_Params.jl"),
                          display_io::IO = stdout )
    CL2D_params = import_json_and_display( params_file, CubicLattice2DParams, display_io )
    temp = CubicLattice2D( CL2D_params, Matrix(undef, num_sites_CL2D(CL2D_params.Lx, CL2D_params.Ly), NN_SQUARE_LATT ) )
    construct_lattice!( temp )
    return temp
end

function CubicLattice2D( lattparams::CubicLattice2DParams )
    temp = CubicLattice2D( lattparams, Matrix(undef, num_sites_CL2D(lattparams), NN_SQUARE_LATT ) )
    construct_lattice!( temp )
    return temp 
end 

CubicLattice2D(Lx::Int, Ly::Int) = CubicLattice2D( CubicLattice2DParams(Lx, Ly) )

nearest_neighbors(latt::CubicLattice2D, site) = view(latt.neighbors, site, :)