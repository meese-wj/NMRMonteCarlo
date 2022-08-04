
using Parameters2JSON

export CubicLattice2D, CubicLattice2DParams

abstract type AbstractCubicLattice <: AbstractLattice end
parameters(latt::AbstractCubicLattice) = latt.params

const NN_SQUARE_LATT = 4
@jsonable struct CubicLattice2DParams
    Lx::Int
    Ly::Int
end
reciprocal_type(::Type{CubicLattice2DParams}) = CubicLattice2D

struct CubicLattice2D <: AbstractCubicLattice
    params::CubicLattice2DParams
    neighbors::Matrix{Int32}
end
@inline num_sites_CL2D(Lx, Ly) = Lx * Ly
@inline num_sites_CL2D(lattparams) = num_sites_CL2D(lattparams.Lx, lattparams.Ly)
@inline num_sites( latt::CubicLattice2D ) = num_sites_CL2D( latt.params )
Base.getindex( latt::CubicLattice2D, site, neighbor ) = latt.neighbors[ site, neighbor ]
Base.setindex!( latt::CubicLattice2D, value, site, neighbor ) = latt.neighbors[ site, neighbor ] = value 
@inline site_index( latt::CubicLattice2D, indices::Tuple{Int, Int} ) = latt.params.Lx * ( indices[2] - 1 ) + indices[1]
@inline indices_from_site( latt::CubicLattice2D, origin_site::Int ) = ( 1 + (origin_site - 1) % latt.params.Lx, 1 + (origin_site - 1) ÷ latt.params.Lx )
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

site_index(latt::CubicLattice2D, origin::Tuple{Int,Int}, Δr::Tuple{Int, Int}) = site_index(latt, 
                                                                                          ( pbc_add(origin[1], Δr[1], latt.params.Lx),
                                                                                            pbc_add(origin[2], Δr[2], latt.params.Ly) ) )

function site_index(latt::CubicLattice2D, origin_site::Int, Δr::Tuple{Int, Int} )
    origin = indices_from_site(latt, origin_site)
    return site_index(latt, origin, Δr)
end

function construct_lattice!( latt::CubicLattice2D )
    for ydx ∈ 1:latt.params.Ly, xdx ∈ 1:latt.params.Lx
        site = site_index(latt, (xdx, ydx))
        
        # Order of neighbors: (x, y-1), (x-1, y), (x+1, y), (x, y+1)
        latt[site, 1] = site_index(latt, (xdx, ydx), ( 0, -1))
        latt[site, 2] = site_index(latt, (xdx, ydx), (-1,  0))
        latt[site, 3] = site_index(latt, (xdx, ydx), ( 1,  0))
        latt[site, 4] = site_index(latt, (xdx, ydx), ( 0,  1))
        # latt[site, 2] = site_index(latt, ( pbc_add(xdx, -1, latt.params.Lx), ydx ) )
        # latt[site, 3] = site_index(latt, ( pbc_add(xdx, 1, latt.params.Lx), ydx ) )
        # latt[site, 4] = site_index(latt, ( xdx, pbc_add(ydx, 1, latt.params.Ly) ) )
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
