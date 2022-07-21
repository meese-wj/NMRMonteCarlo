
abstract type AbstractLattice end
construct_lattice!( latt::AbstractLattice ) = error("No implementation defined for $(typeof(latt)) types.")
num_sites( latt::AbstractLattice ) = error("No implementation defined for $(typeof(latt)) types.")
site_index( latt::AbstractLattice, indices ) = error("No implementation defined for $(typeof(latt)) and $(typeof(indices)) types.")
nearest_neighbors( latt::AbstractLattice, site ) = error("No implementation defined for $(typeof(latt)) and $(typeof(site)) types.")

include("CubicLattices.jl")