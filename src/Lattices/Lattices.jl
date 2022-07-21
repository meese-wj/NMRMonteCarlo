"""
    module Lattices

Main entry point for any simulated model to access an underlying geometrical structure.
Currently comprised of the following subfiles:

- [`AbstractLattices.jl`]:
    - Interface module that defines the [`AbstractLattice](@ref) as the supertype for all lattices
- [`CubicLattices`]: 
    - Defines the [`AbstractCubicLattice`](@ref) subtype
    - Implements the [`CubicLattice2D`](@ref) `<: AbstractCubicLattice` and the [`CubicLattice2DParams`](@ref)  
"""
module Lattices

include("AbstractLattices.jl")
include("CubicLattices.jl")

end # Lattices