"""
    module Lattices

Main entry point for any simulated model to access an underlying geometrical structure.
Currently comprised of the following submodules:

- [`AbstractLattices`](@ref):
    - Interface module for all types of lattices
- [`CubicLattices`](@ref): 
    - Defines the [`AbstractCubicLattice`](@ref) subtype
    - Implements the [`CubicLattice2D`](@ref) `<: AbstractCubicLattice` and the [`CubicLattice2DParams`](@ref)  
"""
module Lattices

using Reexport
include("AbstractLattices.jl")
@reexport using .AbstractLattices
include("CubicLattices.jl")
@reexport using .CubicLattices

end # Lattices