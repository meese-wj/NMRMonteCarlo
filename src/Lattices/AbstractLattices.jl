
import DrWatson

export 
# DrWatson overloads
       savename,
# AbstractLattice interface functions
       construct_lattice!, num_sites, site_index, nearest_neighbors, parameters

"""
    abstract type AbstractLattice end

The `supertype` for all `Type`s of lattices one can dream of.

# Required Interface Methods

One must define the following `methods` for each new (non-`abstract`) `AbstractLattice`:

- [`construct_lattice!`](@ref)
- [`num_sites`](@ref)
- [`site_index`](@ref)
- [`nearest_neighbors`](@ref)
- [`parameters`](@ref)

# Default Interface Methods

The following represent the default `methods` for `AbstractLattice`s. One may still overload
them for any peculiar `subtype`.

- [`DrWatson.savename`](@ref)
"""
abstract type AbstractLattice end

# Required interface methods with MethodError defaults
"""
    construct_lattice!(::AbstractLattice)

Build a lattice and its geometry. Meant to be used in a constructor.
"""
construct_lattice!( latt::AbstractLattice ) = throw(MethodError(construct_lattice!, latt))
"""
    num_sites(::AbstractLattice)

Return the number of sites that an [`AbstractLattice`](@ref) contains.

# Example

```jldoctest
julia> latt = CubicLattice2D(4, 4);

julia> num_sites(latt)
16
```
"""
num_sites( latt::AbstractLattice ) = throw(MethodError(num_sites, latt))
"""
    site_index(::AbstractLattice, indices::itr)

Calculate the flattened index from an iterable set of `indices` in an [`AbstractLattice`](@ref).

# Example

```jldoctest
julia> latt = CubicLattice2D(4, 4);

julia> site_index(latt, (1, 2))
6
```
"""
site_index( latt::AbstractLattice, indices ) = throw(MethodError(site_index, latt, indices))
"""
    nearest_neighbors(::AbstractLattice, site)

Return the set of `nearest_neighbors` for a given `site` in the [`AbstractLattice`](@ref).

# Example

```jldoctest
julia> latt = CubicLattice2D(4, 4);

julia> nearest_neighbors(latt, 1)
4-element view(::Matrix{Int32}, 1, :) with eltype Int32:
 13
  4
  2
  5
```
"""
nearest_neighbors( latt::AbstractLattice, site ) = throw(MethodError(nearest_neighbors, latt, site))
"""
    parameters(::AbstractLattice)

Return the `parameters` used to define the [`AbstractLattice`](@ref).

# Example

```jldoctest
julia> latt = CubicLattice2D(4, 4);

julia> parameters(latt, 1)
CubicLattice2DParams(4, 4)
```
"""
parameters(latt::AbstractLattice) = throw(MethodError(parameters, latt))

# Default interface methods
"""
    DrWatson.savename(::AbstractLattice)

Overload `savename` for a given [`AbstractLattice`](@ref). This assumes that the
[`parameters`](@ref) for that [`AbstractLattice`](@ref) instance return a simple 
enough object for `savename` to act upon.

# Example

```jldoctest
julia> latt = CubicLattice2D(4, 4);

julia> savename(latt)
"Lx=4_Ly=4"
```
"""
DrWatson.savename(latt::AbstractLattice) = DrWatson.savename( parameters(latt) )