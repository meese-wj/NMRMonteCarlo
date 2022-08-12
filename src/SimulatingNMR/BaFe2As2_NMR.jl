
using StaticArrays
# Automatically submodulizes this code
import ..Lattices: CubicLattice2D
using ..Hamiltonians


const hyp_Aaa = 0.66
const hyp_Acc = 0.47
const hyp_Aac = 0.43
const hyp_Aab = 0.33

"""
    const hyp_Aaa

As hyperfine field components of BaFe2As2 in units of T/μB
for an in-plane field. Here ``(a, b, c)`` are the axes of 
the crystallographic 1-Fe unit cell.

# Unique Components:
* [`hyp_Aaa`](@ref) = $(hyp_Aaa)
* [`hyp_Acc`](@ref) = $(hyp_Acc)
* [`hyp_Aac`](@ref) = $(hyp_Aac)
* [`hyp_Aab`](@ref) = $(hyp_Aab)
"""
hyp_Aaa, hyp_Acc, hyp_Aac, hyp_Aab

const hyp_A1 = @SMatrix [ [hyp_Aaa hyp_Aab hyp_Aac] ;
                          [hyp_Aab hyp_Aaa hyp_Aac] ;
                          [hyp_Aac hyp_Aac hyp_Acc] ]

const hyp_A2 = @SMatrix [ [hyp_Aaa -hyp_Aab -hyp_Aac] ;
                          [-hyp_Aab hyp_Aaa  hyp_Aac] ;
                          [-hyp_Aac hyp_Aac  hyp_Acc] ]
                          
const hyp_A3 = @SMatrix [ [hyp_Aaa  -hyp_Aab  hyp_Aac] ;
                          [-hyp_Aab  hyp_Aaa -hyp_Aac] ;
                          [hyp_Aac  -hyp_Aac  hyp_Acc] ]
                          
const hyp_A4 = @SMatrix [ [ hyp_Aaa  hyp_Aab -hyp_Aac] ;
                          [ hyp_Aab  hyp_Aaa -hyp_Aac] ;
                          [-hyp_Aac -hyp_Aac  hyp_Acc] ]

"""
    const rotation_mat

This rotation matrix converts the spin-space coordinates
of the As nuclear moment to the real-space coordinates
where the crystallographic ̂b direction is the one chosen
for the external NMR field.
""" 
const rotation_mat = @SMatrix [ [-1 0 0] ;
                                [ 0 0 1] ;
                                [ 0 1 0] ]

"""
    const hyperfine_Amats

A vector collecting of the hyperfine tensors as

1. `hyp_A1` = $(hyp_A1)
1. `hyp_A2` = $(hyp_A2)
1. `hyp_A3` = $(hyp_A3)
1. `hyp_A4` = $(hyp_A4)

in that specific order.
"""
const hyperfine_Amats = @SVector [ hyp_A1, hyp_A2, hyp_A3, hyp_A4 ]

"""
    @enum As_atoms

There are two As atoms per unit cell, denoted by 
`As_plus` and `As_minus`.
"""
@enum As_atoms As_plus=1 As_minus As_enum_end
# TODO: Destroy this.
const NUM_As_ATOMS = Int(As_enum_end) - Int(As_plus)

"""
    num_atoms(::Type{As_atoms}) → 2

Return the number of atoms relevant for As NMR.
"""
@inline num_atoms(::Type{As_atoms}) = Int(As_enum_end) - Int(As_plus)

"""
    abstract type AbstractNMRConstruct end

Define a supertype for all different model `Type` constructs for NMR modeling.

See also [`mag_vector`](@ref) for specific implementations.
"""
abstract type AbstractNMRConstruct end
"""
    struct Out_of_Plane <: AbstractNMRConstruct

Construct `Type` denoting the Fe spin moments point out of the basal plane.

See also [`mag_vector`](@ref) for specific implementations.
"""
struct Out_of_Plane <: AbstractNMRConstruct end
"""
    struct Easy_Axis_In_Plane <: AbstractNMRConstruct

Construct `Type` denoting the Fe spin moments point along one crystallographic
axis within the basal plane.

See also [`mag_vector`](@ref) for specific implementations.
"""
struct Easy_Axis_In_Plane <: AbstractNMRConstruct end
"""
    struct Spin_Orbit_Coupling <: AbstractNMRConstruct

Construct `Type` denoting that the Fe spins flip which crystallographic
axis within the basal plane they point along depending on the local 
nematic.

See also [`mag_vector`](@ref) for specific implementations.
"""
struct Spin_Orbit_Coupling <: AbstractNMRConstruct end
const mag_vector_types = @SVector [Out_of_Plane, Easy_Axis_In_Plane, Spin_Orbit_Coupling]

"""
    mag_vector(::Type{<: AbstractNMRConstruct}, args...)

Build the magnetic spin vector for a given [`AbstractNMRConstruct`] with respect to 
its expected behavior in the crystallographic basis. Later, one uses the 
[`rotation_mat`](@ref) to calculate the proxy spin-lattice relaxation rate.
"""
mag_vector(ty::Type{T}, args...) where {T <: AbstractNMRConstruct} = throw(MethodError(mag_vector, ty, args...))
"""
    mag_vector(::Type{Out_of_Plane}, state, site, color)

Model the magnetic spins to point out of the basal plane, so along the 
crystallographic ̂c axis.
"""
mag_vector(::Type{Out_of_Plane}, state, site, color) = @SVector [0,0, state[site, color]]
"""
    mag_vector(::Type{Easy_Axis_In_Plane}, state, site, color)

Model the magnetic spins to point in the basal plane along the 
crystallographic ̂a axis.
"""
mag_vector(::Type{Easy_Axis_In_Plane}, state, site, color) = @SVector [state[site, color], 0, 0]
"""
    mag_vector(::Type{Spin_Orbit_Coupling}, state, site, color)

Model the magnetic spins to point in the basal plane along either the 
crystallographic ̂a or ̂b axes, depending on the local nematic.
"""
function mag_vector(::Type{Spin_Orbit_Coupling}, state, site, color)
    if site_Baxter(state, site) == one(eltype(state))
        return @SVector [state[site, color], 0, 0]
    end
    return @SVector [0, state[site, color], 0]
end

"""
    hyperfine_plus_vectors(ty, state, ::CubicLattice2D, site)

Return the set of four hyperfine field vectors for the `As_plus` atom of the [`As_atoms`](@ref) `Enum`.

    # Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `state::AbstractArray`: the field of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_plus_vectors(ty, state, latt::CubicLattice2D, site )
    Atau0     = hyperfine_Amats[4] * mag_vector(ty, state, site, AT_tau)
    Atau0p1   = -hyperfine_Amats[1] * mag_vector(ty, state, site_index(latt, site, (0, 1)), AT_tau)
    Asigma0   = hyperfine_Amats[3] * mag_vector(ty, state, site, AT_sigma)
    Asigma1p0 = -hyperfine_Amats[2] * mag_vector(ty, state, site_index(latt, site, (1, 0)), AT_sigma)
    return @SVector [ Asigma0, Atau0, Asigma1p0, Atau0p1 ]
end
"""
    hyperfine_plus_vectors(ty, state, ::CubicLattice2D, site)

Return the set of four hyperfine field vectors for the `As_plus` atom of the [`As_atoms`](@ref) `Enum`.

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `state::AbstractArray`: the field of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_minus_vectors(ty, state, latt::CubicLattice2D, site )
    Asigma0   = hyperfine_Amats[1] * mag_vector(ty, state, site, AT_sigma)
    Asigma0m1 = -hyperfine_Amats[4] * mag_vector(ty, state, site_index(latt, site, (0, -1)), AT_sigma)
    Atau0     = hyperfine_Amats[2] * mag_vector(ty, state, site, AT_tau)
    Atau1m0   = -hyperfine_Amats[3] * mag_vector(ty, state, site_index(latt, site, (-1, 0)), AT_tau)
    return @SVector [ Asigma0, Asigma0m1, Atau0, Atau1m0 ]
end

"""
    hyperfine_plus(ty, state, ::CubicLattice2D, site)

Sum the set of hyperfine fields and then rotate them from the 
crystallographic axes into the spin-space basis.

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `state::AbstractArray`: the field of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_plus(ty, state, latt::CubicLattice2D, site )
    return rotation_mat * sum( hyperfine_plus_vectors(ty, state, latt, site) )
end

"""
    hyperfine_minus(ty, state, ::CubicLattice2D, site)

Sum the set of hyperfine fields and then rotate them from the 
crystallographic axes into the spin-space basis.

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `state::AbstractArray`: the field of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_minus(ty, state, latt::CubicLattice2D, site )
    return rotation_mat * sum( hyperfine_minus_vectors(ty, state, latt, site ) )
end

"""
    hyperfine_fields(ty, state, ::CubicLattice2D, site)

Return the set of [`hyperfine_plus`](@ref) and [`hyperfine_minus`](@ref) fields acting on
the As atoms within a unit cell. 

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `state::AbstractArray`: the field of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_fields(ty, state, latt::CubicLattice2D, site )
    return @SVector [ hyperfine_plus(ty, state, latt, site), hyperfine_minus(ty, state, latt, site) ]
end

@doc raw"""
    single_hyperfine_fluct( field )

Calculate the instantaneous `field` fluctuations that contribute to the NMR
spin-lattice relaxation rate. For example, if the field is ``\vec{h}`` **in 
spin-space**, then the instantaneous proxy spin-lattice relaxation rate is

```math
\Omega = h_x^2 + h_y^2.
```
"""
single_hyperfine_fluct( field ) = field[1] * field[1] + field[2] * field[2]