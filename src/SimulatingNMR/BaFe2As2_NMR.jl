
import Base: to_index
using StaticArrays
# Automatically submodulizes this code
import ..Lattices: CubicLattice2D, site_index
using ..Hamiltonians
import ..Hamiltonians: site_Baxter

export inst_hyperfine_fluctuations, hyperfine_field_parts_to_save, inst_hyperfine_observables, NMR_OBS_PER_AS, hyperfine_field_susceptibility,
       to_index, As_atom_index, As_atoms, As_plus, As_minus, Easy_Axis_In_Plane, Out_of_Plane, Spin_Orbit_Coupling, Form_Factor_Test

const hyp_Aaa = 0.66
const hyp_Acc = 0.47
const hyp_Aac = 0.43
const hyp_Aab = 0.33

const FeAs_diagram = raw"""
    b
Fe2 : Fe1
    As ---> a
Fe4   Fe3
"""

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
    const hyperfine_Amats

A vector collecting of the hyperfine tensors as

1. `hyp_A1` = $(hyp_A1)
1. `hyp_A2` = $(hyp_A2)
1. `hyp_A3` = $(hyp_A3)
1. `hyp_A4` = $(hyp_A4)

in that specific order. The numbering corresponds to 
the following convention:

$(FeAs_diagram)

given by [Katagawa _et al._ J. Phys. Soc. Jpn. (2008)](https://journals.jps.jp/doi/10.1143/JPSJ.77.114709).
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

@inline to_index(atom::As_atoms) = Int(atom)
@inline As_atom_index(site_idx, atom::As_atoms) = num_atoms(As_atoms) * (site_idx - one(site_idx)) + to_index(atom)

"""
    abstract type AbstractNMRConstruct end

Define a supertype for all different model `Type` constructs for NMR modeling.

See also [`mag_vector`](@ref) for specific implementations.
"""
abstract type AbstractNMRConstruct end
"""
    struct Out_of_Plane <: AbstractNMRConstruct

Construct `Type` denoting the Fe spin moments point out of the basal plane
with the external field applied in the same direction.

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
"""
    struct Form_Factor_Test <: AbstractNMRConstruct

Artificial test where all hyperfine tensors are the identity
matrix to see if they diverge at the transition. In this case,
the ̂c axis is the -̂y axis in spin space, and the external field
is directed along the ̂b axis.
"""
struct Form_Factor_Test <: AbstractNMRConstruct end
const mag_vector_types = @SVector [Out_of_Plane, Easy_Axis_In_Plane, Spin_Orbit_Coupling, Form_Factor_Test]

"""
    spin_space(::Type{<: AbstractNMRConstruct}, component::Type{Val{<: Char}}, hvec)

Define the default spin-space `component`s for a given [`AbstractNMRConstruct`](@ref)
relative to the crystallographic axes of the 1-Fe unit cell.

The default will be taken with the external field applied along the ̂b direction. This 
defines the spin-space ̂z to be along the ̂b direction. Defining ̂a == ̂x, then a 90° rotation
clockwise about the ̂x axis yields the right transformation. This makes the following true:
̂a == ̂x, ̂b == ̂z, ̂c == -̂y.
"""
spin_space(::Type{T}, ::Type{Val{'x'}}, hvec) where T <: AbstractNMRConstruct = hvec[1] 
spin_space(::Type{T}, ::Type{Val{'y'}}, hvec) where T <: AbstractNMRConstruct = hvec[3] 
spin_space(::Type{T}, ::Type{Val{'z'}}, hvec) where T <: AbstractNMRConstruct = -hvec[2] 

"""
    spin_space(::Type{Out_of_Plane}, component::Type{Val{<: Char}}, hvec)

Overwrite the default the spin-space `component`s for the [`Out_of_Plane`](@ref)
NMR construct relative to the crystallographic axes of the 1-Fe unit cell.

In this case, the external field applied along the ̂c direction. This defines the spin-space
z to be along the ̂c direction. Thus, the crystallographic and spin-space axes are coincident.
"""
spin_space(::Type{Out_of_Plane}, ::Type{Val{'x'}}, hvec) = hvec[1] 
spin_space(::Type{Out_of_Plane}, ::Type{Val{'y'}}, hvec) = hvec[2] 
spin_space(::Type{Out_of_Plane}, ::Type{Val{'z'}}, hvec) = hvec[3] 

"""
    mag_vector(::Type{<: AbstractNMRConstruct}, args...)

Build the magnetic spin vector for a given [`AbstractNMRConstruct`] with respect to 
its expected behavior in the crystallographic basis. One must be careful about the 
meaning of the ``x`` and ``y`` components of the internal field in [`spin_space`](@ref)
in calculating the spin-lattice relaxation rate, for example in [`single_hyperfine_fluct`](@ref).
"""
@inline mag_vector(ty::Type{T}, args...) where {T <: AbstractNMRConstruct} = throw(MethodError(mag_vector, ty, args...))
"""
    mag_vector(::Type{Out_of_Plane}, ham, site, color)

Model the magnetic spins to point out of the basal plane, so along the 
crystallographic ̂c axis.
"""
@inline mag_vector(::Type{Out_of_Plane}, ham, site, color) = SVector{3, eltype(ham)}(zero(eltype(ham)), zero(eltype(ham)), ham[site, color])
"""
    mag_vector(::Type{Easy_Axis_In_Plane}, ham, site, color)

Model the magnetic spins to point in the basal plane along the 
crystallographic ̂a axis.
"""
@inline mag_vector(::Type{Easy_Axis_In_Plane}, ham, site, color) = SVector{3, eltype(ham)}(ham[site, color], zero(eltype(ham)), zero(eltype(ham)))
"""
    mag_vector(::Type{Spin_Orbit_Coupling}, ham, site, color)

Model the magnetic spins to point in the basal plane along either the 
crystallographic ̂a or ̂b axes, depending on the local nematic.
"""
function mag_vector(::Type{Spin_Orbit_Coupling}, ham, site, color)
    if site_Baxter(ham, site) == one(eltype(ham))
        return SVector{3, eltype(ham)}(ham[site, color], zero(eltype(ham)), zero(eltype(ham)))
    end
    return SVector{3, eltype(ham)}(zero(eltype(ham)), ham[site, color], zero(eltype(ham)))
end

@inline mag_vector(::Type{Form_Factor_Test}, ham, site, color) = ( S = ham[site, color]; SVector{3, eltype(ham)}(S, S, S) )


"""
    hyperfine_plus_vectors(ty, ham, ::CubicLattice2D, site)

Return the set of four hyperfine field vectors for the `As_plus` atom of the [`As_atoms`](@ref) `Enum`.

    # Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `ham`: the Hamiltonian of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
@inline function hyperfine_plus_vectors(ty, ham, latt::CubicLattice2D, site )
    Atau0     = hyperfine_Amats[4] * mag_vector(ty, ham, site, AT_tau)
    Atau0p1   = -hyperfine_Amats[1] * mag_vector(ty, ham, site_index(latt, site, (0, 1)), AT_tau)
    Asigma0   = hyperfine_Amats[2] * mag_vector(ty, ham, site, AT_sigma)
    Asigma1p0 = -hyperfine_Amats[3] * mag_vector(ty, ham, site_index(latt, site, (1, 0)), AT_sigma)
    return @SVector [ Asigma0, Atau0, Asigma1p0, Atau0p1 ]
end

"""
    hyperfine_minus_vectors(ty, ham, ::CubicLattice2D, site)

Return the set of four hyperfine field vectors for the `As_plus` atom of the [`As_atoms`](@ref) `Enum`.

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `ham`: the Hamiltonian of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
@inline function hyperfine_minus_vectors(ty, ham, latt::CubicLattice2D, site )
    Asigma0   = hyperfine_Amats[1] * mag_vector(ty, ham, site, AT_sigma)
    Asigma0m1 = -hyperfine_Amats[4] * mag_vector(ty, ham, site_index(latt, site, (0, -1)), AT_sigma)
    Atau0     = hyperfine_Amats[3] * mag_vector(ty, ham, site, AT_tau)
    Atau1m0   = -hyperfine_Amats[2] * mag_vector(ty, ham, site_index(latt, site, (-1, 0)), AT_tau)
    return @SVector [ Asigma0, Asigma0m1, Atau0, Atau1m0 ]
end

"""
    hyperfine_plus(ty, ham, ::CubicLattice2D, site)

Sum the set of hyperfine fields in real-space.

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `ham`: the Hamiltonian of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_plus(ty, ham, latt::CubicLattice2D, site )
    return sum( hyperfine_plus_vectors(ty, ham, latt, site) )
end

function hyperfine_plus_spins(::Type{Out_of_Plane}, ham, latt::CubicLattice2D, site)
    S1 = -ham[site_index(latt, site, (0, 1)), AT_tau]
    S2 =  ham[site, AT_sigma]
    S3 = -ham[site_index(latt, site, (1, 0)), AT_sigma]
    S4 =  ham[site, AT_tau]
    return S1, S2, S3, S4
end

hyperfine_plus_spins(::Type{Form_Factor_Test}, ham, latt::CubicLattice2D, site) = hyperfine_plus_spins(Out_of_Plane, ham, latt, site)

function hyperfine_spin_vector(::Type{Out_of_Plane}, S1, S2, S3, S4)
    return ( hyp_Aac * (S1 - S2 + S3 - S4),
             hyp_Aac * (S1 + S2 - S3 - S4),
             hyp_Acc * (S1 + S2 + S3 + S4) )
end

@inline hyperfine_spin_vector(::Type{Form_Factor_Test}, S1, S2, S3, S4) = ( sum = S1 + S2 + S3 + S4; ( sum, sum, sum ) )

function hyperfine_plus(::Type{Out_of_Plane}, ham, latt::CubicLattice2D, site)
    return hyperfine_spin_vector(Out_of_Plane, hyperfine_plus_spins(Out_of_Plane, ham, latt, site)... )
end

function hyperfine_plus(::Type{Form_Factor_Test}, ham, latt::CubicLattice2D, site)
    return hyperfine_spin_vector(Form_Factor_Test, hyperfine_plus_spins(Form_Factor_Test, ham, latt, site)... )
end

"""
    hyperfine_minus(ty, ham, ::CubicLattice2D, site)

Sum the set of hyperfine fields in real-space.

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `ham`: the Hamiltonian of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_minus(ty, ham, latt::CubicLattice2D, site )
    return sum( hyperfine_minus_vectors(ty, ham, latt, site ) )
end

function hyperfine_minus_spins(::Type{Out_of_Plane}, ham, latt::CubicLattice2D, site)
    S1 =  ham[site, AT_sigma]
    S2 = -ham[site_index(latt, site, (-1, 0)), AT_tau]
    S3 =  ham[site, AT_tau]
    S4 = -ham[site_index(latt, site, (0, -1)), AT_sigma]
    return S1, S2, S3, S4
end

hyperfine_minus_spins(::Type{Form_Factor_Test}, ham, latt::CubicLattice2D, site) = hyperfine_minus_spins(Out_of_Plane, ham, latt, site)

function hyperfine_minus(::Type{Out_of_Plane}, ham, latt::CubicLattice2D, site)
    return hyperfine_spin_vector(Out_of_Plane, hyperfine_minus_spins(Out_of_Plane, ham, latt, site)... )
end

function hyperfine_minus(::Type{Form_Factor_Test}, ham, latt::CubicLattice2D, site)
    return hyperfine_spin_vector(Form_Factor_Test, hyperfine_minus_spins(Form_Factor_Test, ham, latt, site)... )
end

"""
    hyperfine_fields(ty, ham, ::CubicLattice2D, site)

Return the set of [`hyperfine_plus`](@ref) and [`hyperfine_minus`](@ref) fields acting on
the As atoms within a unit cell. 

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `ham`: the Hamiltonian of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function hyperfine_fields(ty, ham, latt::CubicLattice2D, site )
    return ( hyperfine_plus(ty, ham, latt, site), hyperfine_minus(ty, ham, latt, site) )
end

const NMR_OBS_PER_AS = 4
"""
    hyperfine_field_parts_to_save(ty, hyp)

Return the set of observables to save from a specific hyperfine field.
The specific order of the observables in the return Tuple are

```math
(|h_x|, |h_y|, h_x^2, h_y^2)
```

where ``x`` and ``y`` are defined in [`spin_space`](@ref).

# Arguments:

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `hyp::SVector{3}`: the hyperfine field in question
"""
function hyperfine_field_parts_to_save(ty, hyp)
    hx, hy = spin_space(ty, Val{'x'}, hyp), spin_space(ty, Val{'y'}, hyp)
    return abs(hx), abs(hy), hx * hx, hy * hy
end
@doc raw"""
    hyperfine_field_susceptibility(hyperfine_tuple)

Compute the susceptibility for the hyperfine field that contributes to 
the spin-lattice relaxation rate. Note that this quantity is a *proxy*
in that it represents proportionality as

```math
\frac{1}{T_1} \propto \sum_{\alpha = x,y} \chi_{\alpha\alpha}^{(h)} \equiv \sum_{\alpha = x,y} \langle h_\alpha^2 \rangle - \langle |h_\alpha |\rangle^2.
```
"""
hyperfine_field_susceptibility(hyperfine_tuple) = ( (abshx, abshy, hx2, hy2) = hyperfine_tuple; (hx2 - abshx * abshx) + (hy2 - abshy * abshy) ) 

"""
    inst_hyperfine_observables(ty, ham, ::CubicLattice2D, site)

Return the pair of hyperfine observables in a single 
unit cell. 

# Arguments

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `ham`: the Hamiltonian of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
@generated function inst_hyperfine_observables(ty, ham, latt::CubicLattice2D, site )
    expr = quote
        fields = hyperfine_fields(ty, ham, latt, site)
        return ( hyperfine_field_parts_to_save(ty, fields[1]), hyperfine_field_parts_to_save(ty, fields[2]) )
    end
    return expr
end
# function inst_hyperfine_observables(ty, ham, latt::CubicLattice2D, site )
#     fields = hyperfine_fields(ty, ham, latt, site)
#     return ( hyperfine_field_parts_to_save(ty, fields[1]), hyperfine_field_parts_to_save(ty, fields[2]) )
# end

@doc raw"""
    single_hyperfine_fluct(::Type{<: AbstractNMRConstruct}, field)

Calculate the instantaneous `field` fluctuations that contribute to the NMR
spin-lattice relaxation rate. For example, if the `field` is ``\vec{h}`` **in 
spin-space**, then the instantaneous proxy spin-lattice relaxation rate is

```math
\Omega = h_x^2 + h_y^2.
```

!!! note 
    By definition, in **spin-space**, the ̂z direction points along the
    external field which doesn't couple to the As nuclear moment's 
    raising and lowering operators. The definition of the ̂z direction
    in spin-space follows from specific [`AbstractNMRConstruct`](@ref)s.
    The default implementation is `x == 1` and `y == 3`. Specific 
    implementations are required for non-default constructs.
"""
function single_hyperfine_fluct(ty::Type{<: AbstractNMRConstruct}, field)
    hx, hy = (spin_space(ty, Val{'x'}, field), spin_space(ty, Val{'y'}, field))
    return hx * hx + hy * hy
end
function single_hyperfine_fluct(::Type{Out_of_Plane}, field) 
    hx, hy = (spin_space(Out_of_Plane, Val{'x'}, field), spin_space(Out_of_Plane, Val{'y'}, field))
    return hx * hx + hy * hy
end
"""
    inst_hyperfine_fluctuations(ty, ham, ::CubicLattice2D, site)

Return the pair of proxy spin-lattice relaxation rates in a single 
unit cell. 

# Arguments

1. `ty::Type{<: AbstractNMRConstruct}`: defines which magnetic spin construct to implement
1. `ham`: the Hamiltonian of spins
1. `::CubicLattice2D`: needed to find the [`nearest_neighbors`](@ref)
1. `site`: which unit cell to operate on
"""
function inst_hyperfine_fluctuations(ty, ham, latt::CubicLattice2D, site )
    fields = hyperfine_fields(ty, ham, latt, site)
    return @SVector [ single_hyperfine_fluct(ty, fields[1]), single_hyperfine_fluct(ty, fields[2]) ]
end