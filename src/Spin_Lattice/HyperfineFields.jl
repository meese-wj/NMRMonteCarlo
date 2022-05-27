using StaticArrays
include("../ParametersJSON.jl")
include("../Lattices/CubicLattice2D.jl")
include("../Ashkin_Teller/AT_Hamiltonian.jl")

# BaFe2As2 values in units of T/μB
# for an in-plane field      
const hyp_Aaa = 0.66
const hyp_Acc = 0.47
const hyp_Aac = 0.43
const hyp_Aab = 0.33

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

# This rotation matrix converts the spin-space coordinates
# of the As nuclear moment to the real-space coordinates
# where the crystallographic ̂b direction is the one chosen
# for the external NMR field. 
const rotation_mat = @SMatrix [ [-1 0 0] ;
                                [ 0 0 1] ;
                                [ 0 1 0] ]

const hyperfine_Amats = @SVector [ hyp_A1, hyp_A2, hyp_A3, hyp_A4 ]


# There are two As atoms per unit cell, denoted by 
# As_plus and As_minus
@enum As_atoms As_plus=1 As_minus As_enum_end
const NUM_As_ATOMS = Int(As_enum_end) - Int(As_plus)

# Define a couple of empty structs to be magnetic traits
struct Out_of_Plane end
struct Easy_Axis_In_Plane end
struct Spin_Orbit_Coupling end
const mag_vector_types = @SVector [Out_of_Plane, Easy_Axis_In_Plane, Spin_Orbit_Coupling]

# Define the mag_vector in the model with respect to 
# the expected behavior in the crystallographic basis.
# Later, use the rotation_mat to calculate the spin-
# lattice relaxation rate.
mag_vector(ty::Type{T}...) where {T} = error("\nNo method defined for mag_vectors with the $(ty) trait.")
mag_vector(::Type{Out_of_Plane}, state, site, color) = @SVector [0,0, state[site, color]]
mag_vector(::Type{Easy_Axis_In_Plane}, state, site, color) = @SVector [state[site, color], 0, 0]
function mag_vector(::Type{Spin_Orbit_Coupling}, state, site, color)
    if site_Baxter(state, site) == one(eltype(state))
        return @SVector [state[site, color], 0, 0]
    end
    return @SVector [0, state[site, color], 0]
end

function hyperfine_plus_vectors(ty, state, latt::CubicLattice2D, site )
    Atau0     = hyperfine_Amats[4] * mag_vector(ty, state, site, AT_tau)
    Atau0p1   = -hyperfine_Amats[1] * mag_vector(ty, state, site_index(latt, site, (0, 1)), AT_tau)
    Asigma0   = hyperfine_Amats[3] * mag_vector(ty, state, site, AT_sigma)
    Asigma1p0 = -hyperfine_Amats[2] * mag_vector(ty, state, site_index(latt, site, (1, 0)), AT_sigma)
    return @SVector [ Asigma0, Atau0, Asigma1p0, Atau0p1 ]
end

function hyperfine_minus_vectors(ty, state, latt::CubicLattice2D, site )
    Asigma0   = hyperfine_Amats[1] * mag_vector(ty, state, site, AT_sigma)
    Asigma0m1 = -hyperfine_Amats[4] * mag_vector(ty, state, site_index(latt, site, (0, -1)), AT_sigma)
    Atau0     = hyperfine_Amats[2] * mag_vector(ty, state, site, AT_tau)
    Atau1m0   = -hyperfine_Amats[3] * mag_vector(ty, state, site_index(latt, site, (-1, 0)), AT_tau)
    return @SVector [ Asigma0, Asigma0m1, Atau0, Atau1m0 ]
end

function hyperfine_plus(ty, state, latt::CubicLattice2D, site )
    return rotation_mat * sum( hyperfine_plus_vectors(ty, state, latt, site) )
end

function hyperfine_minus(ty, state, latt::CubicLattice2D, site )
    return rotation_mat * sum( hyperfine_minus_vectors(ty, state, latt, site ) )
end

function hyperfine_fields(ty, state, latt::CubicLattice2D, site )
    return @SVector [ hyperfine_plus(ty, state, latt, site), hyperfine_minus(ty, state, latt, site) ]
end

single_hyperfine_fluct( field ) = field[1] * field[1] + field[2] * field[2]

# By definition, in spin space, the ̂z direction points along the
# external field which doesn't couple to the As nuclear moment's 
# raising and lowering operators.
function inst_hyperfine_fluctuations(ty, state, latt::CubicLattice2D, site )
    fields = hyperfine_fields(ty, state, latt, site)
    return @SVector [ single_hyperfine_fluct(fields[1]), single_hyperfine_fluct(fields[2]) ]
end

function populate_hyperfine_fluctuations!(ty, flucts, state, latt::CubicLattice2D) where {T}
    if length(flucts) != num_sites_CL2D(latt.params) * NUM_As_ATOMS
        error("\nFluctuations state is of incorrect length: $(length(flucts)) != $(num_sites_CL2D(latt.params) * NUM_As_ATOMS)")
    end

    @inbounds for site ∈ eachindex(1:num_sites_CL2D(latt.params))
        range = @SVector [NUM_As_ATOMS * (site - 1) + 1, NUM_As_ATOMS * (site - 1) + NUM_As_ATOMS]
        flucts[ range[1]:range[2] ] .= inst_hyperfine_fluctuations(ty, state, latt, site)
    end
    return nothing
end

function populate_all_hyperfine_fluctuations!(ty, all_flucts, all_states, latt::CubicLattice2D)
    # all_flucts and all_states have reversed sizes
    # for speed on subsequent all_flucts analyses
    if size(all_flucts) != reverse(size(all_states))
        error("\nFluctuations state and MC state do not match in size: $(size(all_flucts)) != reverse($(size(all_states)))")
    end 
    @inbounds for state_idx ∈ eachindex(1:size(all_flucts)[1])
        flucts = @view all_flucts[state_idx, :]
        state =  @view all_states[:, state_idx]
        populate_hyperfine_fluctuations!(ty, flucts, state, latt)
    end
    return nothing
end

function populate_all_hyperfine_fluctuations(ty, all_states, latt::CubicLattice2D)
    all_fluctuations = zeros(reverse(size(all_states)))
    populate_all_hyperfine_fluctuations!(ty, all_fluctuations, all_states, latt)
    return all_fluctuations
end
 
function analyze_fluctuations!( all_flucts, all_states, latt::CubicLattice2D; analysis = mean )
    fluct_dists = zeros( size(all_flucts)[2], length(mag_vector_types) )
    @inbounds for (tydx, ty) ∈ enumerate(mag_vector_types)
        populate_all_hyperfine_fluctuations!(ty, all_flucts, all_states, latt)
        fluct_dists[:, tydx] .= [ analysis(all_flucts[:, idx]) for idx ∈ 1:size(all_flucts)[2] ]
    end
    return fluct_dists
end



