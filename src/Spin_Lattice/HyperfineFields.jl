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

const hyp_constants = @SVector [hyp_Aaa, hyp_Acc, hyp_Aac, hyp_Aab]
const max_hyp_const = maximum(hyp_constants)

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

const hyperfine_Amats = @SVector [ hyp_A1, hyp_A2, hyp_A3, hyp_A4 ]

# There are two As atoms per unit cell, denoted by 
# As_plus and As_minus
@enum As_atoms As_plus=1 As_minus As_enum_end
const NUM_As_ATOMS = Int(As_enum_end) - Int(As_plus)

struct OutofPlane end
struct SpinOrbitCoupling end
const mag_vector_types = @SVector [OutofPlane, SpinOrbitCoupling]

mag_vector(ty::Type{T}...) where {T} = error("\nNo method defined for mag_vectors with the $(ty) trait.")
mag_vector(::Type{OutofPlane}, state, site, color) = @SVector [0,0, state[site, color]]
function mag_vector(::Type{SpinOrbitCoupling}, state, site, color)
    if site_Baxter(state, site) == one(eltype(state))
        return @SVector [state[site, color], 0, 0]
    end
    return @SVector [0, state[site, color], 0]
end

function hyperfine_plus(ty, state, latt::CubicLattice2D, site )
    hyp1 = hyperfine_Amats[4] * mag_vector(ty, state, site, AT_tau)
    hyp2 = hyperfine_Amats[1] * mag_vector(ty, state, site_index(latt, site, (0, 1)), AT_tau)
    hyp3 = hyperfine_Amats[3] * mag_vector(ty, state, site, AT_sigma)
    hyp4 = hyperfine_Amats[2] * mag_vector(ty, state, site_index(latt, site, (1, 0)), AT_sigma)
    return hyp1 - hyp2 + hyp3 - hyp4
end

function hyperfine_minus(ty, state, latt::CubicLattice2D, site )
    hyp1 = hyperfine_Amats[1] * mag_vector(ty, state, site, AT_sigma)
    hyp2 = hyperfine_Amats[4] * mag_vector(ty, state, site_index(latt, site, (0, -1)), AT_sigma)
    hyp3 = hyperfine_Amats[2] * mag_vector(ty, state, site, AT_tau)
    hyp4 = hyperfine_Amats[3] * mag_vector(ty, state, site_index(latt, site, (-1, 0)), AT_tau)
    return hyp1 - hyp2 + hyp3 - hyp4
end

function hyperfine_fields(ty, state, latt::CubicLattice2D, site )
    return @SVector [ hyperfine_plus(ty, state, latt, site), hyperfine_minus(ty, state, latt, site) ]
end

function inst_hyperfine_fluctuations(ty, state, latt::CubicLattice2D, site )
    fields = hyperfine_fields(ty, state, latt, site)
    return @SVector [fields[1][1] * fields[1][1] + fields[1][3] * fields[1][3],
                     fields[2][1] * fields[2][1] + fields[2][3] * fields[2][3] ]
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

function analyze_fluctuations!( all_flucts, all_states, latt::CubicLattice2D; analysis = mean )
    fluct_dists = zeros( size(all_flucts)[2], length(mag_vector_types) )
    @inbounds for (tydx, ty) ∈ enumerate(mag_vector_types)
        populate_all_hyperfine_fluctuations!(ty, all_flucts, all_states, latt)
        fluct_dists[:, tydx] .= [ analysis(all_flucts[:, idx]) for idx ∈ 1:size(all_flucts)[2] ]
    end
    return fluct_dists
end



