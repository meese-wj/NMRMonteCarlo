using PyPlot, LaTeXStrings, StructTypes, StaticArrays

include("../../Lattices/CubicLattice2D.jl")
include("../../Ashkin_Teller/AT_Hamiltonian.jl")
include("../HyperfineFields.jl")

function figure_3d()
    fig = figure()
    ax = fig[:gca](projection = "3d")
    return fig, ax
end

function pm_one(value)
    if abs(value) != one(value)
        error("\nArgument $value does not have unit magnitude.\n")
    end
    return value == one(value) ? "+1" : "-1"
end 

function state_string( state, latt, site )
    sigma0   = pm_one( state[ site, AT_sigma ] )
    tau0     = pm_one( state[ site, AT_tau ] )
    sigma1p0 = pm_one( state[ site_index(latt, site, (1, 0)), AT_sigma ] )
    tau0p1   = pm_one( state[ site_index(latt, site, (0, 1)), AT_tau ] )
    return "$sigma0, $tau0, $sigma1p0, $tau0p1"
end

plot_dir = joinpath(@__DIR__, "TestData")
mkpath(plot_dir)

################################################################
# Make a 3 × 3 grid and measure everything with respect to the 
# central site located at (2, 2) → 5
################################################################
test_latt_params = CubicLattice2DParams( 3, 3 )
test_latt = CubicLattice2D( test_latt_params )
test_site = site_index(test_latt, (2, 2))

################################################################
# Create a set of Ashkin-Teller spins 
################################################################
test_state = ones( 2 * num_sites(test_latt) )
test_state[test_site, AT_tau] *= -1
# test_state[ site_index(test_latt, test_site, (1, 0)) , AT_sigma] *= -1
test_state[ site_index(test_latt, test_site, (0, 1)) , AT_tau]   *= -1

################################################################
# Create vectors to house all Ω(±) rates for the Out_of_Plane
# and Easy_Axis_In_Plane mag_vector types. There are 16 possible
# values in total for these magnetic model types.
################################################################
test_values = [1, -1]
test_mag_types = [Out_of_Plane, Easy_Axis_In_Plane]

latt_space = 1
arrow_fraction = 0.45
fe_size = 40
bbox_coords = [ (0, 1, 0),  # σ₀
                (0, 0, 0),  # τ₀
                (1, 0, 0),  # σ₍₁,₀₎
                (1, 1, 0) ] # τ₍₀,₁₎

bbox_coords = [ latt_space .* coords for coords ∈ bbox_coords ]

as_center = 0.5 .* (latt_space, latt_space, 0)

fig3d, ax3d = figure_3d()
for coords ∈ bbox_coords
    ax3d[:scatter3D](coords..., color = "red", s=fe_size, edgecolor = "black")
end
ax3d[:plot3D]( [bbox_coords[1][1], bbox_coords[2][1], bbox_coords[3][1], bbox_coords[4][1], bbox_coords[1][1] ],
               [bbox_coords[1][2], bbox_coords[2][2], bbox_coords[3][2], bbox_coords[4][2], bbox_coords[1][2] ],
               [bbox_coords[1][3], bbox_coords[2][3], bbox_coords[3][3], bbox_coords[4][3], bbox_coords[1][3] ],
               color = "red", ls = "dashed", lw = 2)

ax3d[:scatter3D](as_center..., color = "blue", s=6*fe_size, edgecolor = "black", zorder=1000)

hyp_vectors = hyperfine_plus_vectors(Easy_Axis_In_Plane, test_state, test_latt, test_site)
for arrow_prop ∈ zip(bbox_coords, hyp_vectors)
    ax3d[:quiver](arrow_prop[1]..., arrow_prop[2]...,
                  length = arrow_fraction * latt_space, 
                  color = "black", pivot = "middle", zorder = 0)
end

net_hyp_vector = sum(hyp_vectors)
ax3d[:quiver](as_center..., net_hyp_vector..., 
              length = arrow_fraction * latt_space,
              color = "orange", lw = 2, pivot = "middle", zorder = 0)

ylims = ax3d.get_ylim()
ax3d[:quiver](0.5 * latt_space, 0, -latt_space, 0, ylims[2] - ylims[1], 0, color = "purple", lw = 1)
ax3d.text( 0.55 * latt_space, 0.4 * latt_space, -latt_space, 
          L"$\mathbf{H}_0$", color = "purple")

ax3d.set_xlim( latt_space .* (-0.15, 1.15) )
ax3d.set_ylim(ylims)
ax3d.set_zlim(-latt_space, latt_space)
ax3d.set_xlabel(L"a")
ax3d.set_ylabel(L"b")
ax3d.set_zlabel(L"c")
ax3d.set_xticks([])
ax3d.set_yticks([])
ax3d.set_zticks([])
ax3d.view_init(15, -80)
ax3d.set_title("\$\\{$( state_string(test_state, test_latt, test_site) )\\}\$")
ax3d.grid(false)

pygui(false)
fig3d

