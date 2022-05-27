using PyPlot, LaTeXStrings, StructTypes, StaticArrays

include("../../Lattices/CubicLattice2D.jl")
include("../../Ashkin_Teller/AT_Hamiltonian.jl")
include("../HyperfineFields.jl")

function pm_one(value)
    if abs(value) != one(value)
        error("\nArgument $value does not have unit magnitude.\n")
    end
    return value == one(value) ? "+1" : "-1"
end 

function state_substrings( state, latt, site )
    sigma0   = pm_one( state[ site, AT_sigma ] )
    tau0     = pm_one( state[ site, AT_tau ] )
    sigma1p0 = pm_one( state[ site_index(latt, site, (1, 0)), AT_sigma ] )
    tau0p1   = pm_one( state[ site_index(latt, site, (0, 1)), AT_tau ] )
    return (sigma0, tau0, sigma1p0, tau0p1)
end

function state_string(state, latt, site; delim = ", ")
    state_substrs = state_substrings(state, latt, site)
    result = ""
    for (idx, substr) ∈ enumerate(state_substrs)
        result *= substr
        if idx != length(state_substrs)
            result *= delim
        end
    end
    return result
end

function unit_cell_properties(; bbox_coords, fe_size, 
                                latt_space = 1, arrow_fraction = 0.45,
                                as_center = (0.5 * latt_space) .* (1, 1, 0),
                                elev = 15, azim = -80, as_size = 4 )
    properties = Dict(:bbox_coords => bbox_coords,
                        :fe_size => fe_size, 
                        :latt_space => latt_space, 
                        :as_center => as_center,
                        :arrow_fraction => arrow_fraction,
                        :elev => elev,
                        :azim => azim,
                        :as_size => as_size * fe_size)
    return properties
end

function bbox_border( bbox_coords )
    return ( [bbox_coords[1][1], bbox_coords[2][1], bbox_coords[3][1], bbox_coords[4][1], bbox_coords[1][1] ],
             [bbox_coords[1][2], bbox_coords[2][2], bbox_coords[3][2], bbox_coords[4][2], bbox_coords[1][2] ],
             [bbox_coords[1][3], bbox_coords[2][3], bbox_coords[3][3], bbox_coords[4][3], bbox_coords[1][3] ] )
end

function build_unit_cell_3d!( ax3d, properties )

    for coords ∈ properties[:bbox_coords]
        ax3d[:scatter3D](coords..., color = "red", s=properties[:fe_size], edgecolor = "black")
    end
    bbox_outline = bbox_border(properties[:bbox_coords])
    ax3d[:plot3D]( bbox_outline...,
                   color = "red", ls = "dashed", lw = 2)
    
    ax3d[:scatter3D](properties[:as_center]..., color = "blue", 
                     s=properties[:as_size],
                     edgecolor = "black", zorder=1000)
    ax3d.set_xlim( properties[:latt_space] .* (-0.15, 1.15) )
    xlims = ax3d.get_xlim()
    ax3d.set_ylim(xlims)
    ax3d.set_zlim(-properties[:latt_space], properties[:latt_space])
    ax3d.set_xlabel(L"a")
    ax3d.set_ylabel(L"b")
    ax3d.set_zlabel(L"c")
    ax3d.set_xticks([])
    ax3d.set_yticks([])
    ax3d.set_zticks([])
    ax3d.view_init(properties[:elev], properties[:azim])
    ax3d.grid(false)

    return nothing
end

function build_unit_cell_3d( properties )
    fig3d = PyPlot.figure()
    ax3d = fig3d.add_subplot(projection="3d")
    build_unit_cell_3d!(ax3d, properties)
    return (fig3d, ax3d)
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
test_state[ site_index(test_latt, test_site, (1, 0)) , AT_sigma] *= -1
test_state[ site_index(test_latt, test_site, (0, 1)) , AT_tau]   *= -1

################################################################
# Create vectors to house all Ω(±) rates for the Out_of_Plane
# and Easy_Axis_In_Plane mag_vector types. There are 16 possible
# values in total for these magnetic model types.
################################################################
test_values = [1, -1]
test_mag_types = [Out_of_Plane, Easy_Axis_In_Plane]
test_model_type = Easy_Axis_In_Plane

latt_space = 1
arrow_fraction = 0.45
fe_size = 70
bbox_coords = [ (0, 1, 0),  # σ₀
                (0, 0, 0),  # τ₀
                (1, 0, 0),  # σ₍₁,₀₎
                (1, 1, 0) ] # τ₍₀,₁₎

bbox_coords = [ latt_space .* coords for coords ∈ bbox_coords ]

ucprop = unit_cell_properties(; bbox_coords = bbox_coords,
                                fe_size = fe_size, 
                                arrow_fraction = arrow_fraction,
                                latt_space = latt_space)

fig3d = PyPlot.figure(figsize = PyPlot.figaspect(0.5))
mag_ax = fig3d.add_subplot(1,2,1, projection="3d")
hyp_ax = fig3d.add_subplot(1,2,2, projection="3d")
build_unit_cell_3d!(mag_ax, ucprop)
build_unit_cell_3d!(hyp_ax, ucprop)

spin_vectors = @SVector [ mag_vector(test_model_type, test_state, test_site, AT_sigma),
                          mag_vector(test_model_type, test_state, test_site, AT_tau),
                          -mag_vector(test_model_type, test_state, site_index(test_latt, test_site, (1, 0)), AT_sigma),
                          -mag_vector(test_model_type, test_state, site_index(test_latt, test_site, (0, 1)), AT_tau) ]
for (coords, spin) ∈ zip(ucprop[:bbox_coords], spin_vectors)
    mag_ax[:quiver](coords..., spin..., 
                    length = arrow_fraction * ucprop[:latt_space], 
                    color = "blue", pivot = "middle", zorder = 0)
end

hyp_vectors = hyperfine_plus_vectors(test_model_type, test_state, test_latt, test_site)
for arrow_prop ∈ zip(bbox_coords, hyp_vectors)
    hyp_ax[:quiver](arrow_prop[1]..., arrow_prop[2]...,
                  length = arrow_fraction * ucprop[:latt_space], 
                  color = "purple", pivot = "middle", zorder = 0)
end

net_hyp_vector = sum(hyp_vectors)
hyp_ax[:quiver](ucprop[:as_center]..., net_hyp_vector..., 
              length = arrow_fraction * ucprop[:latt_space],
              color = "orange", lw = 2, pivot = "middle", zorder = 0)


axes = @SVector [ mag_ax, hyp_ax ]
for ax ∈ axes
    ylims = ax.get_ylim()
    ax[:quiver](0.5 * ucprop[:latt_space], 0, -ucprop[:latt_space], 0, ylims[2] - ylims[1], 0, color = "black", lw = 1)
    ax.text( 0.55 * ucprop[:latt_space], 0.4 * ucprop[:latt_space], -ucprop[:latt_space], 
             L"$\mathbf{H}_0$", color = "black")
    ax.text( ucprop[:bbox_coords][1][1], ucprop[:bbox_coords][1][2] + 0.15 * ucprop[:latt_space], 0.15 * ucprop[:latt_space],  
             L"$\sigma_0$", color = "black")
    ax.text( ucprop[:bbox_coords][2][1] + 0.05 * ucprop[:latt_space], ucprop[:bbox_coords][2][2] - 0.15 * ucprop[:latt_space], -0.15 * ucprop[:latt_space],  
             L"$\tau_0$", color = "black")
    ax.text( ucprop[:bbox_coords][3][1] + 0.05 * ucprop[:latt_space], ucprop[:bbox_coords][3][2] - 0.15 * ucprop[:latt_space], -0.2 * ucprop[:latt_space],  
             L"$\sigma_{(1,0)}$", color = "black")
    ax.text( ucprop[:bbox_coords][4][1], ucprop[:bbox_coords][4][2] + 0.05 * ucprop[:latt_space], 0.15 * ucprop[:latt_space],  
             L"$\tau_{(0,1)}$", color = "black")
end

fig3d.suptitle("Ashkin-Teller State: \$\\{$( state_string(test_state, test_latt, test_site) )\\}\$")
mag_ax.set_title("Spin Configuration")
hyp_ax.set_title("Hyperfine Field Contributions")

fig3d.tight_layout()

state_str = state_string(test_state, test_latt, test_site; delim = "_" )
fig3d.savefig(joinpath(plot_dir, "$(test_model_type)_$state_str.svg"))
fig3d.savefig(joinpath(plot_dir, "$(test_model_type)_$state_str.png"))

pygui(false)
fig3d

