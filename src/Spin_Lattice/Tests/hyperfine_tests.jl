using Plots, LaTeXStrings, StructTypes

include("../../Lattices/CubicLattice2D.jl")
include("../../Ashkin_Teller/AT_Hamiltonian.jl")
include("../HyperfineFields.jl")

plot_dir = joinpath(@__DIR__, "TestData")
mkpath(plot_dir)

function bar_graph_axes( data_set )
    unique_values = unique(data_set)
    unique_counts = similar(unique_values)
    for (idx, val) ∈ enumerate(unique_values)
        unique_counts[idx] = count(x -> x == val, data_set)
    end
    return (unique_values, unique_counts)
end

function bar_graph!(plot, data_set; kwargs... )
    vals, counts = bar_graph_axes(data_set)
    return bar!(plot, vals, counts; 
                xticks = ( vals, broadcast(x -> "$x", vals) ),
                kwargs...)
end

function bar_graph( data_set; kwargs... )
    plt = Plots.plot()
    bar_graph!(plt, data_set; kwargs...)
    return plt
end

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

################################################################
# Create vectors to house all Ω(±) rates for the Out_of_Plane
# and Easy_Axis_In_Plane mag_vector types. There are 16 possible
# values in total for these magnetic model types.
################################################################
test_values = [1, -1]
test_mag_types = [Out_of_Plane, Easy_Axis_In_Plane]

for test_type ∈ test_mag_types
    test_rate_plus  = zeros(Float64, 0)
    test_rate_minus = zeros(Float64, 0) 

    # Testing the Ω⁺ rates
    for val1 ∈ test_values
        test_state[ site_index(test_latt, test_site, (0, 0)), AT_sigma ] = val1
        for val2 ∈ test_values
            test_state[ site_index(test_latt, test_site, (0, 0)), AT_tau ] = val2
            for val3 ∈ test_values
                test_state[ site_index(test_latt, test_site, (1, 0)), AT_sigma ] = val3 
                for val4 ∈ test_values
                    test_state[ site_index(test_latt, test_site, (0, 1)), AT_tau ] = val4
                    Ωvals = inst_hyperfine_fluctuations(test_type, test_state, test_latt, test_site )
                    push!(test_rate_plus, Ωvals[1])
                    println("{$val1, $val2, $val3, $val4} → $(round(Ωvals[1], digits=4))")
                end
            end
        end
    end

    # Testing the Ω⁻ rates
    for val1 ∈ test_values
        test_state[ site_index(test_latt, test_site, (0, 0)), AT_sigma ] = val1
        for val2 ∈ test_values
            test_state[ site_index(test_latt, test_site, (0, 0)), AT_tau ] = val2
            for val3 ∈ test_values
                test_state[ site_index(test_latt, test_site, (-1, 0)), AT_tau ] = val3 
                for val4 ∈ test_values
                    test_state[ site_index(test_latt, test_site, (0, -1)), AT_sigma ] = val4
                    Ωvals = inst_hyperfine_fluctuations(test_type, test_state, test_latt, test_site )
                    push!(test_rate_minus, Ωvals[2])
                end
            end
        end
    end

    test_rate_plus = round.(test_rate_plus, digits = 4)
    test_rate_minus = round.(test_rate_minus, digits = 4)
    unique(test_rate_plus)
    unique(test_rate_minus)

    @show Symbol(test_type)
    @show mean(test_rate_plus)
    @show mean(test_rate_minus)
    println()

    plt_Ωp = bar_graph(test_rate_plus; 
                    label="$(length(test_rate_plus)) States",
                    ylabel = "Degeneracy",
                    xlabel = L"$\Omega^{(+)}$ $(\mathrm{T}^2/\mu_B^2)$")
                    

    plt_Ωm = bar_graph(test_rate_minus;
                    label="$(length(test_rate_plus)) States",
                    #    ylabel = "Degeneracy",
                    xlabel = L"$\Omega^{(-)}$ $(\mathrm{T}^2/\mu_B^2)$")

    test_plot = plot(plt_Ωp, plt_Ωm; 
                    layout = (1, 2), link=:both,
                    size = (600,500),
                    xrotation = 55, 
                    plot_title = replace(String(Symbol(test_type)), "_" => " "))

    savefig(test_plot, joinpath(plot_dir, "$(String(Symbol(test_type)))_values.svg"))
    test_plot
end

################################################################
# Now create a similar set of plots for the Spin_Orbit_Coupling
# magnetic model. With this model, there are 2⁶ = 64 different
# possible spin states. 
################################################################
soc_test_values_plus  = zeros(Float64, 0)
soc_test_values_minus = zeros(Float64, 0)

for val1 ∈ test_values
    test_state[ site_index(test_latt, test_site, (0, 0)), AT_sigma ] = val1
    for val2 ∈ test_values
        test_state[ site_index(test_latt, test_site, (0, 0)), AT_tau ] = val2
        for val3 ∈ test_values
            test_state[ site_index(test_latt, test_site, (1, 0)), AT_sigma ] = val3 
            for val4 ∈ test_values
                test_state[ site_index(test_latt, test_site, (1, 0)), AT_tau ] = val4 
                for val5 ∈ test_values
                    test_state[ site_index(test_latt, test_site, (0, 1)), AT_sigma ] = val5 
                    for val6 ∈ test_values
                        test_state[ site_index(test_latt, test_site, (0, 1)), AT_tau ] = val6
                        Ωvals = inst_hyperfine_fluctuations(Spin_Orbit_Coupling, test_state, test_latt, test_site )
                        push!(soc_test_values_plus, Ωvals[1])
                    end
                end
            end
        end
    end
end

for val1 ∈ test_values
    test_state[ site_index(test_latt, test_site, (0, 0)), AT_sigma ] = val1
    for val2 ∈ test_values
        test_state[ site_index(test_latt, test_site, (0, 0)), AT_tau ] = val2
        for val3 ∈ test_values
            test_state[ site_index(test_latt, test_site, (0, -1)), AT_sigma ] = val3 
            for val4 ∈ test_values
                test_state[ site_index(test_latt, test_site, (0, -1)), AT_tau ] = val4 
                for val5 ∈ test_values
                    test_state[ site_index(test_latt, test_site, (-1, 0)), AT_sigma ] = val5 
                    for val6 ∈ test_values
                        test_state[ site_index(test_latt, test_site, (-1, 0)), AT_tau ] = val6
                        Ωvals = inst_hyperfine_fluctuations(Spin_Orbit_Coupling, test_state, test_latt, test_site )
                        push!(soc_test_values_minus, Ωvals[2])
                    end
                end
            end
        end
    end
end

soc_test_values_plus = broadcast(x -> round(x, digits = 4), soc_test_values_plus)
soc_test_values_minus = broadcast(x -> round(x, digits = 4), soc_test_values_minus)

@show Symbol(Spin_Orbit_Coupling)
@show mean(soc_test_values_plus)
@show mean(soc_test_values_minus)
println()

soc_plt_p = bar_graph(soc_test_values_plus; 
                      label = "$(length(soc_test_values_plus)) States",
                      ylabel = "Degeneracy",
                      xlabel = L"$\Omega^{(+)}$ $(\mathrm{T}^2/\mu_B^2)$")


soc_plt_m = bar_graph(soc_test_values_minus; 
                      label = "$(length(soc_test_values_plus)) States",
                      xlabel = L"$\Omega^{(-)}$ $(\mathrm{T}^2/\mu_B^2)$")
                      


soc_plot = plot(soc_plt_p, soc_plt_m;
                link = :both,
                xrotation = 55,
                size = (600,500),
                layout = (1,2),
                plot_title = replace(String(Symbol(Spin_Orbit_Coupling)), "_" => " ")
                )

savefig(soc_plot, joinpath(plot_dir, "$(String(Symbol(Spin_Orbit_Coupling)))_values.svg"))
soc_plot
