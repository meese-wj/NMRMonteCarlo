### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4f535963-19e0-46ce-bdcc-35c5210e0034
using DrWatson

# ╔═╡ d1e09ed5-88df-49cf-8905-f0d8dd741831
@quickactivate "NMRMonteCarlo"

# ╔═╡ a786ab69-ef38-44b9-855c-9e602d6c0034
using NMRMonteCarlo

# ╔═╡ 62018f41-adbc-4ca7-890b-7cb458e29f6d
using JLD2, UncertainHistogramming

# ╔═╡ 65fdf1a1-b9a6-4bfa-bfc6-7bc626491edb
using Measurements

# ╔═╡ cb61db54-e330-4520-a356-3efbca769929
using Plots

# ╔═╡ 16517c06-3855-11ed-1e40-e5dedf174ec8
begin
	const Lvalue = 64
	const Kvalue = 0.
	const Tc = critical_temperature(1.0, Kvalue)
	@show const βc = 1 / Tc
	const dTLow = 0.125 * Tc
	const dTHigh = 0.375 * Tc
	@show const betaHigh = 1 / (Tc - dTLow)
	@show const betaLow = 1 / (Tc + dTHigh)
	const sigma_confidence = 1
	const num_replicas = 32
end

# ╔═╡ b8fa1e9f-c9dd-41bb-aa3a-f0c0fe6db28d
const Njobs = 50

# ╔═╡ 0692581e-7534-4453-ae14-d17a41f9b3d8
beta_vals = LinRange(betaLow, betaHigh, Njobs)

# ╔═╡ 9ac73365-2cfa-41c1-89c4-401a31bc26b1
begin
	agate_hyp_χs = []
	conv_betas = Float64[]
	for idx ∈ 1:Njobs
		temp_sim = CleanNMRATMSimulation(; Lx = Lvalue, Kex = Kvalue, 
										   βvalue = beta_vals[idx], Ntherm = 2^20,
										   Nmeas=2^18, Lτ = 2^18)
		filename = savename("threaded_hyperfine_susceptibilites_Out_of_Plane", SimulationParameters(temp_sim) ) * "_#1.jld2"
		if isfile( agatedatadir(filename) )
			push!(conv_betas, beta_vals[idx])
			push!( agate_hyp_χs, JLD2.load_object( agatedatadir( filename ) ) )
		end
	end
end

# ╔═╡ 51338e72-6766-41c6-93d6-19376a7d22ec
temperatures = 1 ./ conv_betas

# ╔═╡ 02d8107b-ba22-442a-acbe-fc975c78ed81
begin
	hyp_χ_ghists = []
	for χvals ∈ agate_hyp_χs
		ghist = GaussianHistogram()
		for χi ∈ χvals
			new_χ = measurement(χi.val, χi.err / sqrt(num_replicas))
			push!(ghist, χi )
		end
		push!(hyp_χ_ghists, ghist)
	end
end

# ╔═╡ 535cd57d-1959-46db-9bfc-3032d289ab63
begin
mean_hyp_χ = similar(temperatures, Measurement{Float64})
for (idx, ghist) ∈ enumerate(hyp_χ_ghists)
	mean_hyp_χ[idx] = measurement( mean(ghist), std(ghist) )
end
end

# ╔═╡ 0a4ef262-a1e8-4c43-b290-42c17751cca7
md"""
## Lattice-Average $1/T_1$
"""

# ╔═╡ 6c7c2778-6ad5-485a-a22e-c20e73c9196b
let
	plt = plot( temperatures, mean_hyp_χ; label = nothing, 
	  	        xlabel = "Temperature \$T\$", ylabel = "\$ \\mathrm{Av}(1 / T_1) \$",
				markershape = :circle,
				title = "\$ L = $Lvalue \$",
				)
	vline!(plt, [Tc]; 
		   label = "\$T_c \\approx $(round(Tc; digits = 3)) \\, J\$", 
		   linestyle = :dash, linewidth = 2.5, z_order = :back)
	savefig(plt, plotsdir("average_W1_L-$(Lvalue)_K-$Kvalue.png"))
	plt
end

# ╔═╡ 44ba1d1d-9cde-440f-9f86-2af893f9cd5f
md"""
## Lattice-Average $1 / T_1T$
"""

# ╔═╡ 89d5a071-d24f-4fd6-b537-f6167fdd3e9b
let
	plt = plot( temperatures, mean_hyp_χ ./ temperatures; label = nothing, 
	  	        xlabel = "Temperature \$T\$", ylabel = "\$ \\mathrm{Av}(1 / T_1T) \$",
				markershape = :circle,
				title = "\$ L = $Lvalue \$")
	vline!(plt, [Tc]; 
		   label = "\$T_c \\approx $(round(Tc; digits = 3)) \\, J\$", 
		   linestyle = :dash, linewidth = 2.5, z_order = :back)
	savefig(plt, plotsdir("average_betaW1_L-$(Lvalue)_K-$Kvalue.png"))
	plt
end

# ╔═╡ d44ca312-2ee2-432d-a454-67d34c630f4b
md"""
## Standard-Deviation of $1/T_1$
"""

# ╔═╡ 997888be-71f3-4e0b-9ee0-98adbbea4ba2
let
	plt = plot( temperatures, std.(hyp_χ_ghists); label = nothing, 
	  	        xlabel = "Temperature \$T\$", ylabel = "\$ \\sigma_1 \$",
				markershape = :circle,
				title = "\$ L = $Lvalue \$")
	vline!(plt, [Tc]; 
		   label = "\$T_c \\approx $(round(Tc; digits = 3)) \\, J\$", 
		   linestyle = :dash, linewidth = 2.5, z_order = :back)
	savefig(plt, plotsdir("average_sigma1_L-$(Lvalue)_K-$Kvalue.png"))
	plt
end

# ╔═╡ 972c61fc-b942-4252-bc5c-1f927835199e
md"""
## Skewness of $1/T_1$
"""

# ╔═╡ 6b901776-ab1e-48f5-94f3-5912e3888224
let
	plt = plot( temperatures, skewness.(hyp_χ_ghists); label = nothing, 
	  	        xlabel = "Temperature \$T\$", ylabel = "\$ \\mathrm{Skewness}(1/T_1) \$",
				markershape = :circle,
				title = "\$ L = $Lvalue \$")
	vline!(plt, [Tc]; 
		   label = "\$T_c \\approx $(round(Tc; digits = 3)) \\, J\$", 
		   linestyle = :dash, linewidth = 2.5, z_order = :back)
	savefig(plt, plotsdir("average_skew1_L-$(Lvalue)_K-$Kvalue.png"))
	plt
end

# ╔═╡ 1bc79c22-e9a6-4f5b-bbbd-b7170b32dec7
md"""
## Kurtosis of $1/T_1$
"""

# ╔═╡ 445d45bc-010d-44cc-b9c7-307e766fcce8
let
	plt = plot( temperatures, kurtosis.(hyp_χ_ghists); label = nothing, 
	  	        xlabel = "Temperature \$T\$", ylabel = "\$ \\mathrm{Kurtosis}(1/T_1) \$",
				markershape = :circle,
				title = "\$ L = $Lvalue \$")
	vline!(plt, [Tc]; 
		   label = "\$T_c \\approx $(round(Tc; digits = 3)) \\, J\$", 
		   linestyle = :dash, linewidth = 2.5, z_order = :back)
	savefig(plt, plotsdir("average_kurt1_L-$(Lvalue)_K-$Kvalue.png"))
	plt
end

# ╔═╡ 5709abbc-2cfd-422a-a097-f839e0d7e5b7
md"""
## $1/T_1$ Distributions
"""

# ╔═╡ cd42356c-399b-4038-b19c-c5aaa22e9722
begin
	TempIndexString = "<input value=\"1\" type=\"range\" min=\"1\" max=\"$(length(temperatures))\"/>"
	@bind TempIndex HTML(TempIndexString)
end

# ╔═╡ 8a5aafa9-6fef-45b2-ba0e-48f799d4a41a
let
Tlabel = round(temperatures[TempIndex]; digits = 3)
Tclabel = round( temperatures[TempIndex] / Tc; digits = 3)
W1vals = LinRange(0., 1.75, 1_000)
color = temperatures[TempIndex] > Tc ? :red : :blue
dist_plt = plot( W1vals, hyp_χ_ghists[TempIndex]; label = false,
				 title = "\$ T = $(Tlabel)\\, J = $(Tclabel)\\, T_c \$",
				 xlabel = "\$ 1/T_1 \\quad \\left(L = $Lvalue, \\, N_{\\mathrm{As}} = $(2 * Lvalue^2) \\right) \$",
			     ylabel = "\$\\mathrm{p}(1/T_1)\$",
				 # ylim = (0, 350),
				 linecolor = color,
				 fillcolor = color,
)
end

# ╔═╡ Cell order:
# ╠═4f535963-19e0-46ce-bdcc-35c5210e0034
# ╠═d1e09ed5-88df-49cf-8905-f0d8dd741831
# ╠═a786ab69-ef38-44b9-855c-9e602d6c0034
# ╠═62018f41-adbc-4ca7-890b-7cb458e29f6d
# ╠═16517c06-3855-11ed-1e40-e5dedf174ec8
# ╠═b8fa1e9f-c9dd-41bb-aa3a-f0c0fe6db28d
# ╠═0692581e-7534-4453-ae14-d17a41f9b3d8
# ╠═9ac73365-2cfa-41c1-89c4-401a31bc26b1
# ╠═51338e72-6766-41c6-93d6-19376a7d22ec
# ╠═02d8107b-ba22-442a-acbe-fc975c78ed81
# ╠═65fdf1a1-b9a6-4bfa-bfc6-7bc626491edb
# ╠═535cd57d-1959-46db-9bfc-3032d289ab63
# ╠═cb61db54-e330-4520-a356-3efbca769929
# ╟─0a4ef262-a1e8-4c43-b290-42c17751cca7
# ╟─6c7c2778-6ad5-485a-a22e-c20e73c9196b
# ╟─44ba1d1d-9cde-440f-9f86-2af893f9cd5f
# ╠═89d5a071-d24f-4fd6-b537-f6167fdd3e9b
# ╟─d44ca312-2ee2-432d-a454-67d34c630f4b
# ╟─997888be-71f3-4e0b-9ee0-98adbbea4ba2
# ╟─972c61fc-b942-4252-bc5c-1f927835199e
# ╟─6b901776-ab1e-48f5-94f3-5912e3888224
# ╟─1bc79c22-e9a6-4f5b-bbbd-b7170b32dec7
# ╟─445d45bc-010d-44cc-b9c7-307e766fcce8
# ╟─5709abbc-2cfd-422a-a097-f839e0d7e5b7
# ╟─cd42356c-399b-4038-b19c-c5aaa22e9722
# ╠═8a5aafa9-6fef-45b2-ba0e-48f799d4a41a
