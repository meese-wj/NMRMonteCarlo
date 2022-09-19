### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 4f535963-19e0-46ce-bdcc-35c5210e0034
using DrWatson

# ╔═╡ d1e09ed5-88df-49cf-8905-f0d8dd741831
@quickactivate "NMRMonteCarlo"

# ╔═╡ a786ab69-ef38-44b9-855c-9e602d6c0034
using NMRMonteCarlo

# ╔═╡ 62018f41-adbc-4ca7-890b-7cb458e29f6d
using JLD2, UncertainHistogramming

# ╔═╡ cb61db54-e330-4520-a356-3efbca769929
using Plots

# ╔═╡ 16517c06-3855-11ed-1e40-e5dedf174ec8
begin
	const Lvalue = 32
	const Kvalue = 0.5
	const Tc = critical_temperature(1.0, Kvalue)
	@show const βc = 1 / Tc
	const dTLow = 0.125
	const dTHigh = 0.375
	@show const betaHigh = 1 / (Tc - dTLow)
	@show const betaLow = 1 / (Tc + dTHigh)
end

# ╔═╡ b8fa1e9f-c9dd-41bb-aa3a-f0c0fe6db28d
const Njobs = 50

# ╔═╡ 0692581e-7534-4453-ae14-d17a41f9b3d8
beta_vals = LinRange(betaLow, betaHigh, Njobs)

# ╔═╡ 51338e72-6766-41c6-93d6-19376a7d22ec
temperatures = 1 ./ beta_vals

# ╔═╡ 9ac73365-2cfa-41c1-89c4-401a31bc26b1
begin
	agate_hyp_χs = []
	for idx ∈ 1:Njobs
		temp_sim = CleanNMRATMSimulation(; Lx = Lvalue, Kex = Kvalue, 
										   βvalue = beta_vals[idx], Ntherm = 2^20,
										   Nmeas=2^18, Lτ = 2^18)
		filename = savename("hyperfine_susceptibilites_Out_of_Plane", SimulationParameters(temp_sim) ) * "_#1.jld2"
		push!( agate_hyp_χs, JLD2.load_object( agatedatadir( filename ) ) )
	end
end

# ╔═╡ 02d8107b-ba22-442a-acbe-fc975c78ed81
begin
	hyp_χ_ghists = []
	for χvals ∈ agate_hyp_χs
		ghist = GaussianHistogram()
		for χi ∈ χvals
			push!(ghist, χi)
		end
		push!(hyp_χ_ghists, ghist)
	end
end

# ╔═╡ 0a4ef262-a1e8-4c43-b290-42c17751cca7
md"""
## Lattice-Average $1/T_1$
"""

# ╔═╡ 6c7c2778-6ad5-485a-a22e-c20e73c9196b
let
	plt = plot( temperatures, mean.(hyp_χ_ghists); label = nothing, 
	  	        xlabel = "Temperature \$T\$", ylabel = "\$ \\mathrm{Av}(1 / T_1) \$",
				markershape = :circle)
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
	plt = plot( temperatures, mean.(hyp_χ_ghists) ./ temperatures; label = nothing, 
	  	        xlabel = "Temperature \$T\$", ylabel = "\$ \\mathrm{Av}(1 / T_1T) \$",
				markershape = :circle)
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
				markershape = :circle)
	vline!(plt, [Tc]; 
		   label = "\$T_c \\approx $(round(Tc; digits = 3)) \\, J\$", 
		   linestyle = :dash, linewidth = 2.5, z_order = :back)
	savefig(plt, plotsdir("average_sigma1_L-$(Lvalue)_K-$Kvalue.png"))
	plt
end

# ╔═╡ Cell order:
# ╠═4f535963-19e0-46ce-bdcc-35c5210e0034
# ╠═d1e09ed5-88df-49cf-8905-f0d8dd741831
# ╠═a786ab69-ef38-44b9-855c-9e602d6c0034
# ╠═62018f41-adbc-4ca7-890b-7cb458e29f6d
# ╠═16517c06-3855-11ed-1e40-e5dedf174ec8
# ╠═b8fa1e9f-c9dd-41bb-aa3a-f0c0fe6db28d
# ╠═0692581e-7534-4453-ae14-d17a41f9b3d8
# ╠═51338e72-6766-41c6-93d6-19376a7d22ec
# ╠═9ac73365-2cfa-41c1-89c4-401a31bc26b1
# ╠═02d8107b-ba22-442a-acbe-fc975c78ed81
# ╠═cb61db54-e330-4520-a356-3efbca769929
# ╟─0a4ef262-a1e8-4c43-b290-42c17751cca7
# ╠═6c7c2778-6ad5-485a-a22e-c20e73c9196b
# ╟─44ba1d1d-9cde-440f-9f86-2af893f9cd5f
# ╠═89d5a071-d24f-4fd6-b537-f6167fdd3e9b
# ╟─d44ca312-2ee2-432d-a454-67d34c630f4b
# ╠═997888be-71f3-4e0b-9ee0-98adbbea4ba2
