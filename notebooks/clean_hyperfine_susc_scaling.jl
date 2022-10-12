### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ b6739750-48d2-11ed-12ad-5bd8f89b0ab6
using DrWatson

# ╔═╡ abc6dc6b-f864-447e-bbf1-40daa2d92594
@quickactivate "NMRMonteCarlo"

# ╔═╡ d9f59e54-ef96-4fdc-8236-f1fa15a08185
using NMRMonteCarlo

# ╔═╡ 4ea5386f-bf4e-4fa2-8023-70d4f9064b76
using JLD2

# ╔═╡ c4c38477-3ab0-4a6c-9a76-dd8f7a712650
using UncertainHistogramming

# ╔═╡ 873e6f73-feba-4c2d-9134-0f9065ca0af9
using Plots

# ╔═╡ 445f6190-b1c7-4f49-be35-6eae76d184c6
begin
	const Lvalues = [32, 64, 128]
	const Kvalue = 0.0
	const Tc = critical_temperature(1.0, Kvalue)
	@show const βc = 1 / Tc
	const dTLow = 0.125 * Tc
	const dTHigh = 0.375 * Tc
	@show const betaHigh = 1 / (Tc - dTLow)
	@show const betaLow = 1 / (Tc + dTHigh)
	const sigma_confidence = 1
	const num_replicas = 32
end

# ╔═╡ 8bdd7dbc-2f91-4db8-9fe8-681ee3f831a7
const Njobs = 50

# ╔═╡ 5833ae1c-d641-409a-bf5f-1b291ca3ed41
beta_vals = LinRange(betaLow, betaHigh, Njobs)

# ╔═╡ 61f806fd-6f17-47fa-8b50-6db4179a00ad
begin
	agate_hyp_χs = Dict()
	conv_betas = Dict()
	for Lvalue ∈ Lvalues
		agate_hyp_χs_L = []
		conv_betas_L = []
		for idx ∈ 1:Njobs
			temp_sim = CleanNMRATMSimulation(; Lx = Lvalue, Kex = Kvalue, 
											   βvalue = beta_vals[idx], Ntherm = 2^20,
											   Nmeas=2^18, Lτ = 2^18)
			filename = savename("threaded_hyperfine_susceptibilites_Out_of_Plane", SimulationParameters(temp_sim) ) * ".jld2"
			if isfile( agatedatadir(filename) )
				push!(conv_betas_L, beta_vals[idx])
				push!( agate_hyp_χs_L, JLD2.load_object( agatedatadir( filename ) ) )
			end
			conv_betas[Lvalue] = conv_betas_L
			agate_hyp_χs[Lvalue] = agate_hyp_χs_L
		end
	end
	for key ∈ keys(conv_betas)
		println("Number completed: $(length(conv_betas[key]))")
	end
end

# ╔═╡ 38dba04a-a68f-4bca-9c9a-a97417624e6c
begin
	hyp_χ_ghists = Dict()
	for (Lvalue, agate_hyp_χs_L) ∈ agate_hyp_χs
		ghists_L = []
		for χvals ∈ agate_hyp_χs_L
			ghist = GaussianHistogram()
			for χi ∈ χvals
				new_χ = measurement(χi.val, χi.err)
				push!(ghist, χi )
			end
			push!(ghists_L, ghist)
		end
		hyp_χ_ghists[Lvalue] = ghists_L
	end
end

# ╔═╡ 4cad2e76-747d-45b8-8fbf-8557f9bfb889
begin
	temperatures = Dict()
	for (Lvalue, betas) ∈ conv_betas
		temperatures[Lvalue] = 1 ./ betas
	end
	temperatures
end

# ╔═╡ 6c5781a7-73b3-41ec-8ae8-bbdd3f1f6ad2
let
plt = plot(; xlabel = "Temperature \$T\$", ylabel = "\$ \\mathrm{Av}(1/T_1) \$")

for (TLval, HistLval) ∈ zip( keys(temperatures), keys(hyp_χ_ghists) )
	plot!(plt, temperatures[TLval], mean.(hyp_χ_ghists[HistLval]);
		  yerror = std.(hyp_χ_ghists[HistLval]),
		  label = "\$ L = $(TLval) \$", markershape = :circle)
end

vline!(plt, [Tc]; ls = :dash, color = :orange, lw = 2.5, label = "\$ T_c \$",
	   z_order = :back)
plt
end

# ╔═╡ eb72f5b4-6a30-4a64-82bf-4596b2599966
let
plt = plot(; xlabel = "Temperature \$T\$", ylabel = "\$ \\sigma_1 \$")

Lmin = minimum(keys(hyp_χ_ghists))	
for (TLval, HistLval) ∈ zip( keys(temperatures), keys(hyp_χ_ghists) )
	plot!(plt, temperatures[TLval], std.(hyp_χ_ghists[HistLval]) ;
		  label = "\$ L = $(TLval) \$", markershape = :circle)
end

vline!(plt, [Tc]; ls = :dash, color = :orange, lw = 2.5, label = "\$ T_c \$",
	   z_order = :back)
plt
end

# ╔═╡ f606caae-2cbb-430b-90fb-44491b326cc0
let
plt = plot(; xlabel = "Temperature \$T\$", ylabel = "\$ \\sigma_1 \\sqrt{L^2} \$")

Lmin = minimum(keys(hyp_χ_ghists))	
for (TLval, HistLval) ∈ zip( keys(temperatures), keys(hyp_χ_ghists) )
	plot!(plt, temperatures[TLval], std.(hyp_χ_ghists[HistLval]) .* HistLval ;
		  label = "\$ L = $(TLval) \$", markershape = :circle)
end

vline!(plt, [Tc]; ls = :dash, color = :orange, lw = 2.5, label = "\$ T_c \$",
	   z_order = :back)
plt
end

# ╔═╡ Cell order:
# ╠═b6739750-48d2-11ed-12ad-5bd8f89b0ab6
# ╠═abc6dc6b-f864-447e-bbf1-40daa2d92594
# ╠═d9f59e54-ef96-4fdc-8236-f1fa15a08185
# ╠═4ea5386f-bf4e-4fa2-8023-70d4f9064b76
# ╠═c4c38477-3ab0-4a6c-9a76-dd8f7a712650
# ╠═445f6190-b1c7-4f49-be35-6eae76d184c6
# ╠═8bdd7dbc-2f91-4db8-9fe8-681ee3f831a7
# ╠═5833ae1c-d641-409a-bf5f-1b291ca3ed41
# ╠═61f806fd-6f17-47fa-8b50-6db4179a00ad
# ╠═38dba04a-a68f-4bca-9c9a-a97417624e6c
# ╠═4cad2e76-747d-45b8-8fbf-8557f9bfb889
# ╠═873e6f73-feba-4c2d-9134-0f9065ca0af9
# ╠═6c5781a7-73b3-41ec-8ae8-bbdd3f1f6ad2
# ╠═eb72f5b4-6a30-4a64-82bf-4596b2599966
# ╠═f606caae-2cbb-430b-90fb-44491b326cc0
