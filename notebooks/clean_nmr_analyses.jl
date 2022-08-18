### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ bf2d7db5-eaa1-4631-879e-c07400d2cd6d
using DrWatson

# ╔═╡ e18e91ba-791b-4bde-825e-502eee274aee
@quickactivate "NMRMonteCarlo"

# ╔═╡ 03f49b47-ef03-4b8a-97a3-9f1e84250ff5
using NMRMonteCarlo

# ╔═╡ f6ac5743-88b6-422d-87d8-5bfc8bae1493
using Plots, UncertainHistogramming, Measurements

# ╔═╡ 69fa0e30-1d83-11ed-1719-6d028032aa3c
md"""
# NMR Analyses of the Clean Ashkin-Teller Model

In this [Pluto Notebook](https://github.com/fonsp/Pluto.jl), we study the NMR spectra in BaFe₂As₂ as modelled by a clean Ashkin-Teller model.

The goals of this Notebook are to understand how the mean and widths of the NMR spectra change as functions of temperature ``T``, linear system size ``L``, and Monte Carlo time ``L_\tau``. It is crucial to understand these specific behaviors before we we introduce random strain into the problem, and model it with the [random-Baxter-field model (RBFM)](https://arxiv.org/abs/2112.05769).
"""

# ╔═╡ a73ee872-a1af-4223-83d5-47a4a7841875
md"""
## Activate and load the `:NMRMonteCarlo` package

Start by loading the [`DrWatson`.jl](https://juliadynamics.github.io/DrWatson.jl/dev/) helper package as well as the NMR simulation code base.

!!! note
	Unlike normal Julia scripts (`*.jl`), one needs to explicitly write out each of the following lines in a separate Pluto cell:

	```julia
	using DrWatson
	@quickactivate "NMRMonteCarlo"
	using NMRMonteCarlo
	```

	!!! tip
		In normal Julia scripts, one only needs to write
	
		```julia
		using DrWatson; @quickactivate :NMRMonteCarlo
		```
	
		to do the same thing.
"""

# ╔═╡ 5f32ec73-2d1d-42bd-9561-8a440bf05f65
md"""
## Create a first NMR model

We will first generate a basic NMR simulation using a clean Ashkin-Teller model while choosing all of the default parameters.
"""

# ╔═╡ 560f5e5b-3d5e-4605-b3fb-d4ba3365d7f8
sim = CleanNMRATMSimulation(; Lx = 24)

# ╔═╡ c9be49e2-0581-46b3-9e3a-a6ff61434a9c
md"""
Now the new `struct` is built and ready to be used. 

One can access its `parameters` and `model` with the help of `SimulationParameters` and `SimulationModel`, respectively.  Similarly, the Monte Carlo updating method can be accessed by `SimulationMethod` as:
"""

# ╔═╡ 45da1ff2-6ab7-4cc6-a40e-e6430b72d016
SimulationMethod(sim)

# ╔═╡ 59e3f6f1-14e0-4b9d-bd29-3c77695b0430
md"""
Now we `simulate!` it.
"""

# ╔═╡ 35ceed2d-af0f-4331-ab4c-f493c547a522
@time simulate!(sim)

# ╔═╡ 0ec7b308-5ba2-43b1-bb21-760ef6e48cbb
md"""
## Now we analyze the simulation data
"""

# ╔═╡ d9c014df-27f5-4b3b-b910-6f720f0323c0
results = analyze(sim)

# ╔═╡ b6183421-ab2a-4b8a-be0e-a1baedb5cef6
md"""
To plot, the resulting NMR proxy spectra, we can simply use `Plots.jl` in collaboration with `UncertainHistogramming.jl`. We choose to use both a `GaussianHistogram` as well as a `UniformHistogram` for the analyses.
"""

# ╔═╡ 92a14544-35d1-476a-bceb-b7e5e15e97d8
start_idx = length(Observables(sim).time_series_obs) + one(Int)

# ╔═╡ 56bc3875-94e6-4607-966a-e8e4f7655861
Ωrange = LinRange(1.5, 2.2, 5000)

# ╔═╡ 4042a83f-a819-4a04-889a-a408491c466b
begin
	ghist = GaussianHistogram();
	for idx ∈ start_idx:length(Observables(sim)) 
		push!(ghist, measurement(results[idx]))
	end
end

# ╔═╡ b75737a2-9842-4ecc-9d1b-bfb300a9bd7a
ghist

# ╔═╡ 5fa7ae24-cfcc-487c-bc49-3fcb2cbaa5b0
gplot = plot( Ωrange, ghist; 
	  		  xlabel = "\$\\Omega\$", 
	  		  ylabel = "\$\\mathcal{H}(\\Omega)\$",
	  		  title = "\$ \\mathtt{GaussianHistogram}ming \$")

# ╔═╡ 9b12b50d-36b4-4f8b-8f28-1ea1559fb744
md"""
Similarly, we can make the Continuous Histogram with a uniform kernel instead of a Gaussian one. 
"""

# ╔═╡ 38657a63-b757-4bba-b765-32c4b77e2545
begin
	uhist = UniformHistogram();
	for idx ∈ start_idx:length(Observables(sim)) 
		push!(uhist, measurement(results[idx]))
	end
end

# ╔═╡ 5980e488-e329-4c94-86cd-a800e026ba00
uhist

# ╔═╡ 7cf992b6-951d-48a9-8eaf-fa4e9185bd01
uplot = plot( Ωrange, uhist; 
	  		  xlabel = "\$\\Omega\$", 
	  		  ylabel = "\$\\mathcal{H}(\\Omega)\$",
	  		  title = "\$ \\mathtt{UniformHistogram}ming \$")

# ╔═╡ d50730bd-9959-4213-90f1-a20d91f0e1d4
md"""
Then, we can plot these distributions next to one another as:
"""

# ╔═╡ e7ab25a6-0bb5-4e67-a7d2-0c88a6128d29
plot( gplot, uplot; link = :both )

# ╔═╡ Cell order:
# ╟─69fa0e30-1d83-11ed-1719-6d028032aa3c
# ╟─a73ee872-a1af-4223-83d5-47a4a7841875
# ╠═bf2d7db5-eaa1-4631-879e-c07400d2cd6d
# ╠═e18e91ba-791b-4bde-825e-502eee274aee
# ╠═03f49b47-ef03-4b8a-97a3-9f1e84250ff5
# ╟─5f32ec73-2d1d-42bd-9561-8a440bf05f65
# ╠═560f5e5b-3d5e-4605-b3fb-d4ba3365d7f8
# ╟─c9be49e2-0581-46b3-9e3a-a6ff61434a9c
# ╠═45da1ff2-6ab7-4cc6-a40e-e6430b72d016
# ╟─59e3f6f1-14e0-4b9d-bd29-3c77695b0430
# ╠═35ceed2d-af0f-4331-ab4c-f493c547a522
# ╟─0ec7b308-5ba2-43b1-bb21-760ef6e48cbb
# ╠═d9c014df-27f5-4b3b-b910-6f720f0323c0
# ╟─b6183421-ab2a-4b8a-be0e-a1baedb5cef6
# ╠═92a14544-35d1-476a-bceb-b7e5e15e97d8
# ╠═56bc3875-94e6-4607-966a-e8e4f7655861
# ╠═f6ac5743-88b6-422d-87d8-5bfc8bae1493
# ╠═4042a83f-a819-4a04-889a-a408491c466b
# ╠═b75737a2-9842-4ecc-9d1b-bfb300a9bd7a
# ╠═5fa7ae24-cfcc-487c-bc49-3fcb2cbaa5b0
# ╟─9b12b50d-36b4-4f8b-8f28-1ea1559fb744
# ╠═38657a63-b757-4bba-b765-32c4b77e2545
# ╠═5980e488-e329-4c94-86cd-a800e026ba00
# ╠═7cf992b6-951d-48a9-8eaf-fa4e9185bd01
# ╟─d50730bd-9959-4213-90f1-a20d91f0e1d4
# ╠═e7ab25a6-0bb5-4e67-a7d2-0c88a6128d29
