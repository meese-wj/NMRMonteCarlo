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
sim = CleanNMRATMSimulation()

# ╔═╡ c9be49e2-0581-46b3-9e3a-a6ff61434a9c
md"""
Now the new `struct` is built and ready to be used. 

One can access its `parameters` and `model` with the help of `SimulationParameters` and `SimulationModel`, respectively.  Similarly, the Monte Carlo updating method can be accessed by `SimulationMethod` as:
"""

# ╔═╡ 45da1ff2-6ab7-4cc6-a40e-e6430b72d016
SimulationMethod(sim)

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
