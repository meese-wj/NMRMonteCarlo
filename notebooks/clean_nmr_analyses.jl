### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ bf2d7db5-eaa1-4631-879e-c07400d2cd6d
using DrWatson

# ╔═╡ e18e91ba-791b-4bde-825e-502eee274aee
@quickactivate :NMRMonteCarlo

# ╔═╡ 69fa0e30-1d83-11ed-1719-6d028032aa3c
md"""
# NMR Analyses of the Clean Ashkin-Teller Model

In this [Pluto Notebook](https://github.com/fonsp/Pluto.jl), we study the NMR spectra in BaFe₂As₂ as modelled by a clean Ashkin-Teller model.

The goals of this Notebook are to understand how the mean and widths of the NMR spectra change as functions of temperature ``T``, linear system size ``L``, and Monte Carlo time ``L_\tau``. It is crucial to understand these specific behaviors before we we introduce random strain into the problem, and model it with the [random-Baxter-field model (RBFM)](https://arxiv.org/abs/2112.05769).
"""

# ╔═╡ a73ee872-a1af-4223-83d5-47a4a7841875
md"""
## Activate and load the `:NMRMonteCarlo` package
"""

# ╔═╡ 5f32ec73-2d1d-42bd-9561-8a440bf05f65
md"""
## Create a first NMR model
"""

# ╔═╡ 560f5e5b-3d5e-4605-b3fb-d4ba3365d7f8
metroparams, model = let
	Temp, Lx, Kex = 2.269, 16, 0.0
	latt = CubicLattice2D(Lx, Lx)
	atparams = AshkinTellerParameters(1.0, Kex)
	ham = BasicAshkinTellerHamiltonian(latt, atparams)
	metroparams = MetropolisParameters{Float64}([1/Temp], 2^16, 2^18, 2^16)

	model = CleanNMRAshkinTellerModel( latt.params.Lx, latt.params.Ly,
									   atparams.Jex, atparams.Kex, 
									   metroparams.total_measurements )
	metroparams, model
end

# ╔═╡ c9be49e2-0581-46b3-9e3a-a6ff61434a9c
md"""
Now the new Metropolis simulation parameters and clean NMR Ashkin-Teller model are available to use in the Main namespace as `metroparams` and `model`, respectively.
"""

# ╔═╡ 86a7273b-bf28-400f-903e-f11fe71691d0


# ╔═╡ 677439a7-3a1e-4e46-b833-4ff9cbdddc80


# ╔═╡ Cell order:
# ╟─69fa0e30-1d83-11ed-1719-6d028032aa3c
# ╟─a73ee872-a1af-4223-83d5-47a4a7841875
# ╠═bf2d7db5-eaa1-4631-879e-c07400d2cd6d
# ╠═e18e91ba-791b-4bde-825e-502eee274aee
# ╟─5f32ec73-2d1d-42bd-9561-8a440bf05f65
# ╠═560f5e5b-3d5e-4605-b3fb-d4ba3365d7f8
# ╟─c9be49e2-0581-46b3-9e3a-a6ff61434a9c
# ╠═86a7273b-bf28-400f-903e-f11fe71691d0
# ╠═677439a7-3a1e-4e46-b833-4ff9cbdddc80
