"""
    module NMRMonteCarlo

Main `src` entry point for the simulation of nuclear magnetic resonance
(NMR) signals in the _Fast Dynamics Limit_ using Monte Carlo Methods.
This limit holds whenever the typical relaxation time is much smaller
than the period of oscillation associated with 
[Larmor precession](https://en.wikipedia.org/wiki/Larmor_precession?oldformat=true).
"""
module NMRMonteCarlo

import Reexport: @reexport
include( joinpath( "Lattices", "Lattices.jl" ) )
include( joinpath( "Hamiltonians", "Hamiltonians.jl" ) )
include( joinpath( "MonteCarloMethods", "MonteCarloMethods.jl" ) )
include( joinpath( "SimulatingNMR", "SimulatingNMR.jl" ) )
include( joinpath( "Models", "Models.jl" ) )
include( joinpath( "Simulations", "Simulations.jl" ) )
include( joinpath( "Utilities", "DrWatsonHelpers.jl" ) )

@reexport using .Lattices
@reexport using .Hamiltonians
@reexport using .MonteCarloMethods
@reexport using .SimulatingNMR
@reexport using .Models
@reexport using .Simulations
@reexport using .DrWatsonHelpers

end # NMRMonteCarlo