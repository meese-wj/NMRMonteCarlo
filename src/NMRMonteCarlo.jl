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

@reexport using .Lattices
@reexport using .Hamiltonians
@reexport using .MonteCarloMethods

end # NMRMonteCarlo