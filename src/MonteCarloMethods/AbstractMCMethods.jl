using Parameters2JSON

export sweeps_per_export, Model

abstract type AbstractMonteCarloParameters end
sweeps_per_export(params::AbstractMonteCarloParameters) = params.measure_sweeps <= params.total_measurements ? 1 : params.measure_sweeps ÷ params.total_measurements

# TODO: Change this to an abstract type
struct Model{L <: AbstractLattice, H <: AbstractHamiltonian}
    latt::L
    ham::H
end