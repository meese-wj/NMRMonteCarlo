using Random, Distributions

abstract type ResamplingAnalyzer end
struct Bootstrap <: ResamplingAnalyzer end
struct Jackknife <: ResamplingAnalyzer end

resample(typ_e::Type{T} ...) where {T} = throw(ArgumentError("\nNo resampling method defined for $(typ_e) types.\n"))``

function resample(record; samples = length(record), sample_size = length(record))
    record_size = length(record)
    # First make the random indices
    indices = rand(DiscreteUniform(1, record_size), sample_size)
    sort!(indices)
    # Then average
end