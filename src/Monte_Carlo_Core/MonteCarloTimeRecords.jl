using Statistics

"""
    var_of_mean(record::AbstractVector)

Calculate the variance of the mean for a given time `record`.

```jldoctest
julia> record = [1, 2, 3, 4, 5];

julia> var_of_mean(record)
0.5

```
"""
var_of_mean( record::AbstractVector ) = var(record; corrected=false) / ( length(record) - 1 )

struct Binner{vec_t <: AbstractVector}
    record_mean::Float64
    record_variance::Float64
    time_record::vec_t
    binned_averages::Vector{vec_t}
    level_variances::vec_t  # store the variances here
end
"""
    bin_size(bin_level) → Int

Number of samples in a bin at a given `bin_level`.

```jldoctest
julia> bin_size(4)
16

```
"""
bin_size(bin_level) = Int(2^bin_level)
"""
    num_bins(num_measurements, bin_level) → Int 

Number of bins at a given `bin_level`. This assumes that each bin is of size `bin_size(bin_level)`.

```jldoctest
julia> num_bins(4096, 4)
256

```
"""
num_bins( num_measurements, bin_level ) = floor(Int, num_measurements / bin_size(bin_level))
num_bins( record::AbstractVector, bin_level ) = num_bins(length(record), bin_level)
num_bins( bins::Binner, bin_level ) = num_bins(bins.time_record, bin_level)

function Binner( record::AbstractVector )
    num_t = length(record)
    max_bin_level = bin_depth(num_t)
    return Binner( mean(record), var_of_mean(record), record,
                   [ zeros( num_bins(num_t, level) ) for level ∈ 1:max_bin_level ],
                   zeros(max_bin_level) )
end

"""
    bin_depth(length_time::Int)

Total number of binning levels attainable from a length of time `length_time`.

```jldoctest
julia> bin_depth(4096)
11

```
"""
bin_depth(length_time::Int) = floor(Int, log2(length_time)) - 1
"""
    bin_depth(time_record::AbstractVector)

Total number of binning levels attainable in the `time_record`.
"""
bin_depth(time_record::AbstractVector) = bin_depth(length(time_record))
"""
    bin_depth(bins::Binner)

Total number of binning levels attainable for a given `Binner` struct.
"""
bin_depth(bins::Binner) = bin_depth(bins.time_record)

naive_variance(bins::Binner) = bins.record_variance

function analyze!(bins::Binner)
    rec_size = length(bins.time_record)
    rec_mean = bins.record_mean
    for level ∈ 1:bin_depth(bins), bindex ∈ 1:num_bins(rec_size, level)
        starter = (bindex - 1) * bin_size(level) + 1
        ender = bindex * bin_size(level)
        bins.binned_averages[level][bindex] = mean( @view bins.time_record[ starter:ender ] )
    end

    for level ∈ 1:bin_depth(bins)
        @show bins.level_variances[level] = var( bins.binned_averages[level]; mean = rec_mean )
    end
    return nothing
end

level_var_of_mean( bins::Binner, level ) = bins.level_variances[level] / num_bins( bins, level )
var_of_mean(bins::Binner) = [ level_var_of_mean(bins, level) for level ∈ 1:bin_depth(bins) ]
Rx(bins::Binner) = begin values = var_of_mean(bins); values ./ naive_variance(bins) end
τeff(bins) = 0.5 .* ( Rx(bins) .- 1 )

# function Binner( time_record::AbstractVector )
#     num_t = length(time_record)
#     num_levels = bin_depth(num_t)
#     return Binner(time_record, [ zeros( num_t ÷ Int(2^ldx) ) for ldx ∈ 1:num_levels ] )
# end

function generate_poisson( τ, signal_length )
    signal::Vector{Float64} = []
    stepsize = floor(Int, -τ * log( 1 - rand() ))
    append!(signal, ones(stepsize))
    while length(signal) < Int(signal_length)
        stepsize = floor(Int, -τ * log( 1 - rand() ))
        if signal[end] == 1
            append!(signal, zeros(stepsize))
        else
            append!(signal, ones(stepsize))
        end
    end
    return signal[1:signal_length]
end

mutable struct BinnedStatistics{T <: Real}
    record_stat::T
    binned_stats::Vector{T}
end
function BinnedStatistics( bins::Binner )
    value_type = eltype(bins.time_record)
    return BinnedStatistics( zero(value_type), zeros(value_type, bin_depth(bins)) )
end

function BinnedStatistics( bins::Binner, statistic::Function )
    bstat = BinnedStatistics(bins)
    bstat.record_stat = statistic(bins.time_record)
    for ldx ∈ 1:bin_depth(bins)
        bstat.binned_stats[ldx] = statistic(bins.binned_averages[ldx])
    end
    return bstat
end

bin_sizes( bstat::BinnedStatistics ) = [ 2^ldx for ldx ∈ length(bstat.binned_stats):-1:1 ]
bin_sizes( bins::Binner ) = [ 2^ldx for ldx ∈ bin_depth(bins):-1:1 ]
bin_variance( values::AbstractVector ) = var(values; corrected=false) / ( length(values) - 1 )
bin_variance( bins::Binner ) = [ bin_variance(bin_values) for bin_values ∈ bins.binned_averages ]



# _Rx_numerator( bstat::BinnedStatistics ) = (2 .* bin_sizes(bstat)) .* bstat.binned_stats
# _Rx_denominator( bins::Binner ) = var( bins.time_record; corrected=false )
# # _Rx_denominator( bins::Binner ) = bin_variance(bins.time_record)
# Rx_values(bstat::BinnedStatistics, bins::Binner) = _Rx_numerator(bstat) ./ _Rx_denominator(bins)
# Rx_values(bins::Binner) = begin bstat = BinnedStatistics( bins, bin_variance ); Rx_value(bstat, bins) end
# reldiff(a, b) = (a - b) / b
# function Rx_limit( rvalues )
#     # Take the last three Rₓ values.
#     # TODO: Include some better test of convergence...
#     # * Note: rvalues[1], rvalues[2] would likely be 
#     #   where the converged values are
#     # if length(rvalues) < 3
#     #     throw(ArgumentError("Tests of convergence require length(Rₓ) ≥ 3."))
#     # end
#     return rvalues[1]
# end
# τ_stat( rvalues ) = 0.5 * ( Rx_limit(rvalues) - 1 )

# function average_halves!(output_avgs::AbstractVector, time_record::AbstractVector)
#     for idx ∈ 1:length(output_avgs)
#         output_avgs[idx] = 0.5 * ( time_record[2 * idx - 1] + time_record[2 * idx] )
#     end
#     return nothing
# end

# function average_halves!(bins::Binner, level::Int)
#     if level == 1
#         average_halves!(bins.binned_averages[level], bins.time_record)
#     else
#         average_halves!(bins.binned_averages[level], bins.binned_averages[level - 1])
#     end
# end

# function average_halves!(bins::Binner)
#     for ldx ∈ 1:bin_depth(bins)
#         average_halves!(bins, ldx)
#     end
# end

# function average_halves(time_record::AbstractVector)
#     num_bins = length(time_record) ÷ 2
#     avgs = zeros(num_bins)
#     average_halves!(avgs, time_record)
#     return avgs
# end

# function bin_record(subrecord::AbstractVector)
#     if length(subrecord) == 1
#         return subrecord[1]
#     end
#     return bin_record( average_halves(subrecord) )
# end