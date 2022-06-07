using Statistics

struct Binner{vec_t <: AbstractVector}
    time_record::vec_t
    binned_averages::Vector{vec_t}
end
bin_depth(length_time::Int) = floor(Int, log2(length_time)) - 1
bin_depth(time_record::AbstractVector) = bin_depth(length(time_record))
bin_depth(bins::Binner) = bin_depth(bins.time_record)

function Binner( time_record::AbstractVector )
    num_t = length(time_record)
    num_levels = bin_depth(num_t)
    return Binner(time_record, [ zeros( num_t ÷ Int(2^ldx) ) for ldx ∈ 1:num_levels ] )
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


function average_halves!(output_avgs::AbstractVector, time_record::AbstractVector)
    for idx ∈ 1:length(output_avgs)
        output_avgs[idx] = 0.5 * ( time_record[2 * idx - 1] + time_record[2 * idx] )
    end
    return nothing
end

function average_halves!(bins::Binner, level::Int)
    if level == 1
        average_halves!(bins.binned_averages[level], bins.time_record)
    else
        average_halves!(bins.binned_averages[level], bins.binned_averages[level - 1])
    end
end

function average_halves!(bins::Binner)
    for ldx ∈ 1:bin_depth(bins)
        average_halves!(bins, ldx)
    end
end

function average_halves(time_record::AbstractVector)
    num_bins = length(time_record) ÷ 2
    avgs = zeros(num_bins)
    average_halves!(avgs, time_record)
    return avgs
end

function bin_record(subrecord::AbstractVector)
    if length(subrecord) == 1
        return subrecord[1]
    end
    return bin_record( average_halves(subrecord) )
end