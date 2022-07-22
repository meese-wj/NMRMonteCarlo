# module MCBinningAnalysis

using Statistics, TelegraphNoise, Plots, LaTeXStrings

# export var_of_mean, bin_size, Binner, analyze!, Rx, τeff

const MIN_NUM_BINS = Int(64)

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
bin_size(bins::Binner) = [ bin_size(level) for level ∈ 1:bin_depth(bins) ]
reliable_bin_size(bins::Binner) = begin arr = bin_size(bins); arr[ length(bins.time_record) ./ arr .>= MIN_NUM_BINS ] end
max_reliable_bin_size(number::Number) = number ÷ MIN_NUM_BINS
max_reliable_bin_size(record) = max_reliable_bin_size(length(record))
max_reliable_bin_size(bins::Binner) = max_reliable_bin_size(bins.time_record)

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
num_bins( bins::Binner ) = [ num_bins(bins.time_record, bin_level) for bin_level ∈ 1:bin_depth(bins) ]

function Binner( record::AbstractVector; analyze = true )
    num_t = length(record)
    max_bin_level = bin_depth(num_t)
    bins = Binner( mean(record), var_of_mean(record), record,
                   [ zeros( num_bins(num_t, level) ) for level ∈ 1:max_bin_level ],
                   zeros(max_bin_level) )
    if analyze
        analyze!(bins)
    end
    return bins
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

function _bin!(output, record)
    @inbounds for idx ∈ 1:length(output)
        output[idx] = 0.5 * ( record[2 * idx - 1] + record[2 * idx] )
    end
    return nothing
end

function _bin!(bins::Binner, level::Int)
    if level > 1
        _bin!(bins.binned_averages[level], bins.binned_averages[level - 1])
        return nothing
    end
    _bin!(bins.binned_averages[1], bins.time_record)
    return nothing
end

_bin!(bins::Binner) = @inbounds for level ∈ 1:bin_depth(bins) _bin!(bins, level) end

function analyze!(bins::Binner)
    _bin!(bins)    
    rec_mean = bins.record_mean
    @inbounds for level ∈ 1:bin_depth(bins)
        bins.level_variances[level] = var( bins.binned_averages[level]; mean = rec_mean )
    end
    return nothing
end

level_var_of_mean( bins::Binner, level ) = bins.level_variances[level] / num_bins( bins, level )
var_of_mean(bins::Binner) = [ level_var_of_mean(bins, level) for level ∈ 1:bin_depth(bins) ]
Rx(bins::Binner) = begin values = var_of_mean(bins); values ./ naive_variance(bins) end
τeff(bins) = 0.5 .* ( Rx(bins) .- 1 )

reldiff(a, b) = abs(a - b) / b
reldiff(a::AbstractVector) = broadcast( (a,b) -> reldiff(a,b), a[2:end], a[1:end-1] )
# end

function bin_plot!( plt1, plt2, bins::Binner )
    plot!(plt1, bin_size(bins), Rx(bins); 
          xscale = :log10, yscale = :log10, 
          xlabel = L"Bin Size $2^\ell$", ylabel = L"$R_X$", 
          markershape = :circle, leg = false )

    plot!(plt2, bin_size(bins)[2:end], reldiff(Rx(bins)); 
          xscale = :log10, yscale = :identity, 
          xlabel = L"Bin Size $2^\ell$", ylabel = L"$\Delta R_X / R_X$", 
          markershape = :circle, label = "\$\\ell_{\\textrm{max}} = $(bin_depth(bins))\$" )
    return nothing
end

plot_max_num_bins( number::Number ) = vline([max_reliable_bin_size(number)]; color = "black", linestyle = :dash, label = "\$N_{\\mathrm{bins}} = $MIN_NUM_BINS\$", alpha = 0.5 )
plot_max_num_bins( record ) = plot_max_num_bins( length(record) )

function bin_plot( record; plot_title = "" )
    plt1 = plot_max_num_bins(record)
    plt2 = plot_max_num_bins(record)
    bin_plot!(plt1, plt2, Binner(record))
    plot(plt1, plt2; layout = (1,2), link = :x, figsize = (800, 450), plot_title = plot_title)
end

function telegraph_plots( dwell_time, minpow2, maxpow2, incr = 2 )
    plt1 = plot_max_num_bins(2^minpow2)
    plt2 = plot_max_num_bins(2^minpow2)
    for pow2 ∈ minpow2:incr:maxpow2
        @show 2^pow2
        pt_binner = Binner( Telegraph(dwell_time, Int(2^pow2)).signal )
        bin_plot!(plt1, plt2, pt_binner)
    end
    plot(plt1, plt2; layout = (1,2), link = :x, figsize = (800, 450))
end