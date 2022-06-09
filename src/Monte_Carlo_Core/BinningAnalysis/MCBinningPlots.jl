using RecipesBase
include("MCBinningAnalysis.jl")

# abstract type BinningPlotType end
# struct RxValues <: BinningPlotType end
# struct τeffValues <: BinningPlotType end

# @recipe f(bp_t::Type{BinningPlotType}, args...) = error("\nNo plotting routine defined for $(typeof(bp_t))") 

@macroexpand @recipe function f(::Type{Val{:RxValues}}, x, record::AbstractVector, z)
    # bins = MCBinningAnalysis.Binner(record; analyze = true)
    # bsizes = MCBinningAnalysis.bin_sizes(bins)
    # seriestype --> :line 
    # markershape --> :circle
    # x := bsizes
    # y := MCBinningAnalysis.Rx(bins)
    # ()
    @series begin
        x := record
        y := record
    end
end

@macroexpand @recipe function f(::Type{Val{:seriespie}}, x, y, z)
    framestyle --> :none
    aspect_ratio --> true
    s = sum(y)
    θ = 0
    for i in eachindex(y)
        θ_new = θ + 2π * y[i] / s
        coords = [(0.0, 0.0); Plots.partialcircle(θ, θ_new, 50)]
        @series begin
            seriestype := :shape
            label --> string(x[i])
            x := first.(coords)
            y := last.(coords)
        end
        θ = θ_new
    end
end

# Plots.plot(bp_t::Type{BinningPlotType}, args...) = error("\nNo plotting routine defined for $(typeof(bp_t))") 

# function Plots.plot(::Type{RxValues}, record; kwargs...)
#     bins = Binner(record; analyze = true)
#     bsizes = bin_size( bins )
#     Plots.plot( bsizes, Rx(bins);
#                 xlabel = L"Bin Size $2^\ell$", 
#                 ylabel = L"$R_X$",
#                 markershape = :circle,
#                 kwargs... )
# end