"""
MagicScan 

a package for genomic QTL scan in connected multiparental populations. 
Export one function: [`magicscan`](@ref). 
"""
module MagicScan

using Pedigrees, MagicBase
using StatsModels, StatsBase, Distributions, LinearAlgebra
using CSV,DataFrames,Dates
using SparseArrays
using ThreadsX
using Plots
using StatsPlots

export
    #private
    # get_haploprob,fit_lm, scan_lm, getpeak,


    #public
    @formula, # StatsModels
    magicscan,plot_manhattan,plot_qqpval

include("scan_lm.jl")
include("scan_plot.jl")
include("magicscan_fun.jl")

end # module

