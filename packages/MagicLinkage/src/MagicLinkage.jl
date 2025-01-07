module MagicLinkage

using Pkg, Distributed, LinearAlgebra
using CSV, DataFrames, SparseArrays, Distributions
using Dates, ProgressMeter, GZip
using Base.Iterators, Random
using LRUCache
using Pedigrees, MagicBase
using MagicPrior, MagicReconstruct, MagicImpute
using MagicLD

export
    read_linkage, plot_linkage,
    # magiclinkage!, 
    magiclinkage

include("pairwise_prior.jl")
include("magiclinkage_fun.jl")

end # module MagicLinkage
