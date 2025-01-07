"""
    MagicMap

a package for genetic map construction in connected multiparental populations. 
Export one function:  [`magicmap`](@ref).
"""
module MagicMap

using CSV, DataFrames, Distributed
using Distributions
using LinearAlgebra, SparseArrays, StatsBase
using JuMP, OSQP, Graphs
using SpectralEmbedding 
# using HMM, Pedigrees, MagicBase, MagicPrior
using MagicBase
using MagicLD, MagicLinkage
using MagicReconstruct, MagicImpute
using GZip, Dates
using ProgressMeter
using Base.Iterators
using Random
using ThreadsX

export
    # private
    
    # public
    magicmap


include("binning.jl")
include("grouping.jl")
include("construct.jl")
include("magicmap_fun.jl")


end # module
