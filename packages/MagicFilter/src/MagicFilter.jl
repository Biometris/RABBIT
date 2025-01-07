"""
    MagicFilter 

a package for filtering markers and individuals in connected multiparental populations. 
Export one function:  [`magicfilter`](@ref).
"""
module MagicFilter

using CSV, DataFrames, LinearAlgebra
using HypothesisTests
using Distributions, Distributed, SparseArrays, StatsBase
using GZip
using ThreadsX
using DataStructures: OrderedDict
using ProgressMeter, LightGraphs
using Plots
import HMM, Pedigrees, MagicPrior
using MagicBase
using MagicReconstruct

export
    # private

    # public
    magicfilter,magicfilter!

include("filter_marker.jl")
include("filter_offspring_missing.jl")
include("filter_offspring_dupe.jl")
include("filter_subpop_size.jl")
include("filter_founder_progeny.jl")
include("filter_loglike_eps.jl")
include("filter_infer_eps.jl")
include("magicfilter_fun.jl")

end # module
