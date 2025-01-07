"""
    MagicReconstruct 

a package for haplotype reconstruction in connected multiparental populations. 
Export two functions:  [`magicreconstruct`](@ref) and [`magicreconstruct!`](@ref).
"""
module MagicReconstruct

using LinearAlgebra, DelimitedFiles, Distributed
using DataStructures,Statistics, DataFrames, SparseArrays, CSV,Dates
using StatsPlots   # no need for `using Plots` as that is reexported here
using SpecialFunctions
using JLD2, GZip
using ProgressMeter
using MagicBase
using HMM, Pedigrees, MagicPrior

export
    # private    

    # public
    magicreconstruct, magicreconstruct!


include("prior_continuous_time.jl")
include("prior_discrete_time.jl")
include("likelihood_ind.jl")
include("likelihood_ind_gbs.jl")
include("likelihood_site.jl")
include("likelihood_precompute.jl")
include("hmm_decode.jl")
include("hmm_decode_permarker.jl")
include("reconstruct.jl")
include("magicreconstruct_fun.jl")
include("plot_recom_inbred.jl")

end # module
