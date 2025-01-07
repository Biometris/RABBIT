
"""
    MagicImpute 

a package for genotype imputation in multiparental populations. 
Export two functions:  [`magicimpute`](@ref) and [`magicimpute!`](@ref).
"""
module MagicImpute

using Base.Iterators, Distributed, Dates
using DataFrames, CSV
using SparseArrays, LinearAlgebra, Distributions, Random
using JLD2, GZip, DelimitedFiles, ProgressMeter
import DataStructures: OrderedDict
using MagicBase
using HMM, Pedigrees, MagicPrior, MagicReconstruct

export
    # private

    # pulic    
    magicmask, magicmask_impute, 
    magicimpute,magicimpute!

include("fix_nonsegregate.jl")
include("fhaploprior.jl")
include("founderimpute.jl")
include("foundercorrect.jl")
include("infer_error_logit.jl")
include("infer_error_logit_viterbi.jl")
include("infer_error_perind.jl")
include("markerdelete.jl")
include("markerspace.jl")
include("markerskeleton.jl")
include("markerorder.jl")
include("magicimpute_fun.jl")
include("magicimpute_founder_chr.jl")
include("magicimpute_founder.jl")
include("magicimpute_offspring.jl")
include("magicmask.jl")
include("magicmask_impute.jl")

# TODO: error message if poscm is missing

end # module
