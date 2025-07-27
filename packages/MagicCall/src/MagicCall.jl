"""
    MagicCall 

a package for single site genotyping calling from sequence data in connected multiparental populations. 
Export one function: [`magiccall`](@ref).
"""
module MagicCall

using GZip
using Distributed
using Distributions
using LinearAlgebra, SparseArrays
using Dates, ThreadsX
using MagicBase
using HMM, Pedigrees, MagicPrior, MagicReconstruct

export 

    magiccall

include("magiccall_fun.jl")
include("shift_outliers.jl")

end # module MagicCall
