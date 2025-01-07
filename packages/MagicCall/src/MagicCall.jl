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
using Dates
using MagicBase
using HMM, Pedigrees, MagicPrior, MagicReconstruct

export 

    magiccall

include("magiccall_fun.jl")

end # module MagicCall
