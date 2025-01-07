"""
    MagicSimulate

a package for simulating genotyping data in multiparental populations. 
See also [`magicsimulate`](@ref).
"""
module MagicSimulate

using DataFrames, DelimitedFiles, StatsBase
using DataStructures, CSV, Distributions
using Dates
using GZip
using Plots
using Plots.Measures
using Pedigrees, MagicBase
using LinearAlgebra: dot
using Random: randperm,shuffle!,shuffle

export
    # priviate
    # get_breakpoints,gridhomologfgl,hist_inbreeding,hist_recomden,
    rand_mateped, rand_magicped,    

    # public
    simfhaplo, magicsimulate

include("simfhaplo.jl")
include("randmate.jl")
include("traitqtl.jl")
include("simfgl.jl")
include("genodrop.jl")
include("simgeno.jl")
include("simpheno.jl")
include("magicsimulate_fun.jl")

end # module
