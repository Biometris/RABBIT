"""
    MagicBase 

a package for basic data structures and functions for genetic analysis in connected multiparental populations. 
"""
module MagicBase

using Pkg, LinearAlgebra, Dates
using CSV, DataFrames, SparseArrays
using DelimitedFiles: readdlm, writedlm
using Random, Distributions
using ThreadsX
using Printf: @sprintf
using Statistics: quantile
using StatsBase
using Combinatorics: combinations
using DataStructures: OrderedDict
using Interpolations
using GZip, Tar, CodecZlib, JLD2
# https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/12
using Plots
using Pedigrees
using Images: load
using Graphs: SimpleGraph, connected_components
using InlineStrings

export
    
    exportall, logit, 
    # for pedinfo
    JuncDist, MateScheme, DesignInfo, MagicPed, 
    parsedesign,
    formmagicped,readmagicped,merge_pedfiles, savemagicped, plotmagicped, generate_magicped,   
    pedfile_designcode2ped,    
    # for genofile
    MagicGeno,
    formmagicgeno,readmagicgeno,savemagicgeno, savegenodata, 
    submagicgeno, submagicgeno!,
    readmarkermap, setmarkermap!, 
    splitby_chromosome!,merge_chromosome!,       
    # for ancestryfile
    MagicAncestry,
    readmagicancestry,savemagicancestry,
    # for inputfiles
    parsebreedped, magicparse, check_infiles, 
    rabbitgeno_mma2jl, rabbitped_mma2jl, 
    vcffilter, merge_vcffiles, 
    vcf_count_markers, vcf_del_samples,vcf_rename_samples,
    arrayfile2vcf, merge_arrayfiles, 
    # for likeparameters
    LikeParameters, parse_likeparameters, 
    ThreshLikeParameters, parse_threshlikeparameters,
    PriorLikeParameters, parse_priorlikeparameters,    
    # for accuracy
    magicaccuracy,magicaccuracy!,mapaccuracy,
    imputeaccuracy, imputeaccuracy!,
    # for visualize
    plotmarkermap,plotmarkererror,
    saveprobplot, plotcondprob, animcondprob, 
    plotmosaic, histmosaic, savemosaic,
    plotrecomheat,
    # other I/O
    memoryuse, getabsfile, printconsole, printpkgst,    
    create_targz, extract_targz,
    loadtarplot, 
    info_argnames,
    # other related functions
    rawgenocall!,rawgenoprob!    


function __init__()    
    ENV["GKSwstype"]="nul"
    # @info string("MagicBase init: ENV[\"GKSwstype\"]=",ENV["GKSwstype"])
end    

include("base/localbrent.jl")
include("base/pathformat.jl")
include("base/basis.jl")
include("base/vcfbase.jl")
include("base/likeparameters.jl")
include("base/getfindexlist.jl")
include("magicped/parsedesign.jl")
include("magicped/magicped.jl")
include("magicped/parsebreedped.jl")
include("magicped/reset_args.jl")
include("magicgeno/parsegeno.jl")
include("magicgeno/parsevcf.jl")
include("magicgeno/magicgeno.jl")
include("magicgeno/get_dosegeno.jl")
include("magicgeno/parse_missing.jl")
include("magicgeno/check_infiles.jl")
include("magicgeno/magicparse.jl")
include("magicgeno/rawgenocall.jl")
include("magicancestry/magicancestry.jl")
include("magicancestry/magicancestry_tocontfgl.jl")
include("convert/convert_rabbit_mma.jl")
include("convert/convert_joinmap_cp.jl")
include("convert/convert_hapmap.jl")
include("convert/convert_snparray.jl")
include("convert/convertto_mergemap.jl")
include("accuracy/magicaccuracy.jl")
include("accuracy/mapaccuracy.jl")
include("accuracy/imputeaccuracy.jl")
include("visualize/plotgenofreq.jl")
include("visualize/plotibd.jl")
include("visualize/plotmap.jl")
include("visualize/ploterror.jl")
include("visualize/plotmosaic.jl")
include("visualize/read_ld_linkage.jl")
include("visualize/plotrecomheat.jl")


end # module MagicBase
