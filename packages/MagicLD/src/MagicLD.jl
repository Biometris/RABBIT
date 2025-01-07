module MagicLD

using Pkg, Distributions,DelimitedFiles, Distributed
using CSV, DataFrames, SparseArrays, DataStructures
using Dates, ProgressMeter, GZip
using Pedigrees, MagicBase

export    
    magicld!, magicld

include("get_subld.jl")
include("magicld_fun.jl")

end # module MagicLD
