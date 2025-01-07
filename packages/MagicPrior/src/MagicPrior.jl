module MagicPrior

using Statistics, LinearAlgebra, SparseArrays

using Pedigrees

#TODO: prior for random mating schemes = magicOrigin, magicPrior

include("memoize_demoize.jl")
include("pedidentity.jl")
include("pedancestry.jl")
include("magicprior_fun.jl")
include("magicorigin.jl")

end # module
