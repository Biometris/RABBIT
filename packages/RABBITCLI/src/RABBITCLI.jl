
module RABBITCLI

using ArgParse

using Pkg, Distributed, Distributions

using HMM, Pedigrees, MagicPrior, SpectralEmbedding

using MagicBase, MagicSimulate, MagicReconstruct, MagicFilter, MagicCall

using MagicImpute, MagicLD, MagicLinkage, MagicMap,  MagicScan

using DataFrames

export 
    rabbitversion

function rabbitversion()
    pkgls = [
        "HMM",
        "Pedigrees",
        "MagicPrior",
        "MagicBase",    
        "MagicSimulate",    
        "MagicReconstruct",
        "MagicCall",
        "MagicFilter",
        "MagicImpute",    
        "SpectralEmbedding",
        "MagicLD",
        "MagicLinkage",
        "MagicMap",
        "MagicScan",
        "RABBITCLI"        
    ]
    vls = [MagicBase.pkgversion_st(pn) for  pn in pkgls]
    msg = string("RABBIT version: ",last(vls), ", same as the version of RABBITCLI")
    @info msg
    DataFrame("pkgname"=>pkgls,"version"=>vls)
end

end # module RABBITCLI
