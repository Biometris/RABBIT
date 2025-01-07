try     
    using Pkg
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
    Pkg.activate()
    for pn in pkgls    
        println("Update ", pn)
        Pkg.update(pn)    
    end
    using RABBITCLI
    @info rabbitversion()
    0
catch err
    @error err
    -1
end

