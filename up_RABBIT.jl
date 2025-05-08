try     
    using Pkg
    pkgls = [
        "HMM",
        "Pedigrees",
        "MagicPrior",
        "MagicBase",    
        "MagicSimulate",    
        "MagicReconstruct",        
        "MagicFilter",
        "MagicImpute",    
        "MagicCall",
        "SpectralEmbedding",
        "MagicLD",
        "MagicLinkage",
        "MagicMap",
        "MagicScan",
        "RABBITCLI"
    ]
    Pkg.activate()
    for pn in pkgls            
        println("-------------update ",pn, "-------------")
        Pkg.update(pn)    
    end
    using RABBITCLI
    @info rabbitversion()
    0
catch err
    @error err
    -1
end

