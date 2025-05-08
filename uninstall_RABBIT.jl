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
    @time for pkg in pkgls
        println("-------------remove ",pkg, "-------------")
        Pkg.rm(pkg)                    
    end
    0
catch err
    @error err
    -1
end