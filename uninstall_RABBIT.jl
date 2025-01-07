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
    @time for pkg in pkgls
        println("Uninstall ", pkg)    
        Pkg.rm(pkg)                    
    end
    0
catch err
    @error err
    -1
end