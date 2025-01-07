
try     
    filedir = abspath(dirname(@__FILE__),"packages")    
    println("Install RABBIT from ",filedir)
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
        println("Install ", pkg)  
        Pkg.develop(PackageSpec(path=joinpath(filedir,pkg)))                        
    end
    0
catch err
    @error err
    -1
end