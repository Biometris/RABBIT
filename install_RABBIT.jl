
try
    using Pkg
    @info string("Start to develop RABBIT packages-----------")    
    filedir = dirname(@__FILE__)
    println("Install RABBIT from ",filedir)    
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
        println("-------------develop ", pkg, "-------------")
        Pkg.develop(PackageSpec(path=joinpath(filedir,"packages", pkg)))                        
    end        
    using RABBITCLI
    @info rabbitversion()
    0
catch err
    @error err
    -1
end
