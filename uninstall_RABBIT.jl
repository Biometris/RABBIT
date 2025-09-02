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
try 
    res = 0
    @time for pkg in pkgls
        println("-------------remove ",pkg, "-------------")
        try 
            Pkg.rm(pkg)   
        catch err     
            res = -1       
            @error err        
        end
    end
    res
catch err
    -1
end

