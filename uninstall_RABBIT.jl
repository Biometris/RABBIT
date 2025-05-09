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
st = 0
@time for pkg in pkgls
    println("-------------remove ",pkg, "-------------")
    try 
        Pkg.rm(pkg)   
    catch err     
        st = -1       
        @error err
    end
end
st
