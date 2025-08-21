
try    
    using Pkg
    pnls = [
        "HMM",
        "Pedigrees",
        "MagicPrior",
        "MagicBase",    
        "MagicSimulate",    
        "MagicReconstruct",        
        "MagicFilter",
        "MagicImpute",    
        "MagicCall",
        "MagicMap",
        "MagicScan",
        "RABBITCLI"
    ]
    Pkg.activate()
    @time for pn in pnls        
        println("-------------test ",pn, "-------------")
        @time Pkg.test(pn)             
    end    
    0
catch err
    rethrow(err)
    -1
end

