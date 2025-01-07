
try    
    using Pkg
    pnls = [
        "HMM",
        "Pedigrees",
        "MagicPrior",
        "MagicBase",    
        "MagicSimulate",    
        "MagicReconstruct",
        "MagicCall",
        "MagicFilter",
        "MagicImpute",    
        "MagicMap",
        "MagicScan",
        "RABBITCLI"
    ]
    Pkg.activate()
    @time for pn in pnls
        println("Test ", pn)    
        @time Pkg.test(pn)             
    end
    0
catch err
    @error err
    -1
end
