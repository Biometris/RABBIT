
try 
    using Pkg
    pkgdeps  = [
        "HMM" => [],
        "Pedigrees" => [],
        "MagicPrior" => ["Pedigrees"],
        "MagicBase" => ["Pedigrees"],    
        "MagicSimulate" => ["Pedigrees","MagicBase"],    
        "MagicReconstruct" => ["HMM","Pedigrees","MagicBase","MagicPrior"],        
        "MagicFilter" => ["HMM","Pedigrees","MagicBase","MagicPrior","MagicReconstruct"],
        "MagicImpute" => ["HMM","Pedigrees","MagicBase","MagicPrior","MagicReconstruct"],    
        "MagicCall" => ["HMM","Pedigrees","MagicBase","MagicPrior","MagicReconstruct"],
        "SpectralEmbedding" => [],
        "MagicLD" => ["Pedigrees","MagicBase"], 
        "MagicLinkage" => ["HMM","Pedigrees","MagicBase","MagicLD","MagicPrior",
            "MagicReconstruct","MagicImpute"],         
        "MagicMap" => ["HMM", "Pedigrees","MagicBase","MagicPrior","MagicReconstruct",
            "MagicImpute","MagicLD","MagicLinkage", "SpectralEmbedding"],
        "MagicScan" => ["Pedigrees","MagicBase"],        
        "RABBITCLI" => ["HMM", "Pedigrees","MagicBase","MagicSimulate", "MagicPrior","MagicReconstruct","MagicFilter", 
            "MagicCall", "MagicImpute","MagicLD","MagicLinkage", "SpectralEmbedding","MagicMap", "MagicScan"]
    ]    
    filedir = abspath(dirname(@__FILE__))
    @time for (pkg,deps) in pkgdeps
        println("-------------update dependencies of ",pkg, "-------------")        
        Pkg.activate(abspath(filedir, "packages",pkg))                       
        for i in deps
            try 
                Pkg.rm(i)   
            catch err            
                @error err
            end
        end   
        for i in deps
            @info string("-------------add ",i, " for ", pkg, "----------------")
            Pkg.add(PackageSpec(path=filedir,subdir=string("packages","/",i)))
        end        
        Pkg.update()                 
    end
    Pkg.activate()         
    0
catch err
    # @error ]err
    rethrow(err)
    -1
finally
    Pkg.activate()
end

