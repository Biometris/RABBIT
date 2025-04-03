function tryusing(pkgname::AbstractString)
    try
        @eval using $(Symbol(pkgname))
    catch
        @eval using Pkg
        Pkg.add(pkgname)
        @eval using $(Symbol(pkgname))
    end
end

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
    # rabbit_url = "https://github.com/Biometris/RABBIT2"
    filedir = abspath(dirname(@__FILE__))
    @time for (pkg,deps) in pkgdeps
        println("-------------update dependencies of ",pkg, "-------------")
        Pkg.activate(abspath(filedir, "packages",pkg))        
        try 
            Pkg.rm.(deps)   
        catch err            
            @error err
        end
        for i in deps
            @info string("-------------add ",i, " for ", pkg, "----------------")
            # Pkg.add(PackageSpec(url=rabbit_url,subdir=string("packages/",i)))
            Pkg.add(PackageSpec(path=filedir,subdir=string("packages/",i)))
        end
        Pkg.update()
        Pkg.activate()          
    end
    0
catch err
    rethrow(err)
    -1
finally
    Pkg.activate()
end


