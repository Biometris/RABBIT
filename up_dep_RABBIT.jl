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
        "MagicCall" => ["HMM","Pedigrees","MagicBase","MagicPrior","MagicReconstruct"],
        "MagicFilter" => ["HMM","Pedigrees","MagicBase","MagicPrior","MagicReconstruct"],
        "MagicImpute" => ["HMM","Pedigrees","MagicBase","MagicPrior","MagicReconstruct"],    
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
    rabbit_url = "https://github.com/Biometris/RABBIT"
    filedir = abspath(dirname(@__FILE__),"packages")    
    @time for (pkg,deps) in pkgdeps
        println("Update dependencies of ",pkg)
        repopath = joinpath(filedir,pkg)
        Pkg.activate(repopath)                
        for i in deps
            try 
                Pkg.rm(i)
            catch err
                @warn err
            end
        end        
        for i in deps
            Pkg.add(PackageSpec(url=rabbit_url,subdir=string("packages/",i)))
        end
        Pkg.update()
        Pkg.activate()          
    end
    0
catch err
    @error err
    -1
finally
    Pkg.activate()
end
