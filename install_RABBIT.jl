function tryusing(pkgname::AbstractString)
    try        
        @eval using $(Symbol(pkgname))
    catch        
        @eval Pkg.add($pkgname)
        @eval using $(Symbol(pkgname))
    end
end

using Pkg

tryusing("ArgParse")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--isdev"
        help = "if true, pakcages will be developped instead of being added"
        arg_type = Bool        
        default = true                
    end
    return parse_args(s)
end

function main()
    println(PROGRAM_FILE)        
    parsed_args = parse_commandline()    
    isdev = parsed_args["isdev"]
    installstr = isdev ? "develop" : "add"
    @info string("Start to ", installstr, "RABBIT packages...")    
    filedir = abspath(dirname(@__FILE__),"packages")    
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
    rabbit_url = "https://github.com/Biometris/RABBIT"
    @time for pkg in pkgls                    
        println("-------------",installstr, " ", pkg, "-------------")
        if isdev            
            Pkg.develop(PackageSpec(path=joinpath(filedir,pkg)))                        
        else
            Pkg.add(PackageSpec(url=rabbit_url,subdir=string("packages/",pkg)))      
        end
    end
    using RABBITCLI
    @info rabbitversion()
    0
end

main()

