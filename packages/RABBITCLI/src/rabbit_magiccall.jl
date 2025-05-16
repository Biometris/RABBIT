function tryusing(pkgname::AbstractString)
    try
        # @eval using Pkg
        # Pkg.update(pkgname)
        @eval using $(Symbol(pkgname))
    catch
        @eval using Pkg
        @eval Pkg.add($pkgname)
        @eval using $(Symbol(pkgname))
    end
end


tryusing("ArgParse")

repodir = abspath(joinpath(dirname(@__FILE__), "..",".."))

try
    @eval using $(Symbol("MagicCall"))
    repodir2 = joinpath(pkgdir(MagicCall),"..")    
    if relpath(repodir2, repodir) == "."
        @eval using $(Symbol("MagicBase")) 
        v1, v2 = MagicBase.pkgversion_st(MagicCall),  MagicBase.pkgversion_toml(MagicCall)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err    
    # installfile = joinpath(repodir, "install_RABBIT.jl")
    # isfile(installfile) || @error string(installfile, "does not exist")
    # include(installfile)
    @eval using Pkg
    @info string("Install MagicCall and its dependencies from ",repodir)
    for pn in ["HMM", "Pedigrees","MagicBase","MagicPrior","MagicReconstruct","MagicCall"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicCall"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Single marker genotype call in connected multiparental populations"
    workdir = pwd()    
    @add_arg_table! s begin
        "--genofile", "-g"
        help = "filename for genotypic data file"
        arg_type = AbstractString
        required = true
        "--pedinfo", "-p"
        help = "pedigree information: filename or stringcode"
        arg_type = AbstractString
        required = true        
        "--model"
        help = "\"depmodel\", \"indepmodel\", or \"jointmodel\" specifies prior dependence of ancestral prior process along two homologous chromosomes within an offspring"
        arg_type = AbstractString
        default = "jointmodel"
        "--likeparameters"
        help = "parameters for genotypic data model. If isinfererror = true, parameters with values being nothing will be inferred. "
        arg_type = AbstractString
        default = "LikeParameters(peroffspringerror=0.0)"   
        "--threshlikeparameters"
        help = "markers with inferred likeparameters values > threshlikeparameters values will be deleted"
        arg_type = AbstractString
        default = "ThreshLikeParameters()"   
        "--priorlikeparameters"
        help = "priors for likelihood parameters"
        arg_type = AbstractString
        default = "PriorLikeParameters(offspringerror=Beta(1.05,9),seqerror=Beta(1.05,9))"   
        "--isfounderinbred"
        help = "if true, founders are inbred, and otherwise outbred"
        arg_type = Bool
        default = true                
        "--threshcall"
        help = "genotye call if maximum posterior probability > threshcall"
        arg_type = AbstractString
        default = "nothing"
        "--israwcall"
        help = "if true, perform raw genotype calling"
        arg_type = Bool
        default = false
        "--minmaf"
        help = "delete markers with minor allele frequency > minmaf"
        arg_type = Float64
        default = 0.05
        "--maxmiss"
        help = "delete markers with genotype missing frequency > maxmiss"
        arg_type = Float64
        default = 0.99
        "--isinfererror"
        help = "if true, infer marker specific error rates. If it is nothing, isinfererror = !israwcall"
        arg_type = AbstractString
        default = "nothing"
        "--nworker"
        help = "number of parallel workers for computing among chromosomes"
        arg_type = Int
        default = 1        
        "--outext"
        help = "extension of output genofile"
        arg_type = AbstractString
        default = ".vcf.gz"
        "--outstem", "-o"
        help = "stem of output filenames"
        arg_type = AbstractString
        default = "outstem"
        "--workdir", "-w"
        help = "directory for reading and writing files"
        arg_type = AbstractString
        default = workdir        
        "--verbose", "-v"
        help = "if true, print messages on console"
        arg_type = Bool
        default = true
    end
    # println(s.description)
    # println(usage_string(s))
    dict = parse_args(s, as_symbols = true)
    for (key,val) in dict
        if isa(val,AbstractString)
            dict[key] = convert(String,strip(dict[key]))
        end
    end
    dict
end

function string2vec(strvec::AbstractString, t::DataType)
    str = strip(strvec)
    @assert str[1] == '[' && str[end] == ']'
    parse.(t, string.(strip.(split(str[2:end-1], ","))))
end

function reset_kwarg_nothing!(parsed_args,kwarg::Symbol,type::DataType)
    a = strip(parsed_args[kwarg])
    if type == Bool
        in(a, ["nothing","true","false"]) || @error string(kwarg, "=",a, " is not in [nothing,true,false]")
    end
    a2 = a == "nothing" ? nothing : parse(type,a)
    delete!(parsed_args, kwarg)
    push!(parsed_args, kwarg => a2)
    parsed_args
end


function main(args::Vector{String})
    parsed_args = parse_commandline()
    verbose = parsed_args[:verbose]
    if verbose
        println("Parsed arguments:")
        for (arg, val) in parsed_args
            println(arg, " => ", val)
        end
    end
    outstem = parsed_args[:outstem]
    logfile = string(outstem, "_magiccall.log")
    genofile = strip(parsed_args[:genofile])
    pedinfo = strip(parsed_args[:pedinfo])
    delete!(parsed_args, :genofile)
    delete!(parsed_args, :pedinfo)
    push!(parsed_args, :logfile => logfile)    
    # setup parallel
    nworker = parsed_args[:nworker]
    isparallel = nworker > 1 
    if isparallel
        tryusing("Distributed")
        nprocs() < nworker+1 && addprocs(nworker+1-nprocs())
        @info string("nworker = ", nworkers())
        @eval @everywhere using MagicCall
    end
    delete!(parsed_args, :nworker)
    push!(parsed_args, :isparallel => isparallel)    
    like = parsed_args[:likeparameters]
    likeparameters = MagicCall.parse_likeparameters(like)
    delete!(parsed_args, :likeparameters)
    maxlike = parsed_args[:threshlikeparameters]
    threshlikeparameters = MagicCall.parse_threshlikeparameters(maxlike)
    delete!(parsed_args, :threshlikeparameters)
    priorlike = parsed_args[:priorlikeparameters]
    priorlikeparameters = MagicCall.parse_priorlikeparameters(priorlike)
    delete!(parsed_args, :priorlikeparameters)
    
    reset_kwarg_nothing!(parsed_args,:threshcall,Float64)    
    if isnothing(parsed_args[:threshcall])
        parsed_args[:threshcall] = parsed_args[:model] == "depmodel" ? 0.95 : 0.9
    end
    reset_kwarg_nothing!(parsed_args,:isinfererror,Bool)    
    if isnothing(parsed_args[:isinfererror]) 
        parsed_args[:isinfererror] = !parsed_args[:israwcall] 
    end
    @time magiccall(genofile, pedinfo; likeparameters, threshlikeparameters, priorlikeparameters, parsed_args...)    
    if isparallel
        rmprocs(workers()...;waitfor=0)
    end
    return 0
end

main(ARGS)


