function tryusing(pkgname::AbstractString)
    try        
        @eval using $(Symbol(pkgname))
    catch
        Pkg.add(pkgname)
        @eval using $(Symbol(pkgname))
    end
end

using Pkg
tryusing("ArgParse")
tryusing("Distributed")

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
        "--likeparam"
        help = "parameters for genotypic data model. If isinfererror = true, parameters with values being nothing will be inferred. "
        arg_type = AbstractString
        default = "LikeParam(offspringerror=0.005)"           
        "--softthreshlikeparam"
        help = "markers with inferred likeparam values > threshlikeparam values will be deleted only if they are also outliers."
        arg_type = AbstractString
        default = "SoftThreshLikeParam()"
        "--threshlikeparam"
        help = "markers with inferred likeparam values > threshlikeparam values will be deleted"
        arg_type = AbstractString
        default = "ThreshLikeParam()"   
        "--priorlikeparam"
        help = "priors for likelihood parameters"
        arg_type = AbstractString
        default = "PriorLikeParam()"   
        "--tukeyfence"
        help = "tukeyfence for detecting outlier error rates"
        arg_type = Float64
        default = 3   
        "--isfounderinbred"
        help = "if true, founders are inbred, and otherwise outbred"
        arg_type = Bool
        default = true 
        "--threshcall"
        help = "genotypes are called if maximum posterior probability > threshcall"
        arg_type = AbstractString
        default = "nothing"        
        "--iscalloffspring"
        help = "if true, offspring genotypes are called"
        arg_type = Bool
        default = true
        "--israwcall"
        help = "if true, perform raw genotype calling"
        arg_type = Bool
        default = false
        "--minmaf"
        help = "delete markers with minor allele frequency > minmaf"
        arg_type = Float64
        default = 0.05
        "--maxfmiss"
        help = "delete markers with founder genotype missing frequency > maxfmiss"
        arg_type = Float64
        default = 1.0
        "--maxomiss"
        help = "delete markers with offspring genotype missing frequency > maxomiss"
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
        nprocs() < nworker+1 && addprocs(nworker+1-nprocs())
        @info string("nworker = ", nworkers())
        @eval @everywhere using MagicCall
    end
    delete!(parsed_args, :nworker)
    push!(parsed_args, :isparallel => isparallel)    

     likedict = Dict()
    for id in ["LikeParam","SoftThreshLikeParam", "ThreshLikeParam","PriorLikeParam"]
        id2 = Symbol(lowercase(id))
        like = parsed_args[id2]
        likeparam = MagicBase.parse_likeparam(like,id)
        push!(likedict, id2 => likeparam)
        delete!(parsed_args, id2)
    end

    for t in [:threshcall]
        reset_kwarg_nothing!(parsed_args,t,Float64)    
        isnothing(parsed_args[t]) && delete!(parsed_args, t)
    end
    
    reset_kwarg_nothing!(parsed_args,:isinfererror,Bool)    
    if isnothing(parsed_args[:isinfererror]) 
        parsed_args[:isinfererror] = !parsed_args[:israwcall] 
    end
    @time magiccall(genofile, pedinfo; likedict..., parsed_args...)    
    if isparallel
        rmprocs(workers()...;waitfor=0)
    end
    return 0
end

main(ARGS)


