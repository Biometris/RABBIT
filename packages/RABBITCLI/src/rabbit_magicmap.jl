function tryusing(pkgname::AbstractString)
    try
        # @eval using Pkg
        # Pkg.update(pkgname)
        @eval using $(Symbol(pkgname))
    catch
        @eval using Pkg
        Pkg.add(pkgname)
        @eval using $(Symbol(pkgname))
    end
end

tryusing("ArgParse")

repodir = abspath(joinpath(dirname(@__FILE__), "..",".."))

try
    @eval using $(Symbol("MagicMap"))
    repodir2 = joinpath(pkgdir(MagicMap),"..")    
    if relpath(repodir2, repodir) == "."
        @eval using $(Symbol("MagicBase"))    
        v1, v2 = MagicBase.pkgversion_st(MagicMap),  MagicBase.pkgversion_toml(MagicMap)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err
    @eval using Pkg
    @info string("Install MagicMap and its dependencies from ",repodir)
    for pn in ["HMM", "Pedigrees","MagicBase","MagicPrior","MagicReconstruct",
        "MagicImpute","MagicLD","MagicLinkage", "SpectralEmbedding"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicMap"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Genetic map construction in connected multiparental populations"
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
        "--ncluster"
        help = "number of linkage groups. If it is nothing, ncluster will be inferred in the range of [minncluster, maxncluster]"
        arg_type = AbstractString
        default = "nothing"
        "--minncluster"
        help = "min number of linkage groups. If it is nothing, minncluster is set to 1 if ncluster = nothing and otherwise it is set to ncluster"
        arg_type = AbstractString
        default = "nothing"
        "--maxncluster"
        help = "max number of linkage groups. If it is nothing, maxncluster is set to 30 if ncluster = nothing and otherwise it is set to ncluster"
        arg_type = AbstractString
        default = "nothing"
        "--minsilhouette"
        help = "delete markers withg silhouette scores < minsilhouette"
        arg_type = Float64
        default = 0.0
        "--formatpriority"
        help = "priorities of genotype formats in a decreasing order"
        arg_type = AbstractString
        default = "[GT,AD]"
        "--model"
        help = "\"depmodel\", \"indepmodel\", or \"jointmodel\" specifies prior dependence of ancestral prior process along two homologous chromosomes within an offspring"
        arg_type = AbstractString
        default = "jointmodel"
        "--likeparameters"
        help = "Set error rate values in the genotypic data model. "
        arg_type = AbstractString        
        default = "LikeParameters()"   
        "--isfounderinbred"
        help = "if true, founders are inbred, and otherwise outbred"
        arg_type = Bool
        default = true                        
        "--snpthin"
        help = "thin markers by taking every snpthin-th markers"
        arg_type = Int
        default = 1
        "--ispermmarker"
        help = "if true, permute input marker ordering"
        arg_type = Bool
        default = true
        "--isdupebinning"
        help = "if ture, bin duplicate marker"
        arg_type = AbstractString
        default = "false"        
        "--mincomponentsize"
        help = "connectecd components of size < mincomponentsize are removed. If it is nothing, it is internally set. "
        arg_type = AbstractString
        default = "nothing"
        "--minlodcluster"
        help = "minimum lod score for clustering. If it is nothing, it is internally set. "
        arg_type = AbstractString
        default = "nothing"
        "--minlodorder"
        help = "minimum lod score for ordering. If it is nothing, it is internally set. "
        arg_type = AbstractString
        default = "nothing"
        "--knncluster"
        help = "number of nearest neighbors for clustering. If -1, it is set to the nearest integer of 0.1*#markers"
        arg_type = Int
        default = -1
        "--knnorder"
        help = "number of nearest neighbors for ordering. If -1, it is set to the nearest integer of sqrt(#markers)"
        arg_type = Int
        default = -1
        "--nworker"
        help = "number of parallel workers for computing"
        arg_type = Int
        default = 1
        "--commentstring"
        help = "rows that begin with commentstring will be ignored"
        arg_type = AbstractString
        default = "##"
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

function reset_priority!(parsed_args)
    # formatpriority
    priority0 = parsed_args[:formatpriority]
    if  occursin("[", priority0) && occursin("]", priority0)
        str = strip(priority0)
        @assert str[1] == '[' && str[end] == ']'
        parsed_args[:formatpriority] = string.(strip.(split(str[2:end-1], ",")))
    else
        @error string("unknown formatpriority: ",formatpriority)
    end
    # snpthin exists in magicmap
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
    logfile = string(outstem, "_magicmap.log")
    genofile = strip(parsed_args[:genofile])
    pedinfo = strip(parsed_args[:pedinfo])
    delete!(parsed_args, :genofile)
    delete!(parsed_args, :pedinfo)
    push!(parsed_args, :logfile => logfile)
    reset_priority!(parsed_args)            
    reset_kwarg_nothing!(parsed_args,:isdupebinning,Bool)
    reset_kwarg_nothing!(parsed_args,:mincomponentsize,Int)    
    reset_kwarg_nothing!(parsed_args,:minlodcluster,Float64)    
    reset_kwarg_nothing!(parsed_args,:minlodorder,Float64)    
    # ncluster, minncluster, maxncluster
    ncluster = tryparse(Int, parsed_args[:ncluster]) 
    minncluster = tryparse(Int, parsed_args[:minncluster])
    if isnothing(minncluster) 
        minncluster = isnothing(ncluster) ? 1 : ncluster    
    end    
    maxncluster = tryparse(Int, parsed_args[:maxncluster])
    if isnothing(maxncluster) 
        maxncluster = isnothing(ncluster) ? 30 : ncluster    
    end
    delete!(parsed_args, :ncluster)
    delete!(parsed_args, :minncluster)
    delete!(parsed_args, :maxncluster)
    # reset knn
    knnccluster0 = parsed_args[:knncluster]
    knnorder0 = parsed_args[:knnorder]
    knncluster = knnccluster0 <= 0 ? nothing : (x->x = knnccluster0)
    knnorder = knnorder0 <= 0 ? nothing : (x->x = knnorder0)
    delete!(parsed_args, :knncluster)
    delete!(parsed_args, :knnorder)
    # setup parallel
    nworker = parsed_args[:nworker]
    isparallel = nworker <= 1 ? false : true
    if isparallel
        tryusing("Distributed")
        nprocs() < nworker+1 && addprocs(nworker+1-nprocs())
        @info string("nworker = ", nworkers())
        @eval @everywhere using MagicMap
    end
    delete!(parsed_args, :nworker)
    push!(parsed_args, :isparallel => isparallel)    
    like = parsed_args[:likeparameters]
    likeparameters = MagicMap.parse_likeparameters(like)
    delete!(parsed_args, :likeparameters)
    @time magicmap(genofile, pedinfo; likeparameters,
        ncluster, minncluster, maxncluster, knncluster, knnorder, parsed_args...)
    if isparallel
        rmprocs(workers()...;waitfor=0)
    end
    return 0
end

main(ARGS)
