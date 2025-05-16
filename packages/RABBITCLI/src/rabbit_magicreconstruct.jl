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
    @eval using $(Symbol("MagicReconstruct"))
    repodir2 = joinpath(pkgdir(MagicReconstruct),"..")    
    if relpath(repodir2, repodir) == "."
        @eval using $(Symbol("MagicBase"))    
        v1, v2 = MagicBase.pkgversion_st(MagicReconstruct),  MagicBase.pkgversion_toml(MagicReconstruct)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err
    @eval using Pkg
    @info string("Install MagicReconstruct and its dependencies from ",repodir)
    for pn in ["HMM", "Pedigrees","MagicBase","MagicPrior","MagicReconstruct"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicReconstruct"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Haplotype reconstruction in connected multiparental populations"
    workdir = pwd()
    tempdirectory = tempdir()
    @add_arg_table! s begin
        "--genofile", "-g"
        help = "filename for genotypic data file"
        arg_type = AbstractString
        required = true
        "--pedinfo", "-p"
        help = "pedigree information: filename or stringcode"
        arg_type = AbstractString
        required = true
        "--formatpriority"
        help = "priorities of genotype formats in a decreasing order"
        arg_type = AbstractString
        default = "[GP,AD,GT]"
        "--isphysmap"
        help = "if true, transform physical map into genetic map using recomrate, and overwrite the exist genetic map. If false, keep input physical and/or genetic map."
        arg_type = Bool
        default = false
        "--recomrate"
        help = "constant recombation rate in cM/Mbp"
        arg_type = Float64
        default = 1.0
        "--model"
        help = "\"depmodel\", \"indepmodel\", or \"jointmodel\" specifies prior dependence of ancestral prior process along two homologous chromosomes within an offspring"
        arg_type = AbstractString
        required = false
        default = "jointmodel"
        "--likeparameters"
        help = "Set error rate values in the genotypic data model. "
        arg_type = AbstractString
        default = "LikeParameters()"   
        "--isfounderinbred"
        help = "if true, founders are inbred, and otherwise outbred"
        arg_type = Bool
        default = true                
        "--chrsubset"
        help = "subset of chromosomes, with nothing denoting all chromosomes,
        e.g, \"[2,10]\" denotes the second and tenth chromosomes"
        arg_type = AbstractString
        default = "nothing"
        "--markerthin"
        help = "subset of markers by taking every markerthin-th markers"
        arg_type = Int
        default = 1
        "--hmmalg"
        help = "HMM alogrithm must be either forwardbackward or viterbi"
        arg_type = AbstractString
        default = "forwardbackward"
        "--isignorephase"
        help = "if true, the phases of offspring genotypes are ignored."
        arg_type = Bool
        default = false
        "--nworker"
        help = "number of parallel workers for computing among chromosomes"
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
        "--tempdirectory", "-t"
        help = "tempdirectory directory for intermediate results"
        arg_type = AbstractString
        default = tempdirectory
        "--nplot_subpop"
        help = "plots for up to nplot_subpop offspring in each subpopulation"
        arg_type = Int
        default = 10
        "--thincm"
        help = "thin ancestry results so that inter-marker distances > thincm"
        arg_type = Float64
        default = 0.0
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
    parse.(t, split(str[2:end-1], ","))
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


function reset_priority_subset!(parsed_args)
    # formatpriority
    priority0 = parsed_args[:formatpriority]
    if  occursin("[", priority0) && occursin("]", priority0)
        str = strip(priority0)
        @assert str[1] == '[' && str[end] == ']'
        parsed_args[:formatpriority] = string.(strip.(split(str[2:end-1], ",")))
    else
        @error string("unknown formatpriority: ",formatpriority)
    end
    # chrsubset
    if haskey(parsed_args,:chrsubset)
        a = parsed_args[:chrsubset]
        if a == "nothing"
            parsed_args[:chrsubset] = nothing
        else
            if  occursin("[", a) && occursin("]", a)
                parsed_args[:chrsubset] = string2vec(a, Int)
            else
                @error string("unknown chrsubset: ", chrsubset)
            end
        end
    end
    if haskey(parsed_args,:markerthin)
        markerthin = parsed_args[:markerthin]
        # assum the max number of markers in a linkage group < 10^6
        snpsubset= markerthin<=1 ? nothing : 1:markerthin:10^6
        delete!(parsed_args, :markerthin)
        push!(parsed_args, :snpsubset => snpsubset)
    end
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
    logfile = string(outstem, "_magicreconstruct.log")
    genofile = strip(parsed_args[:genofile])
    pedinfo = strip(parsed_args[:pedinfo])
    delete!(parsed_args, :genofile)
    delete!(parsed_args, :pedinfo)
    push!(parsed_args, :logfile => logfile)
    reset_priority_subset!(parsed_args)    
    # set up parallel
    nworker = parsed_args[:nworker]
    isparallel = nworker <= 1 ? false : true
    if isparallel
        tryusing("Distributed")
        if nworkers()>1
            pids = workers()
            @info string("remove current workers ", pids)
            rmprocs(pids...; waitfor=0)
        end
        addprocs(nworker) # add worker processes on local machine
        @info string("nworker = ", nworkers())
        @eval @everywhere using MagicReconstruct
    end
    delete!(parsed_args, :nworker)
    push!(parsed_args, :isparallel => isparallel)
    like = parsed_args[:likeparameters]
    likeparameters = MagicReconstruct.parse_likeparameters(like)
    delete!(parsed_args, :likeparameters)
    @time magicreconstruct(genofile, pedinfo; 
        likeparameters, parsed_args...)
    if isparallel
        rmprocs(workers()...; waitfor=0)
    end
    return 0
end

main(ARGS)
