function tryusing(pkgname::AbstractString)
    try 
        @eval using $(Symbol(pkgname))
    catch       
        Pkg.add($pkgname)
        @eval using $(Symbol(pkgname))
    end
end


using Pkg
tryusing("ArgParse")


repodir = abspath(joinpath(dirname(@__FILE__), "..",".."))

try
    @eval using $(Symbol("MagicFilter"))
    repodir2 = joinpath(pkgdir(MagicFilter),"..")    
    if relpath(repodir2, repodir) == "."
        @eval using $(Symbol("MagicBase")) 
        v1, v2 = MagicBase.pkgversion_st(MagicFilter),  MagicBase.pkgversion_toml(MagicFilter)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err    
    @info string("Install MagicFilter and its dependencies from ",repodir)
    for pn in ["HMM", "Pedigrees","MagicBase","MagicPrior","MagicReconstruct","MagicFilter"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicFilter"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Filter markers and individuals in connected multiparental populations"
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
        "--formatpriority"
        help = "priorities of genotype formats in a decreasing order"
        arg_type = AbstractString
        default = "[AD,GT]"
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
        default = "jointmodel"
        "--likeparam"
        help = "Set error rate values in the genotypic data model. "
        arg_type = AbstractString        
        default = "LikeParam()"   
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
        "--isdelmultiallelic"
        help = "delete markers with >2 alleles"
        arg_type = Bool
        default = true        
        "--isdelinconsistent"
        help = "if true, delete markers with inconsistent changes of founder genotypes."
        arg_type = Bool
        default = true
        "--minsubpop"
        help = "delete subpopulaions with size < minsubpop"
        arg_type = Int
        default =1
        "--minnprogeny"
        help = "delete founder and their progeny if the number of progeny < minnprogeny"
        arg_type = Int
        default =1
        "--minmonotest"
        help = "monomorphic test for a subpopulation at a marker is performed if #observed genotypes >= minmonotest and its minor allele frequency >=  minmaf"
        arg_type = Int
        default = 20        
        "--mono2miss"
        help = "if true, offspring genotypes in each monomorphic subpopulation are set to missing, and otherwise only inconsistent offspring genotypes are corrected. And if it is nothing, offspring genotypes are not changed."        
        arg_type = AbstractString
        default = "true"
        "--minmaf"
        help = "test monomorphic for a subpopulation only if its minor allele frequency (maf) < minmaf. And filter for markers with maf >= minmaf."
        arg_type = Float64
        default = 0.05
        "--maxomiss"
        help = "filter for markers with missing fraction in offspring <= maxomiss || missing fraction in founder <  ormaxfmiss"
        arg_type = Float64
        default = 1.0
        "--ormaxfmiss"
        help = "filter for markers with missing fraction in offspring <= maxomiss || missing fraction in founder <  ormaxfmiss"
        arg_type = Float64
        default = 0.0        
        "--offspring_maxmiss"
        help = "delete offspring if its missing fraction > offspring_maxmiss"
        arg_type = Float64
        default = 1.0
        "--isfilterdupe"
        help = "if true, keep only one of duplicated individuals"
        arg_type = Bool
        default = false
        "--offspring_maxcorr"
        help = "two offspring are duplciated if their correlation > offspring_maxcorr"
        arg_type = Float64
        default = 0.99
        "--offspring_cutcorr"
        help = "pairwise offspring correlations are set zeros if they < offspring_cutcorr"
        arg_type = Float64
        default = 0.4        
        "--commentstring"
        help = "rows that begin with commentstring will be ignored"
        arg_type = AbstractString
        default = "##"
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
    logfile = string(outstem, "_magicfilter.log")
    genofile = strip(parsed_args[:genofile])
    pedinfo = strip(parsed_args[:pedinfo])
    delete!(parsed_args, :genofile)
    delete!(parsed_args, :pedinfo)
    push!(parsed_args, :logfile => logfile)
    reset_priority_subset!(parsed_args)    
    reset_kwarg_nothing!(parsed_args,:mono2miss,Bool)       
    # missfilter        
    maxomiss = parsed_args[:maxomiss]
    maxfmiss = parsed_args[:ormaxfmiss]
    missfilter = (maxomiss = maxomiss, ormaxfmiss = maxfmiss)
    delete!(parsed_args, :maxomiss)
    delete!(parsed_args, :ormaxfmiss)
    
    id = "LikeParam"
    id2 = Symbol(lowercase(id))
    likeparam = MagicBase.parse_likeparam(parsed_args[id2],id)    
    delete!(parsed_args, id2)

    @time magicfilter(genofile, pedinfo; likeparam,missfilter, parsed_args...)
    return 0
end

main(ARGS)