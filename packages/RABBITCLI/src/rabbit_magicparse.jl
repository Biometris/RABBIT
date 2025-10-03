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

repodir = abspath(joinpath(dirname(@__FILE__), "..",".."))

try
    @eval using $(Symbol("MagicBase"))
    repodir2 = joinpath(pkgdir(MagicBase),"..")    
    if relpath(repodir2, repodir) == "."
        v1, v2 = MagicBase.pkgversion_st(MagicBase),  MagicBase.pkgversion_toml(MagicBase)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err    
    @info string("Install MagicBase and its dependencies from ",repodir)
    for pn in ["Pedigrees","MagicBase"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicBase"))
end


function parse_commandline()
    s = ArgParseSettings()
    s.description = "Parse genofile and  pedfile"
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
        "--isfounderinbred"
        help = "if true, founders are inbred, and otherwise outbred"
        arg_type = Bool
        default = true
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
        "--isbreedped"
        help = "if true, input pedinfo is a breedped-format file instead of magicped-format pedfile for RABBIT. If it is nothing, it is internally set."        
        arg_type = AbstractString
        default = "nothing"
        "--fixed_nself"
        help = "number of selfing generations for the 3rd col = \"fixed\""
        arg_type = Int
        default = 10
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
    reset_kwarg_nothing!(parsed_args,:isbreedped,Bool)
    outstem = parsed_args[:outstem]
    logfile = string(outstem, "_magicparse.log")
    genofile = strip(parsed_args[:genofile])
    pedinfo = parsed_args[:pedinfo]
    delete!(parsed_args, :genofile)
    delete!(parsed_args, :pedinfo)
    push!(parsed_args, :logfile => logfile)
    reset_priority!(parsed_args)
    @time magicparse(genofile, pedinfo; parsed_args...)
    return 0
end

main(ARGS)
