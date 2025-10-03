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
    s.description = "Reduce ancestry results on a subset of markers such that inter-marker distances <= thincm. "
    workdir = pwd()
    @add_arg_table! s begin
        "--ancestryfile", "-g"
        help = "anncestry file resulting from magicreconstruct"
        arg_type = AbstractString
        required = true        
        "--thincm"
        help = "keep ancestry results on a subseut of markers such that inter-marker distances <= thincm. By default, thincm=0, keeping only the first of markers at the same position. "
        arg_type = Float64
        default = 0.0
        "--outstem", "-o"
        help = "outstem of output files"
        arg_type = AbstractString
        default = "outstem"
        "--workdir", "-w"
        help = "directory for reading and writing files"
        arg_type = AbstractString
        default = workdir
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

function main(args::Vector{String})
    parsed_args = parse_commandline()
    println("Parsed arguments:")
    for (arg, val) in parsed_args
        println(arg, " => ", val)
    end
    ancestryfile = strip(parsed_args[:ancestryfile])
    delete!(parsed_args, :ancestryfile)        
    @time MagicBase.thinmagicancestry(ancestryfile; parsed_args...)
    return 0
end

main(ARGS)
