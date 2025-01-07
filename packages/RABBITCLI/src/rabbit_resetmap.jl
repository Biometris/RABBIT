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
    @eval using Pkg
    @info string("Install MagicBase and its dependencies from ",repodir)
    for pn in ["Pedigrees","MagicBase"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicBase"))
end


function parse_commandline()
    s = ArgParseSettings()
    s.description = "Replace marker map in input genofile with that of mapfile"
    workdir = pwd()
    @add_arg_table! s begin
        "--genofile", "-g"
        help = "vcf genofile"
        arg_type = AbstractString
        required = true        
        "--mapfile"
        help = "file for marker map, it can either be in VCF format or in CSV format. For CSV-format, it must contain at least four columns: marker, chromosome poscm, physposbp. The values for the last three coloumns can be missingstring"
        arg_type = AbstractString
        required = true        
        "--commentstring"
        help = "rows that begin with commentstring will be ignored"
        arg_type = AbstractString
        default = "##"        
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
    genofile = strip(parsed_args[:genofile])
    delete!(parsed_args, :genofile)    
    mapfile = strip(parsed_args[:mapfile])
    delete!(parsed_args, :mapfile)    
    @time MagicBase.resetmap(genofile,mapfile; parsed_args...)
    return 0
end

main(ARGS)
