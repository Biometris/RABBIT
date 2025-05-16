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
    s.description = "Generate RABBIT-format pedigree file"
    workdir = pwd()
    @add_arg_table! s begin
        "--designcodes"
        help = "design codes for each subpopulation"        
        "--founders"
        help = "founders for each subpopulation"        
        "--subpopsizes"
        help = "population sizes for each subpopulation"        
        "--outstem", "-o"
        help = "outstem of output files"        
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

function string2vec(strvec::AbstractString, t::DataType)
    str = strip(strvec)
    @assert str[1] == '[' && str[end] == ']'    
    res = string.(strip.(split(str[2:end-1], ",")))
    t <: AbstractString ? res : parse.(t, res)
end


function main(args::Vector{String})
    parsed_args = parse_commandline()
    println("Parsed arguments:")
    for (arg, val) in parsed_args
        println(arg, " => ", val)
    end
    for keyarg in [:designcodes,:founders,:subpopsizes]
        if isnothing(parsed_args[keyarg])
            delete!(parsed_args, keyarg)
        else
            if keyarg == :subpopsizes
                parsed_args[keyarg] = string2vec(parsed_args[keyarg],Int)
            else
                parsed_args[keyarg] = string2vec(parsed_args[keyarg],String)
            end
        end
    end
    @time generate_magicped(; parsed_args...)
    return 0
end

main(ARGS)
