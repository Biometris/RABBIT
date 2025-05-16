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
    s.description = "Parse breed pedigree file"
    workdir = pwd()
    @add_arg_table! s begin
        "--pedfile", "-p"
        help = "breed-format pedigree file with the first three columns: sample, pedcode, nself"
        arg_type = AbstractString
        required = true
        "--fixed_nself"
        help = "number of selfing generations for the 3rd col = \"fixed\""
        arg_type = Int
        default = 10
        "--commentstring"
        help = "rows that begin with commentstring will be ignored"
        arg_type = AbstractString
        default = "##"
        "--delim"
        help = "delimiter separates each chunk of data within line string"
        arg_type = AbstractChar
        default = ','
        "--outfile", "-o"
        help = "output magicped-format pedigree file"
        arg_type = AbstractString
        default = "magicped.csv"
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
    parse.(t, string.(strip.(split(str[2:end-1], ","))))
end

function main(args::Vector{String})
    parsed_args = parse_commandline()
    println("Parsed arguments:")
    for (arg, val) in parsed_args
        println(arg, " => ", val)
    end
    pedfile = parsed_args[:pedfile]
    delete!(parsed_args, :pedfile)
    @time parsebreedped(pedfile; parsed_args...)
    return 0
end

main(ARGS)
