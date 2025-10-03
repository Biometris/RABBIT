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
    @eval using $(Symbol("MagicSimulate"))
    repodir2 = joinpath(pkgdir(MagicSimulate),"..")    
    if relpath(repodir2, repodir) == "."
        @eval using $(Symbol("MagicBase"))    
        v1, v2 = MagicBase.pkgversion_st(MagicSimulate),  MagicBase.pkgversion_toml(MagicSimulate)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err    
    @info string("Install MagicSimulate and its dependencies from ",repodir)
    for pn in ["Pedigrees","MagicBase","MagicSimulate"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicSimulate"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Simulating founder haplotypes"
    workdir = pwd()
    @add_arg_table! s begin
        "--nsnp"
        help = "total number of SNPs"
        arg_type = Int
        default = 1000
        "--nparent"
        help = "number of parents"
        arg_type = Int
        default = 4
        "--isfounderinbred"
        help = "if true, founders are inbred, and otherwise outbred"
        arg_type = Bool
        default = true
        "--chrlen"
        help = "list of chromosome lengths (in centi-Morgan),
        e.g, \"[100,100]\" denotes 100 cM for each of two chromosomes"
        arg_type = String
        default = "[100]"
        "--outfile", "-o"
        help = "output filename"
        arg_type = AbstractString
        default = "sim_fhaplo.vcf.gz"
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
    parse.(t, split(str[2:end-1], ","))
end

function main(args::Vector{String})
    parsed_args = parse_commandline()
    println("Parsed arguments:")
    for (arg, val) in parsed_args
        println(arg, " => ", val)
    end
    a = parsed_args[:chrlen]
    if  occursin("[", a) && occursin("]", a)
        parsed_args[:chrlen] = string2vec(a, Int)
    else
        @error string("unknown chrlen: ",chrlen)
    end
    @time simfhaplo(; parsed_args...)
    return 0
end

main(ARGS)
