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
    @eval using $(Symbol("MagicScan"))
    repodir2 = joinpath(pkgdir(MagicScan),"..")    
    if relpath(repodir2, repodir) == "."
        @eval using $(Symbol("MagicBase"))    
        v1, v2 = MagicBase.pkgversion_st(MagicScan),  MagicBase.pkgversion_toml(MagicScan)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err
    @eval using Pkg
    @info string("Install MagicScan and its dependencies from ",repodir)
    for pn in ["MagicBase"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicScan"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Genomic scan of QTLs in connected multiparental populations"
    workdir = pwd()
    @add_arg_table! s begin
        "--ancestryfile", "-g"
        help = "ancestry file resulting from `magicreconstruct`"
        arg_type = AbstractString
        required = true
        "--phenofile", "-p"
        help = "phenotypic data file"
        arg_type = AbstractString
        required = true
        "--equation", "-e"
        help = "equation from linear model. If it is nothing, last_colname ~ 1"
        arg_type = AbstractString
        required = false
        default = "nothing"
        "--thresholds"
        help = "list of thresholds for QTL detection. e.g, \"[4.0]\"."
        arg_type = String
        default = "nothing"
        "--nperm"
        help = "number of permutations of phenotypes  for calculating unspecified thresholds"
        arg_type = Int
        default = 200
        "--siglevels"
        help = "significance levels for calculating thresholds by permutations,
        e.g, \"[0.05, 0.10]\"."
        arg_type = String
        default = "[0.05]"
        "--islog10p"
        help = "if islog10p = true, profile refers to -log10 P-value, and otherwise LOD score"
        arg_type = Bool
        default = false
        "--missingstring"
        help = "string denotes a missing phenotypic value"
        arg_type = String
        default = "NA"
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

function string2vec(strvec::AbstractString, t::DataType)
    str = strip(strvec)
    @assert str[1] == '[' && str[end] == ']'
    parse.(t, split(str[2:end-1], ","))
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
    logfile = string(outstem, "_magicscan.log")
    ancestryfile = parsed_args[:ancestryfile]
    phenofile = parsed_args[:phenofile]
    delete!(parsed_args, :ancestryfile)
    delete!(parsed_args, :phenofile)
    push!(parsed_args, :logfile => logfile)
    a = strip(parsed_args[:thresholds])
    if a == "nothing"
        parsed_args[:thresholds] = nothing
    else
        if  occursin("[", a) && occursin("]", a)
            parsed_args[:thresholds] = string2vec(a, Float64)
            @info string("parsed thresholds = ", parsed_args[:thresholds])
        else
            @error string("unknown thresholds: ",thresholds)
        end
    end
    a = parsed_args[:siglevels]
    if  occursin("[", a) && occursin("]", a)
        parsed_args[:siglevels] = string2vec(a, Float64)
        @info string("parsed siglevels = ", parsed_args[:siglevels])
    else
        @error string("unknown siglevels: ",siglevels)
    end
    # TODO parse equation
    if parsed_args[:equation] == "nothing"         
        parsed_args[:equation] = nothing
    else
        # https://discourse.julialang.org/t/build-a-formula-from-a-string/55804/4
        parsed_args[:equation] = @eval(@formula($(Meta.parse(parsed_args[:equation]))))
    end    
    missingstring = [string(strip(parsed_args[:missingstring]))]
    delete!(parsed_args, :missingstring)
    @time magicscan(ancestryfile,phenofile; 
        missingstring, parsed_args...)
    return 0
end

main(ARGS)
