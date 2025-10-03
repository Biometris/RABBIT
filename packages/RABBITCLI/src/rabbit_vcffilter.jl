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
    s.description = "Filter for biallelic markers in vcf file"
    workdir = pwd()
    @add_arg_table! s begin
        "--genofile", "-g"
        help = "vcf genofile"
        arg_type = AbstractString
        required = true
        "--setmarkerid"
        help = "if true, set markerid in format of CHROM_POS or snp{i} with {i} being the marker index if CHROM or POS is missing. If it is nothing, setmarkerid = true only if markerid is missing. "
        arg_type = String
        default = "nothing"
        "--delsamples"
        help = "list of samples to be deleted"
        arg_type = String
        default = "nothing"
        "--deldupe"
        help = "if true, delete sucessive markers that have exactly duplicated genotypes in format of GT"
        arg_type = Bool
        default = false        
        "--isdelmultiallelic"
        help = "delete markers with >2 alleles"
        arg_type = Bool
        default = true
        "--isdelmonomorphic"
        help = "delete markers with only one allele"
        arg_type = Bool
        default = true        
        "--seqstretch"
        help = "delete non-initial markers in a sequence stretch of length <= seqstretch (in bp), assuming marker are ordered by physical positions. If it is not positive, no filtering for short streches."
        arg_type = Int
        default = 0
        "--maxmiss"
        help = "delete markers with missing fraction > maxmiss"
        arg_type = Float64
        default = 0.99
        "--minmaf"
        help = "delete markers with minor allele frequency < minmaf"
        arg_type = Float64
        default = 0.01
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

function main(args::Vector{String})
    parsed_args = parse_commandline()
    println("Parsed arguments:")
    for (arg, val) in parsed_args
        println(arg, " => ", val)
    end
    genofile = strip(parsed_args[:genofile])
    delete!(parsed_args, :genofile)
    reset_kwarg_nothing!(parsed_args, :setmarkerid, Bool)
    a = strip(parsed_args[:delsamples])
    if a == "nothing"
        parsed_args[:delsamples] = nothing
    else
        if  occursin("[", a) && occursin("]", a)
            parsed_args[:delsamples] = strip.(split(a[2:end-1],","))
            @info string("parsed delsamples = ", parsed_args[:delsamples])
        else
            @error string("unknown delsamples: ",delsamples)
        end
    end
    @time vcffilter(genofile; parsed_args...)
    return 0
end

main(ARGS)
