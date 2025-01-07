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
tryusing("Distributions")

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
    @eval using Pkg
    @info string("Install MagicSimulate and its dependencies from ",repodir)
    for pn in ["Pedigrees","MagicBase","MagicSimulate"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicSimulate"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Simulating genotypic data in connected multiparental populations"
    workdir = pwd()
    @add_arg_table! s begin
        "--fhaplofile", "-g"
        help = "filename for founder haplotypes including marker map"
        arg_type = AbstractString
        required = true
        "--pedinfo", "-p"
        help = "pedigree information: filename or stringcode"
        arg_type = AbstractString
        required = true
        "--popsize"
        help = "poulation size, i.e. the number of offspring"
        arg_type = Int
        default = 200
        "--isfounderinbred"
        help = "if true, founders are inbred, and otherwise outbred"
        arg_type = Bool
        default = true
        "--foundererror"
        help = "distribution of founder allelic error rate among markers"
        arg_type = AbstractString
        default = "Beta(1,199)"
        "--offspringerror"
        help = "distribution of offspring allelic error rate among markers"
        arg_type = AbstractString
        default = "Beta(1,199)"
        "--foundermiss"
        help = "distribution of founder genotype msissing fraction among markers"
        arg_type = AbstractString
        default = "Beta(1,9)"
        "--offspringmiss"
        help = "distribution of offspring genotype msissing fraction among markers"
        arg_type = AbstractString
        default = "Beta(1,9)"
        "--seqfrac"
        help = "fraction of markers being genotyped by sequencing, and the rest by SNP array."
        arg_type = Float64
        default = 0.0
        "--seqdepth"
        help = "distribution of mean read depth among markers"
        arg_type = AbstractString
        default = "Gamma(2,5)"
        "--seqerror"
        help = "distribution of sequence read error probability"
        arg_type = AbstractString
        default = "Beta(1,999)"
        "--allelebalancemean"
        help = "distribution of allelic balance mean among markers"
        arg_type = AbstractString
        default = "Beta(10,10)"
        "--allelebalancedisperse"
        help = "distribution of allele balance overdispersion among markers"
        arg_type = AbstractString
        default = "Exponential(0.05)"
        "--ispheno"
        help = "if true, simulate phenotypes"
        arg_type = Bool
        default = false
        "--pheno_nqtl"
        help = "number of QTLs for simulating phenotypes"
        arg_type = Int
        default = 1        
        "--pheno_h2"
        help = "heritablity for simulating phenotypes"
        arg_type = Float64
        default = 0.5
        # "--select_nqtl"
        # help = "number of QTLs in simulating trait for artifical selection"
        # arg_type = Int
        # default = 1
        # "--select_dom"
        # help = "dominant effect for a selecting trait is given by allelic_effect * select_dom"
        # arg_type = Float64
        # default = 0.0
        # "--select_prop"
        # help = "proportion of zygotes selected in artifical selection. By default, no artifical selection"
        # arg_type = Float64
        # default = 1.0
        "--nplot_subpop"
        help = "number of plots per subpoplation"
        arg_type = Int
        default = 10
        "--outstem", "-o"
        help = "stem of output filenames"
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
    # outstem = parsed_args[:outstem]
    # logfile = string(outstem, ".log")
    fhaplofile = parsed_args[:fhaplofile]
    pedinfo = strip(parsed_args[:pedinfo])
    delete!(parsed_args, :fhaplofile)
    delete!(parsed_args, :pedinfo)
    for k in [:foundermiss, :offspringmiss, :seqdepth, :foundererror,:offspringerror,:seqerror,:allelebalancemean, :allelebalancedisperse]
        parsed_args[k] = eval(Meta.parse(parsed_args[k]))
    end
    @time magicsimulate(fhaplofile, pedinfo; parsed_args...)
    # workdir = parsed_args[:workdir]
    # outfiles = filter(x -> occursin(outstem, x), readdir(workdir))
    # println("output files: ", join(outfiles, ","))
    return 0
end

main(ARGS)
