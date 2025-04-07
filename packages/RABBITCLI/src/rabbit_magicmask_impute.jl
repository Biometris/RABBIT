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
    @eval using $(Symbol("MagicImpute"))
    repodir2 = joinpath(pkgdir(MagicImpute),"..")    
    if relpath(repodir2, repodir) == "."
        @eval using $(Symbol("MagicBase")) 
        v1, v2 = MagicBase.pkgversion_st(MagicImpute),  MagicBase.pkgversion_toml(MagicImpute)
        if v1 != v2
            throw(error(string("inconsistent package versions: ", v1, " vs ", v2)))
        end
    else
        throw(error(string("inconsistent package directories: ", repodir, " vs ", repodir2)))
    end
catch err
    @warn err        
    @eval using Pkg
    @info string("Install MagicImpute and its dependencies from ",repodir)
    for pn in ["HMM", "Pedigrees","MagicBase","MagicPrior","MagicReconstruct","MagicImpute"]
        Pkg.develop(path=joinpath(repodir,pn))
    end
    @eval using $(Symbol("MagicBase")) 
    @eval using $(Symbol("MagicImpute"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Genotype masking and imputation in connected multiparental populations"
    workdir = pwd()
    tempdirectory = tempdir()
    byfounder_help = string("alternatively impute founder blocks (a parition of founders). ",
        "If byfounder==-1, impute all founders simulteneously. ",
         "If byfounder==0, partition is based on the founders of each subpopulation. ",
         "If byfounder>0, the size of each block <= byfounder and each block is the subset of a subpopulation's founders.")        
    mapfilehelp = "if it is nothing, use the marker map in the input genofile, and otherwise reset genetic marker map by that in mapfile. "
    mapfilehelp *= "The mapfile can either be in VCF format or in CSV format. For VCF format, genetic map is provided in the \"INFO\" column using keywords \"LINKAGEGROUP\" and \"POSCM\". "
    mapfilehelp *= "For CSV-format, it must contain at least five columns: \"marker\", \"linkagegroup\", \"poscm\", \"physchrom\", and \"physposbp\", where missing values are represented by \"NA\". "
    mapfilehelp *=  "If there exist columns \"binno\" and \"represent\", markers with the same \"binno\" are binned with the represent being the marker with non-zero \"represent\". "
    mapfilehelp *=  "All the rest columns are ignored. "
    phasealg_help = "If phasealg=forwardbackward, the output diplotype probabilities (in format GP), corresonding to the phased genotypes 0|0, 0|1, 1|0, and 1|1, are caculated based on the forward-backward algorithm"
    phasealg_help *= ", and the output phased offspring genotypes (in format GT) are given by those with the largest diplotype probabilities if they are greater than threshcall. "
    phasealg_help *= "If phasealg=viterbi, the output diplotype probabilities (GP) are set to those of phasealg=forwardbackward, and the output phased genotypes (GT) are caculated based on the Viterbi algorithm. "
    phasealg_help *= "If phasealg=unphase, the output genotype probabilities (GP), corresonding to the unphased genotypes 0/0, 0/1, and 1/1, are calculated by transforming the posterior diplotype probabilities of phasealg=forwardbackward"
    phasealg_help *= ", and the output unphased genotypes (GT) are given by those with the largest genotype probabilities if they are greater than threshcall. "
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
        "--foundermask"
        help = "fraction of observed founder genotypes to be masked"
        arg_type = Float64
        default = 0.1
        "--offspringmask"
        help = "fraction of observed offspring genotypes to be masked"
        arg_type = Float64
        default = 0.1
        "--skipmarker"
        help = "skip masking markers with offspring missing fraction >= skipmarker"
        arg_type = Float64
        default = 0.99
        "--minread"
        help = "skip masking genotypes with #reads < minread"
        arg_type = Int
        default = 10
        "--model"
        help = "\"depmodel\", \"indepmodel\", or \"jointmodel\" specifies prior dependence of ancestral prior process along two homologous chromosomes within an offspring"
        arg_type = AbstractString
        default = "jointmodel"
        "--likeparameters"
        help = "parameters for genotypic data model. If isinfererror = true, parameters with values being nothing will be inferred. "
        arg_type = AbstractString
        default = "LikeParameters()"   
        "--threshlikeparameters"
        help = "markers with inferred likeparameters values > threshlikeparameters values will be deleted"
        arg_type = AbstractString
        default = "ThreshLikeParameters()"   
        "--priorlikeparameters"
        help = "priors for likelihood parameters"
        arg_type = AbstractString
        default = "PriorLikeParameters()"   
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
        "--target"
        help = "imputing target: \"all\", \"founder\", or \"offspring\""
        arg_type = AbstractString
        default = "all"
        "--threshimpute"
        help = "impute offspring if maximum posterior probability > threshimpute"
        arg_type = Float64
        default = 0.9
        "--byfounder"
        help = byfounder_help
        arg_type = Int
        default = 0
        "--isallowmissing"
        help = "if true, allow missing in founders during imputation"
        arg_type = Bool
        default = true
        "--isrepeatimpute"
        help = byfounder_help
        arg_type = AbstractString
        default = "nothing"
        "--nrepeatmin"
        help = byfounder_help
        arg_type = Int
        default = 3
        "--nrepeatmax"
        help = byfounder_help
        arg_type = Int
        default = 6
        "--mapfile"
        help = mapfilehelp
        arg_type = AbstractString
        default = "nothing"
        "--iscorrectfounder"
        help = "if true, perform parental error correction. If it is nothing, iscorrectfounder=true if model=depmodel or isinfererror=true or offspring do not have genotypes in AD format"
        arg_type = AbstractString
        default = "nothing"        
        "--phasealg"
        help = phasealg_help
        arg_type = AbstractString
        default = "unphase"        
        "--isdelmarker"
        help = "if true, perform marker deletion"
        arg_type = Bool
        default = true
        # "--delsiglevel"
        # help = "significance level for marker deletion"
        # arg_type = Float64
        # default = 0.01
        "--isinfererror"
        help = "if true, infer marker specific error rates.  If it is nothing, isinfererror=true if model ≠ depmodel or isspacemarker = true. It is necessary to set isinfererror = true for accurately imputing sequence data in heterozygous populations."
        arg_type = AbstractString
        default = "true"   
        "--tukeyfence"
        help = "tukeyfence for detecting outlier error rates"
        arg_type = Float64
        default = 3.0                
        "--isordermarker"
        help = "if true, refine local marker ordering. If it is nothing, isordermarker=true only if mapfile exists."
        arg_type = AbstractString
        default = "nothing"   
        "--inittemperature"
        help = "initial temperature of annealing algorithm for marker ordering. If it is nothing, inittemperature=2.0 if isordermarker and otherwise 0.0"
        arg_type = AbstractString
        default = "nothing"   
        "--coolrate"
        help = "temperature is mutiplied by coolrate after each iteration of annealing agrogrithm"
        arg_type = Float64
        default = 0.85        
        "--isspacemarker"
        help = "if true, refine inter-marker distances. If it is nothing, isspacemarker=true if mapfile exists or isordermarker=true or isphysmap=true."
        arg_type = AbstractString
        default = "nothing"
        "--trimcm"
        help = "remove markers of each segment with distances to the flanking markers > trimcm (cM) "
        arg_type = Float64
        default = 20.0
        "--skeletonsize"
        help = "number of skeleton markers for piecewisely re-scaling inter-marker distances. If it is nothing, skeletonsize is set to the number of marers with distinct positions. "
        arg_type = AbstractString
        default = "nothing"
        "--nworker"
        help = "number of parallel workers for computing among chromosomes"
        arg_type = Int
        default = 1
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
        "--tempdirectory", "-t"
        help = "tempdirectory directory for intermediate results"
        arg_type = AbstractString
        default = tempdirectory
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
    logfile = string(outstem, "_magicimpute.log")
    genofile = strip(parsed_args[:genofile])
    pedinfo = strip(parsed_args[:pedinfo])
    delete!(parsed_args, :genofile)
    delete!(parsed_args, :pedinfo)
    push!(parsed_args, :logfile => logfile)
    mapfile = strip(parsed_args[:mapfile])
    mapfile == "nothing" && (mapfile = nothing)
    delete!(parsed_args, :mapfile)
    reset_priority_subset!(parsed_args)             
    for keyarg in [:iscorrectfounder,:isinfererror]
        reset_kwarg_nothing!(parsed_args,keyarg,Bool) 
        if isnothing(parsed_args[keyarg])
            delete!(parsed_args, keyarg)
        end
    end
    if strip(parsed_args[:skeletonsize]) == "nothing"
        delete!(parsed_args, :skeletonsize)
    else
        reset_kwarg_nothing!(parsed_args,:skeletonsize,Int)    
    end    
    reset_kwarg_nothing!(parsed_args,:isrepeatimpute,Bool)    
    reset_kwarg_nothing!(parsed_args,:isordermarker,Bool)    
    reset_kwarg_nothing!(parsed_args,:isspacemarker,Bool)           
    if isnothing(mapfile)  
        if parsed_args[:isphysmap] == true 
            isnothing(parsed_args[:isspacemarker]) && (parsed_args[:isspacemarker] = true)        
        end 
        for s in [:isordermarker,:isspacemarker]
            isnothing(parsed_args[s]) && (parsed_args[s] = false)
        end
    else        
        for s in [:isordermarker,:isspacemarker]
            isnothing(parsed_args[s]) && (parsed_args[s] = true)
        end
    end
    if parsed_args[:isordermarker] && !parsed_args[:isspacemarker]
        msg = string("isordermarker=",parsed_args[:isordermarker], ", isspacemarker=",parsed_args[:isspacemarker], 
            ". isspacemarker is expected to be true if isordermarker=true!")
        @warn msg
    end
    reset_kwarg_nothing!(parsed_args,:inittemperature,Float64)    
    if isnothing(parsed_args[:inittemperature])
        parsed_args[:inittemperature] = parsed_args[:isordermarker] ? 2.0 : 0.0
    end
    # setup parallel
    nworker = parsed_args[:nworker]
    isparallel = nworker <= 1 ? false : true
    if isparallel
        tryusing("Distributed")
        nprocs() < nworker+1 && addprocs(nworker+1-nprocs())
        @info string("nworker = ", nworkers())
        @eval @everywhere using MagicImpute
    end
    delete!(parsed_args, :nworker)
    push!(parsed_args, :isparallel => isparallel)    
    like = parsed_args[:likeparameters]
    likeparameters = MagicImpute.parse_likeparameters(like)
    delete!(parsed_args, :likeparameters)
    maxlike = parsed_args[:threshlikeparameters]
    threshlikeparameters = MagicImpute.parse_threshlikeparameters(maxlike)
    delete!(parsed_args, :threshlikeparameters)
    priorlike = parsed_args[:priorlikeparameters]
    priorlikeparameters = MagicImpute.parse_priorlikeparameters(priorlike)
    delete!(parsed_args, :priorlikeparameters)
    @time magicmask_impute(genofile, pedinfo; mapfile, likeparameters, threshlikeparameters, priorlikeparameters, parsed_args...)    
    if isparallel
        rmprocs(workers()...;waitfor=0)
    end
    return 0
end

main(ARGS)


