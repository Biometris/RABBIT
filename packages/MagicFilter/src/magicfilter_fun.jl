
"""
    magicfilter(genofile, pedinfo;
        formatpriority, isphysmap, recomrate, commentstring,kwargs...)

filter mrkers and founders/offspring from genofile and pedinfo.

# Positional arguments

`genofile::AbstractString` genotypic data file.

`pedinfo::Union{MagicBase.JuncDist,AbstractString}` specifies pedigree information via a pedigree fille or a string designcode or via a struct juncdist::JuncDist.

# Keyword arguments

See [`formmagicgeno`](@ref) for the arguments (`formatpriority`, `isphysmap`,
`recomrate`,`commentstring`) that are used for formming magicgeno. 
Note that formatpriority=["AD","GT"] by default. 

See [`magicfilter!`](@ref) for kwargs.

# Examples
```julia-repl
julia> magicfilter("geno.vcf.gz","4ril_self3")
```
"""
function magicfilter(genofile::AbstractString,    
    pedinfo::Union{Integer, MagicBase.JuncDist,AbstractString};
    model::AbstractString="jointmodel",    
    likeparameters::LikeParameters=LikeParameters(),  
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,    
    threshcall::Real = 0.95, 
    delmultiallelic::Bool=true,
    del_inconsistent::Bool = false,    
    min_subpop::Integer = 1, 
    min_nprogeny::Integer = 1,
    snp_minmaf::Real=0.05,            
    snp_monosubpop::Integer = max(20,round(Int,1/snp_minmaf)),	    
    snp_mono2miss::Union{Nothing,Bool} = true,	    
    snp_missfilter::Union{NamedTuple,Function}= (fmiss,omiss)-> omiss <= 1.0 || fmiss < 0.0,
    isfilterdupe::Bool=false,
    offspring_maxmiss::Real = 0.99,
    offspring_maxcorr::Real = 0.99,
    offspring_cutcorr::Real = 0.4,    
    commentstring::AbstractString="##",
    outstem::AbstractString= "outstem",
    outext::AbstractString=in(last(MagicBase.split_allext(genofile)),[".vcf",".vcf.gz"]) ? ".vcf.gz" : ".csv.gz",
    logfile::Union{AbstractString,IO} = outstem*"_magicfilter.log",
    isparallel::Bool=false,
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()    
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magicfilter"; verbose,delim="≡")
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(logio,verbose,msg)
        MagicBase.printpkgst(logio,verbose,"MagicFilter")
    end
    MagicBase.check_infiles(genofile,pedinfo; isbreedped=false,
        io = logio, commentstring,workdir,verbose)    
    MagicReconstruct.info_file_arg(genofile, pedinfo, formatpriority,isphysmap, recomrate,
        commentstring,workdir, logio, verbose)
    check_filter_arg(outext) #TODO: check more args
    if isa(snp_missfilter, NamedTuple)
        msg = string("(fmiss,omiss)-> omiss <= ",snp_missfilter.snp_maxomiss, " || fmiss < ", snp_missfilter.or_snp_maxfmiss)        
        msg = string("snp_missfilter=",snp_missfilter, " is transformed into an anonymous function:\n\t",  msg)
        MagicBase.printconsole(logio,verbose,msg)        
        snp_maxomiss = snp_missfilter.snp_maxomiss
        or_snp_maxfmiss = snp_missfilter.or_snp_maxfmiss        
        snp_missfilter = (f,o)-> o <= snp_maxomiss || f < or_snp_maxfmiss            
    end
    tused = @elapsed magicgeno = formmagicgeno(genofile, pedinfo;
        isfounderinbred, isphysmap, recomrate,delmultiallelic,
        formatpriority, commentstring, missingstring="NA", workdir)
    msg = string("tused=", round(tused,digits=1), " seconds by formmagicgeno")
    MagicBase.printconsole(logio,verbose,msg)    
    magicfilter!(magicgeno; model, likeparameters,
        isfounderinbred,  threshcall,
        chrsubset,snpsubset,      
        min_subpop,min_nprogeny,
        snp_monosubpop, del_inconsistent,snp_mono2miss,
        snp_missfilter,snp_minmaf, 
        offspring_maxmiss,
        isfilterdupe, offspring_maxcorr,offspring_cutcorr,                
        outstem, logfile=logio, isparallel,workdir, verbose)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
    MagicBase.info_missing(magicgeno;io=logio,verbose)
    if !isnothing(outstem)
        outstem *= "_magicfilter"
        # save geno
        outfile = outstem*"_geno"*outext
        savegenodata(outfile,magicgeno; commentstring,workdir)
        msg = string("save after-magicfilter genofile: ", outfile)
        MagicBase.printconsole(logio,verbose,msg)
        # save ped
        outfile = outstem*"_ped.csv"
        savemagicped(outfile, magicgeno.magicped;workdir)
        msg = string("save after-magicfilter pedfile: ", outfile)
        MagicBase.printconsole(logio,verbose,msg)
        if size(magicgeno.magicped.founderinfo,1) <= 50
            # TODO: too large pedidgree results in error in plotting
            try 
                fn = string(outstem, "_ped.png")
                gdesign=plotmagicped(magicgeno.magicped;
                    isfounderinbred,outfile = getabsfile(workdir,fn))
                msg = isa(gdesign, AbstractString) ? gdesign : string("save design plot in ", fn)
                MagicBase.printconsole(logio,verbose,msg)
            catch
                @warn string("Could not plot magicmped")
            end
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magicfilter"; verbose,delim="≡")
    magicgeno
end


"""
    magicfilter!(magicgeno::MagicGeno; kwargs...)

removes bad mrkers and founders/offspring from magicgeno.

# Keyword arguments

`model::AbstractString="jointmodel"`: prior dependence of
ancestral prior process along the two homologous chromosomes within an offspring.
It must be "depmodel", "indepmodel", or "jointmodel". 

`likeparameters::LikeParameters=LikeParameters()`: specifies default genotyping error rates. 

`isfounderinbred::Bool=true`: specifies if fouonders are inbred

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosome indices.
  `nothing` denotes all chromosomes. Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of marker indices within each chromosome.
  `nothing` denotes all markers. Marker indices that are larger than the number of markers
  within the chromosome are deleted.

`threshcall::Real = model == "depmodel" ? 0.95 : 0.9`: threshold for genotype calling. The filtering is based on called genotypes. 
  
`min_subpop::Integer = 1`: delete subpopulations with size < min_subpop.

`min_nprogeny::Integer = 1`: delete founder and their progeny if the number of progeny < min_nprogeny.

`snp_monosubpop::Integer = 20`: a subpopulation is tested for monomorphic at a marker only if the number of observed genotypes >= snp\\_mono\\_subpop.

`snp_mono2miss::Union{Nothing,Bool} = true`: if true, all offspring genotypes in a monomorphic subpopulation are set to missing, 
and otherwise only inconsistent offspring genotypes are corrected. And if nothing, offspring genotypes are not changed.

`del_inconsistent::Bool = false`: if true, delete markers with inconsistent changes of founder genotypes. 

`snp_minmaf::Real = 0.05`: keep only markers if maf >= snp\\_min\\_maf; maf denotes minor allele frequency.

`snp_missfilter::Function=(fmiss,omiss)-> omiss <= 1.0 || fmiss < 0.0`: keep only markers if snp_missfilter(fmiss, omiss);
fmiss denotes missing fraction in founders, and omiss for offspring.

`offspring_maxmiss::Real = 0.99`: delete offspring if its missing > offspring\\_max\\_miss

`isfilterdupe::Bool=false`: if true, remove duplicated offspring by their correlations.

`offspring_maxcorr::Real = 0.99`: two offspring are duplciated if their correlation >= offspring_maxcorr,

`offspring_cutcorr::Real = 0.4`: pairwise offspring correlations that < offspring_cutcorr are set to zeros. 

`isparallel::Bool=true`: if true, pefrom parallel multicore computing.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`outstem::Union{Nothing,AbstractString}="outstem"` specifies the stem of output files.

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicreconstruct.log"))`:
log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.

# Examples
```julia-repl
julia> magicgeno = formmagicgeno("geno.vcf.gz","ped.csv")
julia> magicfilter!(magicgeno)
```
"""
function magicfilter!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    likeparameters::LikeParameters=LikeParameters(),   
    isfounderinbred::Bool=true,        
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,        
    threshcall::Real = 0.95, 
    min_subpop::Integer = 1, 
    min_nprogeny::Integer = 1,    
    snp_monosubpop::Integer = 20,	    
    snp_mono2miss::Union{Nothing,Bool} = true,	 
    del_inconsistent::Bool = false,    
    snp_minmaf::Real=0.05,            
    snp_missfilter::Function=(fmiss,omiss)-> omiss <= 1.0 || fmiss < 0.0,
    isfilterdupe::Bool=false,
    offspring_maxmiss::Real = 0.99,        
    offspring_maxcorr::Real = 0.99,    
    offspring_cutcorr::Real = 0.4,    
    outstem::AbstractString= "outstem",    
    logfile::Union{AbstractString,IO} = outstem*"_magicfilter.log",
    isparallel::Bool=false,
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magicfilter!"; verbose,delim="=")
    model = MagicBase.reset_model(magicgeno.magicped,model;io=logio,verbose)    
    submagicgeno!(magicgeno; chrsubset,snpsubset)    
    epso = MagicBase.get_offspringerror(likeparameters)
    seqerror = MagicBase.get_seqerror(likeparameters)
    offspringformats = MagicBase.get_offspringformats(magicgeno)
    targets = model == "depmodel" ? ["founders","offspring"] : ["founders"] 
    MagicBase.rawgenoprob!(magicgeno; targets, seqerror, isfounderinbred)
    MagicBase.rawgenocall!(magicgeno; targets = ["founders"] , callthreshold = 1.0 - 2*epso, isfounderinbred)    
    if model == "depmodel" 
        MagicBase.rawgenocall!(magicgeno; targets = ["offspring"], callthreshold = threshcall, isfounderinbred)    
    end
    if model == "depmodel"         
        printconsole(logio, verbose, string("offspringformats = ",join(offspringformats,",")))
        if in("AD",offspringformats)
            msg = string("AD-format genotypes were transformed into GT for model = dpemodel")
            @warn msg
            printconsole(logio, false, "Warning: "*msg)
        end
        if in("GP",offspringformats)
            msg = string("GP-format genotypes were transformed into GT for model = dpemodel")
            @warn msg
            printconsole(logio, false, "Warning: "*msg)
        end
    end
    outstem *= "_magicfilter"
    nsubpop = length(unique(magicgeno.magicped.offspringinfo[!,:member]))
    if nsubpop > 1
        filter_subpop_size!(magicgeno;
            model, isfounderinbred,
            min_subpop,
            outstem,logfile=logio,workdir, verbose)
        filter_founder_progeny!(magicgeno;
            model, isfounderinbred,
            min_nprogeny,
            outstem,logfile=logio,workdir, verbose)        
    end
    filter_marker!(magicgeno; model, isfounderinbred, threshcall,
        snp_monosubpop, epso,
        del_inconsistent,snp_mono2miss,
        snp_missfilter, snp_minmaf,
        outstem, logfile=logio, 
        isparallel, workdir, verbose)
    filter_offspring_missing!(magicgeno; model, isfounderinbred, offspring_maxmiss,
        outstem,logfile=logio,workdir, verbose)          
    if isfilterdupe
        # population structure will disappear if running filter_offspring_dupe! after test_monomorphic! with snp_mono2miss = true
        filter_offspring_dupe!(magicgeno;
            model, isfounderinbred,
            offspring_maxcorr, offspring_cutcorr,
            outstem,logfile=logio,
            isparallel,workdir, verbose)
    end     
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magicfilter!"; verbose,delim="=")
    magicgeno
end

function check_filter_arg(    
    outext::AbstractString)
    if !in(outext, [".vcf.gz",".csv.gz", ".vcf",".csv"])
        @error string("outext=",outext, ", not in [.vcf.gz, .csv.gz, .vcf, .csv]")
    end
    nothing
end

function call_offgeno(magicgeno::MagicGeno,model::AbstractString, chr::Integer;
    callthreshold::Real=0.95)
    chrid = string(magicgeno.markermap[chr][1,:linkagegroup])
    isautosome=!(lowercase(chrid) in ["chrx"])
    ishaploidls = MagicReconstruct.precompute_ishaploidls(magicgeno.magicped.offspringinfo,model,isautosome)
    formatls = magicgeno.markermap[chr][!,:offspringformat]
    calledgeno = copy(magicgeno.offspringgeno[chr])
    calledgeno = call_offgeno!(calledgeno,formatls, ishaploidls,callthreshold)
    calledgeno
end

function call_offgeno!(calledgeno::AbstractMatrix, formatls::AbstractVector,
    isoffspringinbred::BitVector,callthreshold::Real)
    formatset = unique(formatls)
    if "GT_unphased" in formatset
        if any(isoffspringinbred)
            b = formatls .== "GT_unphased"
            subgeno =view(calledgeno, b, isoffspringinbred)
            subgeno[subgeno .== "12"] == "NN"
            subgeno[subgeno .== "21"] == "NN"
        end
    end
    if formatset != ["GT_unphased"]
        if any(isoffspringinbred)
            b = [i in ["AD","GP"] for i in formatls]
            subgeno =view(calledgeno, b, isoffspringinbred)
            subgeno .= MagicBase.genocallhaplo(subgeno, formatls[b];
                seqerror = 0.001,callthreshold = 0.95) .^ 2
        end
        if !all(isoffspringinbred)
            b = [i in ["AD","GP"] for i in formatls]
            subgeno =view(calledgeno, b, .!isoffspringinbred)
            subgeno .= MagicBase.genocalldiplo(subgeno, formatls[b];
                seqerror = 0.001,callthreshold)
        end
    end
    calledgeno
end
