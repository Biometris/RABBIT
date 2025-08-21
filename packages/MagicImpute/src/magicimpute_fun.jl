
"""
    magicimpute(genofile, pedinfo;
        formatpriority, isphysmap, recomrate, commentstring,kwargs...)

genotype imputation from genofile and pedinfo.

# Positional arguments

`genofile::AbstractString` genotypic data file.

`pedinfo::Union{MagicBase.JuncDist,AbstractString}` specifies pedigree information via a pedigree fille or a string designcode or via a struct juncdist::JuncDist.

# Keyword arguments

See [`formmagicgeno`](@ref) for the arguments (`formatpriority`, `isphysmap`,
`recomrate`,`commentstring`) that are used for formming magicgeno. 
Note that formatpriority=["AD","GT"] by default. 

`mapfile::Union{Nothing, AbstractString}=nothing`: if it is nothing, use the marker map in the input genofile, and otherwise reset genetic marker map by that in mapfile. 
    The mapfile can either be in VCF format or in CSV format. For VCF format, genetic map is provided in the \"INFO\" column using keywords \"LINKAGEGROUP\" and \"POSCM\". 
    For CSV-format, it must contain at least five columns: \"marker\", \"linkagegroup\", \"poscm\", \"physchrom\", and \"physposbp\", where missing values are represented by \"NA\". 
    If there exist columns \"binno\" and \"represent\", markers with the same \"binno\" are binned with the represent being the marker with non-zero \"represent\". 
    All the rest columns are ignored. 

See [`magicimpute!`](@ref) for the other arguments.

# Examples
```julia-repl
julia> magicimpute("geno.vcf.gz","4ril_self3")
```
"""
function magicimpute(genofile::AbstractString,
    pedinfo::Union{MagicBase.JuncDist,AbstractString};
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0, # 1 cM per Mbp    
    model::Union{AbstractString,AbstractVector}="jointmodel",
    mapfile::Union{Nothing, AbstractString}=nothing,   
    likeparam::LikeParam=LikeParam(), 
    softthreshlikeparam::SoftThreshLikeParam=SoftThreshLikeParam(),
    threshlikeparam::ThreshLikeParam=ThreshLikeParam(),
    priorlikeparam::PriorLikeParam=PriorLikeParam(),        
    israndallele::Bool=true, 
    isfounderinbred::Bool=true,        
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    target::AbstractString = "all",        
    threshimpute::Real=0.9,                
    byfounder::Integer=0,
    startbyhalf::Union{Nothing,Integer}=5,     
    isgreedy::Bool=false, 
    threshproposal::Real=0.7,
    isallowmissing::Bool=true,
    isrepeatimpute::Union{Nothing,Bool}=false, 
    nrepeatmin::Integer=3,
    nrepeatmax::Integer=6,         
    isinferjunc::Union{Nothing, Bool} = false,     
    isbinning::Union{Nothing,Bool}=nothing,        
    bincm::Real=0.001, # successive markers with intermarker distance < bincm are binned    
    isinfererror::Union{Nothing, Bool} = true,               
    tukeyfence::Real=2,         
    iscorrectfounder::Union{Nothing, Bool} = true,       
    phasealg::AbstractString="unphase", 
    isdelmarker::Bool= true,         
    delsiglevel::Real = 0.01,    
    isordermarker::Bool = !isnothing(mapfile),        
    isspacemarker::Bool= !isnothing(mapfile) || isordermarker || isphysmap,
    trimcm::Real=20,
	trimfraction::Real=0.025,  #cM    
    skeletonsize::Union{Nothing,Integer} = 100,  
    slidewin_neighbor::Union{Nothing,Integer} = 200,
    slidewin::Union{Nothing,Integer} = nothing,	            
    binriffle::Union{Nothing,Integer} = (!isnothing(mapfile) && isfounderinbred) ? -1 : nothing,  
    orderactions::AbstractVector = ["inverse","permute"],  
    orderactions_neighbor::AbstractVector = ["inverse11","inverse01"],  
    inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.8,
    minaccept::Real=0.15,
    spacebyviterbi::Bool=false,     
    isparallel::Bool=true,    
    isparallelfounder::Bool=true, 
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),
    commentstring::AbstractString="##",
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString= ".vcf.gz",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicimpute.log")),
    verbose::Bool=true,
    more_verbose::Bool=false)
    starttime = time()
    if !isnothing(outstem)         
        # inputoutstem = outstem
        outstem = outstem*"_magicimpute"
    end
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicimpute"; verbose,delim="≡")
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(io,verbose,msg)
        MagicBase.printpkgst(io,verbose,"MagicImpute")
    end
    MagicBase.check_infiles(genofile,pedinfo; isbreedped=false, io, commentstring,workdir,verbose)    
    MagicReconstruct.info_file_arg(genofile, pedinfo, formatpriority,isphysmap, recomrate,
        commentstring,workdir, io, verbose)
    if !isnothing(mapfile)
        isfile(getabsfile(workdir,mapfile)) || error(string("mapfile=",mapfile, " does not exist!"))
    end      
    if !in(phasealg, ["viterbi","forwardbackward"]) && !occursin(r"^unphas",phasealg)
		msg = string("unknown phasealg=",phasealg)
        printconsole(io, verbose, "ERROR: "*msg)
        error(msg)
	end
    MagicReconstruct.check_common_arg(model, 0.005,0.005, 0.001,
        workdir, tempdirectory;chrsubset,snpsubset)
    check_impute_arg(threshimpute,byfounder,
        delsiglevel,trimcm, trimfraction, minaccept,inittemperature, coolrate; target,outext)        
    if isnothing(mapfile)        
        inputneighbor = nothing
        magicgeno=formmagicgeno(genofile,pedinfo; isfounderinbred, 
            formatpriority, isphysmap, recomrate, commentstring, workdir)        
        if bincm > 0
            printconsole(io,verbose,string("bincm=",bincm))
            mapdf = get_binning_mapdf(magicgeno; bincm)
            # isnothing(outstem) || CSV.write(getabsfile(workdir,outstem*"_binning.csv"),mapdf)
            inputbinning = isnothing(mapdf) ? nothing : get_map_binning(mapdf, isbinning; io,verbose)        
        else
            inputbinning = nothing
        end
    else
        genofile2 = MagicBase.resetmap(genofile, mapfile; missingstring = "NA",commentstring, 
            outstem = outstem*"_inputgeno_map", workdir)        
        mapdf = read_construct_map(mapfile; commentstring="##", missingstring="NA",workdir)
        inputneighbor = get_map_neighbor(mapdf)   
        inputbinning = get_map_binning(mapdf, isbinning; io, verbose)        
        magicgeno=formmagicgeno(genofile2,pedinfo; isfounderinbred, 
            formatpriority, isphysmap, recomrate, commentstring, workdir)        
        rm(getabsfile(workdir,genofile2))
    end    
    if isnothing(isbinning) 
        isbinning = !isnothing(inputbinning)
        msg = string("reset isbinning=",isbinning)
        printconsole(io,verbose,msg)
    else
        if isbinning && isnothing(inputbinning)
            isbinning = false
            @warn "reset isbinning=false"
            printconsole(io,false,"Warning: reset isbinning=false")
        else
            printconsole(io,verbose,string("isbinning=",isbinning))
        end
    end        
    if isbinning
        printconsole(io,verbose,string("binriffle=",binriffle))
    end
    magicimpute!(magicgeno;
        target, threshimpute, model, 
        likeparam, softthreshlikeparam, threshlikeparam, priorlikeparam,
        chrsubset, snpsubset,isparallel, isparallelfounder, 
        byfounder,startbyhalf, isgreedy, threshproposal, isallowmissing,
        isrepeatimpute, nrepeatmin, nrepeatmax, 
        isdelmarker, delsiglevel,
        israndallele, isfounderinbred, 
        isinferjunc, iscorrectfounder, phasealg,isinfererror, tukeyfence,          
        inputneighbor, inputbinning, isspacemarker, trimcm, trimfraction, skeletonsize, 
        isordermarker, slidewin_neighbor, slidewin, binriffle, orderactions, orderactions_neighbor, 
        inittemperature, coolrate, minaccept,spacebyviterbi,  
        workdir,tempdirectory,
        outstem,outext,logfile=io, verbose,more_verbose)  
    if !isnothing(outstem) && (isordermarker || isspacemarker)   
        # compare maps        
        mapfiles =[isnothing(i) ? i : getabsfile(workdir,i) for i in  [genofile, mapfile, 
            outstem*"_represent_map.csv",
            outstem*"_represent_enlarged_map.csv",
            outstem*"_map.csv"]]
        maplabels = ["from_genofile(cM)","from_mapfile(cM)","representmap(cM)","represent_enlargedmap(cM)","refinedmap(cM)"]                
        isphysmapls = falses(length(mapfiles))
        markermap1 =  readmarkermap(mapfiles[1]; del_ungrouped = false,
            commentstring,missingstring=["NA","missing"], workdir)        
        if !all(ismissing.(markermap1[!,:physposbp]))
            map1phys  = true
            maplabels[1] = "from_genofile(Mbp)"
        elseif !all(ismissing.(markermap1[!,:poscm]))
            map1phys  = false
        else
            map1phys = nothing
        end        
        if !isnothing(map1phys)
            isphysmapls[1] = map1phys
            b = [isnothing(i) || .!isfile(i) for i in mapfiles]
            deleteat!(mapfiles, b)
            deleteat!(maplabels, b)
            deleteat!(isphysmapls, b)
            if length(mapfiles) > 1
                nchr = length(magicgeno.markermap)
                plotmarkermap(mapfiles...; 
                    maplabels, isphysmap = isphysmapls, 
                    cordigits = nchr <=12 ? 3 : 2,                     
                    outstem= isphysmapls[1] ? outstem*"_compare_physmap" : outstem*"_compare_inputmap",
                    workdir, io, verbose, 
                )
            end
        end      
    end  
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicimpute"; verbose,delim="≡")
    magicgeno
end

"""
    magicimpute!(magicgeno::MagicGeno; kwargs...)

genotype imputation from magicgeno.

# Keyword arguments

`model::Union{AbstractString,AbstractVector}="jointmodel"`:  prior depedence of ancestral prior process
  along the two homologous chromosomes within an offspring. If model is a string, it must be "depmodel",
  "indepmodel", or "jointmodel". If model is a vector, the first element specifies the model 
  for founder imputation and the last element for offspring imputation. 

`likeparam::LikeParam=LikeParam()`: parameters for genotypic data model. 
  If isinfererror = true, parameters with values being nothing will be inferred. 

`softthreshlikeparam::ThreshLikeParam=SoftThreshLikeParam()`: markers with inferred likeparam values > softthreshlikeparam values will be deleted only if they are also outliers. 

`threshlikeparam::ThreshLikeParam=ThreshLikeParam()`: markers with inferred likeparam values > threshlikeparam values will be deleted. 

`priorlikeparam::PriorLikeParam=PriorLikeParam()`: priors for likelihood parameters

`israndallele::Bool=true`: if true, genotyping error model follows the random allelic model, and otherwise the random genotypic model. 

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise they are outbred.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosome indices.
  `nothing` denotes all chromosomes. Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of marker indices within each chromosome.
  `nothing` denotes all markers. Marker indices that are larger than the number of markers
  within the chromosome are deleted.

`target::AbstractString = "all"`: target of imputation. `target=all` imputes
  founders and offspring. `target` must be "all", "founder", "offspring". 

`threshimpute::Real=0.9`: offspring genotypes are imputed if
  the maximum posterior probability > threshimpute.

`byfounder::Integer=0`: alternatively impute each blocks of founders. The founders are partitioned such that the size of each block <= byfounder. 
  If byfounder=-1, impute all founders simulteneously. 
  If byfounder=0, reset to the maximum subpopulation size, and the partition is based on the fouders of each sub-population.

`startbyhalf::Union{Nothing, Integer}=Nothing`: start iteration of obtaining prpoposal founder imputation conditional on the other half of chromosomes. 

`inputneighbor::Union{Nothing,AbstractDict}=nothing`: nearest neighbors for each markers, which is required for neighbor-based marker order refinement. 
  If it is nothing and isordermarker = true, perform only random marker order refinement. 

`inputbinning::Union{Nothing,AbstractDict}=nothing`: a parition of markers into bins.  If it is not nothing, first impute founders for representative markers of each bin 
  and them impute founders for all markers. 

`isinfererror::Bool = true`: if true, infer marker specific likelihood parameters that have values of nothing in likeparam. 

`tukeyfence::Real=2`: tukey fence for detecting outlier error rates (including foundererror, offspringerror, baseerror, and allelicbias). 

`iscorrectfounder::Union{Nothing, Bool} = true`: if true, perform parental error correction.

`phasealg::AbstractString="unphase"`: if phasealg=forwardbackward, the output diplotype probabilities (in format GP), corresonding to the phased genotypes 0|0, 0|1, 1|0, and 1|1, are caculated based on the forward-backward algorithm, 
  and the output phased offspring genotypes (in format GT) are given by those with the largest diplotype probabilities if they are greater than threshcall. 
  If phasealg=viterbi, the output diplotype probabilities (GP) are set to those of phasealg=forwardbackward, and the output phased genotypes (GT) are caculated based on the Viterbi algorithm. 
  If phasealg=unphase, the output genotype probabilities (GP), corresonding to the unphased genotypes 0/0, 0/1, and 1/1, are calculated based on the forward backward algorithm, 
  and the output unphased genotypes (GT) are given by those with the largest genotype probabilities if they are greater than threshcall. 

`isdelmarker::Bool=true`: if true, perform marker deletion.

`delsiglevel::Real = 0.01`: significance level for marker deletion

`isordermarker::Bool=false`: if true, refine local marker ordering.

`isspacemarker::Bool=false`: if true, estimate inter-marker distances.

`trimcm::Real=20`: remove markers of each segment with distances to the flanking markers > trimcm.
  The number of markers of each segment must be less than 5% total number of markers.

`skeletonsize::Union{Nothing,Integer} = 100`: number of skeleton markers for piecewisely re-scaling inter-marker distances. 
  If it is nothing, skeletonsize is set to the number of distint positions in the genetic map before re-scaling. 

`slidewin_neighbor::Union{Nothing,Integer} = 200`: max sliding window size for neighbor-based marker order refinement. 

`slidewin::Union{Nothing,Integer} = nothing`: max sliding window size for random marker order refinement 

`binriffle::Union{Nothing,Integer} = nothing`: valid only in the case of marker binning. Skip magicimpute_founder after replacing representatives with binned markers if binriffle < 0. 
  Keep magicimpute_founder for binned markers but without refinning ordering if 0 <= binriffle <= 1. 
  Keep magicimpute_founder for binned markers if binriffle >=2, and if isordermarker = true set random order refinement with slidewin = binriffle and ignore neighbor-based order refinement. 
    

`orderactions::AbstractVector = ["inverse","permute"]`: update actions for random marker order refinement. 
  It must be a subset of ["permute","inverse","inverse00", "inverse01","inverse10"]. 

`orderactions_neighbor::AbstractVector = ["inverse11","inverse01"]`: update actions for neighbor-based marker order refinement. 

`inittemperature::Real= isordermarker ? 2.0 : 0.0`: initial temperature of annealing algorithm for marker ordering.

`coolrate::Real=0.8`: temperature is mutiplied by coolrate after each iteration of annealing agrogrithm.

`minaccept::Real=0.15`: minimum accept rate for controlling the window size of ordering update.

`isparallel::Bool=true`: if true, parallel multicore computing over chromosomes.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`tempdirectory::AbstractString = tempdir()`: temporary directory for inter-mediate results.

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.

`outext::AbstractString=".vcf.gz"`: extension of output file for imputed geno.

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicimpute.log"))`:
  log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.

# Examples
```julia-repl
julia> magicgeno = formmagicgeno("geno.vcf.gz","ped.csv")
julia> magicimpute!(magicgeno)
```
"""
function magicimpute!(magicgeno::MagicGeno;    
    model::Union{AbstractString,AbstractVector}="jointmodel",
    inputneighbor::Union{Nothing,AbstractDict}=nothing,
    inputbinning::Union{Nothing,AbstractDict}=nothing,            
    likeparam::LikeParam=LikeParam(), 
    softthreshlikeparam::SoftThreshLikeParam=SoftThreshLikeParam(),    
    threshlikeparam::ThreshLikeParam=ThreshLikeParam(),    
    priorlikeparam::PriorLikeParam=PriorLikeParam(),    
    israndallele::Bool=true, 
    isfounderinbred::Bool=true,        
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    target::AbstractString = "all",        
    threshimpute::Real=0.9,        
    byfounder::Integer=0,
    startbyhalf::Union{Nothing,Integer}=5,     
    isgreedy::Bool=false, 
    threshproposal::Real=0.7,
    isallowmissing::Bool=true,
    isrepeatimpute::Union{Nothing,Bool}=false, 
    nrepeatmin::Integer=3,
    nrepeatmax::Integer=6,     
    isinferjunc::Union{Nothing, Bool} = false,     
    isinfererror::Union{Nothing, Bool} = true,                
    tukeyfence::Real=2,        
    iscorrectfounder::Union{Nothing, Bool} = true,     
    phasealg::AbstractString="unphase", 
    isdelmarker::Bool= true,
    delsiglevel::Real = 0.01,                    
    isordermarker::Bool = !isnothing(inputneighbor) ,
    isspacemarker::Bool= !isnothing(inputneighbor) || isordermarker,
    trimcm::Real=20,
	trimfraction::Real=0.025,  #cM    
    skeletonsize::Union{Nothing,Integer} = 100,     
    slidewin_neighbor::Union{Nothing,Integer} = 200,
    slidewin::Union{Nothing,Integer} = nothing,	        
    binriffle::Union{Nothing,Integer} = (!isnothing(inputneighbor) && isfounderinbred) ? -1 : nothing,  
    orderactions::AbstractVector = ["inverse","permute"],  
    orderactions_neighbor::AbstractVector = ["inverse11","inverse01"],  
    inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.8,
    minaccept::Real=0.15,
    spacebyviterbi::Bool=false,     
    isparallel::Bool=true,
    isparallelfounder::Bool=true, 
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString=".vcf.gz",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    verbose::Bool=true,
    more_verbose::Bool=false)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicimpute!"; verbose,delim="=")
    isa(logfile, AbstractString) && printpkgst(io,verbose,"MagicImpute")    
    baseerror = MagicBase.get_likeproperty(likeparam, :baseerror)
	epsf = MagicBase.get_likeproperty(likeparam, :foundererror)
	epso = MagicBase.get_likeproperty(likeparam, :offspringerror)	
    MagicReconstruct.check_common_arg(model, epsf, epso, baseerror,
        workdir, tempdirectory;chrsubset,snpsubset)
    check_impute_arg(threshimpute,byfounder,
        delsiglevel,trimcm, trimfraction,minaccept,inittemperature, coolrate; target,outext)            
    submagicgeno!(magicgeno; chrsubset,snpsubset)    
    msg = string("list of sbuset options: \n",    
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", isnothing(snpsubset) ? "all markers" : snpsubset)
    printconsole(io,verbose,msg)    
    if !isnothing(inputneighbor) && !isnothing(inputbinning)
        issubset(keys(inputbinning),keys(inputneighbor)) || @error "inconsistent inputneighbor and inputbinning"
    end
    MagicBase.check_required_geneticmap(magicgeno; io)
    MagicBase.check_markerorder(magicgeno; io,verbose)    
    model = MagicBase.reset_model(magicgeno.magicped,model;io,verbose)    
    model_founderimpute = isa(model, AbstractVector) ? first(model) : model    
    MagicBase.rawgenoprob!(magicgeno; targets=["founders"] , baseerror, isfounderinbred)
    MagicBase.rawgenocall!(magicgeno; targets=["founders"] , callthreshold = 0.95, isfounderinbred)        
    if target in ["all","founder"]
        # && mean(size.(magicgeno.markermap,1)) > 200      
        isbinning = !isnothing(inputbinning)             
        if isbinning
            if isnothing(binriffle )                
                riffle_prefactor = isnothing(inputneighbor) ? 2 : 1                   
                binriffle = round(Int, riffle_prefactor*mean(length.(values(inputbinning))))                
                binriffle = min(30,max(riffle_prefactor*5, binriffle))         
                printconsole(io,verbose,string("reset binriffle=",binriffle))                         
            end            
            represent_magicgeno, represent_neighbor = get_represent_magicgeno(magicgeno, inputbinning, inputneighbor)            
            # maxiter = ceil(Int,-1*log(2*inittemperature)/log(coolrate))+2 # last temperature ~ 0.5            
            magicimpute_founder!(represent_magicgeno;
                model = model_founderimpute,         
                likeparam, softthreshlikeparam, threshlikeparam, priorlikeparam, 
                israndallele, isfounderinbred,byfounder, startbyhalf, isgreedy, 
                isrepeatimpute, nrepeatmin, nrepeatmax,
                isdelmarker, delsiglevel, iscorrectfounder, threshproposal, isallowmissing,
                isinferjunc, isinfererror, tukeyfence,                  
                inputneighbor=represent_neighbor, 
                isspacemarker, trimcm, trimfraction, skeletonsize, 
                isordermarker, slidewin,slidewin_neighbor, orderactions, orderactions_neighbor, 
                inittemperature, coolrate, minaccept, spacebyviterbi,                
                isparallel=isparallelfounder, logfile=io,workdir,tempdirectory,
                outstem = outstem*"_represent",outext, verbose,more_verbose)
            merge_represent_magicgeno!(magicgeno, represent_magicgeno; inputbinning,isfounderinbred)                                    
            if !isnothing(outstem)
                outputfile= string(outstem,"_represent_enlarged_founder", outext)
                tused = @elapsed savegenodata(outputfile, magicgeno; workdir)
                mem1 = round(Int, memoryuse()/10^6)
                represent_magicgeno = nothing
                represent_neighbor = nothing
                GC.gc()
                GC.gc()
                mem2 = round(Int, memoryuse()/10^6)		
                msg = string("save enlarged genofile: ",outputfile,
                    ", tused=",round(tused,digits=1),"s, mem=",mem1,"|",mem2,"MB")
                printconsole(io,verbose,msg)        
                if isordermarker || isspacemarker || (isnothing(isinfererror) && model != "depmodel") || (!isnothing(isinfererror) && isinfererror)                       
                    outmapfile = savemapfile(magicgeno; workdir, outstem = outstem*"_represent_enlarged")
                    printconsole(io,verbose,string("save enlarged refinedmap file: ", outmapfile))    
                    if binriffle < 0             
                        outmapfile2 = getabsfile(workdir,outstem*"_map.csv")
                        cp(getabsfile(workdir,outmapfile),outmapfile2; force = true)                        
                        printconsole(io,verbose,string("save refinedmap file: ", outmapfile2))    
                    end                    
                end                   
            end            
        end                
        if isbinning && binriffle < 0
            msg = string("skip magicimpute_founder after replacing representatives with binned markers because of binriffle =",binriffle,
                "\n\tTo keep magicimpute_founder for binned markers but without refinning ordering set binriffle = 0 or 1 \n\t",
                "To keep magicimpute_founder for binned markers set binriffle >= 2")
            printconsole(io,verbose,msg)                  
        else
            if isbinning && isordermarker && binriffle <=1 
                msg = string("skip refinning ordering in magicimpute_founder for binned markers for binriffle =",binriffle,
                    "\n\tTo keep refinning ordering in magicimpute_founder for binned markers set binriffle >= 2")
                printconsole(io,verbose,msg)  
            end         
            magicimpute_founder!(magicgeno;
                model = model_founderimpute,
                likeparam, softthreshlikeparam, threshlikeparam,priorlikeparam,
                israndallele, isfounderinbred,byfounder,startbyhalf, isgreedy, 
                isrepeatimpute = isbinning ? false : isrepeatimpute, 
                nrepeatmin, nrepeatmax, 
                isdelmarker, delsiglevel, 
                isinferjunc, iscorrectfounder, threshproposal, isallowmissing,
                isinfererror, tukeyfence,              
                inputneighbor = isbinning ? nothing : inputneighbor,  # if isbinning, neighbor-based order refinement is not performed
                isordermarker = isbinning ? isordermarker && (binriffle > 1) : isordermarker, 
                slidewin = isbinning ? binriffle : slidewin,                 
                inittemperature = isbinning ? 0.0 : inittemperature, 
                coolrate, isspacemarker, trimcm, trimfraction, skeletonsize, 
                slidewin_neighbor, orderactions, orderactions_neighbor, minaccept, spacebyviterbi,                                     
                isparallel=isparallelfounder,logfile=io,workdir,tempdirectory,
                outstem,outext, verbose,more_verbose)
        end     
    end    
    if !isnothing(outstem)
        # plot magicped; magicped will be not changed by magicimpute_offspring
        # outputfile= string(outstem,"_ped.csv")
        # savemagicped(outputfile, magicgeno.magicped; workdir)
        # msg = string("save ", outputfile, " in ",outputfile)
        # printconsole(io,verbose,msg)    
        # TODO: too large pedidgree results in error in plotting
        if size(magicgeno.magicped.founderinfo,1) <= 50
            try 
                fn = string(outstem, "_ped.png")
                gdesign=plotmagicped(magicgeno.magicped;
                    isfounderinbred,outfile = getabsfile(workdir,fn))
                msg = isa(gdesign, AbstractString) ? gdesign : string("save design plot in ", fn)
                MagicBase.printconsole(io,verbose,msg)
            catch
                @warn string("Could not plot magicped")
            end
        end    
    end
    if target in ["all","offspring"]        
        magicimpute_offspring!(magicgeno;
            model = isa(model, AbstractVector) ? last(model) : model, 
            likeparam, threshimpute, 
            israndallele,isfounderinbred,phasealg, 
        	isparallel, workdir,tempdirectory,
            outstem, outext,logfile= io, verbose)
    end    
    MagicBase.info_magicgeno(magicgeno;io,verbose)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicimpute!"; verbose,delim="=")
    magicgeno
end

function check_impute_arg(threshimpute::Real,
    byfounder::Integer,
    delsiglevel::Real,
    trimcm::Real,
    trimfraction::Real, 
    minaccept::Real,
    inittemperature::Real,
    coolrate::Real;
    target::AbstractString="all",
    outext::AbstractString=".vcf.gz")
    if !in(target, ["all","founder","offspring"])
        @error string("target =", target, ", is not in [all, founder, offspring]")
    end
    if !in(outext, [".vcf.gz",".vcf"])
        @error string("outext=",outext, ", not in [.vcf.gz, .vcf]")
    end
    if !(0<=threshimpute<=1)
        @error string("threshimpute=", threshimpute, " is not in [0,1]")
    end
    if byfounder < -1
        @error string("byfounder = ",byfounder, ", is not >= -1")
    end    
    if !(0<delsiglevel<=1)
        @error string("delsiglevel=", delsiglevel, " is not in (0,1]")
    end
    if !(trimcm>0)
        @error string("trimcm=", trimcm, " is not positive")
    end
    if !(0<=trimfraction<1)
        @error string("trimfraction=", trimfraction, " is not in [0,1)")
    end
    if !(0<minaccept<1)
        @error string("minaccept=", minaccept, " is not in (0,1)")
    end
    if !(inittemperature >= 0)
        @error string("inittemperature=", inittemperature, " is not non-negative")
    end
    if !(coolrate >= 0)
        @error string("coolrate=", coolrate, " is not non-negative")
    end
    nothing
end


function read_construct_map(mapfile::AbstractString;
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())
    mapfile2=getabsfile(workdir,mapfile)
    mapdf = CSV.read(mapfile2,DataFrame; delim=',',
        comment=commentstring, missingstring)
    size(mapdf,2) <=4  &&  return nothing
    if issubset([:marker,:linkagegroup],propertynames(mapdf))
        filter!(row ->!ismissing(row[:linkagegroup]),mapdf)
        mapdf[!,:marker] = string.(mapdf[!,:marker])
    else
        @warn string("columns marker and chromosome do not eixst")
    end
    mapdf
end



function get_map_neighbor(mapdf::DataFrame)       
    in(:neighbor,propertynames(mapdf)) || return nothing     
    colnbr = :neighbor    
    ismiss = ismissing.(mapdf[!,colnbr])
    nonmiss = .!ismiss
    snpid = mapdf[nonmiss,:marker]
    neighbor = split.(mapdf[nonmiss,colnbr],"||")
    nsnp = length(snpid)
    rule = Dict(snpid .=> 1:nsnp)   
    iils = Vector{Int}()
    jjls = Vector{Int}()
    for i in eachindex(neighbor)
        for k in neighbor[i]
            nbr = get(rule, k,nothing)
            if !isnothing(nbr) 
                push!(iils, i)
                push!(jjls, nbr)        
            end
        end
    end
    iils, jjls = vcat(iils, jjls), vcat(jjls, iils)
    vvls = ones(Int,length(iils))
    adj = sparse(iils,jjls,vvls,nsnp,nsnp)
    neighbor2 = [first(findnz(adj[:,i])) for i in 1:nsnp]
    res = Dict([snpid[i] => snpid[neighbor2[i]] for i in 1:nsnp])
    if any(ismiss)
        push!(res, [i =>String[] for i in mapdf[ismiss,:marker]]...)
    end
    res
end

function get_map_binning(mapdf::DataFrame,isbinning::Union{Nothing,Bool}; io,verbose)
    if issubset([:binno,:represent],propertynames(mapdf))                
        nsnp = size(mapdf,1)
        nbin = sum(mapdf[!,:represent] .> 0)
        avgbinsize = nsnp/nbin        
        msg = string("#marker=",nsnp,",#bin=",nbin, ",<binsize>=",round(avgbinsize,digits=2))
        printconsole(io,verbose,msg)
        if isnothing(isbinning)             
            isbinning = avgbinsize >= 1.5 
        end
        inputbinning = isbinning ? MagicImpute.get_map_binning(mapdf) : nothing                
    else
        inputbinning = nothing 
    end    
    inputbinning
end

function get_map_binning(mapdf::DataFrame)            
    issubset([:marker,:linkagegroup, :binno,:represent],propertynames(mapdf)) || return nothing
    # isphysmiss = any(ismissing.(mapdf[!,:physposbp])) 
    # isphysmiss = isphysmiss || any(ismissing.(mapdf[!,:physchrom]))
    # chrcol = isphysmiss ? :physchrom : :linkagegroup
    chrcol = :linkagegroup
    binningdf = mapdf[!,[:marker,chrcol,:binno,:represent]]
    sort!(binningdf,[chrcol,:binno])
    gdf = groupby(binningdf,[chrcol,:binno])
    inputbinning = Dict{String,Vector{String}}()
    nmulti_rep = 0
    for df in gdf
        bin = df[:,:marker]
        represents = bin[df[!,:represent] .> 0]
        if length(represents) > 1 # occursin when isdupebinning = true and isrfbinning = true     
            nmulti_rep += 1
            push!(inputbinning,first(represents) => bin)            
        else    
            push!(inputbinning,only(represents) => bin)
        end
    end
    nmulti_rep > 0 && @warn string(nmulti_rep, " out of ", length(gdf), " bins with multiple represents!") 
    all(length.(values(inputbinning)) .== 1) && return nothing
    inputbinning
end

function get_binning_mapdf(magicgeno::MagicGeno; bincm=1e-4)
    dfls = [begin
        nmissing = MagicBase.count_missing(magicgeno.offspringgeno[chr],magicgeno.markermap[chr][!,:offspringformat])
        df = magicgeno.markermap[chr][:,1:5]
        insertcols!(df,6, :nmissing=>nmissing, :binno=>0,:represent=>0)
        bin = 1
        df[1,:binno] = bin
        for i in 2:size(df,1)    
            maxd = maximum(df[i,:poscm] .- df[bin:(i-1),:poscm])
            if 0<= maxd <= bincm
                df[i,:binno] = bin
            else
                bin = i        
                df[i,:binno] = bin
            end
        end
        for subdf in groupby(df,:binno)
            if size(subdf,1)  == 1
                subdf[1,:represent] = subdf[1,:binno]
            else
                pos = argmin(subdf[!,:nmissing])
                subdf[pos,:represent] = subdf[1,:binno]
            end
        end
        df
    end for chr in 1:length(magicgeno.markermap)]
    for i in 2:length(dfls)
        baseno = maximum(dfls[i-1][!,:binno])
        dfls[i][!,:binno] .+= baseno
        b = dfls[i][!,:represent] .> 0
        dfls[i][b,:represent] .+= baseno
    end
    reduce(vcat,dfls)
end


function get_represent_magicgeno(magicgeno::MagicGeno,    
    inputbinning::AbstractDict, inputneighbor::Union{Nothing, AbstractDict})
    represent_markers = keys(inputbinning)
    bls = [[in(i, represent_markers) for i in df[!,:marker]] for df in magicgeno.markermap]
    markermap = map((df,b)->df[b,:], magicgeno.markermap, bls)
    foundergeno=map((A,b)->A[b,:],magicgeno.foundergeno,bls)
    offspringgeno=map((A,b)->A[b,:],magicgeno.offspringgeno,bls)
    represent_magicgeno = MagicGeno(deepcopy(magicgeno.magicped),markermap,foundergeno,
        offspringgeno,deepcopy(magicgeno.misc))
    if isnothing(inputneighbor)
        represent_neighbor2 = nothing
    else
        snpids = reduce(vcat,[i[!,:marker] for i in represent_magicgeno.markermap])
        represent_neighbor = filter(x->in(x.first,snpids),inputneighbor)
        represent_neighbor2 = Dict([j=>intersect(represent_neighbor[j],snpids) for j in snpids])
    end
    represent_magicgeno, represent_neighbor2
end

function merge_represent_magicgeno!(magicgeno::MagicGeno, 
    represent_magicgeno::MagicGeno; 
    inputbinning::AbstractDict,
    isfounderinbred::Bool)
    magicgeno.magicped = represent_magicgeno.magicped
    nchr = length(represent_magicgeno.markermap)
    for chr in 1:nchr
        # update ordering and spacing
        represent_map = represent_magicgeno.markermap[chr]
        ordered_bins = [inputbinning[i] for i in represent_map[!,:marker]]
        ordered_posls = reduce(vcat, map((pos,n) -> repeat([pos],n), represent_map[!,:poscm], length.(ordered_bins)))
        ordered_markers = reduce(vcat,ordered_bins)
        nowdict = Dict(magicgeno.markermap[chr][!,:marker] .=> 1:size(magicgeno.markermap[chr],1))
        oo = [get(nowdict,i,nothing) for i in ordered_markers]
        b = .!isnothing.(oo)
        if !all(b)
            oo = oo[b]
            ordered_posls = ordered_posls[b]
            ordered_markers = ordered_markers[b]
        end
        new_map = magicgeno.markermap[chr][oo,:]
        new_map[!, :marker] == ordered_markers || @error "inconsistent markers"
        new_map[!, :poscm] .= ordered_posls
        magicgeno.markermap[chr] = new_map        
        magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][oo,:]
        magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][oo,:]
        # update foundergeno
        newdict = Dict(new_map[!,:marker] .=> 1:size(new_map,1))
        oo = [get(newdict,i,nothing) for i in represent_map[!,:marker]]
        b = .!isnothing.(oo)
        all(b) || @error "unexpected representive markers"
        magicgeno.foundergeno[chr][oo[b],:] .= represent_magicgeno.foundergeno[chr][b,:]
        magicgeno.markermap[chr][oo[b],:founderformat] .= represent_map[b,:founderformat]
        if isfounderinbred
            magicgeno.foundergeno[chr][oo[b],:] .= represent_magicgeno.foundergeno[chr][b,:]
            magicgeno.markermap[chr][oo[b],:founderformat] .= represent_map[b,:founderformat]
        else
            keepphase = false
            if keepphase 
                magicgeno.foundergeno[chr][oo[b],:] .= represent_magicgeno.foundergeno[chr][b,:]
                magicgeno.markermap[chr][oo[b],:founderformat] .= represent_map[b,:founderformat]
            else
                # ignore phase information for imputed founder genotypes at representtive markers            
                magicgeno.foundergeno[chr][oo[b],:] .= join.(sort.(represent_magicgeno.foundergeno[chr][b,:]))
                unique(represent_map[b,:founderformat]) == ["GT_phased"] || @error string("unexepcted founder format = ",unique(represent_map[b,:founderformat]))
                magicgeno.markermap[chr][oo[b],:founderformat] .= "GT_unphased"
            end
        end
        # update marker-specific error rates
        cols = [:foundererror,:offspringerror,:baseerror,:allelicbias,:allelicoverdispersion,:allelicdropout]
        for col in cols
            ty = eltype(represent_map[!,col])
            ty <: Missing && continue
            ls = Vector{Union{Missing,ty}}(magicgeno.markermap[chr][!,col])
            ls[oo[b]] .= represent_map[b,col]
            magicgeno.markermap[chr][!,col] .= ls
        end
    end
    magicgeno
end
