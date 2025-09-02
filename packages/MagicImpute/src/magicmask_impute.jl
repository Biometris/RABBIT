
function magicmask_impute(genofile::AbstractString,
    pedinfo::Union{MagicBase.JuncDist,AbstractString};
    foundermask::Real = 0.1,
    offspringmask::Real = 0.1,
    skipmarker::Real = 0.99, #skip masking genotypes if offspring missing fraction >= skipmarker    
    minread::Integer = 10,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0, # 1 cM per Mbp    
    model::Union{AbstractString,AbstractVector}="jointmodel",
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
    threshproposal::Real=0.7,
    isallowmissing::Bool=true,
    isrepeatimpute::Union{Nothing,Bool}=false, 
    nrepeatmin::Integer=3,
    nrepeatmax::Integer=6,     
    isinferjunc::Union{Nothing, Bool} = false, 
    mapfile::Union{Nothing, AbstractString}=nothing,       
    isbinning::Union{Nothing,Bool}=nothing,        
    bincm::Real=0.001, # successive markers with intermarker distance < bincm are binned
    isinfererror::Union{Nothing, Bool} = true,        
    tukeyfence::Real=2,           
    iscorrectfounder::Union{Nothing, Bool} = true,    
    phasealg::AbstractString="unphase",
    isdelmarker::Bool= true,         
    delsiglevel::Real = 0.01,    
    isordermarker::Bool = !isnothing(mapfile),        
    isspacemarker::Bool= !isnothing(mapfile) || isordermarker == true,
    trimcm::Real=20,
	trimfraction::Real=0.025,  #cM            
    skeletonsize::Union{Nothing,Integer} =  nothing,      
    slidewin::Union{Nothing,Integer} = nothing,
	slidewin_neighbor::Union{Nothing,Integer} = 200,
    inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.8,
    minaccept::Real=0.15,
    isparallel::Bool=true,    
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),
    commentstring::AbstractString="##",
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString= ".vcf.gz",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicimpute.log")),
    verbose::Bool=true,
    more_verbose::Bool=false)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicmask_impute"; verbose,delim="≣")
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(io,verbose,msg)
        MagicBase.printpkgst(io,verbose,"MagicBase")
        MagicBase.printpkgst(io,verbose,"MagicImpute")
    end
    magicmask(genofile,pedinfo; 
        isfounderinbred, formatpriority, isphysmap=false,recomrate,  
        foundermask, offspringmask, skipmarker, minread, 
        commentstring,outstem, outext, logfile = io, workdir, verbose
    )
    masked_genofile = string(outstem,"_magicmask_geno", outext)    
    tused = @elapsed  magicimpute(masked_genofile,pedinfo;
        formatpriority, isphysmap,recomrate, commentstring, 
        target, threshimpute, isrepeatimpute,nrepeatmin, nrepeatmax, 
        model,  likeparam, softthreshlikeparam,threshlikeparam, priorlikeparam,
        chrsubset, snpsubset,isparallel, byfounder, threshproposal, isallowmissing,
        isdelmarker, delsiglevel,
        israndallele, isfounderinbred, 
        mapfile, isbinning, bincm, isinferjunc, iscorrectfounder, phasealg,isinfererror,
        tukeyfence, isspacemarker, trimcm, trimfraction, skeletonsize, 
        isordermarker, slidewin, slidewin_neighbor, inittemperature, coolrate, minaccept,
        workdir,tempdirectory,
        outstem,outext,logfile=io, verbose,more_verbose)      
    msg  = string("tused = ", round(tused,digits=1), "s in magicimpute", 
        ", see logfile = ",outstem*"_magicimpute.log")        
    printconsole(io, verbose,msg)    
    reversed_genofile = string(outstem,"_magicmask_reversed", outext)
    imputed_genofile = string(outstem,"_magicimpute_geno",outext)
    imputeaccuracy(reversed_genofile,imputed_genofile, pedinfo;
        alignfounder=true, workdir, io, outstem = outstem*"_magicimpute", verbose)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicmask_impute"; verbose,delim="≣")
    nothing
end
