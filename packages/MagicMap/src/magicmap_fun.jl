
"""
    magicmap(genofile, pedinfo; kwargs...)

genetic map construction from genofile and pedinfo.

# Positional arguments

`genofile::AbstractString` genotypic data file.

`pedinfo::Union{MagicBase.JuncDist,AbstractString}` specifies pedigree information via a pedigree fille or a string designcode or via a struct juncdist::JuncDist.

# Keyword arguments

`formatpriority::AbstractVector=["GT","AD"]`: the priority of genotype 
  formats when parasing input vcf genofile.  

`model::AbstractString="jointmodel"`:  prior depedence of ancestral prior process
  along the two homologous chromosomes within an offspring. It must be "depmodel",
  "indepmodel", or "jointmodel". 

`likeparameters::LikeParameters=LikeParameters()`: specifies default genotyping error rates. 

`threshcall::Real = model == "depmodel" ? 0.95 : 0.9`: threshold for genotype calling. The filtering is based on called genotypes. 

`israndallele::Bool=true`: if true, genotyping error model follows the random allelic model, and otherwise the random genotypic model. 

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise they are outbred.

`snpthin::Integer = 1`: take every snpthin-th markers.

`ispermmarker::Bool=true`: if true, permute input marker ordering

`isdupebinning::Union{Nothing,Bool}=false`: if ture, bin duplicate marker. 

`binshare::Real=0.5`: min fraction of shared genotypes between represent marker and each of the rest in a bin. 

`minlodsave::Union{Nothing, Real}=nothing`: results of pairwise analyses are saved only if the LOD score for LD or linkage > minlodsave. 
  If it is nothing, minlodsave increases with number of markers. 

`minldsave::Union{Nothing, Real}=nothing`: results of pairwise LD analyses are saved only if the LD score (squared allelic correlation) > minldsave. 
  If it is nothing, minldsave increases with number of markers. 

`ncluster::Union{Nothing, Integer}=nothing`: number of linkage groups. 
  If it is nothing, ncluster will be inferred in the range of [minncluster, maxncluster]
 
`minncluster::Integer = isnothing(ncluster) ? 1 : ncluster`: min number of linkage groups. 
  If it is nothing, minncluster is set to 1 if ncluster = nothing and otherwise it is set to ncluster

`maxncluster::Integer = isnothing(ncluster) ? 30 : ncluster`: max number of linkage groups. 
  If it is nothing, maxncluster is set to 30 if ncluster = nothing and otherwise it is set to ncluster

`clusteralg::Union{Nothing,AbstractString}=nothing`: clustering algorithm after spectral embedding. 

`minsilhouette::Real=0.0`: delete markers withg silhouette scores < minsilhouette. 

`minlodcluster::Union{Nothing, Real} = nothing`: minimum linkage LOD threshold.
  If it is nothing, estimated internally as the minimum lod keeping the resulting
  graph connected, ignoring the connected components of size < mincomponentsize.

`mincomponentsize::Union{Nothing,Integer} = nothing`: the markers in the graph
  connectecd components of size < mincomponentsize are removed. If it is nothing, it is internally set. 

`maxrf::Union{Nothing,Real} = nothing`: keep pairwise linakge analyses only if recombation fraction <= maxrf. 

`binrf::Union{Nothing,Real}=nothing`: if binrf >= 0, perform linkage-based marker binning such that the recombation fraction 
  for two markers in a bin is always <= binrf, and otherwise not perform linkage-based binning. 
  If it is nothing, binrf is set to 0.001 if #markers > ncluster*2000 and otherwise -1 if ncluster is not nothing, 
  and binrf is set to 0.001 if #markers > maxncluster*2000 and otherwise -1 if ncluster is nothing. 

`alwayskeep::Real=0.99`: neighbors are always kept if its recombation fraction >= alwayskeep, regardless of knncluster or knnorder. 

`minlodcluster::Union{Nothing,Real} = nothing`: min LOD score for clustering. If it is nothing, it is internally set. 

`minlodorder::Union{Nothing,Real} = nothing`: min LOD score for ordering. If it is nothing, it is internally set. 

`maxminlodcluster::Union{Nothing,Real} = nothing,`: if minlodcluster = nothing, minlodcluster is internally estimated with upbound maxminlodcluster. 

`maxminlodorder::Union{Nothing,Real} = nothing`: if minlodorder = nothing, minlodorder is internally estimated with upbound maxminlodorder. 

`knncluster::Union{Nothing,Function} = nothing`: an anonymous function knncluster(x) of #markers x. 
  It returns #nearest neighbors for clustering. If it is nothing, knncluster = x->0.1*x. 

`knnorder::Union{Nothing,Function} = nothing`: an anonymous function knncluster(x) of #markers x in a linkage group. 
  It returns #nearest neighbors for ordering. If it is nothing, knnorder = x->sqrt(x). 

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`commentstring::AbstractString="##"`: rows that begin with commentstring will be ignored.

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.

`outext::AbstractString=".vcf.gz"`: extension of output file for imputed geno.

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicimpute.log"))`:
  log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.

# Example
```julia-repl
julia> magicmap(genofile,pedinfo; ncluster=12)
```
"""
function magicmap(genofile::AbstractString,
    pedinfo::Union{MagicBase.JuncDist,AbstractString};
    formatpriority::AbstractVector=["GT","AD"],    
    model::AbstractString="jointmodel",
    likeparameters::LikeParameters=LikeParameters(),            
    israndallele::Bool=true, 
    threshcall::Real = (model == "depmodel" || !israndallele) ? 0.95 : 0.9, 
    isfounderinbred::Bool = true,        
    snpthin::Integer=1,
    ispermmarker::Bool=true,
    isdupebinning::Union{Nothing,Bool}=false,
    binshare::Real=0.5,    
    byfounder::Integer=0,    
    minlodsave::Union{Nothing, Real}=nothing,
    minldsave::Union{Nothing, Real}=nothing,            
    ncluster::Union{Nothing, Integer}=nothing, 
    minncluster::Integer = isnothing(ncluster) ? 1 : ncluster, 
    maxncluster::Integer = isnothing(ncluster) ? 30 : ncluster,        
    eigselect::AbstractString="eigratio", 
    clusteralg::Union{Nothing,AbstractString}=nothing, 
    minsilhouette::Real=0.0,        
    ncomponent::Union{Nothing,Integer} = nothing,
    mincomponentsize::Union{Nothing,Integer} = nothing,
    maxrf::Union{Nothing,Real} = nothing,
    binrf::Union{Nothing,Real}=nothing, 
    alwayskeep::Real=0.99,            
    maxminlodcluster::Union{Nothing,Real} = nothing,
    maxminlodorder::Union{Nothing,Real} = nothing,
    minlodcluster::Union{Nothing,Real} = nothing,
    minlodorder::Union{Nothing,Real} = nothing,
    knncluster::Union{Nothing,Function} = nothing,
    knnorder::Union{Nothing,Function} = nothing,        
    knnsave::Union{Nothing,Function} = nothing,        
    commentstring::AbstractString="##",
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicmap.log")),
    isparallel::Bool=true,
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicmap"; verbose,delim="=")
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(io,verbose,msg)
        MagicBase.printpkgst(io,verbose,"MagicMap")
    end
    MagicBase.check_infiles(genofile,pedinfo; isbreedped=false, io, commentstring,workdir,verbose)    
    msg = string("list of file options: \n",
            "genofile = ", genofile, "\n",
            "pedinfo = ", pedinfo, "\n",
            "commentstring = ", commentstring, "\n",
            "workdir = ", workdir)
    MagicBase.printconsole(io,verbose,msg)
    outfiles = Vector{String}()
    isnothing(outstem) || (outstem *= "_magicmap")
    # step1 binning
    if isnothing(isdupebinning)
        nmarkers = MagicBase.vcf_count_markers(genofile;commentstring)        
        isdupebinning = nmarkers > 2e4
        printconsole(io,verbose,string("reset isdupebinning=",isdupebinning, " (#markers=",nmarkers, ")"))
    end
    seqerror = MagicBase.get_seqerror(likeparameters)
    if isdupebinning
        binfile = getabsfile(workdir, outstem*"_binning.csv.gz")
        if isfile(binfile)
            msg  = string("skip binning, using existed binfile=",binfile)
            @warn msg
            printconsole(io, false, "Warning: "*msg)        
        else
            startbin = time()
            binfile = first(MagicMap.binning(genofile,pedinfo;
                isdepmodel = model == "depmodel",
                formatpriority, isfounderinbred,
                seqerror, snpthin,binshare,
                isparallel,commentstring, workdir,
                outstem,verbose))
            msg  = string("tused = ", round(time()-startbin,digits=1), "s in binning", 
                ", see logfile = ",outstem*"_binning.log")        
            printconsole(io, verbose,msg)        
        end
    else
        binfile = nothing
    end
    # step2 magicld
    ldfile = getabsfile(workdir, outstem*"_magicld.csv.gz")
    if isfile(ldfile)
        msg  = string("skip magicld, using existed ldfile=",ldfile)
        @warn msg
        printconsole(io, false, "Warning: "*msg)        
    else        
        startld = time()
        if isa(pedinfo, AbstractString)
            if last(splitext(pedinfo))==".csv"
                # pedfile
                magicped = readmagicped(pedinfo; commentstring, workdir)
                skipind = size(magicped.founderinfo,1)
            else
                # designcode
                designinfo = MagicBase.parsedesign(pedinfo)
                skipind = MagicBase.getnfounder(designinfo)
            end
        else
            # juncdist
            skipind = pedinfo.nfounder
        end
        ldfile = magicld(genofile,pedinfo;
            binfile, 
            formatpriority, threshcall,
            isdepmodel = model == "depmodel",
            seqerror, snpthin,minlodsave, minldsave,
            isparallel,commentstring, workdir,
            outstem,verbose
        )
        msg  = string("tused = ", round(time()-startld,digits=1), "s in magicld", 
            ", see logfile = ",outstem*"_magicld.log")        
        printconsole(io, verbose,msg)        
    end
    # step3 magiclinkage
    linkagefile = getabsfile(workdir, outstem*"_magiclinkage.csv.gz")
    if isfile(linkagefile)
        msg  = string("skip magiclinkage, using existed linkagefile=",linkagefile)
        @warn msg
        printconsole(io, false, "Warning: "*msg)   
    else
        startlinage = time()
        linkagefile = magiclinkage(genofile,pedinfo;
            formatpriority, ldfile,  isfounderinbred, 
            model,likeparameters, israndallele, threshcall,
            snpthin,
            byfounder, 
            minlodsave,maxrfsave = 1.0,
            isparallel, commentstring, workdir,
            outstem, verbose
        )
        msg  = string("tused = ", round(time()-startlinage,digits=1), "s in magiclinkage", 
            ", see logfile = ",outstem*"_magiclinkage.log")        
        printconsole(io, verbose,msg)    
    end
    # step4 construct
    if isnothing(maxminlodcluster)
        nmarker = MagicBase.vcf_count_markers(getabsfile(workdir, genofile))
        if nmarker < 2000 || (!isnothing(ncluster) && nmarker < ncluster*200) || (isnothing(ncluster) && nmarker < (minncluster+maxncluster)*100)  
            maxminlodcluster = 5
        else
            magicped = formmagicped(genofile,pedinfo; commentstring, workdir)    
            subpop_weight = get_subpop_weight(magicped)
            maxminlodcluster = round(Int, 10 + 2*log(2,max(1.0,subpop_weight/100)))
        end
        msg = string("reset maxminlodcluster = ", maxminlodcluster)
        printconsole(io, verbose,msg)    
    end    
    mapfile = construct(linkagefile;
        ldfile,ispermmarker,
        ncluster, minncluster, maxncluster, ncomponent, eigselect, minsilhouette, 
        mincomponentsize, maxrf,binrf, alwayskeep, 
        maxminlodcluster, maxminlodorder, minlodcluster, minlodorder,
        knncluster, knnorder, knnsave,
        isparallel,clusteralg,
        workdir, outstem, logfile=io, verbose
    )  
    if !isnothing(outstem)
        plotmarkermap(genofile, mapfile; 
            isphysmap = [true, false], 
            cordigits = isnothing(ncluster) ? 3 : (ncluster<=12 ? 3 : 2), 
            maplabels = ["Physical position(Mbp)", "MagicMap position(cM)"], 
            workdir, io, outstem = outstem*"_construct_compare_physmap", 
            verbose, 
        )        
        plotmarkermap(genofile, mapfile; 
            isphysmap = [false,false],
            cordigits = isnothing(ncluster) ? 3 : (ncluster<=12 ? 3 : 2), 
            maplabels = ["Genetic position(cM)", "MagicMap position(cM)"], 
            workdir, io, outstem = outstem*"_construct_compare_inputmap", 
            verbose, 
        )        
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicmap"; verbose,delim="=")
    outfiles
end

function expandmap(mapfile::AbstractString,
  binfile::Union{Nothing,AbstractString};
  commentstring::AbstractString="##",
  delim::AbstractChar=',',
  workdir::AbstractString = pwd(),
  outstem::AbstractString = "outstem")
  if isnothing(binfile)
      mapfile2 = getabsfile(workdir, mapfile)
      outfile =string(outstem,".csv")
      outfile2 = getabsfile(workdir,outfile)
      if isfile(mapfile2)
          cp(mapfile2,outfile2; force=true)
          return outfile
      else
          @error string(mapfile2, " does not exist")
      end
  end
  mapdf = CSV.read(getabsfile(workdir,mapfile),DataFrame;
      delim,comment=commentstring,missingstring="NA")
  bindf = CSV.read(getabsfile(workdir,binfile),DataFrame;
      delim,comment=commentstring)
  # calculate map from all markers to integers
  gbindf = groupby(bindf,:binno)
  binrule = Dict([begin
      b = gbindf[i][!,:represent] .!= 0
      snps = Vector(gbindf[i][!,:marker])
      only(snps[b]) => snps
  end for i in 1:length(gbindf)])
  issetequal(keys(binrule), mapdf[!,:marker]) || @error "inconsistent markers"
  # expand map
  if :neighbor in propertynames(mapdf)
      nbrrule = Dict(mapdf[!,:marker] .=> [ismissing(i) ? i : split(i,"|") for i in mapdf[!,:neighbor]])
      res = [begin
          snp = mapdf[i,:marker]
          bin = binrule[snp]
          df = reduce(vcat,[DataFrame(mapdf[i,:]) for j in 1:length(bin)])
          df[!,:marker] .= bin
          neighbors = nbrrule[snp]
          for j in 1:size(df,1)
              ismissing(neighbors) && continue
              if df[j,:marker] == snp
                  df[j,:neighbor] = join(neighbors,"|")
              else
                  df[j,:neighbor] = join([rand(binrule[i]) for i in neighbors],"|")
              end
          end
          df
      end for i in 1:size(mapdf,1)]
  else
      res = [begin
          snp = mapdf[i,:marker]
          bin = binrule[snp]
          df = reduce(vcat,[DataFrame(mapdf[i,:]) for j in 1:length(bin)])
          df[!,:marker] .= bin
          df
      end for i in 1:size(mapdf,1)]
  end
  resdf = reduce(vcat, res)
  outfile =string(outstem,".csv")
  outfile2 = getabsfile(workdir,outfile)
  CSV.write(outfile2, resdf; delim=",", header=true, missingstring="NA")
  outfile
end

function get_subpop_weight(magicped::MagicPed)
    offdict = MagicBase.get_subpop2offspring(magicped)
    fdict = MagicBase.get_subpop2founder(magicped)
    subpop_sizels = []
    for popid  in keys(fdict)
        nf = length(fdict[popid])
        noff = length(offdict[popid])
        push!(subpop_sizels, [nf, noff])
    end
    nfls = first.(subpop_sizels)
    noffls = last.(subpop_sizels)
    sum((noffls ./ nfls) .* (noffls ./ sum(noffls)))
end

