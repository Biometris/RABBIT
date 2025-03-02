
function construct(linkagefile::AbstractString;        
    ldfile::Union{Nothing, AbstractString}=nothing,
    ispermmarker::Bool=true,
    ncluster::Union{Nothing, Integer}=nothing,     
    minncluster::Integer = isnothing(ncluster) ? 1 : ncluster, 
    maxncluster::Integer = isnothing(ncluster) ? 30 : ncluster,        
    eigselect::AbstractString="eigratio", 
    ncomponent::Union{Nothing,Integer} = nothing,    
    minsilhouette::Real=0.0,    
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
    clusteralg::Union{Nothing,AbstractString,AbstractVector}=nothing, 
    isparallel::Bool = true,
    workdir::AbstractString=pwd(),
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{AbstractString,IO}= string(outstem,"_construct.log"),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "construct"; verbose,delim="-")   
    msg = string("list of args: \n",
        "linkagefile = ", linkagefile, "\n",
        "ldfile = ", ldfile, "\n",
        "ispermmarker = ", ispermmarker, "\n",
        "ncluster = ", ncluster, "\n",        
        "minncluster = ", minncluster, "\n",
        "maxncluster = ", maxncluster, "\n",
        "eigselect = ", eigselect, "\n",
        "ncomponent = ", ncomponent, "\n",
        "minsilhouette = ", minsilhouette, "\n",
        "mincomponentsize = ", mincomponentsize, "\n",
        "maxrf(scaled from 0 to 1)= ", maxrf, "\n",
        "binrf = ", binrf, "\n",
        "alwayskeep = ", alwayskeep, "\n",           
        "maxminlodcluster = ", maxminlodcluster, "\n",     
        "maxminlodorder = ", maxminlodorder, "\n",     
        "minlodcluster = ", minlodcluster, "\n",
        "minlodorder = ", minlodorder, "\n",        
        "knncluster = ", isnothing(knncluster) ? "nothing" : first(code_lowered(knncluster)), "\n",
        "knnorder = ", isnothing(knnorder) ? "nothing" : first(code_lowered(knnorder)), "\n",
        "knnsave = ", isnothing(knnsave) ? "nothing" : first(code_lowered(knnsave)), "\n",                
        "clusteralg = ", clusteralg, "\n",        
        "isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)
    if !isnothing(ncluster)
        if minncluster > ncluster
            msg = string("minncluster ", minncluster, " > ncluster ", ncluster)
            msg *= string(", reset minncluster=",ncluster)
            printconsole(logio, verbose, msg)
            minncluster = ncluster
        end
        if maxncluster < ncluster
            msg = string("maxncluster ", maxncluster, " < ncluster ", ncluster)
            msg *= string(", reset maxncluster=",ncluster)
            printconsole(logio, verbose, msg)
            maxncluster = ncluster
        end    
    end
    # reading data
    tused = @elapsed inds, markers, physmapdict, nmissingls, recomnonfrac, recomlod = MagicBase.read_linkage(linkagefile; workdir)                
    mem = round(Int, memoryuse()/10^6);
    msg = string("read linkagefile, tused=", round(tused,digits=1),"s, mem=",mem, "MB")
    printconsole(logio, verbose, msg)
    if ispermmarker
        neworder = randperm(length(markers))        
        markers .= markers[neworder]        
        nmissingls .= nmissingls[neworder]
        recomnonfrac .= recomnonfrac[neworder,neworder]
        recomlod .= recomlod[neworder,neworder]        
    end
    if isnothing(ldfile)
        ldlod = nothing      
        dupebindict = nothing        
    else
        tused = @elapsed inds2, markers2, dupebindict, _, ldlod = MagicBase.read_ld(ldfile; workdir)             
        mem = round(Int, memoryuse()/10^6);
        msg = string("read ldfile, tused=", round(tused,digits=1),"s, mem=",mem, "MB")
        printconsole(logio, verbose, msg)
        if !isnothing(dupebindict)
            issetequal(markers2,keys(dupebindict)) ||  @error "inconsistent between markers and dupebin"
        end
        if !issetequal(inds,inds2) 
            msg = string("inconsistent individuals! individuals in linkagefile but not in ldfile: ",setdiff(inds,inds2), 
                "\nindividuals in ldfile but not in linkagefile: ",setdiff(inds2,inds))
            @error msg
            printconsole(logio,verbose,"ERROR: "*msg)
        end
        if !issetequal(markers,markers2) 
            msg = string("inconsistent markers between linkagefile and ldfile")
            @error msg
            printconsole(logio,verbose,"ERROR: "*msg)
        end 
        if markers != markers2 
            if !ispermmarker
                msg = string("inconsistent marker orderings between linkagefile and ldfile")
                @warn msg
                printconsole(logio,verbose,"WARNING: "*msg)
            end
            snprule = Dict(markers2 .=> 1:length(markers2))
            ii = [get(snprule,i,nothing) for i in markers]
            markers2 .= markers2[ii]            
            ldlod .= ldlod[ii,ii]            
        end
        markers == markers2 || @error "inconsistent markers"
    end    
    # return (markers, physposbpls, nmissingls, recomnonfrac, recomlod, ldlod, dupebindict)    
    if isnothing(maxminlodcluster)
        nmarker = length(markers)
        if nmarker < 2000 || (!isnothing(ncluster) && nmarker < ncluster*200) || (isnothing(ncluster) && nmarker < (minncluster+maxncluster)*100)  
            maxminlodcluster = 5
        else
            maxminlodcluster = 10
        end         
        msg = string("reset maxminlodcluster = ", maxminlodcluster)
        printconsole(logio, verbose,msg)    
    end   
    if isnothing(binrf)
        nmarker = length(markers)
        if nmarker < 20000 || (!isnothing(ncluster) && nmarker < ncluster*2000) || (isnothing(ncluster) && nmarker < (minncluster+maxncluster)*1000)  
            binrf = -1.0
        else
            binrf = 1e-3
        end
        msg = string("reset binrf = ", binrf)
        printconsole(logio, verbose,msg)    
    end
    if binrf < 0
        is_cosgregate_binning = false
    else
        is_cosgregate_binning = true
        if isnothing(dupebindict)
            # physmapdict = Dict(markers .=> tuple.(physchromls,physposbpls))
            dupebindict = Dict([i => (i,physmapdict[i]...) for i in convert.(String,markers)])
            isdupebin = false
        else
            isdupebin = true
        end
        nmarkerbef = length(markers)        
        minlod_bin = first(findminlod(recomnonfrac,recomlod,ldlod;
            ncomponent = 1, maxrf=0.7, minminlod = 5.0, maxminlod = 2*maxminlodcluster,
            mincomponentsize=min(20,5+round(Int, length(markers)/5000)), alwayskeep = 1.0))
        # dupebindict is updated to include cosegregation binning
        tused = @elapsed  markers, nmissingls, recomnonfrac,recomlod,ldlod = binning_cosegrate!(markers,nmissingls,recomnonfrac, recomlod, 
            ldlod, dupebindict; binrf,minlod_bin);        
        if isdupebin        
            totmarker = sum(length(split(first(v),"||")) for (k,v) in dupebindict)
            msg = string("#cosegrate_bins = ", length(markers), " for #markers=", totmarker, ", #ld_bins=", nmarkerbef)
        else
            msg = string("#cosegrate_bins = ", length(markers), " for #markers=", nmarkerbef)
        end
        msg *= string(", binrf=",binrf,", minlod_bin=",minlod_bin)
        msg *= string("; tused=",round(tused,digits=1),"s")
        printconsole(logio,verbose,msg)
    end    
    # setup parameters including minlod
    nzlod =  nonzeros(recomlod)
    minlodsave = round(min(unique(nzlod)...),digits=3)    
    msg = string("minlodsave=",minlodsave)
    printconsole(logio,verbose,msg)
    inputmincomponentsize = mincomponentsize
    minminlodcluster = round(Int,minlodsave)
    if isnothing(mincomponentsize)        
        mincomponentsize = min(20,5+round(Int, length(markers)/5000)) 
        minlod = length(markers) < 1e4 ? minminlodcluster : max(5,minminlodcluster)        
        mincomponentsize = findmincomponentsize(mincomponentsize, recomnonfrac,recomlod,ldlod; 
            minlod,maxrf = 1.0, alwayskeep)
        msg = string("set mincomponentsize = ",mincomponentsize)
        printconsole(logio,verbose,msg)
    end         
    if isnothing(maxrf)        
        maxrf = findmaxrf(recomnonfrac,recomlod,ldlod; minlod = minlodsave, 
            mincomponentsize,alwayskeep)
        msg = string("set maxrf=",maxrf, " for #marker=",size(recomnonfrac,1))
        printconsole(logio, verbose,msg)
    end         
    if isnothing(ncomponent)         
        minlod = max(2, minlodsave)           
        nconn0 = first(count_connected_components(recomnonfrac,recomlod,
            ldlod,maxrf, minlod, mincomponentsize, alwayskeep))        
        maxncomponent = min(ceil(Int,minncluster/2),max(2,nconn0))                
        ncomponent = findncomponent(recomnonfrac,recomlod,ldlod; minlodsave, maxminlod = 3, 
            maxncomponent, maxrf, mincomponentsize, alwayskeep)
        msg = string("set ncomponent = ",ncomponent)
        printconsole(logio, verbose, msg)
    end  
    if isnothing(minlodcluster)     
        minlodcluster,nungroup0, _, lodhis = findminlod(recomnonfrac,recomlod,ldlod;
            ncomponent, maxrf, minminlod = minminlodcluster, 
            maxminlod = maxminlodcluster,
            mincomponentsize,alwayskeep)
        msg = string("set minlodcluster = ",minlodcluster)
        msg *= string("; initial #ungrouped=",nungroup0," for minlodsave=",minlodsave, "")
        msg *= string("; minlodcluster_his=", lodhis)
        printconsole(logio,verbose,msg)        
        if isnothing(inputmincomponentsize) && ncomponent>1
            mincomponentsize2 = findmincomponentsize(mincomponentsize, recomnonfrac,recomlod,ldlod; 
                minlod = max(5,minlodcluster), maxrf, alwayskeep)
            if mincomponentsize2 > mincomponentsize
                msg = string("reset mincomponentsize = ",mincomponentsize2)
                printconsole(logio,verbose,msg)
                mincomponentsize = mincomponentsize2
            end       
        end
    end
    if isnothing(maxminlodorder)
        if minlodcluster <= 10
            maxminlodorder = 5
        else
            maxminlodorder = round(Int, minlodcluster/2)            
        end
        msg = string("reset maxminlodorder = ", maxminlodorder)
        printconsole(logio, verbose,msg)    
    end    
    if minlodsave > minlodcluster && minlodcluster>=3
        @warn string("minlodsave is too large! minlodsave(", minlodsave, ") >= ", " minlodcluster(", minlodcluster, ")")
    end
    nconn0, nungroup0, _ = count_connected_components(recomnonfrac,recomlod,ldlod,
        maxrf, minminlodcluster,mincomponentsize, alwayskeep)
    nsnp = length(markers) - nungroup0   
    if isnothing(knncluster)             
        kknncluster = round(Int,nsnp/10)
        kknncluster = min(max(20, kknncluster),nsnp)
        msg = string("set knn(cluster) = ",kknncluster)
        printconsole(logio,verbose,msg)
    else        
        kknncluster = round(Int,knncluster(nsnp))
        kknncluster = min(max(20, kknncluster),nsnp)
    end
    if isnothing(inputmincomponentsize) && ncomponent>1
        mincomponentsize2 = findmincomponentsize(mincomponentsize, recomnonfrac,recomlod,ldlod; 
            minlod = minlodcluster,knn = kknncluster, maxrf, alwayskeep)
        if mincomponentsize2 > mincomponentsize
            msg = string("reset mincomponentsize = ",mincomponentsize2)
            printconsole(logio,verbose,msg)
            mincomponentsize = mincomponentsize2
        end
    end
    # marker grouping        
    if isnothing(clusteralg)         
        clusteralg = (length(markers) < 2e4 && minncluster == maxncluster) ? ["hclust","kmeans"] : "kmeans"               
        msg = string("set clusteralg = ",clusteralg, " for #markers ")
        msg *= string(clusteralg == "kmeans" ? ">=" : "<", "20000")
        printconsole(logio, verbose, msg)
    end      
    eigweightfrac = 0.01
    resgrouping = marker_grouping(recomnonfrac,recomlod, ldlod; 
        ncluster, minncluster,maxncluster,
        maxrf, alwayskeep, minlodcluster, kknncluster, mincomponentsize, minsilhouette, 
        clusteralg, eigselect, eigweightfrac,eigweighttype = "eiggapratio", io=logio,verbose)
    nungroup = length(markers)-sum(length.(resgrouping.clusters))
    msg = string("#ungroupped markers = ", nungroup, " out of ", length(markers))
    printconsole(logio,verbose,msg) 
    if !isnothing(outstem)    
        # eigen
        geigen = SpectralEmbedding.plot_eigen(resgrouping.eigenvals,
            resgrouping.eigenvecs[:,1:length(resgrouping.clusters)], 
            resgrouping.clusters)
        outfile =string(outstem,"_construct_eigen.png")        
        MagicBase.savefig(geigen,getabsfile(workdir,outfile))
        printconsole(logio,verbose,string("eigenplot in ", outfile))
        # silhouette    
        sumsilh = SpectralEmbedding.summary_silhouette(resgrouping.res_silhouettes; eigweightfrac)
        for j in 2:size(sumsilh,2)
            sumsilh[!,j] .= round.(sumsilh[!,j],digits=6)
        end
        if isa(clusteralg, AbstractVector)
            insertcols!(sumsilh,1, :clusteralg => clusteralg)
        end
        CSV.write(logio,sumsilh;header =true, append=true)                
        verbose && @info sumsilh
        outfile = outstem*"_construct_silhouette.csv"
        open(getabsfile(workdir, outfile),"w") do outio
            descripls = ["ncluster, number of clusters", 
                "eigweight, weight based on ncluster-th and (ncluster+1)-th eigenvalues", 
                "eigweight_scaled, eigweight divided by the maximum of all eigweights", 
                "average_silhouette, silhouette score averages over markers && see https://en.wikipedia.org/wiki/Silhouette_(clustering)",                
            ]            
            size(sumsilh,2) == 4 && push!(descripls,"weighted_silhouette, average_silhouette weighted by eigengap")
            if isa(clusteralg, AbstractVector)
                pushfirst!(descripls,"clusteralg, clustering algorithm")
            end
            for i in eachindex(descripls)
                write(outio, string("##col_",i, ", ", descripls[i],"\n"))
            end
            CSV.write(outio,sumsilh; append=true, header=true)
        end
        fig = SpectralEmbedding.plot_silhouette(resgrouping.res_silhouettes;eigweightfrac)
        outfile = outstem*"_construct_silhouette.png"
        MagicBase.savefig(fig,getabsfile(workdir, outfile))
        printconsole(logio,verbose,string("silhouette plot in ",outfile))         
    end        
    linkagegroups = [resgrouping.connectedsnps[i] for i in resgrouping.clusters]    
    # marker ordering and spacing             
    snpmapls = ordering_spacing(recomnonfrac,recomlod,ldlod, linkagegroups;
        knnorder, knnsave, maxrf, alwayskeep, is_cosgregate_binning, 
        minlodorder, minminlodorder=minlodsave, maxminlod=maxminlodorder,
        isparallel, io = logio, verbose)
    # mapdf cols: "marker","linkagegroup","poscm","physposbp", "binno","represent", "neighbor","neighbornonfrac","neighborrecomlod",...
    mapdf = combine2mapdf(markers, physmapdict, snpmapls, resgrouping,dupebindict,ispermmarker)        
    msg = string("#ungrouped markers: ", sum(mapdf[!,:linkagegroup] .== "NA")," out of ", size(mapdf,1))
    printconsole(logio,verbose,msg)
    # exporting
    if !isnothing(outstem)                 
        # mapdf
        magicmapfile = string(outstem,"_construct_map.csv.gz")
        magicmapfile2 = getabsfile(workdir,magicmapfile)
        GZip.open(magicmapfile2,"w") do mapio
            descripls = ["marker, marker ID", 
                "linkagegroup, linkage group ID", 
                "poscm, marker position in centiMorgan",
                "physchrom, physical map chromosome ID for the marker",
                "physposbp, physical map position (in base pair) for the marker",
                "binno, marker bin index",
                "represent, same as binno if the marker is representative and 0 otherwise",
                "neighbor, marker IDs for the neighbors of the marker",
                "neighbornonfrac, 1 - scaled recombination fractions for the neighbors",
                "neighborrecomlod, LOD scores for the neighbors", 
                "silhouette, silhouette scores for each marker"                
            ]
            msg = ""
            for i in eachindex(descripls)
                msg *= string("##col_", i,",", descripls[i],"\n")
            end
            write(mapio, msg)
            neigenvec = size(mapdf,2) - 11
            msg = string("eigenval1 ...eigenval", neigenvec, ", each column gives the eigenvector corresponding to the eigenvalue in the column name")
            msg = string("##col_",length(descripls)+1, "-",length(descripls)+neigenvec, ",", msg, "\n")
            write(mapio, msg)            
            CSV.write(mapio, mapdf; delim=',', missingstring="NA",header=true,append=true)        
        end
        printconsole(logio,verbose,string("construct_map in ", magicmapfile))    
        # heatmap   
        try 
            snpthin = max(1,round(Int,size(mapdf,1)/10^4))        
            plotrecomheat(ldfile,linkagefile,magicmapfile;  snpthin, 
                workdir, iseachlg = false, outstem = outstem*"_construct",io=logio,verbose)
        catch err
            msg = string(err, ". Cound not plot recomheatmap")
            @warn msg
            printconsole(logio, false, "Warning: "*msg)
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"construct"; verbose,delim="-")   
    isnothing(outstem) ? mapdf : magicmapfile
end

function binning_cosegrate!(markers, nmissingls, recomnonfrac, recomlod, ldlod, dupebindict;
    binrf::Union{Nothing,Real}=nothing, # scaled from 0.0 to 1.0        
    minlod_bin::Real=5.0,        
    )
    nmarker= size(recomlod,1)
    length(nmissingls) == nmarker || @error "inconsistent #marekrs"                
    maxrf,minlodcluster = binrf, minlod_bin    
    similarity = getsimilarity(recomnonfrac,recomlod, ldlod,maxrf, minlodcluster,1.0);
    ii,jj, vv = findnz(similarity)
    all(vv .> 1-maxrf) || @error "unexpected similarity"
    adjmtx = sparse(ii,jj, ones(length(ii)),size(similarity)...)
    @time binls = binning_dupe(adjmtx,nmissingls);
    representls = [bin[argmin(nmissingls[bin])] for bin in binls]
    bindf = DataFrame(represent=representls,bin=binls)
    sort!(bindf,:represent)
    # binning results
    for row in eachrow(bindf)
        snp_represent = markers[row[1]]    
        snp_bin = markers[row[2]]
        if length(snp_bin) > 1   
            ls = [dupebindict[s] for s in snp_bin]
            ls1 = join([s[1] for s in ls],"||")
            ls2 = join([s[2] for s in ls],"||")
            ls3 = join([s[3] for s in ls],"||")
            for snp in setdiff(snp_bin,[snp_represent])
                delete!(dupebindict,snp)
            end                
            dupebindict[snp_represent] = (ls1,ls2,ls3)
        end
    end    
    represents = bindf[!,:represent]
    issorted(represents) || @error "represents must be sorted"
    #  keepat! procues an error when markers is ChainVector, not a vector. 
    # keepat!(markers, represents)    
    markers = markers[represents]    
    nmissingls = nmissingls[represents]
    recomnonfrac = recomnonfrac[represents,represents]
    recomlod = recomlod[represents,represents]
    if !isnothing(ldlod)
        ldlod = ldlod[represents,represents]
    end
    markers, nmissingls, recomnonfrac,recomlod,ldlod
end

function ordering_spacing(recomnonfrac::AbstractMatrix,
    recomlod::AbstractMatrix,    
    ldlod::Union{Nothing,AbstractMatrix},    
    linkagegroups::AbstractVector;    
    knnorder::Union{Nothing,Function},
    knnsave::Union{Nothing,Function},    
    maxrf::Real,
    alwayskeep::Real,
    is_cosgregate_binning::Bool,
    minlodorder::Union{Nothing,Real},
    minminlodorder::Real,
    maxminlod::Real,
    isparallel::Bool=true,    
    io::Union{Nothing,IO} = nothing,
    verbose::Bool=false)
    chridls = [string("LG",i) for i in 1:length(linkagegroups)]
    nonfracls = [recomnonfrac[snps, snps] for snps in linkagegroups]
    recomlodls = [recomlod[snps, snps] for snps in linkagegroups]       
    if isnothing(ldlod)        
        ldlodls = repeat([nothing],length(linkagegroups))
    else
        ldlodls = [ldlod[snps, snps] for snps in linkagegroups]
    end    
    tused = @elapsed if isparallel && nworkers()>1
        # snpmapls[i]: msg, orderedsnps, snppos, binnols, representls, nbrsnp, nbrnonfrac, nbrrecomlod,nbrldlod
        snpmapls = pmap((x,y,z,c,chrid)->ordering_spacing_lg(x,y,z,c; 
            chrid, knnorder, knnsave, maxrf, alwayskeep, is_cosgregate_binning,
            minlodorder,minminlodorder, maxminlod,verbose),
            nonfracls,recomlodls, ldlodls,linkagegroups,chridls)        
    else
        snpmapls = map((x,y,z,c,chrid)->ordering_spacing_lg(x,y,z,c; 
            chrid,knnorder, knnsave, maxrf, alwayskeep, is_cosgregate_binning,
            minlodorder,minminlodorder,maxminlod,verbose),
            nonfracls,recomlodls, ldlodls,linkagegroups,chridls)        
    end
    for i in snpmapls
        write(io,i[1],"\n")
    end
    flush(io)
    mem = round(Int, memoryuse()/10^6);
    msg = string("tused=", round(tused,digits=1), "s for all LGs, mem=", mem,"MB")
    printconsole(io, verbose, msg)
    genlen = [round(Int,last(i[3])) for i = snpmapls]
    msg = string(length(genlen), " genetic lengths: ", join(genlen,","), "cM")
    printconsole(io,verbose,msg)
    snpmapls
end

function ordering_spacing_lg(lgnonfrac::AbstractMatrix,
    lgrecomlod::AbstractMatrix,    
    lgldlod::Union{Nothing,AbstractMatrix},     
    snps::AbstractVector;    
    chrid::AbstractString,
    knnorder::Union{Nothing,Function},
    knnsave::Union{Nothing,Function},    
    maxrf::Real,
    alwayskeep::Real,
    is_cosgregate_binning::Bool,
    minlodorder::Union{Nothing,Real},    
    minminlodorder::Real,
    maxminlod::Real,
    verbose::Bool)
    startt = time()
    # ordering         
    nsnp = length(snps)
    nsnp == size(lgnonfrac,1) || @error "inconsistent #markers"
    mincomponentsize = min(5,1+div(nsnp,5000))     
    mincomponentsize = findmincomponentsize(mincomponentsize, lgnonfrac,lgrecomlod,lgldlod; 
        minlod = minminlodorder,maxrf, alwayskeep)      
    if isnothing(minlodorder)                     
        minlodorder = first(findminlod(lgnonfrac,lgrecomlod,lgldlod;
            maxrf, minminlod = minminlodorder, maxminlod,
            ncomponent=1, mincomponentsize,alwayskeep))
    end
    if isnothing(knnorder)     
        kknnorder =  findknn(lgnonfrac,lgrecomlod,lgldlod; 
            maxrf, minlod = minlodorder, is_cosgregate_binning, mincomponentsize,alwayskeep)
    else
        kknnorder = round(Int,knnorder(nsnp))
        kknnorder = min(max(2, kknnorder),nsnp)
    end        
    similarity, connectednodes = calsimilarity_order(lgnonfrac,lgrecomlod,lgldlod;
        kknnorder, maxrf, minlodorder, mincomponentsize, alwayskeep)                  
    oo = SpectralEmbedding.spectralordering(similarity)  
    nneighbor = get_nneighbor(similarity)
    
    # spacing
    similarity2 = view(similarity,oo,oo)
    # set knnpos <=71 to avoid out of memory during snppos_least_square
    snppos = snppos_least_square(similarity2, min(71,kknnorder))
    nrepresent = length(snppos)
    binnols = collect(1:nrepresent)    
    representls = collect(1:nrepresent)    
    # saving neighbors
    orderednodes = connectednodes[oo]
    orderedsnps = snps[orderednodes]
    if isnothing(knnsave)
        knnneighbor = 10
    else
        knnneighbor = round(Int, knnsave(nsnp))
        knnneighbor = min(max(2, knnneighbor),nsnp)
    end
    nbrsnp, nbrnonfrac, nbrrecomlod, nbrldlod = getneighbor4save(lgnonfrac,lgrecomlod,lgldlod,
        orderednodes, orderedsnps; minlod=minlodorder, knnneighbor,alwayskeep)
    ndel = length(snps)-length(orderedsnps)
    # print msg
    mem = round(Int,memoryuse()/10^6)    
    msg = string(chrid, ", size=",length(snps),
        ", l=",round(Int,last(snppos)), "cM",        
        ", minlod=", round(minlodorder,digits=1),
        ", knnorder=", kknnorder, 
        ", #nbr=", nneighbor,               
        ndel > 0 ? string(", #del=", ndel) : "",
        ", t=", round(time()-startt,digits=1),"s", 
        ", mem=",mem,"MB"
    )
    verbose && (@info msg)    
    msg, orderedsnps, snppos, binnols,representls,nbrsnp, nbrnonfrac, nbrrecomlod,nbrldlod
end

function get_nneighbor(similarity::AbstractMatrix)
    round(Int,(nnz(similarity) - size(similarity,1))/size(similarity,1))
end

function calsimilarity_order(lgnonfrac::AbstractMatrix, 
    lgrecomlod::AbstractMatrix,
    lgldlod::Union{Nothing, AbstractMatrix};        
    maxrf::Real, 
    kknnorder::Integer, 
    minlodorder::Real,
    mincomponentsize::Integer,
    alwayskeep::Real)
    similarity = getsimilarity(lgnonfrac,lgrecomlod, lgldlod,maxrf, minlodorder,alwayskeep)
    snps = 1:size(similarity,1)
    similarity, connectednodes, _ = SpectralEmbedding.dropsingletons(similarity;
        mincomponentsize)
    snps = snps[connectednodes]     
    similarity = SpectralEmbedding.toknnsimilarity(similarity,kknnorder;alwayskeep)
    similarity, connectednodes, _= SpectralEmbedding.dropsingletons(similarity;
        mincomponentsize)
    snps = snps[connectednodes]
    similarity, snps
end

function findncomponent(recomnonfrac::AbstractMatrix,
    recomlod::AbstractMatrix,
    ldlod::Union{Nothing,AbstractMatrix};        
    minlodsave::Real,
    maxminlod::Real,
    maxncomponent::Real,
    maxrf::Real,     
    mincomponentsize::Integer,
    alwayskeep::Real)    
    nconn0 = first(count_connected_components(recomnonfrac,recomlod,
        ldlod,maxrf, maxminlod,mincomponentsize, alwayskeep))
    if nconn0 >= maxncomponent 
        return maxncomponent
    elseif nconn0 == 1
        return 1
    else
        minlod = round(minlodsave,digits=1)
        while true
            nconn0 = first(count_connected_components(recomnonfrac,recomlod,
                ldlod,maxrf, minlod,mincomponentsize, alwayskeep))
            if nconn0 <= maxncomponent 
                if minlod < maxminlod
                    minlod += 0.2
                else
                    return nconn0
                end
            else 
                return maxncomponent 
            end
        end   
    end
    
end

function findmincomponentsize(mincomponentsize::Integer, recomnonfrac::AbstractMatrix,
    recomlod::AbstractMatrix,
    ldlod::Union{Nothing,AbstractMatrix};    
    minlod::Real,    
    maxrf::Real, 
    knn::Union{Nothing,Integer}=nothing,
    alwayskeep::Real)
    mincomponentsize2 =mincomponentsize    
    sim = getsimilarity(recomnonfrac,recomlod,ldlod,maxrf, minlod,alwayskeep)
    if !isnothing(knn)
        sim = SpectralEmbedding.toknnsimilarity(sim,knn;alwayskeep)
    end
    lenls = length.(get_connected_components(sim))    
    lenls = sort(lenls[lenls .>= mincomponentsize],rev=true)    
    if length(lenls) > 1
        nsnp = size(recomnonfrac,1)
        mincc = max(20,round(Int,sqrt(nsnp)/5))
        lenls2 = lenls[lenls .< min(mincc,lenls[2]+1)]
        if !isempty(lenls2)             
            mincomponentsize2 = 1+max(lenls2...)
            mincomponentsize2 >= first(lenls) && (mincomponentsize2 -= 1)
            @info string("connnected component sizes: ",lenls,             
                " after removing size <",mincomponentsize,
                "; reset mincomponentsize=",mincomponentsize2)
        end
    end
    mincomponentsize2
end

function findmaxrf(recomnonfrac::AbstractMatrix, recomlod::AbstractMatrix,
    ldlod::Union{Nothing,AbstractMatrix};
    minlod::Real,
    mincomponentsize::Integer,
    alwayskeep::Real)
    # set # inital maxrf from 1 to 0.45, output maxrf+0.05
    nmarker = size(recomnonfrac,1)
    if nmarker <= 1000
        maxrf = 0.85 - 0.01*div(nmarker,100)
    elseif nmarker <= 10000
        maxrf = 0.75 - 0.01*div(nmarker,1000)
    elseif nmarker <= 50000
        maxrf = 0.65 - 0.04*div(nmarker,10^4)
    else
        maxrf = 0.45
    end         
    function func(maxrf::Real)
        nconn = first(count_connected_components(recomnonfrac,recomlod,
            ldlod,maxrf,minlod,mincomponentsize,alwayskeep))
        nconn
    end
    nconn1 = func(1.0)
    while true
        nconn = func(maxrf)        
        nconn == nconn1 && break
        maxrf += 0.05
    end
    round(min(1.0,maxrf+0.05),digits=2)
end

function findknn(recomnonfrac::AbstractMatrix, recomlod::AbstractMatrix,
    ldlod::Union{Nothing,AbstractMatrix};
    maxrf::Real,
    knnmin::Union{Nothing, Integer}=nothing, 
    is_cosgregate_binning::Bool, 
    minlod::Real,
    mincomponentsize::Integer,
    alwayskeep::Real)
    similarity = getsimilarity(recomnonfrac,recomlod,ldlod,maxrf,minlod,alwayskeep)    
    nsnp = size(recomnonfrac,1)    
    if isnothing(knnmin) 
        if nsnp <= 128
            knnpower = 0.5        
        else
            if is_cosgregate_binning
                knnpower =  nsnp >= 512 ? 9/14 : log2(nsnp)/14            
            else
                knnpower =  nsnp >= 1024 ? 10/14 : log2(nsnp)/14            
            end
        end
        knnmin = round(Int, nsnp^knnpower)
    end    
    knnmin = min(max(10, knnmin),nsnp)        
    knn = first(SpectralEmbedding.findknn(similarity;
        alwayskeep,knnmin,mincomponentsize,verbose=false))
    knn
end

function count_connected_components(recomnonfrac::AbstractMatrix, recomlod::AbstractMatrix,
    ldlod::Union{Nothing, AbstractMatrix},
    maxrf::Real,minlod::Real,mincomponentsize::Integer,        
    alwayskeep::Real)
    sim = getsimilarity(recomnonfrac,recomlod,ldlod,maxrf, minlod,alwayskeep)    
    gg = SimpleGraph(sign.(sim))
    cc = connected_components(gg)
    len = length.(cc)
    b = len .>= mincomponentsize
    nconn = sum(b)    
    nungroup = sum(len[.!b])
    if nconn >=1         
        nbrls = reduce(vcat,[degree(first(induced_subgraph(gg,c)))  for c in cc[b]])        
        # 2 includes one edge and self-loop since the diagonal of sim > 0
        nnode_1nbr = sum(nbrls .== 2)
    else
        nnode_1nbr = 0
    end
    nconn, nungroup, nnode_1nbr
end

function findminlod(recomnonfrac::AbstractMatrix, recomlod::AbstractMatrix,
    ldlod::Union{Nothing,AbstractMatrix};
    mincomponentsize::Integer,
    maxrf::Real,
    minminlod::Real = 1.0,    
    maxminlod::Real = Inf,
    ncomponent::Union{Nothing,Integer} = nothing,    
    alwayskeep::Real=0.99)    
    nsnp = size(recomnonfrac,1)
    minlod = minminlod
    maxnconn = ncomponent    
    nconn0, nungroup0, nnode_1nbr0 = count_connected_components(recomnonfrac,recomlod,
        ldlod,maxrf, minlod,mincomponentsize, alwayskeep)
    nconn0 > maxnconn  && return (minlod,nungroup0,nnode_1nbr0, [[minlod,-minlod^2]])
    nsnp = size(recomnonfrac,1)        
    nconn = nungroup = 0
    diffminlod = maxminlod - minminlod    
    if 20 < diffminlod < 200
        step = max(1, round(Int,diffminlod/20))
    else
        step = 1
    end
    # grid = DataFrame(minlod=[minlod], objv=[oldflod],nconn=[nconn0],nungroup=[nungroup0],nnode_1nbr=[nnode_1nbr0])
    grid = DataFrame(minlod=Float64[], objv=Float64[],nconn=Int[],nungroup=Int[],nnode_1nbr=Int[])
    function func(minlod::Real)
        nconn, nungroup, nnode_1nbr = count_connected_components(recomnonfrac,recomlod,
            ldlod,maxrf,minlod,mincomponentsize, alwayskeep)                    
        objv = if nconn > maxnconn || nsnp - nungroup < 20
            -minlod
        else            
            (minlod^2)*nconn/sqrt(nungroup - nungroup0 + nnode_1nbr - nnode_1nbr0 +1)                        
        end           
        # println("minlod=",minlod,",objv=",objv, ",nconn=",nconn,",maxnconn=", maxnconn, ",nsnp=",nsnp, ",nungroup=",nungroup,",nnode_1nbr=",nnode_1nbr)
        push!(grid, (round(minlod,digits=3), round(objv,digits=3), nconn,nungroup,nnode_1nbr))
        objv     
    end    
    oldflod = func(minlod)
    while true
        minlod += step
        flod = func(minlod)        
        # (flod < 0 || (flod < oldflod && minlod > 5 && size(grid,1) >= 3) || minlod>maxminlod) && break
        (flod < 0 || minlod>maxminlod) && break
        oldflod = flod
    end
    # grid[end,2]>0 && return (grid[end-1,1],nungroup0,nnode_1nbr0,grid)
    if size(grid,1)>2 && grid[1,2] > grid[2,2]
        pos = argmax(grid[2:end,2]) + 1
    else
        pos = argmax(grid[!,2])
    end
    xstart = grid[pos,1]
    xmin, xmax = max(minminlod,xstart-step), min(xstart+step,maxminlod)
    x, fx,his = MagicBase.brentMax(func,xmin, xmax;
        xstart, precisiongoal=3,accuracygoal=3)
    # for i in 2:size(his,1)        
    #     push!(grid, (round(his[i,2],digits=3), round(his[i,3],digits=3), -1,-1))
    # end
    unique!(grid)
    floor(x,digits=3), nungroup0, nnode_1nbr0, grid
end

function snppos_least_square(similarity::AbstractMatrix, knnpos::Integer)
    aa = sign.(similarity)
    sim = sparse(similarity .* (triu(aa,1) .* tril(aa, knnpos+1)))    
    ri,rj, rval = findnz(sim)
    n = size(sim,1)-1
    # nonfrac = 1-frac
    # Haldane frac= 1/2 [1-exp(-2dis)]
    # assume frac= [1-exp(-dis)]
    # rdis in cM
    rdis = (-50.0) .* log.(rval)
    model = Model(OSQP.Optimizer)
    # MOI.set(model, MOI.Silent(), true)
    set_silent(model)
    @variable(model, x[1:n] >= 0)
    quadex = QuadExpr(AffExpr(0.0))
    for z in eachindex(rdis)
        i,j, d = ri[z], rj[z], rdis[z]
        # w = orderedlod[i,j]^2
        w = 1.0
        add_to_expression!(quadex, w*(sum(x[i:(j-1)])-d)^2)
    end
    drop_zeros!(quadex)
    @objective(model, Min, quadex)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        optimal_solution = value.(x)
        # optimal_objective = objective_value(model)
        res = accumulate(+,abs.(optimal_solution))
        pushfirst!(res,0)
        res
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        suboptimal_solution = value.(x)
        # suboptimal_objective = objective_value(model)
        res = accumulate(+,abs.(suboptimal_solution))
        pushfirst!(res,0)
        res
    else
        error("The model was not solved correctly.")
    end
    snppos = round.(res, digits=4)
    snppos
end

function getneighbor4save(nonfrac::AbstractMatrix,
    recomlod::AbstractMatrix,
    ldlod::Union{Nothing,AbstractMatrix}, 
    orderednodes::AbstractVector,     
    orderedsnps::AbstractVector;
    minlod::Real,
    knnneighbor::Integer,
    alwayskeep::Real=0.99)
    # nonfrac, recomlod, ldlod are sparse matrices containing only upper triangular elements
    orderednonfrac = nonfrac[orderednodes,orderednodes]    
    orderedrecomlod = recomlod[orderednodes,orderednodes]
    orderedldlod = ldlod[orderednodes,orderednodes]
    if !isnothing(orderedldlod)
        orderednonfrac .*= orderedldlod .>= min(3,minlod)
    end
    orderednonfrac .*= orderedrecomlod .>= min(3,minlod)
    dropzeros!(orderednonfrac)    
    # extract nbr
    nbrsnp = Vector{String}()
    nbrnonfrac = Vector{String}()
    nbrrecomlod = Vector{String}()
    isld = !isnothing(orderedldlod)
    nbrldlod = isld ? Vector{String}() : nothing
    # orderednonfrac, orderedrecomlod, orderedldlod are symmetric
    for j in 1:size(orderednonfrac,2)
        ii,nonfrac = findnz(orderednonfrac[:,j])
        iip = ii[nonfrac .>= alwayskeep]
        if length(iip) < knnneighbor+1
            p = sortperm(nonfrac,rev=true)
            iip = ii[p[1:min(length(p),knnneighbor+1)]]
            if length(p) >=5 && nonfrac[p[5]] >= 0.9
                deleteat!(iip, orderednonfrac[iip,j] .< 0.9)
            end
        end
        setdiff!(iip,[j])
        length(iip) == knnneighbor+1 && pop!(iip)
        rf  = round.(orderednonfrac[iip,j],digits=5)
        push!(nbrnonfrac,join(rf,"||"))
        push!(nbrsnp,join(orderedsnps[iip],"||"))
        push!(nbrrecomlod,join(orderedrecomlod[iip,j],"||"))
        isld && push!(nbrldlod,join(orderedldlod[iip,j],"||"))
    end    
    nbrsnp, nbrnonfrac, nbrrecomlod, nbrldlod
end

function combine2mapdf(markers::AbstractVector, physmapdict::AbstractDict, 
    snpmapls::AbstractVector, resgrouping::NamedTuple,
    dupebindict::Union{Nothing,AbstractDict},ispermmarker::Bool)
    # snpmapls[i]: msg, orderedsnps, snppos, binnols, representls, nbrsnp, nbrnonfrac, nbrrecomlod,nbrldlod
    mapls = [begin
        orderedsnps, snppos, binnols, representls, neighborsnp, nbrnonfrac, nbrrecomlod, nbrldlod = snpmapls[lg][2:end]
        snpids =markers[orderedsnps]
        snpphys = [get(physmapdict, i, nothing) for i in snpids]
        lgcol =  [string("LG", lg) for _ in 1:length(snpids)]
        rule = Dict(string.(orderedsnps) .=> snpids)
        neighbors = [join([get(rule,i,"NA") for i in j],"||") for j in split.(neighborsnp,"||")]
        # nbrldlod is not included
        [snpids lgcol snppos first.(snpphys) last.(snpphys) binnols representls neighbors nbrnonfrac nbrrecomlod]
    end for lg in eachindex(snpmapls)]
    markermap = vcat(mapls...)        
    groupmarkers = markers[resgrouping.connectedsnps]
    dict = Dict(groupmarkers .=> 1:length(groupmarkers))    
    ii = [get(dict,i, nothing) for i in markermap[:,1]]
    res = hcat(markermap,resgrouping.silhouettes[ii], resgrouping.eigenvecs[ii,:])
    ungroup = setdiff(markers,markermap[:,1])
    resungroup = Matrix{Any}(undef,length(ungroup),size(res,2))
    resungroup[:,1] = ungroup
    resungroup[:,2:end] .= "NA"    
    res = vcat(res, resungroup)
    colnames =[string("eigenval",i, "_",
        round(resgrouping.eigenvals[i],digits=5)) for i=1:length(resgrouping.eigenvals)]
    colnames = vcat(["marker","linkagegroup","poscm","physchrom", "physposbp", "binno","represent",
        "neighbor","neighbornonfrac","neighborrecomlod", "silhouette"], colnames)
    mapdf = DataFrame(res,colnames)    
    # put back duplicated markers
    if !isnothing(dupebindict)
        col_eigen1 = size(mapdf,2) - length(resgrouping.eigenvals) + 1
        mapdf = reduce(vcat, [begin 
            df = mapdf[[i],:]
            if df[1,:linkagegroup] != "NA"
                snpid = df[1,:marker]      
                binmap = dupebindict[snpid] 
                binsnps = split(binmap[1],"||")
                binphyschrom = split(binmap[2],"||")
                binphysposbp = split(binmap[3],"||")
                binsize = length(binsnps)
                if ispermmarker && binsize > 1                    
                    oo = sample(1:binsize,binsize,replace=false)
                    binsnps .= binsnps[oo]
                    length(binphyschrom) == binsize && (binphyschrom .= binphyschrom[oo])
                    length(binphysposbp) == binsize && (binphysposbp .= binphysposbp[oo])
                end
                isrepresent =  binsnps .== snpid 
                sum(isrepresent)==1 || @error string("represent=",snpid, " is not in its bin=",binsnps)                
                if binsize != 1                 
                    # set binning
                    df = mapdf[repeat([i],binsize),:]
                    df[.!isrepresent,col_eigen1:end] .= "NA" 
                    df[!,:marker] .= binsnps
                    df[!,:binno] .= i
                    df[!,:represent] .= 0
                    df[isrepresent,:represent] .= i                                            
                    df[!,:physchrom] .= binphyschrom
                    df[!,:physposbp] .= binphysposbp
                    # reset nbr
                    nbrls = split(df[1,:neighbor],"||")
                    nneighbor = length(nbrls)
                    pushfirst!(nbrls,binsnps...)
                    unique!(nbrls)
                    n_diff = length(nbrls) - nneighbor -1 # 1 denotes the marker itself
                    if n_diff >= 1
                        df[!,:neighbor] .= [join(setdiff(nbrls,[i]),"||") for i in binsnps]
                        df[!,:neighbornonfrac] .= join([join(ones(n_diff),"||"),df[1,:neighbornonfrac]],"||")
                        df[!,:neighborrecomlod] .= join([join(100*ones(n_diff),"||"),df[1,:neighborrecomlod]],"||")
                    end
                end
            end
            df
        end for i in 1:size(mapdf,1)])
    end    
    mapdf
end

function tocorrmtx(lodmtx::SparseMatrixCSC)
    nrow, ncol = size(lodmtx)
    @assert nrow == ncol
    sqrtdiag = sqrt.(diag(lodmtx))
    ils, jls, vls  = findnz(lodmtx)
    df = DataFrame(i=ils, j=jls, v=vls)
    gdf = groupby(df,[:i])
    df2 = reduce(vcat, [begin
        d = DataFrame(df)
        d[!,:v] ./= (sqrtdiag[d[1,:i]]*sqrtdiag[d[!,:j]])
        d
    end for df in gdf])
    sparse(df2[!,:i],df2[!,:j], df2[!,:v], nrow,ncol)
end
