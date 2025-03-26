
function marker_grouping(recomnonfrac::AbstractMatrix,
    recomlod::AbstractMatrix,    
    ldlod::Union{Nothing,AbstractMatrix};    
    ncluster::Union{Nothing,Real}, 
    minncluster::Real, 
    maxncluster::Real, 
    maxrf::Real,
    alwayskeep::Real,
    minlodcluster::Real,    
    kknncluster::Integer,
    mincomponentsize::Integer,
    minsilhouette::Union{Nothing,Real},
    clusteralg::Union{AbstractString,AbstractVector},    
    eigselect::AbstractString,
    eigweightfrac::Real,
    eigweighttype::AbstractString,
    io::Union{Nothing,IO} = nothing,
    verbose::Bool=false)    
    # grouping
    similarity, connectedsnps = calsimilarity_cluster(recomnonfrac,recomlod; 
        ldlod = ldlod, ncluster = isnothing(ncluster) ? div(maxncluster, 2) : ncluster, 
        kknncluster, maxrf, alwayskeep, minlodcluster, mincomponentsize, io,verbose)    
    if !issubset(connectedsnps,1:size(recomnonfrac,1)) 
        msg = "wrong marker indices that are the row/column names of similarity"
        @error msg
        printconsole(io, verbose, "Error: "*msg)
    end
    nclusters = collect(minncluster:maxncluster)
    isparallel = in(clusteralg,["kmeans","randswap"]) && (isnothing(ncluster) || ncluster < maxncluster)
    resgrouping = marker_grouping(similarity, nclusters; mincomponentsize, 
        minsilhouette, clusteralg, eigselect,eigweightfrac,eigweighttype, isparallel, io,verbose)
    (resgrouping..., connectedsnps = connectedsnps)
end

function marker_grouping(similarity::AbstractMatrix,        
    nclusters::AbstractVector;
    mincomponentsize::Integer=1,
    minsilhouette::Union{Nothing,Real}=nothing,
    clusteralg::Union{AbstractString,AbstractVector}="kmeans",
    eigselect::AbstractString="eigratio", 
    eigweightfrac::Real=0.01, 
    eigweighttype::AbstractString="eiggapratio",     
    isparallel::Bool=true,    
    io::Union{Nothing,IO} = nothing,
    verbose::Bool=false)
    startt = time()
    if length(nclusters) == 1
        # nev = number of eigenvectors, if nothing, it is determined by eigselect method. 
        nc = only(nclusters)        
        spec_res = SpectralEmbedding.spectralclustering(similarity,nc;
            clusteralg,nev = nothing, eigselect, minclustersize = mincomponentsize, isparallel, io, verbose)
    else
        spec_res = SpectralEmbedding.spectralclustering(similarity,nclusters;
            clusteralg,minclustersize = mincomponentsize, eigweightfrac, eigweighttype, isparallel, io, verbose)
    end    
    clusters, eigenvals, eigenvecs, silhouettes, res_silhouettes =spec_res
    msg = string(length(clusters), " group sizes: ", join(length.(clusters),","))
    printconsole(io,verbose,msg)    
    # remove markers with silhouette < minsilhouette    
    pos = argmax(SpectralEmbedding.weighted_avgsilhls(res_silhouettes.eigweightls, res_silhouettes.avgsilhls; eigweightfrac))
    silhs = res_silhouettes.silhsls[pos] 
    clusters2 = repeat_trunc_clusters_by_silh(clusters,silhs; minsilhouette, io, verbose)  
    clusters2 = clusters2[length.(clusters2) .>= mincomponentsize]        
    res_silhouettes.clustersls[pos] = clusters2
    # print info        
    mem = round(Int, memoryuse()/10^6);    
    msg = string("marker_grouping, tused=", round(time()-startt,digits=1),"s, mem=",mem, "MB")
    printconsole(io,verbose,msg)
    # keep using cluster without deleting markers
    if !issubset(reduce(vcat,clusters2), 1:size(similarity,1))
        msg = "wrong clustering"
        @error msg
        printconsole(io, verbose, "Error: "*msg)        
    end        
    # return 
    (clusters=clusters2, eigenvals=eigenvals, eigenvecs=eigenvecs, silhouettes=silhouettes, res_silhouettes=res_silhouettes)
end

function repeat_trunc_clusters_by_silh(clusters::AbstractVector, silhs::AbstractVector; 
    minsilhouette::Union{Nothing,Real}, io::Union{Nothing, IO},verbose::Bool)
    if isnothing(minsilhouette)
        minsilhouette=0.5
        clusters2, maxdiff, msg = trunc_clusters_by_silh(clusters, silhs; minsilhouette)
        printconsole(io,verbose,msg)    
        while maxdiff > 0.2 && minsilhouette > 0.0
            minsilhouette -= 0.1
            clusters2, maxdiff, msg = trunc_clusters_by_silh(clusters, silhs; minsilhouette)
            printconsole(io,verbose,msg)                
        end        
    else
        clusters2, maxdiff, msg = trunc_clusters_by_silh(clusters, silhs; minsilhouette)
        printconsole(io,verbose,msg)    
        if maxdiff >= 0.2 # Kepp consistent with mapcorndel(minfreq=0.2)
            msg = string("minsilhouette might be too large!")
            if !isnothing(minsilhouette) && minsilhouette > 0
                msg *= string(" Suggest to set minsilhouette = 0")
            end           
            @warn msg
            printconsole(io,verbose,"WARN: "*msg)    
        end    
    end        
    clusters2
end

function trunc_clusters_by_silh(clusters::AbstractVector, silhs::AbstractVector; minsilhouette::Real)
    badii = findall(silhs .< minsilhouette)    
    clusters2 = [setdiff(i,badii) for i in clusters]        
    # check differences 
    nls = length.(clusters)
    nls2 = length.(clusters2)
    maxdiff, lgdiff = findmax((nls .- nls2) ./ nls)        
    msg = string(length(nls2), " group sizes: ", join(nls2,","), "; maxdiff=", round(100*maxdiff,digits=1), "% in LG=", lgdiff)
    msg *= string("; #del_markers = ", sum(nls .- nls2), " with silhouette < ", minsilhouette)
    clusters2, maxdiff, msg
end

function calsimilarity_cluster(recomnonfrac::AbstractMatrix,
    recomlod::AbstractMatrix;
    ldlod::Union{Nothing,AbstractMatrix},
    ncluster::Integer,
    kknncluster::Integer,
    maxrf::Real,
    alwayskeep::Real,
    minlodcluster::Real,
    mincomponentsize::Integer,    
    io::Union{Nothing,IO} = nothing,
    verbose::Bool=false)
    similarity = getsimilarity(recomnonfrac,recomlod, ldlod,maxrf, minlodcluster,alwayskeep)
    cc = get_connected_components(similarity)
    cclenls = sort(length.(cc),rev=true)
    deleteat!(cclenls, cclenls .== 1)
    msg = string("sizes of ", length(cclenls)," components(minlodcluster=",
        minlodcluster, "): ", join(cclenls,","),
        " after removing omponents of size 1")
    printconsole(io,verbose,msg)
    ncc = sum(length.(cc) .>= mincomponentsize)
    if ncc == 0
        msg = string("No connected components")
        @error msg
        printconsole(io, verbose, "Error: "*msg) 
    elseif ncc > ncluster
        msg = string("# connected components = ", ncc, " is larger than ncluster =", ncluster)
        @error msg
        printconsole(io, verbose, "Error: "*msg) 
    end
    snps = 1:size(similarity,1)
    similarity, connectednodes, _ = SpectralEmbedding.dropsingletons(similarity;
        mincomponentsize)
    snps = snps[connectednodes]    
    similarity = SpectralEmbedding.toknnsimilarity(similarity,kknncluster;alwayskeep)       
    similarity, connectednodes, _= SpectralEmbedding.dropsingletons(similarity;
        mincomponentsize)
    snps = snps[connectednodes]
    similarity, snps
end

function get_connected_components(similarity::AbstractMatrix)
    g = SimpleGraph(sign.(similarity))
    connected_components(g)
end

function getsimilarity(recomnonfrac::AbstractMatrix, recomlod::AbstractMatrix,
    ldlod::Union{Nothing, AbstractMatrix},maxrf::Real,minlod::Real,
    alwayskeep::Real)
    maxlod = max(diag(recomlod)...)
    res = (recomlod ./ (maxlod*10000)) .+ recomnonfrac    
    b = @. (recomlod > minlod) && (recomnonfrac >= 1.0-maxrf)  || (recomnonfrac >= alwayskeep)
    b[diagind(b)] .= true
    if !isnothing(ldlod)        
        b .*= (ldlod .- minlod) .> 0
    end        
    sim = dropzeros(res .* b)
    sim
end

function plot_eigen(magicmapfile;
    ispermute::Bool=true,
    workdir::AbstractString=pwd(),
    missingstring = ["NA","missing"])
    magicmapfile2 = getabsfile(workdir,magicmapfile)
    mapdf = CSV.read(magicmapfile2, DataFrame; missingstring)
    eigencol = findfirst(x->occursin(r"^eigenval1_",x), names(mapdf))
    vals = parse.(Float64,last.(split.(names(mapdf)[eigencol:end],"_")))
    b = .!ismissing.(mapdf[!, eigencol])
    vecs = Matrix{Float64}(mapdf[b,eigencol:end])    
    ls = mapdf[b,:linkagegroup]
    clusters = [findall(ls .== i) for i in unique(ls)]
    if ispermute
        clusters = [sample(i,length(i),replace=false) for i in clusters]
    end
    SpectralEmbedding.plot_eigen(vals,vecs,clusters)
end
