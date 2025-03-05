function knn_spectralclustering(similarity::AbstractMatrix, ncluster::Integer;
    nev::Union{Nothing,Integer}=ncluster, 
    ncomponent::Integer=1,
    knnmin::Integer = round(Int,sqrt(size(similarity,1))),
    mincomponentsize::Integer=1,
    laplaciantype::AbstractString= "randomwalk",
    clusteralg::AbstractString="kmeans",
    eigselect::AbstractString="eigratio",
    minclustersize::Integer=5,
    alwayskeep::Union{Nothing,Real}=nothing,
    isparallel::Bool=true,
    io::Union{Nothing,IO} = nothing,
    verbose::Bool=false)
    if !in(clusteralg,["kmeans","kmedoids", "ksilswap", "randswap","hclust"]) 
        msg = string("clusteralg=",clusteralg, ", not in [kmeans,kmedoids,ksilswap,randswap,hclust]")
        printconsole(io, verbose,"ERROR: "*msg)
        @error msg
    end    
    knn, knn_history = findknn(similarity; ncomponent,knnmin,mincomponentsize,
        alwayskeep, io, verbose)
    similarity2=toknnsimilarity(similarity,knn;alwayskeep)
    similarity2, connectednodes, singletons = takecomponents(similarity2; ncomponent)
    clusters, eigenvals, eigenvecs, silhouettes, res_silhouettes= spectralclustering(similarity2,ncluster;
        laplaciantype, clusteralg, nev, minclustersize,eigselect, isparallel,io, verbose)
    if isempty(singletons)
        clusters2 = clusters
    else
        clusters2 = [connectednodes[i] for i=clusters]
    end
    (clusters = clusters2, eigenvals = eigenvals, eigenvecs = eigenvecs, 
        silhouettes=silhouettes,res_silhouettes = res_silhouettes, knn_history = knn_history)
end

function spectralclustering(similarity::AbstractMatrix, ncluster::Integer;    
    laplaciantype::AbstractString= "randomwalk",
    clusteralg::Union{AbstractString,AbstractVector}="kmeans",
    nev::Union{Nothing,Integer}=ncluster, 
    minclustersize::Integer=5,    
    eigselect::AbstractString="eigratio",
    isparallel::Bool=true,
    io::Union{Nothing,IO}=nothing,    
    verbose::Bool=true)
    clusteralgls = isa(clusteralg, AbstractString) ? [clusteralg] : clusteralg
    for alg in clusteralgls
        if !in(alg,["kmeans","kmedoids", "ksilswap", "randswap","hclust"]) 
            msg = string("clusteralg=",alg, ", not in [kmeans,kmedoids,ksilswap,randswap,hclust]")
            printconsole(io, verbose,"ERROR: "*msg)
            @error msg
        end
    end
    n = size(similarity,1)
    ncluster < n || @error string("ncuster =",ncluster, " is not <= #row of similarity!")
    init_nev = isnothing(nev) ? min(n, ceil(Int,1.5*ncluster+1)) : nev + 1
    vals, vecs = spectralembedding(similarity; nev = init_nev, laplaciantype)         
    if isnothing(nev)                   
        eigweightls = get_eigweightls(vals[ncluster:end]; eigweighttype = eigselect)        
        nev = argmax(eigweightls) + ncluster - 1      
        msg = string("ncluster=", ncluster, ", nev=", nev, ", eigenval=",round.(vals,digits=4))
        printconsole(io,verbose,msg)
    end    
    # res_silhouettes nametuple: eigweightls, avgsilhls, silhsls,clustersls
    if isa(clusteralg, AbstractVector)        
        resls = [clustering_silhouettes(vals,vecs, ncluster; nev = nev, clusteralg = alg, 
            minclustersize, eigweighttype = eigselect, logio=io, verbose) for alg in clusteralg]        
        res = first.(resls)
        eigweightls, avgsilhls, silhsls, clustersls = collect.(eachrow(reduce(hcat, res)))    
        res_silhouettes = (eigweightls=eigweightls, avgsilhls=avgsilhls, silhsls=silhsls,clustersls=clustersls)        
        pos = argmax(avgsilhls)    
    else
        res_silhouettes = clustering_silhouettes(vals,vecs, [ncluster]; 
            nevls = [nev], minclustersize, clusteralg, eigweighttype = eigselect, 
            isparallel, logio=io, verbose)        
        pos = 1
    end    
    clusters = res_silhouettes.clustersls[pos]        
    silhouettes = res_silhouettes.silhsls[pos]
    (clusters = clusters, eigenvals = vals, eigenvecs = vecs, silhouettes=silhouettes,res_silhouettes = res_silhouettes)
end

function get_eigweightls(eigvalls; eigweighttype::AbstractString="eiggapratio")
    smallev = 1e-4
    if eigweighttype == "eiggap"
        [eigvalls[i+1]-eigvalls[i] for i in 1:length(eigvalls)-1]    
    elseif eigweighttype == "eigratio"
        [(eigvalls[i+1] + smallev) /(eigvalls[i] + smallev) for i in 1:length(eigvalls)-1]    
    elseif eigweighttype == "eiggapratio"
        [(eigvalls[i+1]-eigvalls[i]) * (eigvalls[i+1] + smallev) /(eigvalls[i] + smallev) for i in 1:length(eigvalls)-1]    
    else
        @error string("unknown eigweighttype=",eigweighttype)
    end
end

function spectralclustering(similarity::AbstractMatrix, nclusters::AbstractVector;
    laplaciantype::AbstractString= "randomwalk",
    clusteralg::AbstractString="kmeans",
    minclustersize::Integer=5,
    eigweightfrac::Real=0.01,
    eigweighttype::AbstractString="eiggapratio",
    isparallel::Bool=true,
    io::Union{Nothing,IO}=nothing,    
    verbose::Bool=true)
    n = size(similarity,1)
    ncls = sort(collect(nclusters))    
    ncmax = last(ncls)
    ncmax<n || @error string("max ncuster =",ncmax, " is not <= #row of similarity!")
    # vecs is a matrix of size (n, nev); each column correponds to a eigenvector
    if ncmax ==1
        n = size(similarity,1)        
        vals = zeros(n)
        vecs = zeros(n,ncmax+1)
    else
        vals, vecs = spectralembedding(similarity; nev = ncmax+1, laplaciantype)        
    end
    # res_silhouettes nametuple: eigweightls, avgsilhls, silhsls,clustersls
    res_silhouettes = clustering_silhouettes(vals,vecs, ncls; 
        nevls = ncls, minclustersize, clusteralg, eigweighttype,
        isparallel, logio=io, verbose)    
    wsilh = weighted_avgsilhls(res_silhouettes.eigweightls, res_silhouettes.avgsilhls; eigweightfrac)
    pos = argmax(wsilh)
    clusters = res_silhouettes.clustersls[pos]    
    silhouettes = res_silhouettes.silhsls[pos]
    (clusters = clusters, eigenvals = vals, eigenvecs = vecs, silhouettes=silhouettes,res_silhouettes = res_silhouettes)
end


function weighted_avgsilhls(eigweightls::AbstractVector, avgsilhls::AbstractVector;
    eigweightfrac::Real=0.0)    
    eigweightfrac <= 0 && return avgsilhls
    if length(eigweightls)==1
        wls  = ones(1)
    else
        wls = eigweightls ./ maximum(eigweightls)
    end
    @. avgsilhls * (one(1.0) + eigweightfrac*wls)
end

function clustering_silhouettes(eigvals::AbstractVector, eigvecs::AbstractMatrix, nclusters::AbstractVector;    
    nevls::AbstractVector = nclusters,
    minclustersize::Integer=5,
    clusteralg::AbstractString="kmeans",
    eigweighttype::AbstractString="eiggapratio",
    isparallel::Bool=true,
    logio::Union{Nothing,IO},    
    verbose::Bool)
    ncls = collect(nclusters)        
    (length(eigvals) < max(ncls...)+1) && @error string("#eigenvals ", length(eigvals), " < #clusters(", max(ncls...), ") +1")    
    resls = Vector(undef,length(ncls))
    if isparallel && nprocs()>1
        try
            resls = pmap((x,y)->clustering_silhouettes(eigvals, eigvecs, x;
                nev = y, minclustersize, clusteralg, eigweighttype, logio=nothing, verbose), ncls,nevls)
        finally
            if !isnothing(logio)
                for msg in last.(resls)
                    write(logio, msg) # msg has \n 
                end
            end
        end        
    else
        resls = map((x,y)->clustering_silhouettes(eigvals, eigvecs, x;
            nev = y, minclustersize, clusteralg, eigweighttype, logio, verbose), ncls, nevls)        
    end
    res = first.(resls)
    eigweightls, avgsilhls, silhsls, clustersls = collect.(eachrow(reduce(hcat, res)))    
    (eigweightls=eigweightls, avgsilhls=avgsilhls, silhsls=silhsls,clustersls=clustersls)
end



function clustering_silhouettes(eigvals::AbstractVector, eigvecs::AbstractMatrix, nc::Integer;    
    nev::Integer = nc,
    minclustersize::Integer=5,        
    clusteralg::AbstractString="kmeans",
    eigweighttype::AbstractString="eiggapratio",
    logio::Union{Nothing,IO}=nothing,    
    verbose::Bool=true)    
    io = isnothing(logio) ? IOBuffer() : logio
    try 
        if nc <= 1
            nobject = size(eigvecs,1)
            sils = zeros(nobject)
            clusters = [1:nobject]            
            eiggap = 0.0
            return ([eiggap, first(sils), sils, clusters],"")    
        end        
        # data2 = view(eigvecs, :, 1:nc)        
        data2 = eigvecs[:,1:nev]'
        for i in eachcol(data2) 
            i .= normalize(i,2) # normalize eachrow
        end      
        # @info string("#rows of kemeans data = ", size(data2,1))  
        tused_cl = @elapsed if clusteralg == "kmeans"
            kclustalg = clusteralg
            assignments, counts = robust_kclust(data2, nc; 
                kclustalg, maxrepeat = 5000,                
                minclustersize,io, verbose)            
        elseif clusteralg == "kmedoids"
            kclustalg = clusteralg
            dists = pairwise(Euclidean(), data2; dims=2)
            assignments, counts = robust_kclust(dists, nc; 
                kclustalg, maxrepeat = 5000,                
                minclustersize,io, verbose)    
        elseif clusteralg == "hclust"                        
            # https://myscale.com/blog/power-cosine-similarity-vs-euclidean-distance-explained/
            # Opt for Cosine Similarity in high-dimensional data or text analysis 
            # where vector magnitude is not critical; select Euclidean Similarity 
            # for lower-dimensional spaces where vector magnitude plays a vital role.
            # Tests for real maize data show that Cosine sometimes results in wrong unbalanced group sizes. 
            dists = pairwise(Euclidean(), data2; dims=2)
            # dists = pairwise(CosineDist(), data2; dims=2)
            # critical option: linkage = :single, see below
            # https://stats.stackexchange.com/questions/195446/choosing-the-right-linkage-method-for-hierarchical-clustering
            assignments, counts = robust_hclust(dists, nc; 
                linkage = :single, 
                minclustersize,io, verbose
            )            
        elseif clusteralg in ["ksilswap", "randswap"]
            kclustalg = clusteralg
            assignments, counts = robust_kclust(data2, nc; 
                kclustalg, maxrepeat = 10,                
                minclustersize,io, verbose)            
        else
            msg = string("unknown clusteralg=",clusteralg)          
            printconsole(io,verbose, "ERROR: "*msg)
            @error msg
        end
        # assignment = 0 for uncluster points whose cluster size < minclustersize
        assigns, clusters = trunc_clusters(assignments, counts; minclustersize)        
        b = assigns .> 0    
        tused_sil = @elapsed if in(clusteralg,["kmeans","ksilswap","randswap"])
            sils = cal_silhouettes(data2, clusters; sil_default=-1.0); # sil_default -1 for uncluster points
        else
            sils = -ones(size(data2,2)) # set sil = -1 for uncluster points
            # sils[b] .= silhouettes(assigns[b], length.(clusters), view(dists,b,b))
            if length(unique(assigns[b]))  > 1
                sils[b] .= silhouettes(assigns[b], view(dists,b,b))
            end
        end    
        avgsil = mean(sils)
        eigweight = first(get_eigweightls(eigvals[nev:nev+1]; eigweighttype))        
        msg = string("ncluster=",nc, ", nev=",nev,
            ", silhouette=", round(avgsil,digits=4),                         
            ", eigweight=", round(eigweight,digits=4),                         
            ", t(cluster|silh)=", join(round.([tused_cl,tused_sil],digits=2),"|"), "s")    
        printconsole(io, verbose, msg)        
        msg = isnothing(logio)  ? String(take!(io)) : ""
        ([eigweight, avgsil, sils, clusters], msg)
    finally
        isnothing(logio) && close(io)
    end
end

function trunc_clusters(assignments::AbstractVector,counts::AbstractVector;
    minclustersize::Integer=5)
    clusters = [findall(assignments .== i) for i in findall(counts .>= minclustersize)]
    oo = sortperm(length.(clusters),rev=true)
    clusters = clusters[oo]
    assigns = zeros(Int, length(assignments))
    counts = length.(clusters)
    for i in eachindex(clusters)
        assigns[clusters[i]] .= i
    end
    assigns, clusters
end

function robust_kclust(data::AbstractMatrix, nc::Integer;  
    kclustalg::AbstractString="kmeans",    
    minclustersize::Integer=5,
    io::Union{Nothing,IO}=nothing,    
    verbose::Bool=true,
    maxrepeat::Integer = (in(kclustalg,["ksilswap","randswap"]) ? 10 : 1000))
    newnc = nc        
    itmax = nc    
    for it = 1:itmax
        # kmeans clustering        
        reskclust = first(bestkclust(data, newnc; kclustalg, maxrepeat))        
        bool =  reskclust.counts .>= minclustersize
        nlarge = sum(bool)
        if nlarge == nc            
            break
        elseif  nlarge < nc
            newnc += (nc-nlarge)
        else
            newnc -= (nlarge-nc)
        end        
        if it == itmax
            msg = string("couldn't find desired ncluster = ", nc, ", return ncluster =", nlarge)            
            verbose && @warn msg
            printconsole(io, false, "Warning: "*msg)
            break
        end
    end            
    tused = @elapsed reskclust, reshis, itlast = bestkclust(data, newnc; kclustalg, maxrepeat) 
    if in(kclustalg,["ksilswap","randswap"]) && verbose 
        println(reskclust.acceptedswap)
    end
    costs = [i[1] for i in reshis]
    lenls = [i[2] for i in reshis]
    nc_est = sum(reskclust.counts .>= minclustersize)
    if length(unique(lenls)) > 1
        msg = string("nc0=", nc,", nc=", nc_est, ", k=",newnc)
        msg *= string(", #repeat_", kclustalg, "=", itlast, ", tused=",round(tused,digits=1),"s")
        msg *= string("; inconsistent clusterings = ", join(join.(lenls,"|"),", "))
        msg *= string("; costs = ", join(round.(costs,digits=2),","))
        verbose && @warn msg
        printconsole(io, false, "Warning: "*msg)
    end 
    reskclust.assignments, reskclust.counts            
end


function robust_hclust(dists::AbstractMatrix, nc::Integer;
    minclustersize::Integer=5,
    linkage = :single,
    io::Union{Nothing,IO}=nothing,    
    verbose::Bool=true)
    res_hl = hclust(dists; linkage)    
    newnc = nc      
    itmax = nc  
    for it in 1:itmax        
        assigns = cutree(res_hl; k=newnc)
        counts = [sum(assigns .== i) for i in 1:newnc]
        bool =  counts .>= minclustersize
        nlarge = sum(bool)
        if nlarge == nc                  
            return assigns, counts
        elseif  nlarge < nc
            newnc += (nc-nlarge)
        else
            newnc -= (nlarge-nc)
        end        
        if it == itmax
            msg = string("couldn't find desired ncluster = ", nc, ", return ncluster =", nlarge)            
            verbose && @warn msg
            printconsole(io, false, "Warning: "*msg)
            return assigns, counts
        end
    end    
end

function bestkclust(data::AbstractMatrix, ncluster::Integer; 
    kclustalg::AbstractString="kmeans",    
    maxrepeat::Integer=500,
    maxstuck = in(kclustalg,["ksilswap", "randswap"]) ? div(maxrepeat,2) : div(maxrepeat,10))
    if kclustalg == "randswap"
        reskclust = randswap(data, ncluster)
    elseif kclustalg == "ksilswap"
        reskclust = ksilswap(data, ncluster)
    elseif kclustalg == "kmeans"
        reskclust = kmeans(data, ncluster)
    elseif kclustalg == "kmedoids"
        reskclust = kmedoids(data, ncluster)
    end    
    mincost =  kclustalg == "ksilswap" ? -1*reskclust.totalsilhouette : reskclust.totalcost
    counts = sort(reskclust.counts,rev=true)
    resls = [[mincost, counts, reskclust]]
    nstuck = 0    
    itlast = [maxrepeat]
    for it  in 1:maxrepeat
        if kclustalg == "randswap"
            reskclust = randswap(data, ncluster)
        elseif kclustalg == "ksilswap"
            reskclust = ksilswap(data, ncluster)
        elseif kclustalg == "kmeans"
            reskclust = kmeans(data, ncluster)
        elseif kclustalg == "kmedoids"
            reskclust = kmedoids(data, ncluster)
        end        
        cost =  kclustalg == "ksilswap" ? -1*reskclust.totalsilhouette : reskclust.totalcost
        counts = sort(reskclust.counts,rev=true)
        if isapprox(cost, mincost, atol=1e-4)
            # push!(resls, [cost, counts, reskclust])
            nstuck += 1
        elseif  cost < mincost
            push!(resls, [cost, counts, reskclust])
            mincost = cost
            nstuck = 0
            # nstuck = sum([isapprox(i[1], mincost,atol=1e-4) for i in resls])
        else            
            continue
        end
        if nstuck == maxstuck                        
            itlast[1] = it
            break
        end
    end
    oo = sortperm([i[1] for i in resls])
    resls = resls[oo]
    resls[1][3], resls, itlast[1]
end


function plot_eigen(eigenvals::AbstractVector, eigenvecs::AbstractMatrix, clusters::AbstractVector)
    gval = scatter(eigenvals,
        xlabel="eigenvalue index",
        ylabel="eigenvalue",
        legend=false,        
    )
    plot!(eigenvals)
    # clusters is a list of vectors; each vector is a list of row indices of eigenvecs
    rowls = vcat(clusters...)
    rowls2 = setdiff(1:size(eigenvecs,1),rowls) 
    rowls2 = vcat(rowls,rowls2)
    eigenvecs2 = view(eigenvecs,rowls2,:)
    gvec = plot(eigenvecs2,
        xlabel="object index",
        ylabel="eigenvector",
        lw = 2,
        legend=false,        
    )
    plot(gval,gvec;
        size = (1200,600), 
        left_margin=20Plots.mm,         
        bottom_margin=15Plots.mm,        
    )
end
