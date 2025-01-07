function cal_silhouettes(data::AbstractMatrix, clusters::AbstractVector; 
    sil_default::Real = -1.0)    
    n = size(data, 2)
    di = zeros(1, n) # di saves distances between row i and each of all rows
    sils = sil_default * ones(n)
    singletons = clusters[length.(clusters) .== 1]
    if !isempty(singletons) 
        sils[reduce(vcat, singletons)] .= 0.0
    end    
    nc = length(clusters)
    for c in eachindex(clusters)    
        cluster = clusters[c]
        c_rest = setdiff(1:nc, [c])
        for i in cluster
            if length(cluster) > 1                
                pairwise!(Euclidean(), di, view(data,:, i:i),data; dims=2)
                ai = (sum(di[1,cluster])-di[1,i])/(length(cluster)-1)
                bi = minimum([mean(di[1,clusters[c2]]) for c2 in c_rest])
                sils[i] = (bi-ai)/max(ai,bi)
            end
        end
    end
    sils
end

function summary_silhouette(res_silhouettes::NamedTuple; eigweightfrac::Real=0.0)    
    ncls = length.(res_silhouettes.clustersls)
    weightls = res_silhouettes.eigweightls
    avgsilhls = res_silhouettes.avgsilhls
    if eigweightfrac <= 0
        DataFrame(ncluster = ncls, eigweight=weightls,average_silhouette=avgsilhls)
    else
        wsilhls = weighted_avgsilhls(weightls,avgsilhls; eigweightfrac)
        DataFrame(ncluster = ncls, eigweight=weightls, eigweight_scaled = weightls ./ maximum(weightls),average_silhouette=avgsilhls, weighted_silhouette=wsilhls)
    end
end

function plot_silhouette(res_silhouettes::NamedTuple; eigweightfrac::Real=0.0)        
    eigweightls, avgsilhls, silhsls,clustersls = res_silhouettes    
    if eigweightfrac > 0
        avgsilhls_w = SpectralEmbedding.weighted_avgsilhls(eigweightls, avgsilhls;eigweightfrac)
        pos = argmax(avgsilhls_w)
    else
        pos = argmax(avgsilhls)
    end    
    kls = length.(clustersls)
    leftpos = max(1,pos-1)
    rightpos = min(length(kls),pos+1)
    posls = range(leftpos,rightpos)    
    gls = [SpectralEmbedding.plotsilhouete(silhsls[pos],clustersls[pos];
        xlabel = "Index",
        ylabel = "Silhouette",
        title=string("k=",kls[pos], ", <silhouette>=",round(mean(silhsls[pos]),digits=4))
    ) for pos in posls]
    if length(gls) == 1        
        plotsize = (1200,400)
    else        
        g1 = plot(kls, avgsilhls;
            legend = eigweightfrac > 0,
            label="average silhouette",
            xlabel = "#clusters k",
            ylabel = "silhouette",        
            color = :blue,             
        )        
        pos = argmax(avgsilhls)
        x = kls[pos]
        y = round(avgsilhls[pos],digits=4)
        scatter!(g1,[x],[y];                     
            label = "", 
            marker = (:star, 6,:blue),            
            series_annotations = text((x,y), :top,:blue),            
        )       
        if eigweightfrac > 0
            plot!(g1, kls, avgsilhls_w;
                label="weighted silhouette",
                xlabel = "#clusters k",
                ylabel = "silhouette",        
                color = :red,                
            )    
            pos = argmax(avgsilhls_w)
            x = kls[pos]
            y = round(avgsilhls_w[pos],digits=4)
            scatter!(g1,[x],[y];                     
                label = "", 
                marker = (:cross, 6,:red),            
                series_annotations = text((x,y), :bottom,:red),            
            )        
        end
        push!(gls,g1)
        plotsize = (1200,1200)
    end
    plot(gls...;    
        layout = (length(gls),1),
        size = plotsize,
        left_margin=25Plots.mm,         
        bottom_margin=15Plots.mm,    
    )       
end

function plotsilhouete(silhouettes::AbstractVector,clusters::AbstractVector; plotargs...)
    sizels = length.(clusters)
    ls = accumulate(+,sizels)
    pushfirst!(ls,0)
    fig = plot(ls[1]+1:ls[2], sort(silhouettes[clusters[1]],rev=true);
        # seriestype = :steppre,
        fillcolor = 1, 
        fillrange = 0,
        legend=false,
        plotargs...
    )
    for i in 2:length(clusters)
        plot!(fig, ls[i]+1:ls[i+1], sort(silhouettes[clusters[i]],rev=true),
            fillcolor = i, 
            fillrange = 0,
            legend=false
        )
    end
    avgsilh = round(mean(silhouettes),digits=4)
    plot!(fig; title=string("k=",length(sizels), ", <silhouette>=",avgsilh))
    # println(plotargs)
    # plot!(fig; plotargs)
end


