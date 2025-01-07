
function plot_manhattan(profile::DataFrame;
    islog10p::Bool=false,
    thresholds::Union{Nothing,AbstractVector} = nothing,
    line_threshold = (1.5,:dash,:black),
    line_boundary = (1.5,:dash,:gray),
    plotkeyargs...)    
    res = groupby(profile,:linkagegroup)
    nchr = length(res)
    xls = [i[:,:poscm] for i in res]
    for chr in 2:nchr
        xls[chr] .+= xls[chr-1][end]-xls[chr][1]
    end    
    col = islog10p ? :neg_log10p : :lod
    xmax = xls[end][end]
    ymin = min(profile[!,col]...)
    ymax = max(profile[!,col]...)
    ymax = isnothing(thresholds) ? 1.05*ymax : 1.05*max(ymax,thresholds...)
    chridls=[i[1,:linkagegroup] for i in res]
    @info chridls
    fig = plot(ylims = [0,ymax])
    if !isnothing(thresholds) 
        for threshold in thresholds
            plot!(fig,[0,xmax],[threshold,threshold]; 
                line=line_threshold)
        end
    end
    for chr=1:nchr
        plot!(fig,xls[chr],res[chr][!,col];   
            legend=false,
            linewidth = 2, 
        )
    end             
    if !isnothing(line_boundary)
        x0 = last.(xls)
        pushfirst!(x0, xls[1][1])
        xy=hcat([[[i,ymin] [i,ymax] [NaN,NaN]] for i=x0]...)
        plot!(xy[1,:],xy[2,:];line = line_boundary,legend=false)
    end
    ylabel = islog10p ? "-log₁₀P" : "LOD"
    plot!(fig; 
        xlabel="linkagegroup",
        ylabel,                        
        bottom_margin=5Plots.mm,         
        xticks = (mean.(xls),string.(chridls)),
        xrotation=15,
        plotkeyargs...       
    )
end


function plot_qqpval(profile::DataFrame)
    obsp = profile[!,:neg_log10p]
    expp= -log10.(rand(length(obsp)))
    xmin,xmax = extrema(expp)
    fig = qqplot(expp,obsp; 
        xlabel="expected -log₁₀P",ylabel="obsereved -log₁₀P")
    plot!(fig,[xmin,xmax],[xmin,xmax])
end


function getpeak(profile::DataFrame;
    CI_probs::AbstractVector = [0.99,0.95,0.9])
    ii = argmax(profile[!,:lod])
    chr = profile[ii,:linkagegroup]
    chrprofile = filter(row->row.linkagegroup==chr,profile)
    ww = diff(chrprofile[!,:poscm])
    ww = vcat(ww[1]/2, (ww[1:end-1] + ww[2:end]) ./ 2,ww[end]/2)
    lod = chrprofile[!,:lod]
    den = normalize(ww .* exp10.(lod),1)
    den = cumsum(den)
    cils = [begin 
        alpha = 1- CI_prob
        lb = findlast(<=(alpha/2.0),den)
        isnothing(lb) && (lb = 1)
        ub = findfirst(>=(1.0-alpha/2.0),den)        
        isnothing(ub) && (ub = length(den))
        # println("[CI_prob,alpha,lb,ub]=",[CI_prob,alpha,lb,ub])
        [max(lb-1,1),min(ub+1,length(den))]
    end for CI_prob in CI_probs]
    resnames = string.(names(chrprofile)[[1:3; 5:6]])
    append!(resnames,[string("CI",i) for i in CI_probs])
    append!(resnames,[string("CI",i,"_marker") for i in CI_probs])
    ii = argmax(chrprofile[!,:lod])
    resvals = Vector{Any}(chrprofile[ii,[1:3; 5:6]])
    resvals[4:5] .= round.(resvals[4:5],digits=4)
    append!(resvals,[round.(chrprofile[ci,:poscm],digits=4) for ci in cils])
    append!(resvals,[chrprofile[ci,:marker] for ci in cils])    
    resvals2 = [isa(i,Vector) ? join(i," ~ ") : string(i) for i in resvals]
    DataFrame(:id=>resnames,:peak=>resvals2)    
end


