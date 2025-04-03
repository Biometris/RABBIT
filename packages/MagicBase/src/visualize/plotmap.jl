

function plotmarkermap(mapxfile::AbstractString,mapyfile::AbstractString,
    moremapfiles::Vararg{AbstractString};
    boundaryline = (1.5,:dot,:gray),    
    isannotate::Bool=true,
    annotatefontsize::Integer=9,
    cordigits::Integer=2,    
    maplabels::Union{Nothing,AbstractVector}=nothing,
    isphysmap::Union{Nothing,AbstractVector}=nothing,
    missingstring=["NA","missing"],
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd(),
    outstem::Union{Nothing, AbstractString} = "outstem",
    io::Union{Nothing, IO}=nothing, 
    verbose::Bool=true,
    plotkeyargs...)
    filels = [mapxfile,mapyfile, moremapfiles...]
    filels = [getabsfile(workdir,mapfile) for mapfile = filels]
    b = isfile.(filels)
    all(b) || error(string("files do not exist: ",filels[.!b]))
    mapls = [readmarkermap(mapfile; del_ungrouped = false,
        commentstring,missingstring, workdir) for mapfile = filels]
    labells = [string("map",i, " position") for i=1:length(mapls)]
    if !isnothing(maplabels)
        n = min(length(maplabels),length(mapls))
        labells[1:n] .= maplabels[1:n]
    end
    n = length(mapls)
    if isnothing(isphysmap) 
        isphysmap = falses(n)        
    else
        length(isphysmap) == length(filels) || @error string("inconsistent number of maps")
    end
    poscol1 = first(isphysmap) ? :physposbp : :poscm
    b = ismissing.(mapls[1][!,poscol1])
    all(b) && return nothing
    nmarker_map1 = size(mapls[1],1)
    markersize = in(:markersize, keys(plotkeyargs)) ? plotkeyargs[:markersize] : (nmarker_map1 < 1000 ? 3.0 : (nmarker_map1 < 10000 ? 2.0 : 1.5))
    res = [plotmarkermap(mapls[1],mapls[j];
        boundaryline,
        markersize, 
        isannotate,annotatefontsize, cordigits,
        isphysmap = isphysmap[[1,j]], 
        maplabels=labells[[1,j]], io, verbose,plotkeyargs...) for j=2:n]    
    res = res[.!isnothing.(res)]
    isempty(res) && return nothing
    chrcol1 = first(isphysmap) ? :physchrom : :linkagegroup
    nchr = length(unique(mapls[1][!,chrcol1]))    
    sizewidth  = max(1200,nchr*60)
    plotsize = in(:size, keys(plotkeyargs)) ? plotkeyargs[:size] : (sizewidth,(n-1)*0.6*sizewidth)
    left_margin = in(:left_margin, keys(plotkeyargs)) ? plotkeyargs[:left_margin] : max(15,5*n)*Plots.mm
    bottom_margin = in(:bottom_margin, keys(plotkeyargs)) ? plotkeyargs[:bottom_margin] : 15*Plots.mm
    fig = plot(res...;
        layout=(n-1,1), 
        left_margin,
        bottom_margin,        
        size=plotsize)
    if !isnothing(outstem)
        outfile = outstem*".png"    
        try 
            savefig(fig, getabsfile(workdir, outfile))
            msg = string("mapcompare in ",outfile)
            printconsole(io,verbose, msg)
        catch err
            msg = string("Could not savefig. ", err)    
        end    
    end
    fig
end


"""
    plotmarkermap(mapx, mapy;
        boundaryline = (1.5,:dot,:black),
        markersize = 1.5,
        isannotate= true,
        maplabels=["mapx(cM)", "mapy(cM)"],
        isphysmap = [false,false],
        plotkeyargs...
    )

plot postions of mapx vs those of mapy.

# Positional arguments

`mapx::Vector{DataFrame}`: marker map for all chromosomes.

`mapy::Vector{DataFrame}`: comparing map for all chromosomes.

# Keyword arguments

`boundaryline = (1.0,:dot,:gray)`: vertical lines for chromosome boundaries.

`markersize::Real= size(mapx,1) <= 1000 ? 3.0 : 1.5`: size of scatter markers.

`isannotate::Bool=true`: if ture,  annotate chromosome ID and kendall correlation.

`mmaplabels::Union{Nothing,AbstractString}=nothing`: labels of comparing marker maps.

`isphysmap::AbstractVector = falses(2)`: specify if mapx and/or map are physical maps. 

`plotkeyargs...`: other Plots.plot keyward arguments. 

"""
function plotmarkermap(mapx::DataFrame,mapy::DataFrame;
    boundaryline = (1.0,:dot,:gray),
    markersize::Real= size(mapx,1) <= 1000 ? 3.0 : 1.5, 
    isannotate::Bool=true,
    annotatefontsize::Integer=9,
    cordigits::Integer=3,
    maplabels::AbstractVector=["mapx(cM)", "mapy(cM)"],
    isphysmap::AbstractVector = falses(2),     
    io::Union{Nothing, IO}=nothing, 
    verbose::Bool=true,
    plotkeyargs...)
    chrcolls =  [i ? :physchrom : :linkagegroup for i in isphysmap]
    poscolls =  [i ? :physposbp : :poscm for i in isphysmap]
    # remove missing in mapx    
    b = ismissing.(mapx[!,chrcolls[1]]) .|| ismissing.(mapx[!,poscolls[1]])
    if any(b) 
        mapx = mapx[.!b,:]
    else
        mapx = copy(mapx)
    end
    isempty(mapx) && return nothing
    # sort groups of physmap by chrid    
    if first(isphysmap)
        physchrls = unique(mapx[!,:physchrom])
        if !all(isa.(physchrls,Integer)) 
            physchrls2 = string.(physchrls)
            ischrint = all(occursin.(r"^[1-9][0-9]{0,6}$",physchrls2))
            if ischrint
                chrdict = Dict(physchrls2 .=> parse.(Int, physchrls2))
                mapx[!,:physchrom] .= [chrdict[i] for i in mapx[!,:physchrom]]          
                gmapx = groupby(mapx,:physchrom; sort=true)
                mapx = reduce(vcat, gmapx)  
            end
        end
    end
    gmapx = groupby(mapx,chrcolls[1]; sort=false)
    isorderedls = get_isorderedls(gmapx; isphysmap = isphysmap[1])
    chridls = [i[1,chrcolls[1]] for i in gmapx]
    badchrs = chridls[ismissing.(isorderedls)]
    if !isempty(badchrs)
        msg = string(maplabels[1], " ", chrcolls[1], " with missing positoins:  ", badchrs)
        printconsole(io, false,"Warning: "*msg)   
        @warn msg
    end
    if !all(skipmissing(isorderedls))        
        badchrs = chridls[isorderedls .=== false]
        if !isempty(badchrs)
            msg = string(maplabels[1], " ", chrcolls[1], " with unsorted positions: ", badchrs)
            printconsole(io, false,"Warning: "*msg)   
            @warn msg
        end        
        mapx = copy(mapx)
        sort!(mapx,[chrcolls[1],poscolls[1]])        
    end
    # seperate missing in mapy      
    b = ismissing.(mapy[!,chrcolls[2]]) .|| ismissing.(mapy[!,poscolls[2]])
    if any(b) 
        mapy_miss = mapy[b,:]
        mapy = mapy[.!b,:]
    else
        mapy_miss = mapy[[],:]
        mapy = copy(mapy)
    end    
    # keep only common snps
    commsnps = intersect(mapx[!,:marker],mapy[!,:marker])
    if isempty(commsnps) 
        msg = string("no common markers")
        @error msg
        printconsole(io, false,msg)
        return nothing
    end
    b = [in(i, commsnps) for i in mapx[!,:marker]]
    all(b) || deleteat!(mapx,.!b)
    isempty(mapx) && return nothing
    b = [in(i, commsnps) for i in mapy[!,:marker]]
    all(b) || deleteat!(mapy,.!b)  
    isempty(mapy) && return nothing          
    # reorder LGs of mapy to match mapx
    snprule = Dict(mapx[!,:marker] .=> 1:size(mapx,1))
    gmapy = groupby(mapy,chrcolls[2]; sort=false)
    lgpos = [median([snprule[i] for i in df[!,:marker]]) for df in gmapy]
    ii = sortperm(lgpos)
    mapy = reduce(vcat,gmapy[ii])    
    # reset poscol
    reset_poscol!(mapx,poscolls[1])
    reset_poscol!(mapy,poscolls[2])        
    gmapy = groupby(mapy,chrcolls[2]; sort=false)
    snprule = Dict(mapy[!,:marker] .=> 1:size(mapy,1))
    gmapx = groupby(mapx,chrcolls[1]; sort=false)
    orderacc = round.(mapcorkendall(gmapx, gmapy; isphysmap), digits = cordigits)
    msg = string("tau=",join(abs.(orderacc),"|"))
    printconsole(io,verbose,msg)
    res=[begin    
        tau = orderacc[lg]
        df = gmapx[lg]
        snpii = [get(snprule,i,nothing) for i in df[!,:marker]]        
        b = .!isnothing.(snpii)
        snpii2 = snpii[b]
        xls = Float32.(df[b,poscolls[1]])        
        yls = Float32.(mapy[snpii2,poscolls[2]])         

        b = [in(i, mapy_miss[!,:marker]) for i in df[!,:marker]]
        xls_ungroup = Vector{Float32}(df[b, poscolls[1]])       
        if tau < 0 
            tau = -tau
            xbegin = Float32(df[begin,poscolls[1]])
            xend = Float32(df[end,poscolls[1]])
            xls .= (xbegin + xend)  .- reverse(xls)
            yls .= reverse(yls)            
            xls_ungroup .= (xbegin + xend) .- reverse(xls_ungroup)
        end         
        [xls,yls,tau,xls_ungroup]
    end for lg in eachindex(orderacc)]    
    # plot chr boundaries
    xmax = mapx[end,poscolls[1]]
    ymax = mapy[end,poscolls[2]]
    xminls = [i[1,poscolls[1]] for i in gmapx]    
    yminls = [i[1,poscolls[2]] for i in gmapy]
    xminls[1] = 0.0
    yminls[1] = 0.0
    push!(xminls, xmax)
    push!(yminls, ymax)
    nchr = length(res)
    fig = plot([0,0]; 
        size = (max(1200,nchr*65),max(1200,nchr*65)),
        dpi = 600,     
        xrange = (-xmax*0.08,xmax*1.05),
        yrange = (-ymax*0.12,ymax*1.05),
        xlabel = maplabels[1],
        ylabel = maplabels[2],
        # xtick = round.(Int,xminls),
        # ytick = round.(Int,yminls),
        legend = false,
        grid = false,        
        plotkeyargs...,
        left_margin=20Plots.mm,         
        bottom_margin=10Plots.mm,
    )
    xy=hcat([[[i,0] [i,ymax] [NaN,NaN]] for i=xminls]...)'
    plot!(fig,xy[:,1],xy[:,2];
        line = boundaryline,
        legend=false,
    )
    chrposls = (xminls[1:end-1] .+ xminls[2:end]) ./ 2    
    chridls = [i[1,chrcolls[1]] for i in gmapx]
    annotate!(fig,[(chrposls[i], ymax*(-0.04), text(chridls[i],annotatefontsize,rotation=30)) for i in eachindex(chrposls)])
    xy=hcat([[[0,i] [xmax, i] [NaN,NaN]] for i=yminls]...)'
    plot!(fig,xy[:,1],xy[:,2];
        line = boundaryline,
        legend=false,
    )
    chrposls = (yminls[1:end-1] .+ yminls[2:end]) ./ 2
    chridls = [string(i[1,chrcolls[2]]) for i in gmapy]
    chridlen = maximum(length.(chridls))
    chridls2 = [lpad(i,chridlen," ") for i in chridls]
    annotate!(fig,[(xmax*(-0.04),chrposls[i], text(chridls2[i],annotatefontsize)) for i in eachindex(chrposls)])
    # scatter        
    marker_colors = distinguishable_colors(nchr, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    for chr=1:nchr
        xls,yls, tau,xls_ungroup = res[chr]
        c = marker_colors[chr]
        scatter!(fig,xls, yls, 
            legend=false,
            marker=(:circle, markersize, 0.5, c,stroke(c)),            
        )
        yls_ungroup = similar(xls_ungroup)
        # yls_ungroup .= rand(Uniform(0,ymax*0.02))
        yls_ungroup .= ymax*(-0.01f0)
        scatter!(fig,xls_ungroup, yls_ungroup, 
            legend=false,
            marker=(:vline, markersize, 0.5, c,stroke(c)),
        )
        if isannotate
            labx=mean(xls[[1,end]])        
            lab0 = chr==1 ? "τ=" : ""
            lab = string(lab0,tau)
            annotate!(fig,[(labx,ymax*1.055,text(lab,:top,annotatefontsize,rotation=30))])
        end
    end
    fig
end

function reset_poscol!(mapdf::AbstractDataFrame,poscol)
    if poscol == :physposbp
        mapdf[!,:physposbp] .*= 1e-6
        gmapdf = groupby(mapdf, :physchrom)
    else
        gmapdf = groupby(mapdf, :linkagegroup)
        for df in gmapdf
            if df[1,:poscm] > 0
                df[!,:poscm] .-= df[1,:poscm] 
            end
        end
    end
    df = gmapdf[1]
    posmin, posmax = extrema(df[!,poscol])
    genolen = posmax
    for g in 2:length(gmapdf)
        df = gmapdf[g]
        posmin, posmax = extrema(df[!,poscol])
        if posmin > genolen 
            genolen = posmax
        else
            df[!,poscol] .+= genolen     
            genolen += posmax
        end        
    end
    mapdf
end

function reset_chrid(chrid::AbstractString)
    if occursin("chr",lowercase(chrid)) 
        replace(chrid,"chr"=>"", "Chr"=>"")
    else
        if occursin("lg",lowercase(chrid))             
            replace(chrid,"lg"=>"", "LG"=>"")
        else
            occursin(r"^[1-9][0-9]{0,6}$",chrid) ? string("LG",chrid) : chrid
        end
    end 
end
