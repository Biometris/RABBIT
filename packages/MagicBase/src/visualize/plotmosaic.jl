
function savemosaic(contfgl::MagicGeno;    
    offspring::AbstractVector=1:3,
    fglcolors::Union{Nothing, AbstractVector}=nothing,
    workdir::AbstractString = pwd(),
    outstem::AbstractString = "outstem",
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    tused = @elapsed begin         
        _, mosaicfile = plotmosaic(contfgl; offspring,fglcolors, workdir,outstem)                    
    end
    msg = string("save in ",mosaicfile, ", tused=",round(tused,digits=1),"s")
    printconsole(io,verbose,msg)
    tused = @elapsed begin         
        _, recomfiles = hist_mosaic(contfgl; workdir,outstem)
    end 
    msg = string("save in ",join(recomfiles,", "), ", tused=",round(tused,digits=1),"s")
    printconsole(io,verbose,msg)
    (mosaicfile,recomfiles...)
end

function plotmosaic(contfgl::MagicGeno;
    offspring::AbstractVector=1:3,
    fglcolors::Union{Nothing, AbstractVector}=nothing,
    workdir::AbstractString = pwd(),
    outstem::AbstractString = "outstem")
    fgls = vec([i[1][2] for i in contfgl.foundergeno[1]])
    fgls = unique(fgls)
    nfgl = length(fgls)
    if isnothing(fglcolors) || length(fglcolors)<nfgl
        if !isnothing(fglcolors)
            length(fglcolors)<nfgl && @error string("length of fglcolors = ", length(fglcolors), " < nfgl, distinguishable_colors will be used.")
        end
        fglcolors = distinguishable_colors(nfgl, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    end
    nind = size(contfgl.magicped.offspringinfo[!,:individual],1)
    indls = offspring[offspring .<= nind]
    nind = length(indls)
    chridls = [i[1,:linkagegroup] for i in contfgl.markermap]
    chrlenls = [i[end,:poscm] - i[1,:poscm] for i in contfgl.markermap]
    if nind > 600
        chrlenls .*= 600/nind
    end
    nchr = length(chrlenls)
    chrgap = mean(chrlenls)*0.1
    chrstartxls = accumulate(+,chrlenls)
    pop!(chrstartxls)
    pushfirst!(chrstartxls,0.0)
    chrstartxls .+= chrgap .* collect(1:nchr)    
    homologheight = 0.8*chrgap
    homologgap = 0.1*chrgap
    indheight = 2*homologheight + homologgap # diploid, two homologs
    indgap = 3*homologgap
    indstartyls = indheight * collect(0:nind-1)
    indstartyls .+= indgap .* collect(1:nind)
    # plot 
    xmax = (sum(chrlenls) + nchr*chrgap)+chrgap
    ymax = nind*(indheight+indgap)+indgap    
    if ymax/xmax < 0.25 
        yscale = 0.25*xmax/ymax
        ymax *= yscale
        homologheight *= yscale
        homologgap *= yscale
        indheight *= yscale
        indstartyls .*= yscale
    end
    # println("(xmax,ymax)=",(xmax,ymax))
    mosaic = plot(size=(1.2*xmax, 1.2*ymax),
        xlims=(-2*chrgap,xmax),
        ylims=(-indheight,ymax),     
        left_margin=20Plots.mm,         
        bottom_margin=15Plots.mm,   
    )    
    for chr in 1:nchr
        chrstartx = chrstartxls[chr]
        for i in eachindex(indls)
            ind = indls[i]
            indfgl = contfgl.offspringgeno[chr][:,ind]
            for j in 1:2
                mosaic_homolog!(mosaic,indfgl[j];
                    fglcolors, chrstartx,
                    indstarty = indstartyls[i], 
                    homologheight, homologgap,firsthomolog = j==1
                )
            end
        end
    end
    fontsize = nind<600 ? 12 : max(3,round(Int,12*(600/nind)))
    y0 = indheight* (nind < 4 ? -0.3 : -1.0)
    for chr in 1:nchr
        annotate!(chrstartxls[chr]+chrlenls[chr]/2, y0, text(chridls[chr],fontsize, "courier"))
    end
    for i in eachindex(indls)
        annotate!(-chrgap, indstartyls[i]+indheight/2, text(string(indls[i]),fontsize,"courier"))
    end
    dpi = min(max(300,round(Int,nind/10)),2000)    
    plot!(showaxis = false,axis=nothing,dpi=dpi)
    mosaicfile = string(outstem,"_mosaic.png")    
    try 
        savefig(getabsfile(workdir,mosaicfile))
    catch err
        @error err
        @error string("Cannot savefig ", mosaicfile)
    end
    mosaic, mosaicfile
end

function hist_mosaic(contfgl::MagicGeno;
    workdir::AbstractString = pwd(),
    outstem::AbstractString="outstem")    
    outfiles = []
    ginbreed,inbreedcoef = hist_inbreedcoef(contfgl)
    gnrecom,nrecom = hist_nrecom(contfgl)
    gseglen1 = hist_seglen(contfgl; ishomologpair=false)    
    gseglen2 = hist_seglen(contfgl; ishomologpair=true)    
    gden1,recomden1 = hist_recomden(contfgl;ishomologpair=false)
    gden2,recomden2 = hist_recomden(contfgl;ishomologpair=true)
    gnone = plot(axis=nothing, showaxis=false)
    ghist = plot(gnrecom,ginbreed,gnone,gnone,
        gden1, gseglen1, gnone,gnone,gden2, gseglen2,
        layout=grid(5,2,
        heights=(0.33,0,0.33,0,0.33)),
        size=(900,900),
        left_margin=10Plots.mm,         
        bottom_margin=10Plots.mm,
    )    
    outfile = outstem*"_recom.png"
    savefig(ghist, getabsfile(workdir,outfile))
    push!(outfiles, outfile)
    # save data
    offls = contfgl.magicped.offspringinfo[!,:individual]
    memls = contfgl.magicped.offspringinfo[!,:member]
    df = DataFrame("individual"=>offls, "member"=>memls,
        "inbreedcoef"=>inbreedcoef, "recomnum"=>nrecom)
    outfile = outstem*"_recom.csv"
    CSV.write(getabsfile(workdir,outfile), df)      
    push!(outfiles,outfile)  
    ghist, outfiles
end

# define a function that returns a Plots.Shape
# (x,y) coordinates of left-bottom
rectangle(x,y,w, h) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

function mosaic_homolog!(mosaic::Plots.Plot,
    homologfgl::AbstractVector; 
    fglcolors::AbstractVector,
    chrstartx::Real, 
    indstarty::Real, 
    homologheight::Real,
    homologgap::Real,
    firsthomolog::Bool)
    y0 = firsthomolog ? indstarty : indstarty + homologheight + homologgap
    h = homologheight
    for s in 1:length(homologfgl)-1
        x0 = homologfgl[s][1] + chrstartx
        w = homologfgl[s+1][1] - homologfgl[s][1] 
        c = fglcolors[round(Int,homologfgl[s][2])]   
        # fill=(0,c)     
        plot!(mosaic, rectangle(x0,y0,w,h), linewidth=0, color=c,legend=false)
    end
end

function get_breakpoints(homolog_fgl::AbstractVector)
    first.(homolog_fgl)
end

function gridhomologfgl(homologfgl::Vector{Vector{T}} where T <: AbstractFloat,
    chrmarkerpos::AbstractVector)
    # @assert chrmarkerpos[end] < homologfgl[end][1]
    res = zeros(Int,length(chrmarkerpos))
    seg = 1
    i=1
    while i<=length(chrmarkerpos)
        if chrmarkerpos[i]<homologfgl[seg+1][1]
            res[i] = homologfgl[seg][2]
            i += 1
        else
            seg += 1
        end
    end
    # pos = findall(abs.(diff(res)) .> 0)
    # [chrmarkerpos[pos] chrmarkerpos[pos .+ 1]]
    # scatter(map((x,y)->(x,y), chrmarkerpos, res))
    res
end

function get_ibdlength(indfgl::AbstractVector)
    breaks = sort(unique(reduce(vcat,get_breakpoints.(indfgl))))
    posls = (breaks[1:end-1] .+ breaks[2:end]) ./ 2
    m = gridhomologfgl(indfgl[1],posls)
    p = gridhomologfgl(indfgl[2],posls)    
    sum(diff(breaks) .* (m .== p))
end

function get_inbreedcoef(contfgl::MagicGeno)
    ibdlen = sum([[get_ibdlength(indfgl) for indfgl in eachcol(chrfgl)] for chrfgl in contfgl.offspringgeno])
    genolen = sum([i[end,:poscm]-i[1,:poscm] for i in contfgl.markermap])
    ibdlen ./ genolen    
end

function hist_inbreedcoef(contfgl::MagicGeno)
    inbreedcoef = get_inbreedcoef(contfgl)
    avg = round(mean(inbreedcoef),digits=3)
    ci = round.(quantile(inbreedcoef,[0.025,0.975]),digits=3)
    title = string(avg, "(",ci[1],"~",ci[2],")")
    ghist = histogram(inbreedcoef;        
        title,
        titlefontsize=10,
        xlabel = "inbreedcoef coefficient",
        ylabel = "#offspring",
        legend = false
    )
    ghist, inbreedcoef
end

function hist_nrecom(contfgl::MagicGeno)
    nrecom2 = sum([begin     
        breaks = get_breakpoints.(a)    
        n2 = [length(unique(reduce(vcat,i))) for i in eachcol(breaks)] .- 2
        n2
    end for a in contfgl.offspringgeno])
    avg = round(mean(nrecom2),digits=2)
    ci = round.(quantile(nrecom2,[0.025,0.975]),digits=3)
    title = string(avg, "(",ci[1],"~",ci[2],")")
    ghist =histogram(nrecom2;        
        title,
        titlefontsize=10,
        xlabel = "#recombination breakpoints",
        ylabel = "#offspring",
        legend = false
    )
    ghist, nrecom2
end

function hist_seglen(contfgl::MagicGeno; ishomologpair::Bool=false)    
    seglenls = [begin 
        breaks = get_breakpoints.(a)  
        if ishomologpair
            seglen = reduce(vcat,[diff(sort(unique(reduce(vcat,i)))) for i in eachcol(breaks)])            
        else
            seglen = reduce(vcat,diff.(breaks))                
        end
    end for a in contfgl.offspringgeno]
    segls2 = reduce(vcat,seglenls)
    avg = round(mean(segls2),digits=2)
    ci = round.(quantile(segls2,[0.025,0.975]),digits=3)
    title = string(avg, "(",ci[1],"~",ci[2],")")
    xlabel = ishomologpair ? "homologous pair" : "homolog"
    xlabel *= " segement length (cM)"
    histogram(segls2;
        title,
        titlefontsize=10,
        xlabel,
        ylabel = "frequency",
        legend = false,
    )
end

function hist_recomden(contfgl::MagicGeno; ishomologpair::Bool=false)    
    recomden = [begin     
        a = contfgl.offspringgeno[chr]
        breaks = get_breakpoints.(a)
        if ishomologpair
            nrecom = [length(unique(reduce(vcat,i))) - 2 for i in eachcol(breaks)]
        else
            nrecom = reduce(vcat,length.(breaks) .- 2)
        end        
        chrlen = a[1,1][end][1]
        nrecom ./ (chrlen/100)
    end for chr in eachindex(contfgl.offspringgeno)]
    recomden2 = reduce(vcat,recomden)
    avg = round(mean(recomden2),digits=2)
    ci = round.(quantile(recomden2,[0.025,0.975]),digits=3)
    title = string(avg, "(",ci[1],"~",ci[2],")")
    xlabel = "#breakpoints per Morgan per"
    xlabel *= ishomologpair ? " homologous pair" : " homolog"
    ylabel = ishomologpair ? " #homologous pairs" : "#homologs"
    ghist = histogram(recomden2;        
        title,
        titlefontsize=10,
        xlabel,
        ylabel,
        legend = false
    )
    ghist, recomden2
end