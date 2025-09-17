
function plotmissing(magicgeno::MagicGeno; 
    targetfounder::Bool = true,
    color = Plots.cgrad([:aqua, :blue], [0,1]),    
    boundaryline = (1.5, :gray),    
    annotefontsize=12,     
    outfile::Union{Nothing, AbstractString}=nothing,
    plotkeyargs...)
    nchr = length(magicgeno.markermap)
    missls = [begin 
        if targetfounder
            genomtx = magicgeno.foundergeno[chr]
            formatvec = magicgeno.markermap[chr][!,:founderformat]
        else
            genomtx = magicgeno.offspringgeno[chr]
            formatvec = magicgeno.markermap[chr][!,:offspringformat]
        end
        MagicBase.getismissing(genomtx, formatvec)
    end for chr in 1:nchr]
    heat = permutedims(reduce(vcat,missls))
    poscmls = [i[!,:poscm] for i in magicgeno.markermap]
    chrlenls = accumulate(+,last.(poscmls))
    pushfirst!(chrlenls,0.0)
    pop!(chrlenls)    
    xposls = map((x,y)->x .+ y, poscmls, chrlenls)    
    xls = reduce(vcat,xposls)
    nrow = size(heat,1)
    @info string("average ", targetfounder ? "founder" : "offspring", " missing = ", mean(heat))
    missmap=heatmap(xls, 1:nrow, heat; 
        color, 
        size = (1500,max(500,nrow*10)), 
        dpi = 300,
        left_margin=20Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=15Plots.mm,
        grid = false,        
        legend = false,            
    )
    scatter!(xls, -0.5*ones(length(xls)), marker=(:vline,6,:black))
    x0 = last.(xposls)
    pushfirst!(x0, xposls[1][1])
    xy=hcat([[[i,0] [i,nrow+0.75] [NaN,NaN]] for i=x0]...)
    chrlabx=xy[1,1:3:end]
    chrlabx=(chrlabx[1:end-1] .+ chrlabx[2:end]) ./ 2
    ypos = nrow + 1.5 
    # ypos = nrow + 1.0
    chridls= [i[1,:linkagegroup] for i=magicgeno.markermap]    
    annotate!(missmap,[(chrlabx[i],ypos,Plots.text(chridls[i],annotefontsize,:black,rotation=15))
        for i=1:length(chrlabx)])    
    indls = targetfounder ? magicgeno.magicped.founderinfo[!,:individual] : magicgeno.magicped.offspringinfo[!,:individual]
    ytick = (collect(1:nrow),indls)
    xmax = last(xls)
    plot!(missmap,xy[1,:],xy[2,:];
        line = boundaryline,
        legend = false,
        yticks = ytick,
        xlabel="Position (cM)",    
        yrange = (-2,nrow+2),
        xrange = (-0.01*xmax, xmax*1.01),
        plotkeyargs...        
    )
    isnothing(outfile) || savefig(missmap,outfile)
    missmap
end
