
function plotrecomheat(ldfile::Union{Nothing,AbstractString}, 
    linkagefile::AbstractString, mapfile::AbstractString;       
    iseachlg::Bool = false, 
    snpthin::Integer=1,
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd(),    
    outstem::AbstractString="outstem",
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    if iseachlg
        try 
            tused = @elapsed begin                                  
                mapfile2 = getabsfile(workdir,mapfile)      
                if !isnothing(ldfile)
                    ldfile2 = getabsfile(workdir,ldfile)      
                    plotrecomheat(ldfile2,mapfile2; 
                        ispairwiseld=true, iseachlg = true, workdir=figdir, missingstring, outstem)                            
                end            
                linkagefile2 = getabsfile(workdir,linkagefile)      
                plotrecomheat(linkagefile2,mapfile2;
                    ispairwiseld=false, iseachlg = true, workdir=figdir, missingstring, outstem)                            
                
            end        
            MagicBase.create_tar(figdir, getabsfile(workdir,tarplotfile))                            
            msg = string("LG heatmap in ",tarplotfile, ", tused = ",round(tused,digits=1), "s")
            printconsole(io,verbose,msg)  
        catch err
            @warn string("Could not plot heatmap. ", err)
        end
        try 
            rm(figdir;recursive=true)
        catch
            @warn string("Could not remove ",figdir)
        end
        tarplotfile
    else
        if !isnothing(ldfile)
            tused = @elapsed outfile = plotrecomheat(ldfile,mapfile;
                ispairwiseld=true,  iseachlg = false, snpthin,workdir, outstem)
            msg = string("LD heatmap in ", outfile, ", tused = ",round(tused,digits=1), "s")
            printconsole(io,verbose,msg)  
        end
        tused = @elapsed outfile = plotrecomheat(linkagefile,mapfile;
            ispairwiseld=false,  iseachlg = false, snpthin, workdir, outstem)
        msg = string("Linkage heatmap in ", outfile, ", tused = ",round(tused,digits=1), "s")
        printconsole(io,verbose,msg)          
    end
end

function plotrecomheat(pairwisefile::AbstractString, mapfile::AbstractString;
    ispairwiseld::Bool=false,    
    iseachlg::Bool = false,
    snpthin::Integer=1,
    color = Plots.cgrad([:aqua, :blue, :pink,:red], ispairwiseld ? [0,0.2,0.4,0.6] : [0,0.4,0.8,0.9]),
    boundaryline = (0.2,:dot,:gray),
    commentstring::AbstractString = "##", 
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd(),    
    outstem::AbstractString="outstem",
    plotkeyargs...)
     # read mapfile
     markercol = 1
     chrcol = 2
     mapfile2 = getabsfile(workdir,mapfile)
     mapdf = CSV.read(mapfile2,DataFrame; delim=",", missingstring, comment=commentstring)    
     missrows = findall(ismissing.(mapdf[!,chrcol]))
     deleteat!(mapdf, missrows)
     select!(mapdf,1:2) 
     if !iseachlg && snpthin > 1
        mapdf = mapdf[1:snpthin:end,:]
     end
    # read pairwisefile
    if ispairwiseld
        # input pairwisefile is a pairwise LD file
        inds, markers, dupebindict, recomnonfrac,recomlod = read_ld(pairwisefile; workdir)        
        outid = "LD_r2"        
    else
        # input pairwisefile is a pairwise linkage file
        inds, markers, physmapdict, nmissingls, recomnonfrac,recomlod = read_linkage(pairwisefile; workdir)
        outid = "1-recomfrac"
    end      
    dict = Dict(markers .=> 1:length(markers))   
    # split map by LG
    gmapls = groupby(mapdf,chrcol)
    iils = [begin 
        ii = [get(dict,i,nothing) for i=gmap[!,markercol]]
        b = isnothing.(ii)
        any(b) && (ii = ii[.!b])
        ii
    end for gmap in gmapls]
    chridls = [MagicBase.reset_chrid(string(gmap[1,chrcol])) for gmap in gmapls]
    outid = ispairwiseld ? "LD" : "linkage"
    if iseachlg        
        outfiles = []        
        for g in eachindex(iils)                        
            chrid = chridls[g]                
            ii = iils[g]   
            for g2 in 1:g
                chrid2 = chridls[g2]                
                ii2 = iils[g2]  
                lglod = recomlod[ii,ii2]
                nzlod =  nonzeros(lglod)
                isempty(nzlod) && continue
                minlodsave = round(min(unique(nzlod)...),digits=2)
                lgnonfrac = recomnonfrac[ii,ii2]
                nnonzero = nnz(lgnonfrac) 
                nnonzero == 0 && continue
                sim = Matrix(lgnonfrac)
                title = string(ispairwiseld ? "LD_r2" : "1-recomfrac"," of ",chrid)     
                g2 < g && (title *= string(" ~ ",chrid2))
                fnz = round(nnonzero/length(lgnonfrac),digits=4)
                title *= string(", freq_nz=",fnz, ", minlod=",minlodsave)
                heatmap(sim; color,title, 
                    xlabel = string("Marker index of ", chrid),
                    ylabel = string("Marker index of ", chrid2),                
                    size = (800,500),
                    dpi = 200,
                    left_margin=20Plots.mm,         
                    right_margin=15Plots.mm,         
                    top_margin=15Plots.mm,
                    bottom_margin=15Plots.mm,
                )                        
                if g2 == g
                    outfile = string(outstem,"_", outid,"_within_", chrid,".png")
                else
                    outfile = string(outstem,"_", outid,"_between_", chrid, "_", chrid2,".png")                    
                end            
                savefig(getabsfile(workdir,outfile))        
                push!(outfiles,outfile)
            end
        end
        outfiles
    else
        sim = Matrix{Float16}(recomnonfrac)
        recomnonfrac = recomlod = nothing
        GC.gc()                
        iiall = reduce(vcat,iils)
        sim = sim[iiall, iiall]
        nsnp = size(sim,1)                
        simmap = heatmap(sim; color,
            size = (1200,1200),
            dpi = 600, 
            xrange = (-nsnp/19,nsnp),
            yrange = (-nsnp/19,nsnp),
            # showaxis = false,
            tick = :none,
            border = :none,            
            left_margin=20Plots.mm,         
            right_margin=15Plots.mm,         
            top_margin=15Plots.mm,
            bottom_margin=15Plots.mm,
            plotkeyargs...
        )         
        chrlenls = accumulate(+,length.(iils))
        pushfirst!(chrlenls,0)
        chrposls = (chrlenls[1:end-1] .+ chrlenls[2:end]) ./ 2
        annotate!(simmap,[(chrposls[i], -nsnp/40, text(chridls[i],10,rotation=15)) for i in eachindex(chrposls)])
        annotate!(simmap,[(-nsnp/30, chrposls[i], text(chridls[i],10)) for i in eachindex(chrposls)])
        xy=hcat([[[i,0] [i,nsnp] [NaN,NaN]] for i=chrlenls]...)'
        plot!(simmap,xy[:,1],xy[:,2];
            line = boundaryline,
            legend=false,
        )
        plot!(simmap,xy[:,2],xy[:,1];
            line = boundaryline,
            legend=false,
        )
        thinid = snpthin==1 ? "" : string("_thin", snpthin)
        outfile = string(outstem,"_", outid, thinid, "_heatmap.png")
        savefig(simmap,getabsfile(workdir,outfile))
        outfile        
    end
end
