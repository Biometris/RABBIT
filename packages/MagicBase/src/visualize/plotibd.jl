"""
    plotcondprob(magicancestry,offspring=nothing,probtype="haploprob",
        colorgradient = cgrad([:white,:blue,:red]),,
        boundaryline = (1.5, :gray),
        truemarker=(:star, 5, 0.5,:gray,stroke(:gray)),
        truefgl=nothing,
        outfile::Union{Nothing, AbstractString}=nothing,
        plotkeyargs...)

plot heatmap for conditional probability.

# Positional arguments

`magicancestry::MagicAncestry`: magicancestry returned from `magicreconstruct`.

# Keyword arguments

`probtype::AbstractString="haploprob"`: specify type of condprob

`offspring::Union{Nothing,Integer}=nothing`: offsprign index. By default, a random offspring index. 

`colorgradient::ColorGradient=cgrad([:white,:blue,:red]),`: color gradient
for heatmap

`boundaryline=(1.5,:gray)`: vertical lines for chromosome boundaries.

`truemarker=(:star, 5, 0.5,:gray,stroke(:gray))`: scatter markers for
true ancestral states.

`truefgl::Union{Nothing,MagicGeno}=nothing`: provides true ancestral origins. 

`outfile::Union{Nothing, AbstractString}=nothing`: if nothing, not save the plot, and otherwise save to outfile.  

`plotkeyargs...`: other Plots.plot keyward arguments. 

"""
function plotcondprob(magicancestry::MagicAncestry;
    offspring::Union{Nothing,Integer}=nothing,
    probtype::AbstractString="haploprob",
    colorgradient::ColorGradient=cgrad([:white,:blue,:red]),
    boundaryline = (1.5, :gray),
    truemarker=(:xcross, 6, 0.5,:gray,stroke(:gray)),
    annotefontsize=12, 
    truefgl::Union{Nothing,MagicGeno}=nothing,
    outfile::Union{Nothing, AbstractString}=nothing,
    plotkeyargs...)
    if !(probtype in ["haploprob", "genoprob", "diploprob"])
        msg = string("unknown probtype: ", probtype)
        msg  = msg*". probtype must be \"haploprob\",\"genoprob\", or \"diploprob\"."
        error(msg)
    end
    if !isnothing(truefgl)       
        # align markers        
        truesnps = [i[!,:marker] for i in truefgl.markermap]
        estsnps = [i[!,:marker] for i in magicancestry.markermap]
        if truesnps != estsnps
            MagicBase.alignchromosome!(truefgl,magicancestry);
            MagicBase.alignmarker!(truefgl,magicancestry);
        end    
    end    
    fglls = magicancestry.statespace["haplotype"]
    nfgl = length(fglls)
    if probtype == "haploprob"
        condprob = magicancestry.haploprob
        ylabel = "Founder"
        titlestem = "Haplotype probablility"
        statespace = fglls
    elseif probtype == "genoprob"
        condprob = magicancestry.genoprob
        ylabel = "Founder origin state"
        titlestem = "Genotype probablility"        
        statespace = ["("*join(fglls[i],",")*")" for i in MagicBase.prior_genoindex(nfgl)]
    elseif probtype == "diploprob"
        condprob = magicancestry.diploprob
        ylabel = "Founder origin state"
        titlestem = "Diplotype probablility"        
        statespace = ["("*join(fglls[i],",")*")" for i in MagicBase.prior_diploindex(nfgl)]
    else
        error("unknown saving probaiblity type: ", probtype)
    end
    isnothing(condprob) && error(string(probtype, " is nothing"))
    # calculate xposls
    poscmls = [i[!,:poscm] for i in magicancestry.markermap]
    chrlenls = accumulate(+,last.(poscmls))
    pushfirst!(chrlenls,0.0)
    pop!(chrlenls)    
    xposls = map((x,y)->x .+ y, poscmls, chrlenls)    
    xls = reduce(vcat,xposls)
    # calculate heat
    noff = length(condprob[1])
    if isnothing(offspring)        
        off = rand(1:noff)
    else
        offspring > noff && error(string("offspring index=",offspring, " > #offspring=",noff))
        off=offspring        
    end    
    offinfo = magicancestry.magicped.offspringinfo    
    if in("isoutlier",names(offinfo))
        isoutlier = ismissing(offinfo[off, :isoutlier]) ? false : offinfo[off, :isoutlier]
    else
        isoutlier = false
    end
    offid = string(offinfo[off,:individual], isoutlier ? ", too many #breakpoints" : "")
    # prob = [hcat([j[off,:] for j=i]...)' for i=condprob]
    prob = [i[off] for i=condprob]
    if !isnothing(truefgl)                 
        # 
        nfgl = length(magicancestry.statespace["haploindex"])             
        diploindex = MagicBase.prior_diploindex(nfgl)   
        truechrls = findtruelg(truefgl.markermap,magicancestry.markermap)        
        offtrue = []
        for chr in eachindex(truechrls)
            true_chr = truechrls[chr]
            estsnpidls = magicancestry.markermap[chr][!,:marker]
            true_ii = [begin
                truesnpidls = truefgl.markermap[c][!,:marker]
                if estsnpidls == truesnpidls
                    1:length(estsnpidls)
                else
                    d = setdiff(estsnpidls,truesnpidls)
                    if !isempty(d)
                        estchrid =  magicancestry.markermap[chr][1,:linkagegroup]
                        msg= string("markers of chr=",estchrid, " in magicancestry but not in truefgl: ",d)
                        @warn msg
                    end
                    rule = Dict(truesnpidls .=> 1:length(truesnpidls))
                    ii = [get(rule, i, missing) for i in estsnpidls]
                    collect(skipmissing(ii))
                end
            end for c in true_chr]
            off_fgl = reduce(vcat, [truefgl.offspringgeno[true_chr[c]][true_ii[c],off] for c in 1:length(true_chr)])
            # for each chromosome, adjust unidentifiable the absoluate phase
            if probtype == "diploprob" 
                offest = diploindex[[argmax(i) for i in eachrow(prob[chr])]]
                nerr0 = sum(off_fgl .!= offest)
                nerr1 = sum(reverse.(off_fgl) .!= offest)
                if nerr0 > nerr1
                    off_fgl .= reverse.(off_fgl)
                end
            end
            append!(offtrue,off_fgl)       
        end        
        # offtrue=reduce(vcat, [i[:,off] for i in truefgl.offspringgeno[chrii]])
        if probtype == "haploprob"
            xy = reduce(vcat,map((x,y)->[[x,i] for i in unique(y)], xls,offtrue))
            xls_haplo = first.(xy) # rewrite xls
            yls = last.(xy)            
        elseif probtype == "genoprob"
            index = MagicBase.prior_genoindex(nfgl)
            # dict=Dict([sort(string.(i)) for i=index] .=> 1:length(index))
            dict=Dict(sort.(index) .=> 1:length(index))
            yls =  [get(dict, sort(i), nothing) for i=offtrue]                         
        else
            probtype == "diploprob" || @error string("unexpected probtype=",probtype)
            index = MagicBase.prior_diploindex(nfgl)
            # dict=Dict([string.(i) for i=index] .=> 1:length(index))
            dict=Dict(index .=> 1:length(index))            
            yls =  [get(dict, i, nothing) for i=offtrue]                                    
        end        
    end    
    # heat= vcat(prob...)' ==> MethodError: no method matching zero(::Type{ColorTypes.RGBA{Float64}})
    heat= Matrix(vcat(prob...)')
    heat[isnan.(heat)] .= 0.0
    nzrows = findall(sum(heat,dims=2)[:,1] .> 0)
    if !isnothing(truefgl)
        nz_yls = setdiff(yls,nzrows)
        if !isempty(nz_yls)
            nzrows = sort(vcat(nzrows,nz_yls))
        end
    end
    if probtype == "diploprob"
        nfgl = length(magicancestry.statespace["haploindex"])        
        diplostates = MagicBase.prior_diploindex(nfgl)        
        dict=Dict(diplostates .=> 1:length(diplostates)) 
        nzrows_index = unique(vcat(diplostates[nzrows],reverse.(diplostates[nzrows])))
        nzrows = sort([dict[i] for i in nzrows_index])
    end    
    ytick = (collect(1:length(nzrows)),statespace[nzrows])    
    heat = heat[nzrows,:]
    isempty(heat) && @info string("empty heat matrix for offspring=",off)        
    ibdmap=heatmap(xls, 1:size(heat,1), heat,
        c=colorgradient, 
        size = (1500,max(500,length(nzrows)*(probtype == "haploprob" ? 40 : 30))), 
        dpi = 300,
        left_margin=20Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=15Plots.mm,
        grid = false,        
    )
    scatter!(xls, 0.08*ones(length(xls)), marker=(:vline,8,:black))
    # len = size.(prob,1)    
    # x0=accumulate(+,len[1:end])
    # pushfirst!(x0,0)
    x0 = last.(xposls)
    pushfirst!(x0, xposls[1][1])
    nstate= size(heat,1)
    xy=hcat([[[i,0] [i,nstate+0.75] [NaN,NaN]] for i=x0]...)
    chrlabx=xy[1,1:3:end]
    chrlabx=(chrlabx[1:end-1] .+ chrlabx[2:end]) ./ 2
    ypos = probtype == "haploprob" ? nstate + 1.5 : nstate+1.5
    # ypos = nstate + 1.0
    chridls= [i[1,:linkagegroup] for i=magicancestry.markermap]
    annotate!(ibdmap,[(chrlabx[i],ypos,Plots.text(chridls[i],annotefontsize,:black,rotation=15))
        for i=1:length(chrlabx)])
    plot!(ibdmap,xy[1,:],xy[2,:];
        line = boundaryline,
        legend=false,
        yticks = ytick,
        xlabel="Position (cM)",
        ylabel=ylabel,
        title=string(titlestem, " for ",off,"-th offspring =",offid),
        titlefont = font("Times", 12),
        plotkeyargs...
    )
    if !isnothing(truefgl)
        dict = Dict(nzrows .=> 1:length(nzrows))                
        yls2=  [get(dict, i, nothing) for i in yls]            
        xls2 = probtype == "haploprob" ? xls_haplo : xls
        (nothing in yls2) && error(string("nothing for offspring=",off))
        scatter!(ibdmap,xls2,yls2,marker=truemarker,legend=false)        
    end
    isnothing(outfile) || savefig(ibdmap,outfile)
    ibdmap
end


"""
    animcondprob(magicancestry,fps=1,outstem="",kewargs...)

animation for plots of conditional probability.

# Positional arguments

`magicancestry::MagicAncestry`: magicancestry returned from `magicreconstruct`.

# Keyword arguments

`fps::Real=1`: number of frames per seconds.

`outfile::AbstractString="condprob.gif"`: output file for saving animation.

see [`plotcondprob`](@ref) for keyargs.

"""
function animcondprob(magicancestry::MagicAncestry;
    fps::Real=1,
    probtype::AbstractString="haploprob",
    colorgradient::ColorGradient=cgrad([:white,:blue,:red]),
    boundaryline = (1.5, :gray),
    truemarker=(:xcross, 3, 0.5,:gray,stroke(:gray)),
    truefgl::Union{Nothing,MagicGeno}=nothing,
    outfile::Union{Nothing,AbstractString}=nothing)
    noff=size(magicancestry.magicped.offspringinfo,1)
    anim = Animation()
    n1 = ceil(Int,fps)
    n2 = fps == 0 ? 1 : ceil(Int,1/fps)
    for off=1:noff
        plotcondprob(magicancestry;
            offspring=off,
            probtype, colorgradient, boundaryline, truemarker, truefgl)
        for i=1:n2
            frame(anim)
        end
    end
    if isnothing(outfile)
        gif(anim,fps=n1)
    else
        gif(anim,outfile,fps=n1)
    end
end

"""
    plotrecombreak(magicancestry,chr=1,
        colorgradient = ColorGradient([:yellow,:blue,:red]),
        truefgl=nothing)

plot recombination breakpoints.

# Positional arguments

`magicancestry::MagicAncestry`: magicancestry returning from `magicreconstruct`.

# Keyword arguments

`chr::Integer=1`: chromosome index.

`colorgradient::ColorGradient=cgrad([:white,:blue,:red])`: color gradient
for heatmap

`truefgl::Union{Nothing,MagicGeno}`: provides true ancestral origins. 

"""
function plotrecombreak(magicancestry::MagicAncestry;
    chr::Integer=1,
    colorgradient::ColorGradient=cgrad([:white,:blue,:red]),
    truefgl::Union{Nothing,MagicGeno}=nothing)
    if isnothing(magicancestry.viterbipath)
        @error string("viterbipath is required")
        return -1
    end
    a=magicancestry.viterbipath[chr]'
    gmarkermap = magicancestry.markermap
    chrid = gmarkermap[chr][1,:linkagegroup]
    g1=heatmap(a,c=colorgradient,
        xlabel=string("SNP index in chr=", chrid),
        ylabel="Offspring index")
    estrecom = calnumrecom(magicancestry.viterbipath)
    g2=histogram(estrecom,xlabel="Estimated #recom in all chrs",
        ylabel="Frequency",legend=false)
    if isnothing(truefgl)
        plot(g1,g2,
            title=["(a)" "(b)"],
            titleloc = :right, titlefont = font(12))
    else
        truerecom = calnumrecom(truefgl.offspringgeno)
        g3=scatter(truerecom,estrecom, smooth = true,
            xlabel="True #recom in all chrs",
            ylabel="Estimated #recom",legend=false)
        plot!(g3,x->x)
        l = @layout [a  [b
                         c]]
        plot(g1,g2,g3,layout=l,
            title=["(a)" "(b)" "(c)"],
            titleloc = :right, titlefont = font(12))
    end
end

function saveprobplot(
    magicancestry::MagicAncestry;    
    nplot_subpop::Integer=10,
    probtype::AbstractString="haploprob",
    colorgradient::ColorGradient = cgrad([:white, :blue, :red]),
    boundaryline = (1.5, :black),
    truemarker = (:star, 5, 0.5, :gray, stroke(:gray)),
    truefgl::Union{Nothing,MagicGeno} = nothing,    
    workdir::AbstractString = pwd(),
    outstem::Union{Nothing,AbstractString}="outstem")    
    figdir = workdir
    if nplot_subpop <= 0 
        @warn string("nplot_subpop = ",nplot_subpop, ", <=0, no probplots produced")
        return nothing
    end
    subpop2off = MagicBase.get_subpop2offspring(magicancestry.magicped; isindex=true)    
    offls = sort(reduce(vcat, [length(i)<nplot_subpop ? i : i[1:nplot_subpop] for i in values(subpop2off)]))
    offinfo = magicancestry.magicped.offspringinfo        
    if in("isoutlier",names(offinfo))
        iiout = findall([ismissing(i) ? false : i for i in offinfo[!, :isoutlier]])
    else
        iiout = []
    end
    iinonout = setdiff(offls, iiout)
    outstem2 = isnothing(outstem) ? "" : outstem
    for k = 1:2
        ii = k == 1 ? iiout : iinonout
        outid = k == 1 ? probtype*"_plot_outlier_" : probtype*"_plot_"
        outid = outstem2*"_"*outid
        isempty(ii) && continue
        for i in ii            
            # offid = magicancestry.magicped.offspringinfo[i, :individual]
            # offid2 = replace(offid,"/"=>"_","\\"=>"_")
            fn = joinpath(
                figdir,
                string(outid, i, "th_offspring.png"),
            )
            plotcondprob(magicancestry; offspring = i,
                probtype,colorgradient,boundaryline,truemarker,truefgl,
                outfile = fn
            )
        end
    end
    figdir
end


function loadtarplot(plottarfile::AbstractString; 
    subfiles::Union{Nothing,AbstractVector}=nothing,
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString=tempdir())
    isnothing(subfiles) || eltype(subfiles) <: Integer || error(string("subfiles = ",subfiles, " is not a list of indices"))
    jltempdir = mktempdir(tempdirectory;
        prefix="jl_readprobplot_", cleanup=true)
    try 
        Tar.extract(getabsfile(workdir,plottarfile),jltempdir)
        filels = readdir(jltempdir;sort=true, join=true)
        if !isnothing(subfiles)
            nf = length(filels)
            ls = filter(x-> 1<=x<=nf, subfiles)
            filels = filels[ls]
        end 
        [load(file) for file in filels]
    finally
        rm.(readdir(jltempdir;join=true);force=true)
        rm(jltempdir; force=true,recursive=true)
    end
end



function plotviterbipath(magicancestry::MagicAncestry;
    offspring::Union{Nothing,Integer}=nothing,    
    boundaryline = (1.5, :gray),    
    viterbimarker=(:circle,10,0.5, :cyan, stroke(:gray)),    
    truemarker = (:xcross, 8, 0.5, :gray30, stroke(:gray30)),    
    truefgl::Union{Nothing,MagicGeno}=nothing,
    outfile::Union{Nothing, AbstractString}=nothing,
    plotkeyargs...)        
    vpath = magicancestry.viterbipath
    isnothing(vpath) && error("viterbi path is nothing")
    fglls = magicancestry.statespace["haplotype"]
    nfgl = length(fglls)
    diplostates = MagicBase.prior_diploindex(nfgl)
    ylabel = "Founder"
    titlestem = "Viterbi path"    
    noff = size(vpath[1],2)
    if isnothing(offspring)        
        off = rand(1:noff)
    else
        offspring > noff && error(string("offspring index=",offspring, " > #offspring=",noff))
        off=offspring        
    end    
    offinfo = magicancestry.magicped.offspringinfo    
    isoutlier = ismissing(offinfo[off, :isoutlier]) ? false : offinfo[off, :isoutlier]
    offid = string(offinfo[off,:individual], isoutlier ? ", too many #breakpoints" : "")
    offpathls = [i[:,off] for i in vpath]
    if !isnothing(truefgl)
        # aligntruegeno!(truefgl,magicancestry)
        truechrls = MagicBase.findtruelg(truefgl.markermap,magicancestry.markermap)        
        offtrue = []
        for chr in eachindex(truechrls)
            true_chr = truechrls[chr]
            estsnpidls = magicancestry.markermap[chr][!,:marker]
            true_ii = [begin
                truesnpidls = truefgl.markermap[c][!,:marker]
                if estsnpidls == truesnpidls
                    1:length(estsnpidls)
                else
                    d = setdiff(estsnpidls,truesnpidls)
                    if !isempty(d)
                        estchrid =  magicancestry.markermap[chr][1,:linkagegroup]
                        msg= string("markers of chr=",estchrid, " in magicancestry but not in truefgl: ",d)
                        @warn msg
                    end
                    rule = Dict(truesnpidls .=> 1:length(truesnpidls))
                    ii = [get(rule, i, missing) for i in estsnpidls]
                    collect(skipmissing(ii))
                end
            end for c in true_chr]
            off_fgl = reduce(vcat, [truefgl.offspringgeno[true_chr[c]][true_ii[c],off] for c in 1:length(true_chr)])
            # for each chromosome, adjust unidentifiable the absoluate phase
            offest = diplostates[offpathls[chr]]
            nerr0 = sum(off_fgl .!= offest)
            nerr1 = sum(reverse.(off_fgl) .!= offest)
            if nerr0 > nerr1
                off_fgl .= reverse.(off_fgl)
            end
            append!(offtrue,off_fgl)       
        end                
        dict=Dict(diplostates .=> 1:length(diplostates))            
        ytruels =  [get(dict, i, nothing) for i=offtrue]                                        
    end   
    # calculate xposls
    poscmls = [i[!,:poscm] for i in magicancestry.markermap]
    chrlenls = accumulate(+,last.(poscmls))
    pushfirst!(chrlenls,0.0)
    pop!(chrlenls)    
    xposls = map((x,y)->x .+ y, poscmls, chrlenls)    
    xls = reduce(vcat,xposls)
    yestls = reduce(vcat, offpathls)
    yset = sort(unique(yestls))
    if !isnothing(truefgl)
        yset = sort(unique(vcat(yset, unique(ytruels))))
    end
    ytick0 = ["("*join(fglls[i],",")*")" for i in diplostates[yset]]
    ytick = (collect(1:length(ytick0)),ytick0)    
    ysetdict = Dict(yset .=> 1:length(yset))
    yestls .= [ysetdict[i] for i in yestls]    
    vpathplot = scatter(xls,yestls; 
        marker = viterbimarker, 
        legend=false,
        size = (1500,max(500,length(yset)*40)), 
    ) 
    if !isnothing(truefgl)
        ytruels .= [ysetdict[i] for i in ytruels]
        scatter!(vpathplot,xls,ytruels,marker=truemarker,legend=false)               
    end
    scatter!(xls, 0.08*ones(length(xls)), marker=(:vline,8,:black))    
    x0 = last.(xposls)
    pushfirst!(x0, xposls[1][1])
    nstate = length(yset)
    xy=hcat([[[i,0] [i,nstate+0.75] [NaN,NaN]] for i=x0]...)
    chrlabx=xy[1,1:3:end]
    chrlabx=(chrlabx[1:end-1] .+ chrlabx[2:end]) ./ 2
    ypos = nstate+0.5
    # ypos = nstate + 1.0
    chridls= [i[1,:linkagegroup] for i=magicancestry.markermap]
    annotate!(vpathplot,[(chrlabx[i],ypos,Plots.text(chridls[i],12,:black,rotation=15))
        for i=1:length(chrlabx)])
    plot!(vpathplot,xy[1,:],xy[2,:];
        line = boundaryline,
        legend=false,
        yticks = ytick,
        xlabel="Position (cM)",
        ylabel=ylabel,
        title=string(titlestem, " for ",off,"-th offspring =",offid),
        titlefont = font("Times", 12),
        plotkeyargs...
    )
    isnothing(outfile) || savefig(vpathplot,outfile)
    vpathplot    
end

function saveviterbiplot(
    magicancestry::MagicAncestry;    
    nplot_subpop::Integer=10,
    boundaryline = (1.5, :black),
    viterbimarker=(:circle,10,0.5, :cyan, stroke(:gray)),    
    truemarker = (:xcross, 8, 0.5, :gray30, stroke(:gray30)),
    truefgl::Union{Nothing,MagicGeno} = nothing,    
    workdir::AbstractString = pwd(),
    outstem::Union{Nothing,AbstractString}="outstem")    
    figdir = workdir
    if nplot_subpop <= 0 
        @warn string("nplot_subpop = ",nplot_subpop, ", <=0, no probplots produced")
        return nothing
    end
    subpop2off = MagicBase.get_subpop2offspring(magicancestry.magicped; isindex=true)    
    offls = sort(reduce(vcat, [length(i)<nplot_subpop ? i : i[1:nplot_subpop] for i in values(subpop2off)]))
    offinfo = magicancestry.magicped.offspringinfo    
    iiout = findall([ismissing(i) ? false : i for i in offinfo[!, :isoutlier]])
    iinonout = setdiff(offls, iiout)
    outstem2 = isnothing(outstem) ? "" : outstem
    for k = 1:2
        ii = k == 1 ? iiout : iinonout
        outid = k == 1 ? "_plot_outlier_" : "_plot_"
        outid = outstem2*outid
        isempty(ii) && continue
        for i in ii            
            # offid = magicancestry.magicped.offspringinfo[i, :individual]
            # offid2 = replace(offid,"/"=>"_","\\"=>"_")
            fn = joinpath(
                figdir,
                string(outid, i, "th_offspring.png"),
            )
            plotviterbipath(magicancestry; offspring = i,
                boundaryline,viterbimarker, truemarker,truefgl,
                outfile = fn
            )
        end
    end
    figdir
end
