

function plot_posterior_recom(recomfile::AbstractString; 
    priorline = (2,:dot,:blue),
    workdir::AbstractString=pwd(),
    outstem::Union{Nothing,AbstractString}="outstem")
    recomdf = CSV.read(getabsfile(workdir,recomfile),DataFrame)
    g1 = MagicReconstruct.plot_posterior_recom(recomdf; isperchrom=false)
    g2 = MagicReconstruct.plot_posterior_recom(recomdf; isperchrom=true,priorline)
    fig = plot(g1,g2,
        layout=(2,1),
        size = (1000,700), 
        left_margin = 10Plots.mm,
        bottom_margin = 10Plots.mm,
    )
    outfile = string(outstem,"_posterior_recom.png")
    outfile2 =getabsfile(workdir, outfile)
    savefig(fig, outfile2)
    outfile
end

function plot_posterior_recom(recomdf::AbstractDataFrame; priorline = (2,:dot,:blue), isperchrom::Bool=false)
    chrlen = 0.01 .* parse.(Float64, split(string(recomdf[1,:chrlen_cM]),"|"))
    genolen = sum(chrlen)
    if isperchrom        
        isoutlier = recomdf[!,:isoutlier]
        recomdf = copy(recomdf[.!isoutlier, :])        
        colnames = names(recomdf)
        colls = occursin.("posterior_recomnum_", colnames)
        nrecom = Matrix(recomdf[:,colls])                
        nrecom_den = nrecom ./ permutedims(chrlen)
        chridls = replace.(colnames[colls], "posterior_recomnum_"=>"")
        # ischrint = all(occursin.(r"^[1-9][0-9]{0,6}$",chridls))
        # if ischrint
        #     chridls = parse.(Int, chridls)
        # end
        xls = repeat(chridls,inner = size(recomdf,1))
        yls = vec(nrecom_den)
        grecom = violin(xls,yls, side=:left, linewidth=0, label="posterior recomden",xdiscrete_values=chridls)
        # dotplot!(grecom, xls,yls, side=:left, linewidth=0, label="")
        yls = vec(nrecom)
        grecom = violin!(grecom,xls,yls, side=:right, linewidth=0, label="posterior recomnum")
        # dotplot!(xls,yls, side=:right, linewidth=0, label="")
        prior_den = mean(recomdf[!,:prior_recomden])
        plot!(grecom, x->prior_den, 
            line = priorline, 
            label="prior density",            
            # legend = :topright
        )
        post_den = mean(recomdf[!,:posterior_recomnum] ./ genolen)
        plot!(grecom; 
            xlabel = "Linkage group", 
            ylabel = "Recombination breakpoints", 
            title = string( "prior density=", round(prior_den,digits=3),", posterior density = ", round(post_den,digits=3), " per Morgan per offspring"), 
            titlefontsize = 10
        )    
        grecom        
    else        
        isoutlier = Bool.(recomdf[!,:isoutlier])
        isnonoutlier = .!isoutlier        
        xls = string.(round.(recomdf[isnonoutlier,:prior_recomden],digits=3))
        yls = recomdf[isnonoutlier,:posterior_recomden]
        grecom = violin(xls, yls; 
            linewidth=0,legend=false,
            size=(800,500),
            xlabel="prior breakpoint density (per Morgan per chromosome-pair)",
            ylabel="posterior breakpoint density",
            left_margin=20Plots.mm,
            bottom_margin=15Plots.mm,                
            bar_width = 0.05,     
            label = ""    
        )
        dotlabel = any(isoutlier) ? "Non-outlier" : ""
        dotplot!(xls, yls; 
            bar_width = 0.05, 
            marker=(:black,:circle,2.5),    
            legend = true, 
            label = dotlabel
        )
        if any(isoutlier)
            xls = string.(round.(recomdf[isoutlier,:prior_recomden],digits=3))            
            yls = recomdf[isoutlier,:posterior_recomden]
            dotlabel = any(isoutlier) ? "Outlier" : ""
            dotplot!(xls, yls; 
                bar_width = 0.05, 
                marker=(:red,:xcross,4),
                label = dotlabel
            )
        end        
        title!(string("genome length = ", round(genolen,digits=2), " Morgan"))
        grecom
    end
end


function save_posterior_recom(magicped::MagicPed; outstem, workdir)
    recomdf = MagicBase.offspringinfo2df(magicped.offspringinfo)
    outfile = string(outstem, "_posterior_recom.csv")
    outfile2 = MagicBase.getabsfile(workdir,outfile)
    CSV.write(outfile2,recomdf;header=true)
    outfile2
end        

# function save_posterior_recom(magicancestry::MagicAncestry;
#     minprob::Real= 0.7,    
#     isautosome::Bool=true,
#     workdir::AbstractString = pwd(),
#     outstem::AbstractString = "outstem",
#     io::Union{Nothing,IO}=nothing,
#     verbose::Bool=true)
#     telapse=@elapsed begin 
#         offinfo = magicancestry.magicped.offspringinfo
#         recomdf = MagicBase.offspringinfo2df(offinfo)
#         outfile = string(outstem, "_posterior_recom.csv")
#         outfile2 = MagicBase.getabsfile(workdir,outfile)
#         open(outfile2,"w") do outio            
#             descripls = ["individual, offspring ID", 
#                 "member, subpopulation ID for the offpsring", 
#                 "prior_inbredcoef, prior inbreeding coefficient for the offpring",
#                 "prior_recomden, prior number of recombination breakpoints per Morgan for the offspring",
#                 "posterior_recomden, posterior number of recombination breakpoints per Morgan for the offspring",
#                 "posterior_recomnum, posterior number of recombination breakpoints for the offspring",
#                 "length_cM, genomic length in cent-Morgan"
#             ]            
#             for i in eachindex(descripls)
#                 write(outio, string("##col_",i, ", ", descripls[i], "\n"))
#             end
#             CSV.write(outio,recomdf; append=true,header=true)
#         end
#     end
# end

function posterior_recom(magicancestry::MagicAncestry;minprob::Real= 0.7)    
    # post_nrecom is a matrix A; A[i,j] = #recombreakpoints for offsprign i in LG j
    post_nrecom = MagicBase.calnumrecom(magicancestry;isperchrom=true,minprob) 
    colid = [string("posterior_recomnum_", i[1,:linkagegroup]) for i in magicancestry.markermap]
    df = DataFrame(post_nrecom,colid)
    insertcols!(df,1, "individual" => magicancestry.magicped.offspringinfo[!,:individual])
    chrlenls = [i[end,:poscm]-i[1,:poscm] for i in magicancestry.markermap]
    chrlen = join(round.(chrlenls,digits=4),"|")
    insertcols!(df,size(df,2)+1,"chrlen_cM"=>chrlen)
    sum_nrecom = sum(post_nrecom,dims=2)[:,1]
    insertcols!(df,size(df,2)+1,"posterior_recomnum"=>sum_nrecom)    
    genolen = 0.01 * sum(chrlenls)
    recomden = sum_nrecom./ genolen
    insertcols!(df,size(df,2)+1,"posterior_recomden"=>recomden)
    df
end

function prior_recom(magicped::MagicPed;
    isfounderinbred::Bool=true,
    isautosome::Bool=true)
    offls = magicped.offspringinfo[!,:individual]
    memls = magicped.offspringinfo[!,:member]
    ishomozygousls = magicped.offspringinfo[!,:ishomozygous]
    design =  magicped.designinfo    
    resdf = DataFrame(:individual=>offls, :prior_inbredcoef =>0.0, :prior_recomden=>0.0)
    if isa(design, Pedigrees.Pedigree)        
        uni_memls = unique(memls)
        res = MagicPrior.magicorigin(design; 
            memberlist=uni_memls,isfounderinbred,
            isautosome,isfglexch = true, isconcise=true,
        )
        # res[1,:] = ["a", "phi12^mp","R^m","R^p","rho^mp","J1112^mp", ...]
        for i in 2:size(res,1)
            b = memls .== res[i,1]            
            ishomozygous = all(ishomozygousls[b])
            inbredcoef = ishomozygous ? 1.0 : 1.0 - res[i,2] 
            recomden = ishomozygous ? (res[i,3]+res[i,4])/2.0 : res[i,5]
            resdf[b,:prior_inbredcoef] .= inbredcoef            
            resdf[b,:prior_recomden] .= recomden            
        end 
    elseif isa(design, Dict{String,DesignInfo})        
        for (popid,subdesign) in design
            if isnothing(subdesign.pedigree)
                juncdist = subdesign.juncdist            
                isnothing(juncdist) && @error string("TODO for subdesign=",subdesign)                
                (; ibd, mapexpansion, j1122, j1211, j1213, j1222,j1232) = juncdist
                b = isnothing.([mapexpansion, j1122, j1211, j1213, j1222,j1232])
                if all(b) 
                    rows = memls .== popid
                    resdf[rows,:prior_inbredcoef] .= ibd
                    @warn string("could not derive recomden from juncdist")                
                else
                    if !isnothing(j1122)
                        Rm = j1122+2*j1222+j1232
                        Rp = j1122+2*j1211+j1213
                        recomden = Rm + Rp - j1122
                    else
                        recomden = isnothing(ibd) ? 2*mapexpansion : (2-ibd)*mapexpansion
                    end
                    rows = memls .== popid
                    resdf[rows,:prior_inbredcoef] .= ibd
                    resdf[rows,:prior_recomden] .= recomden
                end
            else                
                res = MagicPrior.magicorigin(subdesign.pedigree;                     
                    memberlist=[popid], isfounderinbred,
                    isautosome,isfglexch = true, isconcise=true,
                )
                # res[1,:] = ["a", "phi12^mp","R^m","R^p","rho^mp","J1112^mp", ...]
                for i in 2:size(res,1)
                    b = memls .== res[i,1]            
                    ishomozygous = all(ishomozygousls[b])
                    inbredcoef = ishomozygous ? 1.0 : 1.0 - res[i,2] 
                    recomden = ishomozygous ? (res[i,3]+res[i,4])/2.0 : res[i,5]
                    resdf[b,:prior_inbredcoef] .= inbredcoef            
                    resdf[b,:prior_recomden] .= recomden            
                end 
            end
        end
    else
        @error string("unknown typeof designinfo: ", typeof(design))
    end
    resdf
end

function posterior_prior_recom!(magicancestry::MagicAncestry; minprob::Real= 0.7)
    isfounderinbred = MagicBase.get_isfounderinbred(magicancestry)
    prior_nrecom = MagicReconstruct.prior_recom(magicancestry.magicped;isfounderinbred,isautosome=true)
    post_nrecom = MagicReconstruct.posterior_recom(magicancestry; minprob) 
    offinfo = innerjoin(magicancestry.magicped.offspringinfo,post_nrecom, prior_nrecom, on= :individual)
    magicancestry.magicped.offspringinfo = offinfo
end


function detect_outlier!(magicancestry::MagicAncestry;
    tukeyfence::Real=2,istransform::Bool=false)
    offinfo = magicancestry.magicped.offspringinfo 
    postls = offinfo[!,:posterior_recomnum]
    prils = offinfo[!,:prior_recomden]
    outliers = reduce(vcat,[begin         
        offls = findall(prils.== pri)
        nrecomls = postls[offls]        
        anscombe = istransform ? 2 * (sqrt.(nrecomls .+ 3.0/8)) : nrecomls    
        q1,q3 = quantile(anscombe,[0.25,0.75])
        upbound=q3+tukeyfence*(q3-q1)
        offls[anscombe .> upbound]
    end for pri in unique(prils)])    
    isoutlier = falses(size(offinfo,1))
    isoutlier[outliers] .= true
    if in("isoutlier",names(offinfo))
        offinfo[!,:isoutlier] .= isoutlier    
    else
        insertcols!(offinfo,size(offinfo,2)+1,"isoutlier"=>isoutlier)
    end
    offinfo
end

# function detect_outlier!(magicancestry::MagicAncestry;
#     minprob::Real= 0.7,tukeyfence::Real=2,istransform::Bool=false)
#     nrecom= calnumrecom(magicancestry; minprob)
#     anscombe = istransform ? nrecom : 2 * (sqrt.(nrecom .+ 3.0/8))    
#     q1,q3 = quantile(anscombe,[0.25,0.75])
#     upbound=q3+tukeyfence*(q3-q1)
#     outlier=anscombe .> upbound
#     magicancestry.magicped.offspringinfo[!,:isoutlier] = outlier
# end


function plotinbredcoef(magicancestry::MagicAncestry;
    boundaryline = (1.5, :gray,:dash),
    plotsize=(1200,800))
    # scatter plot inbredcoef
    poscmls = [i[!,:poscm] for i in magicancestry.markermap]
    chrlenls = accumulate(+,last.(poscmls))
    pushfirst!(chrlenls,0.0)
    pop!(chrlenls)    
    xposls = map((x,y)->x .+ y, poscmls, chrlenls)
    xls = reduce(vcat,xposls)
    yls = reduce(vcat,[mean(i,dims=2)[:,1] for i in magicancestry.inbredcoef])
    inbredplot = plot(xls,yls;        
        yrange = (0,1.1),                    
    )
    scatter!(inbredplot,xls,yls;                
        marker = (2,:circle),    
    )
    # annotate chromosome
    x0 = last.(xposls)
    pushfirst!(x0, xposls[1][1])
    chrlabx = [(i[1]+i[end])/2 for i in xposls]
    chrlaby = 1.05
    chridls= [i[1,:linkagegroup] for i=magicancestry.markermap]
    annotate!(inbredplot,[(chrlabx[i],chrlaby,Plots.text(chridls[i],12,:black,rotation=15)) for i= eachindex(chrlabx)])
    xy=reduce(hcat,[[[i,0] [i,1.05] [NaN,NaN]] for i=x0])
    plot!(inbredplot, xy[1,:],xy[2,:];
        line = boundaryline,
        legend = false,    
        xlabel="position (cM)",
        ylabel="inbreeding coefficient",    
        # titlefont = font("Times", 12)
    )
    # add theoretical expectation
    isfounderinbred = MagicBase.get_isfounderinbred(magicancestry)
    df = MagicReconstruct.prior_recom(magicancestry.magicped; isfounderinbred,isautosome = true)
    avginbred = mean(df[!,:prior_inbredcoef])
    plot!(inbredplot, xls,avginbred*ones(length(xls)); line=(2.0,:red))    
    gmarker= histogram(yls; 
        legend = false, 
        xrange = (0,1),        
        xlabel = "inbreeding coefficient",
        ylabel = "#markers"    
    )
    scatter!(gmarker, [avginbred],[0.05]; marker = (10,:circle))
    nsnp = length(yls)
    inbredls = sum([sum(i,dims=1)[1,:] for i in magicancestry.inbredcoef]) ./ nsnp
    goff = histogram(inbredls; 
        legend = false, 
        xrange = (0,1),        
        xlabel = "inbreeding coefficient",
        ylabel = "#offpsring"    
    )
    scatter!(goff, [avginbred],[0.05]; marker = (10,:circle))
    l = @layout [a{0.5h}
                grid(1,2)]
    plot(inbredplot,gmarker,goff;     
        layout = l,
        size=plotsize,
        left_margin=20Plots.mm,
        bottom_margin=15Plots.mm,
    )
end


