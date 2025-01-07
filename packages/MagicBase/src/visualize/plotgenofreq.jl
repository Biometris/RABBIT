
function plotgenofreq(magicgeno::MagicGeno;
    boundaryline = (1.5, :gray,:dash),
    plotsize=(1500,1000),
    minmaf_subpop::Real=0.01,
    commentstring::AbstractString = "##", 
    outstem::AbstractString="outstem",    
    workdir::AbstractString=pwd())
    # calculate chridls, xposls and xls    
    poscmls = [i[!,:poscm] for i in magicgeno.markermap]
    chrlenls = accumulate(+,last.(poscmls))
    pushfirst!(chrlenls,0.0)
    pop!(chrlenls)        
    xposls = map((x,y)->x .+ y, poscmls, chrlenls)
    xls = reduce(vcat,xposls)
    # calculate yls
    genofreqfile = calgenofreq(magicgeno; minmaf_subpop,outstem,workdir)
    genofreq = CSV.read(getabsfile(workdir,genofreqfile),DataFrame; comment=commentstring)
    all(genofreq[!,:poscm] .≈ reduce(vcat,poscmls)) || @error "inconsistent poscm"
    heterofreq = genofreq[!,:heterofreq]
    freqA2 = genofreq[!,:freqA2]
    ylabells = ["freq of hetero_genotype","freq of ref_allele"]
    figls = [begin 
        yls = g ==1 ? heterofreq : freqA2
        fig = plot_markerpos(magicgeno.markermap; boundaryline)
        plot(fig, xls,yls;
            yrange = (0,1.1),        
            xlabel="position (cM)",    
            ylabel= ylabells[g],            
        )
        scatter!(xls,yls;
            marker = (2,:circle)
        )   
    end for g in 1:2]    
    gmarker_hetero= histogram(heterofreq;
        legend = false, 
        xrange = (0,1),        
        xlabel= ylabells[1],
        ylabel = "#markers"    
    )
    gmarker_A2= histogram(freqA2; 
        legend = false, 
        xrange = (0,1),        
        xlabel= ylabells[2],
        ylabel = "#markers"    
    )
    l = @layout [a{0.3h}
                 b{0.3h}
                grid(1,2)]
    plot(figls...,gmarker_hetero,gmarker_A2;     
        layout = l,
        size = plotsize,
        left_margin=20Plots.mm,
        bottom_margin=15Plots.mm,
    )
    figfile = getabsfile(workdir,outstem*"_genofreq.png")
    savefig(figfile)
    genofreqfile, figfile
end

function plot_markerpos(markermap::AbstractVector; boundaryline, yminmax=(0,1.05))          
    ymin,ymax = yminmax
    chridls= [i[1,:linkagegroup] for i=markermap]
    poscmls = [i[!,:poscm] for i in markermap]
    chrlenls = accumulate(+,last.(poscmls))
    pushfirst!(chrlenls,0.0)
    pop!(chrlenls)        
    xposls = map((x,y)->x .+ y, poscmls, chrlenls)
    x0 = last.(xposls)
    pushfirst!(x0, xposls[1][1])
    chrlabx = [(i[1]+i[end])/2 for i in xposls]
    chrlaby = 1.05        
    xy=hcat([[[i,ymin] [i,ymax] [NaN,NaN]] for i=x0]...)
    chrplot = plot(xy[1,:],xy[2,:];
         line = boundaryline,
         legend = false)
    annotate!(chrplot,[(chrlabx[i],chrlaby,Plots.text(chridls[i],12,:black,rotation=15)) for i= eachindex(chrlabx)])
    chrplot
end

function get_ppoffls(magicped::MagicPed)
    pop2off = MagicBase.get_subpop2offspring(magicped; isindex = true)
    pop2pp = MagicBase.get_subpop2founder(magicped; isindex = true)
    popls = unique(magicped.offspringinfo[!,:member])
    Dict([i=>[pop2pp[i],pop2off[i]] for i in popls])
end

function calgenofreq(magicgeno::MagicGeno;
    minmaf_subpop::Real=0.01,
    outstem::AbstractString="outstem",    
    workdir::AbstractString=pwd())    
    outfile = getabsfile(workdir, outstem*"_genofreq.csv")
    open(outfile,"w") do outio
        subpopls = unique(magicgeno.magicped.offspringinfo[!,:member])
        descripls = ["chromosome, chromosome ID", "marker, marker ID", 
            "poscm, marker position in centi-Morgan",
            "physposbp, marker position in base pair",
        ]
        for subpop in subpopls
            msg = string(subpop, "_n11|n12|n22", ", counts for genotypes 11|12|22 where 1 = referenc eall && 2 = alternative allele for subpopulation = ",subpop)
            msg *= "&& \"=>0|0|0\" denotes that the subpopulation is monomorphic and it is excluded in caluclating heterofreq and freqA2 "
            push!(descripls,msg)
        end
        push!(descripls, "heterofreq, frequency of heterozygous genotypes")
        push!(descripls, "freqA2, frequency of 2nd allele (alternative allele)")
        for i in eachindex(descripls)
            write(outio, string("##col_",i, ", ", descripls[i], "\n"))
        end
        headstr = string.(subpopls,"_n11|n12|n22") # ignore genotypes 1N,2N
        headstr = string("chromosome,marker,poscm,physposbp,",join(headstr,","),",heterofreq",",freqA2")
        write(outio, headstr,"\n")
        ppoffls = get_ppoffls(magicgeno.magicped)
        for chr in 1:length(magicgeno.markermap)
            fgeno = magicgeno.foundergeno[chr]
            offgeno = magicgeno.offspringgeno[chr]
            nsnp = size(fgeno,1)
            for snp in 1:nsnp                
                nls = [zeros(Int,3) for _ in 1:length(subpopls)] 
                a2freqls = zeros(length(subpopls))
                for k in eachindex(subpopls)                    
                    ppls, offls = ppoffls[subpopls[k]]
                    alleleset = unique(split(join(reduce(vcat, unique(view(fgeno,snp,ppls)))),""))
                    isppmono = alleleset in [["1"],["2"]]
                    isppmono && continue
                    snpsubgeno = join.(view(offgeno,snp,offls))
                    alleleset = unique(split(join(unique(snpsubgeno)),""))
                    setdiff!(alleleset,"N")
                    isoffmono = alleleset in [["1"],["2"],[]]
                    isoffmono && continue
                    n11 = sum(snpsubgeno .== "11")
                    n12 = sum(snpsubgeno .== "12") + sum(snpsubgeno .== "21")
                    n22 = sum(snpsubgeno .== "22")
                    ntot = n11+n12+n22
                    ntot == 0 && continue
                    a2freqls[k] = (0.5*n12 + n22)/ntot
                    nls[k] .= [n11,n12,n22]
                end
                b = [minmaf_subpop <= i <= 1.0-minmaf_subpop for i in a2freqls]
                if any(b)
                    n11,n12,n22 = sum(nls[b])
                    ntot = n11+n12+n22
                    heterofreq = round(n12/ntot,digits=4)
                    a2freq = round((0.5*n12 + n22)/ntot,digits=4)
                else
                    heterofreq = NaN
                    a2freq = NaN
                end
                nls_str = map((x,y)->y ? join(x,"|") : (sum(x)>0 ? string(join(x,"|"),"=>0|0|0") : "0|0|0"),nls,b)                
                res = join(magicgeno.markermap[chr][snp,[:linkagegroup,:marker,:poscm,:physposbp]], ",")
                res = string(res,",",join(nls_str,","),",",heterofreq,",",a2freq)
                write(outio,res,"\n")
            end            
        end
        outfile
    end
end