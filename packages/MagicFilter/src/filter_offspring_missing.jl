
function filter_offspring_missing!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    isfounderinbred::Bool=true,
    threshcall::Real = model == "depmodel" ? 0.95 : 0.9,
    offspring_maxmiss::Real = 0.99,
    outstem::AbstractString= "outstem",
    logfile::Union{AbstractString,IO} = outstem*"_purify_missing.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    offspring_maxmiss > 1.0 && return magicgeno
    logio = MagicBase.set_logfile_begin(logfile, workdir, "filter_offspring_missing!"; verbose)
    model = MagicBase.reset_model(magicgeno.magicped,model;io=logio,verbose)
    msg = string("list of options: \n",
        "model = ", model, "\n",
        "isfounderinbred = ", isfounderinbred, "\n",
        "threshcall =", threshcall, "\n",
        "offspring_maxmiss = ", offspring_maxmiss, "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)
    MagicBase.setunphasedgeno!(magicgeno)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
    MagicBase.info_missing(magicgeno;io=logio,verbose)
    # founder missing
    formatls = reduce(vcat,magicgeno.markermap)[!,:founderformat]
    geno = reduce(vcat,magicgeno.foundergeno)
    missfrac = mean(MagicBase.getismissing(geno, formatls),dims=1)[1,:]
    indls = magicgeno.magicped.founderinfo[!,:individual]
    res = DataFrame(individual=indls,missfrac=missfrac,isfounder=true,keep=true)
    # offspring missing
    formatls = reduce(vcat,magicgeno.markermap)[!,:offspringformat]
    if unique(formatls) == ["GT_unphased"]
        calledgeno = reduce(vcat,magicgeno.offspringgeno)
    else
        calledgeno = reduce(vcat, [call_offgeno(magicgeno, model,chr;callthreshold = threshcall) for chr in 1:length(magicgeno.markermap)])
        formatls = ["GT_unphased" for _ in 1:length(formatls)]

    end
    missfrac = mean(MagicBase.getismissing(calledgeno, formatls),dims=1)[1,:]    
    indls = magicgeno.magicped.offspringinfo[!,:individual]
    offkeep = missfrac .<= offspring_maxmiss
    any(offkeep) || @error "All offspring are deleted!"
    res = vcat(res, DataFrame(individual=indls,missfrac=missfrac,isfounder=false,keep=offkeep))
    outfile = outstem*"_ind_missing.csv"
    outfile2 = joinpath(workdir,outfile)
    open(outfile2,"w") do outio
        commentstring = "##"
        msg = commentstring*"col_1, individual, individual (founder or offsring) ID \n"
        msg *= commentstring*"col_2, missfrac, fraction of missing genotypes in the individual\n"
        msg *= commentstring*"col_3, isfounder, if the individual is a founder\n"        
        msg *= commentstring*"col_4, keep, the individual is removed if keep = false\n"
        write(outio, msg)
        CSV.write(outio,res; header=true, append=true)
    end
    msg = string("save frac_missing_genotype for each individual in ", outfile)
    printconsole(logio, verbose,msg)
    try 
        plot_miss_ind(outfile2;annotate_offspringmiss = offspring_maxmiss)
        savefig(joinpath(workdir,outstem*"_ind_missing.png"))   
    catch
        @warn string("Could not plot missing fractions for each offspring")
    end
    del_offspring!(magicgeno,findall(offkeep);io=logio,verbose)
    MagicBase.info_missing(magicgeno;io=logio,verbose)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"purify_missing"; verbose)
    magicgeno
end

function del_offspring!(magicgeno::MagicGeno,keep_offindices::AbstractVector;
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    offls = magicgeno.magicped.offspringinfo[!,:individual]
    deloffls = setdiff(offls,offls[keep_offindices])
    if isempty(deloffls)
        MagicBase.printconsole(io,verbose,"no offspring deleted")
        return magicgeno
    else
        msg = string("delete ", length(deloffls)," offspring: ",deloffls)
        MagicBase.printconsole(io,verbose,msg)
    end
    for chr in 1:length(magicgeno.offspringgeno)
        magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][:,keep_offindices]
    end
    magicgeno.magicped.offspringinfo = magicgeno.magicped.offspringinfo[keep_offindices,:]
    if isa(magicgeno.magicped.designinfo,Pedigrees.Pedigree)
        # up ped
        newmemls = unique(magicgeno.magicped.offspringinfo[!,:member])
        newped = Pedigrees.getsubped(magicgeno.magicped.designinfo, newmemls)
        magicgeno.magicped.designinfo = newped
        # up founders
        oldfounders = magicgeno.magicped.founderinfo[!,:individual]
        newfounders = Pedigrees.getfounderid(newped)
        if newfounders != oldfounders
            d = setdiff(oldfounders,newfounders)
            if !isempty(d)
                msg = string("delete ", length(d)," founders: ",d)
                MagicBase.printconsole(io,verbose,msg)
            end
            dict = Dict(oldfounders .=> 1:length(oldfounders))
            fkeep = [dict[i] for i in newfounders]
            magicgeno.magicped.founderinfo = magicgeno.magicped.founderinfo[fkeep,:]
            for chr in 1:length(magicgeno.foundergeno)
                magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][:,fkeep]
            end
        end
    end
    magicgeno
end

function plot_miss_ind(missfracfile::AbstractString;
    commentstring::AbstractString = "##", 
    annotate_offspringmiss::Union{Nothing,Real}=nothing)
    missdf = CSV.read(missfracfile, DataFrame; comment=commentstring)
    isfounder = missdf[!,:isfounder]
    fmiss = missdf[isfounder,:missfrac]
    fcdf = ecdf(fmiss)
    g1 = plot(x -> fcdf(x), 0, 1;
        xlabel = "Fraction of missing genotypes per founder",
        ylabel = "Cumulative fraction of founders",
        legend=false,    
        title = "Empirical cumulative distribution for founders"    
    )
    isoff = .!missdf[!,:isfounder]
    offmiss = missdf[isoff,:missfrac]
    offcdf = ecdf(offmiss)
    g2 = plot(x -> offcdf(x), 0, 1;
        xlabel = "Fraction of missing genotypes per offspring",
        ylabel = "Cumulative fraction of offspring",
        legend=false,              
        title = "Empirical cumulative distribution for offspring"            
    )
    if !isnothing(annotate_offspringmiss)
        x = round(annotate_offspringmiss,digits=5)
        y = round(offcdf(x),digits=5)
        scatter!(g2,[x],[y],
            series_annotations = text((x,y), :bottom)
        )
    end
    plot(g1,g2;    
        layout=(2,1),
        linewidth = 2, 
        left_margin=20Plots.mm,
        right_margin=10Plots.mm,
        bottom_margin=15Plots.mm,
        size = (1000,800),                
    )
end
