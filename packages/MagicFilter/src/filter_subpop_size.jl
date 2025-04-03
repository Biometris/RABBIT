
function filter_subpop_size!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    isfounderinbred::Bool=true,
    minsubpop::Real = 1,
    outstem::AbstractString= "outstem",
    logfile::Union{AbstractString,IO} = outstem*"_purify_subpop.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    minsubpop < 1 && return magicgeno
    logio = MagicBase.set_logfile_begin(logfile, workdir, "filter_subpop_size!"; verbose)
    model = MagicBase.reset_model(magicgeno.magicped,model;io=logio,verbose)
    msg = string("list of options: \n",
        "model = ", model, "\n",
        "isfounderinbred = ", isfounderinbred, "\n",
        "minsubpop = ", minsubpop, "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)
    MagicBase.setunphasedgeno!(magicgeno)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
    # MagicBase.info_missing(magicgeno;io=logio,verbose)
    # del small subpops
    offinfo = magicgeno.magicped.offspringinfo
    memls = offinfo[!,:member]
    memset = unique(memls)
    offlsls = [findall(memls .== i) for i in memset]
    subsizels = length.(offlsls)
    b = subsizels .< minsubpop
    resdf = DataFrame(subpop=memset,size=subsizels,keep = .!b)
    outfile = outstem*"_subpop_size.csv"
    outfile2 = getabsfile(workdir,outfile)
    open(outfile2,"w") do outio
        commentstring = "##"
        msg = commentstring*"col_1, subpop, subpopulation ID\n"
        msg *= commentstring*"col_2, size, size of the subpopulation\n"        
        msg *= commentstring*"col_3, keep, the subpopulation is removed if keep = false"
        write(outio, msg, "\n")
        CSV.write(outio,resdf; header=true,append=true)
    end
    msg = string("save sizes of each subpop in ", outfile)
    printconsole(logio, verbose,msg)
    try 
        plot_subpop_size(outfile2;annotate_subpop_size = minsubpop-1)
        savefig(joinpath(workdir,outstem*"_subpop_size.png"))  
    catch
        @warn string("Could not plot sizes of sub-populations")
    end
    del_subpop = memset[b]
    indls = magicgeno.magicped.offspringinfo[!,:individual]
    del_off =  isempty(del_subpop) ? [] : indls[reduce(vcat, offlsls[b])]
    if !isempty(del_subpop)
        msg = string(length(del_subpop), " subpops of size < ", minsubpop,
            ": ", join(del_subpop,","))
        printconsole(logio, verbose, msg)
    end
    offkeep = findall(.![i in del_off for i in indls])
    del_offspring!(magicgeno,offkeep;io=logio,verbose)
    # MagicBase.info_missing(magicgeno;io=logio,verbose)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"purify_subpop"; verbose)
    magicgeno
end

function plot_subpop_size(popsizefile::AbstractString;
    commentstring::AbstractString = "##", 
    annotate_subpop_size::Union{Nothing,Real}=nothing)
    df = CSV.read(popsizefile, DataFrame; comment=commentstring)
    sizels = df[!,:size]
    sizecdf = ecdf(sizels; weights=sizels)
    xmax  = max(sizels...)
    g = plot(x -> sizecdf(x), 0, xmax;
        xlabel = "subpop size",
        ylabel = "cumulative fraction of offspring",
        legend=false,
        xlims = (0, xmax),
        ylims = (-0.05,1.05),
        size = (1000,600), 
        left_margin=20Plots.mm,
        bottom_margin=15Plots.mm,        
        title = "subpop-size among subpopulations"
    )
    if !isnothing(annotate_subpop_size)
        x = annotate_subpop_size
        y = round(sizecdf(x),digits=5)
        scatter!(g,[x],[y],
            series_annotations = text((x,y), :bottom),
        )
    end
    g
end
