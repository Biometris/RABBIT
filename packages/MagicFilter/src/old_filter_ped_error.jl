
function filter_ped_error!(magicgeno::MagicGeno;
    model::AbstractString = "jointmodel",
    seqerror::Real = 0.001,
    isfounderinbred::Bool=true,
    mapavailable::Bool=false,
    thresh_epsf_nprogeny::Integer = 10,
    max_epsf::Real = 0.2,
    max_epso::Real = 0.4,
    accuracy_eps::Real = 0.01,
    maxiteration::Integer = 20,
    outstem::Union{Nothing,AbstractString}= "outstem",
    logfile::Union{Nothing,AbstractString,IO} = isnothing(outstem) ? nothing : outstem*"_purify_ped.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "filter_ped_error!"; verbose,delim="-")
    msg = string("list of options: \n",
        "model = ", model, "\n",
        "seqerror = ", seqerror, "\n",
        "isfounderinbred = ", isfounderinbred, "\n",
        "mapavailable = ", mapavailable, "\n",
        "thresh_epsf_nprogeny = ", thresh_epsf_nprogeny, "\n",
        "max_epsf = ", max_epsf, "\n",
        "max_epso = ", max_epso, "\n",
        "accuracy_eps = ", accuracy_eps, "\n",
        "maxiteration = ", maxiteration, "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)
    MagicBase.setunphasedgeno!(magicgeno)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
	missgeno =  MagicBase.split_missprogeny!(magicgeno)
	nsnpmiss = sum(size.(missgeno.markermap,1))
	if nsnpmiss > 0
		msg = string("remove (and insert afterwards) ", nsnpmiss, " markers with all offspring genotypes being missing")
		printconsole(logio,verbose,msg)
	end
    historyfile = (isnothing(outstem) ? "" : outstem)*"_ind_allelicerror_history.txt"
    historyfile = joinpath(workdir,historyfile)
    epsf_vec,epso_vec = infer_eps(magicgeno; model,snpthin=1,
        isfounderinbred,mapavailable, historyfile,
        accuracy_eps,maxiteration,logio,verbose)
    isnothing(outstem) && rm(historyfile;force=true)
    epsf_df,epso_df = get_eps_df(magicgeno, epsf_vec,epso_vec,
        thresh_epsf_nprogeny,max_epsf,max_epso)
	if nsnpmiss > 0
		MagicBase.merge_missprogeny!(magicgeno,missgeno)
	end
	filter_ped_error!(magicgeno,epsf_df, epso_df;logio,verbose)
    if !isnothing(outstem)
        # save allelic error result
        outfile = outstem*"_ind_allelicerror.csv"
        open(joinpath(workdir,outfile),"w") do io
            write(io, "MagicFilter,founder\n")
            CSV.write(io, epsf_df; header =true,append=true)
            write(io, "MagicFilter,offspring\n")
            CSV.write(io, epso_df; header =true,append=true)
        end
        msg = string("save allelicerror file: ", outfile)
        MagicBase.printconsole(logio,verbose,msg)
        # save allelic error plot
        try 
            figeps = plot_eps(epsf_vec,epso_vec,max_epsf,max_epso)
            # display(figeps)
            outfile = outstem*"_ind_allelicerror.png"
            savefig(figeps, joinpath(workdir,outfile))
            msg = string("save allelicerror plot: ", outfile)
            MagicBase.printconsole(logio,verbose,msg)  
        catch
            @warn string("Could not plot allelic error rates")
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"filter_ped_error!"; verbose,delim="-")
    magicgeno
end

function filter_ped_error!(magicgeno::MagicGeno, epsf_df::DataFrame,epso_df::DataFrame;
    logio::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    nchr = length(magicgeno.markermap)
    # delete offspring
    b_del = epso_df[!,:isdel]
    oldmemls = unique(magicgeno.magicped.offspringinfo[!,:member])
    if any(b_del)
        bad = magicgeno.magicped.offspringinfo[b_del,:individual]
        deleteat!(magicgeno.magicped.offspringinfo,b_del)
        for chr in 1:nchr
            magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][:,.!b_del]
        end
        if !isnothing(logio)
            write(logio, "delete offspring: \n")
            CSV.write(logio,epso_df[b_del,:]; header =true,append=true)
        end
        verbose && @info string("delete offspring: ", epso_df[b_del,:])
    else
        MagicBase.printconsole(logio, verbose, "no bad offspring detected")
    end
    # set founders' genotypes to missing
    b_miss = epsf_df[!,:set2missing]
    if any(b_miss)
        for chr in 1:nchr
            magicgeno.foundergeno[chr][:,b_miss] .= "N"
        end
        bad = magicgeno.magicped.founderinfo[b_miss,:individual]
        if !isnothing(logio)
            write(logio, "remove founders' genotypes: \n")
            CSV.write(logio,epsf_df[b_miss,:]; header =true,append=true)
        end
        verbose && @info string("remove founders' genotypes: ", epsf_df[b_miss,:])
    end
    # delete founders and set new designinfo
    b_del = epsf_df[!,:isdel]
    deletions = epsf_df[b_del,:founder]
    newmemls = unique(magicgeno.magicped.offspringinfo[!,:member])
    if issetequal(newmemls, oldmemls)
        isempty(deletions) || @error string("undeleted founders: ",deletions)
    else
        oldfounders = Pedigrees.getfounderid(magicgeno.magicped.designinfo)
        newped = Pedigrees.getsubped(magicgeno.magicped.designinfo, newmemls)
        newfounders = Pedigrees.getfounderid(newped)
        if isempty(newfounders)
            @error "all founder deleted!"
        else
            if issetequal(oldfounders, newfounders)
                isempty(deletions)|| @error string("untracked founder deletion: ",deletions)
            else
                magicgeno.magicped.designinfo = newped
                del_founders = setdiff(oldfounders, newfounders)
                d = setdiff(deletions,del_founders)
                isempty(d) || @error string("untracked founder deletion: ",d)
                d = setdiff(del_founders, deletions)
                if !isempty(d)
                    b = [in(i,d) for i in epsf_df[!,:founder]]
                    epsf_df[b, :remark] .= "no progeny after other deletions"
                    epsf_df[b, :isdel] .= true
                end
                fdict = Dict(oldfounders .=> 1:length(oldfounders))
                findices = [fdict[i] for i in newfounders]
                magicgeno.magicped.founderinfo = magicgeno.magicped.founderinfo[findices,:]
                for chr in 1:nchr
                    magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][:,findices]
                end
                b_del = epsf_df[!,:isdel]
                if !isnothing(logio)
                    write(logio, "delete founders: \n")
                    CSV.write(logio,epsf_df[b_del,:]; header =true,append=true)
                end
                verbose && @info string("delete founders: ", epsf_df[b_del,:])
            end
        end
    end
    b_del = epsf_df[!,:isdel]
    if !(any(b_miss) || any(b_del))
        MagicBase.printconsole(logio, verbose, "no bad founder detected")
    end
    magicgeno, epsf_df
end

function get_eps_df(magicgeno::MagicGeno,
    epsf_vec::AbstractVector, epso_vec::AbstractVector,
    thresh_epsf_nprogeny::Integer,
    max_epsf::Real,max_epso::Real)
    # founder
    founderls = magicgeno.magicped.founderinfo[!,:individual]
    founder_progeny = MagicBase.get_founder2offspring(magicgeno.magicped)
    nprogeny = [length(founder_progeny[i]) for i in founderls]
    bad_founders = findall(epsf_vec .> max_epsf)
    del_founders = bad_founders[nprogeny[bad_founders] .< thresh_epsf_nprogeny]
    miss_founders = bad_founders[nprogeny[bad_founders] .>= thresh_epsf_nprogeny]
    isdel = falses(length(founderls))
    isdel[del_founders] .= true
    set2missing = falses(length(founderls))
    set2missing[miss_founders] .= true
    remark = repeat(["keep"], length(founderls))
    remark[miss_founders] .= "set genotypes to missing"
    remark[del_founders] .= "del founder and progeny"
    epsf_df = DataFrame(founder=founderls, nprogeny=nprogeny,allelic_error=epsf_vec,
        isdel=isdel,set2missing=set2missing,remark=remark)
    # offspring
    offspringls = magicgeno.magicped.offspringinfo[!,:individual]
    subpop = magicgeno.magicped.offspringinfo[!,:member]
    isdel = epso_vec .> max_epso
    remark =  ["keep" for _ in 1:length(offspringls)]
    remark[isdel] .= string("due to epso>",max_epso)
    f2off = MagicBase.get_founder2offspring(magicgeno.magicped)
    for f in founderls[del_founders]
        off_todel = f2off[f]
        b = [in(i,off_todel) for i in offspringls]
        remark[b] .= [i=="keep" ? string("due to del of ",f) : i*string(" and deletion of ",f)
            for i in remark[b]]
        isdel[b] .= true
    end
    epso_df = DataFrame(offspring=offspringls, subpop=subpop,
        allelic_error=epso_vec,isdel=isdel,remark=remark)
    epsf_df, epso_df
end

# function purifygeno!(magicgeno::MagicGeno,
#     epsf_vec::AbstractVector, epso_vec::AbstractVector,
#     max_epsf::Real,max_epso::Real,
#     logio::Union{Nothing,IO},verbose::Bool)
#     nchr = length(magicgeno.markermap)
#     b = epsf_vec .> max_epsf
#     if any(b)
#         bad = magicgeno.magicped.founderinfo[b,:individual]
#         msg = string("bad founders = ",bad,
#             " with epsf = ",round.(epsf_vec[b],digits=2))
#         msg *= "; set their genotypes to missing"
#         MagicBase.printconsole(logio, verbose, msg)
#         # magicgeno.magicped.founderinfo[b,:individual] .= string.(bad,"_replaced")
#         # if isa(magicgeno.magicped.designinfo,Pedigree)
#         #     # Pedigree is a non-mutable struct
#         #     newfid = magicgeno.magicped.founderinfo[!,:individual]
#         #     magicgeno.magicped.designinfo = setfounderid(magicgeno.magicped.designinfo,newfid)
#         # end
#         for chr in 1:nchr
#             magicgeno.foundergeno[chr][:,b] .= "N"
#         end
#     else
#         MagicBase.printconsole(logio, verbose, "no bad founders detected")
#     end
#     b = epso_vec .> max_epso
#     if any(b)
#         bad = magicgeno.magicped.offspringinfo[b,:individual]
#         msg = string("delete bad offspring = ",bad,
#             " with epsf = ",round.(epso_vec[b],digits=2))
#         MagicBase.printconsole(logio, verbose, msg)
#         deleteat!(magicgeno.magicped.offspringinfo,b)
#         for chr in 1:nchr
#             magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][:,.!b]
#         end
#     else
#         MagicBase.printconsole(logio, verbose, "no bad offspring detected")
#     end
#     magicgeno
# end

function plot_eps(epsf_vec::AbstractVector,epso_vec::AbstractVector,
    max_epsf::Real,max_epso::Real;
    plotmarker=(:star, 3, 0.5,:blue,stroke(:blue)),
    line_threshold = (1.5,:dot,:red))
    figf = plot([0,length(epsf_vec)],[max_epsf,max_epsf];
        line=line_threshold,
        size = (1000,600), 
    )
    plot!(figf,epsf_vec;
        marker=plotmarker,legend=false,color=:gray,
        xlabel = "founder index",
        ylabel = "allelic error (ϵf)",        
    )
    figo = plot([0,length(epso_vec)],[max_epso,max_epso];
        line=line_threshold,
        size = (900,600), 
    )
    plot!(figo,epso_vec;
        marker=plotmarker,legend=false,color=:gray,
        left_margin=20Plots.mm,
        bottom_margin=15Plots.mm,        
        xlabel = "offspring index",
        ylabel = "allelic error  (ϵo)"
    )
    plot(figf,figo,layout=(1,2))
end
