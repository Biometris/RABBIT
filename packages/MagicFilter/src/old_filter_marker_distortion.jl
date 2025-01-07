function filter_distortion(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};
    model::AbstractString="depmodel",
    log10_siglevel::Real = -10.0,
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    commentstring::AbstractString="##",
    outstem::AbstractString= "outstem",
    logfile::Union{AbstractString,IO} = outstem*"_mono.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    if !isfounderinbred
        @error string("filter_distortion required inbred founders!")
        return
    end
    magicgeno = formmagicgeno(genofile, pedinfo;
        isfounderinbred,isphysmap=false,
        formatpriority, commentstring, missingstring="NA", workdir)
    filter_distortion!(magicgeno;
        model,log10_siglevel, outstem, logfile, workdir, verbose)
    outfile = outstem*"_distortion.csv"
    savegenodata(outfile,magicgeno; commentstring,workdir)
    outfile
end


# it works only for homogenous population with inbred parents
function filter_distortion!(magicgeno::MagicGeno;
    model::AbstractString="depmodel",
    log10_siglevel::Real = -10.0,
    outstem::AbstractString="outstem",
    logfile::Union{AbstractString,IO} = outstem*"_mono.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    # model = "depmodel" homogenous population with inbred parents
    if model != "depmodel"
        @error string("TODO for model=",model)
    end
    isfounderinbred = true
    founderformat = unique(reduce(vcat,[unique(i[!,:founderformat]) for i in magicgeno.markermap]))
    offformat = unique(reduce(vcat,[unique(i[!,:offspringformat]) for i in magicgeno.markermap]))
    if founderformat != ["GT_haplo"]
        @error "founders genoformat must be GT_haplo"
        return
    end
    if !in(offformat, ["GT_unphased","GT_phased"] )
        @error "offspring genoformat must be GT_unphased or GT_phased"
        return
    end
    if typeof(logfile) <: AbstractString
        logio=open(MagicBase.getabsfile(workdir,logfile), "w")
        # round(Dates.now(),Dates.Second)
        msg = string("filter_distortion! logfile=", logfile)
        MagicBase.printconsole(logio,verbose,msg)
    else
        logio = logfile
    end
    MagicBase.printconsole(logio,verbose,string(repeat("-",30),"filter_distortion!",repeat("-",29)))
    msg = string("list of options: \n",
        "model = ", model, "\n",
        "log10_siglevel = ", log10_siglevel, "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)
    MagicBase.setunphasedgeno!(magicgeno)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
    nchr = length(magicgeno.markermap)
    for chr in 1:nchr
        chrgeno = magicgeno.offspringgeno[chr]
        gset = unique(chrgeno)
        genodiff = setdiff(unique(chrgeno),["11","12","22","NN","1N","2N"])
        if !isempty(genodiff)
            @error string("unexpected genotypes=",genodiff, " on chr=",chr)
        end
        in("12",gset) && (chrgeno[chrgeno .== "12"] .= "NN")
        in("1N",gset) && (chrgeno[chrgeno .== "1N"] .= "11")
        in("2N",gset) && (chrgeno[chrgeno .== "2N"] .= "22")
    end
    colnames = Symbol.(["linkagegroup","marker","subpop","log10p",
        "genocount","chr_index","marker_index"])
    outfiles = [joinpath(workdir,outstem*i) for i in
        ["_distortion_log.csv","_distortion.csv"]]
    for i in outfiles
        open(i, "w") do io
            write(io,join(colnames[1:5],","),"\n")
        end
    end
    magicprior= MagicReconstruct.calmagicprior(magicgeno,model;isfounderinbred)
    for chr in 1:nchr
        distortions = Vector()
        chrid = magicgeno.markermap[chr][1,:linkagegroup]
        foundergeno = magicgeno.foundergeno[chr]
        offgeno = magicgeno.offspringgeno[chr]
        @time popmakeup = MagicReconstruct.calpopmakeup(magicgeno,chr, model, magicprior;isfounderinbred)
        popidls = collect(keys(popmakeup))
        sizels = [length(popmakeup[popid]["offspring"]) for popid in popidls]
        popidls = popidls[sortperm(sizels,rev=true)]
        for popid in popidls
            founders = popmakeup[popid]["founder"]
            offspring = popmakeup[popid]["offspring"]
            ishaploid = popmakeup[popid]["ishaploid"]
            initprob = popmakeup[popid]["initprob"]
            nzorigin = popmakeup[popid]["nzorigin"]
            fdict = Dict(founders .=> 1:length(founders))
            if ishaploid
                nzorigin2 = [get(fdict, i, nothing) for i in nzorigin ]
            else
                nzorigin2 = [[get(fdict, i, nothing) for i in j ] for j in nzorigin]
            end
            for snp in 1:size(offgeno,1)
                fseq = foundergeno[snp,founders]
                in(unique(fseq),[["1"],["2"]]) && continue
                fseq2 = [get(fgenorule, i, missing) for i in fseq]
                fseq3 = vec(collect(Iterators.product(fseq2...)))
                oseq = offgeno[snp,offspring]
                issubset(unique(oseq), ["NN","N"]) && continue
                if ishaploid
                    gcount = [sum(oseq .== i) for i in ["11","22"]]
                    gfreq = [begin
                        g = [h[i] for i in nzorigin2]
                        sum(initprob .* g)
                    end for h in fseq3]
                    pv = reduce(max,[pvalue(BinomialTest(gcount[2],sum(gcount), i)) for i in gfreq])
                else
                    gcount = [sum(oseq .== i) for i in ["11","12","22"]]
                    gfreq = [begin
                        g = [sum(h[i]) for i in nzorigin2]
                        [sum(initprob[g .== i]) for i in 0:2]
                    end for h in fseq3]
                    # TODO: change into exact test
                    pv = reduce(max,[pvalue(ChisqTest(gcount, i)) for i in gfreq])
                end
                snpid = magicgeno.markermap[chr][snp,:marker]
                push!(distortions, [chrid, snpid, popid, log10(pv), join(gcount,"|"), chr,snp])
            end
        end
        if !isempty(distortions)
            distortions2= permutedims(reduce(hcat,distortions))
            distort_df = DataFrame(distortions2,colnames)
            sort!(distort_df, [:chr_index,:marker_index])
            CSV.write(outfiles[1],distort_df[!,1:5];append=true)
            filter!(row -> row.log10p < log10_siglevel, distort_df)
            isempty(distort_df) || CSV.write(outfiles[2],distort_df[!,1:5];append=true)
        end
        msg = string("chr=", chrid, ", #segregation distortion =", size(distort_final,1))
        MagicBase.printconsole(logio,verbose,msg)
    end
    msg= string("End, tused = ",
        round(time()-starttime), " seconds by filter_distortion!")
    MagicBase.printconsole(logio,verbose,msg)
    MagicBase.printconsole(logio,verbose,repeat("-",76))
    if typeof(logfile) <: AbstractString
        close(logio)
    elseif typeof(logfile) <: IO
        flush(logio)
    end
    magicgeno
end
