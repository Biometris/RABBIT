
function filter_founder_progeny!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    isfounderinbred::Bool=true,
    minnprogeny::Real = 1,
    outstem::AbstractString= "outstem",
    logfile::Union{AbstractString,IO} = outstem*"_purify_founder.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    minnprogeny < 1 && return magicgeno
    logio = MagicBase.set_logfile_begin(logfile, workdir, "filter_founder_progeny!"; verbose)
    model = MagicBase.reset_model(magicgeno.magicped,model;io=logio,verbose)
    msg = string("list of options: \n",
        "model = ", model, "\n",
        "isfounderinbred = ", isfounderinbred, "\n",
        "minnprogeny = ", minnprogeny, "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)
    MagicBase.setunphasedgeno!(magicgeno)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
    # MagicBase.info_missing(magicgeno;io=logio,verbose)
    # del founder with nprogeny < minnprogeny
    del_off = []
    f2off = MagicBase.get_founder2offspring(magicgeno.magicped)
    offls = collect(values(f2off))
    nprogenyls = length.(offls)
    while true
        b = length.(offls) .< minnprogeny
        if any(b)
            append!(del_off,offls[b]...)
            unique!(del_off)
            offls = offls[.!b]
            for i in offls
                setdiff!(i,del_off)
            end
            offls = offls[.!isempty.(offls)]
            isempty(offls) && break
        else
            break
        end
    end
    indls = magicgeno.magicped.offspringinfo[!,:individual]
    offkeep = findall(.![i in del_off for i in indls])
    del_offspring!(magicgeno,offkeep;io=logio,verbose)
    # export
    fls = collect(keys(f2off))
    newf2off = MagicBase.get_founder2offspring(magicgeno.magicped)
    newfls = collect(keys(newf2off))
    b = [i in newfls for i in fls]
    new_nprogenyls = [haskey(newf2off,i) ? length(newf2off[i]) : 0 for i in fls]
    resdf = DataFrame(founder=fls,
        nprogeny=nprogenyls,
        nprogeny_after=new_nprogenyls,
        ndel = new_nprogenyls .- nprogenyls, keep = b)
    outfile = outstem*"_founder_nprogeny.csv"
    outfile2 = joinpath(workdir,outfile)
    open(outfile2,"w") do outio
        commentstring = "##"
        msg = commentstring*"col_1, founder, founder ID\n"
        msg *= commentstring*"col_2, nprogeny, number of progeny for the founder\n"
        msg *= commentstring*"col_3, nprogeny_after, number of progeny for the founder after filtering\n"
        msg *= commentstring*"col_4, ndel, number of progeny removed for the founder\n"
        msg *= commentstring*"col_5, keep, the founder is removed if keep = false"
        write(outio, msg, "\n")
        CSV.write(outio,resdf; header=true, append=true)
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"purify_founder"; verbose)
    magicgeno
end
