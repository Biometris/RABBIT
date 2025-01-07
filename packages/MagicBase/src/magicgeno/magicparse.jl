
function magicparse(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};    
    isbreedped::Union{Nothing,Bool}=nothing,
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0,
    delmultiallelic::Bool=true,
    fixed_nself::Integer = 10,
    commentstring::AbstractString="##",
    missingstring::AbstractString="NA",
    outstem::Union{Nothing,AbstractString} = "outstem",
    outext::AbstractString=in(last(MagicBase.split_allext(genofile)),[".vcf",".vcf.gz"]) ? ".vcf.gz" : ".csv.gz",
    logfile::Union{Nothing, AbstractString,IO} = (isa(outstem, AbstractString) ? outstem : "")*"_magicparse.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magicparse"; verbose,delim="=")
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(logio,verbose,msg)
        MagicBase.printpkgst(logio,verbose,"MagicBase")
    end
    if isnothing(isbreedped)
        fileext = last(split_allext(string(pedinfo)))
        if !isempty(fileext)
            inpedfile = getabsfile(workdir,pedinfo)
            isfile(inpedfile) || error(string("pedfile=",inpedfile, " does not exist"))
            isbreedped = !get_israbbitped(inpedfile)
        end
    end
    MagicBase.check_infiles(genofile,pedinfo; isbreedped,commentstring,workdir,io=logio,verbose)
    msg = string("list of args: \n",
        "genofile = ", genofile,"\n",
        "pedinfo = ", pedinfo,"\n",
        "isfounderinbred = ", isfounderinbred,"\n",
        "formatpriority = ", formatpriority,"\n",
        "isphysmap = ", isphysmap,"\n",
        "recomrate = ", recomrate,"\n",
        "delmultiallelic = ", delmultiallelic, "\n",
        "fixed_nself = ", fixed_nself,"\n",
        "workdir = ",workdir,"\n",
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",
        "outext = ",outext,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)    
    if isa(pedinfo, Integer)
        magicgeno = formmagicgeno(genofile,pedinfo;
            isfounderinbred,formatpriority, isphysmap, recomrate,
            commentstring, missingstring = unique(vcat(missingstring,"missing")), 
            workdir)
    else
        # if pedinfo is a pedfile, it is in format of BASF format
        if isbreedped
            newpedinfo = joinpath(workdir,isnothing(outstem) ? "_magicped.csv" : outstem*"_magicparse_ped.csv")
            newpedinfo = joinpath(workdir,outstem*"_magicparse_ped.csv")
            magicped = parsebreedped(pedinfo; fixed_nself,
                outfile = newpedinfo, commentstring, workdir)
        else
            magicped = readmagicped(pedinfo; commentstring, workdir)
            newpedinfo = pedinfo
        end        
        genodf,commentlines = readgenodf(genofile,newpedinfo; isfounderinbred,formatpriority,
            delmultiallelic,commentstring, missingstring, workdir,
            logfile=logio, verbose)
        isbreedped && isnothing(outstem) && rm(newpedinfo;force=true)
        magicgeno = formmagicgeno!(genodf, magicped;
            isfounderinbred, isphysmap,recomrate)        
        if !isempty(commentlines)
            filecomment = DataFrame(sourcefile=getabsfile(workdir,genofile),
                commentlines=commentlines)
            pushmisc!(magicgeno,"filecomment"=>filecomment)
        end        
    end
    info_magicgeno(magicgeno; io=logio, verbose)
    if !isnothing(outstem)
        outstem *= "_magicparse"
        if isbreedped
            outpedfile = outstem*"_ped.csv" #keep it same as newpedinfo if it is pedifle
            tused = @elapsed savemagicped(outpedfile, magicgeno.magicped; workdir, delim=',')
            msg = string("save pedfile in ",outpedfile, ", tused =",round(tused,digits=1), "s")
            MagicBase.printconsole(logio,verbose,msg)
        end
        if size(magicgeno.magicped.founderinfo,1) <= 50
            # TODO: too large pedidgree results in error in plotting
            fn = string(outstem, "_ped.png")
            gdesign=MagicBase.plotmagicped(magicgeno.magicped;
                isfounderinbred,outfile = joinpath(workdir,fn))
            msg = isa(gdesign, AbstractString) ? gdesign : string("save design plot in ", fn)
            MagicBase.printconsole(logio,verbose,msg)
        end
        outgenofile = outstem*outext
        tused = @elapsed savegenodata(outgenofile, magicgeno;workdir,commentstring,missingstring)
        msg = string("save genofile in ",outgenofile, ", tused =",round(tused,digits=1), "s")
        MagicBase.printconsole(logio,verbose,msg)
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magicparse"; verbose,delim="=")
    magicgeno
end
