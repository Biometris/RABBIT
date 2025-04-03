
function magicld(genofile::AbstractString, 
    pedinfo::Union{Integer,MagicBase.JuncDist,AbstractString}=0;    
    binfile::Union{Nothing,AbstractString}=nothing,
    formatpriority::AbstractVector=["GT","AD"],
    missingstring=["NA","missing"],        
    isdepmodel::Bool=false,        
    threshcall::Real = isdepmodel ? 0.95 : 0.9,   
    seqerror::Real=0.001,
    markerthin::Integer=1,
    minlodsave::Union{Nothing, Real}=nothing,
    minldsave::Union{Nothing, Real}=nothing,
    isparallel::Bool=true,    
    workdir::AbstractString=pwd(),
    commentstring::AbstractString="##",
    outstem::AbstractString="outstem",
    logfile::Union{AbstractString,IO}= string(outstem,"_magicld.log"),
    verbose::Bool=true)    
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicld"; verbose,delim="=")
    if isa(logfile, AbstractString)
        MagicBase.printpkgst(io,verbose,"MagicLD")
    end
    ispath(workdir) || @error string(workdir, " is not a directory")
    msg = string("list of file args/options: \n",
        "genofile = ", genofile, "\n",        
        "pedinfo = ", pedinfo,"\n",
        "isdepmodel = ", isdepmodel,"\n",
        "formatpriority = ", formatpriority,"\n",
        "missingstring = ", missingstring, "\n",        
        "commentstring = ", commentstring, "\n",
        "workdir = ", workdir, "\n",        
    )
    printconsole(io,verbose,msg)   
    tused = @elapsed begin 
        # foundergeno is not used.
        magicgeno=formmagicgeno(genofile,pedinfo;
            isphysmap=false, recomrate=1.0, 
            isfounderinbred=false, formatpriority, commentstring,missingstring, workdir)
        # gsize = round(Base.summarysize(geno)/10^6,digits=3)
        # printconsole(io,verbose,string("genodata summarysize: ", gsize, "MB"))
        mem1 = round(Int, memoryuse()/10^6)
        GC.gc()
        mem2 = round(Int, memoryuse()/10^6)        
    end
    msg = string("formmagicgeno, tused=", round(tused,digits=1), "s, mem=",mem1,"|",mem2,"MB")
    MagicBase.printconsole(io,verbose,msg)  
    res=magicld!(magicgeno; binfile,
        isdepmodel, threshcall,
        seqerror, markerthin, minlodsave,minldsave, 
        isparallel,workdir, outstem,logfile=io,verbose
    )
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicld"; verbose,delim="=")    
    res
end

function magicld!(magicgeno::MagicGeno;
    binfile::Union{Nothing,AbstractString}=nothing,
    isdepmodel::Bool=false,    
    threshcall::Real = isdepmodel ? 0.95 : 0.9,   
    seqerror::Real=0.001,
    markerthin::Integer=1,
    minlodsave::Union{Nothing, Real}=nothing,
    minldsave::Union{Nothing, Real}=nothing,    
    isparallel::Bool=true,
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem",
    logfile::Union{AbstractString,IO}= string(outstem,".log"),
    verbose::Bool=true)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicld"; verbose)
    isa(logfile, AbstractString) && printpkgst(io,verbose,"MagicLD")
    isparallel = isparallel && nprocs() > 1
    ispath(workdir) || @error string(workdir, " is not a directory")
    msg = string("list of options: \n",
        "binfile = ", binfile,"\n",            
        "isdepmodel = ", isdepmodel,"\n",          
        "threshcall = ", threshcall, "\n",          
        "seqerror = ", seqerror,"\n",
        "markerthin = ", markerthin, "\n",
        "minlodsave = ", minlodsave, "\n",
        "minldsave = ", minldsave, "\n",        
        "isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose
    )
    printconsole(io,verbose,msg)
    merge_chromosome!(magicgeno)    # only ignore chrom/linagegroup
    if markerthin>1
        nsnp = sum(size.(magicgeno.markermap,1))
        submagicgeno!(magicgeno; snpsubset = 1:markerthin:nsnp)
    end
    # incorperate binfile    
    if isnothing(binfile)        
        snprepresentls = only(magicgeno.markermap)[:, :marker]
        physchromls0 = only(magicgeno.markermap)[:, :physchrom]         
        physchromls =  [ismissing(i) ? "NA" : i for i in physchromls0]    
        physposbpls0 = only(magicgeno.markermap)[:, :physposbp]                  
        physposbpls =  [ismissing(i) ? "NA" : i for i in physposbpls0]    
        snpbinls = repeat(["NA"],length(snprepresentls)) 
    else
        bindf = CSV.read(getabsfile(workdir,binfile),DataFrame; missingstring="NA",comment="##")
        gdf = groupby(bindf,:binno) # withni each group the order of rows in df is preserved
        dupebinls = Vector{Pair{String,Any}}()        
        col2str(col) = join([ismissing(i) ? "NA" : i for i in col],"||")
        for df in gdf
            represents = df[df[!,:represent] .> 0,:marker]
            isempty(represents) && continue
            if length(represents) > 1 
                msg = string("unexpected binning: ",df)                
                printconsole(io, false, "ERROR: "*msg)    
                error(msg)
            end
            push!(dupebinls, only(represents)=> (col2str(df[!,:marker]), col2str(df[!,:physchrom]), col2str(df[!,:physposbp])))
        end        
        markers = only(magicgeno.markermap)[:, :marker]
        dict = Dict(markers .=> 1:length(markers))        
        snprepresentls = first.(dupebinls)        
        snpsubset = [get(dict,i,missing) for i in snprepresentls]
        b = ismissing.(snpsubset)
        if any(b)
            msg = string("delete markers in binfile but  not genofile: ",snprepresentls[b])
            printconsole(io,verbose,msg)
            deleteat!(snprepresentls,b)
            deleteat!(snpsubset,b)
        end
        submagicgeno!(magicgeno; snpsubset)        
        ls = last.(dupebinls)
        deleteat!(ls,b)
        snpbinls = [i[1] for i in ls]
        physchromls = [i[2] for i in ls]
        physposbpls = [i[3] for i in ls]     
    end
    offinfo = copy(magicgeno.magicped.offspringinfo)
    nsnp = length(snprepresentls)
    noff = size(offinfo,1)
    msg = string("#markers=", nsnp,
        ", #individuals=", noff,
        ", 1st-individual=", offinfo[1,:individual])
    printconsole(io,verbose,msg)
    if isnothing(minlodsave)        
        minlodsave = get_minlodsave(noff, nsnp)        
        msg = string("set minlodsave = ",minlodsave)
        printconsole(io,verbose,msg)
    end
    if isnothing(minldsave)                
        minldsave = get_minldsave(nsnp)                
        msg = string("set minldsave = ",minldsave, "; r=sqrt(ld)=",round(sqrt(minldsave),digits=6))
        printconsole(io,verbose,msg)
    end    
    MagicBase.info_magicgeno(magicgeno;io,verbose)
    offformat = unique(reduce(vcat,[unique(i[!,:offspringformat]) for i=magicgeno.markermap]))    
    setdiff!(offformat,["GT"])
    MagicBase.rawgenoprob!(magicgeno; targets = ["founders","offspring"],
        seqerror, isfounderinbred=false, isoffspringinbred = isdepmodel)        
    MagicBase.rawgenocall!(magicgeno; callthreshold = threshcall, isfounderinbred=false,ishalfcall=true)
    if !isempty(intersect(offformat, ["GP", "AD"]))
        msg = string("offspringformat=",join(offformat,","), "; transformed to GT with threshcall=",threshcall)
        printconsole(io,verbose,msg)
    end
    if in("GT_phased", offformat)
        msg = string("ignore phase information; offspringformat=",join(offformat,","))
        verbose && @warn msg
        printconsole(io,false,"WARN: "*msg)
    end
    _, offspringformat = MagicBase.setunphasedgeno!(magicgeno)
    if in("GT_haplo",offspringformat)
        @error string("offspringformat=",offspringformat, ", offspring have haplotypes, not genotypes")
    end        
    tused = @elapsed  begin 
        dosegeno = MagicBase.get_dosegeno(magicgeno;isdepmodel)    
        dosegeno = permutedims(dosegeno)
        magicgeno = nothing  # save memory usage
        mem1 = round(Int, memoryuse()/10^6)
        GC.gc()
        mem2 = round(Int, memoryuse()/10^6)        
    end
    msg = string("get_dosegeno, tused=", round(tused,digits=1), "s, mem=",mem1,"|",mem2,"MB")
    printconsole(io, verbose, msg)
    nwork = nworkers()
    if isparallel && nwork>1
        if nsnp <= 10*nwork
            nseg = 1
        else
            nseg0 = round(Int, nsnp/(1000*nwork))
            nseg = (nseg0+1)*nwork
        end
    else
        nseg = nsnp <= 1000 ? 1 : div(nsnp,1000)
    end
    snprg = getsnprg(nsnp,nseg)    
    msg = string(length(snprg), " sets of snpranges")
    nwork > 1 && (msg *= string(", distributed over ", nwork, " workers"))
    printconsole(io,verbose,msg)   
    filels = [joinpath(workdir,string(outstem*"_pairwiseld_temporary",i,".txt")) for i=1:length(snprg)]
    startt = time()    
    if isparallel && nwork>1
        show_progress = false
        msgls = pmap((x,y)->get_subld(dosegeno, x, minlodsave,minldsave, show_progress,y,verbose), snprg,filels)
    else
        show_progress = true
        msgls = map((x,y)->get_subld(dosegeno, x, minlodsave,minldsave, show_progress,y,verbose), snprg,filels)
    end
    for msg = msgls
        write(io,string(msg,"\n"))
    end
    flush(io)
    mem = round(Int, memoryuse()/10^6)
    msg = string("tused=",round(time()-startt, digits=1), "s for all markers, mem=",mem, "MB")
    printconsole(io, verbose, msg)
    outfile = string(outstem,"_magicld.csv.gz")    
    outfile2 = joinpath(workdir,outfile)
    tused = @elapsed GZip.open(outfile2,"w") do outio
        initial = "RABBIT"
        delim = ','
        MagicBase.appenddf(outio, offinfo; delim,initial,dfname="offspringinfo")
        markerinfo = DataFrame(markerno = 1:length(snprepresentls), marker_represent=snprepresentls,
            marker_bin = snpbinls,
            physchrom_bin=physchromls,
            physposbp_bin=physposbpls,
            )
        MagicBase.appenddf(outio, markerinfo; delim,initial,dfname="markerinfo")
        write(outio,initial,delim,"pairwiseld\n")
        write(outio, "marker1,marker2,ld_r2,ld_lod\n")       
        for file in reverse(filels)
            if isfile(file)
                write(outio,read(file))
                flush(outio)
            end
        end        
    end
    msg = string("magicld file: ", outfile, ", tused=",round(tused, digits=1), "s in merging outfiles")
    printconsole(io,verbose,msg)
    rm.(filels;force=true)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicld"; verbose)
    outfile
end

function get_minlodsave(nind::Integer, nsnp::Integer)    
    if nsnp < 10000 
        minlodsave = round(1.0+0.1*nsnp/1000,digits=3)        
    else
        slope_save = nind >= 1000 ? 0.5*(1.0+log10(nind/1000)) : 0.5*nind/1000    
        minlodsave = 2.0 + slope_save*log(2,nsnp/10^4)        
        minlodsave = round(min(5.0,minlodsave),digits=3)
    end
    minlodsave
end

function get_minldsave(nsnp::Integer)
    if nsnp <= 10000
        corr = max(0.01,round(0.01*nsnp/1000,digits=3))                
    elseif nsnp <= 640000
        corr = 0.1 + 0.05*log(2,nsnp/10000)
    else
        corr = 0.4
    end
    round(corr^2,digits=6)
end

function getsnprg(nsnp::Integer,nseg::Integer)
    ls = accumulate(+,nsnp-1:-1:1)
    ls2 = round.(Int,(1:(nseg-1)) .* (last(ls)/nseg))
    ls3 = [findfirst(x->x>i, ls) for i=ls2]
    pushfirst!(ls3,1)
    push!(ls3,nsnp)
    ls3 = unique(ls3)
    res = [ls3[i]:(ls3[i+1]-1) for i=1:length(ls3)-1]
    reverse(res)
end
