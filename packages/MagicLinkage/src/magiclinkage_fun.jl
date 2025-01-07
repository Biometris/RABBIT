
function magiclinkage(genofile::AbstractString, pedinfo::Union{MagicBase.JuncDist,AbstractString};
    ldfile::Union{Nothing,AbstractString}=nothing,
    formatpriority::AbstractVector=["GT","AD"],    
    isfounderinbred::Bool=true,
    model::AbstractString="jointmodel",
    likeparameters::LikeParameters=LikeParameters(),   
    israndallele::Bool=true, 
    threshcall::Real = model == "depmodel" ? 0.95 : 0.9,   
    snpthin::Integer=1,    
    byfounder::Integer=0,    
    minlodsave::Union{Nothing, Real}=nothing,
    maxrfsave::Real=1.0, # rf scaled from 0 to 1
    commentstring::AbstractString="##",
    isparallel::Bool=true,
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem",
    logfile::Union{AbstractString,IO}= string(outstem,"_magiclinkage.log"),
    verbose::Bool=true)
    starttime = time()
    ispath(workdir) || string(workdir, " is not a directory")
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magiclinkage"; verbose,delim="=")
    msg = string("list of file args: \n",
        "genofile = ", genofile, "\n",
        "pedinfo = ", pedinfo, "\n",        
        "formatpriority = ", formatpriority, "\n",        
        "isfounderinbred = ", isfounderinbred, "\n",        
        "commentstring = ", commentstring, "\n",
        "workdir = ",workdir,"\n")
    printconsole(logio,verbose,msg)
    tused = @elapsed begin 
        magicgeno=formmagicgeno(genofile, pedinfo;
            formatpriority, isfounderinbred, commentstring,workdir)        
        mem1 = round(Int, memoryuse()/10^6)
        GC.gc()
        mem2 = round(Int, memoryuse()/10^6)        
    end
    msg = string("formmagicgeno, tused=", round(tused,digits=1), 
        "seconds, mem=",mem1,"|",mem2,"MB")
    MagicBase.printconsole(logio,verbose,msg)            
    linkagefile=magiclinkage!(magicgeno;        
        ldfile,model, likeparameters, israndallele, threshcall,
        snpthin, byfounder, minlodsave, maxrfsave,
        isfounderinbred, 
        isparallel,workdir,outstem,
        logfile=logio,verbose
    )
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magiclinkage";verbose,delim="=")
    linkagefile
end


function magiclinkage!(magicgeno::MagicGeno;    
    ldfile::Union{Nothing,AbstractString}=nothing,
    isfounderinbred::Bool=true,
    model::AbstractString="jointmodel",
    likeparameters::LikeParameters=LikeParameters(),   
    israndallele::Bool=true, 
    threshcall::Real = model == "depmodel" ? 0.95 : 0.9,   
    snpthin::Integer=1,    
    byfounder::Integer=0,    
    minlodsave::Union{Nothing, Real}=nothing,
    maxrfsave::Real=1.0, # rf scaled from 0 to 1
    isparallel::Bool=true,
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem",
    logfile::Union{AbstractString,IO}= string(outstem,"_magiclinkage.log"),
    verbose::Bool=true)
    starttime = time()
    ispath(workdir) || string(workdir, " is not a directory")
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magiclinkage";verbose,delim="-")    
    msg = string("list of options: \n",
        "ldfile = ", ldfile, "\n",        
        "isfounderinbred = ", isfounderinbred, "\n",        
        "model = ", model, "\n",
        "likeparameters = ", likeparameters, "\n",        
        "israndallele = ", israndallele, "\n",        
        "threshcall = ", threshcall, "\n",        
        "snpthin = ", snpthin, "\n",        
        "byfounder = ", byfounder, "\n",        
        "minlodsave = ", minlodsave, "\n",
        "maxrfsave = ", maxrfsave, "\n",
        "isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)    
    merge_chromosome!(magicgeno)
    if snpthin > 1
        @info string("take every ",snpthin, "-th marker")
        snpsubset = 1:snpthin:size(only(magicgeno.markermap),1)
        submagicgeno!(magicgeno;snpsubset)
    end
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
    
    model = MagicBase.reset_model(magicgeno.magicped,model; io=logio,verbose)
    MagicBase.reset_juncdist!(magicgeno.magicped,model; io=logio,verbose,isfounderinbred)    
    isdepmodel = model == "depmodel"    
    seqerror = MagicBase.get_seqerror(likeparameters)
    offformat = unique(reduce(vcat,[unique(i[!,:offspringformat]) for i=magicgeno.markermap]))
    setdiff!(offformat,["GT"])
    MagicBase.rawgenoprob!(magicgeno; seqerror,isfounderinbred,isoffspringinbred= isdepmodel)    
    MagicBase.rawgenocall!(magicgeno; callthreshold = threshcall, isfounderinbred, ishalfcall=true)
    if issubset(offformat, ["GP", "AD"])
        msg = string("offspringformat=",join(offformat,","), "; transformed to GT with threshcall=",threshcall)
        printconsole(logio,verbose,msg)
    end
    nsnp = size(only(magicgeno.markermap),1)    
    if isnothing(minlodsave)
        nind = size(magicgeno.magicped.offspringinfo,1)            
        minlodsave = MagicLD.get_minlodsave(nind, nsnp)
        msg = string("set minlodsave = ", round(minlodsave,digits=3))
        printconsole(logio,verbose,msg)
    end
    snppairls = get_snppairls!(magicgeno,ldfile; minlodsave, workdir,io=logio,verbose)    
    popmakeup, recomprior = calpairwiseprior(magicgeno.magicped, model;
        isfounderinbred,isautosome=true)
    if model == "jointmodel"
        for (popid,value) in popmakeup
            msg = string("pop=", popid,
                ", ibd=", round.(value["inbreedingcoef"],digits=4),
                ", jointmodel=>", value["model"])
            printconsole(logio,verbose,msg)
        end
    end
    chr = 1
    fhaploset = calfhaploprior_loci(magicgeno,chr)
    epsf = MagicBase.get_foundererror(likeparameters)
    epso = MagicBase.get_offspringerror(likeparameters)
    condlike=MagicReconstruct.precompute_condike(epsf,epso; isoffphased=false,israndallele)
    # offcode for haplo: dict=Dict(["N"=>1,"1"=>2, "2"=>3])
    # offcode for diplo: dict = Dict(["NN", "1N", "2N", "11", "12", "22"] .=> 1:6)
    # offcode: noffspring x nmarker after permutedims
    missingcode = 1
    offcode = permutedims(MagicReconstruct.precompute_offcode(magicgeno,chr,popmakeup))        
    snpls = only(magicgeno.markermap)[:,:marker]            
    physchromls =  [ismissing(i) ? "NA" : i for i in only(magicgeno.markermap)[:,:physchrom]] 
    physposbpls =  [ismissing(i) ? "NA" : i for i in only(magicgeno.markermap)[:,:physposbp]] 
    genomtx = magicgeno.offspringgeno[chr]
    formatvec = magicgeno.markermap[chr][!,:offspringformat]
    nmissingls = sum(MagicBase.getismissing(genomtx, formatvec),dims=2)[:,1]    
    offinfo = copy(magicgeno.magicped.offspringinfo)
    nsnp = length(snpls)
    noff = size(offinfo,1)
    msg = string("#markers=", nsnp,
        ", #individuals=", noff,
        ", 1st-individual=", offinfo[1,:individual])
    printconsole(logio,verbose,msg)
    magicgeno = nothing
    mem1 = round(Int, memoryuse()/10^6)
    GC.gc()
    mem2 = round(Int, memoryuse()/10^6)
    printconsole(logio,verbose,string("memoryuse=", mem1,"|",mem2,"MB"))
    nwork = nworkers()    
    subpairs = getsubpairs(snppairls,nwork)
    if nwork > 1
        msg = string(length(subpairs), " sets of snppairs, distributed over ", nwork, " workers")
        printconsole(logio,verbose,msg)
    end
    # println("offcode summarysize: ", Base.summarysize(offcode)*10^(-6),"MB")       
    # println("recomprior=",recomprior) 
    # println("popmakeup=",popmakeup)
    startt = time()    
    filels = [joinpath(workdir,string(outstem*"_magiclinkage_temporary",i,".txt"))
        for i=1:length(subpairs)]
    if isparallel && nwork>1
        show_progress = false
        msgls = pmap((x,y)->sublinkage(x,offcode, condlike, fhaploset, popmakeup, recomprior;
            missingcode,byfounder,minlodsave,maxrfsave,outfile=y,show_progress,verbose), subpairs,filels)
        empty!(lru_jointlike)
    else
        show_progress = true
        msgls = map((x,y)->sublinkage(x,offcode, condlike, fhaploset, popmakeup, recomprior;
            missingcode,byfounder,minlodsave,maxrfsave,outfile=y,show_progress,verbose), subpairs,filels)        
        empty!(lru_jointlike)
    end
    for msg = msgls
        write(logio,string(msg,"\n"))
    end
    flush(logio)
    mem1 = round(Int, memoryuse()/10^6)
    GC.gc()
    mem2 = round(Int, memoryuse()/10^6)
    msg = string("tused=", round(time()-startt,digits=1), "s for all markers",
        ", mem=",mem1,"|",mem2, "MB")
    printconsole(logio,verbose,msg)
    outfile = string(outstem,"_magiclinkage.csv.gz")    
    outfile2 = joinpath(workdir,outfile)
    tused = @elapsed GZip.open(outfile2,"w") do outio        
        initial = "RABBIT"
        delim = ','
        MagicBase.appenddf(outio, offinfo; delim,initial,dfname="offspringinfo")
        markerinfo = DataFrame(markerno = 1:length(snpls), marker=snpls,
            physchrom=physchromls, physposbp=physposbpls,nmissing = nmissingls)
        MagicBase.appenddf(outio, markerinfo; delim,initial,dfname="markerinfo")
        write(outio,initial,delim,"pairwiselinkage\n")
        write(outio, "marker1, marker2, linkage_rf, linkage_lod\n")
        for file in filels
            if isfile(file)
                write(outio,read(file))
                flush(outio)
            end
        end
    end
    msg = string("magiclinkage file: ", outfile, ", tused=",round(tused, digits=1), "s in merging outfiles")
    printconsole(logio,verbose,msg)
    rm.(filels;force=true)    
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magiclinkage";verbose,delim="-")
    outfile
end

function get_snppairls!(magicgeno::MagicGeno,
    ldfile::Union{Nothing,AbstractString}=nothing;
    minlodsave::Real=2.0,    
    workdir::AbstractString=pwd(),
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    markers = Vector(only(magicgeno.markermap)[!,:marker])
    nsnp = length(markers)
    if isnothing(ldfile) 
        pairindicator = sparse(triu(trues(nsnp,nsnp),1))        
    else
        pairwiseld = MagicBase.read_ld(ldfile; workdir)
        # only one linagegroup in magicgeno
        dict = Dict(markers .=> 1:length(markers))
        snpsubset = [get(dict,i,nothing) for i in pairwiseld.markers]
        b  = isnothing.(snpsubset)
        if any(b)
            msg = string(sum(b), "markers in ldfile but not in magicgeno")
            printconsole(io, verbose, "WARNING: "*msg)
            @warn msg
        end        
        snpsubset = sort(snpsubset[.!b])
        MagicBase.submagicgeno!(magicgeno;snpsubset)
        # re-order pairindicator
        commonmarkers = markers[snpsubset]
        commonmarkers == only(magicgeno.markermap)[!,:marker] || @error "inconsistent markers"
        markerdiff = setdiff(markers, commonmarkers)
        if !isempty(markerdiff)
            msg = string(length(markerdiff), " markers in genofile but not in ldfile")
            printconsole(io, verbose, "WARNING: "*msg)
            @warn msg            
        end
        dict = Dict(pairwiseld.markers .=> 1:length(pairwiseld.markers))
        ii = [get(dict,i,nothing) for i in commonmarkers]        
        pairindicator = triu(pairwiseld.ldlod[ii,ii] .> minlodsave,1)   
    end
    pairnz = findnz(pairindicator)
    pairls = tuple.(pairnz[1], pairnz[2])
    nsnp = size(pairindicator,1)  
    totpair = div(nsnp*(nsnp-1),2)        
    npair = length(pairls)
    msg = string(npair," out of ", totpair, " pairs, fraction = ",
        round(100*npair/totpair,digits=4), "%")
    printconsole(io,verbose,msg)  
    pairls
end

function sublinkage(subpairls::AbstractVector,offcode::AbstractMatrix,
    condlike::NamedTuple, fhaploset::AbstractVector,
    popmakeup::AbstractDict, recomprior::AbstractDict;
    missingcode,byfounder::Integer,
    minlodsave::Real,maxrfsave::Real, show_progress::Bool,
    outfile::AbstractString,verbose::Bool)
    startt = time()    
    isfounderinbred = !(typeof(fhaploset[1][1][1]) <: AbstractVector)
    open(outfile,"w") do io
        progress = Progress(length(subpairls);
            enabled=show_progress,
            desc=string("worker=",myid()," linkage..."))        
        for (snp1,snp2) in subpairls
            rf,lod = link2loci(snp1,snp2, offcode, condlike, fhaploset, popmakeup,recomprior;
                missingcode,byfounder,isfounderinbred)
            rf = round(rf,digits=5)
            lod = round(lod,digits=2)
            # if snp1 == 1
            #     println("snp1=",snp1,",snp2=",snp2,",rf=",rf, ",lod=",lod)
            # end
            if lod ≥ minlodsave && rf ≤ maxrfsave
                msg = string(join([snp1,snp2],","), ",", join([rf,lod],","))
                write(io, msg, "\n")
            end
            next!(progress)
        end
        flush(io)
    end
    mem1 = round(Int, memoryuse()/10^6)
    GC.gc()
    mem2 = round(Int, memoryuse()/10^6)
    msg = string("worker=",myid(),", #snppairs=",length(subpairls), 
        ", t=", round(time()-startt,digits=1), "s",
        ", mem=",mem1,"|",mem2,"MB")    
    verbose && (@info msg)    
    msg
end

function link2loci(loc1::Integer, loc2::Integer,
    offcode::AbstractMatrix, condlike::NamedTuple,fhaploset::AbstractVector,
    popmakeup::AbstractDict, recomprior::AbstractDict;
    missingcode,byfounder::Integer,isfounderinbred::Bool)
    loglmax = -Inf
    rfmax = 1.0
    lodmax = 0
    # count12 = calcount12(view(offcode,:,loc1), view(offcode,:,loc2), popmakeup; missingcode)
    count12 = calcount12(offcode[:,loc1], offcode[:,loc2], popmakeup; missingcode)
    isempty(count12) && return (rfmax, lodmax)    
    fhaploset12 = fhaploset[[loc1,loc2]]    
    if !isfounderinbred
        # absoluate phase does not affect recom_fraction
        fhaploset12[1] = [unique(sort.(i)) for i in fhaploset12[1]]
    end
    fhaplo12 = randinitfhaplo(fhaploset12)    
    vcatfhaplo12 = isfounderinbred ? fhaplo12 : [reduce(vcat,i) for i in fhaplo12]
    likerecom = Dict([popid=>cal_likerecom_popid(popid, vcatfhaplo12,count12,
        condlike, popmakeup,recomprior) for popid in keys(count12)])
    # calculate findexlist
    subfounders = sort(unique(reduce(vcat,[popmakeup[i]["founder"] for i in keys(count12)])))
    lenls = reduce(vcat,[length.(i[subfounders]) for i in fhaploset12])
    lenls = lenls[lenls .>= 2]
    nphase = length(lenls)>10 ? 2^10 : reduce(*,lenls) # avoid overflow, return result 0!
    nfounder = length(first(fhaplo12))
    if byfounder == -1 
        findexlist = [1:nfounder]        
    else
        if nphase <= 64
            findexlist = [1:nfounder]            
        else
            fmissls = mean([length.(i) .> 1 for i in fhaploset12])
            findexlist = MagicBase.getfindexlist(byfounder,fmissls, popmakeup;defaultby=2)            
        end
    end
    findexlist = [intersect(i,subfounders) for i in findexlist]
    filter!(x->!isempty(x),findexlist)    
    itmax = length(findexlist) ==1 ? 1 : 50    
    for it in 1:itmax
        oldloglmax  = loglmax        
        oldfhaplo12 = deepcopy(fhaplo12)
        for findex in findexlist
            fhaplo12ls = getfhaplo12ls(findex, fhaplo12,fhaploset12)
            if !isinf(loglmax)            
                if length(fhaplo12ls)==1 
                    continue
                else
                    setdiff!(fhaplo12ls,[fhaplo12])  
                end
            end  
            popidls = MagicImpute.popfromfindex(findex,popmakeup)     
            intersect!(popidls, keys(count12))         
            isempty(popidls) && continue   
            res = [begin
                up_likerecom!(likerecom, popidls, fhaplo12ls[i],count12,condlike, popmakeup,recomprior,isfounderinbred)
                inferrecomfrac(likerecom)
            end for i in eachindex(fhaplo12ls)]
            loglls = [i[3] for i in res]
            logl = reduce(max,loglls)
            if logl > loglmax
                pos = rand(findall(loglls .≈ logl))
                fhaplo12 = fhaplo12ls[pos]
                rfmax, lodmax, loglmax= res[pos]
            end
            up_likerecom!(likerecom, popidls, fhaplo12,count12,condlike, popmakeup,recomprior,isfounderinbred)
            # println("it=",it, ",findex=", findex,",nphase=", length(res), ", [rfmax, loglmax]= ",round.([rfmax, loglmax],digits=4))
        end
        if oldfhaplo12 == fhaplo12 || abs(oldloglmax  - loglmax)<0.1 || it == itmax
            # println("it=", it, ", loc1=", loc1, ",loc2=", loc2, ", nphase=",nphase,",findexlist=",findexlist, 
            #     ",[oldlogl,logl]=",[oldloglmax,loglmax])        
            break
        end        
    end
    rfmax, lodmax
end


function inferrecomfrac(likerecom::AbstractDict)
    function f(logr)
        res = zero(0.0)
        r = exp(logr)
        for (like, cc) in values(likerecom)
            ncol = size(like,2)
            rp = r .^ (0:ncol-1)
            res += dot(log.(like * rp),cc)
        end
        res
    end
    loglnull = f(log(1.0))
    logl05 = f(log(0.5))
    logl0 = f(-Inf)    
    est = if logl0 ≈ logl05 ≈ loglnull
        # nit = 1
        [NaN,0, loglnull]
    else
        # r is a scaled recombination fraction [0,1]
        rmin, rmax = 10^(-5), 1.0
        if logl0 > logl05 > loglnull
            rmin, rmax = 10^(-5), 0.5
        elseif loglnull > logl05 > logl0
            rmin, rmax = 0.5, 1.0
        end
        x, logl, _ = MagicBase.brentMax(f,log(rmin), log(rmax);
            xstart=log((rmin+rmax)/2), precisiongoal=3,accuracygoal=3)
        r = exp(x)
        lod = (logl - loglnull) * log(10,ℯ)
        # nit = his[end,1]
        # println("r=",r, ";his=", DataFrame(his,:auto))
        [r,lod,logl]
    end    
    est
end

function up_likerecom!(likerecom::AbstractDict, popidls::AbstractVector,
    fhaplo12::AbstractVector,count12::AbstractDict,condlike::NamedTuple,
    popmakeup::AbstractDict, recomprior::AbstractDict,isfounderinbred::Bool)    
    vcatfhaplo12 = isfounderinbred ? fhaplo12 : [reduce(vcat,i) for i in fhaplo12]
    for popid in popidls
        likerecom[popid] = cal_likerecom_popid(popid, vcatfhaplo12,count12,
            condlike, popmakeup,recomprior)
    end
    return
end

function cal_likerecom_popid(popid::AbstractString, vcatfhaplo12::AbstractVector,
    count12::AbstractDict,condlike::NamedTuple,
    popmakeup::AbstractDict, recomprior::AbstractDict)
    nzorigin = popmakeup[popid]["nzorigin"]
    isnonibd = allunique.(nzorigin)
    ishaploid = popmakeup[popid]["ishaploid"]
    hashcode = popmakeup[popid]["hashcode"]
    fderive1 = calfderive(vcatfhaplo12[1],nzorigin,ishaploid)
    fderive2 = calfderive(vcatfhaplo12[2],nzorigin,ishaploid)
    like12 = Vector{Vector{Float64}}()
    cc12 = Vector{Int}()
    for (g12, cc)  in count12[popid]
        like = caljointlike_lru(g12..., fderive1, fderive2,hashcode, ishaploid,
            isnonibd,condlike,recomprior)
        push!(like12, like)
        push!(cc12, cc)
    end
    reduce(hcat,like12)', cc12
end

# cache speeds up for indepmeodel but not much for depmodel;
# local-brent is time consuming
const lru_jointlike = LRU{UInt, Vector{Float64}}(maxsize=10^6)

function caljointlike_lru(g1::Integer, g2::Integer,
    fderive1::AbstractVector, fderive2::AbstractVector,
    hashcode::Integer, ishaploid::Bool, isnonibd::AbstractVector,
    condlike::NamedTuple,recomprior::AbstractDict)
    # condlike and recomprior are constant
    lru_key = hash([g1,g2,fderive1, fderive2,hashcode, ishaploid, isnonibd])
    get!(lru_jointlike, lru_key) do
        caljointlike(g1,g2, fderive1, fderive2,hashcode, ishaploid,isnonibd,
            condlike,recomprior)
    end
end

function caljointlike(g1::Integer, g2::Integer,
    fderive1::AbstractVector, fderive2::AbstractVector,
    hashcode::Integer, ishaploid::Bool, isnonibd::AbstractVector,
    condlike::NamedTuple,recomprior::AbstractDict)
    like1 = calgenolike(g1, fderive1, ishaploid, isnonibd,condlike)
    like2 = calgenolike(g2, fderive2, ishaploid, isnonibd,condlike)
    [like1' * (coeff * like2) for coeff in recomprior[hashcode]]
end

function calgenolike(gcode::Integer, fderive::AbstractVector,
    ishaploid::Bool,isnonibd::AbstractVector, condlike::NamedTuple)
    if ishaploid
        condlike.haplo[gcode, fderive]
    else
        likediplo = condlike.diplo
        like = zeros(length(isnonibd))
        like[.!isnonibd] .= likediplo.ibd[gcode,fderive[.!isnonibd]]
        like[isnonibd] .= likediplo.nonibd[gcode,fderive[isnonibd]]
        like
    end
end

function calfderive(fhaplo::AbstractVector, nzorigin::AbstractVector, ishaploid)
    # fhaplo contains three possible alleles (integer codes): 0(missing), 1, 2
    if ishaploid
        # [0,1,2] .=> 1:3
        fhaplo[nzorigin] .+1
    else
        # nfgl = length(fhaplo)
        derrule=Dict(["00", "01", "10", "02", "20", "11", "12", "21", "22"] .=> 1:9)
        fderive = [get(derrule,join(fhaplo[i]),-1) for i in nzorigin]
        if -1 in fderive
            error("there exist unknow founder alleles")
        end
        fderive
    end
end


function getfhaplo12ls(findex::AbstractVector,fhaplo12::AbstractVector,
    fhaploset12::AbstractVector)
    haploset1, haploset2 = [begin 
        set1 = vec(collect(Iterators.product(fhaploset12[i][findex]...)))
        haploset1 = [begin
            fhaplo = copy(fhaplo12[i])
            fhaplo[findex] .= seg
            fhaplo
        end for seg in set1]
        haploset1
    end for i in eachindex(fhaploset12)]
    [[i,j] for i in haploset1 for j in haploset2]
end

# function getouterlist(v1::AbstractVector)
#     [[i] for i in v1]
# end

# function getouterlist(v1::AbstractVector,v2::AbstractVector,x...)
#     res=[vcat(i,j) for i=v1 for j=v2]
#     foldl(getouterlist,x,init=res)
# end

# function getfhaplo12ls(findex::AbstractVector,fhaplo12::AbstractVector,
#     fhaploset12::AbstractVector)
#     set1 = getouterlist(fhaploset12[1][findex]...)
#     haploset1 = [begin
#         fhaplo = copy(fhaplo12[1])
#         fhaplo[findex] = seg
#         fhaplo
#     end for seg in set1]
#     set2 = getouterlist(fhaploset12[2][findex]...)
#     haploset2 = [begin
#         fhaplo = copy(fhaplo12[2])
#         fhaplo[findex] = seg
#         fhaplo
#     end for seg in set2]
#     [[i,j] for i in haploset1 for j in haploset2]
# end

function randinitfhaplo(fhaploset12::AbstractVector)
    [rand.(i) for i in fhaploset12]
end

function calcount12(geno1::AbstractVector, geno2::AbstractVector,
    popmakeup::AbstractDict; missingcode="NN")
    res = Dict{String,Any}()
    for popid = keys(popmakeup)
        offls = popmakeup[popid]["offspring"]
        geno12 = [missingcode in [geno1[i],geno2[i]] ? missing : [geno1[i],geno2[i]]
            for i in offls]
        count12 = countmemb(skipmissing(geno12))
        if length(count12)>1
            # || (length(count12)==1 && only(count12).second<10)
            push!(res, popid => count12)
        end
    end
    res
end

function countmemb(itr)
    d = Dict{eltype(itr), Int}()
    for val in itr
        d[val] = get(d, val, 0) + 1
    end
    d
end

function getsubpairs(pairls::AbstractVector,nwork::Integer)
    n = length(pairls)
    if n < 100*nwork
        nseg =1
    else
        nseg0 = round(Int, n/(2*(10^5)*nwork))
        # nseg0 = round(Int, n/(1*(10^5)*nwork))
        nseg = (nseg0+1)*nwork
    end
    l = div(n,nseg)+1
    collect(partition(pairls,l))
end
