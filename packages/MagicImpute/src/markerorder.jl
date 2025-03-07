
# snporder,priorprocess, are modified
function updateorder!(snporder::AbstractVector,
    chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;    
    offspringexcl::AbstractVector, 
    orderactions::AbstractVector, 
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},    
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    israndallele::Bool,     
    issnpGT::AbstractVector, 
    decodetempfile::AbstractString,
    slidewinsize::Union{Integer,Distribution}=10,    
    maxwinsize::Union{Nothing,Integer}=nothing,    
    reversechr::Bool=false,    
    temperature::Real=0)
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end
    nsnp = length(snporder)
    logbook=Vector()
    if isnothing(maxwinsize) 
        maxwinsize = nsnp < 100 ? div(nsnp,2) : max(20,div(nsnp,10))
    end    
    startt=time()
    # actionset = ["permute", "inverse11","inverse00", "inverse01","inverse10"]    
    for proptype in orderactions
        loglhis, accepthis, winsizehis =segementorder!(snporder, proptype,
            chrfhaplo,chroffgeno, popmakeup,priorprocess;
            offspringexcl, slidewinsize, maxwinsize, chrneighbor=nothing,
            epsf, epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
            israndallele, issnpGT,decodetempfile,temperature)
        accept = mean(accepthis)
        tused = round(time()-startt,digits=1)
        temp = [round(accept,digits=4), round(loglhis[end],digits=4)]
        avgwinsize = isa(slidewinsize,Integer) ? slidewinsize : round(mean(winsizehis),digits=1)
        temp = vcat([string("#", proptype,"_Poisson"), avgwinsize,
            length(accepthis)-round(Int,avgwinsize)+1, tused], temp)
        push!(logbook, temp)    
    end
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end
    logbook
end

# snporder,priorprocess, are modified
function updateorder_neighbor!(snporder::AbstractVector,
    chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    offspringexcl::AbstractVector, 
    chrneighbor::AbstractDict,    
    orderactions::AbstractVector, 
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},    
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},  
    israndallele::Bool,     
    issnpGT::AbstractVector, 
    decodetempfile::AbstractString,
    maxwinsize::AbstractVector,
    reversechr::Bool = true,
    temperature::Real=0)
    logbook=Vector()
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end        
    # orderactions = ["inverse11","inverse01"]
    for p in eachindex(orderactions)
        proptype = orderactions[p]        
        startt = time()
        loglhis, accepthis, winsizehis = segementorder!(snporder, proptype,
            chrfhaplo,chroffgeno, popmakeup,priorprocess;
            offspringexcl, slidewinsize=nothing, maxwinsize=maxwinsize[p], chrneighbor,
            epsf, epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
            israndallele, issnpGT,decodetempfile, temperature)
        accept = isempty(winsizehis) ? 0.0 : mean(accepthis)
        meanwinsize = isempty(winsizehis) ? 0.0 : round(mean(winsizehis),digits=1)
        # deltlogl = loglhis[end]-loglhis[1]
        tused = round(time()-startt,digits=1)
        temp = [round(accept,digits=4), round(loglhis[end],digits=4)]
        temp = vcat([string("#", proptype, "_neighbor"), meanwinsize,
            length(accepthis),tused], temp)
        push!(logbook,temp)
    end
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end
    logbook
end

# snporder,priorprocess are modified
function segementorder!(snporder::AbstractVector, proptype::AbstractString,
    chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;    
    offspringexcl::AbstractVector, # less efficient; offpsring are excluded on in last step of calculating log likelihood
    slidewinsize::Union{Nothing,Integer,Distribution}=nothing,    
    maxwinsize::Integer,
    chrneighbor::Union{Nothing, AbstractDict}=nothing,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},    
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    israndallele::Bool,     
    issnpGT::AbstractVector, 
    decodetempfile::AbstractString,
    temperature::Real)
    if isnothing(slidewinsize) && isnothing(chrneighbor)
        error("either slidewinsize or chrneighbor must be provided")
    end        
    nsnp = size(chroffgeno,1)
    epsfls = typeof(epsf) <: Real ? epsf*ones(nsnp) : epsf
    epsols = typeof(epso) <: Real ? epso*ones(nsnp) : epso
    seqerrorls = typeof(seqerror) <: Real ? seqerror*ones(nsnp) : seqerror
    allelebalancemeanls = typeof(allelebalancemean) <: Real ? allelebalancemean*ones(nsnp) : allelebalancemean
    allelebalancedispersels = typeof(allelebalancedisperse) <: Real ? allelebalancedisperse*ones(nsnp) : allelebalancedisperse
    alleledropoutls = typeof(alleledropout) <: Real ? alleledropout*ones(nsnp) : alleledropout
    MagicReconstruct.callogbackward_permarker!(decodetempfile, chrfhaplo,chroffgeno, popmakeup,priorprocess;
        epsf = epsfls, epso = epsols, epso_perind, seqerror = seqerrorls,
        allelebalancemean = allelebalancemeanls, allelebalancedisperse = allelebalancedispersels, 
        alleledropout = alleledropoutls, israndallele,issnpGT, snporder) # results are saved in decodetempfile
    tseq = findall(first(values(priorprocess)).markerincl)    
    dataprobls = MagicReconstruct.init_dataprobls_singlephase(popmakeup)    
    jldopen(decodetempfile,"r") do bwfile
        bwkk = 1
        fwkk = 1
        snp = snporder[tseq[1]]
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls,chrfhaplo[snp,:],
            chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind, 
            seqerror=seqerrorls[snp], allelebalancemean=allelebalancemeanls[snp], 
            allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp],  
            israndallele, issiteGT = issnpGT[snp])
        fw_prob_logl = MagicReconstruct.calinitforward(dataprobls, popmakeup)
        fwls = Vector{Union{Nothing,typeof(fw_prob_logl)}}(nothing,maximum(tseq))
        fwls[tseq[1]] = fw_prob_logl
        logbwprob = bwfile[string("t", tseq[1])]
        bwls = Vector{Union{Nothing,typeof(logbwprob)}}(nothing,maximum(tseq))
        bwls[tseq[1]] = logbwprob        
        logl = sum(calindlogl(fw_prob_logl,logbwprob, popmakeup; offspringexcl))
        loglhis = Vector{Float64}()
        push!(loglhis,logl)
        accepthis = Vector{Bool}()
        winsizehis = Vector{Int}()
        kk = 1        
        kkstep = 1.0+max(0.0,min(temperature,2.0)*log2(nsnp/2000)*0.5)
        stepdiv = trunc(Int, kkstep)
        steprem = kkstep - stepdiv
        while kk < length(tseq)                        
            kkseg = get_kkseg(kk, tseq; slidewinsize,chrneighbor,snporder,maxwinsize)
            if isnothing(kkseg)
                kk += 1
                continue
            end
            kkmin, kkmax = kkseg
            if kkmin >= kkmax || kkmax-kkmin+1 > maxwinsize
                kk += 1
                continue
            end 
            if in(proptype,["inverse11","inverse", "permute"])
                kkmid = nothing
            elseif in(proptype,["inverse00", "inverse01","inverse10"])                                
                if kkmax - kkmin > 1                                     
                    kkmid = rand(kkmin:kkmax-1)
                else
                    kk += 1
                    continue
                end
            else
                @error string("unknown proptype=",proptype)
            end
            if fwkk < kkmin-1
                segtseq, segttranls = getsegprop("original", fwkk+1,kkmin-1, tseq)
                fwseg= calsegforward!(dataprobls, segtseq, segttranls,fwkk+1,tseq, fwls,
                    chrfhaplo, chroffgeno, snporder,popmakeup, priorprocess;
                    epsfls,epsols,epso_perind, seqerrorls,allelebalancemeanls,
                    allelebalancedispersels,alleledropoutls, israndallele, issnpGT)
                fwls[tseq[(fwkk+1):(kkmin-1)]] .= fwseg
                fwkk = kkmin-1
            end
            if bwkk > kkmax                
                bw_t = tseq[bwkk]
                isnothing(bwls[bw_t]) && (bwls[bw_t] = bwfile[string("t", bw_t)])
                updatebwls!(dataprobls, bwls, kkmax, bwkk, tseq, chrfhaplo, chroffgeno,
                    snporder,popmakeup, priorprocess;
                    epsfls,epsols,epso_perind, seqerrorls,allelebalancemeanls,
                    allelebalancedispersels,alleledropoutls, israndallele, issnpGT)
                bwkk = kkmax
            end
            bw_tmax = bwls[tseq[kkmax]]
            logbwprob = isnothing(bw_tmax) ? bwfile[string("t", tseq[kkmax])] : bw_tmax
            if rand() < 1.01                
                segtseq, segttranls = getsegprop("original", kkmin, kkmax, tseq)
                fwsegcheck = calsegforward!(dataprobls, segtseq, segttranls,kkmin,tseq,fwls,
                    chrfhaplo, chroffgeno, snporder,popmakeup, priorprocess;
                    epsfls,epsols,epso_perind, seqerrorls,allelebalancemeanls,
                    allelebalancedispersels,alleledropoutls, israndallele, issnpGT)
                logl2 = sum(calindlogl(last(fwsegcheck),logbwprob, popmakeup;offspringexcl))
                msg = string("inconsisent logl=",logl,",logl2=",logl2, "; action=",proptype, "; [kkmin,kkmax]=",[kkmin,kkmax])                 
                isapprox(logl2,logl;atol=1e-3) || @warn msg maxlog=20                
            end
            # if proptype == "original", segttranls[i] == segtseq[i]
            segtseq, segttranls = getsegprop(proptype, kkmin, kkmax, tseq; kkmid)
            fwsegprop = calsegforward!(dataprobls, segtseq, segttranls,kkmin,tseq,fwls,
                chrfhaplo, chroffgeno, snporder,popmakeup, priorprocess;
                epsfls,epsols,epso_perind, seqerrorls,allelebalancemeanls,
                allelebalancedispersels,alleledropoutls,israndallele, issnpGT)
            proplogl = sum(calindlogl(last(fwsegprop),logbwprob, popmakeup;offspringexcl))
            if temperature <= 0.0
                isaccept = proplogl >= logl
            else
                isaccept = rand() < exp((proplogl - logl)/temperature)
            end
            if isaccept
                logl = proplogl
                set_order_prior!(snporder, priorprocess, tseq, kkmin, kkmax, segtseq, segttranls)
                fwls[tseq[kkmin:kkmax]] .= fwsegprop
                if fwkk>kkmax
                    fwls[tseq[(kkmax+1):fwkk]] .= nothing
                end
                fwkk = kkmax
                if bwkk < kkmax
                    bwls[tseq[min(kkmin,bwkk):kkmax]] .= nothing
                end
                bwkk = min(kkmax+1,length(tseq))
            end
            push!(accepthis,isaccept)
            push!(loglhis,logl)
            push!(winsizehis, kkmax-kkmin+1)
            kk += rand() < 1-steprem ? stepdiv : stepdiv + 1            
        end
        loglhis, accepthis, winsizehis
    end
end

function getsegprop(proptype::AbstractString, kkmin::Integer, kkmax::Integer,
    tseq::AbstractVector; kkmid::Union{Nothing, Integer}=nothing)
    if proptype == "original"
        segtls = tseq[kkmin:kkmax-1]
        segtseq = tseq[kkmin:kkmax]
    elseif in(proptype, ["inverse11","inverse"])  # "inverse" === "inverse11"
        segtls = tseq[kkmin:kkmax-1]
        segtseq = tseq[kkmin:kkmax]
        reverse!(segtls)
        reverse!(segtseq)
    elseif proptype == "permute"          
        segtseq = tseq[kkmin:kkmax]
        segtseq .= sample(segtseq, length(segtseq); replace=false)        
        segtls = copy(segtseq)
        setdiff!(segtls, tseq[[kkmax]])
        # segtls = tseq[kkmin:kkmax-1] # keep the original interdistances        
    elseif proptype == "inverse00"         
        segtseq = vcat(tseq[(kkmid+1):kkmax], tseq[kkmin:kkmid])
        segtls = vcat(tseq[(kkmid+1):(kkmax-1)], tseq[kkmin:kkmid])
    elseif proptype == "inverse01"        
        segtseq = vcat(tseq[(kkmid+1):kkmax], reverse(tseq[kkmin:kkmid]))
        segtls = vcat(tseq[(kkmid+1):(kkmax-1)], reverse(tseq[kkmin:kkmid]))
    elseif proptype == "inverse10"        
        segtseq = vcat(reverse(tseq[(kkmid+1):kkmax]), tseq[kkmin:kkmid])
        segtls = vcat(reverse(tseq[(kkmid+1):(kkmax-1)]), tseq[kkmin:kkmid])
    end
    segtseq, segtls
end

function set_order_prior!(snporder::AbstractVector, priorprocess::AbstractDict,
    tseq::AbstractVector,kkmin::Integer, kkmax::Integer, 
    segtseq::AbstractVector, segttranls::AbstractVector)    
    snporder[tseq[kkmin:kkmax]] .= snporder[segtseq]
    for (strkey, pri) in priorprocess
        oldtls = tseq[kkmin:kkmax-1]
        pri.markerdeltd[oldtls] .= pri.markerdeltd[segttranls]
        pri.tranprobseq[oldtls] .= pri.tranprobseq[segttranls]
        pri.markerid[tseq[kkmin:kkmax]] .= pri.markerid[segtseq]
    end
    snporder, priorprocess
end


function get_kkseg(kk::Integer, tseq::AbstractVector;
    slidewinsize::Union{Nothing, Integer,Distribution} = nothing,
    chrneighbor::Union{Nothing, AbstractDict}=nothing,
    snporder::AbstractVector, 
    maxwinsize::Integer)
    kkmin = kk
    if isnothing(slidewinsize)            
        neighbors = chrneighbor[snporder[kkmin]]
        isempty(neighbors)  && return nothing
        neighbors = shuffle(neighbors)        
        rule = Dict(snporder[tseq] .=> 1:length(tseq))
        kknbr = -1
        for nbr in neighbors            
            kknbr = get(rule, nbr,-1)
            kknbr == -1 && continue
            kkdiff = abs(kknbr - kkmin) 
            (kkdiff <=1  || kkdiff > maxwinsize) && continue
        end
        kknbr == -1 && return nothing
        kkdiff = abs(kknbr - kkmin)
        (kkdiff <=1  || kkdiff > maxwinsize) && return nothing
        kkmax = kknbr > kkmin ? kknbr -1 : kknbr + 1
        if kkmin > kkmax
            kkmin, kkmax = kkmax, kkmin
        end 
        [kkmin, kkmax]            
    else            
        winsize = isa(slidewinsize,Integer) ? slidewinsize : rand(slidewinsize)
        if winsize > maxwinsize 
            @warn string("unexpected winsize = ",winsize, ", larger than maxwinsize=", maxwinsize)
            return nothing
        end
        kkmax = min(kkmin+winsize-1,length(tseq))
        [kkmin, kkmax]
    end
end

function calsegforward!(dataprobls, 
    segtseq::AbstractVector, 
    segttranls::AbstractVector, 
    kkmin::Integer, 
    tseq::AbstractVector,    
    fwls::AbstractVector,
    chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    snporder::AbstractVector, popmakeup::AbstractDict,priorprocess::AbstractDict;    
    epsfls::AbstractVector,
    epsols::AbstractVector,
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerrorls::AbstractVector,
    allelebalancemeanls::AbstractVector,
    allelebalancedispersels::AbstractVector,
    alleledropoutls::AbstractVector,
    israndallele::Bool, 
    issnpGT::AbstractVector)
    fwseg = Vector{eltype(fwls)}()    
    # cal first fw    
    snp = snporder[segtseq[1]]
    MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:],chroffgeno[snp,:],popmakeup;
        epsf=epsfls[snp], epso=epsols[snp],seqerror=seqerrorls[snp], epso_perind, 
        allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
        alleledropout=alleledropoutls[snp], israndallele, issiteGT = issnpGT[snp])
    if kkmin == 1
        fw_prob_logl = MagicReconstruct.calinitforward(dataprobls, popmakeup)
    else        
        fw_prob_logl = deepcopy(fwls[tseq[kkmin-1]])
        MagicReconstruct.calnextforward!(fw_prob_logl,tseq[kkmin-1],dataprobls,popmakeup, priorprocess; ttran = tseq[kkmin-1])
    end
    fwseg = Vector{typeof(fw_prob_logl)}()
    push!(fwseg, fw_prob_logl)
    # cal next fw
    for i in 1:length(segtseq)-1
        snp = snporder[segtseq[i+1]]
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:],
            chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
            allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp], israndallele, issiteGT = issnpGT[snp])
        fw_prob_logl2 =  deepcopy(last(fwseg))
        MagicReconstruct.calnextforward!(fw_prob_logl2,segtseq[i],dataprobls,popmakeup, priorprocess;ttran = segttranls[i])
        push!(fwseg,fw_prob_logl2)
    end
    fwseg
end

# bwfile::JLD2.JLDFile,
function updatebwls!(dataprobls, bwls::AbstractVector,
    kkstart::Integer, kkend::Integer,
    tseq::AbstractVector,    
    chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    snporder::AbstractVector, popmakeup::AbstractDict,priorprocess::AbstractDict;
    epsfls::AbstractVector,
    epsols::AbstractVector,
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerrorls::AbstractVector,
    allelebalancemeanls::AbstractVector,
    allelebalancedispersels::AbstractVector,    
    alleledropoutls::AbstractVector,    
    israndallele::Bool, 
    issnpGT::AbstractVector)
    segtseq = tseq[kkstart:kkend]    
    for i = length(segtseq)-1:-1:1
        tback = segtseq[i+1]
        tnow = segtseq[i]
        snp = snporder[tback]
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
            allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp], israndallele, issiteGT = issnpGT[snp])
        logbwprob = copy.(bwls[tback])
        MagicReconstruct.calnextlogbackward!(logbwprob,tback, tnow, dataprobls,popmakeup,priorprocess)
        bwls[tnow] = logbwprob
    end
    bwls
end
