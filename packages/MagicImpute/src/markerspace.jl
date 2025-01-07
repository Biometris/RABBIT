
function markerspace_chr!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    offspringexcl::AbstractVector, 
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},    
    epso_perind::Union{Nothing,AbstractVector},     
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    snporder::AbstractVector,    
    israndallele::Bool, 
    issnpGT::AbstractVector, 
    decodetempfile::AbstractString,
    reversechr::Bool=false,
    priorlength::Real=10.0,
    temperature::Real=0.0,itmax::Integer=20)
    # in prior, an intermarker distance follows an exponential distribution
    # with mean priorlength in centimorgn
    nsnp = size(chroffgeno,1)
    isnothing(snporder) && (snporder = collect(1:nsnp))    
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
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
    offlogl = jldopen(decodetempfile,"r") do bwfile
        snp = snporder[tseq[1]]
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
            chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso = epsols[snp], epso_perind, 
            seqerror=seqerrorls[snp], allelebalancemean=allelebalancemeanls[snp], 
            allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp], 
            israndallele,issiteGT = issnpGT[snp])
        fw_prob_logl = MagicReconstruct.calinitforward(dataprobls, popmakeup)
        for kk in 1:(length(tseq)-1)
            if kk >= 2
                snp = snporder[tseq[kk]]
                MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
                    chroffgeno[snp,:],popmakeup;
                    epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
                    allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
                    alleledropout=alleledropoutls[snp], israndallele, issiteGT = issnpGT[snp])
                # input fw_prob_logl refers to tseq[kk-1], output fw_prob_logl refers to tseq[kk]
                MagicReconstruct.calnextforward!(fw_prob_logl,tseq[kk-1],dataprobls,popmakeup, priorprocess)
            end            
            snp = snporder[tseq[kk+1]]
            MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
                chroffgeno[snp,:],popmakeup;
                epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
                allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
                alleledropout=alleledropoutls[snp],  israndallele, issiteGT = issnpGT[snp])
            # fw_prob_logl refers to tseq[kk]
            logbwprob = bwfile[string("t", tseq[kk+1])]
            tnowdis = calinterdis(tseq[kk],popmakeup,priorprocess,offspringexcl, 
                dataprobls,fw_prob_logl[1],logbwprob,priorlength, temperature,itmax)
            setdistanceat!(priorprocess,tseq[kk],tnowdis)
        end
        kk = length(tseq)
        snp = snporder[tseq[kk]]
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
            chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
            allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp], israndallele, issiteGT = issnpGT[snp])
        # input fw_prob_logl refers to tseq[kk-1], output fw_prob_logl refers to tseq[kk]
        MagicReconstruct.calnextforward!(fw_prob_logl,tseq[kk-1],dataprobls,popmakeup, priorprocess)
        logbwprob = bwfile[string("t", tseq[kk])]       
    end
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end
    offlogl
end

function calinterdis(tnow::Integer,
    popmakeup::AbstractDict, priorprocess::AbstractDict,
    offspringexcl::AbstractVector, 
    dataprobls::AbstractVector,
    fwprob::AbstractVector, 
    logbwprob::AbstractVector,
    priorlength::Real,    
    temperature::Real=0.0,itmax::Integer=20)
    # Haldane map function
    pri1=first(values(priorprocess))
    markerdeltd = pri1.markerdeltd[1:end-1]
    accuracygoal, precisiongoal = 5, 5
    lowbound,upbound = log(10^(-7)), log(1.0)    
    # priorchrlength in unit of centiMorgan
    function loglfun(x::Real)        
        d = exp(x)
        callogldis(d, popmakeup,priorprocess,dataprobls,fwprob,logbwprob,offspringexcl)-d*100.0/priorlength
    end
    xstart= log(max(exp(lowbound),markerdeltd[tnow]))
    if temperature ≈ 0
        # always use brent
        res= MagicBase.brentMax(loglfun,lowbound,upbound;
            xstart, precisiongoal,accuracygoal,maxiter=itmax)
    else
        constraint = x-> lowbound < x < upbound
        res = last(MagicBase.metroplis1d(loglfun; xstart, constraint, temperature,
            stepsize = log(5.0), nstep = 5))
    end
    tnowdis = round(exp(res[1]),digits=6)    
    tnowdis
end

function callogldis(d::Real,popmakeup::AbstractDict, priorprocess::AbstractDict,
    dataprobls::AbstractVector,fwprob::AbstractVector,logbwprob::AbstractVector,
    offspringexcl::AbstractVector)
    res = 0.0    
    for popid in keys(popmakeup)
        hashcode = popmakeup[popid]["hashcode"]        
        offls = popmakeup[popid]["offspring"]
        offls2 = isempty(offspringexcl) ? offls : setdiff(offls, offspringexcl)                
        tranprob = MagicReconstruct.get_tranprobmatrix(d,priorprocess[hashcode].tranrate)          
        if isnothing(tranprob)
            for off in offls2
                bwprob = exp.(logbwprob[off] .- maximum(logbwprob[off]))
                s = sum(fwprob[off] .* dataprobls[off] .* bwprob)
                res += log(s)
            end
        else
            for off in offls2
                bwprob = exp.(logbwprob[off] .- maximum(logbwprob[off]))
                s = sum((tranprob' * fwprob[off]) .* dataprobls[off] .* bwprob)
                res += log(s)
            end
        end
    end
    res
end


################################################################


function markerspace_chr_viterbi!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    offspringexcl::AbstractVector, 
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    snporder::AbstractVector,    
    israndallele::Bool, 
    issnpGT::AbstractVector, 
    decodetempfile::AbstractString,
    reversechr::Bool=false,
    priorlength::Real=10.0,
    temperature::Real=0.0,itmax::Integer=20)
    # in prior, an intermarker distance follows an exponential distribution
    # with mean priorlength in centimorgn
    nsnp = size(chroffgeno,1)
    isnothing(snporder) && (snporder = collect(1:nsnp))
    nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    
    nfgl == size(chrfhaplo,2) || @error string("inconsistent nfgl")
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end
    # infer hidden origin states; results are saved in decodetempfile        
    MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf,epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
        hmmalg="viterbi", decodetempfile, israndallele, issnpGT,snporder);
    # chrviterbi: nsnp x noff, row i = viterbi-path of input makrer index i
    chrviterbi = MagicReconstruct.get_chr_viterbi(decodetempfile, snporder)
    for popid in keys(popmakeup)        
        # nzcol is concisistent with the results of hmm_decode
        nzcol = MagicReconstruct.get_nzcol(nstate,popmakeup,popid)        
        offls = popmakeup[popid]["offspring"]
        pairs = nzcol .=> 1:length(nzcol)
        replace!(view(chrviterbi,:,offls),pairs...) 
    end
    viterbils = collect(eachrow(chrviterbi))
    # infer distance    
    tseq = findall(first(values(priorprocess)).markerincl)
    originnow = viterbils[snporder[tseq[1]]]
    for kk in 1:(length(tseq)-1)
        originnext = viterbils[snporder[tseq[kk+1]]]
        originls = CartesianIndex.(originnow,originnext)
        tnowdis = calinterdis(tseq[kk], originls,popmakeup, priorprocess,offspringexcl, priorlength, temperature,itmax)
        MagicImpute.setdistanceat!(priorprocess,tseq[kk],tnowdis)
        originnow = originnext
    end
    if reversechr
        reverse!(snporder)
        for val in values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end
    return nothing
end


function calinterdis(tnow::Integer,originls::Vector{CartesianIndex{2}},
    popmakeup::AbstractDict, priorprocess::AbstractDict,
    offspringexcl::AbstractVector,
    priorlength::Real,temperature::Real=0.0,itmax::Integer=20)
    pri1=first(values(priorprocess))
    markerdeltd = pri1.markerdeltd[1:end-1]
    accuracygoal, precisiongoal = 5, 5
    lowbound,upbound = log(10^(-7)), log(1.0)
    # priorchrlength in unit of centiMorgan
    function loglfun(x::Real)        
        d = round(exp(x), digits=7)        
        callogldis(d, originls,popmakeup,priorprocess,offspringexcl)-d*100.0/priorlength
    end
    xstart= log(max(exp(lowbound),markerdeltd[tnow]))
    if temperature == 0
        res= MagicBase.brentMax(loglfun,lowbound,upbound;
            xstart, precisiongoal,accuracygoal,maxiter=itmax)
    else
        constraint = x-> lowbound < x < upbound
        res = last(MagicBase.metroplis1d(loglfun; xstart, constraint, temperature,
            stepsize = log(5.0), nstep = 5))
    end
    tnowdis = round(exp(res[1]),digits=6)    
    tnowdis
end


function  callogldis(d::Real,originls::Vector{CartesianIndex{2}},popmakeup::AbstractDict,priorprocess::AbstractDict,
    offspringexcl::AbstractVector)
    res = 0.0
    for popid in keys(popmakeup)
        hashcode = popmakeup[popid]["hashcode"]
        logtranprob = log.(exp(d .* priorprocess[hashcode].tranrate))
        offls = popmakeup[popid]["offspring"]
        offls2 = isempty(offspringexcl) ? offls : setdiff(offls, offspringexcl)                        
        res += sum(logtranprob[originls[offls2]])
    end
    res
end


################################################################

function setdistanceat!(priorprocess::AbstractDict,tnow::Integer,tnowdis::Real)
    for (_, pri) in priorprocess
        d = tnowdis
        pri.markerdeltd[tnow] = d
        pri.tranprobseq[tnow] = MagicReconstruct.get_tranprobmatrix(d,pri.tranrate)
    end
    priorprocess
end



# trim left ( or right) end with gad to the rest.
function trimchrend!(priorprocess::AbstractDict; trimcm::Real=20,trimfraction::Real=0.05)
    pri1=first(values(priorprocess))
    tseq = findall(pri1.markerincl);    
    deltd = pri1.markerdeltd[tseq[1:end-1]]        
    segls = trimsegls(deltd, trimcm, trimfraction)
    isempty(segls) && return (0,"")
    dells = [tseq[i] for i in segls]        
    trimdeltd = [begin 
        if first(i) == 1 
            segdeltd = deltd[1:last(i)] 
        elseif last(i) == length(deltd) + 1
            segdeltd = deltd[first(i)-1:end]
        else
            segdeltd = deltd[first(i)-1:last(i)]
        end
        round.(segdeltd,digits=3) 
    end for i in segls]   
    msg_trim = string("seg=", segls,"; deltd=", join(trimdeltd,","))
    deltt = reduce(vcat, dells)
    allunique(deltt) || @error string("overlapping segls=",segls)
    isempty(deltt) || MagicReconstruct.setpriorprocess!(priorprocess,deltt)
    ndel = length(deltt)
    ndel,msg_trim
end


function trimsegls(deltd::AbstractVector, trimcm::Real, trimfraction::Real)
    isgaps = deltd .> 0.01*trimcm # cM to Morgan
    gaps = findall(isgaps)
    segls = []
    isempty(gaps) && return segls
    maxseg = max(2,round(Int, trimfraction*(length(deltd)+1)))
    # leftend    
    leftseg = gaps[gaps .<= maxseg]
    if !isempty(leftseg) 
        pos = argmax(deltd[leftseg])
        push!(segls,1:leftseg[pos])
        setdiff!(gaps,leftseg[1:pos])
    end
    # rightend    
    isempty(gaps) && return segls
    maxindex = length(deltd)+1
    rightseg = gaps[gaps .>= maxindex-maxseg]
    if !isempty(rightseg) 
        pos = argmax(deltd[rightseg])
        push!(segls, rightseg[pos]+1:maxindex)
        setdiff!(gaps,rightseg[pos:end])
    end    
    # midseg    
    length(gaps) >= 2 || return segls
    insertpos = isempty(segls) ? 1 : (segls[1][1]==1 ? 2 : 1) 
    pos = argmax(deltd[gaps])        
    if pos == 1
        if gaps[2] - gaps[1] <= maxseg
            insert!(segls,insertpos, gaps[1]+1:gaps[2])
        end
    elseif pos == length(gaps)
        if gaps[pos] - gaps[pos-1] <= maxseg
            insert!(segls,insertpos, gaps[pos-1]+1:gaps[pos])
        end
    else
        if gaps[pos+1] - gaps[pos] < gaps[pos] - gaps[pos-1] 
            pos += 1
        end
        if gaps[pos] - gaps[pos-1] <= maxseg
            insert!(segls,insertpos, gaps[pos-1]+1:gaps[pos])
        end
    end
    # results    
    segls
end


# trim left ( or right) end with gad to the rest.
# also trim middel seg with gad to both sides.
# function trimchrend!(priorprocess::AbstractDict; trimcm::Real=20,trimfraction::Real=0.05)
#     pri1=first(values(priorprocess))
#     tseq = findall(pri1.markerincl);
#     tmax = last(tseq)
#     deltd = pri1.markerdeltd[tseq[1:end-1]]
#     gaps = tseq[1:end-1][deltd .> 0.01*trimcm] # cM to Morgan
#     isempty(gaps) && return (0,"")
#     maxseg = round(Int, trimfraction*length(tseq))
#     segls = trimsegls(gaps, tmax, maxseg)
#     # println("gaps=",gaps,", tmax=",tmax, ",maxseg=",maxseg, ", segls=",segls)
#     isempty(segls) && return (0,"")
#     dells = [intersect(span(i), tseq) for i in segls]
#     trimdeltd = [begin
#         i2 = tmax in i ? i[i .!= tmax] : i
#         round.(pri1.markerdeltd[i2],digits=3)
#     end for i in dells]
#     msg_trim = string("seg=", join(segls,","),
#         "; deltd=", join(trimdeltd,","))
#     for i in 1:length(segls)
#         # segls[i] == 1: left end
#         # segls[i] ~= 1: right end
#         first(segls[i]) !=1 && popfirst!(dells[i])
#     end
#     deltt = reduce(vcat, dells)
#     isempty(deltt) || MagicReconstruct.setpriorprocess!(priorprocess,deltt)
#     ndel = length(deltt)
#     ndel,msg_trim
# end

# function trimsegls(gaps::AbstractVector, tmax::Integer, maxseg::Integer)
#     gaps2 = vcat(gaps,[tmax])
#     segls = Vector{Vector{Int}}()
#     seg = [1]
#     for gap in gaps2
#         if gap - last(seg) < maxseg
#             push!(seg,gap)
#         else
#             push!(segls, seg)
#             seg=[gap]
#         end
#         # println("seg=", seg, ",segls=",segls)
#     end
#     push!(segls, seg)
#     segls[length.(segls) .> 1]
# end
