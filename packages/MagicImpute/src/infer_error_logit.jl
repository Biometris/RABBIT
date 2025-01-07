
function infer_errorls!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    targetls::AbstractVector,
    priorlikeparameters::PriorLikeParameters,
    epsf::Union{Real,AbstractVector}, 
    epso::Union{Real,AbstractVector},    
    epso_perind::Union{Nothing,AbstractVector},     
    seqerror::Union{Real,AbstractVector}, 
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},        
    offspringexcl::AbstractVector, 
    issnpGT::AbstractVector,
    snporder::Union{Nothing,AbstractVector}=nothing,
    decodetempfile::AbstractString,    
    itmax::Integer=20,temperature::Real=0.0)
    nsnp = size(chroffgeno,1)
    nsnp == length(issnpGT) || @error string("inconsistent #markers")
    isnothing(snporder) && (snporder = 1:nsnp)    

    epsfls = typeof(epsf) <: Real ? epsf*ones(nsnp) : epsf
    epsols = typeof(epso) <: Real ? epso*ones(nsnp) : epso
    seqerrorls = typeof(seqerror) <: Real ? seqerror*ones(nsnp) : seqerror
    allelebalancemeanls = typeof(allelebalancemean) <: Real ? allelebalancemean*ones(nsnp) : allelebalancemean
    allelebalancedispersels = typeof(allelebalancedisperse) <: Real ? allelebalancedisperse*ones(nsnp) : allelebalancedisperse  
    alleledropoutls = typeof(alleledropout) <: Real ? alleledropout*ones(nsnp) : alleledropout  
    dataprobls = init_dataprobls_singlephase(popmakeup)    
    for target in targetls
        tseq = findall(first(values(priorprocess)).markerincl)
        if in(target,["seqerror","allelebalancemean","allelebalancedisperse"])
            kmax = findlast(.!issnpGT[snporder[tseq]])
            kmax < length(tseq) && deleteat!(tseq,kmax+1:length(tseq))
        end
        if target == "offspringerror"
            typeof(epso) <: Real && @error string("epso must be a vector")                      
            prior_err = priorlikeparameters.offspringerror 
        elseif target == "foundererror"
            typeof(epsf) <: Real && @error string("espf must be a vector")                
            prior_err = priorlikeparameters.foundererror 
        elseif target == "seqerror"
            typeof(seqerror) <: Real && @error string("seqerror must be a vector")          
            prior_err = priorlikeparameters.seqerror  
        elseif target == "allelebalancemean"
            typeof(allelebalancemean) <: Real && @error string("allelebalancemean must be a vector")        
            prior_err = priorlikeparameters.allelebalancemean  
        elseif target == "allelebalancedisperse"
            typeof(allelebalancedisperse) <: Real && @error string("allelebalancedisperse must be a vector")        
            prior_err = priorlikeparameters.allelebalancedisperse 
        elseif target == "alleledropout"
            typeof(alleledropout) <: Real && @error string("alleledropout must be a vector")                
            prior_err = priorlikeparameters.alleledropout  
        else
            @error string("unknown genoerror type: ",target)
        end        
        # @info string("target=",target, ",prior_err=",prior_err)
        MagicReconstruct.callogbackward_permarker!(decodetempfile, chrfhaplo,chroffgeno, popmakeup,priorprocess;
            epsf, epso,epso_perind, seqerror,allelebalancemean,allelebalancedisperse, alleledropout, 
            issnpGT, snporder) # results are saved in decodetempfile                
        jldopen(decodetempfile,"r") do bwfile
            snp = snporder[tseq[1]]
            MagicReconstruct.calsitedataprob_singlephase!(dataprobls,chrfhaplo[snp,:],
                chroffgeno[snp,:],popmakeup;
                epsf=epsfls[snp], epso = epsols[snp], epso_perind, 
                seqerror = seqerrorls[snp],allelebalancemean=allelebalancemeanls[snp], 
                allelebalancedisperse=allelebalancedispersels[snp],
                alleledropout=alleledropoutls[snp],            
                issiteGT=issnpGT[snp])
            fw_prob_logl = MagicReconstruct.calinitforward(dataprobls, popmakeup)
            # logbwprob = bwfile[string("t", tseq[1])]
            # logl = sum(calindlogl(fw_prob_logl,logbwprob, popmakeup))
            for kk in eachindex(tseq)
                if kk >= 3
                    snp = snporder[tseq[kk-1]]
                    MagicReconstruct.calsitedataprob_singlephase!(dataprobls,chrfhaplo[snp,:],
                        chroffgeno[snp,:],popmakeup;
                        epsf=epsfls[snp], epso=epsols[snp], 
                        seqerror = seqerrorls[snp],allelebalancemean=allelebalancemeanls[snp], 
                        allelebalancedisperse=allelebalancedispersels[snp],
                        alleledropout=alleledropoutls[snp],                    
                        issiteGT = issnpGT[snp])
                    # input fw_prob_logl refers to tseq[kk-2], output fw_prob_logl refers to tseq[kk-1]
                    MagicReconstruct.calnextforward!(fw_prob_logl,tseq[kk-2], dataprobls,popmakeup, priorprocess)
                end
                tpre = kk==1 ? 0 : tseq[kk-1]
                tnow = tseq[kk]
                snp = snporder[tnow]
                issnpGT[snp] && in(target,["seqerror","allelebalancemean","allelebalancedisperse","alleledropout"]) && continue
                logbw_tnow = bwfile[string("t", tnow)]
                errmarker = infer_errmarker!(dataprobls,offspringexcl, tpre, tnow, target, prior_err,
                    epsfls[snp],epsols[snp], epso_perind, seqerrorls[snp],
                    allelebalancemeanls[snp],allelebalancedispersels[snp], alleledropoutls[snp],
                    fw_prob_logl[1], logbw_tnow, 
                    chrfhaplo, chroffgeno, popmakeup, priorprocess;
                    issnpGT, snporder, itmax,temperature)     
                # println("target=", target, ",kk=",kk, ",errmarker=",errmarker)       
                if target == "offspringerror"
                    epsols[snp] = errmarker                
                elseif target == "foundererror"
                    epsfls[snp] = errmarker                
                elseif target == "seqerror"
                    seqerrorls[snp] = errmarker                
                elseif target == "allelebalancemean"
                    allelebalancemeanls[snp] = errmarker    
                elseif target == "allelebalancedisperse"
                    allelebalancedispersels[snp] = errmarker                
                elseif target == "alleledropout"
                    alleledropoutls[snp] = errmarker                
                else
                    @error string("unknown genoerror type: ",target)
                end
            end
        end
    end
    return nothing
end

function fisherz(rho::Real)
    log((1+rho)/(1-rho))/2
end

function infer_errmarker!(dataprobls, offspringexcl::AbstractVector, 
    tpre::Integer,tnow::Integer,
    target::AbstractString,
    prior_err::Union{Nothing,Distribution},
    epsfmarker::Real, 
    epsomarker::Real,
    epso_perind::Union{Nothing,AbstractVector},     
    seqerrormarker::Real, 
    allelebalancemeanmarker::Real,    
    allelebalancedispersemarker::Real,    
    alleledropoutmarker::Real,    
    fwprob::AbstractVector,
    logbw_tnow::AbstractVector,
    chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict, priorprocess::AbstractDict;
    issnpGT::AbstractVector,
    snporder::AbstractVector,    
    itmax::Integer=20,temperature::Real=0.0)
    accuracygoal, precisiongoal = 4, 4
    fwpriorprob = tofwpriorprob(tpre, fwprob,popmakeup, priorprocess)
    bwprobls = [exp.(i .- maximum(i)) for i in logbw_tnow]
    snp = snporder[tnow]
    function loglfun(x::Real)
        epsf2, epso2, seqerror2,allelebalancemean2,allelebalancedisperse2,alleledropout2 = epsfmarker, epsomarker, seqerrormarker,allelebalancemeanmarker,allelebalancedispersemarker,alleledropoutmarker
        if target == "offspringerror"
            epso2 = 1/(1+exp(-x))  # inverse of logit transformation     
        elseif target == "foundererror"
            epsf2 = 1/(1+exp(-x))                      
        elseif target == "seqerror"
            seqerror2 = 1/(1+exp(-x))                     
        elseif target == "allelebalancemean"
            allelebalancemean2 = 1/(1+exp(-x))    
        elseif target == "allelebalancedisperse"
            allelebalancedisperse2 = exp(x)
        elseif target == "alleledropout"
            alleledropout2 = 1/(1+exp(-x))                     
        else
            @error string("unknown genoerror type: ",target)
        end
        if target in ["foundererror", "offspringerror", "seqerror","allelebalancemean","alleledropout"]
            p =  1/(1+exp(-x)) 
            logpri = logpdf(prior_err,p)            
            if target in ["alleledropout"]
                # estimation refers to original parameter space
                logjacobi = 0.0                   
            else
                # estimation refers to transformed space
                logjacobi = x - 2*log(1+exp(x)) # dp/dx = exp(x)/(1+exp(x))^2              
            end
        elseif target == "allelebalancedisperse"    
            p = exp(x)
            logpri = logpdf(prior_err,p)                 
            logjacobi = x   # dp/dx = exp(x) # estimation refers to transformed space       
        end    
        calloglmarker!(dataprobls,offspringexcl, epsf2, epso2, epso_perind, seqerror2,allelebalancemean2,allelebalancedisperse2, alleledropout2, snp, 
            fwpriorprob, bwprobls, chrfhaplo, chroffgeno,popmakeup, issnpGT) + logpri + logjacobi
    end    
    if target == "offspringerror"
        xstart = log(epsomarker/(1-epsomarker))
    elseif target == "foundererror"
        xstart = log(epsfmarker/(1-epsfmarker))
    elseif target == "seqerror"
        xstart = log(seqerrormarker/(1-seqerrormarker))
    elseif target == "allelebalancemean"
        xstart = log(allelebalancemeanmarker/(1-allelebalancemeanmarker))
    elseif target == "allelebalancedisperse"
        xstart = max(-10.0,log(allelebalancedispersemarker))
    elseif target == "alleledropout"
        xstart = log(alleledropoutmarker/(1-alleledropoutmarker))
    else
        @error string("unknown genoerror type: ",target)
    end       
    if target == "allelebalancedisperse"
        lowbound, upbound = max(-12.0,xstart-10.0), xstart+10.0
    else
        errfraction = 1/(1+exp(-xstart))
        fraction_lowbound = max(1e-5,errfraction/20)
        fraction_upbound = max(0.01,min(0.99, 20*errfraction))        
        lowbound, upbound = log(fraction_lowbound), log(fraction_upbound/(1-fraction_upbound))                
    end    
    xstart = min(max(lowbound,xstart),upbound)
    if temperature ≈ 0
        # always using brent
        res= MagicBase.brentMax(loglfun,lowbound,upbound;
            xstart, precisiongoal,accuracygoal,maxiter=itmax)            
        x = res[1]
    else
        constraint = x-> lowbound < x < upbound
        res2 = MagicBase.metroplis1d(loglfun; xstart, constraint, temperature,
            stepsize = log(5.0), nstep = 5)
        x = res2[end][1]        
    end
    est =  target == "allelebalancedisperse" ? exp(x) : 1.0/(1.0+exp(-x)) 
    est
end

function calloglmarker!(dataprobls::AbstractVector,offspringexcl::AbstractVector, 
    epsf::Real,epso::Real, epso_perind::Union{Nothing,AbstractVector},     
    seqerror::Real,
    allelebalancemean::Real,allelebalancedisperse::Real,alleledropout::Real,
    snp::Integer,
    fwpriorprob::AbstractVector,
    bwprobls::AbstractVector,
    chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,
    issnpGT::AbstractVector)
    calsitedataprob_singlephase!(dataprobls,chrfhaplo[snp,:], chroffgeno[snp,:],popmakeup;
        epsf, epso, epso_perind, seqerror,allelebalancemean, allelebalancedisperse,alleledropout, issiteGT = issnpGT[snp])        
    noff = length(dataprobls) 
    offincl = isempty(offspringexcl) ? (1:noff) : setdiff(1:noff, offspringexcl)                
    sum([log(sum(fwpriorprob[off] .* dataprobls[off] .*bwprobls[off])) for off in offincl])    
end

function tofwpriorprob(tpre::Integer,fwprob::AbstractMatrix,
    popmakeup::AbstractDict, priorprocess::AbstractDict)
    fwpriorprob = Vector{Vector{Float64}}(undef, size(fwprob,1))
    if tpre == 0
        for popid in keys(popmakeup)
            initprob = copy(popmakeup[popid]["initprob"])
            offls = popmakeup[popid]["offspring"]
            for off in offls
                fwpriorprob[off] = initprob
            end
        end
    else
        for popid in keys(popmakeup)
            hashcode = popmakeup[popid]["hashcode"]
            tranprob = priorprocess[hashcode].tranprobseq[tpre]            
            offls = popmakeup[popid]["offspring"]
            for off = offls                
                if isnothing(tranprob)
                    fwpriorprob[off] = copy(fwprob[off])
                else
                    fwpriorprob[off] = (tranprob' * fwprob[off])
                end
            end
        end
    end
    fwpriorprob
end
