
function infer_errorls_viterbi!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
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
    isoffphased::Bool=false, 
    israndallele::Bool, 
    snporder::Union{Nothing,AbstractVector}=nothing,
    decodetempfile::AbstractString,    
    itmax::Integer=20,temperature::Real=0.0)
    nsnp = size(chroffgeno,1)
    nsnp == length(issnpGT) || @error string("inconsistent #markers")
    isnothing(snporder) && (snporder = 1:nsnp)    
    nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    
    nfgl == size(chrfhaplo,2) || @error string("inconsistent nfgl")
     # infer hidden origin states; results are saved in decodetempfile        
     MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf,epso, epso_perind,seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
        hmmalg="viterbi", decodetempfile, issnpGT, isoffphased, israndallele, snporder);
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
    
    # infer error        
    epsfls = typeof(epsf) <: Real ? epsf*ones(nsnp) : epsf
    epsols = typeof(epso) <: Real ? epso*ones(nsnp) : epso
    seqerrorls = typeof(seqerror) <: Real ? seqerror*ones(nsnp) : seqerror
    allelebalancemeanls = typeof(allelebalancemean) <: Real ? allelebalancemean*ones(nsnp) : allelebalancemean
    allelebalancedispersels = typeof(allelebalancedisperse) <: Real ? allelebalancedisperse*ones(nsnp) : allelebalancedisperse  
    alleledropoutls = typeof(alleledropout) <: Real ? alleledropout*ones(nsnp) : alleledropout  
    dataprobls = MagicReconstruct.init_dataprobls_singlephase(popmakeup)
    for target in targetls
        tseq = findall(first(values(priorprocess)).markerincl)
        if in(target,["seqerror","allelebalancemean","allelebalancedisperse"])
            kmax = findlast(.!issnpGT[snporder[tseq]])
            !isnothing(kmax) && kmax < length(tseq) && deleteat!(tseq,kmax+1:length(tseq))
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
        for kk in eachindex(tseq)
            # if target == "seqerror" && kk==1
            #     println("target=",target,",prior_err=",prior_err)
            # end
            tnow = tseq[kk]
            snp = snporder[tnow]
            issnpGT[snp] && in(target,["seqerror","allelebalancemean","allelebalancedisperse","alleledropout"]) && continue
            originls = viterbils[snp]        
            errmarker = infer_errmarker_viterbi!(dataprobls,snp, originls, target, prior_err,
                epsfls[snp],epsols[snp], epso_perind, seqerrorls[snp],
                allelebalancemeanls[snp],allelebalancedispersels[snp], alleledropoutls[snp],                            
                offspringexcl, chrfhaplo, chroffgeno, popmakeup;
                isoffphased, israndallele, issnpGT, itmax,temperature)                  
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
    return nothing
end

function infer_errmarker_viterbi!(dataprobls::AbstractVector, snp::Integer,
    originls::AbstractVector,     
    target::AbstractString,
    prior_err::Distribution,
    epsfmarker::Real, 
    epsomarker::Real,
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerrormarker::Real, 
    allelebalancemeanmarker::Real,    
    allelebalancedispersemarker::Real,    
    alleledropoutmarker::Real,        
    offspringexcl::AbstractVector, 
    chrfhaplo::AbstractMatrix,
    chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict;        
    isoffphased::Bool, 
    israndallele::Bool, 
    issnpGT::AbstractVector,        
    itmax::Integer=20,temperature::Real=0.0)
    accuracygoal, precisiongoal = 4, 4        
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
        calloglmarker_viterbi!(dataprobls, originls, offspringexcl, epsf2, epso2, epso_perind, seqerror2,allelebalancemean2,allelebalancedisperse2, alleledropout2, snp, 
            chrfhaplo, chroffgeno,popmakeup, isoffphased, israndallele, issnpGT) + logpri + logjacobi
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

function calloglmarker_viterbi!(dataprobls::AbstractVector,originls::AbstractVector, 
    offspringexcl::AbstractVector, 
    epsf::Real,epso::Real, 
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Real,
    allelebalancemean::Real,allelebalancedisperse::Real,alleledropout::Real,
    snp::Integer,    
    chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,    
    isoffphased::Bool,
    israndallele::Bool,
    issnpGT::AbstractVector)
    MagicReconstruct.calsitedataprob_singlephase!(dataprobls,chrfhaplo[snp,:],chroffgeno[snp,:],popmakeup;
        epsf, epso, epso_perind, seqerror,allelebalancemean, allelebalancedisperse,
        alleledropout, isoffphased, israndallele, issiteGT = issnpGT[snp])  
    noff = length(dataprobls) 
    offincl = isempty(offspringexcl) ? (1:noff) : setdiff(1:noff, offspringexcl)                
    sum(log.(map((x,i)->x[i],dataprobls[offincl],originls[offincl])))    
end

