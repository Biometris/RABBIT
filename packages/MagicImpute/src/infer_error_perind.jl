
function infer_error_perind!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;    
    epsf::Union{Real,AbstractVector}, 
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector}, 
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},        
    prior_peroffspringerror::Distribution, 
    issnpGT::AbstractVector,
    isoffphased::Bool=false,
    israndallele::Bool,  
    snporder::Union{Nothing,AbstractVector}=nothing,    
    itmax::Integer=20,
    temperature::Real=0.0)
    isnothing(epso_perind) && return epso_perind
    fderive = MagicReconstruct.precompute_fderive(chrfhaplo,popmakeup)
    offcode = MagicReconstruct.precompute_offcode(chroffgeno,popmakeup,issnpGT; isoffphased)
    for popid in keys(popmakeup)
        offls = popmakeup[popid]["offspring"]
        for off in offls
            epso_perind[off] = infer_epso_perind(off, popid, fderive,offcode,popmakeup,priorprocess;
                epsf,epso,epso_ind = epso_perind[off], seqerror,allelebalancemean, allelebalancedisperse, alleledropout,
                prior_peroffspringerror,issnpGT,isoffphased,israndallele, snporder, itmax,temperature
            )
        end
    end
    epso_perind
end


function infer_epso_perind(offindex::Integer,popid::AbstractString,
    fderive::NamedTuple,offcode::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_ind::Real, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},    
    prior_peroffspringerror::Distribution, 
    issnpGT::AbstractVector,
    isoffphased::Bool=false,
    israndallele::Bool,
    snporder::Union{Nothing,AbstractVector}=nothing,         
    itmax::Integer=20,temperature::Real=0.0)
    accuracygoal, precisiongoal = 4, 4        
    function loglfun(x::Real)
        epso_ind2 = 1/(1+exp(-x))  # inverse of logit transformation             
        logpri = logpdf(prior_peroffspringerror,epso_ind2)        
        loglike_perind(offindex, popid, fderive,offcode,popmakeup,priorprocess;
            epsf,epso,epso_ind = epso_ind2, seqerror,allelebalancemean, allelebalancedisperse, alleledropout,
            issnpGT,isoffphased,israndallele, snporder) + logpri 
    end    
    if epso_ind ≈ 0.0
        xstart = log(1e-4)
        lowbound, upbound = log(1e-5), 0.0
    else
        xstart = log(epso_ind/(1-epso_ind))
        errfraction = 1/(1+exp(-xstart))
        fraction_lowbound = max(1e-5,errfraction/20)
        fraction_upbound = max(0.01,min(0.99, 20*errfraction))        
        lowbound, upbound = log(fraction_lowbound), log(fraction_upbound/(1-fraction_upbound))            
        xstart = min(max(lowbound,xstart),upbound)
    end
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
    1.0/(1.0+exp(-x))
end


function loglike_perind(offindex::Integer,popid::AbstractString,
    fderive::NamedTuple,offcode::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_ind::Real, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},    
    issnpGT::AbstractVector,
    isoffphased::Bool=false,
    israndallele::Bool,
    snporder::Union{Nothing,AbstractVector}=nothing)
    ishaploidsub = popmakeup[popid]["ishaploid"]
    in(offindex, popmakeup[popid]["offspring"]) && "inconsisent popid and offspring"
    nzstate = popmakeup[popid]["nzstate"]             
    hashcode = popmakeup[popid]["hashcode"]    
    markerincl, initprob, tranprobseq = MagicReconstruct.hashcode2prior(priorprocess,hashcode);
    obsseq = view(offcode, :,offindex);
    dataprobseq = [zeros(MagicReconstruct._float_like,length(nzstate))  for _ in 1:length(obsseq)]
    MagicReconstruct.caldataprobseq!(dataprobseq,obsseq,epsf,epso,epso_ind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout,
        fderive,nzstate,isoffphased,israndallele,issnpGT,ishaploidsub);
    dataprobseq2 = view(dataprobseq, snporder[markerincl])               
    fwscale = last(HMM.forward(initprob,tranprobseq,dataprobseq2))
    loglike = HMM.calloglike(fwscale)
    loglike
end

##########################################################################################################


function infer_error_perind_viterbi!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;    
    epsf::Union{Real,AbstractVector}, 
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector}, 
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},        
    prior_peroffspringerror::Distribution, 
    issnpGT::AbstractVector,
    isoffphased::Bool=false,
    israndallele::Bool,  
    snporder::Union{Nothing,AbstractVector}=nothing,    
    decodetempfile::AbstractString,    
    itmax::Integer=20,
    temperature::Real=0.0)
    isnothing(epso_perind) && return epso_perind            
    nsnp = size(chroffgeno,1)
    nsnp == length(issnpGT) || @error string("inconsistent #markers")
    isnothing(snporder) && (snporder = 1:nsnp)    
    nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    
    nfgl == size(chrfhaplo,2) || @error string("inconsistent nfgl")

    # infer hidden origin states; results are saved in decodetempfile        
    MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf,epso,epso_perind, seqerror,allelebalancemean, allelebalancedisperse, alleledropout,
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

    fderive, offcode  = MagicReconstruct.precompute_chr(chrfhaplo, chroffgeno, popmakeup,isoffphased, issnpGT);
    pri1=first(values(priorprocess))
    snpincl = snporder[pri1.markerincl]		    
    for popid in keys(popmakeup)
        offls = popmakeup[popid]["offspring"]
        for off in offls
            epso_perind[off] = infer_epso_perind_viterbi(chrviterbi, fderive,offcode, off,popid,popmakeup,snpincl; 
                epsf,epso,epso_ind = epso_perind[off], seqerror,allelebalancemean, allelebalancedisperse, alleledropout,                
                prior_peroffspringerror,isoffphased,israndallele, issnpGT, itmax,temperature
            ) 
        end
    end
    epso_perind
end


function infer_epso_perind_viterbi(chrviterbi::AbstractMatrix, 
    fderive::NamedTuple,offcode::AbstractMatrix,
    offindex::Integer,popid::AbstractString,    
    popmakeup::AbstractDict,snpincl::AbstractVector;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_ind::Real, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},    
    prior_peroffspringerror::Distribution, 
    issnpGT::AbstractVector,
    isoffphased::Bool=false,
    israndallele::Bool,        
    itmax::Integer=20,temperature::Real=0.0)
    accuracygoal, precisiongoal = 4, 4        
    function loglfun(x::Real)
        epso_ind2 = 1/(1+exp(-x))  # inverse of logit transformation     
        logpri = logpdf(prior_peroffspringerror,epso_ind2)        
        logllike_perind_viterbi(chrviterbi, fderive,offcode,offindex, popid, popmakeup,snpincl;
            epsf,epso,epso_ind = epso_ind2, seqerror,allelebalancemean, allelebalancedisperse, alleledropout,
            issnpGT,isoffphased,israndallele) + logpri
    end    
    xstart = log(epso_ind/(1-epso_ind))
    errfraction = 1/(1+exp(-xstart))
    fraction_lowbound = max(1e-5,errfraction/20)
    fraction_upbound = max(0.01,min(0.99, 20*errfraction))        
    lowbound, upbound = log(fraction_lowbound), log(fraction_upbound/(1-fraction_upbound))        
    xstart = min(max(lowbound,xstart),upbound)
    if temperature ≈ 0.0
        res= MagicBase.brentMax(loglfun,lowbound,upbound;
            xstart, precisiongoal,accuracygoal,maxiter=itmax)            
        x = res[1]
    else
        constraint = x-> lowbound < x < upbound
        res2 = MagicBase.metroplis1d(loglfun; xstart, constraint, temperature,
            stepsize = log(5.0), nstep = 5)
        x = res2[end][1]        
    end
    1.0/(1.0+exp(-x))
end


function logllike_perind_viterbi(chrviterbi::AbstractMatrix, 
    fderive::NamedTuple,
    offcode::AbstractMatrix,
    offindex::Integer,
    popid::AbstractString,
    popmakeup::AbstractDict,
    snpincl::AbstractVector;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_ind::Real,     
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    isoffphased::Bool,
    israndallele::Bool,
    issnpGT::AbstractVector)        
    ishaploidsub = popmakeup[popid]["ishaploid"]
    in(offindex, popmakeup[popid]["offspring"]) && "inconsisent popid and offspring"
    nzstate = popmakeup[popid]["nzstate"]             
    obsseq = view(offcode, :,offindex);
    dataprobseq = [zeros(MagicReconstruct._float_like,length(nzstate))  for _ in 1:length(obsseq)]
    MagicReconstruct.caldataprobseq!(dataprobseq,obsseq,epsf,epso,epso_ind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout,
        fderive,nzstate,isoffphased,israndallele,issnpGT,ishaploidsub);
    dataprobls = view(dataprobseq, snpincl)
    originls = chrviterbi[snpincl, offindex]
    sum(log.(map((x,i)->x[i],dataprobls,originls)))    
end
