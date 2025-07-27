
function infer_errorls_viterbi!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    targetls::AbstractVector,
    priorlikeparam::PriorLikeParam,
    epsf::Union{Real,AbstractVector}, 
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},
    allelicdropout::Union{Real,AbstractVector},        
    avgerrdict::AbstractDict, 
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
        epsf,epso, epso_perind,baseerror,allelicbias,allelicoverdispersion,allelicdropout,
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
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nsnp) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nsnp) : allelicbias
    allelicoverdispersionls = typeof(allelicoverdispersion) <: Real ? allelicoverdispersion*ones(nsnp) : allelicoverdispersion  
    allelicdropoutls = typeof(allelicdropout) <: Real ? allelicdropout*ones(nsnp) : allelicdropout  
    dataprobls = MagicReconstruct.init_dataprobls_singlephase(popmakeup)    
    for target in targetls
        tseq = findall(first(values(priorprocess)).markerincl)
        if in(target,["baseerror","allelicbias","allelicoverdispersion"])
            kmax = findlast(.!issnpGT[snporder[tseq]])
            !isnothing(kmax) && kmax < length(tseq) && deleteat!(tseq,kmax+1:length(tseq))
        end
        if target == "offspringerror"
            typeof(epso) <: Real && @error string("epso must be a vector")                      
            prior_err = priorlikeparam.offspringerror  
        elseif target == "foundererror"
            typeof(epsf) <: Real && @error string("espf must be a vector")                
            prior_err = priorlikeparam.foundererror 
        elseif target == "baseerror"
            typeof(baseerror) <: Real && @error string("baseerror must be a vector")          
            prior_err = priorlikeparam.baseerror  
        elseif target == "allelicbias"
            typeof(allelicbias) <: Real && @error string("allelicbias must be a vector")        
            prior_err = priorlikeparam.allelicbias  
        elseif target == "allelicoverdispersion"
            typeof(allelicoverdispersion) <: Real && @error string("allelicoverdispersion must be a vector")        
            prior_err = priorlikeparam.allelicoverdispersion  
        elseif target == "allelicdropout"
            typeof(allelicdropout) <: Real && @error string("allelicdropout must be a vector")                
            prior_err = priorlikeparam.allelicdropout  
        else
            @error string("unknown genoerror type: ",target)
        end               
        if target in ["foundererror", "offspringerror", "baseerror","allelicdropout"]
            if isnothing(prior_err)
                prior_err = Beta(1, 1/avgerrdict[target] - 1)
            end
        elseif target in ["allelicoverdispersion"]
            if isnothing(prior_err)
                prior_err = Exponential(max(0.1,avgerrdict[target]))
            end
        elseif target in ["allelicbias"]
            if isnothing(prior_err)
                prior_err = Beta(1.01,1.01)
            end
        end

        # @inf
        for kk in eachindex(tseq)
            # if target == "baseerror" && kk==1
            #     println("target=",target,",prior_err=",prior_err)
            # end
            tnow = tseq[kk]
            snp = snporder[tnow]
            issnpGT[snp] && in(target,["baseerror","allelicbias","allelicoverdispersion","allelicdropout"]) && continue
            originls = viterbils[snp]        
            errmarker = infer_errmarker_viterbi!(dataprobls,snp, originls, target, prior_err,
                epsfls[snp],epsols[snp], epso_perind, baseerrorls[snp],
                allelicbiasls[snp],allelicoverdispersionls[snp], allelicdropoutls[snp],                            
                offspringexcl, chrfhaplo, chroffgeno, popmakeup;
                isoffphased, israndallele, issnpGT, itmax,temperature)                  
            if target == "offspringerror"
                epsols[snp] = errmarker                
            elseif target == "foundererror"
                epsfls[snp] = errmarker                
            elseif target == "baseerror"
                baseerrorls[snp] = errmarker                
            elseif target == "allelicbias"
                allelicbiasls[snp] = errmarker    
            elseif target == "allelicoverdispersion"
                allelicoverdispersionls[snp] = errmarker                
            elseif target == "allelicdropout"
                allelicdropoutls[snp] = errmarker                
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
    baseerrormarker::Real, 
    allelicbiasmarker::Real,    
    allelicoverdispersionmarker::Real,    
    allelicdropoutmarker::Real,        
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
        epsf2, epso2, baseerror2,allelicbias2,allelicoverdispersion2,allelicdropout2 = epsfmarker, epsomarker, baseerrormarker,allelicbiasmarker,allelicoverdispersionmarker,allelicdropoutmarker
        if target == "offspringerror"
            epso2 = 1/(1+exp(-x))  # inverse of logit transformation     
        elseif target == "foundererror"
            epsf2 = 1/(1+exp(-x))                      
        elseif target == "baseerror"
            baseerror2 = 1/(1+exp(-x))                     
        elseif target == "allelicbias"
            allelicbias2 = 1/(1+exp(-x))    
        elseif target == "allelicoverdispersion"
            allelicoverdispersion2 = exp(x)
        elseif target == "allelicdropout"
            allelicdropout2 = 1/(1+exp(-x))                     
        else
            @error string("unknown genoerror type: ",target)
        end        
        if target in ["foundererror", "offspringerror", "baseerror","allelicbias","allelicdropout"]
            p =  1/(1+exp(-x)) 
            logpri = logpdf(prior_err,p)            
            # estimation refers to transformed space
            logjacobi = x - 2*log(1+exp(x)) # dp/dx = exp(x)/(1+exp(x))^2              savefig(outstem*"_fig_error.png")
        elseif target == "allelicoverdispersion"    
            p = exp(x)
            logpri = logpdf(prior_err,p)                 
            logjacobi = x   # dp/dx = exp(x) # estimation refers to transformed space       
        end    
        calloglmarker_viterbi!(dataprobls, originls, offspringexcl, epsf2, epso2, epso_perind, baseerror2,allelicbias2,allelicoverdispersion2, allelicdropout2, snp, 
            chrfhaplo, chroffgeno,popmakeup, isoffphased, israndallele, issnpGT) + logpri + logjacobi
    end    
    if target == "offspringerror"
        xstart = log(epsomarker/(1-epsomarker))
    elseif target == "foundererror"
        xstart = log(epsfmarker/(1-epsfmarker))
    elseif target == "baseerror"
        xstart = log(baseerrormarker/(1-baseerrormarker))
    elseif target == "allelicbias"
        xstart = log(allelicbiasmarker/(1-allelicbiasmarker))
    elseif target == "allelicoverdispersion"
        xstart = max(-10.0,log(allelicoverdispersionmarker))
    elseif target == "allelicdropout"
        xstart = log(allelicdropoutmarker/(1-allelicdropoutmarker))
    else
        @error string("unknown genoerror type: ",target)
    end       
    if target == "allelicoverdispersion"
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
    est =  target == "allelicoverdispersion" ? exp(x) : 1.0/(1.0+exp(-x)) 
    est
end

function calloglmarker_viterbi!(dataprobls::AbstractVector,originls::AbstractVector, 
    offspringexcl::AbstractVector, 
    epsf::Real,epso::Real, 
    epso_perind::Union{Nothing,AbstractVector}, 
    baseerror::Real,
    allelicbias::Real,allelicoverdispersion::Real,allelicdropout::Real,
    snp::Integer,    
    chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,    
    isoffphased::Bool,
    israndallele::Bool,
    issnpGT::AbstractVector)
    MagicReconstruct.calsitedataprob_singlephase!(dataprobls,chrfhaplo[snp,:],chroffgeno[snp,:],popmakeup;
        epsf, epso, epso_perind, baseerror,allelicbias, allelicoverdispersion,
        allelicdropout, isoffphased, israndallele, issiteGT = issnpGT[snp])  
    noff = length(dataprobls) 
    offincl = isempty(offspringexcl) ? (1:noff) : setdiff(1:noff, offspringexcl)                
    sum(log.(map((x,i)->x[i],dataprobls[offincl],originls[offincl])))    
end

