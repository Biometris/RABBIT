# genotypes represented as readcount
function haplolike_depths(allelicdepths::AbstractVector, epsprob::Real,baseerror::Real)
    allelicdepths == [0,0] && return [1.0 1.0]
    n1,n2 = allelicdepths
    logbinom = first(logabsbinomial(n1 + n2, n1))
    logp1 = n1 * log(1-baseerror) +  n2 * log(baseerror) + logbinom
    logp2 = n1 * log(baseerror) + n2 * log(1-baseerror) + logbinom
    logp12 = n1 * log(0.5) + n2 * log(0.5) + logbinom
    logp12 > max(logp1,logp2)
    if logp12 > max(logp1,logp2)
        homo = sum(exp.([logp1,logp2] .- logp12))
        homo /= (1+homo)
        homo <= 0.5 && return [1.0 1.0]
    end
    maxlogp = max(logp1,logp2)
    p1 = exp(logp1 - maxlogp)
    p2 = exp(logp2 - maxlogp)
    pp1 = (1-epsprob)*p1 + epsprob*p2
    pp2 = (1-epsprob)*p2 + epsprob*p1    
    # lowerbound of maxlogp in case of large depths e.g. [1000,2000], exp(-745.0) > 0, exp(-746.0) == 0, assume epsprob > exp(-45.0)
    # might have some effects on inferring error rates
    pscale = exp(max(maxlogp,-700.0))  
    [pp1*pscale pp2*pscale]
end    

# genotypes represented as readcount
function haplolikeGBS(obsseq::AbstractVector, epsprob::Union{Real,AbstractVector},baseerror::Union{Real,AbstractVector})
    res = haplolike_depths.(obsseq,epsprob,baseerror)
    reduce(vcat,res)
end

# genotypes represented as probability vector
function haplolikeGBS(obsseq::AbstractVector,epsprob::Union{Real,AbstractVector})
    p1 = first.(obsseq)
    p2 = last.(obsseq)
    pp1 = @. (1-epsprob)*p1 + epsprob*p2
    pp2 = @. (1-epsprob)*p2 + epsprob*p1
    [pp1 pp2]
end


# genotypes represented as readcount
# allelebalance = probability of read being allele 1 given true genotype being heterozygous
function diplolike_depths(allelicdepths::AbstractVector, epsprob::Real,
    baseerror::Real,allelicbias::Real,allelicoverdispersion::Real,allelicdropout::Real;
    israndallele::Bool)
    allelicdepths == [0,0] && return [1.0 1.0 1.0 1.0]
    n1,n2 = allelicdepths
    logbinom = first(logabsbinomial(n1 + n2, n1))
    logp11 = n1 * log(1-baseerror) +  n2 * log(baseerror) + logbinom
    logp22 = n1 * log(baseerror) + n2 * log(1-baseerror) + logbinom    
    if allelicoverdispersion < 1e-5
        # no overdispersion
        logp12 = n1 * log(allelicbias) + n2 * log(1-allelicbias) + logbinom
    else
        alpha = allelicbias/allelicoverdispersion
        beta = (1.0-allelicbias)/allelicoverdispersion
        logp12 = logbeta(n1+alpha,n2+beta) - logbeta(alpha,beta) + logbinom
    end
    if allelicdropout > 1e-6
        # 0 and 1 inflated 
        ls = [logp11,logp12,logp22] 
        logpmax =maximum(ls)
        @. ls = exp(ls - logpmax)
        p12 = (ls[1] + ls[3])* allelicdropout/2 + ls[2]*(1.0-allelicdropout)
        logp12 = log(p12) + logpmax
    end 
    maxlogp = max(logp11,logp22,logp12)
    p11 = exp(logp11-maxlogp)
    p22 = exp(logp22-maxlogp)
    p12 = exp(logp12-maxlogp)
    p21 = p12
    if israndallele
        # genotyping error model: random allelic model    
        epsls =[(1-epsprob)^2, (1-epsprob)*epsprob, (1-epsprob)*epsprob, epsprob^2]
    else        
        # genotyping error model: random genotypic model    
        epsls =[1-epsprob, epsprob/3, epsprob/3, epsprob/3]
    end
    pp11=sum([p11,p12,p21,p22] .* epsls)
    pp12=sum([p12,p11,p22,p21] .* epsls)
    pp21=sum([p21,p11,p22,p12] .* epsls)
    pp22=sum([p22,p12,p21,p11] .* epsls)
    # lowerbound of maxlogp in case of large depths e.g. [1000,2000], exp(-745.0) > 0, exp(-746.0) == 0, assume epsprob > exp(-45.0)
    # might have some effects on inferring error rates
    pscale = exp(max(maxlogp,-700.0))  
    [pp11 pp12 pp21 pp22] .* pscale
end

function diplolikeGBS(obsseq::AbstractVector,epsprob::Union{Real,AbstractVector},
    baseerror::Union{Real,AbstractVector},
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},
    allelicdropout::Union{Real,AbstractVector};
    israndallele::Bool)
    res = diplolike_depths.(obsseq,epsprob,baseerror,allelicbias,allelicoverdispersion,allelicdropout; israndallele)
    reduce(vcat,res)
end

# genotypes represented as probability vector
function diplolikeGBS(obsseq::AbstractVector,epsprob::Real;
    israndallele::Bool)    
    y=reduce(hcat,obsseq)
    d = size(y,1)
    if d == 2
        # epsprob is regarded as genotyping prob instead of allelic prob
        p11,p22 = [i for i=eachrow(y)]
        pp11 = p11 .* (1-epsprob) + p22 .* epsprob
        pp22 = p22 .* (1-epsprob) + p11 .* epsprob
        n = size(y,2)
        [pp11 zeros(n) zeros(n) pp22]
    elseif d in [3,4]
        if d==3
            p11,p12,p22 = [i for i=eachrow(y)]
            p12 ./= 2
            p21 = p12
        else
            p11,p12,p21,p22 = [i for i=eachrow(y)]
        end
        if israndallele
            # genotyping error model: random allelic model    
            epsls =[(1-epsprob)^2, (1-epsprob)*epsprob, (1-epsprob)*epsprob, epsprob^2]
        else            
            # genotyping error model: random genotypic model    
            epsls =[1-epsprob, epsprob/3, epsprob/3, epsprob/3]
        end
        pp11=sum([p11,p12,p21,p22] .* epsls)
        pp12=sum([p12,p11,p22,p21] .* epsls)
        pp21=sum([p21,p11,p22,p12] .* epsls)
        pp22=sum([p22,p12,p21,p11] .* epsls)
        [pp11 pp12 pp21 pp22]
    else 
        @error string("wrong length of genotype probility vector: ", d)    
    end
end

# genotypes represented as probability vector
function diplolikeGBS(obsseq::AbstractVector,epsprob::AbstractVector;
    israndallele::Bool)           
    y = reduce(hcat,obsseq)    
    d = size(y,1)
    if d == 2
        # epsprob is regarded as genotyping prob instead of allelic prob
        p11,p22 = [i for i=eachrow(y)]
        pp11 = p11 .* (1 .- epsprob) + p22 .* epsprob
        pp22 = p22 .* (1 .- epsprob) + p11 .* epsprob
        n = size(y,2)
        [pp11 zeros(n) zeros(n) pp22]
    elseif d in [3,4]
        if d==3
            p11,p12,p22 = [i for i=eachrow(y)]
            p12 ./= 2
            p21 = p12
        else
            p11,p12,p21,p22 = [i for i=eachrow(y)]
        end
        if israndallele
            # genotyping error model: random allelic model    
            epsls =[(1 .- epsprob) .^ 2, (1 .- epsprob) .* epsprob,
                (1 .- epsprob) .* epsprob, epsprob .^ 2]            
        else
            # genotyping error model: random genotypic model    
            epsls =[1 .- epsprob, epsprob ./ 3, epsprob ./ 3, epsprob ./ 3]
            
        end
        ls = [p11,p12,p21,p22]
        pp11 = sum([ls[i] .* epsls[i] for i in 1:4])
        ls = [p12,p11,p22,p21]
        pp12 = sum([ls[i] .* epsls[i] for i in 1:4])
        ls = [p21,p11,p22,p12]
        pp21 = sum([ls[i] .* epsls[i] for i in 1:4])
        ls = [p22,p12,p21,p11]
        pp22 = sum([ls[i] .* epsls[i] for i in 1:4])
        [pp11 pp12 pp21 pp22]
    else 
        @error string("wrong length of genotype probility vector: ", d)    
    end
end

function callinelikeGBS!(dataprobls::AbstractVector, fderive::AbstractMatrix, nzstate::AbstractVector, obsseq::AbstractVector;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    baseerror::Union{Real,AbstractVector},
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},
    allelicdropout::Union{Real,AbstractVector},
    ishaploid::Bool=false,
    israndallele::Bool)
    if ishaploid
        if eltype(first(obsseq)) <: Integer
            # genofrmat = AD
            # haplotypic data do not depend on allelebalance
            like = haplolikeGBS(obsseq,epso,baseerror)            
        else
            # genofrmat = GP
            like = haplolikeGBS(obsseq,epso)
        end
        if typeof(epsf) <: AbstractVector
            prior = [haploprior(i) for i in epsf]
            likehaplo = reduce(vcat,[like[[i],:] * prior[i] for i in eachindex(prior)])
        else
            prior = haploprior(epsf)
            likehaplo=like*prior
        end
        for snp in eachindex(obsseq)
            dataprobls[snp] .= likehaplo[snp,fderive[snp,nzstate]] 
        end
    else
        nfgl0=sqrt(size(fderive,2))
        if nfgl0 % 1 !=0
            error("unexpected number of columns in fderive")
        end
        nfgl=Int(nfgl0)       
        ibdls = [(i-1)nfgl+i for i in 1:nfgl]
        isnzibd = [i in ibdls for i in nzstate]
        nzibdls = nzstate[isnzibd]
        nznonibdls = nzstate[.!isnzibd]  
        if eltype(first(obsseq)) <: Integer            
            # genofrmat = AD            
            like = diplolikeGBS(obsseq,epso,baseerror,allelicbias,allelicoverdispersion,allelicdropout; israndallele)            
        else
            # genofrmat = GP
            like = diplolikeGBS(obsseq,epso; israndallele)                
        end
        if typeof(epsf) <: AbstractVector
            prior = [diploprior(i) for i in epsf]
            nsnp = length(prior)
            likenonibd = reduce(vcat,[like[[i],:] * prior[i].nonibd for i in 1:nsnp])
            likeibd = reduce(vcat,[like[[i],:] * prior[i].ibd for i in 1:nsnp])
        else
            prior = diploprior(epsf)
            likenonibd=like*prior.nonibd
            likeibd=like*prior.ibd
        end
        issparse = length(nzstate)/(nfgl^2) < 0.5 
        for snp in eachindex(obsseq)
            prob = issparse ? spzeros(_float_like, nfgl^2) : zeros(_float_like,nfgl^2)
            prob[nzibdls] .= likeibd[snp,fderive[snp,nzibdls]]
            prob[nznonibdls] .= likenonibd[snp,fderive[snp,nznonibdls]]
            dataprobls[snp] .= Vector(prob[nzstate])                        
        end 
    end    
end
