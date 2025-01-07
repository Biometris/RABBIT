#################################################

function loglike_ind(epsf_vec::AbstractVector,epso::Real,
    offindex::Integer,subpopid::AbstractString,
    pretuplels::AbstractVector,mapavailable::Bool)
    # epso refers to offindex
    nchr = length(pretuplels)
    loglike = 0.0
    for chr in 1:nchr
        chrissnpGT,chrpopmakeup,chrpriorprocess,chrfderive,chroffcode =pretuplels[chr]
        size(chroffcode,1) <=1 && continue
        ishaploidsub = chrpopmakeup[subpopid]["ishaploid"]
        # nzstate: the set of indices with states being acessible in the subpopulation
        nzstate = chrpopmakeup[subpopid]["nzstate"]
        nzorigin = chrpopmakeup[subpopid]["nzorigin"]
        off = offindex
        obsseq = chroffcode[:,off]
        # println("off=", off, "; loglike_ind before dataprob_given_state")
        dataprob = dataprob_given_state(obsseq, chrfderive, epsf_vec,epso,
            nzstate,nzorigin, chrissnpGT,ishaploidsub)
        # println("off=", off, "; loglike_ind after dataprob_given_state")
        if mapavailable
            hashcode = chrpopmakeup[subpopid]["hashcode"]
            _, initprob, tranprobseq = MagicReconstruct.hashcode2prior(chrpriorprocess,hashcode)
            dataprobseq = Vector{Float632}.(eachrow(dataprob))
            _,fwscale = HMM.forward(initprob,tranprobseq,dataprobseq)
            loglike += HMM.calloglike(fwscale)
        else
            initprob = chrpopmakeup[subpopid]["initprob"]
            loglike += sum([log(sum(initprob .* i)) for i in eachrow(dataprob)])
        end
    end
    loglike
end

function loglike_popls(epsf_vec::AbstractVector,
    epso_vec::AbstractVector,popidls,
    pretuplels::AbstractVector,mapavailable::Bool)
    nchr = length(pretuplels)
    res = zeros(nchr)    
    ThreadsX.foreach(eachindex(pretuplels,res)) do i
        res[i] = loglike_popls_chr(epsf_vec,epso_vec,popidls,pretuplels[i],mapavailable)
    end
    sum(res)
end

function loglike_popls_chr(epsf_vec::AbstractVector,
    epso_vec::AbstractVector,popidls,pretuple::NamedTuple,mapavailable::Bool)
    chrissnpGT,chrpopmakeup,chrpriorprocess,chrfderive,chroffcode = pretuple;
    loglike = spzeros(length(epso_vec))
    size(chroffcode,1) <=1 && return loglike
    for popid in popidls
        ishaploidsub = chrpopmakeup[popid]["ishaploid"]
        offls = chrpopmakeup[popid]["offspring"]
        # nzstate: the set of indices with states being acessible in the subpopulation
        nzstate = chrpopmakeup[popid]["nzstate"]
        nzorigin = chrpopmakeup[popid]["nzorigin"]
        if mapavailable
            hashcode = chrpopmakeup[popid]["hashcode"]
            _, initprob, tranprobseq = MagicReconstruct.hashcode2prior(chrpriorprocess,hashcode)
        else
            initprob = chrpopmakeup[popid]["initprob"]
        end
        for off in offls
            obsseq = chroffcode[:,off]
            # println("off=", off, "; loglike_popls_chr before dataprob_given_state")
            dataprob = dataprob_given_state(obsseq, chrfderive, epsf_vec,epso_vec[off],
                nzstate,nzorigin, chrissnpGT,ishaploidsub)
            # println("off=", off, "; loglike_popls_chr after dataprob_given_state")
            if mapavailable
                dataprobseq = Vector{Float64}.(eachrow(dataprob))
                _,fwscale = HMM.forward(initprob,tranprobseq,dataprobseq)
                loglike[off] = HMM.calloglike(fwscale)
            else
                loglike[off] = sum([log(sum(initprob .* i)) for i in eachrow(dataprob)])
            end
        end
    end
    sum(loglike)
end

function calpretuplels(magicgeno::MagicGeno,
    magicprior::NamedTuple;
    model::AbstractString,
    isfounderinbred::Bool)
    nchr = length(magicgeno.markermap)
    [begin
        nsnp = size(magicgeno.markermap[chr],1)
        popmakeup, priorprocess = MagicReconstruct.calpriorprocess(magicgeno,chr,model,magicprior; isfounderinbred)
        issnpGT = trues(nsnp)
        # [occursin("GT",i) for i in magicgeno.markermap[chr][!,:offspringformat]]
        fderive = MagicReconstruct.precompute_fderive(magicgeno.foundergeno[chr],popmakeup);
        calledgeno = call_offgeno(magicgeno, model,chr)
        offcode = MagicReconstruct.precompute_offcode(calledgeno,popmakeup,issnpGT);
        (issnpGT=issnpGT,popmakeup=popmakeup,priorprocess=priorprocess,fderive=fderive,offcode=offcode)
    end for chr in 1:nchr]
end

#################################################

function dataprob_given_state(obsseq::AbstractVector, fderive::NamedTuple,
    epsf_vec::AbstractVector, epso::Real,
    nzstate::AbstractVector,nzorigin::AbstractVector,
    issnpGT::AbstractVector, ishaploid::Bool)
    # like_true: matrix of size nssnp x nstate
    like_true = like_given_true(obsseq, epso,issnpGT,ishaploid)
    b = sum(like_true,dims=2)[:,1] .== 0.0
    like_true[b,:] .= 1.0
    like_state = like_given_state(like_true, fderive,epsf_vec,
        nzstate, nzorigin, ishaploid)
    like_state
end

function like_given_state(like_true::AbstractMatrix, fderive::NamedTuple,
    epsf_vec::AbstractVector,nzstate::AbstractVector,
    nzorigin::AbstractVector, ishaploid::Bool)
    if ishaploid
        # nstate = nfgl for ishaploid = true
        # length(espf_vec)=nfounder
        # length(epsf_vec)== size(fderive.haplo,2) || @error "dimension mismatch"
        like_state = reduce(hcat,[begin
            true_derive = MagicReconstruct.haploprior(epsf_vec[s])
            true_state = true_derive[:,fderive.haplo[:,s]]'
            sum(like_true .* true_state,dims=2)
        end for s in nzstate])
    else
        like_state = reduce(hcat,[begin
            origin_pair = nzorigin[i]
            epsf_pair = epsf_vec[origin_pair]
            true_derive = MagicReconstruct.diploprior(origin_pair, epsf_pair)
            true_state = true_derive[:,fderive.diplo[:,nzstate[i]]]'
            sum(like_true .* true_state,dims=2)
        end for i in 1:length(nzorigin)])
    end
    like_state
end

# probabily of observed genotypes given true genotypes
function like_given_true(obsseq::AbstractVector, epso::Real,
    issnpGT::AbstractVector,ishaploid::Bool)
    issnpnonGT = .!issnpGT
    nsnp = length(obsseq)
    if ishaploid
        like = zeros(nsnp, 2)
        if any(issnpnonGT)
            like[issnpnonGT,:] = MagicReconstruct.haplolikeGBS(obsseq[issnpnonGT],epso)
        end
        if any(issnpGT)
            haplolike_epso = MagicReconstruct.haplolike(epso)
            like[issnpGT,:] = haplolike_epso[obsseq[issnpGT],:]
        end
    else
        like = zeros(nsnp, 4)
        if any(issnpnonGT)
            like[issnpnonGT,:] = MagicReconstruct.diplolikeGBS(obsseq[issnpnonGT],epso)
        end
        if any(issnpGT)
            diplolike_epso = MagicReconstruct.diplolike(epso)
            like[issnpGT,:] = diplolike_epso[obsseq[issnpGT],:]
        end
    end
    like
end


#################################################
