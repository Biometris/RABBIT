
const _float_tran = Float64
const _float_like = Float64

mutable struct MarkovPrior
    # for a given linkage group
    initprob::Vector{Float64}
    tranrate::Matrix{Float64}
    tranprobseq::Vector{<: Union{Nothing,Matrix{_float_tran}}}
    markerdeltd::Vector{Union{Nothing,Float64}}
    markerincl::BitVector
    markerid::Vector{String}
    function MarkovPrior(initprob,tranrate, tranprobseq,markerdeltd,markerincl,markerid)
        dims = length.([markerdeltd,markerincl,markerid])
        if length(union(dims))!=1
            @error("inconsistent number of markers: ",dims)
        end
        dims = [length(initprob),size(tranrate,1),size(tranrate,2)]
        if length(union(dims))!=1
            @error("inconsistent number of states: ",dims)
        end
        new(initprob,tranrate, tranprobseq,markerdeltd,markerincl,markerid)
    end
end

function reverseprior!(pri::MarkovPrior)
    incl = pri.markerincl
    deltd = circshift(reverse(pri.markerdeltd[incl]),-1)
    tran = circshift(reverse(pri.tranprobseq[incl]),-1)
    reverse!(pri.markerincl) # incl is revised too
    pri.markerdeltd[incl] .= deltd
    pri.markerdeltd[.!incl] .= nothing
    pri.tranprobseq[incl] .= tran
    pri.tranprobseq[.!incl] .= nothing
    reverse!(pri.markerid)
    # check consistency
    # inclsnps = findall(incl)
    # if length(inclsnps) > 1
    #     if any(isnothing.(pri.markerdeltd[inclsnps[1:end-1]]))
    #         @error string("inconsistent markerincl and markerdeltd")
    #     end
    #     if any(isnothing.(pri.tranprobseq[inclsnps[1:end-1]]))
    #         @error string("inconsistent markerincl and tranprobseq")
    #     end
    # end
    pri
end


function uppriorprocess(priorprocess::AbstractDict, popmakeup::AbstractDict)    
    pri1=first(values(priorprocess))
    newpriorprocess = Dict()
    for popid in keys(popmakeup)
        initprob = popmakeup[popid]["initprob"]
        tranrate = popmakeup[popid]["tranrate"]
        hashcode=popmakeup[popid]["hashcode"]
        if !haskey(newpriorprocess, hashcode)            
            # copy markerdeltd, markerincl and markerid for each subpop
            # so that reverseprior will be indepedent for each subpop.
            # round/isapprox is to avoid ERROR: LinearAlgebra.LAPACKException(1)            
            tranprobseq=[get_tranprobmatrix(i,tranrate) for i = pri1.markerdeltd]
            markovprior = MarkovPrior(initprob,tranrate,tranprobseq,
                copy(pri1.markerdeltd), 
                copy(pri1.markerincl),
                copy(pri1.markerid))
            push!(newpriorprocess,hashcode=>markovprior)
        end
    end
    newpriorprocess
end

function calpriorprocess(magicgeno::MagicGeno, chr::Integer, popmakeup::AbstractDict)
    markerdeltd = diff(magicgeno.markermap[chr][!,:poscm]) * 0.01
    markerid = magicgeno.markermap[chr][!,:marker]
    markerincl = trues(length(markerid))
    priorprocess = Dict()
    for popid in keys(popmakeup)
        initprob = popmakeup[popid]["initprob"]
        tranrate = popmakeup[popid]["tranrate"]
        hashcode=popmakeup[popid]["hashcode"]
        if !haskey(priorprocess, hashcode)            
            # copy markerdeltd, markerincl and markerid for each subpop
            # so that reverseprior will be indepedent for each subpop.
            tranprobseq=[get_tranprobmatrix(i,tranrate) for i = markerdeltd]
            markovprior = MarkovPrior(initprob,tranrate,
                vcat(tranprobseq,nothing),
                copy(vcat(markerdeltd,nothing)),
                copy(markerincl),
                copy(markerid))
            push!(priorprocess,hashcode=>markovprior)
        end
    end
    priorprocess
end

function get_tranprobmatrix(d::Union{Nothing,Real},tranrate::AbstractMatrix)
    if isnothing(d)
        nothing
    else        
        if d <= 0.0
            nothing
        elseif d <= 1e-5
            res = Matrix{_float_tran}(d .* tranrate )
            res .+= 0.5 .* (res * res)
            res + I
        else            
            # round digits=6 to avoid ERROR: LinearAlgebra.LAPACKException(1) for matrix exp
            d2 = round(d,digits=6)
            Matrix{_float_tran}(real.(exp(d2 .* tranrate)))                         
        end
    end
end

function calpriorprocess(magicgeno::MagicGeno,chr::Integer,
    model::AbstractString,magicprior::NamedTuple;
    isfounderinbred::Bool=true)
    popmakeup = calpopmakeup(magicgeno,chr, model,magicprior; isfounderinbred)
    priorprocess = calpriorprocess(magicgeno,chr,popmakeup)
    popmakeup,priorprocess
end

function hashcode2prior(priorprocess::AbstractDict, hashcode::Integer)
    ls=priorprocess[hashcode]
    tranprobseq = view(ls.tranprobseq, ls.markerincl)
    (markerincl=ls.markerincl,initprob=ls.initprob, tranprobseq=tranprobseq)
end

function setpriorprocess!(priorprocess::AbstractDict, snporder::AbstractVector,
    delsnps::AbstractVector)
    dict = Dict(snporder .=> 1:length(snporder))
    deltt = [get(dict, snp, nothing) for snp in delsnps]
    setpriorprocess!(priorprocess,deltt)
end

function setpriorprocess!(priorprocess::AbstractDict, deltt::AbstractVector)
    for (_, pri) in priorprocess
         setmarkerincl!(pri,deltt)
    end
    priorprocess
end


function setmarkerincl!(pri::MarkovPrior,deltt::AbstractVector)
    markerincl = pri.markerincl
    markerincl[deltt] .= false
    # calculate new deltd
    deltd = copy(pri.markerdeltd)
    i0 = findfirst(markerincl)
    if isnothing(i0)
        pri.markerdeltd[.!markerincl] .= nothing
        pri.tranprobseq[.!markerincl] .= nothing
    else
        pre = i0
        for i=i0:length(markerincl)-1
            if markerincl[i]
                pre = i
            else
                if !isnothing(deltd[i])
                    deltd[pre] += deltd[i]
                    deltd[i] = nothing
                end
            end
        end
        markerincl[end] || (deltd[pre]=0.0)
        # println("pre=",pre,";i0=",i0,",deltt=",deltt)
        # update pri.markerdeltd
        pos = findall(deltd[1:end-1] .!== pri.markerdeltd[1:end-1])
        pri.markerdeltd[pos] .= deltd[pos]
        pos = pos[markerincl[pos]]
        if !isempty(pos)
            pri.tranprobseq[pos] .= [get_tranprobmatrix(i,pri.tranrate) for i in deltd[pos]]
        end
    end
    pri
end


function check_prior_consistency(magicgeno::MagicGeno,chr::Integer,
    priorprocess::AbstractDict,snporder::AbstractVector)
    markerid_map = magicgeno.markermap[chr][snporder,:marker]
    for (subpop,pri) in priorprocess
        markerid_map == pri.markerid || @error string("inconsistent markerid for subpop=",subpop)
        incl = pri.markerincl
        inclsnps = findall(incl)
        if length(inclsnps) > 1
            if any(isnothing.(pri.markerdeltd[inclsnps[1:end-1]]))
                @error string("inconsistent markerincl and markerdeltd")
            end
            # if any(isnothing.(pri.tranprobseq[inclsnps[1:end-1]]))
            #     @error string("inconsistent markerincl and tranprobseq")
            # end
        end
        deltd = pri.markerdeltd[incl][1:end-1]
        if !isempty(deltd)            
            transeq=[get_tranprobmatrix(i,pri.tranrate) for i = deltd]
            issametran = map((x,y)->any(isnothing.([x,y])) || isapprox(x,y,atol=1e-4), transeq, pri.tranprobseq[incl][1:end-1])
            all(issametran) || @warn string("inconsistent markerdeltd and tranprobseq for hashcode=",subpop)
        end
        prils = values(priorprocess)
        length(unique([pri.markerincl for pri in prils])) == 1 || @error "inconsistent markerincl"
        length(unique([pri.markerid for pri in prils]))==1 || @error "inconsistent markerid"
        # length(unique([pri.markerdeltd for pri in prils]))==1 || @error "inconsistent markerdeltd"
    end
end

function setmagicgeno!(magicgeno::MagicGeno,chr::Integer,
    priorprocess::AbstractDict,snporder::AbstractVector)
    # book deleted markers
    pri1=first(values(priorprocess))
    delsnps = snporder[.!pri1.markerincl]
    if !isempty(delsnps)
        deletiondf = deepcopy(magicgeno.markermap[chr][delsnps,:])
        MagicBase.pushmisc!(magicgeno, "deletion"=>deletiondf)
    end
    # keep un-deleted markers
    keepsnps = snporder[pri1.markerincl]
    magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][keepsnps,:]
    magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][keepsnps,:]
    magicgeno.markermap[chr] = magicgeno.markermap[chr][keepsnps,:]    
    # b = ismissing.(magicgeno.markermap[chr][!,:physposbp])
    # initpos = all(b) ? 0.0 : magicgeno.markermap[chr][1,:poscm] # change initpos will affect merge_missprogeny
    initpos = magicgeno.markermap[chr][1,:poscm] 
    deltd = pri1.markerdeltd[pri1.markerincl][1:end-1]
    snppos = vcat([initpos],initpos .+ 100*accumulate(+,deltd))
    magicgeno.markermap[chr][!,:poscm] .= snppos    
    magicgeno
end
