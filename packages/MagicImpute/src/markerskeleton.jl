
function rescalemap!(priorprocess::AbstractDict,skeletonprior::AbstractDict)
    pri0=first(values(priorprocess))
    incl0 = findall(pri0.markerincl)
    deltd0=deepcopy(pri0.markerdeltd)
    pri1=first(values(skeletonprior))
    incl1 = findall(pri1.markerincl)
    deltd1=pri1.markerdeltd
    issubset(incl1,incl0) || @warn string("unexpected skeleton incl")
    if !in(incl0[1], incl1)
        dlshead = deltd0[1:incl0[1]-1]
        b = [isnothing(i) || i ≈ 0.0 for i in dlshead]
        all(b) || @warn string("first marker is not in skeleton")
    end
    if !in(incl0[end], incl1)
        dlstail = deltd0[incl0[end]:end]
        b = [isnothing(i) || i ≈ 0.0 for i in dlstail]
        all(b) || @warn string("last marker is not in skeleton")
    end
    for i=1:length(incl1)-1
        s=incl1[i]:incl1[i+1]-1        
        d0 = deltd0[s]
        d1 = deltd1[s]
        temp=sum(d0[.!isnothing.(d0)])
        scale = temp ≈ 0 ? 0 : sum(d1[.!isnothing.(d1)])/temp        
        deltd0[s] .= [isnothing(j) ? j : round(j*scale,digits=6) for j in deltd0[s]]
    end
    # sum(deltd0) ≈ sum(deltd1) || @error string("wrong re-scale by skeleton map")
    for (_, pri) in priorprocess
        pri.markerdeltd = deltd0        
        pri.tranprobseq = [MagicReconstruct.get_tranprobmatrix(i,pri.tranrate) for i = deltd0]
    end
    priorprocess
end

function markerskeleton_chr(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    isinfererror::Bool,    
    priorlikeparam::PriorLikeParam,
    liketargetls::AbstractVector,
    likeerrortuple::NamedTuple,     
    offspringexcl::AbstractVector, 
    snporder::AbstractVector,    
    israndallele::Bool,     
    issnpGT::AbstractVector, 
    isinferseq::Bool,
    decodetempfile::AbstractString,        
    skeletonsize::Union{Nothing,Integer},    
    msg_repeatrun::AbstractString, 
    chrid::AbstractString="", 
    startit::Integer=1,
    spacebyviterbi::Bool=false,     
    verbose::Bool=true,
    io::Union{Nothing,IO}=nothing)
    (;epsfls, epsols, epsols_perind, baseerrorls,allelicbiasls,allelicoverdispersionls, allelicdropoutls) = likeerrortuple
    epso_perind = isnothing(epsols_perind) ? nothing : copy(epsols_perind)
    skeleton_epsfls = isa(epsfls,AbstractVector) ? copy(epsfls) : [epsfls for _ in eachindex(snporder)]
    skeleton_epsols = isa(epsols,AbstractVector) ? copy(epsols) : [epsols for _ in eachindex(snporder)]
    skeleton_seqerrls = isa(baseerrorls,AbstractVector) ? copy(baseerrorls) : [baseerrorls for _ in eachindex(snporder)]
    skeleton_allelicbiasls = isa(allelicbiasls,AbstractVector) ? copy(allelicbiasls) : [allelicbiasls for _ in eachindex(snporder)]
    skeleton_allelicoverdispersionls = isa(allelicoverdispersionls,AbstractVector) ? copy(allelicoverdispersionls) : [allelicoverdispersionls for _ in eachindex(snporder)]
    skeleton_allelicdropoutls = isa(allelicdropoutls,AbstractVector) ? copy(allelicdropoutls) : [allelicdropoutls for _ in eachindex(snporder)]
    skeleton_errls = [skeleton_epsols, skeleton_seqerrls, skeleton_allelicbiasls, skeleton_allelicoverdispersionls, skeleton_allelicdropoutls]
    skeleton,nmarkerbin,nmarkerbin2 = getskeleton_seg(skeleton_errls, isinferseq, snporder, priorprocess, skeletonsize)    
    skeletonprior = deepcopy(priorprocess)
    pri1=first(values(skeletonprior))
    if 0 in pri1.markerincl[skeleton]
        @error string("deleted markers are not allowed in skeleton")
    end
    deltt = trues(length(pri1.markerincl))
    deltt[skeleton] .= false
    MagicReconstruct.setpriorprocess!(skeletonprior,deltt)
    
    pri1=first(values(skeletonprior))
    exclsnps = snporder[.!pri1.markerincl]
    skeleton_epsfls[exclsnps] .= NaN
    skeleton_epsols[exclsnps] .= NaN
    skeleton_seqerrls[exclsnps] .= NaN
    skeleton_allelicbiasls[exclsnps] .= NaN
    skeleton_allelicoverdispersionls[exclsnps] .= NaN
    skeleton_allelicdropoutls[exclsnps] .= NaN
    chrlenls = []
    # inclsnps = snporder[pri1.markerincl]
    # df = DataFrame(offerror = skeleton_epsfls[inclsnps],
    #     baseerror = skeleton_seqerrls[inclsnps],
    #     abmean = skeleton_allelicbiasls[inclsnps],
    #     abshape = skeleton_allelicoverdispersionls[inclsnps],
    #     ado = skeleton_allelicdropoutls[inclsnps])
    # println("skeleton,df=",df)
    infer_err_fun = infer_errorls_viterbi! 
    markerspace_fun = spacebyviterbi ? markerspace_chr_viterbi! : markerspace_chr!
    for it in startit:startit+9
        startt = time()
        msg = string("chr=", chrid, msg_repeatrun, ", it=", it, ", #skeleton=",length(skeleton),", #segment=",nmarkerbin)        
        if isinfererror            
            targetls = intersect(["foundererror","offspringerror"],liketargetls)       
            
            pri1=first(values(skeletonprior))
            snpincl = snporder[pri1.markerincl]		
            errnamels = ["foundererror", "offspringerror", "baseerror","allelicoverdispersion","allelicdropout"]            
            avgerrls = [mean(i[snpincl]) for i in [skeleton_epsfls,skeleton_epsols,skeleton_seqerrls,
                skeleton_allelicoverdispersionls,skeleton_allelicdropoutls]]
            avgerrdict = Dict(errnamels .=> avgerrls)   		
            
            infer_err_fun(chrfhaplo,chroffgeno, popmakeup,skeletonprior;
                targetls, priorlikeparam,
                epsf=skeleton_epsfls,epso=skeleton_epsols, epso_perind, 
                baseerror = skeleton_seqerrls,
                allelicbias = skeleton_allelicbiasls,
                allelicoverdispersion = skeleton_allelicoverdispersionls,    
                allelicdropout = skeleton_allelicdropoutls,  avgerrdict,           
                offspringexcl, snporder,decodetempfile, israndallele, issnpGT, temperature = 0.0,itmax=20)        
            # mepsf = mean(skeleton_epsfls[.!isnan.(skeleton_epsfls)])
            mepso = mean(skeleton_epsols[.!isnan.(skeleton_epsols)])            
            msg *= string(", εo=",round(mepso,digits=4))
            # if in("peroffspringerror",liketargetls)		                
            #     infer_error_perind_viterbi!(chrfhaplo, chroffgeno, popmakeup,skeletonprior;
            #         epsf=skeleton_epsfls,epso=skeleton_epsols, epso_perind, 
            #         baseerror = skeleton_seqerrls,
            #         allelicbias = skeleton_allelicbiasls,
            #         allelicoverdispersion = skeleton_allelicoverdispersionls,    
            #         allelicdropout = skeleton_allelicdropoutls,    
            #         prior_peroffspringerror = priorlikeparam.peroffspringerror,         
            #         israndallele, issnpGT,snporder,decodetempfile, temperature = 0.0,itmax=10)	
            #     mepso_perind = mean(epso_perind)
            #     msg *= string(", ξo=", round(mepso_perind,digits=3))			                
            # end
            if isinferseq	
                targetls = setdiff(liketargetls,["foundererror","offspringerror","peroffspringerror"])                
                infer_err_fun(chrfhaplo,chroffgeno, popmakeup,skeletonprior;
                    targetls,  priorlikeparam,
                    epsf=skeleton_epsfls, epso=skeleton_epsols, epso_perind, 
                    baseerror = skeleton_seqerrls,allelicbias = skeleton_allelicbiasls,
                    allelicoverdispersion = skeleton_allelicoverdispersionls,
                    allelicdropout = skeleton_allelicdropoutls, avgerrdict,            
                    offspringexcl,snporder,decodetempfile, israndallele, issnpGT, temperature = 0.0,itmax=20)
                0
                if in("baseerror", targetls)
                    mseqerr = mean(skeleton_seqerrls[.!isnan.(skeleton_seqerrls)])
                    msg *= string(", εbase=",round(mseqerr,digits=4))
                end
                if in("allelicbias", targetls)
                    mbalance = mean(abs.(skeleton_allelicbiasls[.!isnan.(skeleton_allelicbiasls)]))                
                    msg *= string(", bias=",round(mbalance,digits=2))
                end
                if in("allelicoverdispersion", targetls)
                    mdisperse = mean(abs.(skeleton_allelicoverdispersionls[.!isnan.(skeleton_allelicoverdispersionls)]))                
                    msg *= string(", overdisperse=",round(mdisperse,digits=3))
                end
                if in("allelicdropout", targetls)
                    mdropout = mean(abs.(skeleton_allelicdropoutls[.!isnan.(skeleton_allelicdropoutls)]))
                    mdropout > 0 && (msg *= string(", dropout=",round(mdropout,digits=4)))
                end                
            end
        end
        markerspace_fun(chrfhaplo,chroffgeno, popmakeup,skeletonprior;
            epsf=epsfls,epso=skeleton_epsols, epso_perind, 
            baseerror = skeleton_seqerrls,allelicbias = skeleton_allelicbiasls,
            allelicoverdispersion = skeleton_allelicoverdispersionls,
            allelicdropout = skeleton_allelicdropoutls,            
            offspringexcl, snporder,reversechr = iseven(it), israndallele, issnpGT, 
            decodetempfile = decodetempfile,temperature=0.0) 
        chrlen = 100*sum(pri1.markerdeltd[pri1.markerincl][1:end-1])
        push!(chrlenls,chrlen)
        msg *= string(", l=",round(Int,chrlen), "cM", ", t=",round(Int,time()-startt),"s")
        printconsole(io, verbose, msg)            
        if length(chrlenls) >= 3 
            diff_chrlen = min(abs.(chrlen .- chrlenls[end-2:end-1])...)
            if diff_chrlen<=1 || diff_chrlen/chrlen < 0.005
                break
            end
        end        
    end
    skeletonprior,chrlenls
end

function getskeleton_seg(errorls::AbstractVector, isinferseq::Bool,
    snporder::AbstractVector,
    priorprocess::AbstractDict,    
    skeletonsize::Union{Nothing,Integer})
    if isinferseq
        epsols, baseerrorls, allelicbiasls, allelicoverdispersionls, allelicdropoutls = errorls
        tterrorls=[epsols[i] + baseerrorls[i] + 0.001 * (abs(allelicbiasls[i] - 0.5) + allelicoverdispersionls[i] + allelicdropoutls[i]) for i in snporder]
    else
        epsols = first(errorls)
        tterrorls=epsols[snporder]
    end    
    pri1=first(values(priorprocess))
    excl = findall(.!pri1.markerincl)
    deltd= copy(pri1.markerdeltd)
    deltd[isnothing.(deltd)] .= 0.0
    loc =accumulate(+,deltd[1:end-1])
    pushfirst!(loc,0.0)
    bins = MagicBase.splitindex(loc)            
    binls = [setdiff(collect(i),excl) for i=bins]
    binls = binls[length.(binls) .> 0]    
    nmarkerbin = length(binls)
    skeleton=[bin[argmin(tterrorls[bin])] for bin=binls]        
    skeletonerr = [min(tterrorls[bin]...) for bin=binls]  
    # delete bin with too large error    
    skefirst = first(skeleton)
    skelast = last(skeleton)
    b = skeletonerr .<= 0.1
    if 1.0 > mean(b) > 0.7
        skeleton = skeleton[b]
        skeletonerr = skeletonerr[b]
        binls = binls[b]
        b = skeletonerr .<= 0.05
        # no need to update skeletonerr
        if 1.0 > mean(b) > 0.7 && (isnothing(skeletonsize) || sum(b) >= skeletonsize)            
            skeleton = skeleton[b]            
            binls = binls[b]            
        end          
        if !in(skefirst, skeleton) 
            pushfirst!(skeleton, skefirst)
            pushfirst!(binls, [skefirst])
        end
        if !in(skelast, skeleton) 
            push!(skeleton, skelast)
            push!(binls, [skelast])
        end     
    end
    nmarkerbin2 = length(binls)
    if isnothing(skeletonsize) || length(skeleton)<=skeletonsize        
        return skeleton,nmarkerbin,nmarkerbin2
    end
    skeldeltd=deltd[[bin[end] for bin=binls]]
    skeldeltd[end]==0 || @error "unexpected last skeleton deltd"
    if nmarkerbin == nmarkerbin2 
        dsumdiff = abs(sum(skeldeltd) - sum(deltd))
        dsumdiff < 1e-3 || @error string("chrlen_diff = ", dsumdiff, " > 1e-3")
    end
    skelloc=accumulate(+,circshift(skeldeltd,1))
    gridsize=skelloc[end]/(skeletonsize-2)
    grid=gridpartition(skelloc,gridsize)
    grid=grid[length.(grid) .> 0]
    grid=[skeleton[i] for i=grid]
    skeleton2=[bin[argmin(tterrorls[bin])] for bin=grid]
    in(skefirst, skeleton2) && pushfirst!(skeleton2, skefirst)
    in(skelast, skeleton2) && push!(skeleton2, skelast)            
    ndiff = skeletonsize-length(skeleton2)
    if ndiff>0
        ii=setdiff(skeleton,skeleton2)
        s=sortperm(tterrorls[ii])
        length(s) > ndiff && (s=s[1:ndiff])
        skeleton2=sort(vcat(skeleton2,ii[s]))
    end
    skeleton2, nmarkerbin,nmarkerbin2
end

function gridpartition(xls::AbstractVector,gridsize::Real)
    xls .-= xls[1]
    grid=collect(0:gridsize:xls[end])
    grid[end] += gridsize+0.1
    res=Vector{Vector}()
    b=i=1
    bin=[]
    while true
        if grid[b]<=xls[i]<grid[b+1]
            push!(bin,i)
            i+=1
        else
            push!(res,bin)
            b+=1
            bin=[]
        end
        i>length(xls) && break
    end
    push!(res,bin)
    res
end
