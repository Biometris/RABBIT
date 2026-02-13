

function founderimpute_chr!(chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict, priorprocess::AbstractDict, fhaplosetpp::AbstractVector;
    findexlist::AbstractVector,    
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing, AbstractVector}, 
    baseerror::Union{Real,AbstractVector},
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},
    allelicdropout::Union{Real,AbstractVector},
    offspringexcl::AbstractVector,
    snporder::AbstractVector,
    israndallele::Bool,     
    issnpGT::AbstractVector,
    reversechr::Bool, 
    upbyhalf::Bool,    
    alwaysaccept::Bool,         
    threshproposal::Real, 
    imputetempfile::AbstractString)    
    nsnp = size(chrfhaplo,1)
    if nsnp == 1
        @warn string("TODO: nsnp=",nsnp)
    end
    nsnp <=1 && return (0, 0)
    nsnp == length(snporder) || @error "dimension mismatch!"
    loglikels = MagicReconstruct.hmm_loglikels(chrfhaplo,chroffgeno,
        popmakeup,priorprocess;
        epsf,epso,epso_perind, baseerror,
        allelicbias,allelicoverdispersion,allelicdropout, 
        decodetempfile=imputetempfile, israndallele,issnpGT,snporder
    )	
    loglikels[offspringexcl] .= 0.0 
    inputloglike = sum(loglikels)
    
    isfounderinbred = length(fhaplosetpp) == size(chrfhaplo,2)    
    newchrfhaplo = deepcopy(chrfhaplo)
    snpincl = snporder[first(values(priorprocess)).markerincl]				    
    if upbyhalf                    
        reversechrls = reversechr ? [true, false] : [false, true]                 
        for reversechr in reversechrls
            # fixing phasing on one half of the chromosome                      
            nhalf = nsnp ÷ 2                
            fhaplosetpp2 = deepcopy(fhaplosetpp)                
            snpls = reversechr ? snporder[1:nhalf] : snporder[nhalf+1:nsnp]
            for p in eachindex(fhaplosetpp2)                    
                if isfounderinbred
                    for j in snpls
                        fhaplosetpp2[p][j] = [newchrfhaplo[j,p]]                        
                    end
                else
                    for j in snpls
                        fhaplosetpp2[p][j] = [newchrfhaplo[j,[2p-1,2p]]]                        
                    end
                end
            end                
            for findex in findexlist                                            
                founderforwardbackward!(findex,
                    newchrfhaplo,chroffgeno,
                    popmakeup,priorprocess,fhaplosetpp2; 
                    epsf,epso,epso_perind,baseerror,
                    allelicbias,allelicoverdispersion,allelicdropout,
                    offspringexcl, snporder,israndallele, issnpGT,
                    reversechr = !reversechr, forwardtempfile=imputetempfile, 
                    threshproposal, 
                )            
            end 
        end
    else    
        for findex in findexlist                             
            founderforwardbackward!(findex,newchrfhaplo,chroffgeno,
                popmakeup,priorprocess,fhaplosetpp; 
                epsf,epso,epso_perind, baseerror,
                allelicbias,allelicoverdispersion,allelicdropout, 
                offspringexcl, snporder,israndallele, issnpGT,
                reversechr, forwardtempfile=imputetempfile, threshproposal)    
        end                
    end    
    bdiff = get_bdiff(chrfhaplo,newchrfhaplo,snpincl)     
    ndiff = sum(bdiff)
    loglikels = MagicReconstruct.hmm_loglikels(newchrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf,epso,epso_perind, baseerror,
        allelicbias,allelicoverdispersion,allelicdropout, 
        decodetempfile=imputetempfile, israndallele,issnpGT,snporder
    )	
    loglikels[offspringexcl] .= 0.0 
    newloglike = sum(loglikels)
    deltloglike = newloglike - inputloglike    
    isaccept = alwaysaccept ? true : deltloglike > 0.0
    if isaccept
        chrfhaplo .= newchrfhaplo
    end
    deltloglike, ndiff
end

function get_bdiff(oldchrfhaplo, newchrfhaplo, snpincl)    
    oldchrfhaplo2 = oldchrfhaplo[snpincl,:]
    chrfhaplo2 = newchrfhaplo[snpincl,:]
    perm = MagicBase.permfounder(oldchrfhaplo2, chrfhaplo2)
    if perm != 1:length(perm)
        oldchrfhaplo2 .= oldchrfhaplo2[:,perm]
    end		
    bdiff = map(isdiff_allele,oldchrfhaplo2, chrfhaplo2)    
    bdiff
end

function get_fmissls(fhaplosetpp::AbstractVector; isfounderinbred=true)
    nc = isfounderinbred ? 2 : 4 
    [begin 
        len = length.(fhaplosetpp[i])         
        mean(len .>= nc)
    end for i in eachindex(fhaplosetpp)]
end



function popfromfindex(findex::Integer,popmakeup::AbstractDict)
    popfromfindex([findex],popmakeup)
end

# list of all subpopulations whose parnets are contained in the vector findex
function popfromfindex(findex::AbstractVector,popmakeup::AbstractDict)
    res = []
    for (key, val) =popmakeup
        a = intersect(findex,val["founder"])
        isempty(a) || push!(res,key)
    end
    res
end

function calfwphaseprob(prob::AbstractVector)
    origscale = [sum.(p) for p = prob]
    # cal logpgeno
    logpgeno = [sum(log.(p)) for p = origscale] 
    b = isinf.(logpgeno)
    if all(b)
        nphase = length(logpgeno)
        phaseindex = collect(1:nphase)
        phaseprob = ones(nphase) ./ nphase
        origprob = prob
    else
        logpgeno .-= max(logpgeno...)
        phaseindex = findall(logpgeno .- max(logpgeno...) .> log(1e-20))
        phaseprob = normalize(exp.(logpgeno[phaseindex]),1)
        # cal origprob
        origprob = map((x1,x2)->x1./x2, prob[phaseindex],origscale[phaseindex])        
    end
    phaseindex,phaseprob,origprob
end

function founderforward(chroffgeno::AbstractMatrix, popidls::AbstractVector,
    offspringexcl::AbstractVector, 
    popmakeup::AbstractDict, priorprocess::AbstractDict,
    fhaploset::AbstractVector;     
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing, AbstractVector}, 
    baseerror::Union{Real,AbstractVector},
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},    
    allelicdropout::Union{Real,AbstractVector},
    snporder::AbstractVector,
    israndallele::Bool,
    issnpGT::AbstractVector, 
    forwardtempfile::AbstractString)
    nsnp = size(chroffgeno,1)    
    fwphaseprob = Vector{Vector{Float64}}(undef,nsnp)
    fwphaseindex = Vector{Vector{Int}}(undef,nsnp)
    startprob = getstartprob(popmakeup, priorprocess,popidls,offspringexcl)    
    if isempty(startprob)
        error("no non-execluded offspring")
    end
    epsfls = typeof(epsf) <: Real ? epsf*ones(nsnp) : epsf
    epsols = typeof(epso) <: Real ? epso*ones(nsnp) : epso
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nsnp) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nsnp) : allelicbias
    allelicoverdispersionls = typeof(allelicoverdispersion) <: Real ? allelicoverdispersion*ones(nsnp) : allelicoverdispersion
    allelicdropoutls = typeof(allelicdropout) <: Real ? allelicdropout*ones(nsnp) : allelicdropout
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)
    jldopen(forwardtempfile,"w") do file
        snp = snporder[tseq[1]]
        dataprobls = MagicReconstruct.calsitedataprob_multiphase(fhaploset[snp],chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind,baseerror=baseerrorls[snp], 
            allelicbias=allelicbiasls[snp], allelicoverdispersion=allelicoverdispersionls[snp],
            allelicdropout=allelicdropoutls[snp],                    
            popidls, offspringexcl, israndallele,issiteGT = issnpGT[snp])        
        noff =length(startprob)        
        origprob = [[startprob[i] .* dataprob[i] for i=1:noff] for dataprob in dataprobls]
        fwphaseindex[snp], fwphaseprob[snp],fworigprob_pre= calfwphaseprob(origprob)
        # file[string(snp)] = fworigprob_pre
        write(file, string(snp), fworigprob_pre)
        for kk in 2:length(tseq)            
            tranprob = gettranprob(tseq[kk-1], popmakeup, priorprocess,popidls,offspringexcl)
            snp = snporder[tseq[kk]]            
            dataprobls =MagicReconstruct.calsitedataprob_multiphase(fhaploset[snp],chroffgeno[snp,:],popmakeup;
                epsf=epsfls[snp], epso=epsols[snp], epso_perind,baseerror=baseerrorls[snp], 
                allelicbias=allelicbiasls[snp], allelicoverdispersion=allelicoverdispersionls[snp],
                allelicdropout=allelicdropoutls[snp], 
                popidls, offspringexcl, israndallele, issiteGT = issnpGT[snp])
            snppre = snporder[tseq[kk-1]]        
            # println("kk=",kk, ",snppre=",snppre,",fwphaseprob[snppre] =",fwphaseprob[snppre],
            #     "len_dataprobls = ", length(dataprobls))    
            prob = sum(fwphaseprob[snppre] .* fworigprob_pre)            
            prob2=[(prob[i]' * tranprob[i])[1,:] for i=1:noff]
            origprob =  [[prob2[i] .* dataprob[i] for i=1:noff] for dataprob in dataprobls]  
            fwphaseindex[snp], fwphaseprob[snp],fworigprob_pre= calfwphaseprob(origprob)
            write(file, string(snp), fworigprob_pre)
        end
    end
    fwphaseindex, fwphaseprob
end

function posmax(ls::AbstractVector)
    if all(isnan.(ls))
        rand(1:length(ls))
    else 
        v=max(ls...)
        rand(findall(ls .== v))
    end
end

function founderbackwardsample(fwphaseprob::AbstractVector,chroffgeno::AbstractMatrix,
    popidls::AbstractVector, offspringexcl::AbstractVector, popmakeup::AbstractDict,
    priorprocess::AbstractDict,
    snporder::AbstractVector,
    forwardtempfile::AbstractString, 
    isrand::Bool=false)
    nsnp= size(chroffgeno,1)    
    fphaseprob=Vector{Vector{Float64}}(undef,nsnp)
    fphase = zeros(Int, nsnp)
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)
    length(tseq)==1 && @warn string("tseq=",tseq)
    jldopen(forwardtempfile,"r") do file
        snp = snporder[tseq[end]]
        fworigprob = read(file,string(snp))        
        fphaseprob[snp]= normalize(fwphaseprob[snp],1)
        if isrand
            fphase[snp] = rand(Categorical(fphaseprob[snp]))
            orig = [rand(Categorical(i)) for i in fworigprob[fphase[snp]]]                    
        else
            fphase[snp] = posmax(fphaseprob[snp])
            orig = posmax.(fworigprob[fphase[snp]])
        end        
        noff = length(orig)
        for kk in (length(tseq)-1):-1:1        
            tranprob = gettranprob(tseq[kk],popmakeup, priorprocess,popidls,offspringexcl)            
            snp = snporder[tseq[kk]]
            fworigprob = read(file,string(snp))            
            ls = [tranprob[i][:,orig[i]] for i in 1:noff]
            ww = [sum(log.(map(dot, fw,ls))) for fw in fworigprob]
            ww .-= max(ww...)        
            fphaseprob[snp]=normalize(fwphaseprob[snp] .* exp.(ww),1)            
            # println("kk=",kk, ",ww=", ww,",fphaseprob[snp]=",fphaseprob[snp])            
            if isrand
                fphase[snp] = rand(Categorical(fphaseprob[snp]))
                orig = [rand(Categorical(i)) for i in fworigprob[fphase[snp]]]            
            else
                fphase[snp] = posmax(fphaseprob[snp])            
                orig = posmax.(fworigprob[fphase[snp]])
            end            
        end
    end    
    fphaseprob,fphase
end

function calofflogl(logbwprob::AbstractVector,chroffgeno::AbstractMatrix,
    popidls::AbstractVector, offspringexcl::AbstractVector, popmakeup::AbstractDict,
    priorprocess::AbstractDict,initfhaplo::AbstractVector, initsnp::Integer,
    israndallele::Bool,     
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector})
    offls = MagicReconstruct.offspringfrompop(popmakeup,popidls,offspringexcl)
    startprob = getstartprob(popmakeup, priorprocess,popidls)
    snp = initsnp
    initfhaplos = reshape(initfhaplo,1,:)
    epsfmarker = typeof(epsf) <: Real ? epsf : epsf[snp]
    epsomarker = typeof(epso) <: Real ? epso : epso[snp]
    dataprobls =MagicReconstruct.calsitedataprob_multiphase(initfhaplos,chroffgeno[snp,:],popmakeup;
        epsf=epsfmarker, epso=epsomarker, epso_perind, popidls, israndallele,issiteGT = issnpGT[snp])
    dataprob = sum(dataprobls) ./ length(dataprobls)
    fwprob = map((x,y)-> x .* y, startprob, dataprob)
    subofflogl = [begin
        bw = logbwprob[i]
        pmax = max(bw...)
        log(sum(exp.(bw .- pmax) .* fwprob[i])) + pmax
    end for i=eachindex(fwprob)]
    offls, subofflogl
end


# fhaplophase is modified  
function founderforwardbackward!(findex::AbstractVector,
    fhaplophase::AbstractMatrix,
    chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict, priorprocess::AbstractDict,
    fhaplosetpp::AbstractVector;          
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing, AbstractVector}, 
    baseerror::Union{Real,AbstractVector},
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},
    allelicdropout::Union{Real,AbstractVector},
    offspringexcl::AbstractVector,
    snporder::AbstractVector,
    israndallele::Bool,
    issnpGT::AbstractVector, 
    reversechr::Bool,    
    threshproposal::Real,
    forwardtempfile::AbstractString)
    if reversechr
        reverse!(snporder)
        for val=values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end        
    popidls = popfromfindex(findex,popmakeup)
    noffspring = sum([length(popmakeup[popid]["offspring"]) for popid in popidls])  
    # perform garbage collection to reduce memoryuse, while too often GG increases computational time
    isgc = rand() < noffspring/100 
    fhaploset= getfhaploset(findex,fhaplophase,fhaplosetpp)    
    fwphaseindex, fwphaseprob = founderforward(chroffgeno,
        popidls, offspringexcl, popmakeup,priorprocess,fhaploset; 
        epsf,epso,epso_perind, 
        baseerror,allelicbias,allelicoverdispersion,allelicdropout,
        snporder,israndallele,issnpGT,forwardtempfile)
    isgc && GC.gc()
    # fhaploset=map((x,y)->x[y,:],fhaploset,fwphaseindex)
    pri1=first(values(priorprocess))
    snps = snporder[pri1.markerincl]
    for i in snps
        fhaploset[i] = fhaploset[i][fwphaseindex[i],:]
    end     
    fphaseprob, fphase = founderbackwardsample(fwphaseprob,chroffgeno, popidls, offspringexcl, 
        popmakeup,priorprocess,snporder,forwardtempfile)
    isgc && GC.gc()    
    nfounder = length(fhaplosetpp)
    nfgl = size(fhaplophase,2)
    isfounderinbred = nfounder == nfgl    
    if threshproposal > 0.5 
        findex2 = isfounderinbred ? findex : reduce(vcat,[[2*i-1,2*i] for i in findex])
        for i in snps                            
            fhaplo = avg_sitefhaplo(fhaploset[i], fphaseprob[i],findex2; thresh = threshproposal)            
            # b = fhaplo .!= "N"  
            # fhaplophase[i,b] .= fhaplo[b]
            fhaplophase[i, :] .= fhaplo            
        end
    else    
        for i in snps
            fhaplophase[i,:] .= fhaploset[i][fphase[i],:]
        end            
    end            
    if reversechr
        reverse!(snporder)
        for val=values(priorprocess)
            MagicReconstruct.reverseprior!(val)
        end
    end
    fhaplophase
end

function avg_sitefhaplo(sitefhaplo::AbstractMatrix, sitephaseprob::AbstractVector,
    findex::AbstractVector; thresh = 0.9)
    nphase = length(sitephaseprob)
    nphase == size(sitefhaplo,1) || @error "dimension mismatch!"
    nphase == 1 && return sitefhaplo[1,:]
    fhaplo = replace(sitefhaplo[:,findex],"N"=>NaN,"1"=>0.0,"2"=>1.0)
    fhaplo2 = sum(fhaplo .* sitephaseprob,dims=1)[1,:]
    fhaplo3 = [i < 1-thresh ? "1" : (i>thresh ? "2" : "N") for i in fhaplo2]
    res = copy(sitefhaplo[1,:])
    res[findex] .= fhaplo3
    res
end


function getstartprob(popmakeup::AbstractDict, priorprocess::AbstractDict,
    popidls::AbstractVector,offspringexcl::AbstractVector)
    noff =sum([length(i["offspring"]) for i=values(popmakeup)])
    startprob = Vector{Vector{Float64}}(undef,noff)
    for popid = popidls
        hcode = popmakeup[popid]["hashcode"]
        offls = popmakeup[popid]["offspring"]
        initprob = priorprocess[hcode].initprob
        startprob[offls] = repeat([initprob],length(offls))
    end
    offls = MagicReconstruct.offspringfrompop(popmakeup,popidls,offspringexcl)    
    if issetequal(offls, 1:noff)
        startprob
    else
        startprob[offls]
    end
end

function gettranprob(t::Integer, popmakeup::AbstractDict, priorprocess::AbstractDict,
    popidls,offspringexcl::AbstractVector)
    noff =sum([length(i["offspring"]) for i=values(popmakeup)])
    tranprob = Vector(undef,noff)
    for popid = popidls
        hcode = popmakeup[popid]["hashcode"]
        offls = popmakeup[popid]["offspring"]
        prob = priorprocess[hcode].tranprobseq[t]        
        prob2 = isnothing(prob) ? I(size(priorprocess[hcode].tranrate,1)) : prob
        tranprob[offls] = repeat([prob2],length(offls))
    end
    offls = MagicReconstruct.offspringfrompop(popmakeup,popidls,offspringexcl)
    if issetequal(offls, 1:noff)
        tranprob
    else
        tranprob[offls]
    end
end