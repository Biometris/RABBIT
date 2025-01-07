
function markerdelete_chr!(chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    delsiglevel::Real,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    offspringexcl::AbstractVector,
    snporder::AbstractVector,
    decodetempfile::AbstractString,    
    israndallele::Bool,     
    issnpGT::AbstractVector, 
    caldistance::Bool,priorlength::Real)
    deltt = markerdelete_chr(chrfhaplo,chroffgeno, popmakeup, priorprocess;
       delsiglevel,  epsf, epso, epso_perind, seqerror,allelebalancemean,
       allelebalancedisperse,alleledropout, offspringexcl, snporder, decodetempfile,
       issnpGT,israndallele, caldistance,priorlength)
    isempty(deltt) || MagicReconstruct.setpriorprocess!(priorprocess, deltt)
    deltt
end

function markerdelete_chr(chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    delsiglevel::Real,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    offspringexcl::AbstractVector,
    snporder::AbstractVector,
    decodetempfile::AbstractString,    
    israndallele::Bool=true,     
    issnpGT::AbstractVector, 
    caldistance::Bool,priorlength::Real)
    # calculate tseq and initialize res    
    pri1 = first(values(priorprocess))
    tseq = findall(pri1.markerincl)
    kksegls = MagicImpute.get_kksegls(pri1.markerdeltd[tseq];dtol=1e-6)
    length(kksegls) <=2 && (return [])
    res = zeros(length(kksegls))
    # calcualte binlogl and initialize parameters
    # logl for a bin of makrers is given at the first marker of the bin. 
    binlogl = MagicImpute.calbinlogl(chrfhaplo,chroffgeno;
        popmakeup,epsf, epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
        snporder, tseq, kksegls, israndallele,issnpGT)      
    # nfounder = size(chrfhaplo,2)    
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
    dataprobls = MagicReconstruct.init_dataprobls_singlephase(popmakeup)    
    jldopen(decodetempfile,"r") do bwfile
        snp = snporder[tseq[1]]
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:],chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
            allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp], 
            israndallele,issiteGT = issnpGT[snp])
        fw_prob_logl = MagicReconstruct.calinitforward(dataprobls, popmakeup)
        logbwprob = bwfile[string("t", tseq[1])]
        logl = calindlogl(fw_prob_logl,logbwprob, popmakeup; offspringexcl)
        sumlogl = sum(logl)
        for seg in eachindex(kksegls)
            kk,kkmax = kksegls[seg]
            if kk>=3                                
                pre_fw_kk = max(first(kksegls[seg-1])-1,1) # seg > 1 if kk>=3
                for i in pre_fw_kk:(kk-2)
                    snp = snporder[tseq[i+1]]
                    MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:],chroffgeno[snp,:],popmakeup;
                        epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
                        allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
                        alleledropout=alleledropoutls[snp],                     
                        israndallele, issiteGT = issnpGT[snp])
                    # input fw_prob_logl refers to tseq[i], output fw_prob_logl refers to tseq[i+1]
                    MagicReconstruct.calnextforward!(fw_prob_logl,tseq[i], dataprobls,popmakeup, priorprocess)
                end                    
            end
            if kk > 1 && rand() < 0.1    
                # && rand() < 0.1                
                sumlogl2 = sum(calindlogl(fw_prob_logl,bwfile[string("t", tseq[kk-1])], popmakeup; offspringexcl))
                if !isapprox(sumlogl,sumlogl2; atol=1e-2) 
                    @warn string("kk=",kk,", inconsistent [logl,logl2]=",[sumlogl,sumlogl2]) maxlog=20
                end
            end
            if kkmax < length(tseq)
                logbwprob = bwfile[string("t", tseq[kkmax+1])]
            end
            res[seg] = calvuongts!(dataprobls, kk,kkmax, tseq,logl,binlogl,fw_prob_logl, logbwprob,
                chrfhaplo, chroffgeno, snporder,popmakeup, priorprocess;
                epsfls, epsols, epso_perind,seqerrorls,allelebalancemeanls,allelebalancedispersels,
                alleledropoutls, offspringexcl, israndallele,issnpGT,caldistance,priorlength)
        end
        res
    end
    threshold = abs(quantile(Normal(),delsiglevel))
    res2 = res[.!isnan.(res)]    
    isempty(res2) || (threshold = max(threshold,quantile(res2, 0.9))) #delete at most 10%        
    kksegls2 = kksegls[res .> threshold] 
    # @info string("delting kkseg = ",kksegls2) 
    deltt = isempty(kksegls2) ? [] : reduce(vcat, [tseq[range(i...)] for i in kksegls2])    
    deltt
end

function get_kksegls(deltdls::AbstractVector; dtol=1e-6)
    ksegls = []
    k1 =1
    n = length(deltdls)
    for k in 1:n-1
        if deltdls[k] > dtol
            push!(ksegls,(k1,k))
            k1 = k+1    
        end
    end
    k1 <= n || @error string("unexpected k1=",k1, " in get_ksegls")
    push!(ksegls,(k1,n))    
    # issetequal(reduce(vcat, ksegls),1:n) || @error string("unexpected ksegls=",ksegls)
    ksegls
end

function calvuongts!(dataprobls, kk::Integer,kkmax::Integer, tseq::AbstractVector,logl::AbstractVector,
    binlogl::AbstractMatrix,
    fw_prob_logl::AbstractVector,
    logbwprob::AbstractVector,
    chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    snporder::AbstractVector,popmakeup::AbstractDict,priorprocess::AbstractDict;
    epsfls::AbstractVector,epsols::AbstractVector,
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerrorls::AbstractVector,
    allelebalancemeanls::AbstractVector,
    allelebalancedispersels::AbstractVector,
    alleledropoutls::AbstractVector,
    offspringexcl::AbstractVector,
    israndallele::Bool,     
    issnpGT::AbstractVector, 
    caldistance::Bool, priorlength::Real)
    # logl for a bin of makrers is given at the first marker of the bin. 
    isologl = binlogl[snporder[tseq[kk]],:]
    b=.!isnan.(isologl)
    any(b) || return 0.0
    # length(kksegls )>=2 and thus kk != length(tseq)
    if kkmax==length(tseq)
        proplogl0 = fw_prob_logl[2]
        proplogl = proplogl0[b] .+ isologl[b]
    elseif kk==1
        snp = snporder[tseq[kkmax+1]]
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
            chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
            allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp], 
            israndallele,issiteGT = issnpGT[snp])
        prop_prob_logl = MagicReconstruct.calinitforward(dataprobls, popmakeup)
        proplogl0 = calindlogl(prop_prob_logl,logbwprob, popmakeup; offspringexcl)
        proplogl = proplogl0[b] .+ isologl[b]    
    else
        snp = snporder[tseq[kkmax+1]]
        if caldistance
            MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
                chroffgeno[snp,:],popmakeup;
                epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
                allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
                alleledropout=alleledropoutls[snp], 
                israndallele, issiteGT = issnpGT[snp])
            ttrandis = calinterdis(tseq[kk-1],popmakeup,priorprocess,offspringexcl, 
                dataprobls,fw_prob_logl[1],logbwprob,priorlength)
        else
            pri1=first(values(priorprocess))
            ttrandis=sum(pri1.markerdeltd[tseq[kk-1:kkmax]])
        end
        MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
            chroffgeno[snp,:],popmakeup;
            epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
            allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
            alleledropout=alleledropoutls[snp], 
            israndallele, issiteGT = issnpGT[snp])
        prop_prob_logl = deepcopy(fw_prob_logl) # fw_prob_logl refers to tseq[kk-1]
        MagicReconstruct.calnextforward!(prop_prob_logl, tseq[kk-1],dataprobls,popmakeup, priorprocess; ttrandis)
        proplogl0 = calindlogl(prop_prob_logl,logbwprob,popmakeup; offspringexcl)
        proplogl = proplogl0[b] .+ isologl[b]
    end
    logldiff = proplogl .- logl[b]
    dfdiff = length(unique(chrfhaplo[snporder[tseq[kk:kkmax]],:]))==1 ? 0 : -1
    vuongts = calvuongts(logldiff,dfdiff)
    # snp = snporder[tseq[kk]]
    # println("kk=", kk, "; vuongts=", vuongts, "; epsf=", epsfls[snp], "; epso=", epsols[snp], ";logldiff_var=", var(logldiff))
    vuongts
end

function calvuongts(logldiff::AbstractVector,dfdiff::Integer)
    n = length(logldiff)
    v = var(logldiff)
    vuong = sum(logldiff)
    # correction for the difference of number of parameters
    vuong -= log(n)*(dfdiff)/2
    vuong /= sqrt(n*v)
    vuong
end

function calindlogl(fw_prob_logl::AbstractVector,logbwprob::AbstractVector,
    popmakeup::AbstractDict; offspringexcl::AbstractVector)
    fwprob, fwlogl = fw_prob_logl
    noff = length(fwlogl)
    res = zeros(noff)    
    for popid in keys(popmakeup)        
        offls = popmakeup[popid]["offspring"]        
        offls2 = isempty(offspringexcl)  ? offls : setdiff(offls, offspringexcl)        
        for off in offls2
            lpmax = maximum(logbwprob[off])            
            lp = exp.(logbwprob[off] .- lpmax)             
            res[off] = log(dot(lp, fwprob[off])) + lpmax + fwlogl[off]
        end
    end
    res
end

function calindlogl(fw_prob_logl::AbstractVector,decodetempfile::AbstractString,
    popmakeup::AbstractDict,t::Integer)
    logbwprob = JLD2.load(decodetempfile,string("t", t))    
    calindlogl(fw_prob_logl,logbwprob,popmakeup; offspringexcl=[])
end

# chroffgeno[snp, off]: offspring genotype at marker index snp
# chrfhaplo[snp, off]: founder genotype at marker index snp
# logl for a bin of makrers is given at the first marker of the bin. 
function calbinlogl(chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix;    
    popmakeup::AbstractDict,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    snporder::AbstractVector, 
    tseq::AbstractVector,    
    kksegls::AbstractVector,     
    israndallele::Bool,
    issnpGT::AbstractVector)
    isoffphased = false    
    fderive, offcode = MagicReconstruct.precompute_chr(chrfhaplo,
        chroffgeno,popmakeup, isoffphased, issnpGT)
    nsnp, noff = size(chroffgeno)
    binlogl = Matrix{Union{Missing,Float64}}(missing,nsnp,noff)
    isnothing(epso_perind) && (epso_perind = zeros(noff))
    for popid in keys(popmakeup)
        ishaploid = popmakeup[popid]["ishaploid"]
        offls = popmakeup[popid]["offspring"]
        nzstate = popmakeup[popid]["nzstate"]
        initprob  =  popmakeup[popid]["initprob"]
        dataprobseq = [zeros(MagicReconstruct._float_like,length(nzstate))  for _ in 1:nsnp]
        for off = offls
            obsseq = offcode[:,off]                     
            MagicReconstruct.caldataprobseq!(dataprobseq, obsseq,epsf,epso, epso_perind[off], seqerror,
                allelebalancemean,allelebalancedisperse,alleledropout,fderive,nzstate, 
                isoffphased,israndallele,issnpGT,ishaploid)
            for (kk,kkmax) in kksegls
                snpls = snporder[tseq[kk:kkmax]]                
                dataprob = [reduce(.*, i) for i in eachrow(reduce(hcat,dataprobseq[snpls]))]
                binlogl[first(snpls), off] = log(dot(dataprob,initprob))
            end
        end
    end
    binlogl
end



# chroffgeno[snp, off]: offspring genotype at marker index snp
# chrfhaplo[snp, off]: founder genotype at marker index snp
function calsinglelogl(chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix;    
    snpincl= 1:size(chroffgeno,1), 
    popmakeup::AbstractDict,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    israndallele::Bool,
    issnpGT::AbstractVector)
    isoffphased = false        
    if issetequal(snpincl, 1:size(chroffgeno,1)) 
        chroffgeno2 = chroffgeno
        chrfhaplo2 = chrfhaplo        
        issnpGT2 = issnpGT
    else
        chroffgeno2 = view(chroffgeno,snpincl,:)
        chrfhaplo2 = view(chrfhaplo,snpincl,:)                        
        issnpGT2 = issnpGT[snpincl]        
        ls =[isa(x, AbstractVector) ? x[snpincl] : x for x in [epsf, epso,seqerror, allelebalancemean,allelebalancedisperse,alleledropout]]        
        epsf, epso,seqerror, allelebalancemean,allelebalancedisperse,alleledropout = ls
    end    
    fderive, offcode = MagicReconstruct.precompute_chr(chrfhaplo2, chroffgeno2,popmakeup, isoffphased, issnpGT2)
    nsnp, noff = size(chroffgeno2)    
    isnothing(epso_perind) && (epso_perind = zeros(noff))
    singlelogl = Matrix{Union{Missing,Float64}}(missing,nsnp,noff)    
    for popid in keys(popmakeup)
        ishaploid = popmakeup[popid]["ishaploid"]
        offls = popmakeup[popid]["offspring"]
        nzstate = popmakeup[popid]["nzstate"]
        initprob  =  popmakeup[popid]["initprob"]
        dataprobseq = [zeros(MagicReconstruct._float_like,length(nzstate))  for _ in 1:nsnp]
        for off = offls
            obsseq = offcode[:,off]                      
            MagicReconstruct.caldataprobseq!(dataprobseq,obsseq,epsf,epso,epso_perind[off], seqerror,
                allelebalancemean,allelebalancedisperse,alleledropout,fderive,nzstate, 
                isoffphased,israndallele,issnpGT2,ishaploid)
            singlelogl[:,off] .= [log(dot(i,initprob)) for i in dataprobseq]
        end
    end    
    singlelogl
end
