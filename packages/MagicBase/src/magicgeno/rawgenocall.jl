function rawgenocall!(magicgeno::MagicGeno;
    targets::AbstractVector= ["founders","offspring"],
    isfounderinbred::Bool=true,    
    isoffspringphased::Bool=false,
    callthreshold::Real=0.95,
    ishalfcall::Bool=true)    
    chridls = lowercase.([string(i[1,:linkagegroup]) for i=magicgeno.markermap])
    if !issubset(targets, ["founders","offspring"])
        error(string("targets must be a vector of founders or offspring"))
    end
    if in("founders",targets)        
        for chr in eachindex(chridls)
            format = String.(magicgeno.markermap[chr][!,:founderformat])
            snps = findall(format .== "GP")
            if !isempty(snps)
                magicgeno.foundergeno[chr] = Matrix{Any}(magicgeno.foundergeno[chr])
                chrgeno = magicgeno.foundergeno[chr]
                callgeno = callfromprob.(chrgeno[snps,:],callthreshold;
                    isphased=false,ishaplo=isfounderinbred,ishalfcall)
                if chridls[chr] in ["chrx"]
                    ismale = magicgeno.magicped.founderinfo[:gender] .== "male"
                    geno_malex = callfromprob.(chrgeno[snps,ismale],callthreshold;
                        isphased=false,ishaplo=true)
                    callgeno[snps,ismale]=[i*"Y" for i = geno_malex]
                end 
                chrgeno[snps,:] = callgeno            
                format[snps] .= isfounderinbred ? "GT_haplo" : "GT_unphased"
                magicgeno.markermap[chr][!,:founderformat] = format
            end
        end
    end
    if in("offspring",targets)
        for chr in eachindex(chridls)
            format = String.(magicgeno.markermap[chr][!,:offspringformat])
            snps = findall(format .== "GP")
            if !isempty(snps)
                magicgeno.offspringgeno[chr] = Matrix{Any}(magicgeno.offspringgeno[chr])
                chrgeno = magicgeno.offspringgeno[chr]
                calledgeno = callfromprob.(view(chrgeno,snps,:),callthreshold;
                    isphased=isoffspringphased,ishaplo=false,ishalfcall)                
                if chridls[chr] in ["chrx"]
                    ismale = magicgeno.magicped.offspringinfo[:gender] .== "male"
                    geno_malex = callfromprob.(chrgeno[snps,ismale],callthreshold;
                        isphased=isoffspringphased,ishaplo=true)
                    calledgeno[:,ismale]= [i*"Y" for i in geno_malex]
                end
                chrgeno[snps,:] = calledgeno
                format[snps] .= isoffspringphased ? "GT_phased" : "GT_unphased"
                magicgeno.markermap[chr][!,:offspringformat] = format
            end
        end
    end
    magicgeno
end

function rawgenoprob!(magicgeno::MagicGeno;
    targets::AbstractVector= ["founders","offspring"],
    baseerror::Real=0.001,    
    isfounderinbred::Bool=true,
    isoffspringinbred::Bool=false)
    # if isfounderinbred, readprob => genocall for founders
    errdigit0 = split(last(split(@sprintf("%.f",baseerror),".")),"")
    errdigit2=findfirst(x->x!="0",errdigit0)
    errdigit = isnothing(errdigit2) ? 6 : max(6,errdigit2+2)    
    chridls = lowercase.([string(i[1,:linkagegroup]) for i=magicgeno.markermap])
    if in("founders",targets) && !isnothing(magicgeno.foundergeno)
        for chr in eachindex(chridls)
            format = magicgeno.markermap[chr][!,:founderformat]
            snps = findall(format .== "AD")
            if !isempty(snps)
                magicgeno.foundergeno[chr] = Matrix{Any}(magicgeno.foundergeno[chr])
                chrgeno = magicgeno.foundergeno[chr]
                if isfounderinbred
                    prob = [round.(j,digits=errdigit) for j= genoprobhaplo.(chrgeno[snps,:],baseerror)]
                else
                    prob = [round.(j,digits=errdigit) for j= genoprobdiplo.(chrgeno[snps,:],baseerror)]
                end
                if chridls[chr] in ["chrx"]
                    ismale = magicgeno.magicped.founderinfo[:gender] .== "male"
                    prob[:,ismale] = [round.(j,digits=errdigit) for j= genoprobhaplo.(chrgeno[snps,ismale];baseerror)]
                end
                chrgeno[snps,:] .= prob
                format[snps] .= "GP"
            end
        end 
    end
    if in("offspring",targets) && !isnothing(magicgeno.offspringgeno)
        for chr in eachindex(chridls)
            format = magicgeno.markermap[chr][!,:offspringformat]
            snps = findall(format .== "AD")
            if !isempty(snps)
                magicgeno.offspringgeno[chr] = Matrix{Any}(magicgeno.offspringgeno[chr])
                chrgeno = magicgeno.offspringgeno[chr]
                if isoffspringinbred
                    prob = [round.(j,digits=errdigit) for j in genoprobhaplo.(chrgeno[snps,:],baseerror)]
                else
                    prob = [round.(j,digits=errdigit) for j in genoprobdiplo.(chrgeno[snps,:],baseerror)]
                end 
                if chridls[chr] in ["chrx"]
                    ismale = magicgeno.magicped.offspringinfo[:gender] .== "male"
                    prob[:,ismale] = [round.(j,digits=errdigit) for j in genoprobhaplo.(chrgeno[snps,ismale])]
                end
                chrgeno[snps,:] .= prob
                format[snps] .= "GP"
            end
        end
    end
    magicgeno
end

function genocalldiplo(genomtx::AbstractMatrix, formatls::AbstractVector;
    baseerror::Real=0.001,callthreshold::Real=0.95)
    # formatls[i] is the format of genomtx[i,:] for marker i
    calledgeno = copy(genomtx)
    formatset = unique(formatls)
    in("GT_haplo",formatset) && @error string("genocalldiplo does not work for G_haplo, use genocallhaplo instead")
    d = setdiff(formatset,["GT_phased","GT_unphased","AD","GP"])
    isempty(d) || @error string("unknow genotype format: ",d)
    ii = formatls .== "AD"
    subgeno = view(calledgeno,ii,:)
    subgeno .= genocalldiplo(subgeno; baseerror,callthreshold, isphased=false)
    ii = formatls .== "GP"
    subgeno = view(calledgeno,ii,:)
    subgeno .= callfromprob.(subgeno,callthreshold; isphased=false,ishaplo=false)
    if in("GT_phased", formatset)
        @warn string("GT_phased genotypes are transformed into GT_unphased")
        ii = formatls .== "GT_phased"
        subgeno = view(calledgeno,ii,:)
        subgeno .= join.(subgeno)
    end
    calledgeno
end

function genocalldiplo(allelicdepth::AbstractVecOrMat;
    baseerror::Real=0.001,callthreshold::Real=0.95,
    isphased::Bool=false)
    callfromprob.(genoprobdiplo.(allelicdepth, baseerror),callthreshold; isphased,ishaplo=false)
end


function genoprobdiplo_biallelic(reads::AbstractVector,baseerror::Real)
    # reads is a vector of integer
    length(reads)==2 || @error string("read counts for non bi-alleles: ",reads)    
    ploidy = 2
    p=[[dot([1-baseerror,baseerror], j) for j=[[1-i/ploidy,i/ploidy],[i/ploidy,1-i/ploidy]]] for i=0:ploidy]
    # logp matrix size: 2x(ploidy+1)
    # logp[i,j]: log probability of a sampled read being allele i given dosage being j-1
    logp = log.(hcat(p...))
    # ls[j]: log probability of reads given dosage being j-1
    ls0 = reads' * logp
    ls = ls0[1,:]
    w = logsumexp(ls)    
    res = exp.(ls .- w)
    # assume a uniform prior for 11, 12,and 22
    Vector{Float64}(res ./ sum(res))
end

function genoprobdiplo(reads::AbstractVector,baseerror::Real)
    na = length(reads)
    na >=2 || @error string("requires at least two read counts! read counts=",reads)    
    readdepth = sum(reads)
    phetero = 0.5*(1-baseerror+baseerror/(na-1))
    phetero2 = baseerror/(na-1)
    # 2*phetero+phetero2*(na-2) ≈ 1.0
    logpls = [if i == j
            ni = reads[i]
            ni*log(1-baseerror) + (readdepth-ni)*log(baseerror)
        else
            nij = reads[i] + reads[j]
            nij*log(phetero) + (readdepth-nij)*log(phetero2)    
        end for j in 1:na for i in 1:j]        
    wls = MagicBase.logsumexp(logpls)    
    res = exp.(logpls .- wls)
    # assume a equal prior probability for each unphased genotypes 11, 12,and 22
    round.(res, digits=10)
end

function genocallhaplo(genomtx::AbstractMatrix, formatls::AbstractVector;
    baseerror::Real=0.001,callthreshold::Real=0.95)    
    calledgeno = copy(genomtx)
    formatset = unique(formatls)
    if (in("GT_unphased", formatset) || in("GT_phased", formatset))
        @error string("genocallhaplo does not work for GT_unphased/GT_phased, use genocalldiplo instead")
    end
    d = setdiff(formatset,["GT_haplo","AD","GP"])
    isempty(d) || @error string("unexpected genotype format: ",d)
    ii = formatls .== "AD"
    subgeno = view(calledgeno,ii,:)
    subgeno .= genocallhaplo(subgeno; baseerror,callthreshold)
    ii = formatls .== "GP"
    subgeno = view(calledgeno,ii,:)
    subgeno .= callfromprob.(subgeno,callthreshold; ishaplo=true)
    calledgeno
end


function genocallhaplo(allelicdepth::AbstractVecOrMat;
    baseerror::Real=0.001,callthreshold::Real=0.95)
    callfromprob.(genoprobhaplo.(allelicdepth, baseerror),callthreshold; ishaplo=true)
end

function genoprobhaplo(reads::AbstractVector,baseerror::Real)
    # reads is a vector of integer
    length(reads)==2 || @error string("read counts for ≥3 alleles: ",reads)
    if reads[1]>0 && reads[2]==0
        [1.0,0.0]
    elseif reads[1]==0 && reads[2]>0
        [0.0, 1.0]
    else
        p = genoprobdiplo(reads, baseerror)
        p[2] >= 0.5 ? [0.5,0.5] : [p[1]+p[2]/2, p[3]+p[2]/2]
    end
end

function callfromprob(prob::AbstractVector, callthreshold::Real;
    isphased::Bool=false,ishaplo::Bool=false, ishalfcall::Bool=true)
    len = length(prob)
    len > 2 && ishaplo && @error string("wrong prob_vector for ishaplo=true")
    if len==2
        genostr = ishaplo ? ["1", "2"] : ["11","22"]
    elseif len==3
        genostr = ["11", "12", "22"]
    elseif len==4
        if isphased
            genostr =  [["1","1"], ["1","2"], ["2","1"],["2","2"]]
        else
            genostr = ["11", "12", "22"]
            prob = [prob[1],prob[2]+prob[3],prob[4]]
            len = 3
        end
    else
        error("callfromprob: wrong input prob length")
    end
    val, ii= findmax(prob)
    if val > callthreshold
        res = genostr[ii]
    else
         if ishalfcall
            if len==2
                res = ishaplo ? "N" : "NN"
            elseif len == 3
                p1n = prob[1] + prob[2]
                p2n = prob[3] + prob[2]
                if p1n>p2n && p1n > callthreshold
                    res = "1N"
                elseif p2n>p1n && p2n > callthreshold
                    res = "2N"
                else
                    res = "NN"
                end
            else   
                # len == 4             
                p1n = prob[1] + prob[2]+prob[3]
                p2n = prob[4] + prob[2]+prob[3]
                if p1n>p2n && p1n > callthreshold
                    if isphased
                        f1n = prob[1] + prob[2] 
                        fn1 = prob[1] + prob[3] 
                        if f1n > fn1 && f1n > callthreshold
                            res = ["1","N"]
                        elseif fn1 > f1n && fn1 > callthreshold
                            res = ["N","1"]
                        else
                            res = ["N","N"]
                        end
                    else
                        res = "1N"
                    end
                elseif p2n>p1n && p2n > callthreshold
                    if isphased
                        f2n = prob[3] + prob[4] 
                        fn2 = prob[2] + prob[4] 
                        if f2n > fn2 && f2n > callthreshold
                            res = ["2","N"]
                        elseif fn2 > f2n && fn2 > callthreshold
                            res = ["N","2"]
                        else
                            res = ["N","N"]
                        end                    
                    else
                        res = "2N"
                    end
                else
                    res = isphased ? ["N","N"] : "NN"
                end
            end
        else
            res = ishaplo ? "N" : (isphased && len==4 ? ["N","N"] : "NN")
        end
    end
    res
end
