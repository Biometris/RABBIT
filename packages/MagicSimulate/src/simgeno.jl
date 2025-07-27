
function simgeno(truegeno::MagicGeno,magicfgl::MagicGeno;
    foundermiss::Distribution=Beta(1,9),
    offspringmiss::Distribution=Beta(1,9),
    error_randallele::Union{Nothing,Real}=1.0,
    foundererror::Distribution=Beta(1,199),
    offspringerror::Distribution=Beta(1,199),    
    seqfrac::Real=0.5,
    baseerror::Distribution = Beta(1,199),
    allelicbias::Distribution = Beta(5,5),    
    allelicoverdispersion::Distribution = Exponential(0.05),    
    allelicdropout::Distribution = Beta(1,19),    
    seqdepth::Distribution = Gamma(2, 5),
    seqdepth_overdispersion::Distribution = Gamma(1,1))
    # initial obsgeno
    obsgeno = deepcopy(truegeno)
    obsgeno.foundergeno .= [Matrix{Any}(join.(i))  for i in truegeno.foundergeno]
    obsgeno.offspringgeno .= [Matrix{Any}(join.(i))  for i in truegeno.offspringgeno]  
    for m in obsgeno.markermap
        m[!,:founderformat] .= [i=="GT_phased" ? "GT_unphased" : i for i in m[!,:founderformat]]
        m[!,:offspringformat] .= [i=="GT_phased" ? "GT_unphased" : i for i in m[!,:offspringformat]]
    end
    MagicBase.setunphasedgeno!(obsgeno)
    # apply miss
    applymiss!(obsgeno; foundermiss, offspringmiss, targetfounder=true)
    applymiss!(obsgeno; foundermiss, offspringmiss, targetfounder=false)    
    # apply allelic error for SNP array, or read mis-alignment error for for sequencing data
    subpop_randallele = get_subpop_randallele(magicfgl;error_randallele)    
    applyerror!(obsgeno; foundererror,offspringerror,subpop_randallele, targetfounder=true)    
    applyerror!(obsgeno; foundererror,offspringerror,subpop_randallele, targetfounder=false)        
    # apply sequencing error and allelic balance due to unbalance amplification in PCR enrichment    
    simread!(obsgeno; seqfrac, baseerror, allelicbias, allelicoverdispersion, allelicdropout, seqdepth,seqdepth_overdispersion)
    missing2string!(obsgeno)
    # assume obsgeno and true have same markermap
    for chr in eachindex(obsgeno.markermap)        
        if all(obsgeno.markermap[chr][!,:marker] .== truegeno.markermap[chr][!,:marker])
            for col in [:foundererror,:offspringerror,:baseerror,:allelicbias,:allelicoverdispersion,:allelicdropout]
                truegeno.markermap[chr][!,col] .= obsgeno.markermap[chr][!,col] 
                obsgeno.markermap[chr][!,col]  .= missing
            end
        else
            @warn string("inconsistent markers between obsgeno and truegeno in chr=",chr)
        end
    end
    obsgeno
end

function applymiss!(magicgeno::MagicGeno; 
    foundermiss::Distribution, 
    offspringmiss::Distribution, 
    targetfounder::Bool=true)
    nchr = length(magicgeno.markermap)
    misslsls = Vector{Vector{Float64}}(undef, nchr)
    missdist = targetfounder ? foundermiss : offspringmiss
    for chr in 1:nchr 
        chrgeno = targetfounder ? magicgeno.foundergeno[chr] : magicgeno.offspringgeno[chr]
        nsnp,nind = size(chrgeno)
        missls = rand(missdist,nsnp)
        for snp in 1:nsnp
            ismiss = rand(Bernoulli(missls[snp]),nind)
            chrgeno[snp,ismiss] .= missing
        end
        misslsls[chr] = missls
    end
    misslsls
end


function applyerror!(magicgeno::MagicGeno;
    foundererror::Distribution,
    offspringerror::Distribution,    
    subpop_randallele::Union{Nothing,AbstractVector}=nothing,
    targetfounder::Bool=true)
    nchr = length(magicgeno.markermap)    
    errdist = targetfounder ? foundererror : offspringerror
    errcol = targetfounder ? :foundererror : :offspringerror
    for chr in 1:nchr
        chr_randallele = isnothing(subpop_randallele) ? nothing : subpop_randallele[chr]
        chrgeno = targetfounder ? magicgeno.foundergeno[chr] : magicgeno.offspringgeno[chr]                
        nsnp = size(chrgeno,1)
        b = ismissing.(magicgeno.markermap[chr][!,errcol])
        if any(b)
            errls = rand(errdist,nsnp)      
            magicgeno.markermap[chr][!,errcol] .= errls
        else
            errls = magicgeno.markermap[chr][!,errcol] 
        end
        if targetfounder || isnothing(chr_randallele)                 
            randallele = 1.0
            for snp in 1:nsnp
                chrgeno[snp,:] .= trueg2obsg.(chrgeno[snp,:],errls[snp],randallele)
            end
        else            
            for popid in keys(chr_randallele)
                offls,randallele = chr_randallele[popid]            
                subgeno =view(chrgeno,:,offls)
                for snp in 1:nsnp
                    subgeno[snp,:] .= trueg2obsg.(subgeno[snp,:],errls[snp],randallele)
                end
            end
        end        
    end
    magicgeno
end

function trueg2obsg(trueg::Union{Missing,AbstractString},err::Real,randallele::Real)
    ismissing(trueg) && return trueg
    if trueg in ["1","2","N"]
        if trueg == "1"
            wsample(["1","2"], [1-err,err])
        elseif trueg == "2"
            wsample(["2","1"], [1-err,err])
        else
            "N"
        end 
    elseif trueg in ["1N","2N","NN"]
        if trueg == "1N"
            wsample(["1N","2N"], [1-err,err])
        elseif trueg == "2N"
            wsample(["2N","1N"], [1-err,err])
        else
            "NN"
        end 
    else
        if !in(trueg,["11","12","22"])
            @error string("unknown trueg=",trueg)            
        end
        if rand()<randallele
            #rand allele error model
            if trueg == "11"                    
                wsample(["11", "12","22"], [(1-err)^2, 2*err*(1-err),err^2])
            elseif trueg == "22"                    
                wsample(["11", "12","22"], [err^2, 2*err*(1-err),(1-err)^2])
            else
                #  trueg == "12"
                wsample(["11", "12","22"], [err*(1-err),(1-err)^2+err^2, err*(1-err)])
            end       
        else
            #rand genotype error model
            if rand() < 1-err
                # no err occurs
                trueg
            else
                gset = ["11","12","22"]
                setdiff!(gset,[trueg])
                rand(gset)
                # println("trueg=",trueg,",obsg=",obsg)
                # obsg
            end
        end   
    end
end

function randdepth(seqdepth::Distribution,seqdepth_overdispersion::Distribution, 
    nind::Integer, nsnp::Integer)
    meanls = rand(seqdepth,nsnp)
    # NegativeBionomial(r,p), such that mean lam = r(1-p)/p, and variance = r(1-p)/p^2 
    # variance= r(1+lam/r) = r(1+seqdepth_overdispersion)        
    # reduce(hcat,[rand(Poisson(m),nind) for m in meanls]) #no overdispersion
    rpls = [begin 
        overdisperse = rand(seqdepth_overdispersion)
        r =m/overdisperse
        p = 1/(m/r + 1)
        (r,p)
    end for m in meanls]
    reduce(hcat,[rand(NegativeBinomial(r, p),nind) for (r,p) in rpls])
end

function  trueg2AD(trueg::Union{Missing,AbstractString},
    readdepth::Real, readerror::Real,
    readbalancemean::Real,
    readbalancedisperse::Real,
    allle_droupout::Real)
    ismissing(trueg) && return trueg
    readdepth == 0 && return zeros(Int,2)
    if trueg in ["1","11"]
        r1 = rand(Binomial(readdepth, 1-readerror))
    elseif trueg in ["2", "22"]
        r1 = rand(Binomial(readdepth, readerror))
    elseif trueg in ["12","21"]        
        if rand()< allle_droupout
            if rand() < 0.5                
                r1 = rand(Binomial(readdepth, 1-readerror)) # allele 2 is dropped out
            else
                r1 = rand(Binomial(readdepth, readerror)) # allele 1 is dropped out
            end
        else
            alpha = readbalancemean/readbalancedisperse
            beta = (1-readbalancemean)/readbalancedisperse
            prob_allele1 = rand(Beta(alpha,beta))
            r1 = rand(Binomial(readdepth, prob_allele1))
        end
    else
        @error string("unexpected true genotype: ", trueg)
    end
    [r1,readdepth-r1]
end

function simread!(magicgeno::MagicGeno; 
    seqfrac::Real=0.5,
    baseerror::Distribution = Beta(1,199),
    allelicbias::Distribution = Beta(5,5),    
    allelicoverdispersion::Distribution = Exponential(0.05),    
    allelicdropout::Distribution = Beta(1,19),    
    seqdepth::Distribution = Gamma(2, 5),
    seqdepth_overdispersion::Distribution = Gamma(1,1))
    seqfrac ≈ 0.0 && return magicgeno
    for chr in eachindex(magicgeno.markermap)
        isseq = rand(Bernoulli(seqfrac),size(magicgeno.markermap[chr],1))
        if any(isseq) 
            nmarker_seq = sum(isseq)
            col = :baseerror
            b = ismissing.(magicgeno.markermap[chr][isseq,col])
            if any(b)                    
                seqerrls = rand(baseerror, nmarker_seq)
                ls = Vector{Union{Missing, Float32}}(magicgeno.markermap[chr][:,col])
                ls[isseq] .= seqerrls
                magicgeno.markermap[chr][!,col] .= ls
            else
                seqerrls = magicgeno.markermap[chr][isseq,col]
            end

            col = :allelicbias
            b = ismissing.(magicgeno.markermap[chr][isseq,col])
            if any(b)                    
                allelicbiasls = rand(allelicbias, nmarker_seq)
                ls = Vector{Union{Missing, Float32}}(magicgeno.markermap[chr][:,col])
                ls[isseq] .= allelicbiasls
                magicgeno.markermap[chr][!,col] .= ls
            else
                allelicbiasls = magicgeno.markermap[chr][isseq,col]
            end

            col = :allelicoverdispersion
            b = ismissing.(magicgeno.markermap[chr][isseq,col])
            if any(b)                    
                allelicoverdispersionls = rand(allelicoverdispersion, nmarker_seq)
                ls = Vector{Union{Missing, Float32}}(magicgeno.markermap[chr][:,col])
                ls[isseq] .= allelicoverdispersionls
                magicgeno.markermap[chr][!,col] .= ls
            else
                allelicoverdispersionls = magicgeno.markermap[chr][isseq,col]
            end


            col = :allelicdropout
            b = ismissing.(magicgeno.markermap[chr][isseq,col])
            if any(b)                    
                allelicdropoutls = rand(allelicdropout, nmarker_seq)
                ls = Vector{Union{Missing, Float32}}(magicgeno.markermap[chr][:,col])
                ls[isseq] .= allelicdropoutls
                magicgeno.markermap[chr][!,col] .= ls
            else
                allelicdropoutls = magicgeno.markermap[chr][isseq,col]
            end

            for targetfounder in [true, false]   
                if targetfounder
                    chrtrue = permutedims(magicgeno.foundergeno[chr][isseq,:])
                else
                    chrtrue = permutedims(magicgeno.offspringgeno[chr][isseq,:])
                end
                nind = size(chrtrue,1)
                chrdepth = MagicSimulate.randdepth(seqdepth, seqdepth_overdispersion,nind,nmarker_seq)                
                for snp in 1:nmarker_seq
                    chrtrue[:,snp] .= trueg2AD.(chrtrue[:,snp],chrdepth[:,snp],seqerrls[snp],
                        allelicbiasls[snp],allelicoverdispersionls[snp],allelicdropoutls[snp])
                end
                if targetfounder
                    magicgeno.foundergeno[chr][isseq,:] .= permutedims(chrtrue)
                    magicgeno.markermap[chr][isseq,:founderformat] .= "AD"                    
                else
                    magicgeno.offspringgeno[chr][isseq,:].= permutedims(chrtrue)
                    magicgeno.markermap[chr][isseq,:offspringformat] .= "AD"
                end
            end            
        end
    end
    magicgeno
end

function get_subpop_randallele(magicfgl::MagicGeno;error_randallele)    
    if isnothing(error_randallele) 
        pop2off = MagicBase.get_subpop2offspring(magicfgl.magicped;isindex=true)    
        [begin     
            isnonibd = allunique.(magicfgl.offspringgeno[chr])
            Dict([begin         
                offls = pop2off[popid]
                randallele = mean(isnonibd[:,offls])                 
                # println("chr=",chr, ",popid=",popid,", randallele=",randallele)
                popid =>(offspring=offls, randallele=randallele)
            end for popid in keys(pop2off)])
        end for chr in eachindex(magicfgl.offspringgeno)]
    else        
        noff = size(magicfgl.magicped.offspringinfo,1)        
        [Dict(["pop"=>(offspring=1:noff, randallele=error_randallele)]) 
            for chr in eachindex(magicfgl.offspringgeno)]
    end
end

function missing2string!(magicgeno::MagicGeno)
    nchr = length(magicgeno.markermap)
    for targetfounder in [true,false]
        for chr in 1:nchr 
            if targetfounder
                chrgeno = magicgeno.foundergeno[chr]
                formatls = magicgeno.markermap[chr][!,:founderformat]
                subgeno = view(chrgeno,formatls .== "GT_haplo",:)
                subgeno[ismissing.(subgeno)] .= "N"
                subgeno = view(chrgeno,formatls .== "GT_unphased",:)
                subgeno[ismissing.(subgeno)] .= "NN"
            else
                chrgeno = magicgeno.offspringgeno[chr]
                formatls = magicgeno.markermap[chr][!,:offspringformat]                
                subgeno = view(chrgeno,formatls .== "GT_unphased",:)
                subgeno[ismissing.(subgeno)] .= "NN"
            end
        end
    end
    magicgeno
end