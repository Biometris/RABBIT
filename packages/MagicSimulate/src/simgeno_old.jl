
function genedrop(magicfgl::MagicGeno; isfounderinbred::Bool=true, isphased::Bool=true)
    markermap = deepcopy(magicfgl.markermap)
    foundergeno = deepcopy(magicfgl.foundergeno)
    nchr = length(foundergeno)
    offspringgeno = [begin
        fgeno = foundergeno[chr]
        if isfounderinbred
            fhaplo = permutedims(fgeno)
            markermap[chr][!,:founderformat] .= "GT_haplo"
        else
            fhaplo = vcat([hcat(col...) for col = eachcol(fgeno)]...)
            markermap[chr][!,:founderformat] .= "GT_phased"
        end        
        offfgl = permutedims(magicfgl.offspringgeno[chr])
        nsnp = size(offfgl,2)
        if isphased
            trueoff = [[fhaplo[i, snp] for i=offfgl[:, snp]] for snp=1:nsnp]
            markermap[chr][!,:offspringformat] .= "GT_phased"
        else
            trueoff = [[join(fhaplo[i,snp]) for i in offfgl[:,snp]] for snp = 1:nsnp]
            markermap[chr][!,:offspringformat] .= "GT_unphased"
        end
        Matrix{Any}(permutedims(reduce(hcat,trueoff)))
    end for chr=1:nchr]
    MagicGeno(magicfgl.magicped,markermap,foundergeno,offspringgeno,Dict())
end

# function applyerror!(chrgeno::AbstractMatrix,errls::AbstractVector)    
#     nsnp = size(chrgeno,1)
#     @assert nsnp == length(errls)
#     gunique = unique(chrgeno)
#     if gunique ⊆  ["1","2","N"]
#         wwls = [[1-err, err] for err in errls]                
#         gsetls = [["1","2","N"], ["2","1", "N"]]
#         if issetequal(["1","2"],gunique)
#             gsetls = [i[1:2] for i in gsetls]
#         end
#         grule = [Dict(gsetls[1] .=> gsetls[1]),Dict(gsetls[1] .=> gsetls[2])]
#     else
#         wwls = [[(1 - err)^2, err*(1 - err), (1 - err)*err, err^2] for err in errls]
#         gsetls = [["1","2","N", "11", "12", "21", "22", "1N","2N", "N1", "N2", "NN"],
#                   ["2","1","N", "21", "22", "11", "12", "2N","1N", "N1", "N2", "NN"],
#                   ["1","2","N", "12", "11", "22", "21", "1N","2N", "N2", "N1", "NN"],
#                   ["2","1","N", "22", "21", "12", "11", "2N","1N", "N2", "N1", "NN"]]
#         d = setdiff(gunique,gsetls[1])
#         isempty(d) || error(string("unknown genotypes ",d))
#         if gunique ⊆  ["11", "12", "21", "22", "1N","2N", "N1", "N2", "NN"]
#             if issetequal(["11", "12", "21", "22"],gunique)
#                 gsetls = [i[4:7] for i in gsetls]
#             else
#                 gsetls = [i[4:end] for i in gsetls]
#             end
#         end
#         grule = [Dict(gsetls[1] .=> gsetls[1]), Dict(gsetls[1] .=> gsetls[2]),
#                  Dict(gsetls[1] .=> gsetls[3]), Dict(gsetls[1] .=> gsetls[4])]
#     end
#     for snp in 1:nsnp
#         errww= StatsBase.ProbabilityWeights(wwls[snp])
#         chrgeno[snp,:] .=  [get(grule[sample(errww)],g,nothing) for g in chrgeno[snp,:]]
#     end
#     # in(nothing,chrgeno) && @error string("unknow genotypes: ",chrgeno[isnothing.(chrgeno)])
#     chrgeno    
# end


function trueg2obsg(trueg::AbstractString,err::Real,randallele::Real)
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
            if rand() < (1-err)^2
                # no err occurs
                trueg
            else
                if trueg == "11"                    
                    wsample(["12","22"], [1-err,err])
                elseif trueg == "12"                    
                    rand(["11","22"])
                else
                    #  trueg == "22"
                    wsample(["12","11"], [1-err,err])                
                end
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

function applyerror!(chrgeno::AbstractMatrix,errls::AbstractVector;
    chr_randallele::Union{Nothing,AbstractDict}=nothing)        
    gunique = unique(chrgeno)
    if in("21",gunique) 
        chrgeno[chrgeno .== "21"] .= "12"
    end
    if in("N1",gunique) 
        chrgeno[chrgeno .== "N1"] .= "1N"
    end
    if in("N2",gunique) 
        chrgeno[chrgeno .== "N2"] .= "2N"
    end
    setdiff!(gunique,["21","N1","N2"])
    if isnothing(chr_randallele)
        # scenario of outbred founders
        randallele = 1.0
        for snp in 1:size(chrgeno,1)
            chrgeno[snp,:] .= trueg2obsg.(chrgeno[snp,:],errls[snp],randallele)
        end
    else
        for popid in keys(chr_randallele)
            offls,randallele = chr_randallele[popid]            
            subgeno = view(chrgeno,:,offls)
            for snp in 1:size(subgeno,1)
                subgeno[snp,:] .= trueg2obsg.(subgeno[snp,:],errls[snp],randallele)
            end
        end
    end
    chrgeno
end

function applyerror!(geno::AbstractVector,errdist::Distribution;
    subpop_randallele::Union{Nothing,AbstractVector}=nothing)
    for chr in eachindex(geno)
        chrgeno = geno[chr]
        nsnp = size(chrgeno,1)
        errls = rand(errdist,nsnp)
        chr_randallele = isnothing(subpop_randallele) ? nothing : subpop_randallele[chr]
        applyerror!(chrgeno,errls; chr_randallele)
    end 
    geno
end

function applymiss!(chrgeno::AbstractMatrix,
    formatls::AbstractVector,
    missls::AbstractVector)
    nsnp,nind = size(chrgeno)
    nsnp == length(formatls) ==  length(missls) || @error "inconsistent number of markers"
    formatset = unique(formatls)
    for format in formatset
        snpls = findall(formatls .== format)              
        if format == "GT_haplo"
             for snp in snpls 
                ismiss = rand(Bernoulli(missls[snp]),nind)
                chrgeno[snp,ismiss] .= "N"
            end
        elseif format == "GT_unphased"
            for snp in snpls                
                ismiss = rand(Bernoulli(missls[snp]),nind)
                chrgeno[snp,ismiss] .= "NN"
            end
        elseif format == "AD"
            for snp in snpls          
                posls = findall(rand(Bernoulli(missls[snp]),nind))                
                for pos = posls                    
                    chrgeno[snp,pos] = [0,0]
                end
                # ismiss = rand(Bernoulli(missls[snp]),nind)
                # chrgeno[snp,ismiss] .= "NA"
            end
        else
            @error string("unkown format: ",format)
        end       
    end
    chrgeno
end

function applymiss!(genols::AbstractVector,formatlsls::AbstractVector,missdist::Distribution)
    misslsls = [rand(missdist,size(i,1)) for i in genols]    
    applymiss!.(genols,formatlsls,misslsls)
    genols
end


function randdepth(meandepth::Distribution,nsnp::Integer, nind::Integer)
    meanls = rand(meandepth,nsnp)
    hcat([rand(Poisson(i),nind) for i=meanls]...)'
end

function randreads(chrgeno::AbstractMatrix, meandepth::Distribution, seqerror::Real)
    # genotypes 1N, 2N, NN are not allowed
    guniq = unique(chrgeno)
    if issubset(guniq,["1","2"])
        pa1 = [1-seqerror, seqerror]
        genoset = ["1","2"]
    elseif issubset(guniq,["11","12","22"])
        pa1 = [1-seqerror, 0.5,  seqerror]
        genoset = ["11","12","22"]
    else
        error(string("unkown genotype: ", guniq))
    end
    dict = Dict(genoset .=> 1:length(genoset))
    indices = [get(dict, i, nothing) for i=chrgeno]
    if nothing in indices
        error(string("unknown genotypes: ", chrgeno[isnothing.(indices)]))
    end
    pp = pa1[indices]
    nn = randdepth(meandepth,size(chrgeno)...)
    reads = map((x,y)->[rand(Binomial(x,y)),x], nn, pp)
    reads = [[i[1], i[2]-i[1]] for i=reads]
    reads
end

function get_subpop_randallele(magicfgl::MagicGeno;error_randallele)    
    if isnothing(error_randallele) 
        pop2off = MagicBase.get_subpop2offspring(magicfgl.magicped;isindex=true)    
        [begin     
            isnonibd = allunique.(magicfgl.offspringgeno[chr])
            Dict([begin         
                offls = pop2off[popid]
                randallele = mean(isnonibd[:,offls])                 
                println("chr=",chr, ",popid=",popid,", randallele=",randallele)
                popid =>(offspring=offls, randallele=randallele)
            end for popid in keys(pop2off)])
        end for chr in eachindex(magicfgl.offspringgeno)]
    else
        noff = size(magicfgl.magicped.offspringinfo,1)
        [Dict(["subpop"=>(offspring=1:noff,randallele=error_randallele)]) 
            for chr in eachindex(magicfgl.offspringgeno)]
    end
end

function simgeno(truegeno::MagicGeno,magicfgl::MagicGeno;
    error_randallele::Union{Nothing,Real}=1.0,
    foundererror::Distribution=Beta(1,199),
    offspringerror::Distribution=Beta(1,199),
    foundermiss::Distribution=Beta(1,9),
    offspringmiss::Distribution=Beta(1,9),
    seqfrac::Real=0.5,
    seqdepth::Distribution = Gamma(2, 5),
    seqerror::Real = 0.001)        
    obsgeno = deepcopy(truegeno)
    # apply errors
    obsgeno.foundergeno .= [Matrix{Any}(join.(i))  for i in truegeno.foundergeno]
    obsgeno.offspringgeno .= [Matrix{Any}(join.(i))  for i in truegeno.offspringgeno]    
    applyerror!(obsgeno.foundergeno, foundererror)
    subpop_randallele = get_subpop_randallele(magicfgl;error_randallele)        
    applyerror!(obsgeno.offspringgeno, offspringerror;subpop_randallele)
    for chr in 1:length(obsgeno.markermap)        
        # sim sequencing
        isseq = rand(Bernoulli(seqfrac),size(obsgeno.markermap[chr],1))
        isnonseq = .!isseq
        fformat = [i=="GT_phased" ? "GT_unphased" : "GT_haplo" 
            for i in obsgeno.markermap[chr][isnonseq,:founderformat]]
        obsgeno.markermap[chr][isnonseq,:founderformat] .= fformat
        obsgeno.markermap[chr][isnonseq,:offspringformat] .= "GT_unphased"        
        if any(isseq)
            obsgeno.markermap[chr][isseq,:founderformat] .= "AD"
            obsgeno.markermap[chr][isseq,:offspringformat] .= "AD"
            alleles = obsgeno.foundergeno[chr][isseq,:]
            obsgeno.foundergeno[chr][isseq,:] .= randreads(alleles, seqdepth, seqerror) 
            alleles = obsgeno.offspringgeno[chr][isseq,:]
            obsgeno.offspringgeno[chr][isseq,:].= randreads(alleles, seqdepth, seqerror) 
        end        
    end
    formatlsls = [i[!,:founderformat] for i  in obsgeno.markermap]
    applymiss!(obsgeno.foundergeno, formatlsls, foundermiss)
    formatlsls = [i[!,:offspringformat] for i  in obsgeno.markermap]
    applymiss!(obsgeno.offspringgeno, formatlsls, offspringmiss)
    obsgeno
end
