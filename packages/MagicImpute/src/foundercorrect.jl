
function foundercorrect_chr!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,    
    popmakeup::AbstractDict,priorprocess::AbstractDict, priorspace::AbstractDict, fhaplosetpp::AbstractVector;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    snporder::AbstractVector,    
    isallowmissing::Bool, 
    decodetempfile::AbstractString,    
    founder2progeny::AbstractVector,
    ismalexls::AbstractVector,    
    israndallele::Bool,     
    offspringexcl::AbstractVector, 
    issnpGT::AbstractVector, 
    itmax::Integer=3)  
    isfounderinbred = length(founder2progeny) == size(chrfhaplo,2)	    
    rescorrect = []
    for it in 1:itmax        
        calledchroffgeno = singlesite_genocall(chrfhaplo,chroffgeno; ismalexls, popmakeup,
            epsf,epso, epso_perind, seqerror, allelebalancemean, allelebalancedisperse,alleledropout, 
            israndallele,issnpGT,callthreshold = 0.7,tempjld2file = decodetempfile)
        correctdf = last(foundererror_chr(chrfhaplo, chroffgeno, calledchroffgeno, 
            founder2progeny, popmakeup,priorprocess,priorspace;
            epsf,epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
            snporder, decodetempfile,offspringexcl, issnpGT,israndallele))
        if !isempty(correctdf)
            setcorrectdf!(chrfhaplo,correctdf,fhaplosetpp; isfounderinbred,isallowmissing)
            push!(rescorrect,correctdf)
        end
        isempty(correctdf) && break
    end
    res = isempty(rescorrect) ? rescorrect : reduce(vcat, rescorrect)
    res
end

function foundererror_chr(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,    
    calledchroffgeno::AbstractMatrix,        
    founder2progeny::AbstractVector,
    popmakeup::AbstractDict,priorprocess::AbstractDict, priorspace::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    snporder::AbstractVector,    
    decodetempfile::AbstractString,    
    israndallele::Bool,
    offspringexcl::AbstractVector, 
    issnpGT::AbstractVector)
    loglike = first(MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf,epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
        hmmalg="viterbi", decodetempfile, israndallele, issnpGT,snporder));
    # calculate bestgeno from viterbi path    
    nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)	
    if nstate == nfgl
        ancestralstate = [[i,i] for i in priorspace["haploindex"]]
    else
        nstate == nfgl^2 || @error "inconsistent number of states"
        ancestralstate = MagicBase.prior_diploindex(nfgl)
    end 
    chrviterbi = MagicReconstruct.get_chr_viterbi(decodetempfile, snporder)
    bestgeno = get_chr_bestgeno(chrfhaplo, chrviterbi,ancestralstate)
    chrbadsnp = getchrbadsnp(calledchroffgeno, bestgeno, popmakeup; minnerr = 2)
    # max 50% markers are labled as chrbadsnp
    wwls = [max(i[2]...) for i in chrbadsnp]
    q0=1.0-0.50*length(snporder)/length(wwls)
    if q0 > 0
        q = quantile(wwls, q0)
        chrbadsnp = chrbadsnp[wwls.>q]
    end
    nselsnp = length(chrbadsnp)
    correctdf = getchrchange(chrfhaplo,calledchroffgeno,chrviterbi,chrbadsnp, 
        decodetempfile, ancestralstate,founder2progeny; offspringexcl, minnerrdiff = 2)
    loglike, nselsnp, correctdf    
end


function singlesite_genocall(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix;
    popmakeup::AbstractDict,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    ismalexls::AbstractVector,
    israndallele::Bool, 
    issnpGT::AbstractVector,    
    callthreshold::Real,
    tempjld2file::AbstractString,
    )
    singlesite_genoprob(chrfhaplo,chroffgeno;  popmakeup,
        epsf,epso, epso_perind, seqerror, allelebalancemean, allelebalancedisperse,
        alleledropout, israndallele,issnpGT,outjld2file= tempjld2file)
    res = copy(chroffgeno)    
    singlesite_callfromprob!(res,tempjld2file;ismalexls,issnpGT,callthreshold)
    res
end

function singlesite_callfromprob!(res, probjld2file::AbstractString; ismalexls::AbstractVector, issnpGT::AbstractVector, callthreshold::Real)
    jldopen(probjld2file,"r") do file                
        noffspring = file["noffspring"]         
        nonissnpGT = .!issnpGT
        for off in 1:noffspring
            probls = file[string(off)]
            res[nonissnpGT,off] .= [MagicBase.callfromprob(i,callthreshold; isphased=false) for i in probls[nonissnpGT]]
        end
         # gintersect = intersect(["21","N1","N2"],unique(calledchroffgeno))
        # isempty(gintersect) || @error string("unexpected called offspring genotypes: ",gintersect)
        size(res,2)==length(ismalexls) || @error "inconsistent size" 
        if any(ismalexls)
            # TODO: to check the consistency of calledgeno for male x
            # male: N, 1, 2 => NN, 11, 22
            # male: NY, 1Y, 2Y => NN, 11, 22
            res[:,ismalexls]=[i[1]*i[1] for i= res[:,ismalexls]]
        end
        res
    end
end



function singlesite_genoprob(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix;    
    popmakeup::AbstractDict,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    israndallele::Bool,     
    issnpGT::AbstractVector, 
    outjld2file::AbstractString)
    nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    
    nsnp,nfgl2 = size(chrfhaplo)
    nfgl == nfgl2 || @error "inconsistent nfgl"
    f(x) = isa(x,AbstractVector) ? x : x*ones(nsnp)
    epsfls,epsols, seqerrorls = f(epsf), f(epso), f(seqerror)
    allelebalancemeanls, allelebalancedispersels, alleledropoutls  = f(allelebalancemean), f(allelebalancedisperse), f(alleledropout)
    isepsf_equal = allequal(epsfls)	    
    jldopen(outjld2file,"w") do file
        write(file,"nmarker",nsnp)
        write(file,"noffspring",size(chroffgeno,2))
        for popid in keys(popmakeup)        
            offls = popmakeup[popid]["offspring"]
            nzstate =  popmakeup[popid]["nzstate"]
            ishaploid =  popmakeup[popid]["ishaploid"]		
            initprob = spzeros(nstate)
            initprob[nzstate] .= popmakeup[popid]["initprob"]      
            res = Matrix(undef, nsnp, length(offls))  
            if ishaploid                       
                fderivehaplo = MagicReconstruct.calfderive(chrfhaplo; nzstate,ishaploid)		
                if isepsf_equal
                    prior = MagicReconstruct.haploprior(first(epsfls))
                end			
                dict = Dict("11"=>[1.0,0.0],"22"=>[0.0,1.0],"1N"=>[1.0,0.0],"2N"=>[0.0,1.0],"12"=>[0.0,0.0],"NN"=>[0.5,0.5])
                for snp in 1:nsnp                    
                    if issnpGT[snp]
                        sitegeno = view(chroffgeno, snp, offls)
                        for i in eachindex(offls)
                            res[snp,i] = dict[sitegeno[i]]
                        end
                    else
                        ismiss = [(i==[0,0] || i==[0.5,0.5]) for i in view(chroffgeno, snp, offls)]
                        isnonmiss  = .!ismiss
                        offls2 = offls[isnonmiss]
                        sitegeno = view(chroffgeno, snp, offls2)
                        if !isempty(sitegeno)                                                        
                            if isnothing(epso_perind) || all(epso_perind .< 2e-5)               
                                epso2 = epsols[snp] 
                            else
                                epso2 = [i + epsols[snp]  - i * epsols[snp] for i in epso_perind[offls2]]       
                            end
                            if eltype(first(sitegeno)) <: Integer
                                # format = AD
                                like = MagicReconstruct.haplolikeGBS(sitegeno,epso2,seqerrorls[snp])
                            else
                                # format = GP
                                like = MagicReconstruct.haplolikeGBS(sitegeno,epso2)
                            end
                            if !isepsf_equal
                                prior = MagicReconstruct.haploprior(epsfls[snp])
                            end
                            ls = prior[:,fderivehaplo[snp,nzstate]] * initprob[nzstate]
                            post = ls' .* like
                            res[snp,isnonmiss] .= [round.(normalize(i,1),digits=6) for i in eachrow(post)]   
                        end
                        res[snp,ismiss] .= [[0.5,0.5] for _ in 1:sum(ismiss)]
                    end                                     
                end
            else
                nstate == nfgl^2 || error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))					
                zzstate = setdiff(1:nstate, nzstate)					
                ibdbool= falses(nfgl^2)
                ibdbool[[(i-1)nfgl+i for i=1:nfgl]] .= true
                nonibdbool = .!ibdbool
                ibdbool[zzstate] .= false
                nonibdbool[zzstate] .= false		
                initprob_ibd = initprob[ibdbool]
                initprob_nonibd = initprob[nonibdbool]
                fderivediplo = MagicReconstruct.calfderive(chrfhaplo; nzstate,ishaploid)		            
                if isepsf_equal
                    prior = MagicReconstruct.diploprior(first(epsfls))
                end            
                dict = Dict("11"=>[1.0,0.0,0.0],"22"=>[0.0,0.0,1.0],"12"=>[0.0,1.0,0.0],
                    "1N"=>[0.333,0.666,0.0],"2N"=>[0.0,0.666,0.333],"NN"=>[0.25,0.5,0.25])                 
                for snp in 1:nsnp
                    if issnpGT[snp]                        
                        sitegeno = view(chroffgeno, snp, offls)
                        for i in eachindex(offls)
                            res[snp,i] = dict[sitegeno[i]]
                        end
                    else                        
                        ismiss = [(i==[0,0] || i==[0.25,0.5,0.25]) for i in view(chroffgeno, snp, offls)]
                        isnonmiss  = .!ismiss                  
                        offls2 = offls[isnonmiss]      
                        sitegeno = view(chroffgeno, snp, offls2)
                        if !isempty(sitegeno)                 
                            if isnothing(epso_perind) || all(epso_perind .< 2e-5)               
                                epso2 = epsols[snp] 
                            else
                                epso2 = [i + epsols[snp]  - i * epsols[snp] for i in epso_perind[offls2]]       
                            end
                            if eltype(first(sitegeno)) <: Integer
                                # format = AD
                                like = MagicReconstruct.diplolikeGBS(sitegeno,epso2,seqerrorls[snp],
                                    allelebalancemeanls[snp],allelebalancedispersels[snp],alleledropoutls[snp]; israndallele)
                            else
                                # format = GP
                                like = MagicReconstruct.diplolikeGBS(sitegeno,epso2;israndallele)
                            end
                            if !isepsf_equal
                                prior = MagicReconstruct.diploprior(epsfls[snp])
                            end
                            ls = prior.nonibd[:,fderivediplo[snp,nonibdbool]] * initprob_nonibd
                            ls .+= prior.ibd[:,fderivediplo[snp,ibdbool]] * initprob_ibd 
                            post = ls' .* like
                            res[snp,isnonmiss] .= [begin 
                                v = normalize(i,1)
                                round.([v[1],v[2]+v[3],v[4]],digits=6) 
                            end for i in eachrow(post)]
                        end
                        res[snp,ismiss] .= [[0.25,0.5,0.25] for _ in 1:sum(ismiss)]
                    end
                end
            end            
            for i in eachindex(offls)
                write(file, string(offls[i]), collect(res[:,i]))
            end            
        end
        outjld2file
    end
end


function getismalexls(magicped::MagicPed,chrid::AbstractString)
    isautosome=!(chrid in ["chrx"])
    if isautosome
        noff = size(magicped.offspringinfo,1)
        ismalexls = falses(noff)
    else
        ismalexls = magicped.offspringinfo[!, :gender] .== "male"
    end
end

function setchrpostprob!(chrpostprob::AbstractVector,delsnps::AbstractVector)
    isempty(delsnps) && return nothing
    nsnp = size(chrpostprob[1],1)
    isnondel = trues(nsnp)
    isnondel[delsnps] .= false
    for i in eachindex(chrpostprob)
        chrpostprob[i]=chrpostprob[i][isnondel,:]
    end
end

function setcorrectdf!(chrfhaplo::AbstractMatrix, correctdf::DataFrame, fhaplosetpp::AbstractVector; isfounderinbred::Bool,isallowmissing)
    for i=1:size(correctdf,1)
        snp, p, newg = correctdf[i,[1,2,4]]        
        if isfounderinbred        
            chrfhaplo[snp,p]=newg
            if isallowmissing
                fhaplosetpp[p][snp] = ["1", "2","N"]            
            else
                fhaplosetpp[p][snp] = ["1", "2"]
            end
            
        else
            chrfhaplo[snp,[2*p-1, 2*p]] .= newg            
            if isallowmissing
                fhaplosetpp[p][snp] = [["1","1"],["1","2"],["2","1"],["2","2"],["1","N"],["2","N"],["N","1"],["N","2"],["N","N"]]        
            else
                fhaplosetpp[p][snp] = [["1", "1"], ["1", "2"], ["2", "1"],["2", "2"]]
            end            
        end
    end
end

function transferchangedf(correctdf::DataFrame,chrmarkermap::DataFrame,
    founderinfo::DataFrame)
    if isempty(correctdf)
        res=[]
    else
        res = [begin
            # [snp, p, oldfg,newfg,oldnerr,newnerr,nnonmiss]
            snp, p,oldg,newg = correctdf[i,1:4]
            snpid = chrmarkermap[snp,:marker]
            chrid = chrmarkermap[snp,:linkagegroup]
            pos =  chrmarkermap[snp,:poscm]
            pid = founderinfo[p,:individual]
            newg2 = join(newg,"|")
            oldg2 = join(oldg,"|")
            vcat([chrid,pos,snpid,pid,oldg2,newg2],Vector(correctdf[i,5:7]))
        end for i=1:size(correctdf,1)]
    end
    pushfirst!(res,["linkagegroup", "position", "marker", "parent",
        "old_genotype", "new_genotype", "old_nerr","new_nerr","new_nmiss"])
    res2=permutedims(hcat(res...))
    df =DataFrame(res2[2:end,:],Symbol.(res2[1,:]))
    for i=[1,3,4,5,6]
        df[!,i]=String.(df[!,i])
    end
    df[!,2]=Float64.(df[!,2])
    for i=[7,8,9]
        df[!,i]=Integer.(df[!,i])
    end
    df
end

function getchrchange(chrfhaplo::AbstractMatrix, calledchroffgeno::AbstractMatrix, 
    beststate::AbstractMatrix,
    chrbadsnp::AbstractVector, decodetempfile::AbstractString, 
    ancestralstate::AbstractVector,founder2progeny::AbstractVector;    
    offspringexcl::AbstractVector, 
    minnerrdiff::Integer=3)
    res=[]    
    isfounderinbred = length(founder2progeny) == size(chrfhaplo,2)	
	jldopen(decodetempfile, "r") do file
	    for (snp,ww) = chrbadsnp
            if isfounderinbred
                isnonmiss = chrfhaplo[snp,:] .!= "N"
            else
                isnonmiss = [i != ["N","N"] for i in partition(chrfhaplo[snp,:],2)]
            end
            pp = findall((ww .>0) .&& isnonmiss)
	        change = [begin            
                if isfounderinbred
                    genoset = ["1","2"]
                    oldfg=chrfhaplo[snp,p]              
                else
                    genoset  = [[i,j] for i=["1","2"] for j=["1","2"]]
                    oldfg=chrfhaplo[snp,[2p-1,2p]]          
                end
                deleteat!(genoset,[i == oldfg for i in genoset])
                pushfirst!(genoset,oldfg)
                snps = snp*ones(Int,length(genoset))
                newfhaplo=Matrix{Any}(chrfhaplo[snps,:])
                if isfounderinbred
                    newfhaplo[:,p] .= genoset
                else
                    newfhaplo[:,[2*p-1,2*p]] .= permutedims(reduce(hcat,genoset))
                end
                progeny = founder2progeny[p]                
                if !isempty(offspringexcl)
                    progeny = setdiff(progeny, offspringexcl)                
                end
                newbeststate = beststate[snps,progeny]
                newbestgeno = best_state2geno(newfhaplo,newbeststate,ancestralstate)
                newoffgeno = calledchroffgeno[snps,progeny]
                mismatch,nonmissing = calgenomismatch(newoffgeno,newbestgeno)
                nerr = sum(mismatch,dims=2)[:,1]
                nnonmiss = sum(nonmissing,dims=2)[:,1]
                newnerr,index = findmin(nerr)
                newfg = genoset[index]
                oldnerr = nerr[1]
                [snp, p, oldfg,newfg,oldnerr,newnerr,nnonmiss[index]]
	        end for p in pp]
			change2=permutedims(hcat(change...))
	        val,ii=findmax(change2[:,5]-change2[:,6])
	        # reduced errors by an amound of at least minnerrdiff
	        if val >= minnerrdiff
	            push!(res,change2[ii,:])
	        end
	    end
    end
    colid = Symbol.(["marker","parent","oldg", "newg", "olderr","newerr", "newnonmiss"])
    if isempty(res)
        res2 = reshape(res, :, 7)
    else
        res2 = permutedims(hcat(res...))
    end
    correctdf=DataFrame(res2,colid)
    correctdf
end


function getchrbadsnp(calledchroffgeno::AbstractMatrix,
    bestgeno::AbstractMatrix, popmakeup::AbstractDict; minnerr::Integer=3)
    # mismatch[m, o] =1 if mismatch between posteiror_called geno and raw/input geno
    # nonmissing[m,o] if neither of two comparing genotypes is NN
    snpls = findall(sum(bestgeno,dims=2)[:,1] .> 0)
    mismatch,nonmissing = calgenomismatch(calledchroffgeno[snpls,:],bestgeno[snpls,:])
    ressnp=[]
    resww=[]
    nfounder  = length(unique(vcat([i["founder"] for i=values(popmakeup)]...)))
    for popid in keys(popmakeup)
        ppls = popmakeup[popid]["founder"]
        offls = popmakeup[popid]["offspring"]
        ppbool = zeros(nfounder)
        ppbool[ppls] .= 1
        popmismatch = view(mismatch, :, offls)
        popnonmiss = view(nonmissing, :, offls)
        # nerr[i]: number of mismatches at marker i
        nerr = sum(popmismatch,dims=2)[:,1]
        nsize = sum(popnonmiss,dims=2)[:,1]
        pos = findall(nerr .>= minnerr)
        ww = nerr[pos] ./ nsize[pos]
        minw = max(0.02,0.2/length(ppls))
        bool = ww .> minw
        snps = snpls[pos[bool]]
        ww2 = ww[bool]
        push!(ressnp,snps...)
        for i in eachindex(ww2)
            push!(resww,ww2[i]*ppbool)
        end
    end
    if isempty(ressnp)
        res  = ressnp
    else
        unisnp = sort(unique(ressnp))
        res = [(i,sum(resww[ressnp .== i])) for i=unisnp]
    end
    res
end


function calgenomismatch(calledchroffgeno::AbstractMatrix, bestgeno::AbstractMatrix)    
    # bestgeno: 1:6 .=> ["NN","1N","2N", "11", "12", "22"]    
    codes= [["NN","1N","2N", "11", "12", "22"],
        ["NN","1N","11","12"],
        ["NN","2N","22","12"],
        ["NN","1N","11"],
        ["NN","1N","2N","12"],
        ["NN","2N","22"]]
    postgenoset = codes[bestgeno]
    mismatch =  .!map(in,calledchroffgeno,postgenoset)
    postnonmissing = bestgeno .> 1
    callnonmissing = calledchroffgeno .!= "NN"
    nonmissing = postnonmissing .* callnonmissing
    mismatch, nonmissing
end

function get_chr_bestgeno(chrfhaplo::AbstractMatrix,chrviterbi::AbstractMatrix, 
    ancestralstate::AbstractVector)    
    size(chrfhaplo, 1) == size(chrviterbi,1)|| @error  "inconsistent #markers"    
    derdict=Dict(["NN","1N","2N", "11", "12", "22"] .=> 1:6)    
    bestgeno  = zeros(Int, size(chrviterbi)...)       
    snpls = findall(sum(chrviterbi,dims=2)[:,1] .> 0)
    chrfhaplo2 = permutedims(chrfhaplo)
    beststate2 = chrviterbi'        
    bestgeno2 = bestgeno'        
    for snp in snpls
        snpfhaplo = view(chrfhaplo2, :,snp)
        snpstate = view(beststate2, :,snp)
        snpgeno = view(bestgeno2,:,snp)
        ss = unique(snpstate)            
        dict = Dict([begin     
            g = join(sort(snpfhaplo[ancestralstate[s]]))
            g2 = get(derdict,g,missing)
            s => g2
        end for s in ss])
        snpgeno .= [get(dict,i,missing) for i in snpstate]
    end
    bestgeno   
end



function best_state2geno(newfhaplo::AbstractMatrix,newbeststate::AbstractMatrix,
	ancestralstate::AbstractVector)
    derdict=Dict(["NN","1N","2N", "11", "12", "22"] .=> 1:6)
    bestate2 = ancestralstate[newbeststate']
    newfhaplo2 = permutedims(newfhaplo)
    res = [begin 
        snphaplo = newfhaplo2[:,snp]
        gls = [join(sort(snphaplo[i])) for i in bestate2[:,snp]]
        [get(derdict,g,missing) for g in gls]   
    end for snp in 1:size(newfhaplo2,2)]
    reduce(hcat,res)'
end

