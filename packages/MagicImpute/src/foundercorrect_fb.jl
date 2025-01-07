
function foundercorrect_chr!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,    
    ismalexls::AbstractVector,
    founder2progeny::AbstractVector,
    popmakeup::AbstractDict,priorprocess::AbstractDict, priorspace::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector},     
    snporder::AbstractVector,
    hmmalg::AbstractString="forwardbackward",
    decodetempfile::AbstractString,    
    israndallele::Bool,    
    issnpGT::AbstractVector, 
    itmax::Integer=1)
    rescorrect = []
    for it in 1:itmax
        correctdf = last(foundererror_chr(chrfhaplo, chroffgeno, 
            ismalexls,founder2progeny, popmakeup,priorprocess,priorspace;
            epsf,epso, epso_perind, snporder, hmmalg,decodetempfile,israndallele,issnpGT))
        if !isempty(correctdf)
            setcorrectdf!(chrfhaplo,correctdf)
            push!(rescorrect,correctdf)
        end
        isempty(correctdf) && break
    end
    res = isempty(rescorrect) ? rescorrect : reduce(vcat, rescorrect)
    res
end

function foundererror_chr(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,    
    ismalexls::AbstractVector,
    founder2progeny::AbstractVector,
    popmakeup::AbstractDict,priorprocess::AbstractDict, priorspace::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector},     
    snporder::AbstractVector,
    hmmalg::AbstractString,
    decodetempfile::AbstractString,    
    israndallele::Bool,
    issnpGT::AbstractVector)
    loglike, _ = MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,
        popmakeup,priorprocess;
        epsf,epso, epso_perind, hmmalg, decodetempfile,
        israndallele, issnpGT,snporder)    
    MagicReconstruct.tomarker_major_order!(decodetempfile; hmmalg)        
	nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)	
	if nstate == nfgl
		ancestralstate = [[i,i] for i in priorspace["haploindex"]]
	else
		nstate == nfgl^2 || @error "inconsistent number of states"
        ancestralstate = MagicBase.prior_diploindex(nfgl)
	end    
    postcalled = callpostsnpprob(chrfhaplo,decodetempfile, snporder,ancestralstate; callthreshold=0.55)    
    # calculuate called offpsring genotypes
    calledchroffgeno = copy(chroffgeno)
    if !all(issnpGT)
        readcallthreshold = 0.95
        calledchroffgeno[.!issnpGT,:] .= MagicBase.callfromprob.(chroffgeno[.!issnpGT,:],
            readcallthreshold; isphased=false)
    end    
    gintersect = intersect(["21","N1","N2"],unique(calledchroffgeno))
    isempty(gintersect) || @error string("unexpected called offspring genotypes: ",gintersect)
    # TODO: to check the consistency of genotype calling for male x
    size(calledchroffgeno,2)==length(ismalexls) || @error "inconsistent size"
    if any(ismalexls)
        # male: N, 1, 2 => NN, 11, 22
        # male: NY, 1Y, 2Y => NN, 11, 22
        calledchroffgeno[:,ismalexls]=[i[1]*i[1] for i= calledchroffgeno[:,ismalexls]]
    end
    chrbadsnp = getchrbadsnp(calledchroffgeno,popmakeup, postcalled; minnerr = 1)
    # max 50% markers are labled as chrbadsnp
    wwls = [max(i[2]...) for i in chrbadsnp]
    q0=1.0-0.50*size(calledchroffgeno,1)/length(wwls)
    if q0 > 0
        q = quantile(wwls, q0)
        chrbadsnp = chrbadsnp[wwls.>q]
    end
    nselsnp = length(chrbadsnp)
    # getchrchange is time consuming
    # correctdf = getchrchange(chrfhaplo,calledchroffgeno,chrbadsnp, decodetempfile, 
    #     snporder, ancestralstate,founder2progeny; minnerrdiff = 1, callthreshold=0.55)
    # println("tused in getchrchange = ",tused)
    # loglike, nselsnp, correctdf
    chrfhaplo,calledchroffgeno,chrbadsnp, decodetempfile, snporder, ancestralstate
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

function setcorrectdf!(chrfhaplo::AbstractMatrix, correctdf::DataFrame)
    for i=1:size(correctdf,1)
        snp, p, newg = correctdf[i,[1,2,4]]
        chrfhaplo[snp,p]=newg
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
    chrbadsnp::AbstractVector, decodetempfile::AbstractString, snporder::AbstractVector,
    ancestralstate::AbstractVector,founder2progeny::AbstractVector;    
    minnerrdiff::Integer=3, callthreshold::Real=0.55)
    res=[]
    # println("#badsnp=",length(chrbadsnp), " out of #snp=",size(chrfhaplo,1))
    isfounderinbred = length(founder2progeny) == size(chrfhaplo,2)
	snpdict = Dict(snporder .=> 1:length(snporder))
	jldopen(decodetempfile, "r") do file
	    for (snp,ww) = chrbadsnp
			t = snpdict[snp]
			condprob = file[string(t)]
	        # pp = sortperm(ww,rev=true)
	        # pp =  pp[ww[pp] .> 0]
	        pp = findall((ww .>0) .&& (chrfhaplo[snp,:] .!= "N"))
	        change = [begin
                if isfounderinbred
	                genoset = ["1","2"]
                    oldfg=chrfhaplo[snp,p]                
	            else
	                genoset  = [[i,j] for i=["1","2"] for j=["1","2"]]
                    oldfg=chrfhaplo[snp,[2p-1,2p]]                
	            end
	            deleteat!(genoset, genoset .== oldfg)
	            pushfirst!(genoset,oldfg)
				snps = [snp for _ in 1:size(genoset,1)]
	            newfgeno=Matrix{Any}(chrfhaplo[snps,:])
	            newfgeno[:,p] .= genoset
	            progeny = founder2progeny[p]
				newcondprob = [condprob[progeny,:] for _ in 1:length(snps)]
	            newpostcalled = callpostsnpprob(newfgeno,newcondprob,ancestralstate; callthreshold)
	            newoffgeno = calledchroffgeno[snps,progeny]
	            # mismatch[m, o] =1 if mismatch between posteiror_called geno and raw/input geno
	            # nonmissing[m,o] if neither of two comparing genotypes is NN
				# println("newoffgeno = ")
				# display(newoffgeno)
	            mismatch,nonmissing = calgenomismatch(newoffgeno,newpostcalled)
	            nerr = sum(mismatch,dims=2)[:,1]
	            nnonmiss = sum(nonmissing,dims=2)[:,1]
	            newnerr,index = findmin(nerr)
	            newfg = genoset[index]
	            oldnerr = nerr[1]
	            [snp, p, oldfg,newfg,oldnerr,newnerr,nnonmiss[index]]
	        end for p in pp]
			change2=permutedims(hcat(change...))
	        val,ii=findmax(change2[:,5]-change2[:,6])
	        # reduced errors by an amound of at least minnerrdiff=3
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
    popmakeup::AbstractDict,postcalled::AbstractMatrix;
    minnerr::Integer=3)
    # mismatch[m, o] =1 if mismatch between posteiror_called geno and raw/input geno
    # nonmissing[m,o] if neither of two comparing genotypes is NN
    mismatch,nonmissing = calgenomismatch(calledchroffgeno,postcalled)
    ressnp=[]
    resww=[]
    nfounder  = length(unique(vcat([i["founder"] for i=values(popmakeup)]...)))
    for popid in keys(popmakeup)
        ppls = popmakeup[popid]["founder"]
        offls = popmakeup[popid]["offspring"]
        ppbool = zeros(nfounder)
        ppbool[ppls] .= 1
        popmismatch = mismatch[:, offls]
        popnonmiss = nonmissing[:,offls]
        # nerr[i]: number of mismatches at marker i
        nerr = sum(popmismatch,dims=2)[:,1]
        nsize = sum(popnonmiss,dims=2)[:,1]
        pos = findall(nerr .>= minnerr)
        ww = nerr[pos] ./ nsize[pos]
        minw = max(0.02,0.2/length(ppls))
        bool = ww .> minw
        snps = pos[bool]
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


function calgenomismatch(calledchroffgeno::AbstractMatrix, postcalled::AbstractMatrix)    
    # postcalled: 1:6 .=> ["NN","1N","2N", "11", "12", "22"]    
    codes= [["NN","1N","2N", "11", "12", "22"],
        ["NN","1N","11","12"],
        ["NN","2N","22","12"],
        ["NN","1N","11"],
        ["NN","1N","2N","12"],
        ["NN","2N","22"]]
    postgenoset = codes[postcalled]
    mismatch =  .!map(in,calledchroffgeno,postgenoset)
    postnonmissing = postcalled .> 1
    callnonmissing = calledchroffgeno .!= "NN"
    nonmissing = postnonmissing .* callnonmissing
    mismatch, nonmissing
end

function callpostsnpprob(prob::AbstractVector; callthreshold::Real)
    val,index=findmax(prob[4:6])
    if val>callthreshold
        res=index+3
    else
        val,index=findmax(prob[1:3])
        res = val>callthreshold ? index : 1
    end
    res
end

function callpostsnpprob(chrfhaplo::AbstractMatrix,decodetempfile::AbstractString,
	snporder::AbstractVector, ancestralstate::AbstractVector; callthreshold::Real)
    derdict=Dict(["NN","1N","2N", "11", "12", "22"] .=> 1:6)
	jldopen(decodetempfile,"r") do file
		nsnp = size(chrfhaplo,1)
		noff = size(file[string(1)],1)
		rescalled  = ones(Int, nsnp, noff)        
		for t in 1:nsnp            
			snp = snporder[t]
			condprob = file[string(t)]
			if issparse(condprob)
				offls, states,vals = findnz(condprob)
			else
				pos = findall(condprob .> 1e-8)
				pos2 = Tuple.(pos)
				offls = first.(pos2)
				states = last.(pos2)
				vals = condprob[pos]
				# nr,nc = size(condprob)
				# offls = repeat(1:nr,outer=nc)
				# states = repeat(1:nc,inner=nr)
				# vals = vec(condprob)
			end            
			# sort(["N","1"]) == ["1","N"]
			d0 = [join(sort(chrfhaplo[snp,s])) for s in ancestralstate[states]]
	        dd=[get(derdict,i,missing) for i in d0]
            ressnp = [zeros(6) for _ in 1:noff]
            for i in eachindex(offls)
                ressnp[offls[i]][dd[i]] += vals[i]
            end
            for p in ressnp
                p[2] += p[4]+p[5]
                p[3] += p[5]+p[6]
            end
            rescalled[snp,:] .=  callpostsnpprob.(ressnp; callthreshold)
		end		
		rescalled
	end
end

function callpostsnpprob(chrfhaplo::AbstractMatrix,chrcondprob::AbstractVector,
	ancestralstate::AbstractVector; callthreshold::Real)
    derdict=Dict(["NN","1N","2N", "11", "12", "22"] .=> 1:6)
    nsnp = size(chrfhaplo,1)
	noff = size(first(chrcondprob),1)	
	rescalled  = ones(Int, nsnp, noff)        
	for snp in 1:nsnp
		if issparse(chrcondprob[snp])
			offls, states,vals = findnz(chrcondprob[snp])
		else
			pos = findall(chrcondprob[snp] .> 1e-8)
			pos2 = Tuple.(pos)
			offls = first.(pos2)
			states = last.(pos2)
			vals = chrcondprob[snp][pos]
			# nr,nc = size(chrcondprob[snp])
			# offls = repeat(1:nr,outer=nc)
			# states = repeat(1:nc,inner=nr)
			# vals = vec(chrcondprob[snp])
		end
		d0 = [join(sort(chrfhaplo[snp,s])) for s in ancestralstate[states]]
		dd=[get(derdict,i,missing) for i in d0]        
        ressnp = [zeros(6) for _ in 1:noff]
        for i in eachindex(offls)
            ressnp[offls[i]][dd[i]] += vals[i]
        end
        for p in ressnp
            p[2] += p[4]+p[5]
            p[3] += p[5]+p[6]
        end
        rescalled[snp,:] .=  callpostsnpprob.(ressnp; callthreshold)
    end
    rescalled
end
