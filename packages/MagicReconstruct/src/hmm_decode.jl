
function hmmdecode_chr(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,Real,AbstractVector},
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    hmmalg::AbstractString,
    posteriordigits::Integer = 4, 
    decodetempfile::Union{Nothing,AbstractString},    
    resetvirtual::Bool=false,
    issnpGT::AbstractVector,
    isoffphased::Bool=false,
    israndallele::Bool, 
    popidls::AbstractVector=collect(keys(popmakeup)),
    snporder::Union{Nothing,AbstractVector}=nothing)
    if isnothing(decodetempfile)
        decodetempfile = tempname(tempdir(); cleanup=true)
    end    
    fderive, offcode  = precompute_chr(chrfhaplo, chroffgeno, popmakeup,isoffphased, issnpGT)
    nstate, nfgl = hmm_nstate_nfgl(popmakeup)    
    nfgl == size(chrfhaplo,2) || @error string("inconsistent nfgl")
    nsnp,noff=size(chroffgeno)
    isnothing(snporder) && (snporder = 1:nsnp)    
    if isnothing(epso_perind) 
        epso_perind = zeros(noff)
    elseif isa(epso_perind,Real)
        epso_perind = epso_perind .* ones(noff)
    else
        length(epso_perind) == noff || @error string("length of epso_perind is not equal to #offspring=",noff)
    end
    loglike=zeros(noff)           
    jldopen(decodetempfile,"w") do file
        write(file, "hmmalg", hmmalg)
        write(file, "snporder", snporder)
        write(file, "noffspring", noff)
        write(file, "tlength", nsnp)
        write(file, "nstate", nstate)        
        old_hmm_nstate = length(popmakeup[popidls[1]]["nzstate"])
        dataprobseq = [zeros(_float_like,old_hmm_nstate)  for _ in 1:nsnp]
        for popid in popidls
            # startt = time()            
            ishaploidsub = popmakeup[popid]["ishaploid"]
            offls = popmakeup[popid]["offspring"]            
            virtualdict = popmakeup[popid]["virtualdict"]            
            # nzstate: the set of indices with states being acessible in the subpopulation
            nzstate = popmakeup[popid]["nzstate"]         
            nzcol = get_nzcol(nstate,popmakeup,popid)           
            hashcode = popmakeup[popid]["hashcode"]                  
            markerincl, initprob, tranprobseq = hashcode2prior(priorprocess,hashcode)
            # length(initprob) == length(nzstate)         
            if length(nzstate) != old_hmm_nstate
                old_hmm_nstate = length(nzstate)
                dataprobseq = [zeros(_float_like,old_hmm_nstate)  for _ in 1:nsnp]
            end            
            for off in offls
                if hmmalg in ["forwardbackward","logforward","logbackward"]                    
                    off_decode = length(nzcol)/nstate > 0.5 ? zeros(nsnp,nstate) : spzeros(nsnp,nstate)
                else
                    # hmmalg == "viterbi"
                    off_decode = zeros(Int,nsnp)
                end
                if resetvirtual && haskey(virtualdict, off) 
                    # copy founders'geno to virtual offspring for intermeidate parents 
                    virtual_f = virtualdict[off]
                    if isfounderinbred 
                        newgeno = chrfhaplo[:,virtual_f] .^ 2
                    else
                        newgeno = join.(eachrow(chrfhaplo[:,[2*virtual_f-1,2*virtual_f]]))
                    end
                    isoffphased2 = false
                    issnpGT2 = trues(length(issnpGT))
                    obsseq = caloffcode(newgeno; ishaploid=ishaploidsub, isoffphased=isoffphased2)                    
                    caldataprobseq!(dataprobseq, obsseq,epsf,0.005,0.0, 0.001, 0.5, 0.0, 0.0, 
                        fderive,nzstate,isoffphased2,israndallele, issnpGT2,ishaploidsub)
                else
                    obsseq = offcode[:,off]
                    caldataprobseq!(dataprobseq, obsseq,epsf,epso, epso_perind[off], seqerror, allelebalancemean,allelebalancedisperse,alleledropout,
                        fderive,nzstate,isoffphased,israndallele, issnpGT,ishaploidsub)
                end
                dataprobseq2 = view(dataprobseq, snporder[markerincl])    
                if hmmalg == "forwardbackward"                    
                    loglike[off],prob=HMM.posteriordecode(initprob,tranprobseq,dataprobseq2; digits=posteriordigits) # larger digits  results in more memoryuse!!
                    off_decode[markerincl,nzcol] .= reduce(hcat,prob)' 
                    issparse(off_decode) && dropzeros!(off_decode)
                elseif hmmalg == "calloglike"                    
                    loglike[off]= HMM.calloglike(initprob,tranprobseq,dataprobseq2)
                elseif hmmalg == "logforward"
                    prob = HMM.logforward(initprob,tranprobseq,dataprobseq2)
                    off_decode[markerincl,nzcol] .= reduce(hcat,prob)'
                elseif hmmalg == "logbackward"
                    isempty(tranprobseq) && @error "empty tranprobseq"
                    prob = HMM.logbackward(tranprobseq,dataprobseq2)
                    off_decode[markerincl,nzcol] .= reduce(hcat,prob)'
                elseif hmmalg == "viterbi"
                    loglike[off],vpath = HMM.viterbidecode(initprob,tranprobseq,dataprobseq2)                    
                    off_decode[markerincl] .= nzcol[vpath]        
                else
                    @error string("unknown hmm algorithm: ",hmmalg)
                end
                if hmmalg != "forward"                    
                    if !issparse(off_decode) && hmmalg == "forwardbackward"                     
                        write(file, string(off), (nzcol,sparse(off_decode)))
                    else
                        write(file, string(off), (nzcol,off_decode))
                    end
                end
            end       
            # @info string("popid=", popid, ", tsused=", round(time()-startt,digits=1), "s, mem=", round(Int, memoryuse()/10^6),"MB")     
        end
    end
    if hmmalg in ["logforward","logbackward"]
        decodetempfile
    else
        # "forwardbackward" or "viterbi"
        (loglike, decodetempfile)
    end
end


function read_decodefile(decodefile::AbstractString)
    jldopen(decodefile, "r") do file
        hmmalg = file["hmmalg"]
        if !in(hmmalg,["viterbi", "forwardbackward","logforward","logbackward"]) 
            error(string("tomarker_major_order! does not work for hmmalg=",hmmalg))
        end
        snporder = file["snporder"]
        noffspring = file["noffspring"]
        decode = [haskey(file,string(off)) ? file[string(off)] : nothing for off in 1:noffspring]
        hmmalg, snporder, decode
    end    
end


function caldataprobseq!(dataprobseq::AbstractVector, 
    obsseq::AbstractVector,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_ind::Real,     
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    fderive::NamedTuple,
    nzstate::AbstractVector,
    isoffphased::Bool,    
    israndallele::Bool,
    issnpGT::AbstractVector,
    ishaploid::Bool)    
    @assert length(obsseq) == length(issnpGT)
    if epso_ind > 0
        if isa(epso,AbstractVector)
            epso = [i + epso_ind - i * epso_ind for i in epso]
        else
            epso = epso + epso_ind - epso * epso_ind
        end
    end    
    if !all(issnpGT)
        issnpnonGT = .!issnpGT
        fderive2 = ishaploid ? view(fderive.haplo, issnpnonGT,:) : view(fderive.diplo,issnpnonGT,:)
        if isa(epsf,AbstractVector) || isa(epso,AbstractVector) || isa(seqerror,AbstractVector) || isa(allelebalancemean,AbstractVector) || isa(allelebalancedisperse,AbstractVector) || isa(alleledropout,AbstractVector) 
            f(x) = isa(x,AbstractVector) ? x[issnpnonGT] : x
            epsf2, epso2,seqerror2,allelebalancemean2,allelebalancedisperse2,alleledropout2 = f(epsf),f(epso),f(seqerror),f(allelebalancemean),f(allelebalancedisperse), f(alleledropout)
        else
            epsf2, epso2,seqerror2,allelebalancemean2,allelebalancedisperse2,alleledropout2 = epsf, epso,seqerror,allelebalancemean,allelebalancedisperse,alleledropout
        end        
        callinelikeGBS!(view(dataprobseq,issnpnonGT),fderive2,nzstate, view(obsseq,issnpnonGT);
            epsf=epsf2,epso=epso2,seqerror = seqerror2, 
            allelebalancemean = allelebalancemean2, 
            allelebalancedisperse = allelebalancedisperse2, 
            alleledropout = alleledropout2, ishaploid,israndallele)
    end
    if any(issnpGT)
        # condlike is a function of epsf, epso        
        condlike = precompute_condike(epsf,epso; isoffphased,israndallele)        
		if isa(epsf,AbstractVector) || isa(epso,AbstractVector)
        	condlike2 = ishaploid ? condlike.haplo[issnpGT] : condlike.diplo[issnpGT]
		else
			condlike2 = ishaploid ? condlike.haplo : condlike.diplo
		end        
        fderive2 = ishaploid ? view(fderive.haplo, issnpGT,:) : view(fderive.diplo,issnpGT,:)
        callinelikesnp!(view(dataprobseq,issnpGT), condlike2,fderive2,nzstate, view(obsseq,issnpGT); ishaploid);
    end
    dataprobseq    
end

function hmm_loglikels(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,        
    popmakeup::AbstractDict,priorprocess::AbstractDict;    
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    snporder::AbstractVector,    
    decodetempfile::AbstractString,    
    hmmalg::AbstractString="calloglike",
    israndallele::Bool,    
    issnpGT::AbstractVector)
    isnothing(snporder) && (snporder = collect(1:size(chroffgeno,1)))    
    loglikels = first(MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf,epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
        hmmalg, decodetempfile, israndallele,issnpGT,snporder))
    loglikels
end

function get_chr_viterbi(viterbitempfile::AbstractString, snporder::AbstractVector)
    jldopen(viterbitempfile,"r") do file
        noff = file["noffspring"]
        nsnp = file["tlength"] 
        length(snporder) == nsnp || @error  "inconsistent #markers"        
        beststate  = zeros(Int, nsnp,noff)        
        for off in 1:noff
            _,off_decode = file[string(off)]
            beststate[snporder,off] .= off_decode
        end
        beststate
    end
end
