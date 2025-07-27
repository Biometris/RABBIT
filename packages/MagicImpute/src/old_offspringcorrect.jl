
function offspringcorrect_chr!(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,    
    popmakeup::AbstractDict,priorprocess::AbstractDict, priorspace::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector},
    baseerror::Union{Real,AbstractVector},
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},
    allelicdropout::Union{Real,AbstractVector},    
    snporder::AbstractVector,    
    correctalg::AbstractString,
    decodetempfile::AbstractString,    
    ismalexls::AbstractVector,    
    israndallele::Bool,     
    correctthesh::Real,
    threshlikeparam::ThreshLikeParam)    
    mismatch = get_mismatch(chrfhaplo,chroffgeno,popmakeup,priorprocess,priorspace;
        epsf,epso, epso_perind, baseerror,allelicbias,allelicoverdispersion,allelicdropout,        
        snporder, correctalg, decodetempfile,ismalexls,israndallele, issnpGT,correctthesh
    )
    nerrorls = sum(mismatch,dims=2)[:,1]
    noff = size(chroffgeno,2)
    delsnps_off = findall(nerrorls .> 2*threshlikeparam.offspringerror*noff)        
    isempty(delsnps_off) || MagicReconstruct.setpriorprocess!(priorprocess,snporder, delsnps_off)    
    # set offspring errorneous genotypes to missing
    mismatch[delsnps_off,:] .= false    
    chroffgeno[mismatch] .= [begin 
        et =eltype(i)
        if et <: AbstractChar 
            "NN" 
        elseif et <: Integer
            et[0,0]
        elseif et <: AbstractFloat
            length(i)==3 ? [0.25,0.5,0.25] : ones(length(i))/length(i)
        else
            @error string("unexpected offspring genotype = ",i)
            nothing
        end
    end for i in chroffgeno[mismatch]]    
    length(delsnps_off), sum(mismatch)
end


function get_mismatch(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix, 
    popmakeup::AbstractDict,priorprocess::AbstractDict, priorspace::AbstractDict;
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector},
    baseerror::Union{Real,AbstractVector},
    allelicbias::Union{Real,AbstractVector},
    allelicoverdispersion::Union{Real,AbstractVector},
    allelicdropout::Union{Real,AbstractVector},
    snporder::AbstractVector,        
    issnpGT::AbstractVector, 
    correctalg::AbstractString,
    decodetempfile::AbstractString,    
    ismalexls::AbstractVector,    
    israndallele::Bool,     
    correctthesh::Real)  
    mismatch = falses(size(chroffgeno)...)    
    tempjld2file = tempname(dirname(abspath(decodetempfile)); cleanup = true)		
    MagicImpute.singlesite_genoprob(chrfhaplo,chroffgeno; popmakeup,
        epsf,epso, epso_perind, baseerror, allelicbias, allelicoverdispersion,
        allelicdropout, issnpGT,outjld2file = tempjld2file)
    pri1=first(values(priorprocess))
	snpincl = snporder[pri1.markerincl]
    if correctalg in ["viterbi","viterbi_forwardbackward"]
        MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
            epsf,epso, epso_perind, baseerror,allelicbias,allelicoverdispersion,allelicdropout,        
            hmmalg="viterbi", decodetempfile, israndallele,issnpGT,snporder)       
        chrviterbi = MagicReconstruct.get_chr_viterbi(decodetempfile, snporder)
        nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)	
        if nstate == nfgl
            ancestralstate = [[i,i] for i in priorspace["haploindex"]]
        else
            nstate == nfgl^2 || @error "inconsistent number of states"
            ancestralstate = MagicBase.prior_diploindex(nfgl)
        end 
        bestgeno = get_chr_bestgeno(chrfhaplo, chrviterbi,ancestralstate)        
        snpls = findall(sum(bestgeno,dims=2)[:,1] .> 0)        
        intersect!(snpls,snpincl)
        calledchroffgeno = copy(chroffgeno)
        singlesite_callfromprob!(calledchroffgeno,tempjld2file;ismalexls,issnpGT,callthreshold=correctthesh)
        mismatch[snpls,:] .= first(calgenomismatch(view(calledchroffgeno, snpls,:), view(bestgeno,snpls,:)))    
        # println("mismatch df=", DataFrame(offgeno=chroffgeno[mismatch],calledchroffgeno=calledchroffgeno[mismatch],pos=findall(mismatch)))            
    end
    if correctalg in ["forwardbackward","viterbi_forwardbackward"]
        MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
            epsf,epso, epso_perind, baseerror,allelicbias,allelicoverdispersion,allelicdropout,        
            hmmalg="forwardbackward", decodetempfile, israndallele, issnpGT,snporder)                   
        MagicImpute.condprob2postprob!(chrfhaplo,decodetempfile,tempjld2file,popmakeup; epsf)         
        MagicImpute.singlesite_genoprob(chrfhaplo,chroffgeno; popmakeup,
            epsf,epso, epso_perind, baseerror, allelicbias, allelicoverdispersion,
            allelicdropout, issnpGT,outjld2file = decodetempfile)              
        mismatch2 = MagicImpute.get_mismatch(decodetempfile, tempjld2file; snpincl, threshold=correctthesh)   
        # println("nsnp=",size(chrfhaplo,1),",nsnpincl=",length(snpincl)) 
        # println("mismatch df=", DataFrame(offgeno=chroffgeno[mismatch2],os=findall(mismatch2)))            
        mismatch .= mismatch .|| mismatch2               
    end   
    rm(tempjld2file;force=true)
    mismatch
end
function get_mismatch(postprobfile::AbstractString, rawprobfile::AbstractString; 
    snpincl::AbstractVector, threshold::Real)
    jldopen(postprobfile,"r") do probfile
        jldopen(rawprobfile,"r") do file
            nmarker = probfile["nmarker"]
            noffspring = probfile["noffspring"]
            file["nmarker"] == nmarker || @error "inconsistent #marker"
            file["noffspring"] == noffspring || @error "inconsistent #offspring"            
            mismatch = falses(nmarker, noffspring)                        
            for off in 1:noffspring
                postprobls = probfile[string(off)]
                rawprobls = file[string(off)]
                for snp in snpincl
                    if isa(rawprobls[snp],AbstractVector)       
                        p, i = findmax(postprobls[snp])                 
                        p2, i2 = findmax(rawprobls[snp])                 
                        mismatch[snp,off] = i != i2 && min(p,p2)> threshold
                        # if mismatch[snp,off] 
                        #     println("snp=",snp," ;postprob=",postprobls[snp],",rawprob=",rawprobls[snp])
                        # end
                    end
                end
            end
            mismatch
        end
    end
end
