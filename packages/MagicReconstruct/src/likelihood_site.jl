

#######################called genotype and GBS prob vector######################
function calsitedataprob_multiphase(sitefhaplo::AbstractMatrix,sitegeno::AbstractVector,
    popmakeup::AbstractDict;
    epsf::Real,epso::Real, epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Real,allelebalancemean::Real,allelebalancedisperse::Real,alleledropout::Real,
    popidls=keys(popmakeup),
    offspringexcl::AbstractVector, 
    isoffphased::Bool=false,
    israndallele::Bool,
    issiteGT::Bool)
    nfgl = size(sitefhaplo,2)
    ishaploidls = [i["ishaploid"] for i in values(popmakeup)]
    nstate = false in ishaploidls ? nfgl^2 : nfgl
    noff = length(sitegeno)
    noff2 =sum([length(i["offspring"]) for i=values(popmakeup)])
    noff == noff2 || @error "inconsistent #offspring"    
    nphase = size(sitefhaplo,1)
    # dataprobls[i][j]: a vector for the i-th fhaplo and j_th offspring, length = nzstate for j_th offspring
    dataprobls = [Vector{Vector{_float_like}}(undef,noff) for i in 1:nphase]
    for popid = popidls
        offls = popmakeup[popid]["offspring"]
        nzstate =  popmakeup[popid]["nzstate"]
        ishaploid =  popmakeup[popid]["ishaploid"]           
        for i in 1:nphase
            for o in offls
                dataprobls[i][o] = zeros(_float_like, length(nzstate))
            end
        end 
        if issiteGT          
            sitelikelihood!(dataprobls, offls, sitefhaplo,sitegeno, epsf,epso,epso_perind,
                nzstate,nstate, ishaploid,isoffphased,israndallele)        
        else
            sitelikelihoodGBS!(dataprobls,offls,sitefhaplo,sitegeno,nzstate,nstate,
                ishaploid,epsf,epso,epso_perind,
                seqerror,allelebalancemean,allelebalancedisperse,alleledropout,israndallele)
        end
    end
    offls = offspringfrompop(popmakeup,popidls,offspringexcl)
    if issetequal(offls,1:noff)
        dataprobls
    else
        [i[offls] for i=dataprobls]
    end
end

function init_dataprobls_singlephase(popmakeup::AbstractDict)
    noff =sum([length(i["offspring"]) for i=values(popmakeup)])
    dataprobls = Vector{Vector{_float_like}}(undef,noff) 
    for popid = keys(popmakeup)        
        nzstate =  popmakeup[popid]["nzstate"]                
        offls =  popmakeup[popid]["offspring"]                
        for o in offls
            dataprobls[o] = zeros(_float_like, length(nzstate))
        end
    end
    dataprobls
end

function calsitedataprob_singlephase(sitefhaplo::AbstractVector,sitegeno::AbstractVector,
    popmakeup::AbstractDict;
    epsf::Real,epso::Real, epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Real,allelebalancemean::Real,allelebalancedisperse::Real,alleledropout::Real,    
    isoffphased::Bool=false,
    israndallele::Bool,
    issiteGT::Bool)    
    dataprobls = init_dataprobls_singlephase(popmakeup)    
    calsitedataprob_singlephase!(dataprobls,sitefhaplo,sitegeno,popmakeup;
        epsf, epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
        isoffphased,israndallele,issiteGT
    )
    dataprobls    
end



function calsitedataprob_singlephase!(dataprobls::AbstractVector, sitefhaplo::AbstractVector,sitegeno::AbstractVector,
    popmakeup::AbstractDict;
    epsf::Real,epso::Real, epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Real,allelebalancemean::Real,allelebalancedisperse::Real,alleledropout::Real,    
    isoffphased::Bool=false,
    israndallele::Bool,
    issiteGT::Bool)
    nstate, nfgl = hmm_nstate_nfgl(popmakeup)    
    nfgl == length(sitefhaplo) || @error string("inconsistent nfgl")
    if rand() < 0.01
        noff = length(sitegeno)
        noff2 =sum([length(i["offspring"]) for i=values(popmakeup)])
        noff == noff2 || @error "inconsistent #offspring"
    end    
    sitefhaplo2 = permutedims(sitefhaplo)
    dataprobls2 = [dataprobls]
    for popid = keys(popmakeup)
        offls = popmakeup[popid]["offspring"]
        nzstate =  popmakeup[popid]["nzstate"]
        ishaploid =  popmakeup[popid]["ishaploid"]
        if issiteGT          
            sitelikelihood!(dataprobls2, offls, sitefhaplo2,sitegeno, epsf,epso,epso_perind,
                nzstate,nstate, ishaploid,isoffphased,israndallele)        
        else
            sitelikelihoodGBS!(dataprobls2,offls,sitefhaplo2,sitegeno,nzstate,nstate,
                ishaploid,epsf,epso,epso_perind,
                seqerror,allelebalancemean,allelebalancedisperse,alleledropout,israndallele)
        end
    end    
    dataprobls    
end

############################Called genotype##################################

function sitelikelihood!(dataprobls::AbstractVector,offls::AbstractVector,
    sitefhaplo::AbstractMatrix,sitegeno::AbstractVector,
    epsf::Real,epso::Real, 
    epso_perind::Union{Nothing,AbstractVector}, 
    nzstate::AbstractVector, nstate::Integer,
    ishaploid::Bool, isoffphased::Bool,israndallele::Bool)
    nhaplo, nfgl=size(sitefhaplo)
    noff = length(offls)    
    if ishaploid                
        fderivehaplo = calfderive(sitefhaplo; ishaploid) # fderivehaplo does not depends on nzstate if ishaploid=true
        offcodehaplo = caloffcode(sitegeno[offls]; ishaploid,isoffphased)        
        if nstate == nfgl
            haplostate = 1:nfgl
            nzcol = nzstate
        elseif nstate == nfgl^2
            haplostate = [(i-1)nfgl+i for i=1:nfgl]
            nzcol = haplostate[nzstate]
        else
            error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))
        end
        if isnothing(epso_perind) || all(epso_perind .< 2e-5)
            # likehaplo = haplolike(epsf,epso)
            # likehaplo: size 3 obsvered haplotypes x 2 true haplotypes
            likehaplo=haplolike(epsf,epso)            
            for i=1:nhaplo
                prob = spzeros(noff,nstate)
                prob[:, haplostate] .= likehaplo[offcodehaplo,fderivehaplo[i,:]]                      
                for o in eachindex(offls)
                    dataprobls[i][offls[o]] .= prob[o,nzcol]
                end                
            end
        else        
            offprob = spzeros(nstate)    
            likehaplo = haplolike(epsf,epso)            
            for i=1:nhaplo                
                for o in eachindex(offls)                    
                    epso2 = epso_perind[offls[o]]
                    epso2 = epso + epso2 - epso * epso2
                    likehaplo .= haplolike(epsf,epso2)            
                    offprob[haplostate] .= likehaplo[offcodehaplo[o],fderivehaplo[i,:]]            
                    dataprobls[i][offls[o]] .= offprob[nzcol]
                end
            end
        end
    else
        nstate == nfgl^2 || error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))
        fderivediplo = calfderive(sitefhaplo; nzstate, ishaploid)
        offcodediplo= caloffcode(sitegeno[offls]; ishaploid,isoffphased)
        ibdbool= falses(nfgl^2)
        ibdbool[[(i-1)nfgl+i for i=1:nfgl]] .= true
        nonibdbool = .!ibdbool
        # missing pattern is the same for each row        
        # missing beeing zero in the sparse matrix of fderive        
        zzstate = setdiff(1:nstate, nzstate)
        ibdbool[zzstate] .= false
        nonibdbool[zzstate] .= false
        if isnothing(epso_perind) || all(epso_perind .< 2e-5)
            likediplo=diplolike(epsf,epso; isoffphased,israndallele)            
            for i=1:nhaplo
                prob = spzeros(noff,nstate)
                prob[:, ibdbool] .= likediplo.ibd[offcodediplo,fderivediplo[i,ibdbool]]
                prob[:, nonibdbool] .= likediplo.nonibd[offcodediplo,fderivediplo[i,nonibdbool]]                
                for o in eachindex(offls)
                    dataprobls[i][offls[o]] .= prob[o,nzstate]
                end
            end
        else
            offprob = spzeros(nstate)            
            for i=1:nhaplo
                for o in eachindex(offls)                    
                    epso2 = epso_perind[offls[o]]
                    epso2 = epso + epso2 - epso * epso2
                    likediplo = diplolike(epsf,epso2; isoffphased,israndallele)            
                    offprob[ibdbool] .= likediplo.ibd[offcodediplo[o],fderivediplo[i,ibdbool]]
                    offprob[nonibdbool] .= likediplo.nonibd[offcodediplo[o],fderivediplo[i,nonibdbool]]
                    dataprobls[i][offls[o]] .= offprob[nzstate]
                end
            end
        end
    end
    dataprobls
end


##########################GBS=>probability vector##############################

function sitelikelihoodGBS!(dataprobls::AbstractVector,offls::AbstractVector,
    sitefhaplo::AbstractMatrix,sitegeno::AbstractVector,
    nzstate::AbstractVector,nstate::Integer,
    ishaploid::Bool, epsf::Real,epso::Real, epso_perind::Union{Nothing,AbstractVector}, 
    seqerror::Real,allelebalancemean::Real,allelebalancedisperse::Real,alleledropout::Real,
    israndallele::Bool)
    nhaplo, nfgl=size(sitefhaplo)
    noff = length(offls)
    if ishaploid        
        # fderivehaplo size: (nphase,nfgl) for ishaploid = true
        fderivehaplo = calfderive(sitefhaplo; ishaploid) # fderivehaplo does not depends on nzstate if ishaploid=true
        if nstate == nfgl
            haplostate = 1:nfgl
            nzcol = nzstate
        elseif nstate == nfgl^2
            haplostate = [(i-1)nfgl+i for i=1:nfgl]
            nzcol = haplostate[nzstate]
        else
            error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))
        end
        offgenols = sitegeno[offls]
        dataformat = typeof(first(offgenols)) 
        if dataformat <: AbstractVector            
            # jldopen("test.jld2","w") do file
            #     file["breakpoint"] = [offls,sitefhaplo,sitegeno,nzstate,nstate,ishaploid,epsf,epso,seqerror]
            # end
            # error("TODO: GBS for haploid offspring")
            datatype = eltype(first(offgenols)) 
            if isnothing(epso_perind) || all(epso_perind .< 2e-5)
                epsols = epso
            else
                epsols = [epso + i -  epso * i for i in epso_perind[offls]]
            end
            if datatype <: Integer
                like = MagicReconstruct.haplolikeGBS(offgenols,epsols,seqerror)  # format = AD
            elseif datatype <: AbstractFloat                
                like = MagicReconstruct.haplolikeGBS(offgenols,epsols) # format = GP
            else
                @error string("TODO for datatype=", datatype, " for genotypes=",offgenols)            
            end
            prior = MagicReconstruct.haploprior(epsf)
            likehaplo=like*prior
            for i=1:nhaplo
                prob = spzeros(noff,nstate)
                prob[:,haplostate] .= likehaplo[:,fderivehaplo[i,:]]
                for o in eachindex(offls)
                    dataprobls[i][offls[o]] .= prob[o,nzcol]
                end
            end
        elseif dataformat <: AbstractString
            # called genotypes for the subpopulation 
            genoset = unique(offgenols)
            isinbred = issubset(genoset,["NN","11","22"])
            ismalex = issubset(genoset,["N","1","2"])
            if !(isinbred || ismalex)
                @error(string("unexpected sitegeno: ", genoset))
            end            
            offcodehaplo =caloffcode(offgenols; ishaploid, isoffphased=false)            
            if isnothing(epso_perind) || all(epso_perind .< 2e-5)                
                # likehaplo: size 3 obsvered haplotypes x 2 true haplotypes
                likehaplo=haplolike(epsf,epso)            
                for i=1:nhaplo
                    prob = spzeros(noff,nstate)
                    prob[:, haplostate] .= likehaplo[offcodehaplo,fderivehaplo[i,:]]                                
                    for o in eachindex(offls)
                        dataprobls[i][offls[o]] .= prob[o,nzcol]
                    end
                end
            else        
                offprob = spzeros(nstate)    
                likehaplo = haplolike(epsf,epso)            
                for i=1:nhaplo                
                    for o in eachindex(offls)                        
                        epso2 = epso_perind[offls[o]]
                        epso2 = epso + epso2 - epso * epso2
                        likehaplo .= haplolike(epsf,epso2)            
                        offprob[haplostate] .= likehaplo[offcodehaplo[o],fderivehaplo[i,:]]            
                        dataprobls[i][offls[o]] .= offprob[nzcol]
                    end
                end
            end
        else
            @error string("TODO for dataformat=", dataformat, " for genotypes=",offgenols)            
        end
    else
        nstate == nfgl^2 || error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))
        zzstate = setdiff(1:nstate, nzstate)
        fderivediplo = calfderive(sitefhaplo; nzstate,ishaploid)
        sitegeno_off = sitegeno[offls]
        if isnothing(epso_perind)  || all(epso_perind .< 2e-5)
            epsols = epso
        else
            epsols = [epso + i -  epso * i for i in epso_perind[offls]]
        end
        if eltype(first(sitegeno_off)) <: Integer
            like = diplolikeGBS(sitegeno_off,epsols,seqerror,allelebalancemean,allelebalancedisperse,alleledropout; israndallele)
        else
            like = diplolikeGBS(sitegeno_off,epsols; israndallele)
        end
        prior = diploprior(epsf)
        likenonibd=like*prior.nonibd
        likeibd=like*prior.ibd
        ibdbool= falses(nfgl^2)
        ibdbool[[(i-1)nfgl+i for i=1:nfgl]] .= true
        nonibdbool = .!ibdbool
        ibdbool[zzstate] .= false
        nonibdbool[zzstate] .= false
        for i=1:nhaplo
            prob = spzeros(noff,nstate)
            prob[:,ibdbool] .= likeibd[:, fderivediplo[i,ibdbool]]
            prob[:,nonibdbool] .= likenonibd[:, fderivediplo[i,nonibdbool]]            
            for o in eachindex(offls)
                dataprobls[i][offls[o]] .= prob[o,nzstate]
            end
        end
    end
    dataprobls
end

function offspringfrompop(popmakeup::AbstractDict,popidls,offspringexcl::AbstractVector)
    offls = sort(reduce(vcat,[popmakeup[i]["offspring"] for i=popidls]))
    setdiff!(offls,offspringexcl)
    offls
end
