

function calpairwiseprior(magicped::MagicPed, model::AbstractString;
    isfounderinbred::Bool=true, isautosome::Bool=true)        
    # offspringinfo
    memdf = unique(magicped.offspringinfo[!,[:member,:ishomozygous,:isfglexch]])
    memberls = memdf[!,:member]    
    isfglexch_dict = Dict(memberls .=> memdf[!,:isfglexch])
    ishomozygous_dict = Dict(memberls .=> memdf[!,:ishomozygous])
    inputmodel=lowercase(model)
    model_dict= Dict([mem => (ishomozygous_dict[mem] ? "depmodel" : inputmodel) for mem in memberls])
    # designinfo 
    priorspace=MagicReconstruct.getpriorstatespace(magicped; isfounderinbred)    
    design=magicped.designinfo
    designtype=typeof(design)        
    isconcise = model!="jointmodel"
    recomprior = Dict{UInt,Any}()        
    if isnothing(designtype)
      error("design info is required, use setDesign! to provide design information")
    elseif designtype <: Dict{String,MagicBase.DesignInfo}      
        allfounders = magicped.founderinfo[!,:individual]    
        founderindexdict = Dict(allfounders.=> 1:length(allfounders))
        popmakeup = Dict([begin          
            subfounders = subdesign.founders            
            subped = subdesign.pedigree            
            if isnothing(subped)
                juncdist = subdesign.juncdist            
                isnothing(juncdist) && @error string("TODO for subdesign=",subdesign)                
                (; ibd,j1122,j1211,j1213,j1222,j1232) = juncdist                        
                phi12 = 1-ibd
                ismalex = false
                fglindicators = MagicReconstruct.get_fglindicators(allfounders,subfounders;isfounderinbred)
                dict=MagicPrior.identityprior(fglindicators,ismalex,phi12,j1122,j1211,j1213,
                    j1222,j1232,isconcise)
            else
                in(popid,subped.member) || @error string("popid =",popid, " is not in pedigree members=",ped.member)            
                dict = MagicPrior.magicsubprior(subped,allfounders; member = popid, isfounderinbred, 
                    isautosome,isfglexch=isfglexch_dict[popid], isconcise)                                      
                if model_dict[popid] == "jointmodel"                             
                    initprob = first(dict[popid].jointmodel)
                    nfgl = round(Int, sqrt(length(initprob)))
                    nfgl == length(priorspace["haploindex"]) || @error string("inconsistent nfgl")
                    ibd = sum(diag(reshape(initprob, nfgl,nfgl)))
                elseif model_dict[popid] == "depmodel"                
                    ibd = 1.0                    
                elseif model_dict[popid] == "indepmodel"      
                    ibd =  0.0                    
                else
                    @error string("unknown model=",model_dict[popid] )                 
                end
            end 
            prior = only(values(dict))        
            if model_dict[popid] == "jointmodel"
                model2 = ibd >= 0.85 ? "depmodel" : "indepmodel"
            else
                model2 =  model_dict[popid] 
            end
            nzstate,initprob,recomcoefs = contmarkov2prior(prior.contmarkov,model2)
            offspring=findall(magicped.offspringinfo[!,:member] .== popid)            
            ishaploid =  model2=="depmodel" || ishomozygous_dict[popid] 
            if ishaploid
                states = priorspace["haploindex"]
            else
                nfgl = length(priorspace["haploindex"])
                states = MagicBase.prior_diploindex(nfgl)
            end
            subfounder_indices = [founderindexdict[i] for i in subfounders]
            hashcode=hash(recomcoefs)
            haskey(recomprior, hashcode) || push!(recomprior, hashcode=>recomcoefs)
            popid => Dict(["founder"=>subfounder_indices,"offspring"=>offspring,
                "nzstate"=>nzstate, "ishaploid"=> ishaploid,
                "nzorigin"=>states[nzstate],
                "inbreedingcoef"=>ibd, "model"=>model2,
                "initprob"=>initprob,"hashcode"=>hashcode])
        end for (popid,subdesign) in design])
    elseif designtype <: Pedigree           
        isfglexchls = unique(values(isfglexch_dict))
        if all(isfglexchls)
            res = MagicPrior.magicprior(design; memberlist=memberls,
                isfounderinbred,isautosome,isfglexch=true,isconcise)
        elseif all(.!isfglexchls)
            res = MagicPrior.magicprior(design; memberlist=memberls,
                isfounderinbred,isautosome,isfglexch=false,isconcise)
        else
            res1 = MagicPrior.magicprior(design; memberlist=memberls[isfglexchls],
                isfounderinbred,isautosome,isfglexch=true,isconcise)
            res2 = MagicPrior.magicprior(design; memberlist=memberls[.!isfglexchls],
                isfounderinbred,isautosome,isfglexch=false,isconcise)
            res = merge(res1,res2)
        end
        popmakeup = Dict([begin      
            if model_dict[popid] == "jointmodel"                             
                initprob = first(prior.jointmodel)
                nfgl = round(Int, sqrt(length(initprob)))
                nfgl == length(priorspace["haploindex"]) || @error string("inconsistent nfgl")
                ibd = sum(diag(reshape(initprob, nfgl,nfgl)))
                model2 = ibd >= 0.85 ? "depmodel" : "indepmodel"
            elseif model_dict[popid] == "depmodel"                
                ibd = 1.0
                model2 = "depmodel"                                
            elseif model_dict[popid] == "indepmodel"      
                ibd =  0.0
                model2 = "indepmodel"               
            else
                @error string("unknown model=",model_dict[popid] )                 
            end
            nzstate, initprob, recomcoefs = contmarkov2prior(prior.contmarkov,model2)
            ishaploid =  model2=="depmodel" || ishomozygous_dict[popid] 
            if ishaploid
                states = priorspace["haploindex"]
            else
                nfgl = length(priorspace["haploindex"])
                states = MagicBase.prior_diploindex(nfgl)
            end             
            founders=sort(unique(vcat(states[nzstate]...)))
            founders = isfounderinbred ? founders : unique(div.(founders .+ 1,2))
            hashcode=hash(recomcoefs)
            haskey(recomprior, hashcode) || push!(recomprior, hashcode=>recomcoefs)
            offspring=findall(magicped.offspringinfo[!,:member] .== popid)            
            popid=>Dict(["founder"=>founders,"offspring"=>offspring,
                "nzstate"=>nzstate,"ishaploid"=> ishaploid,
                "nzorigin"=>states[nzstate],
                "inbreedingcoef"=>ibd, "model"=>model2,
                "initprob"=>initprob,"hashcode"=>hashcode])
        end for (popid,prior) in res])
    end
    popmakeup, recomprior
end


function contmarkov2prior(contmarkov::AbstractVector, model::AbstractString)
    if model == "depmodel"
        init  = (contmarkov[1][1]+contmarkov[2][1]) ./ 2
        init ./= sum(init)
        n = length(init)
        coef1 = [i==j ? init[i]*(-1+ init[i]) : init[i]*init[j] for i in 1:n, j in 1:n]
        coef0 = diagm(init)
        coefs = [coef0,coef1]
        nzstate = findall(init .!= 0)
        coefs = [i[nzstate, nzstate] for i in coefs]
        initprob = init[nzstate]
    else
        initm, initp = contmarkov[1][1], contmarkov[2][1]
        initm ./= sum(initm)
        initp ./= sum(initp)
        initmp = kron(initm, initp)
        coefm, coefp = [[i==j ? init[i]*(-1 + init[i]) : init[i]*init[j]
            for i in 1:length(init), j in 1:length(init)]
            for init in [initm, initp]]
        coef1 = kron(coefm, diagm(initp)) + kron(diagm(initm), coefp)
        # kron(diagm(1:3),diagm(2:4)) == diagm(kron(1:3,2:4))
        coef0 = diagm(initmp)
        coef2 = kron(coefm, coefp)
        coefs = [coef0,coef1,coef2]
        nzstate = findall(initmp .!= 0)
        coefs = [i[nzstate, nzstate] for i in coefs]
        initprob = initmp[nzstate]
    end
    # recomprior = r -> coefs[1] + sum(coefs[i] * r^(i-1) for i in 2:length(coefs))
    recomcoefs = coefs
    [nzstate, initprob, recomcoefs]
end

function calfhaploprior_loci(magicgeno::MagicGeno,chr::Integer)
    fgeno = magicgeno.foundergeno[chr]
    founderkind = only(unique(magicgeno.markermap[chr][!,:founderformat]))
    if founderkind=="GT_haplo"
        dict = Dict(["1"=>[1],"2"=>[2], "N"=>[1,2]])        
    elseif founderkind=="GT_unphased"
        dict = Dict(["11"=>[[1,1]],"12"=>[[1,2],[2,1]],"22"=>[[2,2]],"NN"=>[[1,1],[1,2],[2,1],[2,2]]])
        gset = unique(fgeno)
        d = setdiff(gset, keys(dict))
        if !isempty(d)
            dict2 = Dict(["1N"=>[[1,1],[1,2],[2,1]],"N1"=>[[1,1],[1,2],[2,1]],"21"=>[[1,2],[2,1]],
                "2N"=>[[2,2],[1,2],[2,1]],"N2"=>[[2,2],[1,2],[2,1]]])
            merge!(dict,dict2)
            d = setdiff(gset, keys(dict))
            isempty(d) || @error string("unexpected genotypes: ",d)
            filter!(x->in(x.first,gset),dict)
        end
    else
        # TODO founderkind=="unphasedgeno", "probability" (outbred)
        @error string("unknow kind of founder geno: ",founderkind)        
    end
    fhaploset0 = permutedims([get(dict, i, missing) for i in fgeno])
    fhaploset = [fhaploset0[:,i] for i in 1:size(fhaploset0,2)]
    fhaploset
end
