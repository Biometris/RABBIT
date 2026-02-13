# function getfhaploweight(findex::AbstractVector, fhaploweightpp::AbstractVector)
#     if length(findex)==1
#         deepcopy(fhaploweightpp[findex[1]])
#     else
#         weight = hcat(fhaploweightpp[findex]...)
#         [kron(i...) for i=eachrow(weight)]
#     end
# end

function getfhaploset(findex::AbstractVector,fhaplophase::AbstractMatrix,
    fhaplosetpp::AbstractVector)
    nsnp = size(fhaplophase,1)
    hh = permutedims(fhaplophase)
    ff= permutedims(reduce(hcat,fhaplosetpp[findex]))
    isfounderinbred = length(fhaplosetpp) ==size(fhaplophase,2)
    if isfounderinbred                
        bool = length(findex)==1        
        fhaploset = [begin            
            a = bool ? [[i] for i= ff[1,m]] : split.(kron(ff[:,m]...),"")
            b = reduce(hcat,repeat([hh[:,m]],length(a)))
            b[findex,:] = reduce(hcat,a)
            permutedims(b)            
        end for m=1:nsnp]
    else
        # outbred parents
        findex2 = reduce(vcat,[[2*i-1,2*i] for i in findex])
        fhaploset = [begin 
            a = [[i...] for i in product(ff[:,m]...)]
            b = reduce(hcat,repeat([hh[:,m]],length(a)))
            a2 = reduce(hcat,[reduce(vcat,i) for i in vec(a)])
            b[findex2,:] .= a2
            permutedims(b)
        end for m = 1:nsnp]
    end
    fhaploset

end

function get_initfhaplo(fhaplosetpp::AbstractVector; isimputefounder, isfounderinbred,minfmiss)    
    fmissls = get_fmissls(fhaplosetpp; isfounderinbred) 	   
    if isnothing(isimputefounder)
        isimputels = fmissls .>= minfmiss    
        isimputefounder = any(isimputels)        
    else
        if isimputefounder
            isimputels = fmissls .>= minfmiss    
        else
            isimputels = falses(length(fmissls))
        end
    end
    gmiss = isa(fhaplosetpp[1][1][1],AbstractVector) ? ["N","N"] : "N"
    fhaplols = reduce(hcat, [begin 
        ff_haplo = fhaplosetpp[f]
        if isimputels[f]
            [rand(i) for i in ff_haplo]
        else
            [length(i)==1 ? only(i) : gmiss for i in ff_haplo]
        end    
    end for f in eachindex(fhaplosetpp)])
    fhaplo = permutedims(reduce(hcat,[reduce(vcat,i) for i in eachrow(fhaplols)]))    
    fhaplo,fmissls, isimputefounder
end

function randfhaplophase(fhaplosetpp::AbstractVector)    
    ls = [[rand(i) for i in h] for h in fhaplosetpp]
    fhaplo0 = reduce(hcat,ls)
    fhaplo = permutedims(reduce(hcat,[reduce(vcat,i) for i in eachrow(fhaplo0)]))
    fhaplo
end

function calfhaploprior(magicgeno::MagicGeno,chr::Integer)
    fgeno = magicgeno.foundergeno[chr]
    fhaploset = Matrix{Any}(undef, size(fgeno)...)
    formatls = magicgeno.markermap[chr][!,:founderformat]
    formatset = unique(formatls)
    d = setdiff(formatset,["GT_haplo","GT_phased","GT_unphased","GT"])
    isempty(d) || @error string("calfhaploprior does work for founder format=",d)
    for founderkind in formatset        
        b = formatls .== founderkind
        gset = unique(view(fgeno,b,:))
        if founderkind=="GT_haplo"
            gdiff = setdiff(gset,["1","2","N"])
            isempty(gdiff) || @error string("unexpected genotypes: ",gdiff)
            grule = Dict(["1"=>["1"],"2"=>["2"], "N"=>["1","2"]])            
            fhaploset[b,:] .= [grule[i] for i in view(fgeno,b,:)]
        else           
            grule_phased = Dict([["1","N"] =>[["1","1"],["1","2"]],["2","N"] =>[["2","1"],["2","2"]],
                ["N","1"] =>[["1","1"],["2","1"]],["N", "2"] =>[["1","2"],["2","2"]],
                ["N","N"] =>[["1","1"],["1","2"],["2","1"],["2","2"]],
                ["1","1"]  => [["1","1"]], ["2","2"]  => [["2","2"]], 
                ["1","2"]  => [["1","2"]],["2","1"]  => [["2","1"]]]
            )
            # grule_unphased
            grule_unphased = Dict(["11"=>["11"],"12"=>["12","21"],"22"=>["22"], "NN"=>["11","12","21","22"]])
            gmissrule_unphased = Dict(["1N"=>["11","12","21"],"N1"=>["11","12","21"],
                "2N"=>["22","12","21"],"N2"=>["22","12","21"]])
            merge!(grule_unphased, gmissrule_unphased)        
            grule_unphased = Dict([i =>split.(grule_unphased[i],"") for i in keys(grule_unphased)])

            if founderkind=="GT_unphased"                
                gdiff = setdiff(gset,keys(grule_unphased))      
                isempty(gdiff) || @error string("unexpected genotypes: ",gdiff, ", founderformat=",founderkind)                                  
                grule2 = Dict([i =>grule_unphased[i] for i in gset])
                fhaploset[b,:] .= [grule2[i] for i in view(fgeno,b,:)]
            elseif founderkind=="GT_phased"                
                gdiff = setdiff(gset,keys(grule_phased))        
                isempty(gdiff) || @error string("unexpected genotypes: ",gdiff, ", founderformat=",founderkind)                                  
                grule2 = Dict([i =>grule_phased[i] for i in gset])
                fhaploset[b,:] .= [grule2[i] for i in view(fgeno,b,:)]
            elseif founderkind=="GT"
                grule = merge(grule_unphased, grule_phased)        
                gdiff = setdiff(gset,keys(grule))        
                isempty(gdiff) || @error string("unexpected genotypes: ",gdiff, ", founderformat=",founderkind)                                  
                grule2 = Dict([i =>grule[i] for i in gset])
                fhaploset[b,:] .= [grule2[i] for i in view(fgeno,b,:)]
            else        
                @error string("unknown kind of foundergeno: ",founderkind)
            end
        end
    end    
    collect.(eachcol(fhaploset))
end

# old codes
function calfweight(hhset::AbstractVector,obshh::AbstractVector, epsf::Real,
    isfounderinbred::Bool,maxnerror::Integer)
    dict = Dict(["1"=>1,"2"=>2,"N"=>missing])
    hhset2 = [[get(dict,i,missing) for i= h] for h=hhset]
    obshh2=[get(dict,i,missing) for i= obshh]
    if !isfounderinbred
        for i = 1:length(hhset2)
            h=sort(reshape(hhset2[i],2,:),dims=1)
            hhset2[i] =reshape(h,length(h))
        end
        h=sort(reshape(obshh2,2,:),dims=1)
        obshh2=reshape(h,length(h))
    end
    bool = ismissing.(obshh2)
    if all(bool)
        ww= ones(length(hhset)) ./ length(hhset)
        [hhset ww]
    else
        nmissing = count(bool)
        nerror=[sum(abs.(collect(skipmissing(i .- obshh2)))) for i = hhset2]
        bool = nerror .≤ maxnerror
        n1=nerror[bool]
        n2=(length(obshh)-nmissing) .- n1
        ww = @. (epsf^n1)*((1-epsf)^n2)
        ww ./= sum(ww)
        [hhset[bool] ww]
    end
end

function calfhaploprior(obsf::AbstractVector,epsf::Real,
    isfmalex::AbstractVector,isfounderinbred::Bool,
    founderkind::AbstractString, maxnerror::Integer)
    if isfounderinbred
        # for inbred founders, genotypes are already called  from founder GBS data
        if epsf ≈ 0 || maxnerror == 0
            hh=collect(product([i=="N" ? ["1","2"] : [i] for i = obsf]...))
            hhset=[vcat(i...) for i=reshape(hh,length(hh))]
            ww= ones(length(hhset)) ./ length(hhset)
            fhaplo = [hhset ww]
            fhaplo
        else
            # to check for repeat in hh
            hh=collect(product(Iterators.repeated(["1","2"],length(obsf))...))
            hhset=[vcat(i...) for i=reshape(hh,length(hh))]
            fhaplo=calfweight(hhset,obsf,epsf,isfounderinbred,maxnerror)
        end
        fhaplo
    else
        if (epsf ≈ 0 || maxnerror == 0) && founderkind !== "readcount"
            dict = Dict(["11"=>["11"],"22"=>["22"],
                "12"=>["12","21"],"21"=>["12","21"],
                "1N"=>["11","12","21"],"2N"=>["12", "21","22"],
                "NN"=>["11","12","21","22"]])
            fgeno1 = [get(dict,i,missing) for i= obsf]
        else
            fgeno1=Iterators.repeated(["11","12","21","22"],length(obsf))
        end
        fgeno2=collect(product(fgeno1...))
        fgeno3=reshape(fgeno2,length(fgeno2))
        fgeno4=[vcat(i...) for i=fgeno3]
        if any(isfmalex)
            for i=1:length(fgeno4)
                # For male XY, set Y same as X; genotypic data are for X not Y
                fgeno4[i][isfmalex]=replace(fgeno4[i][isfmalex],"12"=>"11","21"=>"22")
            end
        end
        hhset=[string.(vcat([split(i,"") for i = j]...)) for j=fgeno4]
        if founderkind === "readcount"
            obsww=[[i[1],i[2]/2,i[2]/2,i[3]] for i=obsf]
            obsww2=collect(product(obsww...))
            obsww3=reshape(obsww2,length(obsww2))
            obsww4=prod.(obsww3)
            obsww4=round.(obsww4 ./ sum(obsww4),digits=10)
            b=obsww4 .> 0
            ww = obsww4[b]
            ww ./= sum(ww)
            fhaplo = [hhset[b] ww]
        else
            if epsf ≈ 0 || maxnerror == 0
                ww= ones(length(hhset)) ./ length(hhset)
                fhaplo=[hhset ww]
            else
                obshh=string.(vcat([split(i,"") for i=obsf]...))
                fhaplo=calfweight(hhset,obshh,epsf,isfounderinbred,maxnerror)
            end
        end
        fhaplo
    end
end

# TODO: revise for founderkind of probability vector
function calfhaploprior(magicgeno::MagicGeno,chr::Integer, epsf::Real,baseerror::Real,
    isfounderinbred::Bool,founderkind::AbstractString)
    maxnerror=0
    if founderkind=="readcount"
        if isfounderinbred
            chrfgeno=MagicBase.rawgenocallhaplo.(magicgeno.foundergeno[chr])
        else
            chrfgeno=MagicBase.rawgenoprobdiplo.(magicgeno.foundergeno[chr],baseerror)
        end
    else
        chrfgeno= magicgeno.foundergeno[chr]
    end
    chridls = lowercase.([i[1,:linkagegroup] for i=magicgeno.markermap])
    chrid=chridls[chr]
    isautosome=!(chrid in ["chrx","x"])
    nfounder = size(magicgeno.magicped.founderinfo,1)
    isfmalex =isautosome ? falses(nfounder) : magicgeno.magicped.founderinfo[:gender] .=="male"
    [calfhaploprior(chrfgeno[i,:],epsf,isfmalex,
        isfounderinbred,founderkind,maxnerror) for i =1:size(chrfgeno,1)]
end
