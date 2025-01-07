
###################################################################
# ferive = (haplo=fderivehaplo,diplo=fderivediplo)
# fderivehaplo: nsnp x nstate, nstate=nfgl
# fderivediplo: nsnp x nstate, nstate = nfgl^2,
# fderivediplo[i,j] = derived phased genotype_code for marker i and hidden state j
###################################################################
# offcode: nsnp x noffspring
# offcode[i,j] = phased/unphased genotype code for marker i and offspring j
# observed/derived haplo_code ["N"=>1,"1"=>2, "2"=>3]; TODO for outbred male founders
# unphased genotype_code ["NN", "1N"/"N1", "2N"/"N2", "11", "12"/"21", "22"] .=> 1:6
# phased genotype_code ["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"] .=> 1:9
###################################################################
# condlike = (haplo=condlikehaplo,diplo=condlikediplo)
# condlikehaplo: 3 obsvered haplo_code x 3 derived haplo_code
# condlike.diplo = (:nonibd,:ibd): 6x9 for unphased offspring, 9x9 for phased offspring
# condlike.diplo: likelihood of obsered_genotypes given true derived genotypes (9 possible)
###################################################################

function precompute_chr(chrfhaplo::AbstractMatrix,
    chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,    
    isoffphased::Bool,    
    issnpGT::AbstractVector)
    fderive = precompute_fderive(chrfhaplo,popmakeup)
    offcode = precompute_offcode(chroffgeno,popmakeup,issnpGT; isoffphased)    
    fderive, offcode
end


function precompute_chr(chrfhaplo::AbstractMatrix,
    chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
    issnpGT::AbstractVector;
    isoffphased::Bool=false,
    israndallele::Bool)
    fderive = precompute_fderive(chrfhaplo,popmakeup)
    offcode = precompute_offcode(chroffgeno,popmakeup,issnpGT; isoffphased)
    if any(issnpGT)
        condlike = precompute_condike(epsf,epso; isoffphased,israndallele)
    else
        condlike = nothing
    end
    fderive, offcode, condlike
end


function precompute_chr_mma(magicgeno::MagicGeno,chr::Integer,model::AbstractString, isoffphased::Bool,isfounderinbred::Bool)
    chrid = magicgeno.markermap[chr][1,:linkagegroup]
    isautosome=!(lowercase(chrid) in ["chrx"])
    issnpGT = [occursin("GT",i) for i in magicgeno.markermap[chr][!,:offspringformat]]
    ishaploidls = precompute_ishaploidls(magicgeno.magicped.offspringinfo,model,isautosome)
    if isfounderinbred
        chrfhaplo = magicgeno.foundergeno[chr]
    else
        chrfhaplo = permutedims(reduce(hcat,[reduce(vcat,i) for i in eachrow(magicgeno.foundergeno[chr])]))
    end
    fderive = precompute_fderive(chrfhaplo,ishaploidls)
    offcode = precompute_offcode(magicgeno.offspringgeno[chr],ishaploidls,issnpGT; isoffphased)    
    isautosome, issnpGT, ishaploidls, fderive, offcode
end

function precompute_ishaploidls(offspringinfo::DataFrame,model::AbstractString,isautosome::Bool)
    noff=size(offspringinfo,1)
    res = falses(noff)
    memls = offspringinfo[!,:member]
    for popid in memls
        offs=findall(memls .== popid)
        ishomozygous = unique(offspringinfo[offs,:ishomozygous])
        res[offs] .= lowercase(model)=="depmodel" || only(ishomozygous)
        if !isautosome
            genderls = unique(offspringinfo[offs,:gender])
            if length(genderls) ==1
                res[offs] .= only(genderls) == "male"
            else
                @error string("different genders in subpoulation = ", popid, ": ", genderls)
            end
        end
    end
    res
end

function precompute_ishaploidls(popmakeup::AbstractDict)
    noff = sum([length(val["offspring"]) for val in values(popmakeup)])
    ishaploidls = falses(noff)
    for (_,val) in popmakeup
        ishaploidls[val["offspring"]] .= val["ishaploid"]
    end
    ishaploidls
end

function precompute_fderive(chrfhaplo::AbstractMatrix,popmakeup::AbstractDict)
    nzstate= sort(unique(reduce(vcat,[val["ishaploid"] ? [] : val["nzstate"] for val = values(popmakeup)])))
    ishaploidls = unique([val["ishaploid"] for val in values(popmakeup)])
    precompute_fderive(chrfhaplo,ishaploidls;nzstate)
end

function precompute_fderive(chrfhaplo::AbstractMatrix,ishaploidls::AbstractVector;
    nzstate::Union{Nothing,AbstractVector}=nothing)
    if true in ishaploidls
        fderivehaplo=calfderive(chrfhaplo; nzstate=nothing,ishaploid=true)
    else
        fderivehaplo=nothing
    end
    if false in ishaploidls
        fderivediplo=calfderive(chrfhaplo; nzstate,ishaploid=false)
    else
        fderivediplo=nothing
    end
    (haplo=fderivehaplo,diplo=fderivediplo)
end

function precompute_offcode(magicgeno::MagicGeno,chr::Integer,
    popmakeup::AbstractDict; isoffphased::Bool = false)
    issnpGT = [occursin("GT",i) for i in magicgeno.markermap[chr][!,:offspringformat]]
    offcode = precompute_offcode(magicgeno.offspringgeno[chr],popmakeup,issnpGT; isoffphased)
    offcode 
end

function precompute_offcode(chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,issnpGT::AbstractVector;
    isoffphased::Bool = false)
    ishaploidls = precompute_ishaploidls(popmakeup)
    precompute_offcode(chroffgeno,ishaploidls, issnpGT; isoffphased)
end

function precompute_offcode(chroffgeno::AbstractMatrix,
    ishaploidls::AbstractVector, issnpGT::AbstractVector;
    isoffphased::Bool = false)
    #TODO: for outbred phased parents
    isdiploidls = .!ishaploidls
    if all(issnpGT)
        if all(ishaploidls)
            offcode = caloffcode(chroffgeno; ishaploid=true, isoffphased)
        elseif all(isdiploidls)
            offcode = caloffcode(chroffgeno; ishaploid=false,isoffphased)
        else
            offcode = Matrix{Union{Missing,Int8}}(undef,size(chroffgeno))
            offcode[:,ishaploidls] .= caloffcode(view(chroffgeno,:,ishaploidls);
                ishaploid=true, isoffphased)
            offcode[:,isdiploidls] .= caloffcode(view(chroffgeno,:,isdiploidls);
                ishaploid=false,isoffphased)
        end
    else
        offcode = Matrix{Any}(undef,size(chroffgeno))
        offcode[.!issnpGT,:] = chroffgeno[.!issnpGT,:]
        offcode[issnpGT,ishaploidls] .= caloffcode(view(chroffgeno,issnpGT,ishaploidls);
            ishaploid=true, isoffphased)
        offcode[issnpGT,isdiploidls] .= caloffcode(view(chroffgeno,issnpGT,isdiploidls);
            ishaploid=false,isoffphased)
    end    
    offcode
end

function precompute_condike(epsf::Real,epso::Real;
    isoffphased::Bool=false,israndallele::Bool)
    likehaplo=haplolike(epsf,epso)
    likediplo=diplolike(epsf,epso; isoffphased,israndallele)
    (haplo=likehaplo,diplo=likediplo)
end

function precompute_condike(epsf::Real, epsols::AbstractVector;
    isoffphased::Bool=false,israndallele::Bool)
    epsfls = repeat([epsf],length(epsols))
    precompute_condike(epsfls,epsols;isoffphased,israndallele)
end

function precompute_condike(epsfls::AbstractVector, epso::Real;
    isoffphased::Bool=false,israndallele::Bool)
    epsols = repeat([epso],length(epsfls))
    precompute_condike(epsfls,epsols;isoffphased,israndallele)
end

function precompute_condike(epsfls::AbstractVector, epsols::AbstractVector;
    isoffphased::Bool=false,israndallele::Bool)
    length(epsfls) == length(epsols) || @error string("inconsistent length between epsfls and epsols")
    nsnp = length(epsols)
    condlikehaplo= [haplolike(epsfls[i],epsols[i]) for i in 1:nsnp]
    condlikediplo= [diplolike(epsfls[i],epsols[i]; isoffphased,israndallele) for i in 1:nsnp]
    (haplo=condlikehaplo,diplo=condlikediplo)
end

