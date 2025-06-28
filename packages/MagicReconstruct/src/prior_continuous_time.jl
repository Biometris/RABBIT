
function getpriorstatespace(magicped::MagicPed;     
    isfounderinbred::Bool=true)
    fgls = string.(magicped.founderinfo[!,:individual])
    isfounderinbred || (fgls = [string(i,"_",j) for i=fgls for j=1:2])
    nfgl=length(fgls)
    haplotypes=fgls
    haploindex=collect(1:nfgl)
    Dict("haplotype"=>haplotypes,"haploindex"=>haploindex)
    # genoindex = [[i,j] for i=1:nfgl for j=1:i]
    # genotypes=permutedims(hcat([fgls[i] for i=genoindex]...))
    # diploindex = [[i,j] for i=1:nfgl for j=1:nfgl]
    # diplotypes=permutedims(hcat([fgls[i] for i=diploindex]...))
    # Dict("haplotype"=>haplotypes,"genotype"=>genotypes,"diplotype"=>diplotypes,
    #     "haploindex"=>haploindex,"genoindex"=>genoindex,"diploindex"=>diploindex)
end

function calpopmakeup(magicgeno::MagicGeno,chr::Integer,
    model::AbstractString,magicprior::NamedTuple;
    isfounderinbred::Bool=true)
    chrid=string(magicgeno.markermap[chr][1,:linkagegroup])
    isautosome=!(chrid in ["chrx"])
    # chrmagicprior: dict of popid=>[nzstate,initiprob,tranrate]
    chrmagicprior = isautosome ? magicprior.autosome : magicprior.allosome
    calpopmakeup(magicgeno.magicped,model, chrmagicprior; isfounderinbred,isautosome)
end

function calpopmakeup(magicped::MagicPed,model::AbstractString,
    chrmagicprior::AbstractDict;    
    isfounderinbred::Bool=true,
    isautosome::Bool=true)
    priorspace=getpriorstatespace(magicped; isfounderinbred)            
    isdepmodel = lowercase(model) == "depmodel"    
    offidls = magicped.offspringinfo[!,:individual]
    founderidls = magicped.founderinfo[!,:individual]
    popmakeup=Dict([begin
        offs=findall(magicped.offspringinfo[!,:member] .== popid)
        ishomozygous = unique(magicped.offspringinfo[offs,:ishomozygous])
        length(ishomozygous)> 1 && @error string("DH and non-DH are mixed in subpoulation = ", popid)
        ishaploid =  isdepmodel || only(ishomozygous)
        if !isautosome
            genderls = unique(offspringinfo[offs,:gender])
            if length(genderls) ==1
                ishaploid = only(genderls) == "male"
            else
                @error string("different genders in subpoulation = ", popid, ": ", genderls)
            end
        end
        fglset = priorspace["haploindex"]
        if ishaploid
            states= fglset
        else
            nfgl = length(fglset)
            states = MagicBase.prior_diploindex(nfgl)
        end
        nzstate,initprob,tranrate= chrmagicprior[popid]
        hashcode=hash([initprob,tranrate])
        founders=sort(unique(reduce(vcat,states[nzstate])))
        founders = isfounderinbred ? founders : unique(div.(founders .+ 1,2))
        # calculate virtualdict
        # maping duplicated id_virtualoffspring to original id (act as founder)
        virtual_offls = string.(offidls[offs])
        keepat!(virtual_offls, occursin.(r"_virtualoffspring$", virtual_offls))
        virtual_fls = replace.(virtual_offls, r"_virtualoffspring$"=>"")
        b = [in(i,founderidls) for i in virtual_fls]
        keepat!(virtual_offls,b)
        keepat!(virtual_fls,b)
        if isempty(virtual_offls)
            virtualdict = Dict{Int,Int}()
        else
            virtual_offls2 = [findfirst(==(i), offidls) for i in virtual_offls]
            virtual_fls2 = [findfirst(==(i), founderidls) for i in virtual_fls]
            virtualdict = Dict(virtual_offls2 .=> virtual_fls2)
            if !issubset(virtual_offls2, offs) 
                @error string("unexpected virtual offspring=", virtual_offls)
            end
        end
        # 
        popid=>Dict(["founder"=>founders, 
            "isfounderinbred"=>isfounderinbred, "fglset"=>fglset, 
            "offspring"=>offs, "ishaploid"=>ishaploid,
            "virtualdict"=>virtualdict,
            "nzstate"=>nzstate, "nzorigin"=>states[nzstate],
            "initprob"=>initprob,"tranrate"=>tranrate,"hashcode"=>hashcode])
    end for popid in keys(chrmagicprior)])
    popmakeup
end

function hmm_nstate_nfgl(popmakeup::AbstractDict)
    popid = first(keys(popmakeup))    
    nfgl = length(popmakeup[popid]["fglset"])
	ishaploidsubls = [i["ishaploid"] for i in values(popmakeup)]
	nstate= all(ishaploidsubls) ? nfgl : nfgl^2
	nstate, nfgl
end

function get_nzcol(nstate::Integer, popmakeup::AbstractDict,popid::AbstractString)
    ishaploid = popmakeup[popid]["ishaploid"]
    nzstate = popmakeup[popid]["nzstate"]
    if ishaploid 
        nfgl = length(popmakeup[popid]["fglset"])
        if nstate == nfgl^2
            nzcol = [(i-1)*nfgl+i for i in nzstate]
        elseif nstate == nfgl
            nzcol = nzstate
        else
            @error string("inconsistent nfgl=",nfgl, ", and nstate=",nstate)
        end
    else
        nzcol = nzstate
    end 
    nzcol
end        

function initprob2tranrate(initprob::AbstractVector,tranrate::AbstractMatrix)
    d = abs.(diag(tranrate))
    tranrate2 = permutedims(reduce(hcat, [begin 
        ls = copy(initprob)            
        ls[i] = zero(eltype(ls))
        normalize!(ls,1)
        ls .*= d[i] 
        ls[i] = -1*d[i]
        ls
    end for i in eachindex(d)]))
    tranrate2
end


function calmagicprior(magicgeno::MagicGeno,model::AbstractString;
    isfounderinbred::Bool=true)
    # nfounder = size(magicgeno.magicped.founderinfo,1)
    # nfgl = isfounderinbred ? nfounder : 2*nfounder
    # smodel=Symbol(lowercase(model))    
    chrls=[string(i[1,:linkagegroup]) for i = magicgeno.markermap]
    isautosomels=[!in(chr,["chrx"]) for chr = chrls]
    if true in isautosomels
        prioraa = calprior(magicgeno.magicped,model;
            isfounderinbred,isautosome=true)        
    end
    if false in isautosomels
        priorxx = calprior(magicgeno.magicped,model;
            isfounderinbred,isautosome=false)        
    end
    if issubset([true,false],isautosomels)
        (autosome=prioraa,allosome=priorxx)
    elseif true in isautosomels
        (autosome=prioraa,)
    else
        (allosome=priorxx,)
    end
end

function calprior(magicped::MagicPed,model::AbstractString;
    isfounderinbred::Bool=true, 
    isautosome::Bool=true)
    # offspringinfo
    memdf = unique(magicped.offspringinfo[!,[:member,:ishomozygous,:isfglexch]])
    memberls = memdf[!,:member]    
    isfglexch_dict = Dict(memberls .=> memdf[!,:isfglexch])
    ishomozygous_dict = Dict(memberls .=> memdf[!,:ishomozygous])
    inputmodel=Symbol(lowercase(model))
    smodel_dict= Dict([mem => (ishomozygous_dict[mem] ? :depmodel : inputmodel) for mem in memberls])      
    # designinfo
    design=magicped.designinfo
    designtype=typeof(design)
    isconcise = lowercase(model)=="jointmodel" ? false : true    
    if isnothing(designtype)
      @error "designinfo is required"
    elseif designtype <: Dict{String,DesignInfo}      
        memdf = unique(magicped.offspringinfo[!,[:member,:ishomozygous,:isfglexch]])
        allfounders = magicped.founderinfo[!,:individual]                  
        Dict([begin 
            subfounders = subdesign.founders            
            subped = subdesign.pedigree            
            if isnothing(subped)
                juncdist = subdesign.juncdist            
                isnothing(juncdist) && @error string("TODO for subdesign=",subdesign)                
                (; ibd,j1122,j1211,j1213,j1222,j1232) = juncdist                        
                phi12 = 1-ibd
                ismalex = false
                fglindicators = get_fglindicators(allfounders,subfounders;isfounderinbred)
                dict=MagicPrior.identityprior(fglindicators,ismalex,phi12,j1122,j1211,j1213,
                    j1222,j1232,isconcise)
            else
                in(popid,subped.member) || @error string("popid =",popid, " is not in pedigree members=",ped.member)            
                dict = MagicPrior.magicsubprior(subped,allfounders; member = popid, isfounderinbred, 
                    isautosome,isfglexch=isfglexch_dict[popid], isconcise)         
            end            
            val = only(values(dict))
            smodel = smodel_dict[popid]
            popid => tranformmarkov(val[smodel]...)
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
        Dict([begin 
            prior = priorls[smodel_dict[mem]]
            mem => tranformmarkov(prior...)
        end for (mem,priorls) in res])
    else
        @error "TODO"
    end
end

function get_fglindicators(allfounders::AbstractVector,founders::AbstractVector;
    isfounderinbred::Bool=true)
    issubset(founders,allfounders) || @error string("founders=",founders, ", not a subset of ",allfounders)
    nfounder = length(allfounders)
    dict = Dict(allfounders .=> 1:nfounder)
    founder_indices = [dict[i] for i in founders]
    if isfounderinbred
        fglindicators = falses(nfounder)        
        fglindicators[founder_indices] .= true
    else
        fglindicators = falses(2*nfounder)
        fgl_indices = reduce(vcat,[[2*i-1,2*i] for i in founder_indices])
        fglindicators[fgl_indices] .= true
    end
    fglindicators
end

function tranformmarkov(initprob::AbstractVector,rate::AbstractMatrix)
    length(initprob) == size(rate,1) == size(rate,2) || @error string("inconsistent dimensions!")
    zstate = initprob .== 0.0
    nzstate= .!zstate
    atol = sum(zstate)*length(initprob)*eps(Float64)
    if !isapprox(sum(rate[zstate,:]), 0; atol)
        @error("inconsistent initial and rate matrix")
    end
    if !isapprox(sum(rate[:,zstate]), 0; atol)
        @error("inconsistent initial and rate matrix")
    end
    findall(nzstate), Vector(initprob[nzstate]), Matrix(rate[nzstate,nzstate])
end

# """
#     pedfile_design2ped(pedfile; commentstring="##",workdir=pwd())

# convert a pedfile from `designinfo` being `designcode` to `pedigree`.

# # Keyword arguments

# `commentstring::AbstractString="##"`: the lines beginning with commentstring are ignored
# in pedfile.

# `workdir::AbstractString=pwd()`: directory for reading pedfile.

# """
# function pedfile_design2ped(pedfile::AbstractString;
#     outfile = "oustem.csv", 
#     commentstring::AbstractString="##",
#     workdir::AbstractString=pwd())
#     magicped = MagicBase.readmagicped(pedfile;commentstring,workdir)
#     MagicBase.savemagicped(MagicBase.getabsfile(workdir,outfile),magicped)
# end
