
"""
    JuncDist

struct that stores junction information, that is, the prior distribution of recombination breakpoints. 

"""
Base.@kwdef struct JuncDist  
    nfounder::Integer      
    ibd::Union{Nothing,Float64} = nothing        
    mapexpansion::Union{Nothing, Float64} = nothing    
    # j1112=j1211
    # j1121=j1222
    # Rm=j1122+j1121+j1222+j1232
    # Rp=j1122+j1121+j1211+j1213
    j1122::Union{Nothing,Float64} = nothing    
    j1211::Union{Nothing,Float64} = nothing    
    j1213::Union{Nothing,Float64} = nothing    
    j1222::Union{Nothing,Float64} = nothing    
    j1232::Union{Nothing,Float64} = nothing             
end

function check_juncdist(juncdist::JuncDist)    
    (;ibd, mapexpansion,j1122,j1211,j1213,j1222,j1232) = juncdist
    if !isnothing(ibd)        
        1>=ibd >= 0 || @error string("ibd = ", ibd, " is not in [0,1]")
    end    
    jidls = ["j1122","j1211","j1213","j1222","j1232"]
    jls = [j1122,j1211,j1213,j1222,j1232]
    b = isnothing.(jls)                
    if any(b)
        if !all(b)                
            @error string("junction densities ",join(jidls[b],","), " are not provided")
        end        
    else
        isneg = jls .< 0         
        if any(isneg)
            @error string("junction densities ", jidls[isneg], " = ",  jls[isneg], " are negative")                
        end        
        if !isnothing(mapexpansion)
            mapexpansion > 0 || @error string("mapexpansion = ", mapexpansion, " is not positive")
            Rm = j1122+2*j1222+j1232
            Rp = j1122+2*j1211+j1213
            mapexpansion0 = (Rm + Rp)/2            
            if !isapprox(mapexpansion, mapexpansion0; atol = 1e-4)
                @error string("inconsistent between mapexpansion = ",mapexpansion, ", and junction density derived mapexpansion = ", mapexpansion0)                    
            end
        end
    end                
end   

"""
    MateScheme

struct that stores mating scheme information for pedinfo. 

"""
struct MateScheme
    nfounder::Integer
    matings::Vector{String}
    generations::Vector{Int}    
    function MateScheme(nfounder, matings, generations)
        nfounder >=1 || @error string("nfounder=",nfounder, " is not >=1")
        length(matings) == length(generations) || @error string("inconsistent size between matings and generations")
        if "DH" in matings[1:end-1]
            @error string("matings=",matings,", DH must be in the last generation")
        end
        if "DH" == last(matings) && last(generations) !== 1
            @error string("DH repeated generation=",last(generations),",  ≠ 1")
        end

        b = check_mating.(matings)    
        all(b) || @error string("unknown matings: ", matings[.!b], 
            ", neither in [Pairing, DH, Selfing,Sibling,HalfDiallel,FullDiallel,RM1E, RM2E, CPME, MAIE], ",
            " nor match any [RM1NE_n, RM2NE_n, WFNE_n] where n is popsize.")
        b = generations .> 0 
        all(b) || @error string("numbers of genrations for ", matings[.!b], " = ", generations[.!b], " are not positive")        
        new(nfounder, matings, generations)
    end    
end

function check_mating(mating::AbstractString)
    matingsets = ["Pairing","DH", "Selfing","Sibling","HalfDiallel","FullDiallel", "RM1E","RM2E","CPME","MAIE"]
    if mating in matingsets
        true
    elseif occursin(r"^(RM1NE|RM2NE|WFNE)_[1-9][0-9]{0,}",mating)  # with population size
        true
    else
        false
    end
end

"""
    DesignInfo

mutable struct that stores design information for a subpopulation. See also [`parsedesign`](@ref). 

# Fields    

`designtype::Symbol`: type of the subpopulation design. It must be :commoncross, :breedcross, :juncdist, or :matescheme. 

`founders::Union{Nothing,AbstractVector}`: founders for the subpopulation. 

`designcode::Union{Nothing,AbstractString}`: string code for the design. 

`pedigree::Union{Nothing, Pedigree}`: pedigree for the design 

`matescheme::Union{Nothing, MateScheme}`: mate schemes for the design

`juncdist::Union{Nothing, JuncDist}`: junctdist for the design. 

"""
mutable struct DesignInfo        
    designtype::Symbol 
    founders::Union{Nothing,AbstractVector}
    designcode::Union{Nothing,AbstractString}
    pedigree::Union{Nothing, Pedigree}
    matescheme::Union{Nothing, MateScheme}
    juncdist::Union{Nothing, JuncDist}
    function DesignInfo(designtype,founders,designcode,pedigree,matescheme,juncdist)   
        if !isnothing(founders)
            nfounder = length(founders)
            if designtype in [:commoncross,:breedcross]
                isnothing(pedigree) && @error string("pedigree is not specified for commoncross design")
                nfounder  == pedigree.nfounder || @error string("inconsistent nfounder between founders and pedigree")
                pedigree.member[1:nfounder] == founders || @error string("inconsistent between founders and pedigree")  
            elseif designtype == :juncdist
                isnothing(juncdist) && @error string("juncdist is not specified for juncdist design")
                nfounder  == juncdist.nfounder || @error string("inconsistent nfounder between founders and juncdist")
            elseif designtype == :matescheme
                isnothing(matescheme) && @error string("matescheme is not specified for matescheme design")
                nfounder  == matescheme.nfounder || @error string("inconsistent nfounder between founders and matescheme")
            else
                @error string("unknown designtype=",designtype, "; designtype must be in [:commoncross, :breedcross, :matescheme, :juncdist]")
            end
        end
        new(designtype,founders,designcode,pedigree,matescheme,juncdist)
    end
end

function DesignInfo(;
    designtype::Symbol,
    founders::Union{Nothing,AbstractVector} = nothing,
    designcode::Union{Nothing,AbstractString} = nothing,
    pedigree::Union{Nothing, Pedigree} = nothing,     
    matescheme::Union{Nothing, MateScheme} = nothing,       
    juncdist::Union{Nothing, JuncDist} = nothing)
    DesignInfo(designtype,founders,designcode,pedigree,matescheme,juncdist)
end

function getnfounder(designinfo::DesignInfo)
    if isnothing(designinfo.founders)        
        isnothing(designinfo.pedigree) || return designinfo.pedigree.nfounder
        isnothing(designinfo.juncdist) || return designinfo.juncdist.nfounder
        isnothing(designinfo.matescheme) || return designinfo.matescheme.nfounder        
    else
        return length(designinfo.founders)
    end    
end

function setfounders!(designinfo::DesignInfo,founders::Union{Nothing,AbstractVector})
    isnothing(founders) && return designinfo
    nfounder = getnfounder(designinfo)    
    if nfounder != length(founders)
        @error "inconsistent nfounder"    
    end
    designinfo.founders = copy(founders)
    if !isnothing(designinfo.pedigree)
        designinfo.pedigree = Pedigrees.setfounderid(designinfo.pedigree,founders)
    end
    designinfo    
end


"""
    parsedesign(designcode; kwargs...)

parse string designcode into [`DesignInfo`](@ref)

# Keyword Argument

`founders::Union{Nothing,AbstractVector}=nothing`: a list of founders. 

`popid="pop"`: population id. 

`fixed_nself::Integer=20`: number of selfing generation for designcode in form of pedcode=>FIXED. 

"""
function parsedesign(designcode::AbstractString;
    founders::Union{Nothing,AbstractVector}=nothing,
    popid="pop",
    fixed_nself::Integer=20)            
    designtype = parse_designtype(designcode)
    isnothing(designtype) && return nothing
    nfounder = isnothing(founders) ? nothing : length(founders)
    if designtype == :commoncross
        res = parse_commoncross(designcode;nfounder,popid)
        setfounders!(res,founders)
    elseif designtype == :breedcross
        res = parse_breedcross(designcode; fixed_nself, popid)
    elseif designtype == :juncdist        
        res = parse_juncdist(designcode;nfounder)               
        setfounders!(res,founders)              
    elseif designtype == :matescheme
        msg = string("TODO for matescheme designcode=",designcode)        
        error(msg) 
        # res = parse_matescheme(designcode;nfounder)     
        return nothing       
    else
        msg = string("Could not parse designcode=",designcode, "; designtype=",designtype)      
        error(msg)
        return nothing         
    end        
    res
end

function parse_designtype(designcode::AbstractString)    
    # characters/substrings are not allowed in parent IDs: /, =, =>, 
    designcode = strip(designcode)
    if occursin(r"^bc[1-9][0-9]{0,6}-",designcode) || 
        occursin(r"^[1-9][0-9]{0,6}ril-",designcode) || 
        occursin(r"^[1-9][0-9]{0,6}star-",designcode)
        res = :commoncross
    elseif occursin("/",designcode) && occursin("=>",designcode) && !occursin("||",designcode) 
        ls = strip.(split(designcode,"=>"))
        length(ls) == 2 || @error string("Could not parse designcode=",designcode, " for breedcross; too many separators =>")        
        res = :breedcross
    elseif  occursin("=>",designcode) 
        matingnames = ["nfounder", "Pairing","DH", "Selfing","Sibling","HalfDiallel","FullDiallel", "RM1E","RM2E","CPME","MAIE"]            
        namels = [string(strip(first(split(i,"_")))) for i in split(designcode,"||")] # remove popsize in the namels
        d = setdiff(namels,matingnames)
        if isempty(d)
            res = :matemcheme                                
        else
            if length(d) < length(namels)
                @error string("Could not parse designcode=",designcode, " for matescheme. Unknown matings = ", d)
            end      
            @info string("parsed fieldnames=",namels, " for designcode=",designcode)              
            return nothing
        end
    elseif occursin("=",designcode)         
        namels = [string(strip(first(split(i,"=")))) for i in split(designcode,"||")]
        juncnames = ["nfounder", "ibd", "mapexpansion", "j1122","j1211","j1213","j1222","j1232"]
        d = setdiff(namels,juncnames)
        if isempty(d)
            res = :juncdist           
        else
            if length(d) < length(namels)
                @error string("Could not parse designcode=",designcode, 
                    " for juncdist. Unknown junctions = ", d, ", fieldnames=",namels)
            end
            @info string("parsed fieldnames=",namels, " for designcode=",designcode)
            return nothing
        end
    else
        occursin("=",designcode)  && @warn string("designcode for juncidst must contain nfounder")
        @error string("Could not parse designcode=",designcode)        
        return nothing
    end
    res 
end


function parse_juncdist(designcode::AbstractString;nfounder::Union{Nothing,Integer}=nothing)    
    juncdist1 = split(designcode,"||")
    juncdist2 = [strip.(split(i,"=")) for i in juncdist1]
    all(length.(juncdist2) .== 2) || @error string("Could not parse juncdist=",juncdist1,"; each field must be separated by "||" and specified by =")
    juncdist3 = [i[1]=>parse(Float64,i[2]) for i in juncdist2]
    b = first.(juncdist3) .== "nfounder"
    if any(b)
        nfounderls = last.(juncdist3[b])
        if !isnothing(nfounder)
            only(nfounderls) == nfounder || @error string("nfounder=",nfounder, ", inconsistent with designcode=",designcode)        
        end
    else
        if !isnothing(nfounder)
            pushfirst!(juncdist3,"nfounder"=>nfounder)
        end
    end
    juncdist = JuncDist(; (Symbol(k) => v for (k,v) in juncdist3)...)
    check_juncdist(juncdist)    
    DesignInfo(; designtype=:juncdist, designcode, juncdist)
end

# see also examples of R/qtl2: https://kbroman.org/qtl2/assets/vignettes/developer_guide.html
"""
    parse_commoncross(designcode::AbstractString)

return a pedigree::Pedigree from a designcode. See [`Pedigree`](@ref) for Pedigree
struct.

`designcode` has four possible formats

* "nril-selfm" where "nril" denotes n-way recombination inbred lines (RIL), and
   "selfm" denotes m generations of inbreeding by selfing, e.g. 8ril-self4.

* "nril-sibm" where "nril" denotes n-way RIL, and "sibm" denoes m generations
  of inbreeding by sibling, e.g. 8ril-sib20.

* "bcg-selfm" where "bcg" denotes g generations of back cross with the first parent,
  and "selfm" denoes m generations of inbreeding by selfing, e.g. bc2-self6.
  "bcg" is the same as "bcg-self0", e.g. "bc2".

* "nstar-selfm" where "nstar" denotes star-like crosses between 1st parent and
  each of the other parents, and "selfm" denotes m generations of inbreeding by
  selfing, e.g. 11star-self6.

"""
function parse_commoncross(designcode::AbstractString; 
    nfounder::Union{Nothing,Integer}=nothing,
    popid=nothing)    
    designcode = strip(designcode)
    inputnfounder = nfounder
    if occursin(r"^[1-9][0-9]{0,6}ril-self[0-9]{1,6}",designcode)
        nfounder,nself=parse.(Int,split(replace(replace(designcode,"ril"=>""),"self"=>""),"-"))        
        p=Int(log2(nfounder))
        if p % 1 != 0.0
            error("wrong number of founders for rils")
        end
        n = sum(2^i for i=0:p) + nself
        mem = collect(1:n)
        m=Int((2-1/2^(p-1))*nfounder)
        mot=[zeros(Int,nfounder); 1:2:m; m+1:m+nself]
        fat=[zeros(Int,nfounder); 2:2:m; m+1:m+nself]
        if !isnothing(popid)
            mem = Vector{Any}(mem)
            mem[end] = popid
        end
        df=DataFrame([mem,mot,fat],Symbol.(["member","mother","father"]))
        pedigree=Pedigree(df)
    elseif occursin(r"^[1-9][0-9]{0,6}ril-sib[0-9]{1,6}",designcode)
        nfounder,nsib=parse.(Int,split(replace(replace(designcode,"ril"=>""),"sib"=>""),"-"))
        p=Int(log2(nfounder)-1)
        if p % 1 != 0.0
            @error "wrong number of founders for rils"
        end
        n = sum(2^i for i=1:p+1) + 2*nsib
        mem = collect(1:n)
        gender = repeat(["female","male"],outer=n ÷ 2)
        m=Int((2-1/2^(p-1))*nfounder)
        mot=[zeros(Int,nfounder); 1:2:m; repeat(m+1:2:m+2*nsib,inner=2)]
        fat=[zeros(Int,nfounder); 2:2:m; repeat(m+2:2:m+2*nsib,inner=2)]
        if !isnothing(popid)            
            mem = Vector{Any}(mem)
            mem[end-1:end] .= [string(popid,"_",i) for i in 1:2]
        end
        df=DataFrame([mem,mot,fat,gender],Symbol.(["member","mother","father","gender"]))
        pedigree=Pedigree(df)
    elseif occursin(r"^bc[1-9][0-9]{0,6}",designcode)
        if occursin(r"^bc[1-9][0-9]{0,6}$",designcode)
            nbc=parse(Int,replace(designcode,"bc"=>""))
            nfounder = 2
            nself=0
        elseif occursin(r"^bc[1-9][0-9]{0,6}-self[0-9]{1,6}",designcode)
            nfounder = 2
            nbc,nself=parse.(Int,split(replace(replace(designcode,"bc"=>""),"self"=>""),"-"))
        else
            error(string("unknown population design: ", designcode))
        end
        # @info "back cross with the first parent"
        n=3+nbc+nself
        mem = collect(1:n)
        mot=[zeros(Int,2); repeat([1],nbc+1); nbc+3:nbc+2+nself]
        fat=[zeros(Int,2); 2:nbc+2+nself]
        if !isnothing(popid)
            mem = Vector{Any}(mem)
            mem[end] = popid
        end
        df=DataFrame([mem,mot,fat],Symbol.(["member","mother","father"]))
        pedigree=Pedigree(df)
    elseif occursin(r"^[1-9][0-9]{0,6}star-self[0-9]{1,6}",designcode)
        nfounder,nself=parse.(Int,split(replace(replace(designcode,"star"=>""),"self"=>""),"-"))
        mot=[zeros(Int,nfounder); ones(Int,nfounder-1)]
        fat=[zeros(Int,nfounder); 2:nfounder]
        generation=[zeros(Int,nfounder); ones(Int,nfounder-1)]
        for i=1:nself
            premem = length(mot)-(nfounder-1)
            self = premem+1:premem+(nfounder-1)
            gg = (i+1)*ones(Int, length(self))
            push!(mot,self...)
            push!(fat,self...)
            push!(generation, gg...)
        end
        mem = 1:length(mot)
        if !isnothing(popid)
            posls = findall(generation .== last(generation))
            mem = Vector{Any}(mem)
            mem[posls] .= [string(popid,"_",i) for i in 1:length(posls)]
        end
        df=DataFrame([mem,mot,fat],Symbol.(["member","mother","father"]))
        pedigree=Pedigree(df)
    else
        error(string("unknown designcode for pedinfo: ", designcode))
    end
    if !isnothing(inputnfounder) 
        nfounder == inputnfounder || @error string("nfounder=",inputnfounder, ", inconsistent with designcode=",designcode)        
    end
    pedigree = Pedigrees.setfounderid(pedigree,string.("P", 1:nfounder))
    DesignInfo(; designtype=:commoncross, designcode, pedigree)
end

function parse_breedcross(designcode::AbstractString; 
    fixed_nself::Integer = 20,
    popid=nothing)     
    designcode = strip(designcode)
    peddf = parsebreedcode(designcode; fixed_nself)
    if !isnothing(popid)        
        peddf[end,1] = popid
    end
    pedigree = Pedigree(peddf)
    nf = pedigree.nfounder
    founders = copy(pedigree.member[1:nf])
    DesignInfo(; designtype=:breedcross, designcode, founders, pedigree)
end

function parsebreedcode(breedcode::AbstractString; fixed_nself::Integer = 20)    
    pedcode,ninbred_str = string.(strip.(split(breedcode,"=>")))    
    # set breeding pedigree without selfing
    occursin(";", pedcode) && @error string("; is not allowed in pedcode: ", pedcode)    
    code0 = replace(pedcode,r"/([1-9][0-9]*)/"=>s";\1;")
    code0 = replace(code0,r"/\s*/"=>";2;")
    code0 = replace(code0,"/"=>";1;")
    occursin("/",code0) && @error string("unxpected delim / in transformed pedcode: ",code0, "; character / is not allowed in founder IDs")
    code1 = strip.(split(code0,";"))
    # check parent id
    id = code1[1:2:end]
    id = intersect(id,string.(1:length(code1)))
    isempty(id) || @warn string("unxpected parent ids: ",id)
    # check generation delim
    gen = code1[2:2:end]
    gen_uniq = unique(gen)
    gen_set = string.(1:length(gen_uniq))
    gen_diff = setdiff(gen_set,gen_uniq)
    isempty(gen_diff) || @error string("missing generations: ", gen_diff, ",pedcode=",pedcode)
    # parse from code1 to ped
    # add parents/founders to ped
    ped = [[i,"0","0"] for i in unique(code1[1:2:end])]
    # add offspring to ped
    for gen in gen_set
        while true
            pos = findfirst(==(gen),code1)
            isnothing(pos) && break
            delim = gen == "1" ? "/" : (gen == "2" ? "//" : string("/",gen,"/"))
            mother = code1[pos-1]
            father  = code1[pos+1]            
            ind = string(mother,delim,father)
            if parse(Int,gen) >= 2 && in(ind, code1)
                k = 2
                ind2 = string(ind, "_",k)
                while in(ind2, code1)
                    k += 1
                    ind2 = string(ind, "_",k)
                end
                ind = ind2
            end
            push!(ped,[ind, mother, father])
            code1[pos-1] = ind
            deleteat!(code1,[pos,pos+1])
        end
    end
    # set inbreeding stage for breeding pedigee
    ninbred_parsed = tryparse(Int,ninbred_str)
    ninbred = isnothing(ninbred_parsed) ? ninbred_str : ninbred_parsed
    memcode = ""    
    if isa(ninbred, AbstractString)
        if uppercase(ninbred) == "DH"
            ninbred = 1
            memcode = "_DH"
        elseif uppercase(ninbred) == "FIXED"
            ninbred = fixed_nself
            memcode = "_FIXED"
        else
            @error string("unknow ninbred : ",ninbred,", breedcode=",breedcode)
        end
    else
        ninbred >= 0 || @error string("ninbreed =", ninbred, " is negative")    
    end
    p00 = ped[end][1]
    ped[end][1] = string(p00, "_self0")
    for i in 1:ninbred
        p = ped[end][1]
        o = string(p00,"_self",i)
        push!(ped, [o, p,p])
    end    
    ped[end][1] = string(ped[end][1], memcode)
    df =DataFrame(permutedims(reduce(hcat,ped)),[:member,:mother, :father])
    unique(df)
end

function todesigncode(juncdist::JuncDist)
    res = ""
    for id in propertynames(juncdist)
        val = getfield(juncdist,id)
        if !isnothing(val) 
            isempty(res) || (res *= "||")
            res *= string(id,"=",val)
        end
    end
    res
end
