"""
    Pedigree{T}

immutable struct that stores information of a pedigree. See also [`readped`](@ref).

The Pedigree fields include `nfounder`, `generation`, `member`, `gender`,
`mother`, and `father`. The mothers and fathers of founders must be denoted by 0.

    Pedigree(nfounder,member,mother,father,gender,generation)

inner constructor. The gender of each member must be
"notapplicable", "female", or "male".

    Pedigree(df::DataFrame)

extenral constructor. The dataframe must have columns `member`,
`mother`, and `father`. A gender column is optional. The constructor calculates
the generation for each member, which is the max path length from founders to
the member. The constructor orders pedigree by calling [`orderped`](@ref).

"""
struct Pedigree{T}
    nfounder::Int
    member::Vector{T}
    mother::Vector{T}
    father::Vector{T}
    gender::Vector{T}
    generation::Vector{Int}
    function Pedigree(nfounder,member,mother,father,gender,generation)
        if !allunique(member)
            error("pedigree members are not unique. ", member)
        end
        new{typeof(first(member))}(nfounder,member,mother,father,gender,generation)
    end
end

function Pedigree(df::DataFrame)
    colls = propertynames(df)
    if issubset([:member, :mother,:father],colls)
        if in(:gender, colls)
            df = df[:,[:member, :mother,:father,:gender]]
        else
            df = df[:,[:member, :mother,:father]]
        end
    else
        @error string("pedigree dataframe must have columns member, mother, father")
    end
    n=size(df,1)
    if :gender in propertynames(df)
        gender=string.(df[!,:gender])
        unknown=setdiff(union(gender),["notapplicable","female","male"])
        if !isempty(unknown)
            @error("unknown genders = ", unknown, "; gender must be female, male, or notapplicable")
        end
    else
        gender =repeat(["notapplicable"],n)
    end
    member,mother,father=[String.(strip.(string.(df[!,i]))) for i=[:member, :mother,:father]]
    allunique(member) || @error("pedigree members are not unique")
    ism0= mother .== "0"
    isf0= father .== "0"
    isfounder= ism0 .* isf0
    isnonfounder = .!isfounder
    d=setdiff(mother[isnonfounder],member)
    isempty(d) || @error("nonfounders have unknown mothers  ", d)
    d=setdiff(father[isnonfounder],member)
    isempty(d) || @error("nonfounders have unknown father  ", d)
    nfounder=Int(sum(isfounder))
    generation=Int.(isnonfounder)
    ped=Pedigree(nfounder,member,mother,father,gender,generation)
    orderped(ped)
end

"""
    orderped(ped::Pedigree)

order pedigree members so that parents always come first.

"""
function orderped(ped::Pedigree)
    member, mother, father = ped.member, ped.mother, ped.father
    gen = [findall(ped.generation .== 0)]
    n = length(member)
    while true
        indexbef =vcat(gen...)
        indexaft=setdiff(1:n,indexbef)
        indexaft ==[] && break
        popbef=member[indexbef]
        b=[(mother[i] in popbef) && (father[i] in popbef)  for i= indexaft]
        push!(gen,indexaft[b])
    end
    ii=vcat(gen...)
    generation=vcat([(i-1)*ones(Int,length(gen[i])) for i=1:size(gen,1)]...)
    Pedigree(ped.nfounder,member[ii],mother[ii],father[ii],ped.gender[ii],generation)
end

"""
    readped(pedfile::AbstractString, delim=',', commentstring="##")

read a pedigree from pedfile, ignoring the lines beginning with commentstring.

The pedigree file  must contain columns: member, mother, and father. The mothers
and fathers of founders must be denoted by 0. The value in the optional gender
column must  be "notapplicable", "female" or "male".

"""
function readped(pedfile::AbstractString, delim=',', commentstring="##")
    df=CSV.File(pedfile, delim=delim,comment=string(commentstring)) |> DataFrame
    Pedigree(df)
end

"""
    saveped(filename::AbstractString,ped::Pedigree)

save pedigree into a csv file. The file contains columns: member,
mother, father, gender, generation.

"""
function saveped(filename::AbstractString,ped::Pedigree)
    df = ped2df(ped)
    CSV.write(filename,df)
end

function ped2df(ped::Pedigree)    
    # first field is nfounder, a scalar, not vector
    title=[fieldnames(typeof(ped))[2:end]...]
    df=DataFrame([getfield(ped,i) for i=title],title)
    df
end

"""
    setfounderid(ped::Pedigree,founders::AbstractVector)

change founder IDs in a pedigree

"""
function setfounderid(ped::Pedigree{String},founders::AbstractVector)
    isfounder =  @. (ped.mother =="0") * (ped.father == "0")
    if sum(isfounder) != ped.nfounder
        error("founders does not match nfounder")
    end
    if length(founders) != ped.nfounder
        error("the number of input founders must be  ",ped.nfounder)
    end
    if !allunique(founders)
        error("input founders are not unique: ",founders)
    end
    founders=string.(founders)
    ls = intersect(founders,ped.member[.!isfounder])
    if !isempty(ls)
        error("input founders overlap with non-founders: ",ls)
    end
    old= ped.member[isfounder]
    fdict = Dict(map((x,y)->x => convert(String,y),old, founders))
    member=[get(fdict,i,i) for i in ped.member]
    mother=[get(fdict,i,i) for i in ped.mother]
    father=[get(fdict,i,i) for i in ped.father]
    Pedigree(ped.nfounder,member,mother,father,ped.gender,ped.generation)
end

function getfounderid(ped::Pedigree)
    nf = ped.nfounder
    ped.member[1:nf]
end

"""
    toindexped(ped::Pedigree)

return a pedigree where members are coded as integers staring from 1.
The input pedigree must be ordered. See [`orderped`](@ref).

"""
function toindexped(ped::Pedigree)
    genderset = union(ped.gender)
    gendertype=promote_type(typeof.(genderset)...)
    if gendertype <: AbstractString
        unknown=setdiff(genderset,["notapplicable","female","male"])
        if isempty(unknown)
            sexrule=Dict(["notapplicable"=>0,"female"=>1,"male"=>2])
            genderno=[get(sexrule,i,-1) for i = ped.gender]
        else
            error("unknown genders = ", unknown, "; gender must be notapplicable, female, or male")
        end
    elseif gendertype <: Integer
        unknown=setdiff(genderset,[0,1,2])
        if isempty(unknown)
            genderno= ped.gender
        else
            error("unknown genders = ", unknown, "; gender must be  0=notapplicable, 1=female or 2=male")
        end
    else
        error("unexpected type of gender: ", gendertype)
    end
    member=ped.member
    mother=ped.mother
    father=ped.father
    dict=Dict(member .=> collect(1:size(member,1)))
    memberno= [get(dict,i,0) for i in member]
    motherno= [get(dict,i,0) for i in mother]
    fatherno= [get(dict,i,0) for i in father]
    Pedigree(ped.nfounder,memberno,motherno,fatherno,genderno,ped.generation)
end


function getparents(peddict::Dict,memls::AbstractVector)
    res = memls
    while true
        memls = unique(vcat([get(peddict,i,nothing) for i=memls]...))
        memls = memls[.!isnothing.(memls)]
        isempty(memls) && break
        push!(res,memls...)
    end
    res
end

function getsubped(pedigree::Pedigree{T},member::T) where {T <: Union{String,Integer}}
    nlast = findlast(x->x==member,pedigree.member)
    isnothing(nlast) && error(string(member, " is not in the pedigree"))
    nf = pedigree.nfounder
    dict = Dict([pedigree.member[i] => [pedigree.mother[i],pedigree.father[i]] for i=nlast:-1:nf+1])
    parents = getparents(dict,[member])
    b = [in(i,parents) for i=pedigree.member]
    mother, father = pedigree.mother[b],pedigree.father[b]
    gender, generation = pedigree.gender[b], pedigree.generation[b]
    # cculcate new nf
    isintped = eltype(member) <: Integer
    typezero = isintped ? 0 : "0"
    ism0= mother .==  typezero
    isf0= father .== typezero
    isfounder= ism0 .* isf0
    Pedigree(sum(isfounder),pedigree.member[b],mother,father,gender,generation)
end

function getsubped(pedigree::Pedigree{T},memberls::Vector{T}) where {T <: Union{String,Integer}}
    issubset(memberls, pedigree.member) || @error string("members not in pedigree: ",
        setdiff(memberls,pedigree.member))
    nf = pedigree.nfounder
    nlast = length(pedigree.member)
    dict = Dict([pedigree.member[i] => [pedigree.mother[i],pedigree.father[i]] for i=nlast:-1:nf+1])
    parents = Pedigrees.getparents(dict, memberls)
    b = [in(i,parents) for i in pedigree.member]
    mother, father = pedigree.mother[b],pedigree.father[b]
    gender, generation = pedigree.gender[b], pedigree.generation[b]
    member = pedigree.member[b]
    # cculcate new nf
    isintped = eltype(member) <: Integer
    typezero = isintped ? 0 : "0"
    ism0= mother .==  typezero
    isf0= father .== typezero
    isfounder= ism0 .* isf0
    Pedigree(sum(isfounder),member,mother,father,gender,generation)
end

function pedadjmtx(pedigree::Pedigree)
    iped=Pedigrees.toindexped(pedigree)
    n = length(iped.member)
    adj = zeros(Int,n,n)
    nf=iped.nfounder
    for i=nf+1:n
        adj[i, iped.mother[i]] = 1
        adj[i, iped.father[i]] = 1
    end
    for i=1:n
        adj[1:i-1,i] = adj[i,1:i-1]
    end
    adj
end

function splitvec(f::Function,A::AbstractVector)
    res=Vector{typeof(A)}()
    i0=1
    for i=2:size(A,1)
        if !f(A[i-1],A[i])
            push!(res,A[i0:i-1])
            i0=i
        end
    end
    push!(res,A[i0:end])
    res
end

function splitvec(A::AbstractVector)
    f(x,y)= x == y
    splitvec(f,A)
end



"""
    plot(ped::Pedigree)

plot a pedigree. 

"""
function plotped(pedigree::Pedigree;
    curves = true,
    nodename = nothing,
    nodecolor = nothing,
    plotsize=nothing,
    markersize = nothing,
    fontsize = 10,
    nodesize = 0.1,
    edgecolor = :black,
    edgestyle = :solid,    
    plotargs... 
    )
    adj = pedadjmtx(pedigree)
    n=size(adj,1)
    if issubset(["female","male"], unique(pedigree.gender))
        nodeshape = [i=="male" ? :rect : :circle for i= pedigree.gender]
    else
        nodeshape = :circle
    end
    nf  = pedigree.nfounder
    if isnothing(nodecolor)    
        founder_colors = distinguishable_colors(nf, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
        nodecolor = vcat(founder_colors,[RGB(0.9,0.9,0.9) for i=nf+1:n])
    end
    nodename = isnothing(nodename) ? string.(pedigree.member) : string.(nodename)    
    xls, yls = get_ped_xy(pedigree)    
    xmin, xmax = extrema(xls)
    if isnothing(plotsize)
        plotsize = (120*(xmax-xmin),100*first(yls))
    end

    if isnothing(markersize)
        lenls = length.(splitvec(pedigree.generation))
        maxlen = maximum(lenls)
        markersize = min(maxlen<=4 ? 2 : 5, max(length(lenls),maxlen)/3.0)
    end
    fig = graphplot(adj;
        x=xls, y=yls,
        nodeshape,
        nodecolor,        
        curves,
        markersize, 
        nodesize, 
        edgecolor, 
        edgestyle,
        size = plotsize, 
        plotargs...
    )
    labls = [(xls[i],yls[i],Plots.text(nodename[i],fontsize,:black)) for i=eachindex(xls, yls, nodename)]
    annotate!(fig,labls)
    fig
end

function get_ped_xy(ped::Pedigree; minxspace=nothing)
    gen = ped.generation
    gen2 = splitvec(gen)        
    gensize = reduce(vcat,[length(i)*ones(Int,length(i)) for i in gen2])      
    x_gen=reduce(vcat,[begin
        c = length(i)
        xx = (1:c) .- ((c+1)/2)        
    end for i=gen2])
    nf = ped.nfounder
    iped = toindexped(ped)    
    x = zeros(length(iped.member))    
    x[1:nf] .= (1:nf) .- ((nf+1)/2)
    for i in nf+1:length(iped.member)                
        x[i] = (x[iped.mother[i]] + x[iped.father[i]])/2
    end
    b = gensize .> 1
    @. x[b] = (x[b] + x_gen[b])/2    
    maxlen  = max(length.(gen2)...)
    xscale = 1.5*maxlen / length(gen2[1])
    x .*= xscale
    y = (last(gen)+1.0) .- gen    
    xmin,xmax = extrema(x)
    yscale = max(1.15, 1.5*(xmax-xmin)/first(y))
    y .*= yscale
    # enlarge founder populations
    b = gen .== 0
    x[b] .*= 1.1
    y[b] .+= 1.0

    minxspace = isnothing(minxspace) ? xscale : 1
    # set minxspace 
    iils = splitindex(gen)
    for ii in iils        
        subx = view(x, ii)
        oo = sortperm(subx)
        d = diff(subx[oo])
        d .= [i < minxspace ? minxspace : i for i in d]
        d .= accumulate(+, d)
        pushfirst!(d,0)
        start = subx[oo[1]] - (d[end] - (subx[oo[end]] - subx[oo[1]]))/2
        subx[oo] .= d .+ start        
    end

    x, y
end

function splitindex(f::Function,A::AbstractVector)
    size(A,1)==1 && return [1:1]
    res=Vector{typeof(1:1)}()
    i0=1
    for i=2:size(A,1)
        if !f(A[i-1],A[i])
            push!(res,i0:i-1)
            i0=i
        end
    end
    push!(res,i0:size(A,1))
    res
end

function splitindex(A::AbstractVector)
    f(x,y)= if ismissing(x) 
        ismissing(y)
    else
        if ismissing(y)
            false
        else
            x==y
        end
    end        
    splitindex(f,A)
end


function ped_nodename(pedigree::Pedigree)
    isstr_ped = eltype(pedigree.member) <: AbstractString
    if isstr_ped
        nf = pedigree.nfounder
        fls = pedigree.member[1:nf]
        offls = nf+1:length(pedigree.member)
        offls2 = [in(i,fls) ? string("O_",i) : i for i in offls]
        nodename = vcat(fls,offls2)
    else
        nodename = pedigree.member
    end
    nodename
end

function rename_nonfounder(ped::Pedigree;
    suffix::AbstractString=Random.randstring(),
    keepmembers::AbstractVector=[])
    b = ped.generation .> 0
    old=ped.member[b]
    if !isempty(keepmembers)   
        d =setdiff(keepmembers,ped.member)
        isempty(d) || @warn string(d," are not members of pedigree")
    end
    setdiff!(old,keepmembers)
    new = [string(i,"_",suffix) for i in old]
    fdict=Dict(old .=> new)
    member=[get(fdict,i,i) for i in ped.member]
    mother=[get(fdict,i,i) for i in ped.mother]
    father=[get(fdict,i,i) for i in ped.father]
    Pedigree(ped.nfounder,member,mother,father,ped.gender,ped.generation)
end

function mergeped(pedigrees::Vector{Pedigree{String}}; 
    isfounderinbred::Bool=true,
    keeplast::Bool=true)
    pedls = [Pedigrees.rename_nonfounder(ped; 
        suffix=string(i),
        keepmembers = (keeplast ? [ped.member[end]] : [])) for (i,ped) in enumerate(pedigrees)]
    peddf = unique(reduce(vcat,Pedigrees.ped2df.(pedls)))
    sort!(peddf,:generation)
    if isfounderinbred
        f1df = peddf[peddf[!,:generation] .== 1,:]
        dict = Dict()
        for row in eachrow(f1df)            
            push!(dict,row[:member]=>join(row[[:mother,:father]],"/"))
        end
        peddf[!,[:member,:mother,:father]] .= [get(dict,i,i) for i in Matrix(peddf[!,[:member,:mother,:father]])]
    end
    unique!(peddf)
    # sort founders
    b = @. peddf[!,:mother] == "0" && peddf[!,:father]  == "0"
    peddf[b,:] .= sort(peddf[b,:], :member)
    Pedigree(peddf)
end
