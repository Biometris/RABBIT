function rand_magicped(matescheme::MateScheme; popsize::Union{Nothing,Integer}=nothing)    
    designinfo = rand_mateped(matescheme)
    ishomozygous = last(matescheme.matings) == "DH"
    formmagicped(designinfo,popsize;ishomozygous)
end

function rand_mateped(matescheme::MateScheme)    
    popls = rand_popls(matescheme)
    popls2ped(popls)
end

function popls2ped(popls::AbstractVector)
    popls2 = permutedims(reduce(hcat,reduce(vcat,popls)))
    df = DataFrame(popls2[:,2:end],[:founder, :mother, :father])
    Pedigree(df)
end

function rand_popls(matescheme::MateScheme)
    popls = init_popls(matescheme)
    for i in eachindex(matescheme.matings)
        mating = matescheme.matings[i]
        ng = matescheme.generations[i]
        if mating == "Pairing"
            n = length(last(popls))
            ngmax = log(2,n)
            rem(ngmax,1) == 0 || @error string("current popsize=",n, ", is not a power of 2")
            ng <= ngmax || @error string("repeat #geneations =",ng, ", is not <= log(2,",n,")")           
        end
        for _ in 1:ng
            rand_popnext!(popls,mating)
            del_non_ancestor!(popls)
        end
    end
    popls
end

function init_popls(matescheme::MateScheme)
    nf = matescheme.nfounder
    ng = 0
    popfounders = [[ng, string("P",i), 0,0] for i in 1:nf]
    [popfounders]
end

function del_non_ancestor!(popls::AbstractVector)
    ng = length(popls)    
    for g in ng:-1:2
        parents_next = unique(reduce(vcat, [i[[3,4]] for i in popls[g]]))
        b = .![i[2] in parents_next for i in popls[g-1]]
        if any(b)
            deleteat!(popls[g-1],b)
        else
            break
        end
    end
    popls
end

function rand_popnext!(popls::AbstractVector,mating::AbstractString)    
    popnow = last(popls)
    n = length(popnow)
    parents = rand_mating_parents(n, mating)
    ng = popnow[1][1]+1    
    idls = [i[2] for i in popnow]
    popnext = [[ng, string("mem_",ng,"_",i), idls[parents[i]]...] for i in eachindex(parents)]
    push!(popls,popnext)    
end

function rand_mating_parents(n::Integer,mating::AbstractString)    
    # n = number of parents in current population
    # simulate parents of each offspring in the next generation    
    if mating == "Pairing"
        iseven(n) || @warn string("current popsize=",n,"; it must be even for ",mating)
        [[i,i+1] for i in 1:2:n]
    elseif mating == "DH" 
        # DH is same as Selfing, except that DH is alway in the last genration, and ishomozygous = true for its offspring
        [[i,i] for i in 1:n]
    elseif mating == "Selfing" 
        [[i,i] for i in 1:n]    
    elseif mating == "Sibling"
        iseven(n) || @warn string("current popsize=",n,"; it must be even for ",mating)
        [[i,i+1] for i in 1:2:n for _ in 1:2]
    elseif mating == "HalfDiallel"
        [[i,j] for i in 1:n-1 for j in i+1:n]
    elseif mating == "FullDiallel"
        res = [[i,j] for i in 1:n for j in 1:n]
        res[allunique.(res)]
    elseif mating == "RM1E"        
        ls = rand_derangement(n)
        res = [[i,ls[i]] for i in 1:n]
        shuffle!(res)        
    elseif mating == "RM2E"
        mother = shuffle(collect(1:2:n))
        father = shuffle(collect(2:2:n))
        res = collect.(Pair.(mother, father))
        shuffle(vcat(res,res))
    elseif mating == "CPME"
        res = [[i,i+1] for i in 1:2:n]
        res = vcat(res,res)
        circshift(res,1)
    elseif mating == "MAIE"
        res = [[i,i+1] for i in 1:2:n]
        vcat(res,res)
    elseif occursin(r"^WFNE_[1-9][0-9]{0,}",mating) 
        popsize = tryparse(Int,replace(mating,"WFNE_"=>""))
        isnothing(popsize) && @error string("cannot extract popsize from ",mating)
        [sample(1:n,2,replace=true) for _ in 1:popsize]
    elseif occursin(r"^RM1NE_[1-9][0-9]{0,}",mating) 
        popsize = tryparse(Int,replace(mating,"RM1NE_"=>""))
        isnothing(popsize) && @error string("cannot extract popsize from ",mating)
        [sample(1:n,2,replace=false) for _ in 1:popsize]
    elseif occursin(r"^RM2NE_[1-9][0-9]{0,}",mating) 
        popsize = tryparse(Int,replace(mating,"RM2NE_"=>""))
        isnothing(popsize) && @error string("cannot extract popsize from ",mating)
        iseven(n) || @warn string("current popsize=",n,"; it must be even for ",mating)
        nmate = round(Int,popsize/2)
        mother = sample(1:2:n,nmate,replace=true)
        father = sample(2:2:n,nmate,replace=true)
        mates = collect.(Pair.(mother, father))
        shuffle(vcat(mates,mates))
    else
        @error string("unknown mating=",mating)
    end    
end

function rand_derangement(n::Integer)
    res = randperm(n)
    while any(res .== 1:n)
        res = randperm(n)
    end
    res
end

