
function getfindexlist(byfounder::Integer,fmissls::AbstractVector, 
    popmakeup::AbstractDict; defaultby::Integer,minfmiss=1e-4)    
    if byfounder <= -1  || (byfounder == 0 && defaultby <= -1)               
        founders = sortperm(fmissls,rev=true)        
        deleteat!(founders,fmissls[founders] .<= minfmiss)
        res = [founders]
    elseif byfounder == 1 || (byfounder == 0 && defaultby==1)
        founders = sortperm(fmissls,rev=true)        
        deleteat!(founders,fmissls[founders] .<= minfmiss)
        res = [[i] for i in founders]
    else
        res = get_founder_partition(byfounder, fmissls,popmakeup; defaultby,minfmiss)
    end
    res
end


function get_partition_defaultby(popmakeup::AbstractDict, isfounderinbred::Bool, isfounderphased::Bool)
    nfls = [length(i["founder"]) for i=values(popmakeup)]
    maxnf = maximum(nfls)
    if isfounderinbred         
        defautby = maxnf <= 9 ? 8 : 4
    else
        defautby = isfounderphased ? 4 : 2        
    end    
end

function get_founder_partition(byfounder::Integer,fmissls::AbstractVector, 
    popmakeup::AbstractDict;defaultby::Integer=4,minfmiss=1e-4)    
    byfounder >= 0 || @error string("unexpected byfounder=",byfounder)
    fprogenyls = get_fprogenyls(popmakeup)
    founders_nomiss = findall(fmissls .<= minfmiss)
    popfounderls = unique([i["founder"] for i=values(popmakeup)])
    popfounderls = [setdiff(i, founders_nomiss) for i in popfounderls]
    deleteat!(popfounderls,isempty.(popfounderls))
    popls = deepcopy(popfounderls)
    res = Vector{Vector{Int}}()    
    isempty(popls) && return res
    while true             
        fblock = get_fblock(byfounder, popls, fprogenyls; defaultby)    
        push!(res,fblock)        
        for pop in popls
            setdiff!(pop,fblock)
        end        
        deleteat!(popls, isempty.(popls))
        isempty(popls) && break        
    end
    res
end

function get_fprogenyls(popmakeup::AbstractDict)
    pop_ff_offs=[[popmakeup[popid]["founder"],popmakeup[popid]["offspring"]] for popid in keys(popmakeup)]
    fls = sort(unique((reduce(vcat,first.(pop_ff_offs)))))
    fls == 1:length(fls) || "wrong founder indcies"
    fprogenyls = [Vector{Int}() for _ in fls]
    for (ff,offs) in pop_ff_offs
        for f in ff
            append!(fprogenyls[f],offs)
        end
    end
    fprogenyls
end

function get_fblock(byfounder::Integer, popls::AbstractVector,     
    fprogenyls::AbstractVector;
    defaultby::Integer=4)        
    nprogenies = [length(reduce(union,fprogenyls[i])) for i in popls]
    weights = ProbabilityWeights(nprogenies)
    subpop = sample(popls,weights)    
    if byfounder == 0   
        byfounder = defaultby             
    end
    if length(subpop) <= byfounder+1
        fblock = subpop
    else        
        fblock = sample_fblock(subpop,byfounder)
    end        
    sort(fblock)
end

function sample_fblock(subpop::AbstractVector,blocksize::Integer)
    sample(subpop, blocksize;replace=false)            
end


# function get_fblock(byfounder::Integer, popls::AbstractVector, 
#     fmissls::AbstractVector,
#     fprogenyls::AbstractVector;
#     defaultby::Integer=4)
#     fls = unique(reduce(vcat,popls))    
#     f = fls[argmax(fmissls[fls])]    
#     subpops = popls[[f in i for i in popls]]
#     nprogenies = [length(reduce(union,fprogenyls[i])) for i in subpops]
#     subpop = sort(subpops[argmax(nprogenies)])        
#     subpop = subpop[sortperm(fmissls[subpop],rev=true)]    
#     if byfounder == 0   
#         byfounder = defaultby             
#     end
#     if length(subpop) <= byfounder+iseven(byfounder)
#         fblock = subpop
#     else        
#         fblock = sample_fblock(subpop,byfounder)
#     end        
#     sort(fblock)
# end

# function sample_fblock(subpop::AbstractVector,fmiss_subpop::AbstractVector,blocksize::Integer)
#     islargemiss = fmiss_subpop .> 0.9
#     nlargemiss = sum(islargemiss)
#     if nlargemiss == 0 || nlargemiss == length(subpop)            
#         weights = ProbabilityWeights(fmiss_subpop .+ 1e-5)
#         fblock = sample(subpop, blocksize;replace=false)            
#     elseif nlargemiss >= blocksize                
#         fblock = sample(subpop[islargemiss], weights, blocksize;replace=false)            
#     else
#         0 < nlargemiss < blocksize || @error string("unexpected nlargemiss=",nlargemiss, ", blocksize=",blocksize)
#         b = .!islargemiss
#         weights = ProbabilityWeights(fmiss_subpop[b] .+ 1e-5)
#         fblock = vcat(subpop[islargemiss], sample(subpop[b], weights, blocksize-nlargemiss;replace=false))
#     end
#     fblock
# end