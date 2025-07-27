function mapaccuracy(mapxfile::AbstractString,mapyfile::AbstractString;    
    isphysmap::AbstractVector=[false,false],
    missingstring=["NA","missing"],    
    commentstring::AbstractString="##",
    isgroupacc::Bool=true, 
    workdir::AbstractString = pwd())    
    mapfilels = [mapxfile,mapyfile]
    length(isphysmap) == 2 || @error string("isphysmap must be a vector of length 2")
    mapx, mapy = [begin
        col = isphysmap[i] ? :physchrom : :linkagegroup
        mapdf = readmarkermap(mapfilels[i]; del_ungrouped=true, commentstring,missingstring, workdir)
        DataFrame.(collect(groupby(mapdf,col)))
    end for i in eachindex(mapfilels)]    
    mapaccuracy(mapx,mapy; isgroupacc,isphysmap)
end

function mapaccuracy(mapx::Vector{DataFrame}, mapy::Vector{DataFrame}; 
    isgroupacc=true,isphysmap::AbstractVector=[false,false])    
    groupres = grouping_accuracy(mapx,mapy; isgroupacc)
    orderacc = mapcorkendall(mapx,mapy; isphysmap)
    map1len = maplength(mapx; isphysmap=isphysmap[1])
    map2len = maplength(mapy; isphysmap=isphysmap[2])
    map1nsnp = sum(size.(mapx,1))
    map2nsnp = sum(size.(mapy,1))
    snpx = reduce(vcat, [i[!,:marker] for i in mapx])
    snpy = reduce(vcat, [i[!,:marker] for i in mapy])
    snpdiff12 = setdiff(snpx,snpy)
    snpdiff21 = setdiff(snpy,snpx)
    (groupres..., orderacc=orderacc,
        snpdiff12 = snpdiff12, snpdiff21 = snpdiff21,
        map1len = map1len, map2len = map2len,
        map1nsnp = map1nsnp, map2nsnp=map2nsnp)
end

function maplength(markermap::Vector{DataFrame}; isphysmap::Bool)
    col = isphysmap ? :physposbp : :poscm
    [i[end,col]-i[1,col] for i in markermap]
end

function grouping_accuracy(mapx::Vector{DataFrame}, mapy::Vector{DataFrame}; isgroupacc)
    mapylg = MagicBase.findtruelg(mapy,mapx;verbose=false)    
    mapy = [reduce(vcat,mapy[i]) for i in mapylg]    #
    snpx = [i[!,:marker] for i in mapx]
    snpy = [isempty(i) ? [] : i[!,:marker] for i in mapy]
    snpls = intersect(reduce(vcat,snpx),reduce(vcat,snpy))
    dict = Dict(snpls .=> 1:length(snpls))
    snpx2 = [filter(x->!isnothing(x),[get(dict,j,nothing) for j in i]) for i in snpx]
    snpy2 = [filter(x->!isnothing(x),[get(dict,j,nothing) for j in i]) for i in snpy]    
    nsnpinconsist = sum(map((x,y)->length(setdiff(x,y)),snpx2,snpy2))
    nsnpcommon = length(dict)
    if isgroupacc
        groupacc = paircounting_F1score(snpx2,snpy2)
    else
        groupacc = nothing
    end
    (nsnpinconsist=nsnpinconsist,nsnpcommon=nsnpcommon, groupacc=groupacc)
end

function paircounting_F1score(partition1::AbstractVector, partition2::AbstractVector)
    (partition1 == partition2) && return 1.0
    comb1 = [combinations(sort(i),2) for i in partition1]
    comb2 = [combinations(sort(i),2) for i in partition2]
    a = length(reduce(vcat,[reduce(vcat, [intersect(i,j) for j in comb2]) for i in comb1]))
    b = length(reduce(vcat,[reduce(intersect,[setdiff(i,j) for j in comb2]) for i in comb1]))
    c = length(reduce(vcat,[reduce(intersect,[setdiff(i,j) for j in comb1]) for i in comb2]))
    2a/(2a + b + c)
end

function mapcorkendall(mapx, mapy; isphysmap = [false,false], minfreq = 0.3)
    mapylg = MagicBase.findtruelg(mapy,mapx; minfreq, isphysmap = reverse(isphysmap))    
    mapy2 = [reduce(vcat,mapy[i]) for i in mapylg]
    map(MagicBase.calcorkendall, mapx, mapy2)
end

# function mapcorkendall(mapx::Vector{DataFrame}, mapy::Vector{DataFrame})    
#     allmapy = reduce(vcat,mapy)
#     [calcorkendall(df, allmapy) for df in mapx]    
# end

function calcorkendall(chrmapx::AbstractDataFrame, chrmapy::AbstractDataFrame)
    (isempty(chrmapx) || isempty(chrmapy)) && return 0.0
    commonsnps = intersect(chrmapx[!,:marker],chrmapy[!,:marker])
    df=filter(x->in(x[1],commonsnps),chrmapx)
    rule = Dict(chrmapy[!,:marker] .=> 1:size(chrmapy,1))
    pos = [get(rule,i,missing) for i=df[!,:marker]]
    isnonmiss= .!ismissing.(pos)
    posnon = Vector{Int}(pos[isnonmiss])
    corkendall(posnon,collect(1:length(posnon)))
end
