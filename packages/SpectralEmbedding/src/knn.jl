function dropsingletons(similarity::AbstractMatrix;
    mincomponentsize::Integer=2)
    cc = get_connected_components(similarity)
    b = length.(cc) .<  mincomponentsize
    singletons = sort(vcat(cc[b]...))
    connectednodes = setdiff(1:size(similarity,1),singletons)
    similarity[connectednodes, connectednodes], connectednodes, singletons
end

function takecomponents(similarity::AbstractMatrix;
    ncomponent::Integer=1)
    cc = get_connected_components(similarity)
    p = sortperm(length.(cc),rev=true)
    ncomponent2 = min(ncomponent, length(p))
    connectednodes = vcat(cc[p[1:ncomponent2]]...)
    singletons = setdiff(1:size(similarity,1),connectednodes)
    similarity[connectednodes, connectednodes], connectednodes, singletons
end


function findknn(similarity::AbstractMatrix;
    ncomponent::Union{Nothing,Integer}=nothing,
    knnmin::Integer = round(Int,sqrt(size(similarity,1))),
    mincomponentsize::Integer= 2,
    isordered::Bool=false,
    alwayskeep::Union{Nothing,Real}=nothing,
    io::Union{Nothing,IO} = nothing,
    verbose::Bool=false)
    if isnothing(ncomponent)
        lenls = length.(get_connected_components(similarity))
        ncomponent = sum(lenls .>= mincomponentsize)
        if ncomponent > 1
            @info string("connected component sizes =",sort(lenls[lenls .>= mincomponentsize],rev=true),
                " after removing size < ", mincomponentsize)
        end
    end
    knnmax = size(similarity,1)-1
    his = Vector{Vector{Int}}()
    knn = knnmin
    while true
        cc = first(SpectralEmbedding.knn_connected_components(similarity, knn, isordered;alwayskeep))
        nc = sum(length.(cc) .>=  mincomponentsize)
        push!(his, [knn, nc])
        if nc <= ncomponent || knn == knnmax
            break
        end
        knn = min(2*knn, knnmax)
    end
    printconsole(io,verbose,string("knn_history=",his))
    his[1][2] <= ncomponent && return (knn, his)
    his[end][2] > ncomponent && return (his[end][1], his)
    ressim = similarity
    while true
        knn = min(2*his[end-1][1],knnmax)
        if knn >= his[end][1]
            knn = floor(Int,(his[end-1][1]+his[end][1])/2)
        end
        cc, sim = knn_connected_components(ressim, knn, isordered;alwayskeep)
        nc = sum(length.(cc) .>=  mincomponentsize)
        if nc == his[end][2] || nc <= ncomponent
            # push!(his, [knn,nc])
            his[end] = [knn,nc]
            ressim = sim
        else
            insert!(his, length(his),[knn,nc])
        end
        printconsole(io,verbose,string("knn_history=",his))
        if his[end][1] == his[end-1][1]+1
            break
        end
    end
    his[end][1], his
end

function printconsole(io::Union{Nothing,IO},verbose::Bool,msg::AbstractString)
    verbose && @info(msg)
    if !isnothing(io)
        write(io,string(msg,"\n"))
        flush(io)
    end
end

function knn_connected_components(similarity::AbstractMatrix,
    knn::Integer, isordered::Bool;
    alwayskeep::Union{Nothing,Real}=nothing)
    sim = SpectralEmbedding.toknnsimilarity(similarity,knn;alwayskeep)
    aa = sign.(sim)
    if isordered
        bb = triu(aa,1) .* tril(aa, knn+1)
        bb .+= bb'
    else
        bb = aa
    end
    get_connected_components(bb), sim
end

function get_connected_components(similarity::AbstractMatrix)
    g = SimpleGraph(sign.(similarity))
    connected_components(g)
end

function toknnsimilarity(similarity::SparseMatrixCSC,knn::Integer;
    alwayskeep::Union{Nothing,Real}=nothing)
    iils,jjls, vvls = findnz(similarity)
    df = DataFrame(i=iils,j=jjls,v=vvls)
    gd = groupby(df,:i)    
    df = reduce(vcat, [begin 
        di = sort(gd[i],:v,rev=true)        
        if knn < size(di,1)
            if isnothing(alwayskeep)
                di = di[1:knn,:]
            else
                lastkeep = findlast(x -> x >= alwayskeep, di[!,:v])
                if isnothing(lastkeep)
                    di = di[1:knn,:]
                else                    
                    lastcol = max(knn, lastkeep)
                    di = di[1:lastcol,:]
                end
            end
        end
        di
    end for i in eachindex(gd)])    
    dftran = DataFrame(i=df[!,2],j=df[!,1],v=df[!,3])
    df = vcat(df, dftran)
    gd = groupby(df,[:i,:j])
    df2 = combine(gd, :v => (x->max(x...)) => :vmax)
    sparse(df2[!,:i],df2[!,:j],df2[!,:vmax],size(similarity)...)
end

function toknnsimilarity2(similarity::SparseMatrixCSC,knn::Integer)
    nn = size(similarity,1)
    iils = Vector{Int}()
    jjls = Vector{Int}()
    vvls = Vector{eltype(similarity)}()
    for i = 1:nn
        jj, vv = findnz(similarity[i,:])
        k = min(length(vv),knn)
        pos = partialsortperm(vv,1:k,rev=true)
        newjj = jj[pos]
        newvv = vv[pos]
        newii= repeat([i],length(newjj))
        append!(iils, newii)
        append!(jjls, newjj)
        append!(vvls, newvv)
    end
    df = DataFrame(i=vcat(iils, jjls), j=vcat(jjls, iils), v=vcat(vvls, vvls))
    gd = groupby(df,[:i,:j])
    df2 = combine(gd, :v => (x->max(x...)) => :vmax)
    sparse(df2[!,:i],df2[!,:j],df2[!,:vmax],nn,nn)
end

function toknnsimilarity3(similarity::SparseMatrixCSC,knn::Integer)
    @assert 2<=knn<=size(similarity,1)
    res = spzeros(size(similarity)...)
    for i = 1:size(res,1)
        jj, vv = findnz(similarity[i,:])
        k = min(length(vv),knn)
        pos = partialsortperm(vv,1:k,rev=true)
        pos = jj[pos]
        res[i,pos] .= similarity[i,pos]
    end
    sparse(map(max,res,res'))
end

function toknnsimilarity(similarity::Matrix{T},knn::Integer) where {T<:Real}
    @assert 2<=knn<=size(similarity,1)
    res = zeros(size(similarity)...)
    for i = 1:size(res,1)
        pos = partialsortperm(similarity[i,:],1:knn,rev=true)
        res[i,pos] .= similarity[i,pos]
    end
    map(max,res,res')
end
