# ksilswap algorithm

#### Interface

# C is the type of centers, an (abstract) matrix of size (d x k)
# D is the type of pairwise distance computation from points to cluster centers
# WC is the type of cluster weights, either Int (in the case where points are
# unweighted) or eltype(weights) (in the case where points are weighted).
"""
    KsilResult{C,D<:Real,WC<:Real} 

The output of [`ksil`](@ref) or [`ksilswap`](@ref)

# Type parameters
 * `C<:AbstractMatrix{<:AbstractFloat}`: type of the `centers` matrix
 * `D<:Real`: type of the assignment silhouette
 * `WC<:Real`: type of the cluster weight
"""
struct KsilResult{C<:AbstractMatrix{<:AbstractFloat},D<:Real,WC<:Real} <: ClusteringResult
    centers::C                 # cluster centers (d x k)
    assignments::Vector{Int}   # assignments (n)
    silhouettes::Vector{D}     # silhouette of the assignments (n)
    counts::Vector{Int}        # number of points assigned to each cluster (k)
    wcounts::Vector{WC}        # cluster weights (k)
    totalsilhouette::D         # total silhouette (i.e. objective)    
    iteration::Int             # number of elapsed ksil iterations
    acceptedswap::DataFrame    # list of iterations in which swap was accepted
end

const _ksil_default_init = :kmpp
const _ksil_default_maxiter = 10
const _ksil_default_maxstuck = 2
const _ksil_default_tol = 1.0e-3

"""
ksilswap(X, k; weights, init, maxiter, maxstuck, tol, distance,rng) -> KsilResult

Kmeans-based ksil clustering of the ``d×n`` data matrix `X` (each column of `X`
is a ``d``-dimensional data point) into `k` clusters.

# Arguments

`maxstuck=100`: max number of sequential iterations without accepting random swap. 

See ksil for other keyargs. 

"""
function ksilswap(X::AbstractMatrix{<:Real},                # in: data matrix (d x n) columns = obs
                k::Integer;                               # in: number of centers
                weights::Union{Nothing, AbstractVector{<:Real}}=nothing, # in: data point weights (n)
                init::Union{Symbol, SeedingAlgorithm, AbstractVector{<:Integer}}=
                    _ksil_default_init,               # in: initialization algorithm
                maxiter::Integer=_ksil_default_maxiter, # in: maximum number of iterations                
                maxstuck::Int=_ksil_default_maxstuck,   # in: nuber of iterations without changes
                tol::Real=_ksil_default_tol,              # in: tolerance  of change at convergence
                distance::SemiMetric=SqEuclidean(),         # in: function to calculate distance with                
                rng::AbstractRNG=Random.GLOBAL_RNG)         # in: RNG object
    
    # ksil    
    centers, assignments, dmat, smat, silhouettes, counts, wcounts, to_update, unused = ksil_initialize(X,k; 
        weights,init,distance,rng)
    for _ in 1:3
        # counts,to_update, unused are alway re-initialzed in each update            
        ksil_update!(X, weights, centers, assignments, dmat, smat,silhouettes, counts, wcounts, 
            to_update, unused,distance,rng)
    end    

    # ksil        
    objv = weights === nothing ? sum(silhouettes) : dot(weights, silhouettes)    
    acceptedswap = DataFrame(iteration=[0],delcluster=[-1], newcluster=[-1],silhouette=[objv])    
    k == 1 && return KsilResult(centers, assignments, silhouettes, counts,
        wcounts, objv, 0, acceptedswap)    
    ksil_state = deepcopy([centers, assignments, dmat, smat,silhouettes, wcounts])
    t = 0
    nstuck = 0    
    while t < maxiter
        t += 1        
        acceptedswap2 = ksilswap_update!(X, weights, centers, assignments, dmat, smat,silhouettes, counts, wcounts, 
            to_update, unused, ksil_state, distance,rng)
        if isempty(acceptedswap2)
            nstuck += 1            
        else
            nstuck = 0
            acceptedswap2[!,1] .= t
            append!(acceptedswap, acceptedswap2)            
        end      
        for _ in 1:3                        
            ksil_update!(X, weights, centers, assignments, dmat, smat, silhouettes, counts, wcounts, 
                to_update, unused,distance,rng)        
        end
        ksil_state .= deepcopy([centers, assignments, dmat, smat, silhouettes, wcounts])  
        newobjv = weights === nothing ? sum(silhouettes) : dot(weights, silhouettes)    
        objv_change = objv - newobjv   
        objv = newobjv
        if objv_change > tol 
            @warn("The clustering silhouette decreased at ksil iteration #$t")
        end 
        abs(objv_change) < tol && nstuck >= maxstuck && break       
    end
    push!(acceptedswap, (t,-1, -1,round(objv,digits=4)))           
    n = size(X,2)
    acceptedswap[!,4] ./= n
    return KsilResult(centers, assignments, silhouettes, counts,
        wcounts, objv, t, acceptedswap)
end

# one inerations random-swap udpate
function ksilswap_update!(X::AbstractMatrix{<:Real},          # in: data matrix (d x n)
    weights::Union{Nothing, Vector{T}},       # in: data point weights (n)
    centers::AbstractMatrix{<:AbstractFloat}, # in/out: matrix of centers (d x k)
    assignments::Vector{Int},                 # in/out: assignment vector (n)       
    dmat::Matrix{<:Real},                     # in/out:  distance matrix (k x n)      
    smat::Matrix{<:Real},                     # in/out:  silhouette matrix (k x n)      
    silhouettes::Vector{<:Real},              # in/out: silhouettes of the resultant assignment (n)
    counts::Vector{Int},                      # in/out: # of points assigned to each cluster (k)
    wcounts::Vector{T},                       # in/out: updated cluster weights (k)
    to_update::Vector{Bool},                  # in/out: whether a center needs update (k)                             
    unused::Vector{Int},                      # in/out: list of centers with no points assigned        
    ksil_state::AbstractVector,               # in/out: store all variable value of ksil. 
    distance::SemiMetric,                     # in: function to calculate distance,    
                                              # in: initialization algorithm    
    rng::AbstractRNG) where T<:Real           # in: RNG object

    objv = weights === nothing ? sum(silhouettes) : dot(weights, silhouettes)    
    d, k = size(centers)        
    acceptedswap = DataFrame(iteration=Int[],delcluster=Int[], newcluster=Int[],silhouette=Float64[])    
    for oldc in 1:k 
        for newc in 1:k 
            # oldc == newc && continue
            # proposal for deleting a cluster and creating a new cluster
            newc_jls = findall(assignments .== newc)
            # newc_j = wsample(rng, newc_jls, 1.0 .- silhouettes[newc_jls])        
            newc_j = rand(newc_jls)
            centers[:,oldc] .= X[:,newc_j]
            fill!(to_update, false)
            to_update[oldc] = true

            # update distance matrix dmat            
            pairwise!(view(dmat, to_update, :), distance, view(centers, :, to_update), X, dims=2)
            # tused += @elapsed 
            ksil_update_smat!(smat, dmat)

            # three rounds of keans udpate
            for _ in 1:3                        
                ksil_update!(X, weights, centers, assignments, dmat, smat, silhouettes, counts, wcounts, 
                    to_update, unused,distance,rng)
            end        
            # To maximize objv by accept/reject        
            objv_prop = weights === nothing ? sum(silhouettes) : dot(weights, silhouettes)            
            # println("oldc=",oldc,",newc=",newc,",objv=",objv,",objv_prop=",objv_prop,",accept=",objv_prop>objv)
            if objv_prop > objv
                objv = objv_prop
                push!(acceptedswap, (0,oldc, newc,round(objv,digits=4)))            
                ksil_state .= deepcopy([centers, assignments, dmat, smat, silhouettes, wcounts])
            else            
                # change back 
                centers .= ksil_state[1]
                assignments .= ksil_state[2]
                dmat .= ksil_state[3]
                smat .= ksil_state[4]
                silhouettes .= ksil_state[5]
                wcounts .= ksil_state[6]
            end
        end
    end    
    # println("update smat, tused=",tused,"s")
    return acceptedswap
end    

"""
ksil(X, k; weights, init, maxiter, tol, distance,rng) -> KsilResult

Kmeans-based ksil clustering of the ``d×n`` data matrix `X` (each column of `X`
is a ``d``-dimensional data point) into `k` clusters.

"""
function ksil(X::AbstractMatrix{<:Real},                # in: data matrix (d x n) columns = obs
                k::Integer;                               # in: number of centers
                weights::Union{Nothing, AbstractVector{<:Real}}=nothing, # in: data point weights (n)
                init::Union{Symbol, SeedingAlgorithm, AbstractVector{<:Integer}}=
                    _ksil_default_init,               # in: initialization algorithm
                maxiter::Integer=_ksil_default_maxiter, # in: maximum number of iterations                                
                tol::Real=_ksil_default_tol,              # in: tolerance  of change at convergence
                distance::SemiMetric=SqEuclidean(),         # in: function to calculate distance with                
                rng::AbstractRNG=Random.GLOBAL_RNG)         # in: RNG object
    
    # ksil    
    centers, assignments, dmat, smat, silhouettes, counts, wcounts, to_update, unused = ksil_initialize(X,k; 
        weights,init,distance,rng)    
    objv = weights === nothing ? sum(silhouettes) : dot(weights, silhouettes)    
    t = 0    
    println("t=",t, ",objv=",objv)
    while t < maxiter
        t += 1       
        ksil_update!(X, weights, centers, assignments, dmat, smat, silhouettes, counts, wcounts, 
            to_update, unused,distance,rng)             
        newobjv = weights === nothing ? sum(silhouettes) : dot(weights, silhouettes)    
        println("t=",t, ",newobjv=",newobjv)
        objv_change = objv - newobjv   
        objv = newobjv
        if objv_change > tol 
            @warn("The clustering silhouette decreased at ksil iteration #$t")
        else
            abs(objv_change) < tol && break
        end        
    end    
    acceptedswap = DataFrame(iteration=Int[],delcluster=Int[], newcluster=Int[],silhouette=Float64[])
    return KsilResult(centers, assignments, silhouettes, counts,
        wcounts, objv, t, acceptedswap)
end

# ksil initialization
function ksil_initialize(X::AbstractMatrix{<:Real},           # in: data matrix (d x n)
    k::Integer;                                                 # in: number of centers
    weights::Union{Nothing, Vector{<:Real}},                   # in: data point weights (n)
    init::Union{Symbol, SeedingAlgorithm, AbstractVector{<:Integer}}, 
                                                                # in: initialization algorithm    
    distance::SemiMetric,                                       # in: function to calculate distance
    rng::AbstractRNG)                                           # in: RNG object

    d, n = size(X)
    (1 <= k <= n) || throw(ArgumentError("k must be from 1:n (n=$n), k=$k given."))

    # initialize the centers using a type wide enough so that the updates
    # centers[i, cj] += X[i, j] * wj will occur without loss of precision through rounding
    T = float(weights === nothing ? eltype(X) : promote_type(eltype(X), eltype(weights)))
    iseeds = Clustering.initseeds(init, X, k, rng=rng)
    centers = Clustering.copyseeds!(Matrix{T}(undef, d, k), X, iseeds)
    to_update = Vector(trues(k)) # whether a center needs to be updated
    unused = Vector{Int}()

    assignments = Vector{Int}(undef, n)
    counts = Vector{Int}(undef, k)  # number of data points assigned to each cluster    
    WC = (weights === nothing) ? Int : eltype(weights)
    wcounts = Vector{WC}(undef, k)
    
    dmat = pairwise(distance, centers, X, dims=2) # compute pairwise distances    
    D = typeof(one(eltype(dmat)) * one(WC))
    silhouettes = Vector{D}(undef, n)
    smat = similar(dmat)
    # update smat from dmat
    ksil_update_smat!(smat, dmat)
    ksil_update_assignments!(smat, true, assignments, silhouettes, counts, to_update, unused)

    centers, assignments, dmat, smat, silhouettes, counts, wcounts, to_update, unused    
end

function ksil_re_initialize!(counts::Vector{Int},  # in/out: # of points assigned to each cluster (k)    
    to_update::Vector{Bool},                  # in/out: whether a center needs update (k)                             
    unused::Vector{Int})                      # in/out: list of centers with no points assigned
    fill!(counts, 0)
    fill!(to_update, false)
    isempty(unused) || empty!(unused)
    return nothing
end

# one-iteration ksil update
function ksil_update!(X::AbstractMatrix{<:Real},          # in: data matrix (d x n)
                  weights::Union{Nothing, Vector{T}},       # in: data point weights (n)
                  centers::AbstractMatrix{<:AbstractFloat}, # in/out: matrix of centers (d x k)
                  assignments::Vector{Int},                 # in/out: assignment vector (n)            
                  dmat::Matrix{<:Real},                     # in/out:  distance matrix (k x n)      
                  smat::Matrix{<:Real},                     # in/out:  silhouette matrix (k x n)      
                  silhouettes::Vector{<:Real},                    # in/out: silhouettes of the resultant assignment (n)
                  counts::Vector{Int},                      # in/out: # of points assigned to each cluster (k)
                  wcounts::Vector{T},                       # in/out: updated cluster weights (k)
                  to_update::Vector{Bool},                  # in/out: whether a center needs update (k)                             
                  unused::Vector{Int},                       # in/out: list of centers with no points assigned
                  distance::SemiMetric,                     # in: function to calculate distance
                  rng::AbstractRNG) where T<:Real           # in: RNG object    
    
    
    # re_initialization and update assignments        
    ksil_re_initialize!(counts, to_update, unused)
    ksil_update_assignments!(smat, false, assignments, silhouettes, counts, to_update, unused)

    # update (affected) centers    
    ksil_update_centers!(X, weights, assignments, to_update, centers, wcounts)
    if !isempty(unused)
        ksil_update_unused_centers!(X, silhouettes, centers, unused, distance, rng)
        to_update[unused] .= true
    end    
    if sum(to_update) > 0.75 * size(centers,2)
        pairwise!(dmat, distance, centers, X, dims=2)
    else
        # if only a small subset is affected, only compute for that subset        
        pairwise!(view(dmat, to_update, :), distance,
                    view(centers, :, to_update), X, dims=2)
    end     
    # update smat from dmat
    ksil_update_smat!(smat, dmat)
    # return
    centers, assignments, dmat, smat, silhouettes, counts, wcounts, to_update, unused
end

#
# Calculate simplified sihouette: see https://en.wikipedia.org/wiki/Silhouette_(clustering)
#
function ksil_update_smat!(smat::Matrix{<:Real},     # in/out:  silhouette matrix (k x n)
    dmat::Matrix{<:Real},                            # in:  silhouette matrix (k x n)
    )
    k,n = size(smat)            
    @inbounds for j = 1:n
        dls = view(dmat, :,j)
        i1,i2 = partialsortperm(dls,1:2)
        for i = 1:k
            a = dls[i]
            b = i==i1 ? dls[i2] : dls[i1] 
            smat[i,j] = (b-a)/max(a,b)
        end
    end
end    

#
# The following is modefied from Clustering.jl/kmeans.jl
#
function ksil_update_assignments!(smat::Matrix{<:Real},     # in:  silhouette matrix (k x n)
                             is_init::Bool,            # in:  whether it is the initial run
                             assignments::Vector{Int}, # out: assignment vector (n)
                             silhouettes::Vector{<:Real},    # out: silhouettes of the resultant assignment (n)
                             counts::Vector{Int},      # out: # of points assigned to each cluster (k)
                             to_update::Vector{Bool},  # out: whether a center needs update (k)
                             unused::Vector{Int}       # out: list of centers with no points assigned
                             )
    k, n = size(smat)
    # process each point
    @inbounds for j = 1:n
        # find the closest cluster to the i-th point. Note that a
        # is necessarily between 1 and size(smat, 1) === k as a result
        # and can thus be used as an index in an `inbounds` environment
        c, a = findmax(view(smat, :, j))

        # set/update the assignment
        if is_init
            assignments[j] = a
        else  # update
            pa = assignments[j]
            if pa != a
                # if assignment changes,
                # both old and new centers need to be updated
                assignments[j] = a
                to_update[a] = true
                to_update[pa] = true
            end
        end

        # set silhouettes and counts accordingly
        silhouettes[j] = c
        counts[a] += 1
    end

    # look for centers that have no assigned points
    for i = 1:k
        if counts[i] == 0
            push!(unused, i)
            to_update[i] = false # this is handled using different mechanism
        end
    end
end

#
#  Update centers based on updated assignments
#
#  (specific to the case where points are not weighted)
#
function ksil_update_centers!(X::AbstractMatrix{<:Real},        # in: data matrix (d x n)
                         weights::Nothing,                 # in: point weights
                         assignments::Vector{Int},         # in: assignments (n)
                         to_update::Vector{Bool},          # in: whether a center needs update (k)
                         centers::AbstractMatrix{<:AbstractFloat}, # out: updated centers (d x k)
                         wcounts::Vector{Int})             # out: updated cluster weights (k)
    d, n = size(X)
    k = size(centers, 2)

    # initialize center weights
    wcounts[to_update] .= 0

    # accumulate columns
    @inbounds for j in 1:n
        # skip points assigned to a center that doesn't need to be updated
        cj = assignments[j]
        if to_update[cj]
            if wcounts[cj] > 0
                for i in 1:d
                    centers[i, cj] += X[i, j]
                end
            else
                for i in 1:d
                    centers[i, cj] = X[i, j]
                end
            end
            wcounts[cj] += 1
        end
    end

    # sum ==> mean
    @inbounds for j in 1:k
        if to_update[j]
            cj = wcounts[j]
            for i in 1:d
                centers[i, j] /= cj
            end
        end
    end
end

#
#  Update centers based on updated assignments
#
#  (specific to the case where points are weighted)
#
function ksil_update_centers!(X::AbstractMatrix{<:Real}, # in: data matrix (d x n)
                         weights::Vector{W},        # in: point weights (n)
                         assignments::Vector{Int},  # in: assignments (n)
                         to_update::Vector{Bool},   # in: whether a center needs update (k)
                         centers::AbstractMatrix{<:Real}, # out: updated centers (d x k)
                         wcounts::Vector{W}         # out: updated cluster weights (k)
                         ) where W<:Real
    d, n = size(X)
    k = size(centers, 2)

    # initialize center weights
    wcounts[to_update] .= 0

    # accumulate columns
    @inbounds for j in 1:n
        # skip points with negative weights or assigned to a center
        # that doesn't need to be updated
        wj = weights[j]
        cj = assignments[j]
        if wj > 0 && to_update[cj]
            if wcounts[cj] > 0
                for i in 1:d
                    centers[i, cj] += X[i, j] * wj
                end
            else
                for i in 1:d
                    centers[i, cj] = X[i, j] * wj
                end
            end
            wcounts[cj] += wj
        end
    end

    # sum ==> mean
    @inbounds for j in 1:k
        if to_update[j]
            cj = wcounts[j]
            for i in 1:d
                centers[i, j] /= cj
            end
        end
    end
end


#
#  Re-picks centers that have no points assigned to them.
#
function ksil_update_unused_centers!(X::AbstractMatrix{<:Real}, # in: the data matrix (d x n)
                               silhouettes::Vector{<:Real},     # in: the current assignment silhouettes (n)
                               centers::AbstractMatrix{<:AbstractFloat}, # out: the centers (d x k)
                               unused::Vector{Int},       # in: indices of centers to be updated
                               distance::SemiMetric,      # in: function to calculate the distance with
                               rng::AbstractRNG)          # in: RNG object
    # pick new centers using a scheme like ksil++
    ds = similar(silhouettes)
    tcosts = eltype(silhouettes)(1) .- silhouettes
    n = size(X, 2)

    for i in unused
        j = wsample(rng, 1:n, tcosts)
        tcosts[j] = 0
        v = view(X, :, j)
        centers[:, i] = v
        colwise!(ds, distance, v, X)
        tcosts = min(tcosts, ds)
    end
end

