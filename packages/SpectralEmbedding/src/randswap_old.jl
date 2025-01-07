# Randswap algorithm

#### Interface

# C is the type of centers, an (abstract) matrix of size (d x k)
# D is the type of pairwise distance computation from points to cluster centers
# WC is the type of cluster weights, either Int (in the case where points are
# unweighted) or eltype(weights) (in the case where points are weighted).
"""
    RandswapResult{C,D<:Real,WC<:Real} <: ClusteringResult

The output of [`kmeans`](@ref) and [`kmeans!`](@ref).

# Type parameters
 * `C<:AbstractMatrix{<:AbstractFloat}`: type of the `centers` matrix
 * `D<:Real`: type of the assignment cost
 * `WC<:Real`: type of the cluster weight
"""
struct RandswapResult{C<:AbstractMatrix{<:AbstractFloat},D<:Real,WC<:Real} <: ClusteringResult
    centers::C                 # cluster centers (d x k)
    assignments::Vector{Int}   # assignments (n)
    costs::Vector{D}           # cost of the assignments (n)
    counts::Vector{Int}        # number of points assigned to each cluster (k)
    wcounts::Vector{WC}        # cluster weights (k)
    totalcost::D               # total cost (i.e. objective)    
    iteration::Int             # number of elapsed randswap iterations
    accepted_iterations::Vector{Int}   # list of iterations in which random swap was accepted
end

const _randswap_default_init = :kmpp
const _randswap_default_maxiter = 10000
const _randswap_default_maxstuck = 2000
const _randswap_default_tol = 1.0e-6

"""
    randswap(X, k; weights, init, maxiter, maxstuck, tol, distance,rng) -> RandswapResult

Kmeans-based randswap clustering of the ``d×n`` data matrix `X` (each column of `X`
is a ``d``-dimensional data point) into `k` clusters.

# Arguments

`maxstuck=1000`: max number of sequential iterations without accepting random swap. 

See kmeans for other keyargs. 

"""
function randswap(X::AbstractMatrix{<:Real},                # in: data matrix (d x n) columns = obs
                k::Integer;                               # in: number of centers
                weights::Union{Nothing, AbstractVector{<:Real}}=nothing, # in: data point weights (n)
                init::Union{Symbol, SeedingAlgorithm, AbstractVector{<:Integer}}=
                    _randswap_default_init,               # in: initialization algorithm
                maxiter::Integer=_randswap_default_maxiter, # in: maximum number of iterations                
                maxstuck::Int=_randswap_default_maxstuck,   # in: nuber of iterations without changes
                tol::Real=_randswap_default_tol,              # in: tolerance  of change at convergence
                distance::SemiMetric=SqEuclidean(),         # in: function to calculate distance with
                rng::AbstractRNG=Random.GLOBAL_RNG)         # in: RNG object
    
    # kmeans    
    centers, assignments, dmat, costs, counts, wcounts, to_update, unused = kmeans_initialize(X,k; 
        weights,init,distance,rng)
    for _ in 1:3
        # counts,to_update, unused are alway re-initialzed in each update            
        kmeans_update!(X, weights, centers, assignments, dmat, costs, counts, wcounts, 
            to_update, unused,distance,rng)
    end    

    # randswap    
    oldobjv = weights === nothing ? sum(costs) : dot(weights, costs)    
    accepted_iterations = Vector{Int}()
    k == 1 && return RandswapResult(centers, assignments, costs, counts,
        wcounts, oldobjv, 0, accepted_iterations)
    randswap_state = deepcopy([centers, assignments, dmat, costs, wcounts])
    t = 0
    nstuck = 0    
    while t < maxiter
        t += 1
        isaccept, objv = randswap_update!(X, weights, centers, assignments, dmat, costs, counts, wcounts, 
            to_update, unused, oldobjv, randswap_state, distance,rng)
        if isaccept 
            nstuck = 0
            push!(accepted_iterations, t)
        else
            nstuck += 1
        end        
        objv_change = objv - oldobjv
        oldobjv = objv
        if objv_change > tol 
            @warn("The clustering cost increased at randswap iteration #$t")
        else
            abs(objv_change) < tol && nstuck >= maxstuck && break
        end        
    end
    return RandswapResult(centers, assignments, costs, counts,
        wcounts, oldobjv, t, accepted_iterations)
end

# one inerations random-swap udpate
function randswap_update!(X::AbstractMatrix{<:Real},          # in: data matrix (d x n)
    weights::Union{Nothing, Vector{T}},       # in: data point weights (n)
    centers::AbstractMatrix{<:AbstractFloat}, # in/out: matrix of centers (d x k)
    assignments::Vector{Int},                 # in/out: assignment vector (n)       
    dmat::Matrix{<:Real},                     # in/out:  distance matrix (k x n)      
    costs::Vector{<:Real},                    # in/out: costs of the resultant assignment (n)
    counts::Vector{Int},                      # in/out: # of points assigned to each cluster (k)
    wcounts::Vector{T},                       # in/out: updated cluster weights (k)
    to_update::Vector{Bool},                  # in/out: whether a center needs update (k)                             
    unused::Vector{Int},                      # in/out: list of centers with no points assigned    
    objv::Real,                               # in: object target function value
    randswap_state::AbstractVector,           # in/out: store all variable value of randswap. 
    distance::SemiMetric,                     # in: function to calculate distance
    rng::AbstractRNG) where T<:Real           # in: RNG object

    # random swap prototype, and keep current state for recovering
    delcluster = randswap_prototype!(X,centers, rng)    

    # update distance matrix dmat
    fill!(to_update, false)
    to_update[delcluster] = true
    pairwise!(view(dmat, to_update, :), distance, view(centers, :, to_update), X, dims=2)

    # three rounds of keans udpate
    for _ in 1:3                        
        kmeans_update!(X, weights, centers, assignments, dmat, costs, counts, wcounts, 
            to_update, unused,distance,rng)
    end

    #  check if accept 
    # To minimize objv
    objv_prop = weights === nothing ? sum(costs) : dot(weights, costs)    
    if objv_prop < objv
        isaccept = true
        objv = objv_prop
        randswap_state .= deepcopy([centers, assignments, dmat, costs, wcounts])
    else
        isaccept = false
        # change back 
        centers .= randswap_state[1]
        assignments .= randswap_state[2]
        dmat .= randswap_state[3]
        costs .= randswap_state[4]
        wcounts .= randswap_state[5]
    end
    return isaccept, objv
end    

# random swap prototype
function randswap_prototype!(X::AbstractMatrix{<:Real},                # in: data matrix (d x n)
    centers::AbstractMatrix{<:AbstractFloat},                 # in/out: matrix of centers (d x k)        
    rng::AbstractRNG)                                         # in: RNG object    
    delcluster = rand(rng,1:size(centers,2))
    prop_Xindex = rand(rng,1:size(X,2))    
    centers[:,delcluster] .= X[:,prop_Xindex]
    return delcluster
end

# kmeans initialization
function kmeans_initialize(X::AbstractMatrix{<:Real},           # in: data matrix (d x n)
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
    costs = Vector{D}(undef, n)
    update_assignments!(dmat, true, assignments, costs, counts, to_update, unused)

    centers, assignments, dmat, costs, counts, wcounts, to_update, unused    
end

function re_initialize!(counts::Vector{Int},  # in/out: # of points assigned to each cluster (k)    
    to_update::Vector{Bool},                  # in/out: whether a center needs update (k)                             
    unused::Vector{Int})                      # in/out: list of centers with no points assigned
    fill!(counts, 0)
    fill!(to_update, false)
    isempty(unused) || empty!(unused)
    return nothing
end

# one-iteration kmeans update
function kmeans_update!(X::AbstractMatrix{<:Real},          # in: data matrix (d x n)
                  weights::Union{Nothing, Vector{T}},       # in: data point weights (n)
                  centers::AbstractMatrix{<:AbstractFloat}, # in/out: matrix of centers (d x k)
                  assignments::Vector{Int},                 # in/out: assignment vector (n)            
                  dmat::Matrix{<:Real},                     # in/out:  distance matrix (k x n)      
                  costs::Vector{<:Real},                    # in/out: costs of the resultant assignment (n)
                  counts::Vector{Int},                      # in/out: # of points assigned to each cluster (k)
                  wcounts::Vector{T},                       # in/out: updated cluster weights (k)
                  to_update::Vector{Bool},                  # in/out: whether a center needs update (k)                             
                  unused::Vector{Int},                       # in/out: list of centers with no points assigned
                  distance::SemiMetric,                     # in: function to calculate distance
                  rng::AbstractRNG) where T<:Real           # in: RNG object    
    
    
    # re_initialization and update assignments        
    re_initialize!(counts, to_update, unused)
    update_assignments!(dmat, false, assignments, costs, counts, to_update, unused)

    # update (affected) centers    
    update_centers!(X, weights, assignments, to_update, centers, wcounts)
    if !isempty(unused)
        repick_unused_centers!(X, costs, centers, unused, distance, rng)
        to_update[unused] .= true
    end    
    if sum(to_update) > 0.75 * size(centers,2)
        pairwise!(dmat, distance, centers, X, dims=2)
    else
        # if only a small subset is affected, only compute for that subset        
        pairwise!(view(dmat, to_update, :), distance,
                    view(centers, :, to_update), X, dims=2)
    end        
    # return
    centers, assignments, dmat, costs, counts, wcounts, to_update, unused
end


# The folowing is copied from Clustering.jl/kmeans.jl 
#
#  Updates assignments, costs, and counts based on
#  an updated (squared) distance matrix
#
function update_assignments!(dmat::Matrix{<:Real},     # in:  distance matrix (k x n)
                             is_init::Bool,            # in:  whether it is the initial run
                             assignments::Vector{Int}, # out: assignment vector (n)
                             costs::Vector{<:Real},    # out: costs of the resultant assignment (n)
                             counts::Vector{Int},      # out: # of points assigned to each cluster (k)
                             to_update::Vector{Bool},  # out: whether a center needs update (k)
                             unused::Vector{Int}       # out: list of centers with no points assigned
                             )
    k, n = size(dmat)
    # process each point
    @inbounds for j = 1:n
        # find the closest cluster to the i-th point. Note that a
        # is necessarily between 1 and size(dmat, 1) === k as a result
        # and can thus be used as an index in an `inbounds` environment
        c, a = findmin(view(dmat, :, j))

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

        # set costs and counts accordingly
        costs[j] = c
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
function update_centers!(X::AbstractMatrix{<:Real},        # in: data matrix (d x n)
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
function update_centers!(X::AbstractMatrix{<:Real}, # in: data matrix (d x n)
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
function repick_unused_centers!(X::AbstractMatrix{<:Real}, # in: the data matrix (d x n)
                               costs::Vector{<:Real},     # in: the current assignment costs (n)
                               centers::AbstractMatrix{<:AbstractFloat}, # out: the centers (d x k)
                               unused::Vector{Int},       # in: indices of centers to be updated
                               distance::SemiMetric,      # in: function to calculate the distance with
                               rng::AbstractRNG)          # in: RNG object
    # pick new centers using a scheme like kmeans++
    ds = similar(costs)
    tcosts = copy(costs)
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

