module SpectralEmbedding

using Distributed
using StatsBase, LinearAlgebra, SparseArrays
using Arpack
using Distances  # calculating pariwise distance for hclus
using Clustering
using LightGraphs, DataFrames
using Plots
using Random

export    
    # plot_silhouette, plot_eigen,    
    ksil_initialize, ksil_update!,ksil_re_initialize!,ksil_update_assignments!,
    spectralembedding, randswap,
    knn_spectralclustering, spectralclustering, 
    knn_spectralordering, spectralordering

include("randswap.jl")    
include("ksilswap.jl")    
include("knn.jl")
include("spectralembedding_fun.jl")
include("silhouette.jl")
include("spectralclustering.jl")
include("spectralordering.jl")

end # module SpectralEmbedding
