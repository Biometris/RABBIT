# using Distributed
# nprocs() < 4 && addprocs(4-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using SpectralOrdering

using SpectralEmbedding
cd(@__DIR__)
pwd()


using Distributions, SparseArrays, Random, LinearAlgebra
using Plots
using Distances, Clustering

# Ctrl-S pauses the REPL terminal, so you can either put that keybinding
# into the whitelisted keybindings setting (under julia-client) or
# press Ctrl-Q to resume.

xy = [begin
    angle = rand(Uniform(0,2pi),1000)
    r = rand(Uniform(s-0.5,s+0.5),1000)
    r .* [cos.(angle) sin.(angle)]
end for s = [10,4,1]]
points = collect(eachrow(vcat(xy...)))
shuffle!(points)
xs = [i[1] for i=points]
ys = [i[2] for i=points]
scatter(xs,ys)

similarity = sparse([round(exp(-norm(i .- j)),digits=5) for i = points, j=points])
# heatmap(Matrix(similarity))

ncluster = 3
knn, knn_history = SpectralEmbedding.findknn(similarity; 
    ncomponent=1,
    knnmin = round(Int,sqrt(size(similarity,1))),
    mincomponentsize=1,
    io=nothing, verbose=true)
similarity2= SpectralEmbedding.toknnsimilarity(similarity,knn);
similarity2, connectednodes, singletons = SpectralEmbedding.takecomponents(similarity2; 
    ncomponent=1
);



ncluster = 3
ncls = 1:2ncluster
@time clusters, vals, vecs, sils, res_silhouettes = spectralclustering(similarity2,ncls;    
    clusteralg = "kmeans", 
    # verbose=true,
);


using Plots
SpectralEmbedding.summary_silhouette(res_silhouettes)
SpectralEmbedding.plot_silhouette(res_silhouettes)

keys(res_silhouettes)


silhouettes = res_silhouettes[:silhsls][3]
clusters = res_silhouettes[:clustersls][3]


scatter(silhouettes)

propertynames(res_silhouettes)
clusters = res_silhouettes.clustersls[3];
ncluster = length(clusters)
scatter(vecs[:,1:ncluster])
plot(plot(vecs[:,1:ncluster]),scatter(vals))
xys = tuple.(xs,ys)
scatter(xys[clusters[1]])
scatter!(xys[clusters[2]])
scatter!(xys[clusters[3]])



# @time vals, vecs = SpectralEmbedding.spectralembedding(similarity2; 
#     nev = 6, laplaciantype="randomwalk"
# )        

@time clusters, vals, vecs, res_silhouettes = spectralclustering(similarity2,3;
    # verbose=true,
);

